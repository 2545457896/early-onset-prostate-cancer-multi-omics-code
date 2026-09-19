#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01_run_infercnv_Tref.py
=======================
PCa 上皮细胞 inferCNV —— 全样本合并 + T 细胞参考


运行方式:
  /root/miniconda3/envs/scanpy2/bin/python 01_run_infercnv_Tref.py
"""

import os
import sys
import glob
import json as _json
import random
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import infercnvpy as cnv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from datetime import datetime

today = datetime.now().date().strftime("%Y%m%d")

# ======================== 随机种子 ========================
SEED = 20260815
random.seed(SEED)
np.random.seed(SEED)

# ======================== 参数配置 ========================
WD = "/mnt/02.PCa-zhu/01.code"

# 输入
_h5ad_files = sorted(glob.glob(os.path.join(WD, "01_input", "PCa_T_epi_merged_*.h5ad")))
assert len(_h5ad_files) > 0, "01_input 下未找到 PCa_T_epi_merged_*.h5ad, 请先运行 01_export_T_epi_h5ad.R"
H5AD_PATH = _h5ad_files[-1]
GENEINFO_PATH = os.path.join(WD, "01_input", "infercnvpy_geneInfor_20240606.csv")

# 输出
OUT_DIR = os.path.join(WD, "01_output")
os.makedirs(OUT_DIR, exist_ok=True)

# ---- inferCNV 参数 ----
REFERENCE_KEY = "celltype"
REFERENCE_CAT = ["T_cells"]
WINDOW_SIZES = [100, 200]
INFERCNV_N_JOBS = 4
INFERCNV_CHUNKSIZE = 5000

# ---- 肿瘤簇判定参数  ----
NORMAL_LABEL = "normal"   # lesion 列中正常来源标签 (对应 RCC 的 NR)
NR_MAX = 0.25             # 条件①: 正常来源占比上限
SD_MULT = 1.0             # 条件②: 超基线 s.d. 倍数 (基线 = T_cells)

# ---- 染色体热图子采样 ----
HEATMAP_MAX_PER_CLUST = 300

PREFIX = "PCa_Epithelial_Tref"

print("=" * 70)
print(f"[{today}] 01: PCa Epithelial inferCNV (T cell reference, merged)")
print("=" * 70)
print(f"  h5ad:     {H5AD_PATH}")
print(f"  geneInfor:{GENEINFO_PATH}")
print(f"  results:  {OUT_DIR}")
print(f"  window:   {WINDOW_SIZES}")
print(f"  reference:{REFERENCE_CAT}")
print(f"  判定: normal_frac < {NR_MAX:.0%} AND mean cnv_score >= T_base + {SD_MULT:.0f} SD")

# ======================== Step 1: 读取 h5ad ========================
print("\n" + "=" * 70)
print("Step 1: 读取合并 h5ad")
print("=" * 70)

adata = sc.read_h5ad(H5AD_PATH)
print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

# X 必须为 raw counts (整数)
_X = adata.X
_Xd = _X.data if hasattr(_X, "data") else np.asarray(_X).ravel()
_frac_int = np.mean(_Xd == np.round(_Xd))
print(f"  X 整数比例: {_frac_int:.4f} | max = {_Xd.max()}")
assert _frac_int > 0.999, "X 疑似非 raw counts, 请检查 01_export_T_epi_h5ad.R"

print(f"  obs columns: {list(adata.obs.columns)}")
for _c in ["celltype", "cell_source", "tissue", "lesion"]:
    assert _c in adata.obs.columns, f"obs 缺少列: {_c}"
print(f"  celltype 分布:\n{adata.obs['celltype'].value_counts().to_string()}")
print(f"  lesion 分布:\n{adata.obs['lesion'].value_counts().to_string()}")

# barcode 唯一性
assert adata.obs_names.is_unique, "barcode 不唯一"

# ======================== Step 2: 注释基因染色体坐标 ========================
print("\n" + "=" * 70)
print("Step 2: 注释基因坐标 (geneInfor)")
print("=" * 70)

geneInfor = pd.read_csv(GENEINFO_PATH)
print(f"  geneInfor: {geneInfor.shape[0]} genes")

df_var = adata.var.copy()
df_var["features"] = df_var.index
df_var = pd.merge(df_var, geneInfor, how="left", left_on="features", right_on="SYMBOL")
df_var = df_var.drop_duplicates(["features"])

adata.var["chromosome"] = df_var["chr"].values
adata.var["start"] = df_var["start"].values
adata.var["end"] = df_var["end"].values

n_before = adata.n_vars
adata = adata[:, adata.var.chromosome.notna()].copy()
print(f"  基因: {n_before} -> {adata.n_vars} (去除 {n_before - adata.n_vars} 个无坐标基因)")

chr_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
chr_map = {c: i for i, c in enumerate(chr_order)}
adata.var["chr_rank"] = adata.var["chromosome"].astype(str).map(chr_map)

n_before = adata.n_vars
adata = adata[:, adata.var["chr_rank"].notna()].copy()
print(f"  仅保留 chr1-22/X/Y: {n_before} -> {adata.n_vars} genes")

gene_order = adata.var.sort_values(["chr_rank", "start"]).index.tolist()
adata = adata[:, gene_order].copy()
print(f"  已按 chr + start 排序")

# ======================== Step 3: 保存上皮表达 UMAP ========================
print("\n" + "=" * 70)
print("Step 3: 上皮表达 UMAP (epi_UMAP_1/2)")
print("=" * 70)

if {"epi_UMAP_1", "epi_UMAP_2"}.issubset(adata.obs.columns):
    orig_umap = adata.obs[["epi_UMAP_1", "epi_UMAP_2"]].to_numpy(dtype=float)
    orig_barcodes = adata.obs_names.copy()
    print(f"  epi UMAP 可用 (NA 行数: {int(np.isnan(orig_umap).any(axis=1).sum())}, 为 T_ref)")
else:
    orig_umap, orig_barcodes = None, None
    print("  无 epi_UMAP 列, 后续跳过表达 UMAP 面板")

del _X, _Xd


# ======================== 工具函数 ========================
def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)


def _ax_umap(ax, title):
    _style(ax)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP 1", fontsize=12)
    ax.set_ylabel("UMAP 2", fontsize=12)
    ax.set_title(title, fontsize=14)


def plot_expr_umap(df_epi, w):
    """上皮表达 UMAP 5 面板 (09 号 Panel 风格):
    cluster_call / cnv_score / tissue / lesion / expr_cluster
    """
    plot_df = df_epi.dropna(subset=["epi_UMAP_1", "epi_UMAP_2"]).copy()
    if len(plot_df) > 100000:
        plot_df = plot_df.sample(100000, random_state=SEED)
    x = plot_df["epi_UMAP_1"].values
    y = plot_df["epi_UMAP_2"].values

    fig, axes = plt.subplots(2, 3, figsize=(24, 14), gridspec_kw=dict(wspace=0.35))
    cmap_tab20 = plt.get_cmap("tab20")

    # --- Panel 1: cluster_call ---
    ax = axes[0, 0]
    cc = {"Cancer": "#d62728", "Non-cancer": "#7f9cb3"}
    for call, color in cc.items():
        m = (plot_df["cluster_call"] == call).values
        if m.sum() == 0:
            continue
        ax.scatter(x[m], y[m], c=color, s=2, alpha=0.4, rasterized=True, label=call)
    _ax_umap(ax, "Cancer cluster call (M3)")
    ax.legend(fontsize=11, frameon=False, loc="best", markerscale=4)

    # --- Panel 2: cnv_score 连续色 ---
    ax = axes[0, 1]
    vals = plot_df["cnv_score"].values.astype(float)
    vmax = np.nanpercentile(vals, 98)
    sc_ = ax.scatter(x, y, c=vals, cmap="RdBu_r", s=2, alpha=0.5,
                     rasterized=True, vmin=0, vmax=vmax)
    plt.colorbar(sc_, ax=ax, fraction=0.046, pad=0.04).ax.tick_params(labelsize=10)
    _ax_umap(ax, f"CNV score (w{w})")

    # --- Panel 3: tissue ---
    ax = axes[0, 2]
    tc = {"EOPC": "#e41a1c", "LOPC": "#377eb8",
          "old_normal": "#4daf4a", "young_normal": "#984ea3"}
    for k, color in tc.items():
        m = (plot_df["tissue"] == k).values
        if m.sum() == 0:
            continue
        ax.scatter(x[m], y[m], c=color, s=2, alpha=0.4, rasterized=True, label=k)
    _ax_umap(ax, "Tissue")
    ax.legend(fontsize=11, frameon=False, loc="best", markerscale=4)

    # --- Panel 4: lesion ---
    ax = axes[1, 0]
    lc = {"normal": "#4daf4a", "low grade": "#ff7f00",
          "high grade": "#d62728", "BONE": "#377eb8"}
    for k, color in lc.items():
        m = (plot_df["lesion"] == k).values
        if m.sum() == 0:
            continue
        ax.scatter(x[m], y[m], c=color, s=2, alpha=0.4, rasterized=True, label=k)
    _ax_umap(ax, "Lesion")
    ax.legend(fontsize=11, frameon=False, loc="best", markerscale=4)

    # --- Panel 5: expr_cluster (tab20, 簇心标号) ---
    ax = axes[1, 1]
    ld = plot_df["expr_cluster"].astype(str)
    cats = sorted([v for v in ld.unique() if v != "NA"],
                  key=lambda v: int(v) if v.isdigit() else 999)
    for i, cat in enumerate(cats):
        m = (ld == cat).values
        ax.scatter(x[m], y[m], c=[cmap_tab20(i % 20)], s=2, alpha=0.4, rasterized=True)
        if m.sum() >= 20:
            ax.text(x[m].mean(), y[m].mean(), cat, fontsize=9, ha="center",
                    va="center", fontweight="bold", color="white",
                    path_effects=[pe.withStroke(linewidth=2.2, foreground="black")])
    _ax_umap(ax, "Expression cluster (seurat_clusters)")

    # --- Panel 6: cnv_leiden ---
    ax = axes[1, 2]
    cl = plot_df["cnv_leiden"].astype(str)
    cats = sorted(cl.unique(), key=lambda v: int(v) if v.isdigit() else 999)
    for i, cat in enumerate(cats):
        m = (cl == cat).values
        ax.scatter(x[m], y[m], c=[cmap_tab20(i % 20)], s=2, alpha=0.4, rasterized=True)
        if m.sum() >= 20:
            ax.text(x[m].mean(), y[m].mean(), cat, fontsize=9, ha="center",
                    va="center", fontweight="bold", color="white",
                    path_effects=[pe.withStroke(linewidth=2.2, foreground="black")])
    _ax_umap(ax, f"cnv_leiden (w{w})")

    plt.suptitle(
        f"Cancer cluster calling — PCa Tref (w{w}) | "
        f"normal<{NR_MAX:.0%} & CNA>=T_base+{SD_MULT:.0f}SD",
        fontsize=16, y=0.99)
    fig.subplots_adjust(top=0.93, bottom=0.05, hspace=0.25)
    path = os.path.join(OUT_DIR, f"01_expr_umap_Tref_w{w}_{today}.png")
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close("all")
    print(f"    -> {path}")


# ======================== Step 4: 逐窗口运行 ========================
print("\n" + "=" * 70)
print("Step 4: 逐窗口 inferCNV + 肿瘤簇判定")
print("=" * 70)

for w in WINDOW_SIZES:
    print(f"\n{'─' * 60}\n  Window size = {w}\n{'─' * 60}")
    t0 = datetime.now()

    adata_w = adata.copy()

    # ---- 4.1 inferCNV ----
    present_cats = adata_w.obs[REFERENCE_KEY].astype(str).unique().tolist()
    valid_refs = sorted(set(REFERENCE_CAT) & set(present_cats))
    assert len(valid_refs) == len(REFERENCE_CAT), \
        f"reference_cat 缺失: {set(REFERENCE_CAT) - set(valid_refs)}"
    print(f"  [4.1] cnv.tl.infercnv (reference_cat={valid_refs}) ...")
    cnv.tl.infercnv(
        adata_w,
        reference_key=REFERENCE_KEY,
        reference_cat=valid_refs,
        window_size=w,
        n_jobs=INFERCNV_N_JOBS,
        chunksize=INFERCNV_CHUNKSIZE,
    )

    # ---- 4.2 CNV PCA / Leiden / UMAP / Score ----
    print(f"  [4.2] CNV PCA + neighbors + Leiden + UMAP + Score ...")
    cnv.tl.pca(adata_w)
    cnv.pp.neighbors(adata_w)
    cnv.tl.leiden(adata_w)
    cnv.tl.umap(adata_w)
    cnv.tl.cnv_score(adata_w)
    print(f"    cnv_leiden clusters: {adata_w.obs['cnv_leiden'].nunique()}")

    # ---- 4.3 T 细胞基线 + 簇级判定 ----
    print(f"  [4.3] T 细胞基线 + cnv_leiden 簇级双条件判定 ...")
    obs = adata_w.obs.copy()
    if "X_cnv_umap" in adata_w.obsm:
        obs["cnv_UMAP_1"] = adata_w.obsm["X_cnv_umap"][:, 0]
        obs["cnv_UMAP_2"] = adata_w.obsm["X_cnv_umap"][:, 1]

    # 基线: T 细胞 cnv_score (ddof=0, 与 RCC 09 号 diploid 一致)
    t_scores = obs.loc[obs["celltype"] == "T_cells", "cnv_score"].dropna()
    t_mean, t_sd = t_scores.mean(), t_scores.std(ddof=0)
    thr = t_mean + SD_MULT * t_sd
    print(f"    T_cells baseline (n={len(t_scores):,}): "
          f"mean={t_mean:.6f} sd={t_sd:.6f} | thr = mean+{SD_MULT:.0f}sd = {thr:.6f}")

    epi_mask = obs["celltype"] == "Epithelium"
    rows = []
    for cl, g in obs[epi_mask].groupby("cnv_leiden"):
        n = len(g)
        n_normal = int((g["lesion"].astype(str) == NORMAL_LABEL).sum())
        normal_frac = n_normal / n if n else np.nan
        mean_cnv = g["cnv_score"].mean()
        sd_from_T = (mean_cnv - t_mean) / (t_sd + 1e-12)
        n_T = int((obs.loc[obs["cnv_leiden"] == cl, "celltype"] == "T_cells").sum())

        ld_vals = g["expr_cluster"].astype(str)
        ld_vals = ld_vals[ld_vals != "NA"]
        if len(ld_vals):
            ld_vc = ld_vals.value_counts()
            dom_ld = str(ld_vc.idxmax())
            top3 = ", ".join(f"C{k}:{v}" for k, v in ld_vc.head(3).items())
        else:
            dom_ld, top3 = "NA", ""
        dom_tissue = g["tissue"].value_counts().idxmax() if n else "NA"
        dom_lesion = g["lesion"].value_counts().idxmax() if n else "NA"

        cond_nr = normal_frac < NR_MAX
        cond_cnv = mean_cnv >= thr
        is_cancer = bool(cond_nr and cond_cnv)
        rows.append({
            "cnv_leiden": cl,
            "cluster_label": f"C{cl}",
            "n_epi": n, "n_T_ref": n_T,
            "n_normal": n_normal,
            "n_tumor_sample": n - n_normal,
            "normal_frac": normal_frac,
            "cnv_score_mean": mean_cnv,
            "cnv_score_median": g["cnv_score"].median(),
            "cnv_sd_from_Tbase": sd_from_T,
            "dominant_tissue": dom_tissue,
            "dominant_lesion": dom_lesion,
            "dominant_expr_cluster": dom_ld,
            "top3_expr_clusters": top3,
            "pass_normal_lt25pct": cond_nr,
            "pass_CNA_ge1sd": cond_cnv,
            "cancer_cluster": is_cancer,
            "call": "Cancer" if is_cancer else "Non-cancer",
        })

    summ = pd.DataFrame(rows)
    summ["cl_num"] = pd.to_numeric(summ["cnv_leiden"], errors="coerce").fillna(99).astype(int)
    summ = summ.sort_values("cl_num").reset_index(drop=True).drop(columns=["cl_num"])
    summ["T_base_mean"] = t_mean
    summ["T_base_sd"] = t_sd
    summ["threshold_cnv"] = thr

    # cluster_call: 上皮继承簇判定, T 细胞标记 T_ref
    call_map = dict(zip(summ["cnv_leiden"], summ["call"]))
    obs["cluster_call"] = obs["cnv_leiden"].map(call_map)
    obs.loc[obs["celltype"] == "T_cells", "cluster_call"] = "T_ref"

    n_cancer_cl = int(summ["cancer_cluster"].sum())
    n_cancer_cell = int((obs.loc[epi_mask, "cluster_call"] == "Cancer").sum())
    n_epi_total = int(epi_mask.sum())
    print(f"\n    -> {n_cancer_cl}/{len(summ)} cnv_leiden clusters called Cancer")
    print(f"    -> {n_cancer_cell:,}/{n_epi_total:,} epithelial cells in cancer "
          f"clusters ({n_cancer_cell / n_epi_total * 100:.1f}%)")
    print("\n    Cluster calls:")
    print(summ[["cluster_label", "n_epi", "n_T_ref", "normal_frac",
                "cnv_score_mean", "cnv_sd_from_Tbase",
                "dominant_expr_cluster", "call"]].to_string(index=False))

    # ---- 4.4 细胞级 CSV ----
    cell_cols = [c for c in [
        "orig.ident", "study", "tissue", "lesion", "celltype", "cell_source",
        "expr_cluster", "cnv_leiden", "cnv_score",
        "cnv_UMAP_1", "cnv_UMAP_2", "epi_UMAP_1", "epi_UMAP_2", "cluster_call",
    ] if c in obs.columns]
    df_cell = obs[cell_cols].copy()
    df_cell.insert(0, "barcode", df_cell.index)
    cell_path = os.path.join(OUT_DIR, f"01_cell_annotation_Tref_w{w}_{today}.csv")
    df_cell.to_csv(cell_path, index=False)
    print(f"\n  [4.4] -> {cell_path} ({len(df_cell):,} rows)")

    # ---- 4.5 簇级 CSV ----
    summ_path = os.path.join(OUT_DIR, f"01_cancer_cluster_call_Tref_w{w}_{today}.csv")
    summ.to_csv(summ_path, index=False)
    print(f"  [4.5] -> {summ_path} ({len(summ)} clusters)")

    # ---- 4.6 染色体 heatmap (每 cnv_leiden 簇子采样) ----
    print(f"  [4.6] chromosome heatmap (每簇 max {HEATMAP_MAX_PER_CLUST} cells) ...")
    rng = np.random.default_rng(SEED)
    idx_keep = []
    cl_all = adata_w.obs["cnv_leiden"].astype(str).values
    for cl in pd.unique(cl_all):
        idx_cl = np.where(cl_all == cl)[0]
        if len(idx_cl) > HEATMAP_MAX_PER_CLUST:
            idx_cl = np.sort(rng.choice(idx_cl, HEATMAP_MAX_PER_CLUST, replace=False))
        idx_keep.extend(idx_cl)
    idx_keep = np.array(sorted(idx_keep))

    # X=None 避免 X 列数与 var 行数不一致; heatmap 只用 obsm['X_cnv']
    adata_sub = ad.AnnData(
        X=None,
        obs=adata_w.obs.iloc[idx_keep][["cnv_leiden", "celltype"]].copy(),
        var=adata_w.var[["chromosome"]].copy(),
    )
    adata_sub.obsm["X_cnv"] = adata_w.obsm["X_cnv"][idx_keep]
    adata_sub.uns["cnv"] = adata_w.uns["cnv"]

    X_cnv_dense = (adata_sub.obsm["X_cnv"].toarray()
                   if hasattr(adata_sub.obsm["X_cnv"], "toarray")
                   else adata_sub.obsm["X_cnv"])
    vmin = float(np.percentile(X_cnv_dense, 1))
    vmax = float(np.percentile(X_cnv_dense, 99))
    del X_cnv_dense

    heatmap_path = os.path.join(OUT_DIR, f"01_heatmap_Tref_w{w}_{today}.png")
    try:
        plt.figure(figsize=(16, 30))
        cnv.pl.chromosome_heatmap(adata_sub, groupby="cnv_leiden",
                                  dendrogram=True, vmin=vmin, vmax=vmax)
        plt.savefig(heatmap_path, dpi=200, bbox_inches="tight")
    except Exception as e:
        print(f"    dendrogram heatmap 失败 ({e}), 降级为无树状图重画")
        plt.close("all")
        plt.figure(figsize=(16, 30))
        cnv.pl.chromosome_heatmap(adata_sub, groupby="cnv_leiden",
                                  dendrogram=False, vmin=vmin, vmax=vmax)
        plt.savefig(heatmap_path, dpi=200, bbox_inches="tight")
    plt.close("all")
    del adata_sub
    print(f"    -> {heatmap_path}")

    # ---- 4.7 CNV UMAP 4 面板 ----
    print(f"  [4.7] CNV UMAP plots ...")
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 12))
    cnv.pl.umap(adata_w, color="cnv_leiden", legend_loc="on data",
                legend_fontoutline=2, ax=ax1, show=False)
    cnv.pl.umap(adata_w, color="cnv_score", ax=ax2, show=False)
    cnv.pl.umap(adata_w, color="celltype", ax=ax3, show=False)
    # cluster_call 列注入后画第 4 面板
    adata_w.obs["cluster_call"] = obs["cluster_call"].reindex(adata_w.obs_names).values
    cnv.pl.umap(adata_w, color="cluster_call", ax=ax4, show=False)
    plt.suptitle(f"{PREFIX} - CNV UMAP (window={w})", fontsize=14)
    umap_path = os.path.join(OUT_DIR, f"01_cnv_umap_Tref_w{w}_{today}.png")
    plt.savefig(umap_path, dpi=300, bbox_inches="tight")
    plt.close("all")
    print(f"    -> {umap_path}")

    # ---- 4.8 上皮表达 UMAP 5+1 面板 ----
    print(f"  [4.8] epithelial expression UMAP panels ...")
    df_epi = obs[epi_mask].copy()
    if orig_umap is not None:
        umap_idx = orig_barcodes.get_indexer(df_epi.index)
        assert (umap_idx >= 0).all(), "epi UMAP barcode 对齐失败"
        df_epi["epi_UMAP_1"] = orig_umap[umap_idx, 0]
        df_epi["epi_UMAP_2"] = orig_umap[umap_idx, 1]
        plot_expr_umap(df_epi, w)
    else:
        print("    SKIP: 无 epi_UMAP 坐标")

    # ---- 4.9 X_cnv parquet + bin/chr 注释 ----
    print(f"  [4.9] X_cnv parquet + bin 注释 ...")
    X_cnv = adata_w.obsm["X_cnv"]
    if hasattr(X_cnv, "toarray"):
        X_cnv = X_cnv.toarray()
    n_bins = X_cnv.shape[1]
    df_cnv = pd.DataFrame(X_cnv.astype(np.float32),
                          index=adata_w.obs_names,
                          columns=[f"bin_{i:05d}" for i in range(n_bins)])
    df_cnv.index.name = "barcode"
    parquet_path = os.path.join(OUT_DIR, f"01_X_cnv_matrix_Tref_w{w}_{today}.parquet")
    df_cnv.reset_index().to_parquet(parquet_path, engine="pyarrow", compression="snappy")
    print(f"    -> {parquet_path} ({df_cnv.shape[0]:,} x {df_cnv.shape[1]:,})")
    del X_cnv, df_cnv

    chr_pos = (dict(adata_w.uns["cnv"]["chr_pos"])
               if "cnv" in adata_w.uns and "chr_pos" in adata_w.uns["cnv"] else {})
    chr_pos_sorted = sorted(chr_pos.items(), key=lambda kv: kv[1])
    bin_chr = [None] * n_bins
    for i, (c, start) in enumerate(chr_pos_sorted):
        end = chr_pos_sorted[i + 1][1] if i + 1 < len(chr_pos_sorted) else n_bins
        for j in range(start, end):
            bin_chr[j] = c
    bin_df = pd.DataFrame({"bin_idx": range(n_bins), "chromosome": bin_chr})
    bin_path = os.path.join(OUT_DIR, f"01_bin_annotation_Tref_w{w}_{today}.csv")
    bin_df.to_csv(bin_path, index=False)
    print(f"    -> {bin_path}")

    chr_pos_path = os.path.join(OUT_DIR, f"01_chr_pos_Tref_w{w}_{today}.json")
    with open(chr_pos_path, "w") as f_:
        _json.dump(chr_pos, f_, indent=2, default=int)
    print(f"    -> {chr_pos_path}")

    # cnv_leiden × celltype 交叉表
    crosstab = pd.crosstab(adata_w.obs["cnv_leiden"], adata_w.obs["celltype"])
    crosstab_path = os.path.join(OUT_DIR, f"01_cnv_leiden_x_celltype_Tref_w{w}_{today}.csv")
    crosstab.to_csv(crosstab_path)
    print(f"    -> {crosstab_path}")

    # cnv_score 统计
    score_summary = adata_w.obs.groupby("celltype")["cnv_score"].agg(
        ["mean", "median", "std", "count"])
    score_path = os.path.join(OUT_DIR, f"01_cnv_score_summary_Tref_w{w}_{today}.csv")
    score_summary.to_csv(score_path)
    print(f"    -> {score_path}")
    print(score_summary.to_string())

    # ---- 4.10 保存 h5ad ----
    if adata_w.raw is not None and "_index" in adata_w.raw.var.columns:
        print(f"  [4.10] 检测到 raw.var['_index'] 保留列, 清空 adata_w.raw 后再保存")
        adata_w.raw = None
    h5ad_out = os.path.join(OUT_DIR, f"01_{PREFIX}_w{w}_{today}.h5ad")
    adata_w.write(h5ad_out)
    print(f"  [4.10] -> {h5ad_out} "
          f"({os.path.getsize(h5ad_out) / 1024 ** 3:.2f} GB)")

    del adata_w
    print(f"\n  Window {w} done in "
          f"{(datetime.now() - t0).total_seconds() / 60:.1f} min")

print(f"\nAll done. 输出目录: {OUT_DIR}")
