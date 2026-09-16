# =========================
# 0. 安装环境
# 第一次运行需要安装
# =========================

# 如果你已经有 squidpy 环境，可以跳过这个 cell
# 在 Jupyter 中运行：
# !pip install squidpy scanpy anndata pandas numpy scipy matplotlib seaborn statsmodels openpyxl

# =========================
# 1. 导入包
# =========================

import os
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")

# =========================
# 2. 路径和样本信息
# =========================

base_dir = Path("/mnt/DATA/home/zqy1234560915/EOPC/RCTD")
input_dir = base_dir / "CytoSPACE_input"
out_dir = base_dir / "Squidpy_CytoSPACE_7samples_nhood_enrichment"
out_dir.mkdir(parents=True, exist_ok=True)

sample_info = pd.DataFrame({
    "sample": [
        "NEADT_11",
        "NEADT_3",
        "NEADT_5",
        "NEADT_7",
        "NEADT_8",
        "patient 11 treatment-naive",
        "patient 3 treatment-naive",
    ],
    "group": [
        "LOPC",
        "EOPC",
        "EOPC",
        "LOPC",
        "EOPC",
        "EOPC",
        "LOPC",
    ]
})

print(sample_info)

# =========================
# 3. 参数设置
# =========================

cluster_key = "CellType_simplified"

# permutation 次数
# 先测试可以设 500，正式结果建议 5000 或 10000
n_perms = 5000

# 如果细胞太多，可抽样；不想抽样设为 None
max_cells_per_sample = None

# 目标互作
target_a_alias = [
    "SPP1_TAM",
    "SPP1_TAMs",
    "SPP1_plus_TAM",
    "SPP1_macro",
    "SPP1_macrophage"
]

target_b_alias = [
    "malignant_epi_0",
    "malignant_epi0",
    "malignant_epithelium_0",
    "malignant_epi_0",
    "malignant_epithelium",
    "malignant_epi",
    "malignant"
]

# =========================
# 4. 细胞类型名称清洗
# =========================

def safe_name(x):
    x = str(x)
    x = re.sub(r"[/\+\-\s]+", "_", x)
    x = re.sub(r"[^A-Za-z0-9_]", "_", x)
    x = re.sub(r"_+", "_", x)
    x = re.sub(r"^_|_$", "", x)
    return x

def simplify_celltype(x):
    x0 = safe_name(x)

    # 这里根据你之前的 CytoSPACE label 做兼容
    low = x0.lower()

    if low in [s.lower() for s in target_a_alias]:
        return "SPP1_TAM"

    if low in [s.lower() for s in target_b_alias]:
        return "malignant_epi_0"

    return x0

# =========================
# 5. 读取一个样本的 CytoSPACE assigned_locations
# =========================

def read_one_assigned(sample_name, group_name):
    f = input_dir / sample_name / "cytospace_results" / "assigned_locations.csv"

    if not f.exists():
        raise FileNotFoundError(f"Cannot find: {f}")

    df = pd.read_csv(f)

    print("\n", sample_name)
    print(df.columns.tolist())

    required = {"CellType", "SpotID"}
    if not required.issubset(df.columns):
        raise ValueError(
            f"{sample_name} assigned_locations.csv must contain CellType and SpotID. "
            f"Current columns: {df.columns.tolist()}"
        )

    # CytoSPACE 输出通常包含 row / col，可直接作为空间坐标
    if {"row", "col"}.issubset(df.columns):
        x_col, y_col = "col", "row"
    elif {"x", "y"}.issubset(df.columns):
        x_col, y_col = "x", "y"
    else:
        raise ValueError(
            f"{sample_name} has no row/col or x/y coordinates in assigned_locations.csv."
        )

    df = df.copy()
    df["sample"] = sample_name
    df["group"] = group_name
    df["CellType_raw"] = df["CellType"].astype(str)
    df[cluster_key] = df["CellType_raw"].map(simplify_celltype)

    df["x"] = pd.to_numeric(df[x_col], errors="coerce")
    df["y"] = pd.to_numeric(df[y_col], errors="coerce")

    df = df.dropna(subset=["x", "y", cluster_key])
    df = df[df[cluster_key] != ""]

    # 每个 mapped cell 作为一个 observation
    if "UniqueCID" in df.columns:
        obs_names = df["UniqueCID"].astype(str).values
    elif "OriginalCID" in df.columns:
        obs_names = (
            sample_name + "_" + df["OriginalCID"].astype(str)
        ).values
    else:
        obs_names = (
            sample_name + "_cell_" + np.arange(df.shape[0]).astype(str)
        )

    obs_names = pd.Index(obs_names).astype(str)
    obs_names = obs_names.where(~obs_names.duplicated(), obs_names + "_" + np.arange(len(obs_names)).astype(str))

    obs = df[[
        "sample",
        "group",
        "SpotID",
        "CellType_raw",
        cluster_key,
        "x",
        "y"
    ]].copy()

    obs.index = obs_names

    # Squidpy 做 nhood enrichment 只需要 obs + spatial coordinates
    X = np.zeros((obs.shape[0], 1), dtype=np.float32)

    adata = ad.AnnData(X=X, obs=obs)
    adata.obsm["spatial"] = obs[["x", "y"]].to_numpy(dtype=np.float32)

    adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

    if max_cells_per_sample is not None and adata.n_obs > max_cells_per_sample:
        np.random.seed(123)
        idx = np.random.choice(adata.n_obs, max_cells_per_sample, replace=False)
        adata = adata[idx, :].copy()

    print("n mapped cells:", adata.n_obs)
    print(adata.obs[cluster_key].value_counts())

    return adata

# =========================
# 6. 自动估计 immediate neighbor 半径
# =========================

def estimate_radius(adata):
    from sklearn.neighbors import NearestNeighbors

    xy = adata.obsm["spatial"]
    nn = NearestNeighbors(n_neighbors=2)
    nn.fit(xy)
    dists, _ = nn.kneighbors(xy)

    # 第 2 列是最近邻距离，第 1 列是自己
    radius = np.nanmedian(dists[:, 1]) * 1.5

    if not np.isfinite(radius) or radius <= 0:
        radius = None

    return radius

# =========================
# 7. 对单个样本跑 Squidpy neighborhood enrichment
# =========================

def run_squidpy_one_sample(adata, sample_name):
    print(f"\nRunning Squidpy nhood enrichment: {sample_name}")

    radius = estimate_radius(adata)
    print("estimated radius:", radius)

    if radius is not None:
        sq.gr.spatial_neighbors(
            adata,
            spatial_key="spatial",
            coord_type="generic",
            radius=radius
        )
    else:
        sq.gr.spatial_neighbors(
            adata,
            spatial_key="spatial",
            coord_type="generic",
            n_neighs=6
        )

    sq.gr.nhood_enrichment(
        adata,
        cluster_key=cluster_key,
        n_perms=n_perms,
        seed=123,
        n_jobs=8
    )

    return adata

# =========================
# 8. 提取 enrichment matrix
# =========================

def extract_nhood_matrix(adata):
    key = f"{cluster_key}_nhood_enrichment"

    cats = list(adata.obs[cluster_key].cat.categories)

    z = adata.uns[key]["zscore"]
    c = adata.uns[key]["count"]

    z_df = pd.DataFrame(z, index=cats, columns=cats)
    c_df = pd.DataFrame(c, index=cats, columns=cats)

    return z_df, c_df

def get_pair_value(mat, a="SPP1_TAM", b="malignant_epi_0"):
    if a not in mat.index or b not in mat.columns:
        return np.nan

    v1 = mat.loc[a, b]
    v2 = mat.loc[b, a] if (b in mat.index and a in mat.columns) else np.nan

    return np.nanmean([v1, v2])

# =========================
# 9. 批量运行所有样本
# =========================

adata_dict = {}
z_dict = {}
count_dict = {}
pair_records = []

for _, row in sample_info.iterrows():
    sample_name = row["sample"]
    group_name = row["group"]

    adata = read_one_assigned(sample_name, group_name)
    adata = run_squidpy_one_sample(adata, sample_name)

    z_df, c_df = extract_nhood_matrix(adata)

    adata_dict[sample_name] = adata
    z_dict[sample_name] = z_df
    count_dict[sample_name] = c_df

    pair_records.append({
        "sample": sample_name,
        "group": group_name,
        "pair": "SPP1_TAM__malignant_epi_0",
        "zscore_mean_bidirectional": get_pair_value(z_df),
        "count_mean_bidirectional": get_pair_value(c_df),
        "n_cells": adata.n_obs,
        "n_SPP1_TAM": int((adata.obs[cluster_key] == "SPP1_TAM").sum()) if "SPP1_TAM" in adata.obs[cluster_key].cat.categories else 0,
        "n_malignant_epi_0": int((adata.obs[cluster_key] == "malignant_epi_0").sum()) if "malignant_epi_0" in adata.obs[cluster_key].cat.categories else 0,
    })

pair_df = pd.DataFrame(pair_records)
print(pair_df)

# =========================
# 10. 汇总 EOPC / LOPC enrichment matrix
# =========================

def mean_matrix_by_group(z_dict, sample_info, group_name):
    samples = sample_info.loc[sample_info["group"] == group_name, "sample"].tolist()
    mats = [z_dict[s] for s in samples if s in z_dict]

    all_rows = sorted(set().union(*[set(m.index) for m in mats]))
    all_cols = sorted(set().union(*[set(m.columns) for m in mats]))

    aligned = []
    for m in mats:
        tmp = pd.DataFrame(np.nan, index=all_rows, columns=all_cols)
        tmp.loc[m.index, m.columns] = m
        aligned.append(tmp)

    mean_mat = sum(aligned) / len(aligned)
    return mean_mat

z_eopc = mean_matrix_by_group(z_dict, sample_info, "EOPC")
z_lopc = mean_matrix_by_group(z_dict, sample_info, "LOPC")

all_types = sorted(set(z_eopc.index).union(set(z_lopc.index)))

z_eopc = z_eopc.reindex(index=all_types, columns=all_types)
z_lopc = z_lopc.reindex(index=all_types, columns=all_types)

z_delta = z_eopc - z_lopc

# =========================
# 11. 统计 SPP1_TAM - malignant_epi_0 的 EOPC vs LOPC
# =========================

eopc_values = pair_df.loc[pair_df["group"] == "EOPC", "zscore_mean_bidirectional"].dropna()
lopc_values = pair_df.loc[pair_df["group"] == "LOPC", "zscore_mean_bidirectional"].dropna()

if len(eopc_values) >= 1 and len(lopc_values) >= 1:
    p_pair = mannwhitneyu(eopc_values, lopc_values, alternative="two-sided").pvalue
else:
    p_pair = np.nan

print("SPP1_TAM - malignant_epi_0 enrichment z-score")
print("EOPC:", eopc_values.tolist())
print("LOPC:", lopc_values.tolist())
print("Mann-Whitney p:", p_pair)

# =========================
# 12. 作图函数
# =========================

def plot_heatmap(mat, title, cmap="Reds", center=None, vmin=None, vmax=None):
    plt.figure(figsize=(9, 7))
    sns.heatmap(
        mat,
        cmap=cmap,
        center=center,
        vmin=vmin,
        vmax=vmax,
        square=True,
        linewidths=0.2,
        linecolor="white",
        cbar_kws={"label": "Neighborhood enrichment z-score"}
    )
    plt.title(title, fontsize=14, fontweight="bold")
    plt.xlabel("Neighbor / Predictor cell type")
    plt.ylabel("Target cell type")
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()

def save_heatmap(mat, title, filename, cmap="Reds", center=None, vmin=None, vmax=None):
    plt.figure(figsize=(9, 7))
    sns.heatmap(
        mat,
        cmap=cmap,
        center=center,
        vmin=vmin,
        vmax=vmax,
        square=True,
        linewidths=0.2,
        linecolor="white",
        cbar_kws={"label": "Neighborhood enrichment z-score"}
    )
    plt.title(title, fontsize=14, fontweight="bold")
    plt.xlabel("Neighbor / Predictor cell type")
    plt.ylabel("Target cell type")
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(out_dir / filename, format="pdf")
    plt.close()

# =========================
# 13. 保存 EOPC / LOPC / Delta 热图
# =========================

max_abs_delta = np.nanmax(np.abs(z_delta.to_numpy()))

save_heatmap(
    z_eopc,
    "Squidpy neighborhood enrichment: EOPC",
    "Squidpy_nhood_enrichment_EOPC_mean_zscore.pdf",
    cmap="Reds"
)

save_heatmap(
    z_lopc,
    "Squidpy neighborhood enrichment: LOPC",
    "Squidpy_nhood_enrichment_LOPC_mean_zscore.pdf",
    cmap="Reds"
)

save_heatmap(
    z_delta,
    "Squidpy neighborhood enrichment: EOPC minus LOPC",
    "Squidpy_nhood_enrichment_EOPC_minus_LOPC_delta_zscore.pdf",
    cmap="vlag",
    center=0,
    vmin=-max_abs_delta,
    vmax=max_abs_delta
)

# =========================
# 14. SPP1_TAM - malignant_epi_0 pair boxplot
# =========================

plt.figure(figsize=(4.2, 4.2))

sns.boxplot(
    data=pair_df,
    x="group",
    y="zscore_mean_bidirectional",
    order=["LOPC", "EOPC"],
    palette={"LOPC": "#4DBBD5", "EOPC": "#E64B35"},
    width=0.45,
    showfliers=False
)

sns.stripplot(
    data=pair_df,
    x="group",
    y="zscore_mean_bidirectional",
    order=["LOPC", "EOPC"],
    color="black",
    size=5,
    jitter=0.08
)

plt.title("SPP1+ TAM - malignant epi-0\nneighborhood enrichment", fontweight="bold")
plt.xlabel("")
plt.ylabel("Mean bidirectional enrichment z-score")
plt.text(
    0.5,
    pair_df["zscore_mean_bidirectional"].max() * 1.05,
    f"p = {p_pair:.3g}" if np.isfinite(p_pair) else "p = NA",
    ha="center"
)
plt.tight_layout()
plt.savefig(out_dir / "Squidpy_SPP1_TAM_malignant_epi0_EOPC_vs_LOPC_boxplot.pdf")
plt.close()

# =========================
# 15. 单样本热图也保存
# =========================

sample_heatmap_dir = out_dir / "per_sample_heatmaps"
sample_heatmap_dir.mkdir(exist_ok=True)

for sample_name, z_df in z_dict.items():
    save_heatmap(
        z_df,
        f"Squidpy neighborhood enrichment: {sample_name}",
        f"per_sample_heatmaps/{sample_name}_nhood_enrichment_zscore.pdf",
        cmap="Reds"
    )

# =========================
# 16. 保存结果表
# =========================

pair_df.to_csv(out_dir / "Squidpy_SPP1_TAM_malignant_epi0_pair_zscore_by_sample.csv", index=False)
z_eopc.to_csv(out_dir / "Squidpy_EOPC_mean_nhood_enrichment_zscore_matrix.csv")
z_lopc.to_csv(out_dir / "Squidpy_LOPC_mean_nhood_enrichment_zscore_matrix.csv")
z_delta.to_csv(out_dir / "Squidpy_EOPC_minus_LOPC_delta_nhood_enrichment_zscore_matrix.csv")

print("Finished.")
print("Output directory:", out_dir)
print("Key result:")
print(pair_df)