suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(reticulate)
  library(harmony)
  library(viridis)
  library(tibble)
  library(FNN)
})

setwd("~/EOPC/malignant")
epi <- readRDS("epithelium.rds")

epi$celltype_epi <- case_when(
  epi$RNA_snn_res.0.5 %in% c(18,7) ~ "basal",
  epi$RNA_snn_res.0.5 %in% c(4) ~ "club",
  epi$RNA_snn_res.0.5 %in% c(3,21) ~ "KRT13+_basal_Hillock",
  TRUE ~ "luminal"  # 默认值
)
epi <- subset(epi, subset = celltype_epi %in% c("club","luminal"))
setwd("~/EOPC/malignant")
epi<-NormalizeData(epi,verbose = T) 
epi<-FindVariableFeatures(epi,selection.method = "vst", nfeatures = 2000)
epi<-ScaleData(epi,verbose = FALSE)
epi <- RunPCA(epi,verbose = T,npcs = 50)
epi <- RunHarmony(epi, group.by.vars = c("study","orig.ident"),plot_convergence = TRUE)
epi <- epi %>%
  RunUMAP(reduction = "harmony", dims = 1:15)%>%#
  FindNeighbors(reduction = "harmony", dims = 1:15) 


# =========================
# 0. 路径设置
# =========================

setwd("~/EOPC/malignant")

epi_file <- "epithelium.rds"
out_dir <- file.path(getwd(), "sctour_epi_output_reverse")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

reverse_sctour_direction <- TRUE

if (!file.exists(epi_file)) {
  stop("没有找到 epithelium.rds，请检查 setwd() 或 epi_file 路径。")
}

# =========================
# 1. 使用 scTour 环境
# =========================

env_name <- "sctour_env"
reticulate::use_condaenv(env_name, required = TRUE)

message("Python used:")
print(reticulate::py_config())

# =========================
# 2. 读取 epi
# =========================


DefaultAssay(epi) <- "RNA"

if (inherits(epi[["RNA"]], "Assay5")) {
  epi[["RNA"]] <- JoinLayers(epi[["RNA"]])
}

if (!"cell" %in% colnames(epi@meta.data)) {
  epi$cell <- rownames(epi@meta.data)
}

if (!"umap" %in% names(epi@reductions)) {
  stop("epi 里没有 umap reduction，请确认 reduction 名称。")
}

# =========================
# 3. 获取 counts
# =========================

counts <- GetAssayData(
  epi,
  assay = "RNA",
  layer = "counts"
)

if (!inherits(counts, "dgCMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

gene_keep <- Matrix::rowSums(counts) > 0
cell_keep <- Matrix::colSums(counts) > 0

counts <- counts[gene_keep, cell_keep, drop = FALSE]
epi <- subset(epi, cells = colnames(counts))

message("Cells after nonzero filtering: ", ncol(epi))
message("Genes after nonzero filtering: ", nrow(counts))

# =========================
# 4. 选择 scTour 使用基因
# =========================

if (length(VariableFeatures(epi)) < 500) {
  epi <- NormalizeData(
    epi,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  
  epi <- FindVariableFeatures(
    epi,
    assay = "RNA",
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
}

genes_use <- VariableFeatures(epi)
genes_use <- intersect(genes_use, rownames(counts))

if (length(genes_use) > 2000) {
  genes_use <- genes_use[1:2000]
}

if (length(genes_use) < 500) {
  gene_counts <- Matrix::rowSums(counts)
  genes_use <- names(sort(gene_counts, decreasing = TRUE))[1:min(2000, length(gene_counts))]
}

counts <- counts[genes_use, colnames(epi), drop = FALSE]

message("Genes used for scTour: ", nrow(counts))
message("Cells used for scTour: ", ncol(counts))

# =========================
# 5. 构建 AnnData，X 强制 float32
# =========================

counts_t <- Matrix::t(counts)

anndata <- reticulate::import("anndata", convert = FALSE)
scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)
np <- reticulate::import("numpy", convert = FALSE)

X_py <- scipy_sparse$csr_matrix(
  reticulate::r_to_py(counts_t)
)

X_py <- X_py$astype("float32")$copy()

adata <- anndata$AnnData(
  X = X_py
)

adata$obs_names <- reticulate::r_to_py(rownames(counts_t))
adata$var_names <- reticulate::r_to_py(colnames(counts_t))

# =========================
# 6. 写入 QC columns
# =========================

qc_total_counts <- as.numeric(Matrix::rowSums(counts_t))
qc_n_genes_by_counts <- as.numeric(Matrix::rowSums(counts_t > 0))

adata$obs$`__setitem__`("total_counts", reticulate::r_to_py(qc_total_counts))
adata$obs$`__setitem__`("n_genes_by_counts", reticulate::r_to_py(qc_n_genes_by_counts))
adata$obs$`__setitem__`("log1p_total_counts", reticulate::r_to_py(log1p(qc_total_counts)))
adata$obs$`__setitem__`("log1p_n_genes_by_counts", reticulate::r_to_py(log1p(qc_n_genes_by_counts)))

# =========================
# 7. 写入 metadata 和 UMAP
# =========================

meta <- epi@meta.data
meta <- meta[colnames(epi), , drop = FALSE]

for (cn in colnames(meta)) {
  value <- meta[[cn]]
  
  if (is.factor(value)) value <- as.character(value)
  if (is.logical(value)) value <- as.character(value)
  
  if (is.character(value)) {
    adata$obs$`__setitem__`(cn, reticulate::r_to_py(value))
  } else if (is.numeric(value) || is.integer(value)) {
    adata$obs$`__setitem__`(cn, reticulate::r_to_py(as.numeric(value)))
  } else {
    adata$obs$`__setitem__`(cn, reticulate::r_to_py(as.character(value)))
  }
}

umap_mat <- Embeddings(epi, reduction = "umap")
umap_mat <- umap_mat[colnames(epi), , drop = FALSE]
colnames(umap_mat)[1:2] <- c("UMAP_1", "UMAP_2")

adata$obsm$`__setitem__`(
  "X_umap_seurat",
  reticulate::r_to_py(as.matrix(umap_mat[, 1:2]))
)

message("adata.X dtype:")
print(reticulate::py_to_r(adata$X$dtype$name))

# =========================
# 8. 训练 scTour
# =========================

sct <- reticulate::import("sctour", convert = FALSE)

percent_train <- ifelse(ncol(epi) > 10000, 0.3, 0.8)

set.seed(123)

tnode <- sct$train$Trainer(
  adata,
  loss_mode = "nb",
  alpha_recon_lec = 0.5,
  alpha_recon_lode = 0.5,
  percent = percent_train,
  random_state = 123L
)

tnode$train()

# =========================
# 9. 提取 pseudotime 和 latent space
# =========================

ptime_original <- reticulate::py_to_r(tnode$get_time())
names(ptime_original) <- colnames(epi)

if (reverse_sctour_direction) {
  ptime_np <- np$array(
    reticulate::r_to_py(as.numeric(ptime_original)),
    dtype = "float32"
  )
  
  ptime_reversed <- reticulate::py_to_r(
    sct$train$reverse_time(ptime_np)
  )
  
  names(ptime_reversed) <- colnames(epi)
  
  epi$sctour_pseudotime_original <- ptime_original[colnames(epi)]
  epi$sctour_pseudotime <- ptime_reversed[colnames(epi)]
  
  message("scTour pseudotime has been reversed by sct.train.reverse_time().")
} else {
  epi$sctour_pseudotime_original <- ptime_original[colnames(epi)]
  epi$sctour_pseudotime <- ptime_original[colnames(epi)]
}

latent_out_py <- tnode$get_latentsp(
  alpha_z = 0.5,
  alpha_predz = 0.5
)

latent_out_r <- reticulate::py_to_r(latent_out_py)

mix_zs <- latent_out_r[[1]]
zs <- latent_out_r[[2]]

if (length(latent_out_r) >= 3) {
  pred_zs <- latent_out_r[[3]]
} else {
  pred_zs <- NULL
}

rownames(mix_zs) <- colnames(epi)
rownames(zs) <- colnames(epi)

colnames(mix_zs) <- paste0("scTour_", seq_len(ncol(mix_zs)))
colnames(zs) <- paste0("scTour_z_", seq_len(ncol(zs)))

if (!is.null(pred_zs)) {
  rownames(pred_zs) <- colnames(epi)
  colnames(pred_zs) <- paste0("scTour_predz_", seq_len(ncol(pred_zs)))
}

epi[["sctour"]] <- CreateDimReducObject(
  embeddings = mix_zs[colnames(epi), , drop = FALSE],
  key = "sctour_",
  assay = DefaultAssay(epi)
)

# =========================
# 10. 计算 vector field
# 重要：如果 pseudotime 反向，则画图使用 -vf_mat
# =========================

ptime_for_vf_np <- np$array(
  reticulate::r_to_py(as.numeric(ptime_original)),
  dtype = "float32"
)

mix_zs_np <- np$array(
  reticulate::r_to_py(as.matrix(mix_zs)),
  dtype = "float32"
)

vf <- tnode$get_vector_field(
  ptime_for_vf_np,
  mix_zs_np
)

vf_mat_original <- reticulate::py_to_r(vf)
rownames(vf_mat_original) <- rownames(mix_zs)

if (reverse_sctour_direction) {
  vf_mat <- -vf_mat_original
  message("Vector field has been reversed by multiplying -1, equivalent to reverse=True in scTour plotting.")
} else {
  vf_mat <- vf_mat_original
}

# =========================
# 11. 保存 AnnData 和 Seurat object
# =========================

adata$obs$`__setitem__`(
  "ptime_original",
  reticulate::r_to_py(epi$sctour_pseudotime_original)
)

adata$obs$`__setitem__`(
  "ptime",
  reticulate::r_to_py(epi$sctour_pseudotime)
)

adata$obsm$`__setitem__`(
  "X_TNODE",
  reticulate::r_to_py(mix_zs)
)

adata$obsm$`__setitem__`(
  "X_VF_original",
  reticulate::r_to_py(vf_mat_original)
)

adata$obsm$`__setitem__`(
  "X_VF",
  reticulate::r_to_py(vf_mat)
)

if (!is.null(pred_zs)) {
  adata$obsm$`__setitem__`(
    "X_pred_zs",
    reticulate::r_to_py(pred_zs)
  )
}

adata_file <- file.path(out_dir, "epi_sctour_result_reverse_adjusted.h5ad")
adata$write_h5ad(adata_file)

saveRDS(
  epi,
  file = file.path(out_dir, "epithelium_with_sctour_reverse_adjusted.rds")
)

# =========================
# 12. Seurat UMAP pseudotime
# =========================

p_ptime <- FeaturePlot(
  epi,
  features = "sctour_pseudotime",
  reduction = "umap",
  pt.size = 0.1
) +
  scale_color_viridis_c(option = "plasma") +
  ggtitle("scTour pseudotime on Seurat UMAP") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

pdf(
  file.path(out_dir, "epi_sctour_pseudotime_umap_reverse_adjusted.pdf"),
  width = 5,
  height = 4.5
)
print(p_ptime)
dev.off()

# =========================
# 13. 原始 vs 反向 pseudotime 对比
# =========================

p_ptime_original <- FeaturePlot(
  epi,
  features = "sctour_pseudotime_original",
  reduction = "umap",
  pt.size = 0.1
) +
  scale_color_viridis_c(option = "plasma") +
  ggtitle("Original scTour pseudotime") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_ptime_reversed <- FeaturePlot(
  epi,
  features = "sctour_pseudotime",
  reduction = "umap",
  pt.size = 0.1
) +
  scale_color_viridis_c(option = "plasma") +
  ggtitle("Reverse-adjusted scTour pseudotime") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

pdf(
  file.path(out_dir, "epi_sctour_pseudotime_original_vs_reverse.pdf"),
  width = 10,
  height = 4.5
)
print(p_ptime_original + p_ptime_reversed)
dev.off()

# =========================
# 14. lesion 横向 pseudotime
# =========================

if ("lesion" %in% colnames(epi@meta.data)) {
  epi$lesion <- factor(
    epi$lesion,
    levels = c("normal", "low grade", "high grade", "BONE")
  )
  
  df_plot <- Embeddings(epi, "umap") %>%
    as.data.frame()
  
  colnames(df_plot)[1:2] <- c("umap_1", "umap_2")
  df_plot$cell_id <- rownames(df_plot)
  
  meta_plot <- epi@meta.data %>%
    as.data.frame()
  
  meta_plot$cell_id <- rownames(meta_plot)
  
  meta_plot <- meta_plot %>%
    dplyr::select(
      cell_id,
      lesion,
      sctour_pseudotime
    )
  
  df_plot <- df_plot %>%
    dplyr::left_join(meta_plot, by = "cell_id") %>%
    dplyr::filter(
      !is.na(lesion),
      !is.na(sctour_pseudotime)
    )
  
  p_ptime_lesion <- ggplot(
    df_plot,
    aes(
      x = umap_1,
      y = umap_2,
      color = sctour_pseudotime
    )
  ) +
    geom_point(size = 0.1, alpha = 0.9) +
    facet_wrap(~ lesion, nrow = 1) +
    scale_color_viridis_c(option = "plasma") +
    labs(
      x = "umap_1",
      y = "umap_2",
      color = "scTour\npseudotime"
    ) +
    theme_classic(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      axis.text = element_text(color = "black", size = 7),
      axis.title = element_text(color = "black", size = 9),
      axis.line = element_line(color = "black", linewidth = 0.4),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
  
  pdf(
    file.path(out_dir, "epi_sctour_pseudotime_by_lesion_horizontal_reverse_adjusted.pdf"),
    width = 12,
    height = 3.6
  )
  print(p_ptime_lesion)
  dev.off()
}

# =========================
# 15. scTour latent space vector field
# =========================

latent_vector_df <- data.frame(
  cell_id = rownames(mix_zs),
  x = mix_zs[, 1],
  y = mix_zs[, 2],
  vx = vf_mat[, 1],
  vy = vf_mat[, 2],
  pseudotime = epi$sctour_pseudotime[match(rownames(mix_zs), colnames(epi))]
)

latent_vector_df <- latent_vector_df %>%
  dplyr::filter(
    is.finite(x),
    is.finite(y),
    is.finite(vx),
    is.finite(vy),
    is.finite(pseudotime)
  )

latent_arrow_scale <- 0.08

latent_vector_df <- latent_vector_df %>%
  dplyr::mutate(
    xend = x + vx * latent_arrow_scale,
    yend = y + vy * latent_arrow_scale
  )

set.seed(123)

latent_arrow_df <- latent_vector_df %>%
  dplyr::mutate(random_order = runif(dplyr::n())) %>%
  dplyr::arrange(random_order) %>%
  dplyr::slice_head(n = 1200) %>%
  dplyr::select(-random_order)

p_vector_latent <- ggplot(latent_vector_df, aes(x = x, y = y)) +
  geom_point(
    aes(color = pseudotime),
    size = 0.22,
    alpha = 0.75
  ) +
  geom_segment(
    data = latent_arrow_df,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    inherit.aes = FALSE,
    arrow = arrow(
      length = unit(0.025, "inches"),
      type = "closed"
    ),
    linewidth = 0.13,
    color = "black",
    alpha = 0.45
  ) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "scTour predicted vector field",
    subtitle = ifelse(reverse_sctour_direction, "reverse-adjusted vector field", "original vector field"),
    x = "scTour latent 1",
    y = "scTour latent 2",
    color = "Pseudotime"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

pdf(
  file.path(out_dir, "epi_sctour_predicted_vector_field_latent_space_reverse_adjusted.pdf"),
  width = 5,
  height = 4.5
)
print(p_vector_latent)
dev.off()

# =========================
# 16. Seurat UMAP 上 grid-smoothed vector field
# =========================

common_cells <- intersect(
  rownames(mix_zs),
  colnames(epi)
)

umap_mat <- Embeddings(epi, reduction = "umap")
umap_mat <- umap_mat[common_cells, , drop = FALSE]
colnames(umap_mat)[1:2] <- c("umap_1", "umap_2")

latent_mat <- as.matrix(mix_zs)
latent_mat <- latent_mat[common_cells, , drop = FALSE]

vf_mat2 <- as.matrix(vf_mat)
vf_mat2 <- vf_mat2[common_cells, , drop = FALSE]

pseudotime_vec <- epi$sctour_pseudotime[match(common_cells, colnames(epi))]

latent_step <- 0.05

latent_endpoint <- latent_mat + vf_mat2 * latent_step

nn <- FNN::get.knnx(
  data = latent_mat,
  query = latent_endpoint,
  k = 8
)

start_umap <- umap_mat

end_umap <- matrix(
  NA_real_,
  nrow = nrow(umap_mat),
  ncol = 2
)

for (i in seq_len(nrow(umap_mat))) {
  neighbor_idx <- nn$nn.index[i, ]
  
  candidate_end <- umap_mat[neighbor_idx, , drop = FALSE]
  
  d_umap <- sqrt(
    (candidate_end[, 1] - start_umap[i, 1])^2 +
      (candidate_end[, 2] - start_umap[i, 2])^2
  )
  
  local_keep <- d_umap <= quantile(d_umap, 0.75, na.rm = TRUE)
  
  if (sum(local_keep) >= 2) {
    end_umap[i, ] <- colMeans(candidate_end[local_keep, , drop = FALSE])
  } else {
    end_umap[i, ] <- candidate_end[which.min(d_umap), ]
  }
}

vector_cell_df <- data.frame(
  cell_id = common_cells,
  umap_1 = start_umap[, 1],
  umap_2 = start_umap[, 2],
  umap_1_end = end_umap[, 1],
  umap_2_end = end_umap[, 2],
  pseudotime = pseudotime_vec
)

vector_cell_df <- vector_cell_df %>%
  dplyr::mutate(
    dx = umap_1_end - umap_1,
    dy = umap_2_end - umap_2,
    arrow_len = sqrt(dx^2 + dy^2)
  ) %>%
  dplyr::filter(
    is.finite(umap_1),
    is.finite(umap_2),
    is.finite(dx),
    is.finite(dy),
    is.finite(pseudotime),
    arrow_len > 0
  ) %>%
  dplyr::filter(
    arrow_len <= quantile(arrow_len, 0.65, na.rm = TRUE)
  )

grid_n <- 35

vector_cell_df <- vector_cell_df %>%
  dplyr::mutate(
    grid_x = cut(umap_1, breaks = grid_n, labels = FALSE),
    grid_y = cut(umap_2, breaks = grid_n, labels = FALSE)
  )

grid_df <- vector_cell_df %>%
  dplyr::group_by(grid_x, grid_y) %>%
  dplyr::summarise(
    x = median(umap_1, na.rm = TRUE),
    y = median(umap_2, na.rm = TRUE),
    dx = median(dx, na.rm = TRUE),
    dy = median(dy, na.rm = TRUE),
    n_cells = dplyr::n(),
    pt = median(pseudotime, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    n_cells >= 15,
    is.finite(dx),
    is.finite(dy)
  ) %>%
  dplyr::mutate(
    len = sqrt(dx^2 + dy^2)
  ) %>%
  dplyr::filter(len > 0)

arrow_len_fixed <- 0.35

grid_df <- grid_df %>%
  dplyr::mutate(
    dx_unit = dx / len,
    dy_unit = dy / len,
    xend = x + dx_unit * arrow_len_fixed,
    yend = y + dy_unit * arrow_len_fixed
  )

plot_point_df <- data.frame(
  cell_id = common_cells,
  umap_1 = umap_mat[, 1],
  umap_2 = umap_mat[, 2],
  pseudotime = pseudotime_vec
)

p_vector_umap_grid <- ggplot(
  plot_point_df,
  aes(x = umap_1, y = umap_2)
) +
  geom_point(
    aes(color = pseudotime),
    size = 0.12,
    alpha = 0.8
  ) +
  geom_segment(
    data = grid_df,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    inherit.aes = FALSE,
    arrow = arrow(
      length = unit(0.035, "inches"),
      type = "closed"
    ),
    linewidth = 0.28,
    color = "black",
    alpha = 0.75
  ) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "scTour vector field on Seurat UMAP",
    subtitle = "Grid-smoothed local direction, reverse-adjusted",
    x = "umap_1",
    y = "umap_2",
    color = "Pseudotime"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

pdf(
  file.path(out_dir, "epi_sctour_vector_field_on_umap_grid_smoothed_reverse_adjusted.pdf"),
  width = 5,
  height = 4.5
)
print(p_vector_umap_grid)
dev.off()

# =========================
# 17. 保存最终对象
# =========================

saveRDS(
  epi,
  file = file.path(out_dir, "epithelium_with_sctour_reverse_adjusted_final.rds")
)

message("Finished.")
message("Output directory: ", out_dir)
message("Saved Seurat object: ", file.path(out_dir, "epithelium_with_sctour_reverse_adjusted_final.rds"))
message("Saved AnnData: ", adata_file)