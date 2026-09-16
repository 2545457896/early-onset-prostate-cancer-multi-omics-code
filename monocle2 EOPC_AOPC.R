suppressMessages({
  library(Seurat)
  library(monocle)
  library(Biobase)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

set.seed(123)

setwd("/home/rstudio/single cell EOPC/epi")

epi <- readRDS("epithelium.rds")

epi$celltype_epi <- dplyr::case_when(
  epi$RNA_snn_res.0.5 %in% c(18, 7) ~ "basal",
  epi$RNA_snn_res.0.5 %in% c(4) ~ "club",
  epi$RNA_snn_res.0.5 %in% c(3, 21) ~ "KRT13+_basal_Hillock",
  TRUE ~ "luminal"
)

epi$celltype_epi_new <- paste0("c0", epi$RNA_snn_res.0.5, "_", epi$celltype_epi)

epi_use <- subset(epi, subset = celltype_epi != "basal" & celltype_epi != "KRT13+_basal_Hillock")

table(epi_use$celltype_epi)
table(epi_use$tissue)

outdir <- "/home/rstudio/single cell EOPC/epi/monocle2"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

get_root_state <- function(cds, root_celltype = "club") {
  state_tab <- table(pData(cds)$State, pData(cds)$celltype_epi)
  if (!root_celltype %in% colnames(state_tab)) {
    stop("Root cell type not found: ", root_celltype)
  }
  root_counts <- state_tab[, root_celltype]
  root_state <- as.numeric(names(root_counts)[which.max(root_counts)])
  return(root_state)
}

run_monocle2 <- function(seu, tissue_keep, prefix, root_celltype = "club", nfeatures = 2000) {
  message("========================================")
  message("Running: ", prefix)
  message("Tissues: ", paste(tissue_keep, collapse = ", "))
  message("========================================")
  
  seu_sub <- subset(seu, subset = tissue %in% tissue_keep)
  DefaultAssay(seu_sub) <- "RNA"
  
  message("Cell number: ", ncol(seu_sub))
  print(table(seu_sub$tissue))
  print(table(seu_sub$celltype_epi))
  
  if (!root_celltype %in% unique(seu_sub$celltype_epi)) {
    stop("No ", root_celltype, " cells found in ", prefix)
  }
  
  seu_sub <- NormalizeData(seu_sub, verbose = FALSE)
  seu_sub <- FindVariableFeatures(
    seu_sub,
    selection.method = "vst",
    nfeatures = nfeatures,
    verbose = FALSE
  )
  
  ordering_genes <- VariableFeatures(seu_sub)
  ordering_genes <- head(ordering_genes, nfeatures)
  
  counts <- GetAssayData(seu_sub, assay = "RNA", layer = "counts")
  counts <- as(counts, "dgCMatrix")
  
  pd <- seu_sub@meta.data
  pd$cell_id <- rownames(pd)
  pd <- pd[colnames(counts), , drop = FALSE]
  
  fd <- data.frame(
    gene_short_name = rownames(counts),
    row.names = rownames(counts),
    stringsAsFactors = FALSE
  )
  
  pd_annotated <- new("AnnotatedDataFrame", data = pd)
  fd_annotated <- new("AnnotatedDataFrame", data = fd)
  
  cds <- newCellDataSet(
    counts,
    phenoData = pd_annotated,
    featureData = fd_annotated,
    expressionFamily = negbinomial.size(),
    lowerDetectionLimit = 0.5
  )
  
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  ordering_genes <- intersect(ordering_genes, rownames(cds))
  cds <- setOrderingFilter(cds, ordering_genes)
  
  message("Ordering genes used: ", length(ordering_genes))
  
  pdf(paste0(prefix, "_ordering_genes.pdf"), width = 6, height = 5)
  print(plot_ordering_genes(cds))
  dev.off()
  
  cds <- reduceDimension(
    cds,
    max_components = 2,
    method = "DDRTree",
    norm_method = "log",
    pseudo_expr = 1,
    relative_expr = TRUE,
    auto_param_selection = TRUE,
    verbose = TRUE
  )
  
  cds <- orderCells(cds)
  
  root_state <- get_root_state(
    cds,
    root_celltype = root_celltype
  )
  
  message("Root state for ", prefix, ": ", root_state)
  
  cds <- orderCells(
    cds,
    root_state = root_state
  )
  
  saveRDS(
    cds,
    file = paste0(prefix, "_monocle2.rds")
  )
  
  root_info <- data.frame(
    Cohort = prefix,
    Root_celltype = root_celltype,
    Root_state = root_state,
    N_cells = ncol(cds),
    N_ordering_genes = length(ordering_genes)
  )
  
  write.csv(
    root_info,
    paste0(prefix, "_root_state.csv"),
    row.names = FALSE
  )
  
  p_state <- plot_cell_trajectory(
    cds,
    color_by = "State"
  ) +
    ggtitle(paste0(prefix, " - State"))
  
  p_pseudotime <- plot_cell_trajectory(
    cds,
    color_by = "Pseudotime"
  ) +
    ggtitle(paste0(prefix, " - Pseudotime"))
  
  p_celltype <- plot_cell_trajectory(
    cds,
    color_by = "celltype_epi"
  ) +
    ggtitle(paste0(prefix, " - Cell type"))
  
  p_tissue <- plot_cell_trajectory(
    cds,
    color_by = "tissue"
  ) +
    ggtitle(paste0(prefix, " - Tissue"))
  
  ggsave(
    paste0(prefix, "_trajectory_State.pdf"),
    p_state,
    width = 6,
    height = 5
  )
  
  ggsave(
    paste0(prefix, "_trajectory_Pseudotime.pdf"),
    p_pseudotime,
    width = 6,
    height = 5
  )
  
  ggsave(
    paste0(prefix, "_trajectory_Celltype.pdf"),
    p_celltype,
    width = 6,
    height = 5
  )
  
  ggsave(
    paste0(prefix, "_trajectory_Tissue.pdf"),
    p_tissue,
    width = 6,
    height = 5
  )
  
  p_all <- (p_state | p_pseudotime) / (p_celltype | p_tissue)
  
  ggsave(
    paste0(prefix, "_trajectory_all.pdf"),
    p_all,
    width = 12,
    height = 10
  )
  
  write.csv(
    pData(cds),
    paste0(prefix, "_cell_pseudotime_metadata.csv"),
    row.names = TRUE
  )
  
  message("Finished: ", prefix)
  message("")
  
  return(cds)
}

monocle_EOPC <- run_monocle2(
  seu = epi_use,
  tissue_keep = c("EOPC", "young_normal"),
  prefix = "EOPC",
  root_celltype = "club",
  nfeatures = 2000
)

monocle_LOPC <- run_monocle2(
  seu = epi_use,
  tissue_keep = c("LOPC", "old_normal"),
  prefix = "LOPC",
  root_celltype = "club",
  nfeatures = 2000
)

saveRDS(
  monocle_EOPC,
  "monocle_EOPC_final.rds"
)

saveRDS(
  monocle_LOPC,
  "monocle_LOPC_final.rds"
)