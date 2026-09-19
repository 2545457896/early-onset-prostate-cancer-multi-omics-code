# ==============================================================================
# EOPC / Prostate Cancer
# scATOMIC annotation pipeline
# ==============================================================================


# ==============================================================================
# 0. 加载包
# ==============================================================================

library(scATOMIC)
library(Seurat)
library(Matrix)
library(dplyr)
library(openxlsx)
library(data.table)
library(randomForest)
library(caret)
library(parallel)
library(reticulate)
library(Rmagic)
library(agrmt)
library(cutoff.scATOMIC)
library(copykat)
library(ggplot2)


# ==============================================================================
# 1. 路径设置
# ==============================================================================

setwd("/mnt/DATA/home/zqy1234560915/EOPC")

input_file <- "/mnt/DATA/home/zqy1234560915/EOPC/scRNA_integrated_clustered.rds"

output_dir <- "/mnt/DATA/home/zqy1234560915/EOPC/scATOMIC_PRAD"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

cat("当前工作目录：", getwd(), "\n")
cat("输入文件：", input_file, "\n")
cat("输出目录：", output_dir, "\n")


# ==============================================================================
# 2. Python / MAGIC 环境
# ==============================================================================

Sys.setenv(
  RETICULATE_PYTHON =
    "/mnt/DATA/home/zqy1234560915/miniconda3/envs/scatomic/bin/python"
)

library(reticulate)

pymagic <- import("magic")

cat("Python：", py_config()$python, "\n")


# ==============================================================================
# 3. 检查 create_summary_matrix
#
# 必须使用你前面修改后的版本：
# create_summary_matrix(..., known_cancer_type = NULL)
# ==============================================================================

if (!exists("create_summary_matrix", mode = "function")) {
  
  stop(
    "\n没有找到 create_summary_matrix()。\n",
    "请先运行你之前修改后的 create_summary_matrix() 函数定义部分，",
    "然后再运行本段代码。\n"
  )
}

fun_args <- names(formals(create_summary_matrix))

if (!"known_cancer_type" %in% fun_args) {
  
  stop(
    "\n当前 create_summary_matrix() 没有 known_cancer_type 参数。\n",
    "请使用我们之前修改后的版本。\n"
  )
}

cat("create_summary_matrix 检查通过。\n")


# ==============================================================================
# 4. 读取 Seurat 对象
# ==============================================================================

cat("\n========== 读取 Seurat 对象 ==========\n")

seurat_object <- readRDS(input_file)

cat("对象读取完成\n")
cat("细胞数：", ncol(seurat_object), "\n")
cat("基因数：", nrow(seurat_object), "\n")
cat("Assays：", paste(Assays(seurat_object), collapse = ", "), "\n")


# ==============================================================================
# 5. 检查 RNA assay
# ==============================================================================

if (!"RNA" %in% Assays(seurat_object)) {
  
  stop("Seurat 对象中不存在 RNA assay。")
}

DefaultAssay(seurat_object) <- "RNA"

cat("\nRNA assay class：")
print(class(seurat_object[["RNA"]]))

if (inherits(seurat_object[["RNA"]], "Assay5")) {
  
  cat("\nRNA layers：\n")
  
  print(
    Layers(seurat_object[["RNA"]])
  )
}


# ==============================================================================
# 6. 检查 metadata
# ==============================================================================

cat("\n========== Metadata ==========\n")

print(
  colnames(seurat_object@meta.data)
)

cat("\n前几行 metadata：\n")

print(
  head(seurat_object@meta.data)
)


# ==============================================================================
# 7. 检查 orig.ident
# ==============================================================================

if (!"orig.ident" %in% colnames(seurat_object@meta.data)) {
  
  stop(
    "\nmetadata 中没有 orig.ident。\n",
    "请检查哪个变量代表独立患者/样本，然后把 split.by 改成对应变量。\n"
  )
}


cat("\n========== 每个样本细胞数量 ==========\n")

sample_table <- sort(
  table(seurat_object$orig.ident),
  decreasing = TRUE
)

print(sample_table)

cat("\n样本总数：", length(sample_table), "\n")


# 保存样本统计
write.xlsx(
  data.frame(
    Sample_ID = names(sample_table),
    Cell_number = as.numeric(sample_table)
  ),
  file.path(
    output_dir,
    "PRAD_sample_cell_numbers.xlsx"
  ),
  rowNames = FALSE
)


# ==============================================================================
# 8. 基因名简单检查
# ==============================================================================

gene_names <- rownames(seurat_object)

ensg_prop <- mean(
  grepl("^ENSG", gene_names)
)

cat(
  "\nENSG 格式基因比例：",
  round(ensg_prop * 100, 2),
  "%\n"
)

if (ensg_prop > 0.5) {
  
  warning(
    "超过一半基因看起来是 Ensembl ID。",
    "scATOMIC 通常需要 gene symbol，请确认基因名。"
  )
}

cat("\n前20个基因：\n")
print(head(gene_names, 20))


# ==============================================================================
# 9. 按样本拆分
#
# 非常重要：
# scATOMIC 不建议把不同患者直接一起进行 malignant/normal 判断
# ==============================================================================

cat("\n========== 按 orig.ident 拆分 ==========\n")

sample_list <- SplitObject(
  seurat_object,
  split.by = "orig.ident"
)

cat(
  "成功拆成",
  length(sample_list),
  "个样本\n"
)


# ==============================================================================
# 10. Seurat v4 / v5 通用 counts 提取函数
# ==============================================================================

get_raw_counts <- function(obj) {
  
  DefaultAssay(obj) <- "RNA"
  
  # -------------------------
  # Seurat v5
  # -------------------------
  
  if (inherits(obj[["RNA"]], "Assay5")) {
    
    current_layers <- Layers(obj[["RNA"]])
    
    count_layers <- grep(
      "^counts",
      current_layers,
      value = TRUE
    )
    
    if (length(count_layers) == 0) {
      
      stop(
        "RNA assay 中没有 counts layer。"
      )
    }
    
    # 如果有多个 counts.xxx layer
    # 先 JoinLayers
    if (length(count_layers) > 1) {
      
      cat(
        "检测到多个 counts layers：",
        paste(count_layers, collapse = ", "),
        "\n"
      )
      
      obj <- JoinLayers(
        obj,
        assay = "RNA"
      )
    }
    
    counts <- LayerData(
      obj,
      assay = "RNA",
      layer = "counts"
    )
    
  } else {
    
    # -------------------------
    # Seurat v4
    # -------------------------
    
    counts <- GetAssayData(
      obj,
      assay = "RNA",
      slot = "counts"
    )
  }
  
  
  counts <- as(
    counts,
    "dgCMatrix"
  )
  
  
  return(
    list(
      object = obj,
      counts = counts
    )
  )
}


# ==============================================================================
# 11. 创建保存目录
# ==============================================================================

sample_result_dir <- file.path(
  output_dir,
  "each_sample"
)

dir.create(
  sample_result_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ==============================================================================
# 12. 逐个样本运行 scATOMIC
# ==============================================================================

sample_results_list <- list()

failed_samples <- character()

skipped_samples <- character()


for (sample_name in names(sample_list)) {
  
  
  cat("\n")
  cat("============================================================\n")
  cat("样本：", sample_name, "\n")
  cat("============================================================\n")
  
  
  # --------------------------------------------------------------------------
  # 12.1 当前样本
  # --------------------------------------------------------------------------
  
  current_sample_obj <- sample_list[[sample_name]]
  
  cat(
    "原始细胞数：",
    ncol(current_sample_obj),
    "\n"
  )
  
  
  # --------------------------------------------------------------------------
  # 12.2 获取 raw counts
  # --------------------------------------------------------------------------
  
  counts_res <- get_raw_counts(
    current_sample_obj
  )
  
  current_sample_obj <- counts_res$object
  
  sparse_matrix <- counts_res$counts
  
  rm(counts_res)
  
  gc()
  
  
  cat(
    "Count matrix：",
    nrow(sparse_matrix),
    "genes ×",
    ncol(sparse_matrix),
    "cells\n"
  )
  
  
  # --------------------------------------------------------------------------
  # 12.3 细胞太少跳过
  # --------------------------------------------------------------------------
  
  if (ncol(sparse_matrix) < 30) {
    
    cat(
      "⚠️ 样本 ",
      sample_name,
      " 只有 ",
      ncol(sparse_matrix),
      " 个细胞，跳过。\n",
      sep = ""
    )
    
    skipped_samples <- c(
      skipped_samples,
      sample_name
    )
    
    next
  }
  
  
  # --------------------------------------------------------------------------
  # 12.4 去掉完全不表达的基因
  # --------------------------------------------------------------------------
  
  keep_gene <- Matrix::rowSums(
    sparse_matrix
  ) > 0
  
  sparse_matrix <- sparse_matrix[
    keep_gene,
    ,
    drop = FALSE
  ]
  
  cat(
    "去除 zero genes 后：",
    nrow(sparse_matrix),
    "genes\n"
  )
  
  
  # --------------------------------------------------------------------------
  # 12.5 scATOMIC
  # --------------------------------------------------------------------------
  
  cat("\n>>> 开始 run_scATOMIC\n")
  
  
  cell_predictions <- tryCatch(
    
    {
      
      run_scATOMIC(
        sparse_matrix
      )
      
    },
    
    error = function(e) {
      
      cat(
        "\n❌ run_scATOMIC 失败：",
        conditionMessage(e),
        "\n"
      )
      
      return(NULL)
    }
    
  )
  
  
  if (is.null(cell_predictions)) {
    
    failed_samples <- c(
      failed_samples,
      sample_name
    )
    
    next
  }
  
  
  cat(">>> run_scATOMIC 完成\n")
  
  
  # --------------------------------------------------------------------------
  # 12.6 保存原始 prediction
  #
  # 防止后面 summary 阶段失败后必须重新跑 random forest
  # --------------------------------------------------------------------------
  
  saveRDS(
    cell_predictions,
    file.path(
      sample_result_dir,
      paste0(
        sample_name,
        "_raw_scATOMIC_prediction.rds"
      )
    )
  )
  
  
  # --------------------------------------------------------------------------
  # 12.7 create_summary_matrix
  #
  # 这里明确告诉算法：
  #
  #    Prostate Cancer Cell
  #
  # 也就是 PRAD
  # --------------------------------------------------------------------------
  
  cat("\n>>> 开始 create_summary_matrix\n")
  
  
  sample_res <- tryCatch(
    
    {
      
      create_summary_matrix(
        
        raw_counts = sparse_matrix,
        
        prediction_list = cell_predictions,
        
        use_CNVs = FALSE,
        
        modify_results = TRUE,
        
        mc.cores = 1,
        
        min_prop = 0.5,
        
        breast_mode = FALSE,
        
        fine_grained_T = TRUE,
        
        confidence_cutoff = TRUE,
        
        pan_cancer = FALSE,
        
        cancer_confidence = "default",
        
        normal_tissue = FALSE,
        
        low_res_mode = FALSE,
        
        known_cancer_type = "Prostate Cancer Cell"
        
      )
      
    },
    
    error = function(e) {
      
      cat(
        "\n❌ create_summary_matrix 失败：",
        conditionMessage(e),
        "\n"
      )
      
      return(NULL)
    }
    
  )
  
  
  if (is.null(sample_res)) {
    
    failed_samples <- c(
      failed_samples,
      sample_name
    )
    
    next
  }
  
  
  cat(">>> create_summary_matrix 完成\n")
  
  
  # --------------------------------------------------------------------------
  # 12.8 添加样本信息
  # --------------------------------------------------------------------------
  
  sample_res$Cell_Barcode <- rownames(
    sample_res
  )
  
  sample_res$Sample_ID <- sample_name
  
  
  # --------------------------------------------------------------------------
  # 12.9 保存当前样本结果
  # --------------------------------------------------------------------------
  
  saveRDS(
    sample_res,
    file.path(
      sample_result_dir,
      paste0(
        sample_name,
        "_scATOMIC_PRAD_result.rds"
      )
    )
  )
  
  
  write.xlsx(
    sample_res,
    file.path(
      sample_result_dir,
      paste0(
        sample_name,
        "_scATOMIC_PRAD_result.xlsx"
      )
    ),
    rowNames = FALSE
  )
  
  
  # --------------------------------------------------------------------------
  # 12.10 打印分类结果
  # --------------------------------------------------------------------------
  
  cat(
    "\n当前样本最终 annotation：\n"
  )
  
  print(
    sort(
      table(sample_res$scATOMIC_pred),
      decreasing = TRUE
    )
  )
  
  
  # --------------------------------------------------------------------------
  # 12.11 放进总列表
  # --------------------------------------------------------------------------
  
  sample_results_list[[sample_name]] <- sample_res
  
  
  # --------------------------------------------------------------------------
  # 12.12 清内存
  # --------------------------------------------------------------------------
  
  rm(
    current_sample_obj,
    sparse_matrix,
    cell_predictions,
    sample_res
  )
  
  gc()
}


# ==============================================================================
# 13. 检查成功样本
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("scATOMIC 全部样本运行结束\n")
cat("============================================================\n")

cat(
  "成功样本数：",
  length(sample_results_list),
  "\n"
)

cat(
  "跳过样本数：",
  length(skipped_samples),
  "\n"
)

cat(
  "失败样本数：",
  length(unique(failed_samples)),
  "\n"
)


if (length(sample_results_list) == 0) {
  
  stop(
    "没有任何样本成功生成结果，请检查上面的报错。"
  )
}


# ==============================================================================
# 14. 不同样本结果列补齐
# ==============================================================================

all_cols <- unique(
  unlist(
    lapply(
      sample_results_list,
      colnames
    )
  )
)


sample_results_list_padded <- lapply(
  
  sample_results_list,
  
  function(df) {
    
    missing_cols <- setdiff(
      all_cols,
      colnames(df)
    )
    
    if (length(missing_cols) > 0) {
      
      df[missing_cols] <- NA
    }
    
    
    # 保持相同列顺序
    df <- df[
      ,
      all_cols,
      drop = FALSE
    ]
    
    
    return(df)
  }
  
)


# ==============================================================================
# 15. 合并所有样本
# ==============================================================================

combined_results <- do.call(
  rbind,
  sample_results_list_padded
)


rownames(combined_results) <- combined_results$Cell_Barcode


cat(
  "\n总 annotation 细胞数：",
  nrow(combined_results),
  "\n"
)


# ==============================================================================
# 16. 查看最终 scATOMIC 分类
# ==============================================================================

cat("\n========== 最终 cell type ==========\n")

final_table <- sort(
  table(combined_results$scATOMIC_pred),
  decreasing = TRUE
)

print(final_table)


# 保存统计
final_table_df <- data.frame(
  
  Cell_type = names(final_table),
  
  Cell_number = as.numeric(
    final_table
  ),
  
  Proportion = as.numeric(final_table) /
    sum(final_table)
  
)


write.xlsx(
  final_table_df,
  file.path(
    output_dir,
    "PRAD_scATOMIC_celltype_summary.xlsx"
  ),
  rowNames = FALSE
)


# ==============================================================================
# 17. 保存完整结果
# ==============================================================================

saveRDS(
  combined_results,
  file.path(
    output_dir,
    "PRAD_scATOMIC_all_cells.rds"
  )
)


# Excel 最大大约 104 万行
if (nrow(combined_results) < 1000000) {
  
  write.xlsx(
    combined_results,
    file.path(
      output_dir,
      "PRAD_scATOMIC_all_cells.xlsx"
    ),
    rowNames = FALSE
  )
  
} else {
  
  cat(
    "\n细胞数接近 Excel 行数限制，",
    "不输出总 Excel，仅保存 RDS。\n"
  )
}


# ==============================================================================
# 18. 将 scATOMIC annotation 加回原始 Seurat 对象
# ==============================================================================

cat("\n========== 将结果写回 Seurat ==========\n")


# 检查 barcode 是否匹配
matched_cells <- intersect(
  colnames(seurat_object),
  combined_results$Cell_Barcode
)


cat(
  "原始 Seurat 细胞数：",
  ncol(seurat_object),
  "\n"
)

cat(
  "scATOMIC 成功 annotation：",
  nrow(combined_results),
  "\n"
)

cat(
  "成功匹配回 Seurat：",
  length(matched_cells),
  "\n"
)


# ==============================================================================
# 19. 创建 metadata
# ==============================================================================

annotation_meta <- combined_results[
  matched_cells,
  ,
  drop = FALSE
]


# 只选比较重要的列写回去
wanted_cols <- c(
  
  "scATOMIC_pred",
  
  "classification_confidence",
  
  "layer_1",
  "layer_2",
  "layer_3",
  "layer_4",
  "layer_5",
  "layer_6",
  
  "pan_cancer_cluster"
  
)


wanted_cols <- intersect(
  wanted_cols,
  colnames(annotation_meta)
)


annotation_meta <- annotation_meta[
  ,
  wanted_cols,
  drop = FALSE
]


# 为了防止与你已有 metadata 重名
colnames(annotation_meta) <- paste0(
  "scATOMIC_",
  colnames(annotation_meta)
)


# ==============================================================================
# 20. AddMetaData
# ==============================================================================

seurat_object <- AddMetaData(
  seurat_object,
  metadata = annotation_meta
)


# ==============================================================================
# 21. 看结果
# ==============================================================================

cat("\n新增 metadata：\n")

print(
  grep(
    "^scATOMIC_",
    colnames(seurat_object@meta.data),
    value = TRUE
  )
)


if ("scATOMIC_scATOMIC_pred" %in%
    colnames(seurat_object@meta.data)) {
  
  cat("\n最终 annotation：\n")
  
  print(
    sort(
      table(
        seurat_object$scATOMIC_scATOMIC_pred,
        useNA = "ifany"
      ),
      decreasing = TRUE
    )
  )
}


# ==============================================================================
# 22. 为了方便使用，额外创建一个简单变量：major_celltype
# ==============================================================================

seurat_object$major_celltype <-
  seurat_object$scATOMIC_scATOMIC_pred


# 没有成功跑 scATOMIC 的细胞保留 NA
table(
  seurat_object$major_celltype,
  useNA = "ifany"
)


# ==============================================================================
# 23. 保存最终 Seurat
# ==============================================================================

final_rds <- file.path(
  output_dir,
  "scRNA_integrated_clustered_scATOMIC_PRAD.rds"
)


saveRDS(
  seurat_object,
  final_rds
)


cat("\n")
cat("============================================================\n")
cat("全部完成\n")
cat("============================================================\n")

cat(
  "最终 Seurat：",
  final_rds,
  "\n"
)

cat(
  "完整 annotation：",
  file.path(
    output_dir,
    "PRAD_scATOMIC_all_cells.rds"
  ),
  "\n"
)


# ==============================================================================
# 24. 简单检查前列腺癌细胞数量
# ==============================================================================

cat("\n========== Prostate Cancer Cell ==========\n")

print(
  table(
    seurat_object$major_celltype ==
      "Prostate Cancer Cell",
    useNA = "ifany"
  )
)


cat("\n每个样本中的 Prostate Cancer Cell 数量：\n")

prostate_count <- seurat_object@meta.data %>%
  
  mutate(
    Cell = rownames(seurat_object@meta.data)
  ) %>%
  
  filter(
    major_celltype ==
      "Prostate Cancer Cell"
  ) %>%
  
  count(
    orig.ident,
    name = "Prostate_Cancer_Cell_number"
  ) %>%
  
  arrange(
    desc(Prostate_Cancer_Cell_number)
  )


print(prostate_count)


write.xlsx(
  prostate_count,
  file.path(
    output_dir,
    "Prostate_Cancer_Cell_number_per_sample.xlsx"
  ),
  rowNames = FALSE
)