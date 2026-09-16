memory.limit(9999999999999)
rm(list=ls())
set.seed(61234)
library(stringr)
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
memory.limit(9999999999999)
set.seed(61234)
library(data.table)
library(openxlsx)
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(monocle)
library(ggsci)
library(harmony)
library(presto)
library(ggpubr)
library(msigdbr)
library(clusterProfiler)
library(reticulate)
library(CytoTRACE2)
library(cowplot)
#library(dittoSeq)
library(viridis)
#library(NEPAL)
npg_colors <- pal_npg("nrc")(10)  
npg_palette_extended <- colorRampPalette(npg_colors)(13)  
library(openxlsx)
library(igraph)

##之前的整合方式


##现在的整合方式
##NC
setwd("/home/rstudio/single cell EOPC")
NC <- readRDS("EOPC_NC.rds")
NC$study <- "NC"

##WCH
WCH <- readRDS("WCH_BONE.rds")
WCH$study <- "WCH"


##NG
NG <- readRDS("nature aging.rds")
NG$study <- "NG"



##合并数据
scRNAlist <- merge(NC,list(WCH,NG))
scRNAlist <- JoinLayers(scRNAlist)
sample <- read.xlsx("sample.xlsx",sheet=2)
metadata <- scRNAlist@meta.data
metadata$cell <- rownames(metadata)
metadata <- metadata %>% 
  left_join(sample, by = c("patient" = "Patient.ID"))
rownames(metadata) <- metadata$cell
scRNAlist <- AddMetaData(scRNAlist,metadata = metadata)
scRNAlist <- subset(scRNAlist,subset=inclusion==1)

scRNAlist <- NormalizeData(scRNAlist)
scRNAlist <- FindVariableFeatures(scRNAlist,selection.method = 'vst',nfeatures = 2000)
scRNA.integrated <- scRNAlist
scRNA.integrated <-ScaleData(scRNA.integrated,verbose = FALSE)
scRNA.integrated <- RunPCA(scRNA.integrated,verbose = T,npcs = 50)
save(scRNA.integrated, file = "scRNA_integrated.Rdata")
scRNA.integrated <- RunHarmony(scRNA.integrated, group.by.vars = c("study","orig.ident"), plot_convergence = TRUE)
scRNA.integrated <- scRNA.integrated %>%
  RunUMAP(reduction = "harmony", dims = 1:50)%>%#
  FindNeighbors(reduction = "harmony", dims = 1:50) 
#saveRDS(scRNA.integrated,'scRNA_integrated_new.rds')
#save.image('20250421.rdata')
saveRDS(scRNA.integrated,'scRNA_integrated_new.rds')
scRNA.integrated <- readRDS("scRNA_integrated_new.rds")
#RunTSNE(reduction = "harmony", dims = 1:50)
for (res in c(0.5,0.8,1)) {
  scRNA.integrated=FindClusters(scRNA.integrated, #graph.name = "un_inte",
                                resolution = res, algorithm = 1)
}
Idents(scRNA.integrated) <- scRNA.integrated$RNA_snn_res.0.5

##添加average onset
scRNA.integrated$celltype <- case_when(
  scRNA.integrated$RNA_snn_res.0.5 %in% c(10, 2, 20, 22,24, 25, 26, 30,31, 32, 33, 35, 8, 9) ~ "Epithelium",
  scRNA.integrated$RNA_snn_res.0.5 %in% c(19, 27, 3) ~ "Endothelium",
  scRNA.integrated$RNA_snn_res.0.5 == 7 ~ "B cells",
  scRNA.integrated$RNA_snn_res.0.5 %in% c(0, 1, 12, 16, 21, 29, 36, 37, 38, 39, 40, 41, 42) ~ "T cells",
  scRNA.integrated$RNA_snn_res.0.5 == 13 ~ "NK cells",
  scRNA.integrated$RNA_snn_res.0.5 == 11 ~ "Mast cells",
  scRNA.integrated$RNA_snn_res.0.5 == 5 ~ "Perivascular",
  scRNA.integrated$RNA_snn_res.0.5 == 28 ~ "DC cells",
  scRNA.integrated$RNA_snn_res.0.5 %in% c(4, 34) ~ "Mono-macro",
  scRNA.integrated$RNA_snn_res.0.5 %in% c(6, 17) ~ "Myeloid",
  scRNA.integrated$RNA_snn_res.0.5 %in% c(15, 23) ~ "RBC",
  scRNA.integrated$RNA_snn_res.0.5 == 18 ~ "Plasma",
  scRNA.integrated$RNA_snn_res.0.5 == 14 ~ "Fibroblast",
  TRUE ~ "others"
)

#scRNA.integrated <- total
scRNA.integrated$tissue <- ifelse(scRNA.integrated$orig.ident=="SCG-PCA7-T-HG"|scRNA.integrated$orig.ident=="SCG-PCA10-T-LG"|scRNA.integrated$orig.ident=="SCG-PCA12-T-LG"|scRNA.integrated$orig.ident=="EOPC1"|scRNA.integrated$orig.ident=="EOPC2"|scRNA.integrated$orig.ident=="EOPC3"|scRNA.integrated$orig.ident=="EOPC4"|scRNA.integrated$orig.ident=="Rep-EOPC1"|scRNA.integrated$orig.ident=="Rep-EOPC2"|scRNA.integrated$orig.ident=="Rep-EOPC3"|scRNA.integrated$orig.ident=="Rep-EOPC4"|scRNA.integrated$orig.ident=="Rep-EOPC5"|scRNA.integrated$orig.ident=="MHSPC01BONE"|scRNA.integrated$orig.ident=="MHSPC02BONE","EOPC",ifelse(scRNA.integrated$orig.ident=="HP3"|scRNA.integrated$orig.ident=="HP3"|scRNA.integrated$orig.ident=="SCG-PCA12-N-LG","young_normal",ifelse(scRNA.integrated$orig.ident=="Healthy-PC1"|scRNA.integrated$orig.ident=="HP1"|scRNA.integrated$orig.ident=="HP2"|scRNA.integrated$orig.ident=="HP4"|scRNA.integrated$orig.ident=="SCG-PCA11-N-LG"|scRNA.integrated$orig.ident=="SCG-PCA12-N-LG"|scRNA.integrated$orig.ident=="SCG-PCA15-N-HG"|scRNA.integrated$orig.ident=="SCG-PCA17-N-LG"|scRNA.integrated$orig.ident=="SCG-PCA18-N-LG"|scRNA.integrated$orig.ident=="SCG-PCA19-N-HG"|scRNA.integrated$orig.ident=="SCG-PCA20-N-LG"|scRNA.integrated$orig.ident=="SCG-PCA21-N-LG"|scRNA.integrated$orig.ident=="SCG-PCA22-N-HG"|scRNA.integrated$orig.ident=="SCG-PCA24-N-CplusD"|scRNA.integrated$orig.ident=="SCG-PCA24-N-Rocky"|scRNA.integrated$orig.ident=="SCG-PCA3-N-LG"|scRNA.integrated$orig.ident=="SCG-PCA4-N-HG"|scRNA.integrated$orig.ident=="SCG-PCA5-N-LG"|scRNA.integrated$orig.ident=="SCG-PCA6-N-HG"|scRNA.integrated$orig.ident=="SCG-PCA9-N-LG","old_normal","LOPC")))
scRNA.integrated$lesion <- ifelse(scRNA.integrated$orig.ident=="MHSPC01BONE"|scRNA.integrated$orig.ident=="MHSPC02BONE"|scRNA.integrated$orig.ident=="MHSPC03BONE"|scRNA.integrated$orig.ident=="MHSPC04BONE"|scRNA.integrated$orig.ident=="MHSPC05BONE"|scRNA.integrated$orig.ident=="MHSPC06BONE","BONE",ifelse(scRNA.integrated$tissue=="young_normal"|scRNA.integrated$tissue=="old_normal","normal",ifelse(scRNA.integrated$orig.ident=="Rep-EOPC3"|scRNA.integrated$orig.ident=="Rep-EOPC1"|scRNA.integrated$orig.ident=="Rep-EOPC5"|scRNA.integrated$orig.ident=="Rep-LOPC1"|scRNA.integrated$orig.ident=="Rep-LOPC5"|scRNA.integrated$orig.ident=="SCG-PCA3-T-LG"|scRNA.integrated$orig.ident=="SCG-PCA5-T-LG"|scRNA.integrated$orig.ident=="SCG-PCA9-T-LG"|scRNA.integrated$orig.ident=="SCG-PCA10-T-LG"|scRNA.integrated$orig.ident=="SCG-PCA11-T-LG"|scRNA.integrated$orig.ident=="SCG-PCA12-T-LG"|scRNA.integrated$orig.ident=="SCG-PCA15-T-LG"|scRNA.integrated$orig.ident=="SCG-PCA17-T-LG"|scRNA.integrated$orig.ident=="SCG-PCA18-T-LG"|scRNA.integrated$orig.ident=="SCG-PCA20-T-LG"|scRNA.integrated$orig.ident=="SCG-PCA21-T-LG"|scRNA.integrated$orig.ident=="SCG-PCA24-T-LG","low grade","high grade")))


scRNA.integrated$cell <- rownames(scRNA.integrated@meta.data)
saveRDS(scRNA.integrated,"scRNA_integrated_clustered.rds")
