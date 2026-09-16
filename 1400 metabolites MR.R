##1400 metabolites code

#BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
#BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
library(MRPRESSO)
library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(readxl)
library(writexl)
library(ieugwasr)
library(plinkbinr)
library(coloc)
library(locuscomparer)
library("MungeSumstats")
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library("BSgenome.Hsapiens.1000genomes.hs37d5")
setwd("E:\\华大单细胞课题\\1400代谢物\\1400代谢物数据1e-5")
FileNames <-list.files(paste0(getwd()),pattern=".csv")
exp_dat_ids <- FileNames
exps <- FileNames
# 结局文件放到工作路径下
setwd("D:\\孟德尔随机化法\\血清代谢物课程（壹隅_）\\血清数据新新新\\血清数据新新新\\cancer\\泌尿系统肿瘤及其亚型")
FileNames2 <-list.files("D:\\孟德尔随机化法\\血清代谢物课程（壹隅_）\\血清数据新新新\\血清数据新新新\\cancer\\泌尿系统肿瘤及其亚型",pattern=".csv")
FileNames2 <- FileNames2[13]
outcomeids <- FileNames2
out_comes <- FileNames2
qaqR <- 1
for (qaqR in 1:length(FileNames2)) {{ 
  setwd("D:\\孟德尔随机化法\\血清代谢物课程（壹隅_）\\血清数据新新新\\血清数据新新新\\cancer\\泌尿系统肿瘤及其亚型")
  out<- fread(outcomeids[qaqR],header = T)
  #out<- read.csv("prostate_cancer.csv")
  setwd("D:\\孟德尔随机化法\\smr-1.3.1-win-x86_64\\smr-1.3.1-win-x86_64")
  #######以下为循环代码，不需要进行更改######################
  qaq <- 1
  for (qaq in 1:length(exp_dat_ids)) { 
    setwd("E:\\华大单细胞课题\\1400代谢物\\1400代谢物数据1e-5")
    # setwd("D:\\孟德尔随机化法\\smr-1.3.1-win-x86_64\\smr-1.3.1-win-x86_64")
    d3 <- fread(paste0(getwd(), "/", FileNames[qaq]),quote = "\"") 
    if (ncol(d3) == 12) {
      original_colnames <- as.character(colnames(d3))  # 当前列名（合并后是11列）
      merged_col <- paste0(d3[[10]], d3[[11]])
      d3[, (10) := merged_col]  # 替换第10列
      # 删除第11列
      d3 <- d3[,-11]
      # 更新列名：改为原第2列到第12列的列名（跳过原第1列）
      new_colnames <- original_colnames[2:12]  # 原第2列到第12列的列名
      setnames(d3, new_colnames)  # 直接赋予新列名
    }
    exp <- exps[qaq]
    d3$SNP <- str_replace(d3$SNP, "\\*$", "") 
    # d3<- try(fread(paste0(getwd(),"/",FileNanmes[qaq]),fill=TRUE),silent = T)
    d3<-subset(d3,d3$SNP!="")
    d3 <- subset(d3,d3$pval<5e-6)
    d3$name <- FileNames[qaq]
    #d3$N <- 70000
    #rm(d1)
    #d3<-d2[,c(1,2,3,4,8,9,10,16)]
    #names(d3)[names(d3) == 'SNP'] <- 'SNP'
    d3$N <- d3$samplesize
    #d3$N <- 31684
    d3 <- as.data.frame(d3)
    d3<-format_data(d3,
                    type="exposure",
                    phenotype_col = "Phenotype",
                    snp_col = "SNP",
                    beta_col = "beta",
                    se_col = "se",
                    pval_col = "pval",
                    samplesize_col = "N",
                    eaf_col = "eaf",
                    effect_allele_col = "effect_allele",
                    other_allele_col = "other_allele")
    
    #exp_data <- clump_data(d3,clump_kb = 500,clump_r2 = 0.01)
    
    library(ieugwasr) 
    
    if (T) {
      
      #补充ld_clump需要的三列
      d3$id <- d3$id.exposure #id列用于区分是单个或多个来源的gwas数据，clump按照单个gwas数据进行
      d3$rsid <- d3$SNP
      d3$pval <- d3$pval.exposure
      
      # 执行ld_clump
      # bfile指定参考文件的路径
      bfile <- "D:\\孟德尔随机化法\\血清代谢物课程（壹隅_）\\血清数据新新新\\g1000_eur\\g1000_eur"
      d4 <- ld_clump(d3,
                     # get_plink_exe()
                     plink_bin = get_plink_exe(),
                     bfile = bfile ,
                     clump_kb = 10000, clump_r2 = 0.001)
      #  运行结果：                 
      # Removing 69 of 92 variants due to LD with other variants or absence from LD reference panel
      
      
      # 移除重复的列
      #d4$rsid <- NULL
      #d4$pval <- NULL
    }
    
    exp_data<-subset(d3,SNP %in% d4$rsid) 
    # exp_data$beta.exposure <- exp_data$beta.exposure*-1
    #exp_data <- d3
    #size1 <- c(402195,412181,342438,396410)
    #size1 <- c(316386,315137,314296,382887,134299,163095,314419,140254,3301) 
    
    outcomeid <- outcomeids[qaqR]
    out_come <- out_comes[qaqR]
    
    if(length(exp_data[,1])>2){
      size1 <- 208631
      outcomeid <- out
      outcome_dat<-merge(exp_data,outcomeid,by.x = "SNP",by.y = "ID")#by.y=rsids
      outcome_dat$samplesize.outcome <- as.numeric(size1)
      outcome_dat$trait <- FileNames2[qaqR]
      
      #outcome_dat$beta = log(outcome_dat$Odds_ratio)
      #outcome_dat$sebeta <- (-outcome_dat$beta + outcome_dat$ci_upper)/1.96
      outcome_dat$beta_outcome=log(outcome_dat$OR)
      outcome_dat$se=sqrt(((outcome_dat$beta_outcome^2)/qchisq(outcome_dat$P,1,lower.tail=F)))
      write.csv(outcome_dat,file = "d.csv")
      outcome_dat <- subset(outcome_dat,outcome_dat$P>5e-6)
      
      out_data <- read_outcome_data(
        snps = exp_data$SNP,
        phenotype_col = "trait",
        filename = "d.csv",
        sep = ",",
        eaf_col = "A1_FREQ",
        samplesize_col = "samplesize.outcome",
        snp_col = "SNP",
        beta_col = "beta_outcome",
        se_col = "se",
        effect_allele_col = "ALT",
        other_allele_col = "REF",
        pval_col = "P")
      
      if(length(out_data[,1])>0){  
        dat <- TwoSampleMR::harmonise_data(
          exposure_dat = exp_data,
          outcome_dat = out_data)
        
        ####回文的直接去除
        dat <-subset(dat,dat$mr_keep==TRUE)
        
        
        #计算F值和R2
        get_f<-function(dat,F_value=10){
          log<-is.na(dat$eaf.exposure)
          log<-unique(log)
          if(length(log)==1)
          {if(log==TRUE){
            print("数据不包含eaf，无法计算F统计量")
            return(dat)}
          }
          if(is.null(dat$beta.exposure[1])==T || is.na(dat$beta.exposure[1])==T){print("数据不包含beta，无法计算F统计量")
            return(dat)}
          if(is.null(dat$se.exposure[1])==T || is.na(dat$se.exposure[1])==T){print("数据不包含se，无法计算F统计量")
            return(dat)}
          if(is.null(dat$samplesize.exposure[1])==T || is.na(dat$samplesize.exposure[1])==T){print("数据不包含samplesize(样本量)，无法计算F统计量")
            return(dat)}
          
          
          if("FALSE"%in%log && is.null(dat$beta.exposure[1])==F && is.na(dat$beta.exposure[1])==F && is.null(dat$se.exposure[1])==F && is.na(dat$se.exposure[1])==F && is.null(dat$samplesize.exposure[1])==F && is.na(dat$samplesize.exposure[1])==F){
            R2<-(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))/((2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))+(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$se.exposure^2)*dat$samplesize.exposure))
            F<- (dat$samplesize.exposure-2)*R2/(1-R2)
            dat$R2<-R2
            dat$F<-F
            dat<-subset(dat,F>F_value)
            return(dat)
          }
        }
        
        dat <- get_f(dat, F_value = 10)
        
        
        res=TwoSampleMR::mr(dat,method_list= c("mr_ivw" ,
                                               "mr_weighted_median" ,
                                               "mr_egger_regression",
                                               "mr_simple_mode",
                                               "mr_weighted_mode",
                                               "mr_wald_ratio"))
        ##贝叶斯孟德尔随机化
        # myBWMR <- BWMR(gammahat = dat$beta.exposure,                 
        # Gammahat = dat$beta.outcome,                 
        #  sigmaX = dat$se.exposure,                 
        # sigmaY = dat$se.outcome) 
        
        
        print(paste0(out_come,"_SNP数_",res$nsnp[1]))
        
        try(results <- TwoSampleMR::generate_odds_ratios(res),silent = TRUE)
        setwd("E:\\华大单细胞课题\\1400代谢物")
        
        try(results$estimate <- paste0(
          format(round(results$or, 2), nsmall = 2), " (", 
          format(round(results$or_lci95, 2), nsmall = 2), "-",
          format(round(results$or_uci95, 2), nsmall = 2), ")"),silent = TRUE)
        
        resdata <- dat
        dir.create(path = outcomeids[qaqR])
        openxlsx::write.xlsx(dat,file = paste0(outcomeids[qaqR],"/",exp,"-dat.xlsx"), rowNames = FALSE)
        
        names(resdata)
        Assumption13 <- subset(resdata,mr_keep==TRUE,
                               select = c("SNP","pval.exposure",
                                          "pval.outcome", #"F_statistic",
                                          "mr_keep"))
        
        try( openxlsx::write.xlsx(x = list(
          "main"=results,
          "Assumption13"=Assumption13),
          overwrite = TRUE,
          paste0(outcomeids[qaqR],"/",exp,"-res.xlsx")),silent = TRUE)
        
      }}
    if(length(dat[,1])>2){
      res_hete <- TwoSampleMR::mr_heterogeneity(dat)
      res_plei <- TwoSampleMR::mr_pleiotropy_test(dat)
      try(res_leaveone <- mr_leaveoneout(dat),silent = TRUE)  # 
      
      ######steiger检验######
      dat$r.exposure <- get_r_from_bsen(b = dat$beta.exposure,
                                        dat$se.exposure,
                                        dat$samplesize.exposure)
      dat$r.outcome <- get_r_from_bsen(b = dat$beta.outcome,
                                       dat$se.outcome,
                                       dat$samplesize.outcome)
      try(res_steiger <- mr_steiger(
        p_exp = dat$pval.exposure,
        p_out = dat$pval.outcome,
        n_exp = dat$samplesize.exposure,
        n_out = dat$samplesize.outcome,
        r_exp = dat$r.exposure,
        r_out = dat$r.outcome
      ),silent = TRUE)
      try(res_steiger <- directionality_test(dat),silent = TRUE)
      
      
      
      
      #若需运行MR—presso,将最左侧7个#删去即可 
      #res_presso <- TwoSampleMR::run_mr_presso(dat,
      #NbDistribution = 100)
      # [["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
      #sink(paste0("77种糖基化结局/",out_come,"_PRESSO.txt"),
      #append=FALSE,split = FALSE) 
      #print(res_presso)
      #sink()
      # print(res_presso)
      
      
      
      setwd("E:\\华大单细胞课题\\1400代谢物")
      p1 <- mr_scatter_plot(res, dat)
      try(p1[[1]],silent = TRUE)
      pdf(paste0(outcomeids[qaqR],"/",exp,"_scatter.pdf"))
      try(print(p1[[1]]),silent = TRUE)
      dev.off()
      
      try(res_single <- mr_singlesnp(dat),silent = TRUE)
      try(p2 <- mr_forest_plot(res_single),silent = TRUE)
      pdf(paste0(outcomeids[qaqR],"/",exp,"_forest.pdf"))
      try(print(p2[[1]]),silent = TRUE)
      dev.off()
      
      try(p3 <- mr_funnel_plot(res_single),silent = TRUE)
      pdf(paste0(outcomeids[qaqR],"/",exp,"_funnel.pdf"))
      try(print(p3[[1]]),silent=TRUE)
      dev.off()
      
      try(res_loo <- mr_leaveoneout(dat),silent = TRUE)
      pdf(paste0(outcomeids[qaqR],"/",exp,"_leave_one_out.pdf"))
      try(print(mr_leaveoneout_plot(res_loo)),silent = TRUE)
      dev.off()
      
      
      library(magrittr)
      res3 <- res[1:3,]
      
      
      
      # 转换成论文格式
      library(magrittr)
      # Main result 
      try(res4 <- tidyr::pivot_wider(
        res3,names_from ="method",names_vary = "slowest",
        values_from = c("b","se","pval") ),silent = TRUE)
      # Heterogeneity statistics
      try(res_hete2 <- tidyr::pivot_wider(
        res_hete,names_from ="method",names_vary = "slowest",
        values_from = c("Q","Q_df","Q_pval") ) %>% 
          dplyr::select( -id.exposure,-id.outcome,-outcome,-exposure),silent = TRUE)
      # Horizontal pleiotropy
      try(res_plei2 <- dplyr::select(res_plei,
                                     egger_intercept,se,pval),silent = TRUE)
      
      ##
      try(res_steiger2 <- dplyr::select(res_steiger,
                                        correct_causal_direction,steiger_pval),silent = TRUE)
      
      
      # Merge
      res_ALL <- cbind(res4, res_hete2, res_plei2,res_steiger2)
      
      write.csv(res_ALL,file = paste0(outcomeids[qaqR],"/",exp,".csv"), row.names = FALSE)
      
    }}
  
}}
#导出合并的结果
library(openxlsx)
setwd("E:\\华大单细胞课题\\1400代谢物\\EOPC.csv")
Filenames3 <- list.files("E:\\华大单细胞课题\\1400代谢物\\EOPC.csv",pattern = "res.xlsx")
df <- read.xlsx("GCST90199621_buildGRCh38.csv-res.xlsx")
df$pheno <- Filenames3[1]
for (i in 2:length(Filenames3)) {{
  df1 <- read.xlsx(Filenames3[i])
  df1$pheno <- Filenames3[i]
  df <- rbind(df1,df)
}}

library(data.table)
library(openxlsx)
setwd("E:\\华大单细胞课题\\1400代谢物")
ID <- read.xlsx("accessionId_reportedTrait.xlsx")
df <- read.xlsx("all associations.xlsx",sheet = 2)

df2 <- merge(df,ID,by.x="pheno",by.y="accessionId")
df2$reportedTrait <- sub(" levels$", "", df2$reportedTrait)
subpathway <- read.xlsx("subpathway.xlsx",sheet=2)
subpathway$Metabolites <- gsub("\\*", "", subpathway$Metabolites)
write.xlsx(subpathway,"subpathway2.xlsx")


library(dplyr)
library(ggplot2)
library(patchwork)
library(readxl)
library(ggthemes)

df2_temp <- df2
subpathway_temp <- subpathway

# 将所有键转换为小写（或大写）
df2_temp$reportedTrait_lower <- tolower(df2_temp$reportedTrait)
subpathway_temp$Metabolites_lower <- tolower(subpathway_temp$Metabolites)

# 使用转换后的键进行合并
df3 <- merge(df2_temp, subpathway_temp, 
             by.x = "reportedTrait_lower", 
             by.y = "Metabolites_lower")
df3 <- df3[!duplicated(df3),]
# 移除临时创建的列（如果需要）
df3$reportedTrait_lower <- NULL
df3$Metabolites_lower <- NULL
write.xlsx(df3,"df3.xlsx")

setwd("D:\\孟德尔随机化法\\血清代谢物课程（壹隅_）")
id <- read.csv("代谢id.csv")
df <- merge(df,id,by.x="exposure",by.y="Metabolite.ID",all.y = TRUE)
write.xlsx(df,"结局.xlsx")

##
##画一个森林图
library(openxlsx)
library(data.table)
library(ggplot2)
library(ggthemes)
setwd("E:\\华大单细胞课题\\1400代谢物")
dataset <- read.xlsx("all associations.xlsx",sheet=4)
dataset$index <- factor(dataset$index, levels = unique(dataset$index))
dataset <- subset(dataset,dataset$Estimate>1)
b1 <- ggplot(dataset, aes(Estimate, index)) +
  geom_point(size = 1, aes(color = comp, shape = comp)) +
  geom_errorbarh(
    aes(xmax = upper.limit, xmin = lower.limit, color = comp),
    size = 1, height = 0
  ) +
  scale_shape_manual(values =  c("circle", "square", "triangle", "diamond", "cross", "star", "asterisk")) +
  scale_color_manual(values = c(
    "#FA7F6F", "#8ECFC9", "#FFBE7A",
    "#82B0D2", "#BEB8DC", "#E7DAD2","#F60000"
  )) +
  scale_x_continuous(limits = c(0, 4)) +
  geom_vline(
    aes(xintercept = 1),
    color = "gray", linetype = "dashed", size = 0.8
  ) +
  geom_text(aes(label = sig, x = 0.5)) +
  geom_text(aes(label = comp, x = 3.5)) +
  labs(x = "ln RR (%)", y = NULL) +
  theme_few() +
  theme(
    axis.text.x = element_text(
      size = 10, 
      color = "black",
      margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
    ),
    axis.text.y = element_text(
      size = 10, 
      color = "black",
      margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
    ),
    title = element_text(size = 10),
    panel.grid = element_blank(),
    axis.ticks.length = unit(0.1, "cm"),
    legend.position = "none"
  )
b1
pdf("forest plot.pdf",width = 8,height=9)
print(b1)
dev.off()
