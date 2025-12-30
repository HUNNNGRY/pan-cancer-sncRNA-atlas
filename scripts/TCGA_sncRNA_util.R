setwd("/lulabdata3/baopengfei/projects/motif-RBP-EDA/")
options(scipen = 4,stringsAsFactors = F,digits = 5)
options(bedtools.path = "/lulabdata3/baopengfei/biosoft")
library(bedtoolsr)
library(plyr)      # Load plyr BEFORE dplyr to avoid conflicts
library(tidyr)
library(dplyr)
#library(magrittr)
library(conflicted)
#unloadNamespace("conflicted")
conflict_prefer("colnames", "base")
conflict_prefer("count", "dplyr")

# load all func.
source("/lulabdata3/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/util.R")
source("/lulabdata3/baopengfei/projects/motif-RBP-EDA/scripts/lulab/sncRNA_utils.R")
source("/lulabdata3/baopengfei/shared_utils/stats.R")

# read sncRNA/peak tab/mat 
dst <- "TCGA_small16" # TCGA_small16
dedup <- "_all"
diffpre <- "/lulabdata3/baopengfei/projects/motif-RBP-EDA"
pre="/lulabdata3/baopengfei/projects/WCHSU-FTC"
outpre="/lulabdata3/baopengfei/projects/WCHSU-FTC" 
dir.create(paste0(diffpre,"/output/",dst),recursive = T,showWarnings = F) # ,"/",disease.label,"vs",normal.label


## table1: read sample meta
sample.table0 <- read.csv( paste0(pre,"/exSeek-dev/data/",dst,"/sample_table.txt"),check.names = F,sep = "\t",header = T,encoding="UTF-8")
rownames(sample.table0) <- sample.table0$sample
table(sample.table0$group)
table(sample.table0$cell)


# replace new clinica from TCGAbiolinks
tmp <- data.table::fread("/lulabdata3/baopengfei/projects/tcga/data/TCGA-all_clinical/TCGAbiolinks_clinical.tsv",data.table = F,sep = '\t',header = T,stringsAsFactors = F,check.names = F)
colnames(tmp)[1] <- "patient"
tmp <- tmp[,c("patient","project","age","gender","race","ajcc_stage","ajcc_n","prior_treatment","cigarettes_per_day","years_smoked","days_to_death")]

sample.table0$age <- NULL
sample.table0$gender <- NULL
sample.table0$race <- NULL
sample.table0 <- dplyr::left_join(sample.table0, tmp, by="patient")
table(sample.table0$ajcc_stage)
sample.table0$tumor_stage <- sample.table0$ajcc_stage

tmp <- read.table(paste0("/lulabdata3/baopengfei/projects/tcga/data/TCGA-all_clinical/","UCEC","_clinical/","UCEC","_clinical.tsv"), quote = "",sep = "\t", header = T) # HNSC_gdc_sample_sheet.2024-01-09.tsv
tmp <- tmp[,c("submitter_id","figo_stage")]
colnames(tmp) <- c("patient","figo_stage")
sample.table0 <- dplyr::left_join(sample.table0, tmp, by="patient")
sample.table0$tumor_stage[sample.table0$cell %in% c("TCGA-UCEC_small")] <- sample.table0$figo_stage[sample.table0$cell %in% c("TCGA-UCEC_small")]

tmp <- read.table(paste0("/lulabdata3/baopengfei/projects/tcga/data/TCGA-all_clinical/","HNSC","_clinical/","HNSC","_clinical.tsv"), quote = "",sep = "\t", header = T) # HNSC_gdc_sample_sheet.2024-01-09.tsv
table(tmp$ajcc_clinical_stage)
tmp <- tmp[,c("submitter_id","ajcc_clinical_stage")]
colnames(tmp) <- c("patient","ajcc_clinical_stage")
sample.table0 <- dplyr::left_join(sample.table0, tmp, by="patient")
sample.table0$tumor_stage[sample.table0$cell %in% c("TCGA-HNSC_small")] <- sample.table0$ajcc_clinical_stage[sample.table0$cell %in% c("TCGA-HNSC_small")]

sample.table0$tumor_stage <- tolower(sample.table0$tumor_stage)
sample.table0$tumor_stage <- gsub("stage ia|stage ib|stage ic","stage i",sample.table0$tumor_stage, perl = T)
sample.table0$tumor_stage <- gsub("stage iia|stage iib|stage iic","stage ii",sample.table0$tumor_stage, perl = T)
sample.table0$tumor_stage <- gsub("stage iiic1|stage iiic2|stage iii1|stage iii2|stage iiia|stage iiib|stage iiic","stage iii",sample.table0$tumor_stage, perl = T)
sample.table0$tumor_stage <- gsub("stage iva|stage ivb|stage ivc","stage iv",sample.table0$tumor_stage, perl = T)
sample.table0$tumor_stage <- gsub("stage x","",sample.table0$tumor_stage, perl = T)
sample.table0$tumor_stage <- gsub("t2a|t2b|t2c|t2d","stage ii",sample.table0$tumor_stage, perl = T)
sample.table0$tumor_stage <- gsub("t3a|t3b|t3c|t3d","stage iii",sample.table0$tumor_stage, perl = T)
sample.table0$tumor_stage <- gsub("t4a|t4b|t4c|t4d|t4","stage iv",sample.table0$tumor_stage, perl = T)
sample.table0$tumor_stage <- gsub("n0","stage i",sample.table0$tumor_stage, perl = T)
sample.table0$tumor_stage <- gsub("n1","stage i",sample.table0$tumor_stage, perl = T)
sample.table0$tumor_stage <- gsub("n2a|n2b|n2c|n2d|n2","stage iii",sample.table0$tumor_stage, perl = T)
sample.table0$tumor_stage <- gsub("n3a|n3b|n3c|n3d|n3","stage iii",sample.table0$tumor_stage, perl = T)
sample.table0$tumor_stage <- gsub("n4a|n4b|n4c|n4d|n4","stage iv",sample.table0$tumor_stage, perl = T)
sample.table0$tumor_stage <- gsub("nx","",sample.table0$tumor_stage, perl = T)
sample.table0$group <- gsub("01A","Tumor",sample.table0$group)
sample.table0$group <- gsub("11A","Normal",sample.table0$group)
sample.table0$tumor_stage[sample.table0$group=="Normal"] <- "-"
table(sample.table0$tumor_stage,sample.table0$group)
##longer pattern to forward !!
table(sample.table0$tumor_stage)
sample.table0$tumor_stage[sample.table0$tumor_stage==""] <- "NA"

## table2: read count mat
#count <- read.table( paste0(outpre,"/output/",dst,"/call_peak",dedup,"/count_matrix/cfpeakCNN_b5_d50_p1.txt"),header = T,sep = "\t",stringsAsFactors = F,check.names = F)
count0 <- data.table::fread( paste0(outpre,"/output/",dst,"/call_peak",dedup,"/count_matrix/cfpeakCNN_b5_d50_p1.txt"), header = T,sep = "\t",stringsAsFactors = F,check.names = F, data.table = F)
#count0 <- data.table::fread( paste0(outpre,"/output/",dst,"/call_peak",dedup,"/count_matrix/cfpeakCNN_b5_d50_p1_1042_rmPatient.txt"), header = T,sep = "\t",stringsAsFactors = F,check.names = F, data.table = F)
#count0 <- data.table::fread( paste0(outpre,"/output/",dst,"/call_peak",dedup,"/count_matrix/cfpeakCNN_b5_d50_p1_1041_rmDate.txt"), header = T,sep = "\t",stringsAsFactors = F,check.names = F, data.table = F)
rownames(count0) <- count0$feature
count0$feature <- NULL


## define common vec/var
cmp.df <- data.frame(row.names=c("state","tumorigenesis","metastasis","relapse"),"col"=c("group","group1","group2","group3"),"pos.grp"=c("Tumor","Primary","Yes","Relapse"),"neg.grp"=c("Normal","Normal","No","Primary"))
sample.table0$group <- gsub("01A","Tumor",sample.table0$group, perl = T)
sample.table0$group <- gsub("11A","Normal",sample.table0$group, perl = T)
sample.table0$group2 <- ""
sample.table0$group2[sample.table0$group=="Tumor" ] <- "NA"
sample.table0$group2[sample.table0$tumor_stage %in% c("stage i","stage ia","stage ib","stage ii","stage iia","stage iib") & sample.table0$group=="Tumor" ] <- "No"
sample.table0$group2[sample.table0$tumor_stage %in% c("stage iii","stage iiia","stage iiib","stage iiic","stage iv") & sample.table0$group=="Tumor" ] <- "Yes"
sample.table0$group2[ sample.table0$group=="Normal" ] <- "-"
table(sample.table0$cell,sample.table0$group2)
table(sample.table0$cell,sample.table0$tumor_stage)

sample.table0$group1 <- ""
sample.table0$group1[ sample.table0$group=="Normal" ] <- "Normal"
sample.table0$group1[ sample.table0$group=="Tumor" & sample.table0$group2 %in% c("No","NA") ] <- "Primary"
sample.table0$group1[ sample.table0$group=="Tumor" & sample.table0$group2=="Yes" ] <- "Metastasis"
table(sample.table0$cell,sample.table0$group1)


## update sample tables by latest peak matrix
sample.table0 <- sample.table0[sample.table0$sample %in% colnames(count0),]
sample.table0$cellShort <- gsub("TCGA-|_small","",sample.table0$cell,perl = T)
sample.table0$cell <- sample.table0$cellShort

table(sample.table0$race)
sample.table0$race[sample.table0$race=="american indian or alaska native"] <- "alaska_indian"
sample.table0$race[sample.table0$race=="native hawaiian or other pacific islander"] <- "hawaiian"
sample.table0$race[sample.table0$race=="black or african american"] <- "black"
sample.table0$race[sample.table0$race=="not reported"] <- "NA"
sample.table0$race[sample.table0$race==""] <- "NA"

cell.table <- read.table("/lulabdata3/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/TCGA_small16/TCGA_superclass.txt",header = T,check.names = F,sep = "\t")
#table(sample.table0$cell %in% cell.table$Cancer)
cell.table <- cell.table[cell.table$Cancer %in% unique(sample.table0$cell),]
colnames(cell.table)[c(1,4)] <- c("cell","lineage")
cell.table$cellShort <- cell.table$cell
sample.table0$lineage <- cell.table$lineage[match(sample.table0$cell,cell.table$cell)]
unique(sample.table0$lineage )
#table(is.na(cell.table$lineage))
rownames(cell.table) <- cell.table$cell

#cellShort <- unique(sample.table0$cellShort)
#table(sample.table0$group,sample.table0$cellShort)
colorList <- list()
# cellShort <- c("STAD","COAD","LIHC","HNSC","LUAD","PRAD","BRCA",
#                "THCA","LUSC","ESCA","BLCA","CHOL","UCEC","KIRC","KIRP","KICH")
cellShort <- c("ESCA","STAD","COAD","LIHC","CHOL",   
               "HNSC",    
               "LUAD","LUSC",
               "BLCA","KIRC","KIRP","KICH",
               "THCA",
               "UCEC",
               "PRAD","BRCA"
)
#table(cellShort %in% cell.table$cell)
cellColor <- c(#"NC"="grey60","Normal"="grey60",
  # "Pancreatic disease"="grey30","Lung disease"="grey30","Colorectal disease"="grey30","Stomach disease"="grey30","Liver disease"="grey30","Prostate disease"="grey30","Immune disease"="grey30","Other disease"="grey30",
  "COAD"="orange3",#"CRC"="orange3",
  "STAD"="firebrick3", # "GC"="firebrick3",
  "ESCA"="#EFD493FF",
  # "PACA"="salmon","PAAD"="salmon","PDAC"="salmon",
  "CHOL"="#DBBC92FF","LIHC"="#FBBC92FF", #"Cholangiocarcinoma"="#DBBC92FF","HCC"="#FBBC92FF",
  "HNSC"="#FDF192FF",
  "THCA"="#97B792",
  
  "LUAD"="#7876B1FF","LUSC"="#9876B1FF",# "NSCLC"="#6876B1FF","LUCA"="#7876B1FF",
  "BRCA"="#EE4C97FF",
  "PRAD"="steelblue1", # "PRCA"="steelblue1",
  
  # "GBM"="seagreen1",
  # "LAML"="purple",
  
  "BLCA"="#91C5DE",
  # "RECA"="#B8862FFF",
  "UCEC"="#7A91B8", # "URCA"="#71853DFF",
  "KIRC"="seagreen4","KIRP"="seagreen1","KICH"="lightgreen"
  # "OVCA"="#F48194FF","ENCA"="#FB12FF",
  # "Sarcoma"="#7777B0FF","Lymphoma"="#9972B9FF","Melanoma"="#764F5CFF","Myeloma"="#6A7B71FF"
) #pal_nejm_adaptive()(10)
colorList[["cell"]] <- cellColor
colorList[["cellShort"]] <- cellColor

#unique(cell.table$lineage)
lineageShort <- c("Digest","Head&Neck","Thoracic","Urinary","Thyroid","Reproductive","Prostate","Breast")
lineageColor <- c(
  "Digest"="orange",
  "Head&Neck"="#FDF192FF",
  "Thoracic"="#7856B1FF",
  "Urinary"="steelblue3",
  "Thyroid"="salmon",
  "Reproductive"="yellow4",
  "Prostate"="steelblue1",
  "Breast"="#EE4C97FF"
  # "Brain","Eye",
  # "Skin&SoftTissue","Hematopoietic",
  # "Pleura","AdrenalGland","Thyroid"
)

cell.table$lineage <- factor(cell.table$lineage,levels = lineageShort)
cell.table$cellShort <- factor(cell.table$cellShort,levels = cellShort)
table(cell.table$lineage)
cell.table <- cell.table[order(cell.table$lineage,cell.table$cellShort,decreasing = F),]
colorList[["lineage"]] <- lineageColor

groupShort <- c("Normal","Tumor")
groupColor <- c("Normal"="grey60","Tumor"="firebrick2")
group1Short <- c("Normal","Primary")
group1Color <- c("Normal"="grey60","Primary"="firebrick2")
group2Short <- c("No","Yes")
group2Color <- c("No"="salmon","Yes"="firebrick4")
group3 <- c("Normal","Primary","Metastasis")
group3Short <- group3

colorList[["group"]] <- groupColor
colorList[["group1"]] <- c(group1Color,"Metastasis"="firebrick4")
colorList[["group2"]] <- group2Color
groupAllShort <- c("Normal"="grey60","Primary"="firebrick2","Metastasis"="firebrick4","No"="salmon","Yes"="firebrick4")

#DEG <- c("Tumorigenesis only","Metastasis only","Both significant","Not significant")
DEGfactor <- c("Not significant", "Tumorigenesis only", "Metastasis only", "Both significant")
DEGColor <- c("Tumorigenesis only"="#0072B5FF","Metastasis only"="#BC3C29FF","Both significant"="#E18727FF","Not significant"="#cccccc")
colorList[["DEG"]] <- DEGColor

for (i in c("gender","race","tumor_stage")){
  colorList[[i]] <- pal_npg_adaptive()( length(unique(sample.table0[[i]])) )
  names(colorList[[i]]) <- unique(sample.table0[[i]])
}

colorList[['response']] <- c("nonresponse"="grey","response"="firebrick")

colorList[['species']] <- RNA_colors
colorList[['precursor']] <- RNA_colors


## append ESTIMATE purity et al
tmp <- data.table::fread(file = "/lulabdata3/baopengfei/projects/tcga/ref/TCGAbiolinks_ESTIMATE.txt",sep = "\t",header = T, check.names = F, stringsAsFactors = F,data.table = F)
tmp <- tmp[,c(1,3:ncol(tmp))]
colnames(tmp)[1] <- "tissue_id"
sample.table0 <- dplyr::left_join(sample.table0,tmp,by="tissue_id")
rownames(sample.table0) <- sample.table0$sample
# table(is.na(sample.table0$CPE))
# cor.test(sample.table0$CPE,sample.table0$LUMP)
# cor.test(sample.table0$CPE,sample.table0$IHC)

## append purity metainfo
# priList <- list()
# for (i in paste0("/lulabdata3/wangtaiwei/project/RNA/SLE/output/ImmReg/pan_cancer/",cellShort,"_purity_scores.txt") ){
#   purity <- read.table(i,header = T,row.names = 1,check.names = F,stringsAsFactors = F,sep = "\t")
#   # table(duplicated(purity$sample))
#   purity$tissue_id <- substr(purity$sample,1,16)
#   priList[[i]] <- purity
# }
# priList <- as.data.frame(do.call(rbind,priList))
# data.table::fwrite(priList,paste0(pre,"/exSeek-dev/data/",dst,"/purity.txt"),quote = F, sep = '\t',row.names = F,col.names = T )
tmp <- data.table::fread(file = paste0(pre,"/exSeek-dev/data/",dst,"/purity.txt"),sep = "\t",header = T, check.names = F, stringsAsFactors = F,data.table = F)
tmp <- tmp[tmp$tissue_id %in% sample.table0$tissue_id,]
table(sample.table0$tissue_id %in% tmp$tissue_id)

tmp$sample <- NULL
table(duplicated(tmp$tissue_id))
tmp <- tmp[!duplicated(tmp$tissue_id),]
sample.table0 <- dplyr::left_join(sample.table0,tmp,by="tissue_id")
sample.table0$purity <- as.numeric(sample.table0$purity)
#table(is.na(sample.table0$purity))
#tmp2 <- sample.table0[is.na(sample.table0$purity),]
#data.table::fwrite(tmp2,paste0(pre,"/exSeek-dev/data/",dst,"/purity_NA.txt"),quote = F, sep = '\t',row.names = F,col.names = T ) # not matched mRNA-seq ?
#data.table::fwrite(x = sample.table0, paste0(pre,"/exSeek-dev/data/",dst,"/sample_table2.txt"), sep = "\t", quote = F, row.names = F, col.names = T )
rownames(sample.table0) <- sample.table0$sample
#some samples lost for NA purity
sample.table0 <- sample.table0[!is.na(sample.table0$purity),]


## get libsize
libSizePath <- paste0(outpre,"/output/",dst,"/call_peak",dedup,"/count_matrix/EM.libSize")
tmp <- read.table(libSizePath,header = F,sep = "\t",check.names = F,row.names = 1)
# tmp <- tmp[colnames(count),]
tmp <- tmp[sample.table0$sample,,drop=F]
externalLibSize <- tmp$V2 #[,"star_hg38_v38_dedup"]
#names(externalLibSize) <- colnames(count) # [2:ncol(count)]
names(externalLibSize) <- sample.table0$sample


# ## optional table7: read peak diff table
# #diff.sum.df <- read.table(paste0(diffpre,"/output/",dst,"/cfPeakCNN_smallDomain_diff.summary"), quote = "",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
# diff.sum.df <- data.table::fread(paste0(diffpre,"/output/",dst,"/cfPeakCNN_smallDomain_diff.summary"), quote = "",sep = "\t",header = T,check.names = F,stringsAsFactors = F,data.table=F)
# #diff.sum.df.sig2 <- read.table(paste0(diffpre,"/output/",dst,"/cfPeakCNN_smallDomain_diffSig.summary"), quote = "",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
# diff.sum.df.sig2 <- data.table::fread(paste0(diffpre,"/output/",dst,"/cfPeakCNN_smallDomain_diffSig.summary"), quote = "",sep = "\t",header = T,check.names = F,stringsAsFactors = F,data.table=F)
# #diff.sum.df.cellFreq <- read.table(paste0(diffpre,"/output/",dst,"/cell_diff.summary"), quote = "",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
# diff.sum.df.cellFreq <- data.table::fread(paste0(diffpre,"/output/",dst,"/cell_diff.summary"), quote = "",sep = "\t",header = T,check.names = F,stringsAsFactors = F,data.table=F)




# ## rm batch in TPM/CPM
logcpm0 <- data.table::fread( paste0(outpre,"/output/",dst,"/call_peak",dedup,"/count_matrix/cfpeakCNN_b5_d50_p1_CPM_rmPurity.txt"), sep = "\t",check.names = F,header = T,data.table = F,stringsAsFactors = F )
rownames(logcpm0) <- logcpm0$feature
logcpm0$feature <- NULL
summary(logcpm0[,2])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       1       1      28       1  470051 
logcpm0 <- log2(logcpm0+1)
dim(logcpm0) # 54354  3436





# # re-assign sncRNA types
# #not change 11RNA types! , only re-assign 4 DNA regions!
# #check sequential_assign_peak.R
# # append sncRNA types
# # sncRNA.table <- data.table::fread( paste0(outpre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN/gnAssign/snc_assign.txt.bak"), header = T,check.names = F,stringsAsFactors = F,sep = "\t")
# sncRNA.table <- data.table::fread( paste0(outpre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN/gnAssign/snc_assign.txt"), header = T,check.names = F,stringsAsFactors = F,sep = "\t")
# sncRNA.table <- sncRNA.table[order(as.integer(gsub("peak_","",sncRNA.table$name)),decreasing = F),]
# dim(sncRNA.table)
# #sncRNA.table[1:13,]  # filterd [10,200]
# rownames(logcpm0) <- gsub("_rev|_for","",rownames(logcpm0),perl = T)
# dim(logcpm0)
# rownames(logcpm0)[1:2]
# tmp <- strsplit(rownames(logcpm0),"|",fixed=T)
# tmp.all <- unlist(sapply(tmp,"[",1))
# tmp.region <- unlist(sapply(tmp,"[",2))
# #knitr::kable(table(tmp.region))
# #tmp.tx <- unlist(sapply(tmp,"[",3))
# tmp.name <- unlist(sapply(tmp,"[",4))
# tmp.tx <- unlist(sapply(tmp,"[",5))
# tmp.start <- unlist(sapply(tmp,"[",6))
# tmp.end <- unlist(sapply(tmp,"[",7))
# #tmp.name[1:12]
# #all(gsub("peak_","",tmp.name)==1:length(tmp.name))
# tmp.snc <- sncRNA.table$snc_type[match(tmp.name,sncRNA.table$name)]
# table(is.na(tmp.snc))
# tmp.snc[is.na(tmp.snc)] <- "exclude"
# table(tmp.region %in% rna)
# tmp.snc[tmp.region %in% rna] <- tmp.region[tmp.region %in% rna] # not change 11RNA types! , only re-assign 4 DNA regions!
# knitr::kable(table(tmp.snc))
# table(tmp.region,tmp.snc)
# # all(tmp.name==sncRNA.table$name)
# #sncRNA.table$name[1:3]
# length(tmp.snc)
# #paste0(c("1","2","3"),collapse = "|")
# tmp.merge <- lapply(seq_along(tmp.snc), function(i) {
#   paste(tmp.all[i], tmp.region[i], tmp.tx[i], tmp.name[i], tmp.tx[i], tmp.start[i], tmp.end[i], tmp.snc[i], sep = "|")
# })
# tmp.merge <- unlist(tmp.merge)
# tmp.merge[1:2]
# rownames(logcpm0) <- tmp.merge # paste0(tmp.all,"|",tmp.snc)
# rownames(logcpm0)[1:2]
# logcpm0 <- 2^logcpm0-1
# summary(logcpm0[,2])
# rownames(logcpm0)[1:2]
# data.table::fwrite(cbind("feature"=rownames(logcpm0),logcpm0),
#                    paste0(outpre,"/output/",dst,"/call_peak",dedup,"/count_matrix/cfpeakCNN_b5_d50_p1_CPM_rmPurity.txt"), sep = "\t",quote = F,row.names = F,col.names = T )
# 
# rownames(count0) <- rownames(logcpm0)
# rownames(count0)[1:2]
# data.table::fwrite(cbind("feature"=rownames(count0),count0),
#                    paste0(outpre,"/output/",dst,"/call_peak",dedup,"/count_matrix/cfpeakCNN_b5_d50_p1.txt"), sep = "\t",quote = F,row.names = F,col.names = T )
# 
# real.logcpm0 <- data.table::fread( paste0(outpre,"/output/",dst,"/call_peak",dedup,"/count_matrix/cfpeakCNN_b5_d50_p1_log2CPM_rmPurity.txt"), sep = "\t",check.names = F,header = T,data.table = F,stringsAsFactors = F )
# rownames(real.logcpm0) <- real.logcpm0$feature
# real.logcpm0$feature <- NULL
# rownames(real.logcpm0) <- rownames(logcpm0)
# #rownames(real.logcpm0)[1:2]
# data.table::fwrite(cbind("feature"=rownames(real.logcpm0),real.logcpm0),
#                    paste0(outpre,"/output/",dst,"/call_peak",dedup,"/count_matrix/cfpeakCNN_b5_d50_p1_log2CPM_rmPurity.txt"), sep = "\t",quote = F,row.names = F,col.names = T )
# 
# 
# #plot bar (real gn assign)
# dat <- as.data.frame(table(tmp.snc))
# colnames(dat) <- c("RNA","ratio")
# #str(dat)
# dat$RNA <- factor(dat$RNA, levels=c(rna2,dna))
# #table(dat$RNA)
# ggplot(dat,aes(x=RNA,y=ratio,fill=RNA))+
#   geom_bar(position = "stack",width = 0.8,size=0.3, stat = "sum",color="black")+
#   xlab("")+
#   ylab(paste0("Peak Assginment"))+
#   scale_y_log10()+
#   geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
#   scale_fill_manual(name="Precursor",values = c("pri_miRNA"="#5E576FFF","other"="grey50",RNA_colors)) + # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
#   theme_bw() +
#   theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
#         axis.title = element_text(size = 24,color ="black"),
#         axis.text = element_text(size= 20,color = "black"),
#         #panel.grid=element_blank(),
#         panel.grid.major.x=element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
#         #panel.grid.minor.y = element_blank(),
#         panel.border = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), #   , color = c(rep("steelblue",length(rna)),rep("firebrick",length(dna)))
#         legend.position = "right",#c(.25,.6),
#         legend.text = element_text(size= 16),
#         legend.title= element_text(size= 16))
# ggsave(paste0("./figure/",dst,"_","_assign_num.png"),width = 11,height = 7,units = "in",dpi = 320)







# filter test Samples after batch removal ------------------------------
## table6: read peak meta
# peak.table <- read.table(paste0(outpre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN/b5_d50_p1.summary"),header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "")
peak.table <- data.table::fread(paste0(outpre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN/b5_d50_p1.summary"),header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",data.table = F)
table(peak.table$gn_type,peak.table$tx_type)
rownames(peak.table) <- peak.table$name


## filter
testID <- read.table( paste0("/lulabdata3/baopengfei/projects/WCHSU-FTC/exSeek-dev/data/",dst,"/sample_ids_test.txt"), header = F,check.names = F,stringsAsFactors = F,sep = "\t")$V1
length(testID) # 2011  2068
dim(sample.table0) # 3507  3436
dim(sample.table0)
sample.table0 <- sample.table0[sample.table0$sample %in% testID,] # 2011
table(sample.table0$cell,sample.table0$group1)
table(is.na(sample.table0$purity))
logcpm0 <- logcpm0[,sample.table0$sample]
count0 <- count0[,sample.table0$sample]
dim(logcpm0)

