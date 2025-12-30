# TCGA long 16
setwd("/lulabdata3/baopengfei/projects/motif-RBP-EDA")
dst <- "TCGA_long16"
diffpre <- "/lulabdata3/baopengfei/projects/motif-RBP-EDA"

# library(tidyverse)
library(ggplot2)
library(ggrepel)
library(corrplot)
library(pheatmap)
library(RColorBrewer)
outSurfix="_rmPurity"


# prep/update var -------
mRNA_table <- read.table("output/TCGA_long16/mRNA_meta_intersect_rmDupName.txt",header = T,sep = '\t')
rownames(mRNA_table) <- mRNA_table$symbol
mRNA_table$RNA <- mRNA_table$group
table(mRNA_table$RNA)
# -       RBP RBP-RNase 
# 16766      1697       161 
RNA_colors <- c(RNA_colors,"-"="grey90","RBP"="orange","RBP-RNase"="firebrick")

logcpm0 <- data.table::fread("output/TCGA_long16/log2TPM_mRNA.txt",sep = '\t',check.names = F,header = T, data.table = F)
rownames(logcpm0) <- logcpm0$feature
logcpm0$feature <- NULL
#


# diff TCGA_long16 EDA (by cell/tissue type, grp1+2, 20240802) -------

# plot all PCA by var
sample.table <- sample.table0 #[(sample.table0$Threatment=="-") & (sample.table0$treatmentnaive %in% c("-","0")),] # 376
logcpm <- logcpm0[,sample.table$sample]
dim(logcpm)
all(sample.table$sample==colnames(logcpm)) # T !
#table(sample.table0$group)
table(sample.table0$cell,sample.table0$group1)

for(this.type in c("cell","group1")){ # "gender","race","tumor_stage",
  # this.type <- "cell"
  message(this.type)
  sample.table$Color <- sample.table[[this.type]]
  plotPCA(sample.table = sample.table, logcpm = logcpm, topVarNum = 10000, plotUMAP = "TRUE", outFile = paste0("./output/",dst,"/PCA/cpm.PCA.by",this.type,".png"), 
          plotGrpColor = colorList[[this.type]], plotDotSize = 2, plotDotAlpha = 0.99,legend.position = "right", figWid = 10, figHigh = 8) #  plotGrpFill = cellColor, plotGrpShape = c("21"=21),
}


# plot all PCA by each tissue/cell
for(cellID in cellShort){
  # cellID <- "HNSC"
  message(cellID)
  sample.table <- sample.table0[sample.table0$cell==cellID,]
  logcpm <- logcpm0[,sample.table$sample]
  for(this.type in c("cell","group1")){   # "gender","race","tumor_stage",
    message(this.type)
    sample.table$Color <- sample.table[[this.type]]
    plotPCA(sample.table = sample.table, logcpm = logcpm, topVarNum = 10000, plotUMAP = "TRUE", outFile = paste0("./output/",dst,"/PCA/",cellID,".cpm.PCA.by",this.type,".png"), 
            plotGrpColor = colorList[[this.type]], plotDotSize = 5, plotDotAlpha = 0.99,legend.position = "right", figWid = 10, figHigh = 8) #  plotGrpFill = cellColor, plotGrpShape = c("21"=21),
  }
}


## run diff
outDir = paste0(diffpre,"/output/",dst,"/")
dir.create(outDir,recursive = T,showWarnings = F)
#cmp.df
outSurfix="_rmPurity" # _rmPurity
#cellShort <- cellShort[cellShort!="EV"]
gL <- data.table::fread("/lulabdata3/baopengfei/projects/tcga/ref/TCGAbiolinks_geneInfoHT.txt")
gL <- gL[!duplicated(gL$ENSG),]
#rownames(gL) <- gL$ENSG
tmp <- as.numeric(gL$geneLength)
names(tmp) <- gL$ENSG
gL <- tmp
for(cellID in cellShort){
  # cellID <- "STAD"
  message(cellID)
  
  for(i in c("tumorigenesis","metastasis")){ # "tumorigenesis","metastasis"
    # i <- "tumorigenesis"
    message(i)
    
    pos.grp <- cmp.df[i,"pos.grp"]
    neg.grp <- cmp.df[i,"neg.grp"]
    col.grp <- cmp.df[i,"col"]
    
    method="edger_glmlrt"
    norm_method="TMM"
    
    
    sample.table <- sample.table0[sample.table0$cell==cellID,]
    
    longMat <- data.table::fread(file =  paste0("/lulabdata3/baopengfei/projects/tcga/data/","",cellID,"_long/","count.txt"), sep = "\t", header = T, check.names = F, stringsAsFactors = F,data.table = F)
    rownames(longMat) <- longMat$feature
    longMat$feature <- NULL
    longMat <- longMat[,!duplicated(colnames(longMat))] # duplicated col names !
    colnames(longMat) <- substr(colnames(longMat),1,16)
    longMat <- longMat[,colnames(longMat) %in% sample.table$tissue_id]
    longMat <- longMat[,sample.table$tissue_id]
    longMat <- longMat[!duplicated(unlist(sapply(strsplit(rownames(longMat),".",fixed = T),"[",1))),]
    rownames(longMat) <- unlist(sapply(strsplit(rownames(longMat),".",fixed = T),"[",1))
    longMat <- longMat[(rowSums(longMat)>0),]
    longMat <- na.omit(longMat)
    count <- longMat
    colnames(count) <- sample.table$sample[match(colnames(count),sample.table$tissue_id)]
    
    
    sample.table[["group"]] <- sample.table[[col.grp]] 
    #table(sample.table$group1)
    positive_samples <- sample.table[sample.table$group==pos.grp,"sample"]  # 
    negative_samples <- sample.table[sample.table$group==neg.grp,"sample"] #
    if(length(positive_samples)==0 | length(negative_samples)==0 ) {
      message(cellID,i,"skip for zero smp")
      next
    }
    samples <- c(positive_samples, negative_samples)
    group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
    sample.table <- sample.table[samples,]
    nrow(sample.table)==length(samples)
    sample.table$group <- factor(sample.table$group,levels = c(neg.grp,pos.grp))
    # sample.table$Month <- factor(sample.table$Month)
    
    mat <- count[,samples]
    useExternalLibSize <- "FALSE"
    #str(sample.table)
    
    form <- ~ group + purity # + Month + Operator
    res <- diff.v3.dcb(mat = mat, form = form,sample.table = sample.table,samples = samples, group = group,
                       method = method, norm_method = norm_method, filterType = "NULL", featureType = "gene", 
                       getCountGt1SmpFrac = TRUE, getPercExpByGrp = TRUE, useExternalLibSize = useExternalLibSize,geneLength = gL[rownames(mat)])

    # write.table(res[["diffTable"]],paste0(outDir,"/diff_",cellID,"_",pos.grp,"vs",neg.grp,"",outSurfix,".diff"),quote = F,sep = "\t",row.names = T,col.names = T) # 
    data.table::fwrite(cbind(feature=rownames(res[["diffTable"]]),res[["diffTable"]]),paste0(outDir,"/diff_",cellID,"_",pos.grp,"vs",neg.grp,"",outSurfix,".diff"),quote = F,sep = "\t",row.names = F,col.names = T) # 
    
  }
}


#











# diff TCGA_long16 summary (top/sig heatmap, volcano scatter) --------------
outDir = paste0(diffpre,"/output/",dst,"/")
dir.create(outDir,recursive = T,showWarnings = F)
pval.cutoff1 <- 0.01
padj.cutoff1 <- 0.01
pval.cutoff2 <- 0.01 
padj.cutoff2 <- 0.01 # set to 0.1 not change much
fc.cutoff <- 1

## opt1: by sigNum
#top.cutoff <- 1000 #200
diff.sum.list <- list()
for(cellID in cellShort){
  # cellID <- "STAD"
  message(cellID)
  for(i in c("tumorigenesis","metastasis")){ # "state","activity"
    # i <- "tumorigenesis"
    message(i)
    
    disease.label <- cmp.df[i,"pos.grp"]
    normal.label <- cmp.df[i,"neg.grp"]
    col.grp <- cmp.df[i,"col"]    
    x <- paste0("./output/",dst,"/diff_",cellID,"_",disease.label,"vs",normal.label,outSurfix,".diff") # cfPeakCNN_smallDomain_diff_
    # }
    if(!file.exists(x)){
      message("skip")
      next
    }
    tmp <- data.table::fread(x, header = T, sep="\t",check.names = F,stringsAsFactors = F, data.table = F) # _extLibSize
    tmp <- tmp[tmp$feature %in% mRNA_table$ENSG,]
    rownames(tmp) <- mRNA_table$symbol[match(tmp$feature,mRNA_table$ENSG)]
    tmp$feature <- NULL
    tmp$disease <- disease.label
    tmp$normal <- normal.label
    tmp$cell <- cellID  
    
    if(i=="tumorigenesis"){
      padj.cutoff <- padj.cutoff1
      pval.cutoff <- pval.cutoff1
    }else{
      padj.cutoff <- padj.cutoff2
      pval.cutoff <- pval.cutoff2
    }
    tmp2 <- filterDiffTabByP(diff = tmp, pvalue = pval.cutoff, padj = padj.cutoff, abs.log2FC = fc.cutoff)
    
    # tmp$diff <- i
    tmp$trend <- "insig"
    tmp$trend[rownames(tmp) %in% rownames(tmp2)[tmp2$log2FoldChange>0]] <- "up"
    tmp$trend[rownames(tmp) %in% rownames(tmp2)[tmp2$log2FoldChange<0]] <- "down"
    table(tmp$trend)
    diff.sum.list[[paste0(cellID,"|",disease.label)]] <- tmp
  }
}

diff.sum.df <- as.data.frame(do.call(cbind,unname(lapply(diff.sum.list, function(x) x[mRNA_table$symbol, ]))))
dim(diff.sum.df)
tmp <- diff.sum.df[,grepl("trend",colnames(diff.sum.df))]
hist(rowSums(tmp!="insig"))
tmp <- tmp[rowSums(tmp!="insig")>=1,]
#tmp
data.table::fwrite(cbind(feature=rownames(tmp),tmp), paste0("./output/",dst,"/diff",outSurfix,".diffSig"),quote = F,row.names = F,col.names = T,sep = '\t')
#only keep sig diff

tmp.l <- list()
for(i in 1:length(diff.sum.list)){
  # i <- 1
  diffTab <- diff.sum.list[[i]]
  diffTab <- as.data.frame(table(as.data.frame(diffTab$trend)))
  diffTab$name <- names(diff.sum.list)[i]
  tmp.l[[i]] <- diffTab
}
tmp.df <- as.data.frame(do.call(rbind,tmp.l))
table(tmp.df$Var1)
tmp.df <- tmp.df[tmp.df$Var1 %in% c("up","down"),]
tmp.df$cell <- unlist(sapply(strsplit(tmp.df$name,"|",fixed=T),"[",1))
tmp.df$diff <- unlist(sapply(strsplit(tmp.df$name,"|",fixed=T),"[",2))
table(tmp.df$diff)
tmp.df$diff2 <- tmp.df$diff
tmp.df$diff2 <- gsub("Primary","tumorigenesis",tmp.df$diff2)
tmp.df$diff2 <- gsub("Yes","metastasis",tmp.df$diff2)
tmp.df <- tmp.df[,c("cell","diff2","Var1","Freq")]
colnames(tmp.df) <- c("cell","diff","trend","number")
tmp.df <- reshape2::acast(data = tmp.df, formula = cell~diff+trend,value.var="number")
table(is.na(tmp.df))
tmp.df[is.na(tmp.df)] <- 0
#tmp.df <- tmp.df[,c("state_up","state_down","activity_up","activity_down")]
tmp.df <- tmp.df[,c(
  "tumorigenesis_up","tumorigenesis_down", "metastasis_up","metastasis_down"
)]
tmp.df <- tmp.df[cellShort[cellShort %in% rownames(tmp.df)],]
data.table::fwrite(cbind(cell=rownames(tmp.df),tmp.df), paste0("./output/",dst,"/diff",outSurfix,".diffSig.summary"),quote = F,row.names = F,col.names = T,sep = '\t')




## opt2: by topNum
outSurfix="_rmPurity" # _rmMonthOperator
for(cellID in cellShort){
  # cellID <- "STAD"
  message(cellID)
  
  diff <- list()
  peak.list <- list()
  seqs <- c(50,500) # 100,200,500,1000
  for (top in seqs){
    peak.list[[paste0("top",top)]] <- list()
  }
  
  for(i in c("tumorigenesis","metastasis")){
    # i <- "tumorigenesis"
    message(i)
    
    disease.label <- cmp.df[i,"pos.grp"]
    normal.label <- cmp.df[i,"neg.grp"]
    col.grp <- cmp.df[i,"col"]
    
    ## select peak by top
    x <- paste0("./output/",dst,"/diff_",cellID,"_",disease.label,"vs",normal.label,outSurfix,".diff")
    if(!file.exists(x)){
      message(x,"not exist, skip next")
      next
    }
    tmp <- data.table::fread(x, header = T, sep="\t",check.names = F,stringsAsFactors = F,data.table = F)
    tmp <- tmp[tmp$feature %in% mRNA_table$ENSG,]
    rownames(tmp) <- mRNA_table$symbol[match(tmp$feature,mRNA_table$ENSG)]
    tmp <- tmp[mRNA_table$symbol,]
    tmp <- na.omit(tmp)
    tmp$feature <- NULL
    diff[[disease.label]] <- tmp
    tmp <- tmp[order(as.numeric(tmp$pvalue),decreasing = F),]
    
    # op1: fixed cutoff
    #tmp.sig <- tmp[tmp$pvalue<0.00001,]
    
    # op2: fixed top
    for (i in seqs){
      peak.id <- rownames(tmp)[1:i]
      # print(length(peak.id))
      peak.list[[paste0("top",i)]][[disease.label]] <- peak.id
    }
  }
  
  peak.ids <- list()
  for (top in seqs){
    peak.ids[[paste0("top",top)]] <- unique(do.call("c",peak.list[[paste0("top",top)]]))
    print(length( peak.ids[[paste0("top",top)]] ))
  }
  
}
all(rownames(diff[[disease.label]])==rownames(logcpm0)) # T !!!
#



## heatmap 3grps
seqs <- c(500) # c(50,500,1000)
for(cellID in "LIHC"){
  # cellID <- "LIHC"
  message(cellID)
  
  for(top in paste0("top",seqs)){
    # top <- "top500"
    print(top)
    
    sample.table <- sample.table0[sample.table0$cellShort==cellID,] # 60
    #table(sample.table$group1)
    sample.table$group <- sample.table[["group1"]] # col.grp
    sample.table$group <- factor(sample.table$group,levels = c("Normal","Primary","Metastasis")) # "DC",
    sample.table <- sample.table[order(sample.table$group),]
    logcpm <- logcpm0[peak.ids[[top]],colnames(logcpm0) %in% sample.table$sample] # rownames(diff)[diff$pvalue<=0.01]
    logcpm <- logcpm[,sample.table$sample]
    
    table(sample.table$sample %in% colnames(logcpm))
    
    annotation_col <- sample.table[colnames(logcpm),]
    annotation_col <- annotation_col[,c("cell","group","tumor_stage","gender","race"),drop=F] # "group1","group2","age","sex","ExprDate","Month","Operator","MonthOperator"
    
    annotation_row <- data.frame( RNA = mRNA_table$RNA[match(rownames(logcpm),mRNA_table$symbol)], row.names = rownames(logcpm))
    annotation_row$trend <- ""
    
    #names(diff)
    a.pos <- "Primary"
    b.pos <- "Yes"
    annotation_row$trend[diff[[a.pos]][rownames(logcpm),]$log2FoldChange>0 & diff[[b.pos]][rownames(logcpm),]$log2FoldChange>0 ] <- "Both up"
    annotation_row$trend[diff[[a.pos]][rownames(logcpm),]$log2FoldChange>0 & diff[[b.pos]][rownames(logcpm),]$log2FoldChange<0 ] <- "up&down"
    annotation_row$trend[diff[[a.pos]][rownames(logcpm),]$log2FoldChange<0 & diff[[b.pos]][rownames(logcpm),]$log2FoldChange>0 ] <- "down&up"
    annotation_row$trend[diff[[a.pos]][rownames(logcpm),]$log2FoldChange<0 & diff[[b.pos]][rownames(logcpm),]$log2FoldChange<0 ] <- "Both down"
    annotation_row$trend <- factor(annotation_row$trend,levels = c("Both down","up&down","down&up","Both up"))
    table(annotation_row$trend)
    
    tmp <- as.data.frame(table(annotation_row$RNA))
    tmp$lab <- paste0(tmp$Var1," (n=",tmp$Freq,")")
    annotation_row$RNA.lab <-tmp$lab[match(annotation_row$RNA,tmp$Var1)]
    annotation_row$RNA <- factor(annotation_row$RNA, levels = c("-","RBP","RBP-RNase"))
    annotation_row <- annotation_row[order(annotation_row$RNA),]
    table(annotation_row$RNA)
    
    RNA_colors_tmp <- RNA_colors
    names(RNA_colors_tmp) <- tmp$lab[match(names(RNA_colors_tmp),tmp$Var1)]
    RNA_colors_tmp <- RNA_colors_tmp[!is.na(names(RNA_colors_tmp))]
    # annotation_row[annotation_row==""] <- "-"
    ann_colors <- list("Group"=groupAllShort[unique(as.character(sample.table$group1))], # "Group1"=c(group1Color,group2Color), "Group2"=c(group1Color,group2Color), 
                       "Tissue"=cellColor[unique(as.character(sample.table$cell))],
                       "Stage"=colorList[["tumor_stage"]],
                       "Race"=colorList[["race"]],
                       "Gender"=colorList[["gender"]],
                       "Regulation"=c("Both up"="salmon", "Both down"="steelblue4","up&down"="seagreen","down&up"="chocolate"), # , 
                       "Precursor"=RNA_colors_tmp
                       # "RNA"=c('-'='white','RNase'='grey30'))
                       ) # ,"salmon",seagreen,orange,chocolate

    annotation_row <- annotation_row[order(annotation_row$trend,annotation_row$RNA,decreasing = F),]
    annotation_row$RNA <- NULL
    logcpm <- logcpm[rownames(annotation_row),] #positive_samples
    colnames(annotation_row) <- c("Regulation","Precursor")
    colnames(annotation_col)[1:5] <- c("Tissue","Group","Stage","Gender","Race")
    
    
    dir.create(paste0("./output/",dst,"/heatmap/"),recursive = T,showWarnings = F)
    pheatmap::pheatmap(mat = logcpm, 
                       annotation_col = annotation_col, #data frame that specifies the annotations shown on left side of the heatmap.
                       annotation_row = annotation_row,
                       annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
                       border_color="NA",
                       # labels_col = coldata$age
                       scale = "row", # row,none
                       #labels_col = 3, labels_row = 6,
                       cluster_cols = F,cluster_rows = F,
                       # gaps_row = sum(annotation_row$trend=="up"), # invalid if cluster row/col
                       #cutree_cols = 2,cutree_rows = 3,
                       show_colnames=F, show_rownames=F,
                       fontsize = 12,
                       height = 14, width = 14,
                       color = rev( colorRampPalette( RColorBrewer::brewer.pal(n = 10, name = "RdBu"))(100) ),  # steelblue/dodgerblue4-firebrick2, viridis::viridis_pal()(100), #
                       #fontsize_row = 5,
                       filename = paste0("./output/",dst,"/heatmap/",cellID,"_",top,"_diff_",disease.label,"vs",normal.label,".cpm.heatmap2.png")
    )
  }
  
  
  
  ## heatmap 2grps
    
  ## plot diff volcano 
  # my_theme <-   theme(aspect.ratio=1,
  #                     plot.title = element_text(size = 28,color="black",hjust = 0.5),
  #                     axis.title = element_text(size = 28,color ="black"), 
  #                     axis.text = element_text(size= 24,color = "black"), #,face="bold
  #                     panel.grid.minor.y = element_blank(),
  #                     panel.grid.minor.x = element_blank(),
  #                     axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
  #                     axis.text.y = element_text( hjust = 0 ), # angle = 45,
  #                     panel.grid=element_blank(),
  #                     legend.position = "right",#c(0.5,0.3),
  #                     legend.text = element_text(size= 24,color = "black"),
  #                     legend.title= element_text(size= 26,color = "black"))
  # 
  for(i in c("tumorigenesis","metastasis")){
    # i <- "tumorigenesis"
    message(i)
    # i <- "tumorigenesis"
    disease.label <- cmp.df[i,"pos.grp"]
    normal.label <- cmp.df[i,"neg.grp"]
    col.grp <- cmp.df[i,"col"]
    
    diff.plot <- list()
    x <- paste0("./output/",dst,"/diff_",cellID,"_",disease.label,"vs",normal.label,outSurfix,".diff")
    if(!file.exists(x)){
      message(x,"not exist")
      next
    }
    tmp <- read.table(x,sep = "\t",header = T,check.names = F,stringsAsFactors = F)
    # rownames(tmp) <- tmp$feature
    tmp <- tmp[tmp$feature %in% mRNA_table$ENSG,]
    rownames(tmp) <- mRNA_table$symbol[match(tmp$feature,mRNA_table$ENSG)]
    tmp <- tmp[mRNA_table$symbol,]
    tmp$feature <- NULL
    
    # tmp <- tmp[diff_genes,]
    tmp <- na.omit(tmp)
    
    tmp$minus.log10.pvalue <- -log10(tmp$pvalue)
    tmp$log10.baseMean <- log10(tmp$baseMean)
    
    tmp$alpha <- maxmin.normalize.vec(tmp$minus.log10.pvalue)
    
    tmp$RNA <- mRNA_table$RNA[match(rownames(tmp),mRNA_table$symbol)] #unlist(sapply(strsplit(tmp$feature,"|",fixed=T),"[",8))
    tmp$RNA <- factor(tmp$RNA,levels=c("-","RBP","RBP-RNase"))

    tmp <- tmp[order(tmp$RNA,decreasing = F),]
    ## color by RNA
    diff.plot[[disease.label]] <- ggplot(tmp, aes(x=log10.baseMean,y=log2FoldChange, fill=RNA)) +  # minus.log10.pvalue
      geom_point(aes(alpha=0.8),shape=21,size=3,color="white") + 
      labs(title = paste0(disease.label," vs. ",normal.label)) +
      # ggraph::scale_fill_viridis(direction = -1) + 
      scale_fill_manual(values = RNA_colors)+
      geom_hline(yintercept = 0,linetype="dashed",color="salmon") + 
      theme_bw()+
      guides(fill = guide_legend(title = "-log10P" ), alpha = "none") +
      my_theme_volcano
    ggsave(plot = diff.plot[[disease.label]], paste0("./output/",dst,"/heatmap/",cellID,"_",top,"_diff_",disease.label,"vs",normal.label,"_diff_basemean2.png"), device = "png",dpi = 150, width=9, height=9) # "_",sample
  }
}
#








# diff TCGA_long16 stat: log2FC cor scatter --------------------
diffFCcor.list <- list()
outSurfix="_rmPurity"

for(cellID in cellShort){ # cellShort
  # cellID <- "LIHC"
  message(cellID)

  disease.label <- "Yes" # "HDA","Inactive"
  normal.label <- "No" # "Inactive","HC"
  diff1 <- read.table(paste0("./output/",dst,"/diff_",cellID,"_",disease.label,"vs",normal.label,outSurfix,".diff"), header = T, sep="\t",check.names = F)
  diff1 <- diff1[diff1$feature %in% mRNA_table$ENSG,]
  rownames(diff1) <- mRNA_table$symbol[match(diff1$feature,mRNA_table$ENSG)]
  diff1$feature <- NULL
  
  disease.label <- "Primary" # "HDA",""
  normal.label <- "Normal" # "Inactive","HC"
  diff2 <- read.table(paste0("./output/",dst,"/diff_",cellID,"_",disease.label,"vs",normal.label,outSurfix,".diff"), header = T, sep="\t",check.names = F)
  diff2 <- diff2[diff2$feature %in% mRNA_table$ENSG,]
  rownames(diff2) <- mRNA_table$symbol[match(diff2$feature,mRNA_table$ENSG)]
  diff2$feature <- NULL
  
  all(rownames(diff1)==rownames(diff2))
  diffFCcor.list[[cellID]] <- data.frame("HighvsLow"=diff1$log2FoldChange,"SLEvsHC"=diff2$log2FoldChange,
                                         "HighvsLow.padj"=diff1$padj,"SLEvsHC.padj"=diff2$padj,
                                         cell=cellID,"feature"=rownames(diff1))
}
all(rownames(diff1)==rownames(diff2)) # all T
diffFCcor.df <- as.data.frame(do.call(rbind,diffFCcor.list))
diffFCcor.df$cell <- factor(diffFCcor.df$cell, levels = cellShort)

#str(diffFCcor.df)
p <- 
  ggplot(data=diffFCcor.df,aes(x=HighvsLow, y=SLEvsHC))+
  geom_point(aes(color=cell), alpha=0.6, size=0.5 ) + # fill="grey",color="grey",
  scale_color_nejm_adaptive() +
  geom_smooth(method='lm', formula= y~x,color="grey30") +
  geom_vline(xintercept = 0, color="grey", alpha=0.8, linetype="dashed")+
  geom_hline(yintercept = 0, color="grey", alpha=0.8, linetype="dashed")+
  xlab("Metastasis log2FC") + ylab("Tumorigenesis log2FC") +
  theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 1,
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20,angle=90),
    axis.text.x = element_text(size = 20), #,hjust = 0,vjust = 1
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "right", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16)) #+
# facet_grid(.~cell)
dir.create("output/TCGA_long16/diff")
ggsave(plot = p, filename = paste0("./output/TCGA_long16/diff/","logFCcor_all",".png"),device = "png", width = 12,height = 10, dpi = 200 ) # , limitsize = FALSE

# for(cellID in cellShort){
#   # cellID <- "pDC"
#   message(cellID)
p <- ggplot(data=diffFCcor.df,aes(x=HighvsLow, y=SLEvsHC))+ # [diffFCcor.df$cell==cellID,]
  geom_point(alpha=0.5, color="white",aes(fill=cell), size=2, shape=21) + # fill="grey",color="grey",
  scale_fill_nejm_adaptive() +
  geom_smooth(method='lm', formula= y~x,color="grey30") +
  geom_vline(xintercept = 0, color="grey", alpha=0.8, linetype="dashed")+
  geom_hline(yintercept = 0, color="grey", alpha=0.8, linetype="dashed")+
  xlab("Metastasis log2FC") + ylab("Tumorigenesis log2FC") +
  theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 1,
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20,angle=90),
    axis.text.x = element_text(size = 20), #,hjust = 0,vjust = 1
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "none", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16)) +
  facet_wrap(facets = "cell", nrow = 4)
# ggsave(plot = p, filename = paste0("./figure/","logFCcor_byCell",".png"), device = "png", width = 20,height = 20, dpi = 200, limitsize = FALSE)
ggsave(plot = p, filename = paste0("./output/TCGA_long16/diff/","logFCcor_byCell",".png"),device = "png", width = 20,height = 20, dpi = 200, limitsize = FALSE ) # , limitsize = FALSE
#


# for visualization logFC lim =5,-5
dim(diffFCcor.df)
this.cells2 <- c("LIHC") #
res_sum_max5 = diffFCcor.df[diffFCcor.df$cell %in% this.cells2,]
colnames(res_sum_max5)[1:2] <- c("metastasis_logFC","tumorigenesis_logFC") # "HDAvsLDA" "LDAvsHC"  "cell"
log2FC.cutoff <- 10
res_sum_max5$tumorigenesis_logFC[res_sum_max5$tumorigenesis_logFC> log2FC.cutoff] = log2FC.cutoff
res_sum_max5$tumorigenesis_logFC[res_sum_max5$tumorigenesis_logFC< -log2FC.cutoff] = -log2FC.cutoff
res_sum_max5$metastasis_logFC[res_sum_max5$metastasis_logFC> log2FC.cutoff] = log2FC.cutoff
res_sum_max5$metastasis_logFC[res_sum_max5$metastasis_logFC< -log2FC.cutoff] = -log2FC.cutoff
padj.cutoff <- 0.05
res_sum_max5$DEGtype <- "NA"
res_sum_max5$DEGtype[(res_sum_max5$HighvsLow.padj > padj.cutoff) & (res_sum_max5$SLEvsHC.padj > padj.cutoff)] <- "Not significant"
res_sum_max5$DEGtype[(res_sum_max5$HighvsLow.padj <= padj.cutoff) & (res_sum_max5$SLEvsHC.padj > padj.cutoff)] <- "Metastasis only"
res_sum_max5$DEGtype[(res_sum_max5$HighvsLow.padj > padj.cutoff) & (res_sum_max5$SLEvsHC.padj <= padj.cutoff)] <- "Tumorigenesis only"
res_sum_max5$DEGtype[(res_sum_max5$HighvsLow.padj <= padj.cutoff) & (res_sum_max5$SLEvsHC.padj <= padj.cutoff)] <- "Both significant"
res_sum_max5$DEGtype <- factor(res_sum_max5$DEGtype,levels = c("Not significant","Metastasis only","Tumorigenesis only","Both significant"))
res_sum_max5$cell <- factor(res_sum_max5$cell, levels = this.cells2)
res_sum_max5$siggroup005 <- res_sum_max5$DEGtype
res_sum_max5 = res_sum_max5[order(res_sum_max5$siggroup005,decreasing = F),]
res_sum_max5 <- left_join(res_sum_max5, splitPeakFeatureID(unique(res_sum_max5$feature)), by="feature")
table(res_sum_max5$cell,res_sum_max5$DEGtype)
#label_group = c("Not significant","Tumorigenesis only","Metastasis only","Both significant")
# 9x 3 ver
p =  ggplot(res_sum_max5,aes(x=`tumorigenesis_logFC`,y=`metastasis_logFC`))+
  geom_point(size=0.2,aes(color=factor(siggroup005)))+
  theme_classic()+
  scale_color_manual(values = c("#cccccc","#BC3C29FF","#0072B5FF","#E18727FF"))+
  facet_wrap(~ cell, ncol = 3) +
  theme(aspect.ratio = 1, axis.text.x=element_text(colour="black",size=14),
        axis.text.y=element_text(colour="black",size=14),
        axis.title.x=element_text(colour="black",size=14),
        axis.title.y=element_text(colour="black",size=14),
        plot.title=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(colour="black",size=14),
        panel.background =element_rect(fill="transparent",colour="black",size=1),
        panel.grid = element_blank(),
        legend.position="right",
        legend.title=element_blank(),
        legend.text=element_text(colour="black",size=14))+
  guides(colour=guide_legend(override.aes=list(size=3)))+
  scale_x_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))+
  scale_y_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))+
  labs(x=paste0("Primary vs Normal (logFC)"),y=paste0("Metastasis vs Primary (logFC)"))
# p
# ggsave(plot = p, filename = paste0(diffpre,"/figure/","logFCcor_byCell_withClass",".png"),height=6,width=10,device = "png")
ggsave(plot = p, filename = paste0("./output/TCGA_long16/diff/","logFCcor_byCell_withClass",".png"),device = "png", height=6,width=10) # , limitsize = FALSE
#
# tmp <- diffFCcor.df[diffFCcor.df$cell=="Th2",] # 
# this.cells
#






# diff TCGA_long16 long diff candidates (merge diffTable and MetaTable) ---------------------
this.cells <- c("LIHC") # "pDC"
dir.create(paste0(diffpre, "/output/",dst,"/candidate/"),recursive = T,showWarnings = F)
dir.create(paste0("/lulabdata3/baopengfei/smb/TCGA_longRNA_candidate"),recursive = T,showWarnings = F)
# library(xlsx)
#install.packages("xlsx")
for(cellID in this.cells){
  # cellID <- "LIHC"
  message(cellID)
  
  for(i in c("tumorigenesis","metastasis")){
    # i <- "tumorigenesis"
    message(i)
    
    disease.label <- cmp.df[i,"pos.grp"]
    normal.label <- cmp.df[i,"neg.grp"]
    col.grp <- cmp.df[i,"col"]
    
    x <- paste0("./output/",dst,"/diff_",cellID,"_",disease.label,"vs",normal.label,outSurfix,".diff")
    if(!file.exists(x)){
      message("skip ",x)
      next
    } 
    tmp.diff <- read.table( x, header = T, sep = "\t", check.names = F) # row.names=1, 
    tmp.diff <- tmp.diff[tmp.diff$feature %in% mRNA_table$ENSG,]
    rownames(tmp.diff) <- mRNA_table$symbol[match(tmp.diff$feature,mRNA_table$ENSG)]
    tmp.diff$feature <- NULL
    
    tmp.diff2 <- filterDiffTabByP(diff = tmp.diff, pvalue = 0.01, padj = 0.01, abs.log2FC = 1) # 466/32139
    tmp.diff2 <- filterDiffTabByTopN(diff = tmp.diff, N = 5000, rankBy = "pvalue", countBy = "trend")
    # tmp.diff2 <- tmp.diff2[,grepl("logCPM|log2FoldChange|pvalue|padj|baseMean", colnames(tmp.diff2), perl = T)]
    
    colnames(tmp.diff2) <- paste0(cellID,"_", colnames(tmp.diff2)) # 304 --> 914 (PNK)
    
    consensus.diff <- tmp.diff2 # as.data.frame(cbind(plasma.diff2[g,],tcga.diff2[g,],tmp.diff2[g,]))
    consensus.diff$feature <- rownames(consensus.diff)
    consensus.diff <- consensus.diff[,c(ncol(consensus.diff),1:(ncol(consensus.diff)-1))]
    
    # consensus.diff <- merge.data.frame(splitPeakFeatureID(rownames(consensus.diff)),consensus.diff, by = "feature")
    
    # only same trend kept
    #consensus.diff <- consensus.diff[(consensus.diff$plasmaDT2_log2FoldChange * consensus.diff$plasmaMIR_log2FoldChange) > 0,] # 19 left
    
    # append peak meta
    consensus.diff <- consensus.diff[,grepl("logCPM|log2FoldChange|pvalue|padj|baseMean", colnames(consensus.diff), perl = T)] # PeakID
    # colnames(consensus.diff)[1] <- "name"
    consensus.diff$name <- rownames(consensus.diff)
    mRNA_table_tmp <- mRNA_table
    mRNA_table_tmp$name <- mRNA_table_tmp$symbol
    consensus.diff <- dplyr::left_join(consensus.diff,mRNA_table_tmp,by="name")
    
    consensus.diff <- consensus.diff[consensus.diff$RNA!="-",]
    consensus.diff[,c("RNA","symbol")] <- NULL
    # write.table(tmp, paste0(diffpre, "/output/",dst,"/candidate/",cellID,"_",disease.label,"vs",normal.label,"_candidate_diff_summary.txt") , sep = "\t",row.names = F,col.names = T, quote = F)
    write.table(consensus.diff, paste0(diffpre, "/output/",dst,"/candidate/",cellID,"_",disease.label,"vs",normal.label,"_candidate_diff_summary.txt") , sep = "\t",row.names = F,col.names = T)
    
    # # export simplified to multiple xlsx sheets
    # tmp <- tmp[,!grepl("site_overlap_base|site_overlap_ratio|gn_thickStart|gn_thickEnd|gn_itemRgb|entrez|biotype|description|MFE_structure|intraPair|peptide|coding_label",colnames(tmp),perl = T)]
    # xlsx::write.xlsx(tmp, file=paste0("/lulabdata3/baopengfei/smb/SLE_sncRNA_candidate/candidate_summary.xlsx"), sheetName=paste0(cellID,"_",disease.label,"vs",normal.label), append=TRUE, row.names=FALSE)
    # 
    # write.table(tmp, paste0("/lulabdata3/baopengfei/smb/SLE_sncRNA_candidate/",cellID,"_",disease.label,"vs",normal.label,"_candidate_diff_summary.txt") , sep = "\t",row.names = F,col.names = T, quote = F)
  }
}
#







# diff TCGA_long16 long stats boxplot among diff signature grps -----
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)

# 假设基因名为gene_name, 可以指定该基因的logTPM数据进行绘图
cand.mRNA <- list()
for(cellID in this.cells){
  # cellID <- "LIHC"
  message(cellID)
  
  for(i in c("tumorigenesis","metastasis")){
    # i <- "tumorigenesis"
    message(i)
    
    disease.label <- cmp.df[i,"pos.grp"]
    normal.label <- cmp.df[i,"neg.grp"]
    col.grp <- cmp.df[i,"col"]
    cand.mRNA[[paste0(cellID,"|",i)]] <- read.table(paste0(diffpre, "/output/",dst,"/candidate/",cellID,"_",disease.label,"vs",normal.label,"_candidate_diff_summary.txt") , sep = "\t",check.names = F, header = T)
  }
}
# names(cand.mRNA)
cand.mRNA <- as.data.frame(do.call(rbind,cand.mRNA))
cand.mRNA <- unique(cand.mRNA$name)
# cand.mRNA
table(goi %in% cand.mRNA)
goi.cand <- goi[goi %in% cand.mRNA]



# tname <- rnase_proteins # filter diff
# this.peak.ids <- unique(mRNA_table$ENSG[mRNA_table$symbol %in% tname])  # 替换为你要绘制的基因名
this.peak.ids <- c("IGF2BP1","PEG10","EEF1A2")  #goi.cand
# this.peak.ids <- rnase_genes
this.cells <- "LIHC" # c("LUAD","HNSC","LIHC","BRCA")
sample.table <- sample.table0 # 1042
sample.table$group <- sample.table$group1
table(sample.table$group)
sample.table$group <- factor(sample.table$group,levels = group3)
logcpm <- logcpm0 #logcpm0[,sample.table$sample]
dim(logcpm)

library(dplyr)
#gene.mat <- read.table("./output/rna-matrix/rna-matrix/CPM-TMM_matrix_gene.txt",sep = "\t",header = T, row.names = 1,check.names = F)
# tmp.feature <- 

gene.mat <- as.data.frame(t(logcpm[(rownames(logcpm)) %in% this.peak.ids,]))
#colnames(gene.mat) <- tmp.feature$[match(colnames(gene.mat),tmp.feature$feature)]
gene.mat.t <- gene.mat
gene.mat.t$sample <- rownames(gene.mat.t)
#gene.mat.t[1:4,]
gene.mat.t <- left_join(gene.mat.t,sample.table,by="sample")
dim(gene.mat.t)
## plot expr&pathway expr
dir.create(paste0(diffpre,"/output/",dst,"/candidate/"),showWarnings = F,recursive = T)
# dir.create("./output/candidate/")
for (g in this.peak.ids){
  # g <- "DROSHA"
  print(g)
  # this.peak.name <- mRNA_table$symbol[mRNA_table$ENSG %in% g][0]
  
  print(g %in% colnames(gene.mat.t))
  gene.mat.df <- reshape2::melt(gene.mat.t,id.var=c("group","cell"),measure.var=g)
  df <- gene.mat.df[gene.mat.df$variable==g,]
  
  table(df$group)
  df <- df[df$group %in% group3,]
  df$group <- factor(df$group,levels = group3)
  df <- df[df$cell %in% this.cells,]
  df$cell <- factor(df$cell,levels = this.cells)
  #str(df)
  df$value <- as.numeric(df$value)
  # str(df)
  set.seed(1)
  b <- runif(nrow(df), -0.2, 0.2)
  library(ggplot2)
  p <- ggplot(df,aes(x=group,y=value,fill=group))+  # ,label=label
    geom_boxplot(aes(x=(group),y=value,fill=group),color="grey30",alpha=0.6,outlier.size = 0) + # 
    geom_point(aes(x=(as.numeric(group)+b),y=value,fill=group),
               shape=21,size=2,color="white") +
    scale_fill_manual(values = alpha(groupAllShort,0.6)) +
    # scale_x_discrete(label=c("CT","CN","PT","PN","RT","RN"))+
    #ylim(c(0,1)) +
    labs(title = paste0(g,""))+
    #ylab("ssGSEA score")+
    ylab("log2(TPM+1)")+
    xlab("") +
    ggpubr::stat_compare_means(
      comparisons = list( c("Normal","Primary"),c("Primary","Metastasis"),c("Normal","Metastasis") ),
      label.x.npc = 0.5,hide.ns=F,size =8, vjust = 0.5, step.increase = 0.06,#paired = TRUE,
      aes(label = p.signif), # p.signif, group = Group,
      symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,1),symbols = c( "****","***", "**", "*", "ns")),
      method = "wilcox.test") +
    theme_minimal(base_size=12) +
    theme(#axis.ticks.x=element_blank(),
      # panel.border = element_rect(colour="white", fill = NA), #for panel.border # 
      # panel.grid.major.y = element_line(color = "grey"), #size= 1,
      strip.text.y = element_blank(),
      strip.text.x = element_text(family="sans",size=20), # serif
      axis.title.x = element_text(family="sans",size=20),
      axis.title.y = element_text(family="sans",size=20),
      axis.text.x = element_text(size = 20,family="sans",angle = 45,hjust=1,vjust=1),
      axis.text.y = element_text(size = 20,family="sans"),
      plot.title = element_text(size=20),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) +
    facet_wrap(facets = "cell", ncol = length(this.cells))
  # print(p)
  ggsave(plot = p,filename = paste0("./output/",dst,"/candidate/","TCGA_long","_mRNA_box_",g,".pdf"), device = "pdf", height = 6,width =8)
  # ggsave(plot = p,filename = paste0("./figure/",dst,"_candidate_box_",g,".",candidate.class,".pdf"), device = "pdf", height = 12,width =3)
}
#

