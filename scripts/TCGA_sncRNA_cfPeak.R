
# Paths and global options come from scripts/config.R; the sample table,
# count matrix, peak summary and colour palettes come from
# scripts/TCGA_sncRNA_util.R. Run this script from the repository root.
source("scripts/config.R")
source("scripts/TCGA_sncRNA_util.R")
setwd(PROJECT_DIR)

suppressPackageStartupMessages(library(bedtoolsr))
library(forcats)
library(plyr)      # Load plyr BEFORE dplyr to avoid conflicts
library(tidyr)
library(dplyr)
# library(tidyverse)  # REMOVE THIS LINE - you don't need it
library(ggplot2)
library(ggpubr)
library(ggsci)
library(hrbrthemes)
library(ComplexHeatmap)
library(circlize)
library(GetoptLong)
library(conflicted)

conflict_prefer("arrange", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("which", "base")

# func
getMethod <- function(method){
  if (grepl("expeak|cfpeak",method,perl = TRUE)){
    method.name <- "cfPeak"
  }else if  (grepl("localmax",method)){
    method.name <- "LocalMax"
  }else if  (grepl("clam",method)){
    method.name <- "CLAM"
  }else if  (grepl("clipper",method)){
    method.name <- "CLIPper"
  }else if (grepl("piranha",method)){
    method.name <- "Piranha"
  }
  return(method.name)
}



# peak number heatmap -----
# contingency table of genome region x feature type
row_order <- (paste0('chr',c(1:22,'X','Y','M')))
tmp <- peak.table[peak.table$gn %in% row_order,]

cross_table <- table(tmp$gn, tmp$gn_type)
cross_matrix <- as.matrix(cross_table)
cross_matrix <- log10(cross_matrix+1)

# reorder rows and columns in the specified order
cross_matrix <- cross_matrix[row_order, rev(c(c("pri_miRNA","piRNA","rRNA","lncRNA","mRNA","snoRNA","snRNA","srpRNA","tRNA",
                                             "tucpRNA","Y_RNA","ktRNA","Rfam"),
                                           dna))]

# colour palette
colors <- colorRampPalette(c("white", "firebrick"))(100)

pdf("./figure/tumor_chr_peak_num.pdf",width = 4,height = 4)
# heatmap without clustering so that the specified order is kept
heatmap(t(cross_matrix), 
        scale = "none",
        Rowv = NA,        # no row clustering
        Colv = NA,        # no column clustering
        col = colors,
        margins = c(8, 8),
        cexRow = 0.8, 
        cexCol = 0.8)
dev.off()



# diff TCGA_small16 EDA (EdgeR by cell/tissue type, grp1+2, 20240802) ----------------
setwd(PROJECT_DIR)

outDir <- paste0("./output/",dst)
all(sample.table0$sample==colnames(logcpm0))


# plot all PCA by var
sample.table <- sample.table0[(sample.table0$sample %in% testID),] #  & sample.table0$group!="Normal"
table(sample.table$group,sample.table$group1)
logcpm <- logcpm0[,sample.table$sample]
dim(logcpm)
all(sample.table$sample==colnames(logcpm))
#table(sample.table0$group)
knitr::kable(table(sample.table0$cell,sample.table0$group1))


# plot all PCA by RNA (Fig4)
peakType <- unlist(sapply(strsplit(rownames(logcpm0),"|",fixed = TRUE),"[",8)) #2
peakType <- gsub("_for|_rev","",peakType,perl = TRUE)
table(peakType)
for(this.type in c("cell")){ 
  # this.type <- "cell"
  message(this.type)
  sample.table$Color <- sample.table[[this.type]]
  for (j in unique(peakType)){
    # j <- "repeats"
    message(j)
    # table(peakType==j)
    plotPCA(sample.table = sample.table, logcpm = logcpm[peakType==j,], topVarNum = 10000, plotUMAP = "TRUE", outFile = paste0("./output/",dst,"/PCA/RNA/cpm.PCA.by",this.type,".by",j,".png"), 
            plotGrpColor = colorList[[this.type]], plotDotSize = 2, plotDotAlpha = 0.99,legend.position = "right", figWid = 7, figHigh = 5.6) #  plotGrpFill = cellColor, plotGrpShape = c("21"=21),
  }
}
#


# plot all PCA by each tissue/cell (normal + tumor)
for(cellID in cellShort){
  # cellID <- "LUAD"
  message(cellID)
  sample.table <- sample.table0[(sample.table0$cell==cellID) & (sample.table0$sample %in% testID),]
  logcpm <- logcpm0[,sample.table$sample]
  for(this.type in c("cell","gender","race","tumor_stage","group1")){ 
    message(this.type)
    sample.table$Color <- sample.table[[this.type]]
    #sample.table$Fill <- sample.table$cellShort
    #sample.table$Shape <- factor(21,levels = 1:21) #sample.table$group
    plotPCA(sample.table = sample.table, logcpm = logcpm, topVarNum = 10000, plotUMAP = "TRUE", outFile = paste0("./output/",dst,"/PCA/",cellID,".cpm.PCA.by",this.type,".png"), 
            plotGrpColor = colorList[[this.type]], plotDotSize = 5, plotDotAlpha = 0.99,legend.position = "right", figWid = 10, figHigh = 8) #  plotGrpFill = cellColor, plotGrpShape = c("21"=21),
  }
}

# plot all PCA by each tissue/cell (only tumor, prep res for parVar)
for(cellID in cellShort){
  # cellID <- "LUAD"
  message(cellID)
  sample.table <- sample.table0[(sample.table0$cell==cellID) & (sample.table0$sample %in% testID)  & (sample.table0$group!="Normal"),]
  logcpm <- logcpm0[,sample.table$sample]
  for(this.type in c("cell","gender","race","tumor_stage","group1")){ 
    message(this.type)
    sample.table$Color <- sample.table[[this.type]]
    plotPCA(sample.table = sample.table, logcpm = logcpm, topVarNum = 10000, plotUMAP = "TRUE", outFile = paste0("./output/",dst,"/PCA/",cellID,".cpm.PCA.by",this.type,".tumor.png"), 
            plotGrpColor = colorList[[this.type]], plotDotSize = 5, plotDotAlpha = 0.99,legend.position = "right", figWid = 10, figHigh = 8) #  plotGrpFill = cellColor, plotGrpShape = c("21"=21),
  }
}


# plot all PCA by RNA by each tissue/cell
for(cellID in cellShort){
  # cellID <- "HNSC"
  message(cellID)
  sample.table <- sample.table0[(sample.table0$cell==cellID) & (sample.table0$sample %in% testID),]
  logcpm <- logcpm0[,sample.table$sample]
  # table(sample.table0$group1)
  for(this.type in c("group1")){ 
    # this.type <- "group1"
    message(this.type)
    sample.table$Color <- sample.table[[this.type]]
    for (j in unique(peakType)){
      # j <- "enhancer"
      message(j)
      plotPCA(sample.table = sample.table, logcpm = logcpm[peakType==j,], topVarNum = 10000, plotUMAP = "TRUE", outFile = paste0("./output/",dst,"/PCA/cell_RNA/",cellID,".cpm.PCA.by",this.type,".by",j,".png"), 
              plotGrpColor = colorList[[this.type]], plotDotSize = 2, plotDotAlpha = 0.99,legend.position = "right", figWid = 10, figHigh = 8) #  plotGrpFill = cellColor, plotGrpShape = c("21"=21),
    }
  }
}
#

## run diff
outDir = paste0(diffpre,"/output/",dst,"/")
dir.create(outDir,recursive = TRUE,showWarnings = FALSE)
#cmp.df
outSurfix="_rmPurity" # _rmMonthOperato
#cellShort <- cellShort[cellShort!="EV"]
for(cellID in cellShort){
  # cellID <- "STAD"
  message(cellID)
  
  
  for(i in c("tumorigenesis","metastasis")){ # "tumorigenesis","metastasis"
    # i <- "metastasis"
    message(i)
    
    pos.grp <- cmp.df[i,"pos.grp"]
    neg.grp <- cmp.df[i,"neg.grp"]
    col.grp <- cmp.df[i,"col"]
    
    method="edger_glmlrt"
    norm_method="TMM"
    # outDir = paste0(diffpre,"/output/",dst,"/")
    # outSurfix="_rmMonthOperator"
    
    sample.table <- sample.table0[ (sample.table0$cellShort==cellID & !(is.na(sample.table0$purity)) & (sample.table0$sample %in% testID)) ,] # 1042
    sample.table[["group"]] <- sample.table[[col.grp]] 
    #sample.table$purityN <- 1-sample.table$purity
    
    positive_samples <- sample.table[sample.table$group==pos.grp,"sample"]  # 
    negative_samples <- sample.table[sample.table$group==neg.grp,"sample"] #
    if(length(positive_samples)==0 | length(negative_samples)==0 ) {
      message(cellID,i,"skip for zero smp")
      next
    }
    samples <- c(positive_samples, negative_samples)
    group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
    sample.table <- sample.table[samples,]
    sample.table$group <- factor(sample.table$group,levels = c(neg.grp,pos.grp))
    # sample.table$Month <- factor(sample.table$Month)
    
    mat <- count0[,samples]
    useExternalLibSize <- externalLibSize[samples]
    #str(sample.table)
    
    form <- ~ group + purity # + Month + Operator
    res <- diff.v3.dcb(mat = mat, form = form,sample.table = sample.table,samples = samples, group = group,
                       method = method, norm_method = norm_method, filterType = "NULL", featureType = "domain", 
                       getCountGt1SmpFrac = TRUE, getPercExpByGrp = TRUE, useExternalLibSize = useExternalLibSize)

    # write.table(res[["diffTable"]],paste0(outDir,"/diff_",cellID,"_",pos.grp,"vs",neg.grp,"",outSurfix,".diff"),quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE) # 
    data.table::fwrite(cbind(feature=rownames(res[["diffTable"]]),res[["diffTable"]]),paste0(outDir,"/diff_",cellID,"_",pos.grp,"vs",neg.grp,"",outSurfix,".diff"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE) # 
    
  }
}
#



# basic stat (QC): read map ratio, number ------------
## tx map ratio
# the QC module below reads and writes inside the cfPeak workspace
pre    <- DATA_DIR
outpre <- DATA_DIR
setwd(outpre)

# small/SE map ratio (main QC module)
read.SE.map <- function(x){
  print(x)
  sample <- unlist(sapply(strsplit(x,"/",fixed = TRUE),"[",4))
  rna <- gsub(".log","",basename(x))  
  map0 <- readLines(x)
  in.num <- as.numeric(gsub(" reads; of these:","",map0[grepl("reads; of these:",map0)]))
  not.map.num <- as.numeric(unlist(sapply(strsplit(map0[grepl("aligned 0 times",map0)]," (", fixed = TRUE), "[" ,1)))
  uniq.map.num <- as.numeric(unlist(sapply(strsplit(map0[grepl("aligned exactly 1 time",map0)]," (", fixed = TRUE), "[" ,1)))
  multi.map.num <- as.numeric(unlist(sapply(strsplit(map0[grepl("aligned >1 time",map0)]," (", fixed = TRUE), "[" ,1)))
  map <- as.numeric(gsub("% overall alignment rate","",map0[6]))
  map <- data.frame(sample=sample,type=rna,overall.align.ratio=map,in.num=in.num, uniq.map.num=uniq.map.num, multi.map.num=multi.map.num,not.map.num=not.map.num)
  return(map)
}


map.df <- Sys.glob( paste0("output/",dst,"/tbam/*/log/*.log") )
length(map.df)
#map.df[1:3]
this.smps <- unlist(sapply(strsplit(map.df,"/",fixed=TRUE),"[",4)) # map.df[1:3]
head(this.smps)
sample.table <- sample.table0[sample.table0$sample %in% testID,]
map.df <- map.df[this.smps %in% sample.table$sample]
length(map.df)
map.df2 <- parallel::mclapply(map.df,read.SE.map,mc.cores = 10)
map.df2 <- do.call(rbind,map.df2)
map.df2$overall.align.ratio <- map.df2$overall.align.ratio/100


library(tidyverse)
dat <- dplyr::as_tibble(map.df2)
dat <- dat %>% 
  dplyr::group_by(sample) %>% 
  dplyr::mutate(map.num=uniq.map.num+multi.map.num, 
                all.num=max(in.num),
                unmap.num=min(not.map.num),
                map.ratio=map.num/all.num,
                unmap.ratio=unmap.num/all.num) %>% 
  dplyr::select(sample,type,all.num,map.ratio,unmap.ratio)
dat.wide <- dat %>% tidyr::pivot_wider(names_from = "type", values_from = "map.ratio")
colnames(dat.wide)[1:3] <- c("sample","read_num","unmap")
colnames(dat.wide) <- gsub(".","_",colnames(dat.wide),fixed = TRUE)
colnames(dat.wide) <- gsub("miRNA","pri_miRNA",colnames(dat.wide),fixed = TRUE)
dat.wide <-  dat.wide[,c("sample","read_num","unmap","rRNA","spikein","univec","pri_miRNA","lncRNA", "mRNA", "piRNA", 'snoRNA', "snRNA", 'srpRNA', 'tRNA', 'tucpRNA', 'Y_RNA', "intron_for", "intron_rev", "promoter_for", "promoter_rev", "enhancer_for", "enhancer_rev", "repeats_for", "repeats_rev")]
rowSums(dat.wide[,3:ncol(dat.wide)])
dat <- dat.wide %>% tidyr::pivot_longer(cols = 3:ncol(dat.wide), names_to = "type", values_to = "ratio")

dat$group <- sample.table0$group[match(dat$sample,sample.table0$sample)]
dat$cell <- sample.table0$cellShort[match(dat$sample,sample.table0$sample)]

tmp <- data.frame(cell=dat$cell,group=dat$group)
tmp <- tmp[!(duplicated(tmp)),]
tmp$cell <- factor(tmp$cell ) # levels = c("tongue")
tmp$group <- factor(tmp$group) # , levels = c("N","Y")
tmp <- tmp[order(tmp$cell,tmp$group),]
tmp$lib <- paste(tmp$cell,tmp$group)

dat$lib <- paste(dat$cell,dat$group)
# dat$lib <- paste(dat$cell,dat$group)
dat$lib <- factor(dat$lib, levels = unique(tmp$lib))
table(dat$lib)
dat$sample.group <- dat$lib
table(dat$lib)
dat$type <- gsub("miRNA","pri_miRNA",dat$type)
dat$type <- gsub("pri_pri_miRNA","pri_miRNA",dat$type)
dat$RNA <- gsub("\\.for|\\.rev|_for|_rev|\\.ratio","",dat$type,perl=TRUE)
dat$RNA <- gsub("spikein_small|spikein_long","spikein",dat$RNA,perl=TRUE)
other <- c("spikein","univec")
dat$RNA <- factor(dat$RNA,levels = c(other,rna,dna,"unmap"))
dat$source <- ifelse(as.character(dat$RNA) %in% rna,"RNA","DNA")
dat$source[as.character(dat$RNA) %in% other] <- "other"
dat$source[as.character(dat$RNA) %in% "unmap"] <- "unmap"
dat$source <- factor(dat$source,levels = c("other","RNA","DNA","unmap"))
summary(dat$RNA)

### plot stack bar
library(ggplot2)

dat <- dat[order(dat$cell,dat$group,dat$sample,decreasing = FALSE),] 
#dat[1:3,]
dat$sample <- factor(dat$sample,levels = unique(dat$sample))

#with unmap
dat.sum <- dat %>% 
  # dplyr::filter(RNA!="unmap") %>% 
  dplyr::group_by(RNA,sample.group) %>% 
  dplyr::summarize(mean.ratio=mean(ratio,trim=0.1,na.rm=TRUE)) # %>% 
#str(dat.sum)
table(dat.sum$sample.group )

combinations <- outer(cellShort, c("Normal","Tumor"), FUN = function(x, y) paste0(x," ",y))
combinations <- as.vector(t(combinations))
dat.sum$sample.group <- factor(dat.sum$sample.group, levels = combinations) 
p1 <- ggplot(dat.sum,aes(x=sample.group,y=mean.ratio,fill=RNA))+ 
  geom_bar(position = "fill",width = 0.8,size=0.3, stat = "identity", color="black") +
  #geom_errorbar()+
  xlab("")+
  ylab(paste0("Reads mapped ratio"))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+  
  scale_fill_manual(name="Species",values = c("unmap"="grey",RNA_colors) )+ # colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(length(c(other,rna,dna,"unmap")))
  theme_bw() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.text = element_text(size= 20,color = "black"),
        axis.text.x = element_text(size= 20,color = "black",angle = 90,vjust = 0.5,hjust = 1),
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        legend.position = "right",#c(.25,.6),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))
ggsave(plot = p1, filename =paste0("./figure/",dst,"_group_small_map_ratio_withUnmap.png"), width = 14,height = 8)


#rm unmap
dat.sum <- dat %>% 
  dplyr::filter(RNA!="unmap") %>% 
  dplyr::group_by(RNA,sample.group) %>% 
  dplyr::summarize(mean.ratio=mean(ratio,trim=0.1,na.rm=TRUE)) # %>% 
p1 <- ggplot(dat.sum,aes(x=sample.group,y=mean.ratio,fill=RNA))+ 
  geom_bar(position = "fill",width = 0.8,size=0.3, stat = "identity", color="black") +
  #geom_errorbar()+
  xlab("")+
  ylab(paste0("Reads mapped ratio"))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+  
  scale_fill_manual(name="Species",values = RNA_colors )+ # 
  theme_bw() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.text = element_text(size= 20,color = "black"),
        axis.text.x = element_text(size= 20,color = "black",angle = 90,vjust = 0.5,hjust = 1),
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(), 
        legend.position = "right",#c(.25,.6),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)
  )



#table(dat$RNA)
dat.sum <- dat %>% 
  dplyr::filter(RNA!="unmap") %>% 
  dplyr::group_by(RNA,sample.group) %>% 
  dplyr::summarize(mean.number=mean(read_num,trim=0.1,na.rm=TRUE)) %>% 
  dplyr::group_by(sample.group) %>% 
  dplyr::summarize(sum.number=sum(mean.number/1000000,na.rm=TRUE))
# dplyr::group_by(sample.group) %>% 
# dplyr::mutate(sum.ratio=sum(mean.ratio),norm.mean.ratio=mean.ratio/sum.ratio) # 
p2 <- ggplot(dat.sum,aes(x=sample.group,y=sum.number ))+  # fill=RNA
  geom_bar(width = 0.8,size=0.3, stat = "identity", color="grey10",fill="grey50") + # position = "fill",
  #geom_errorbar()+
  xlab("")+
  ylab(paste0("Clean reads number (M)"))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  scale_fill_manual(name="Species",values = RNA_colors )+ # colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(length(c(other,rna,dna,"unmap")))
  theme_bw() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.text = element_text(size= 20,color = "black"),
        axis.text.x = element_blank(), #element_text(size= 20,color = "black",angle = 90,vjust = 0.5,hjust = 1), # element_blank(), #
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = c(rep("grey50",length(other)),rep("steelblue",length(rna)),rep("firebrick",length(dna)))), #  
        legend.position = "right",#c(.25,.6),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16),
        plot.margin = unit(x=c(10,0,-20,0),units="pt") ) # top, right, bottom, left

p <- ggpubr::ggarrange(plotlist = list(p2,p1),ncol = 1,nrow=2,common.legend = TRUE,align = "v",heights=c(0.16,0.2),legend = "right")
ggsave(plot = p, filename =  paste0("./figure/",dst,"_group_small_map_merge_noUnmap.png"), width = 14,height = 11)
#








# basic stat (QC): sample tx/peak stats num/ratio stackbar -------------------------------------------------------
## read ref
ref <- data.table::fread(file.path(DATA_DIR, "exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt"),header = TRUE,sep = "\t", stringsAsFactors = FALSE)
ref <- ref[!is.na(ref$tx.length),]
colnames(ref)

## read sample table
### plot detected tx number
readDomain <- function(sample){
  x <- paste0(outpre,"/output/",dst,"/call_peak_all/cfpeak_by_sample/b5_d50_p1/",sample,".bed") # not use CNN 
  if(!file.exists(x)){
    x <- paste0(outpre,"/output/",dst,"/call_peak_dedup/cfpeak_by_sample/b5_d50_p1/",sample,".bed") # # not use CNN 
  }
  
  peak <- data.table::fread(x,data.table = FALSE,sep = '\t',header = FALSE,stringsAsFactors = FALSE)
  peak <- as.data.frame(peak)
  peak <- peak[,1:6]
  colnames(peak) <- c("chr","start","end","name","score","strand") # cfpeak
  peak$sample <- sample
  return(peak)
}


#sample.table <- sample.table0[sample.table0$sample %in% testID,]
#tmp <- lapply(sample.table0$sample,readDomain)
tmp0 <- parallel::mclapply(sample.table0$sample,readDomain,mc.cores = 10)
x <- do.call("rbind",tmp0)

x$tx.length <- ref$tx.length[match(x$chr,ref$transcript_id)]
x$RNA <- ref$transcript_type[match(x$chr,ref$transcript_id)]
x$RNA <- gsub("\\.for|\\.rev|_for|_rev","",x$RNA,perl = TRUE)
x$RNA <- factor(x$RNA,levels = c(rna,dna))
#table(x$RNA)

x <- as.data.frame(x)
x$length <- as.numeric(x$end-x$start) #/x$tx.length # /x$length
x$length.ratio <- x$length/x$tx.length
x[1:3,]
dat <- x
rm(list = "x")
dat$group <- sample.table0$group[match(dat$sample,sample.table0$sample)]
dat$cell <- sample.table0$cellShort[match(dat$sample,sample.table0$sample)]

tmp <- data.frame(cell=sample.table0$cellShort,group=sample.table0$group)
tmp <- tmp[!(duplicated(tmp)),]
tmp$cell <- factor(tmp$cell, levels = cellShort) # 
tmp$group <- factor(tmp$group, levels = groupShort) # 
tmp <- tmp[order(tmp$cell,tmp$group),]
tmp$lib <- paste(tmp$cell,tmp$group)

dat$lib <- paste(dat$cell,dat$group)
# dat$lib <- paste(dat$cell,dat$group)
dat$lib <- factor(dat$lib, levels = unique(tmp$lib))
table(dat$lib)
dat$sample.group <- dat$lib
table(dat$lib)
dat$group.lib <- dat$lib

my_theme <- theme_bw() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.text = element_text(size= 20,color = "black"),
        axis.text.x = element_text(size= 16,color = "black",angle = 90,vjust = 0.5,hjust = 1),
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = c(rep("grey50",length(other)),rep("steelblue",length(rna)),rep("firebrick",length(dna)))), #  
        legend.position = "right",#c(.25,.6),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))



### plot peak,tx num & ratio stackbar
dim(dat)
dat[1:3,]
library(tidyverse)
df <- dplyr::as_tibble(dat) %>% 
  dplyr::mutate(total.distinct.tx.num=n_distinct(chr)) %>% 
  dplyr::group_by(sample) %>% 
  dplyr::mutate(sample.tx.num=n_distinct(chr), sample.peak.num=n_distinct(name)) %>% 
  dplyr::group_by(sample,RNA) %>% 
  dplyr::mutate(sample.RNA.tx.num=n_distinct(chr), sample.RNA.peak.num=n_distinct(name), 
                RNA.tx.ratio=sample.RNA.tx.num/sample.tx.num, RNA.peak.ratio=sample.RNA.peak.num/sample.peak.num,
                total.number=dplyr::n()) %>% 
  dplyr::select(chr,name,sample,group,lib,group.lib,RNA,sample.tx.num,sample.peak.num,sample.RNA.tx.num,sample.RNA.peak.num,RNA.tx.ratio,RNA.peak.ratio,total.number,total.distinct.tx.num) # no name & chr ? (add name & chr at 2403)

dim(df)
df[1:3,]
#table(duplicated(df))
df <- df[!duplicated(paste0(df$name,df$sample)),]


#### peak chr location
ref2 <- data.table::fread(file.path(DATA_DIR, "exSeek-dev/genome/hg38/transcript_table/all_newTxID.txt"),header = TRUE,sep = "\t", stringsAsFactors = FALSE)
ref2 <- ref2[!is.na(ref2$tx.length),]
colnames(ref2)
ref2[1:3,]
ref2 <- ref2[ref2$transcript_id %in% df$chr,]
table(df.tmp$chr %in% ref2$transcript_id)
table(ref2$chrom)
ref2$gn.chr <- ref2$chrom

gene.ref2 <- data.table::fread(file.path(REF_DIR, "geneset/encodev27_gene_v2.txt"),header = TRUE,sep = "\t", stringsAsFactors = FALSE)
table(df.tmp$chr %in% gene.ref2$name)
ref2$gn.chr[ref2$transcript_type=='rRNA'] <- gene.ref2$chr[match(ref2$gene_name[ref2$transcript_type=='rRNA'], gene.ref2$name)]
chr.seqs <- paste0("chr",c(1:22,"X","Y","M"))
table(ref2$gn.chr %in%  chr.seqs)
table(ref2$gn.chr)

df.tmp <- df
df.tmp[1:4,]
df.tmp$gn.chr <- ref2$gn.chr[match(df.tmp$chr,ref2$transcript_id)]
table(df.tmp$gn.chr)
df.tmp <- df.tmp[df.tmp$gn.chr %in% chr.seqs,]
#

#### heatmap of peak stat by gn.chr or RNA
df.sum <- df.tmp %>% 
  dplyr::group_by(RNA,group.lib,gn.chr) %>% 
  dplyr::summarize(mean.ratio=mean(RNA.peak.ratio,trim=0.1,na.rm=TRUE),mean.total.number=mean(total.number,na.rm=TRUE,trim=0.1)) %>%  # average biol rep within same group.lib grp
  dplyr::group_by(group.lib,gn.chr) %>% 
  dplyr::mutate(sum.ratio=sum(mean.ratio),norm.mean.ratio=mean.ratio/sum.ratio)  #%>% 

table(df.sum$RNA)
df.sum <- df.sum[grepl("Tumor",df.sum$group.lib),] # Tumor, Normal

#opt1: by gn.chr
df.sum2 <- df.sum %>% 
  # dplyr::filter(RNA=="tRNA") %>% 
  dplyr::select(group.lib, RNA, gn.chr, mean.total.number) %>% 
  tidyr::pivot_wider(names_from = "gn.chr",id_cols = "group.lib",
                     values_from = "mean.total.number",
                     values_fn = sum
                     )
for(i in chr.seqs){
  if(!(i %in% colnames(df.sum2))){
    message("add ",i)
    df.sum2[[i]] <- 0
  }
}
df.sum2 <- as.data.frame(df.sum2)
rownames(df.sum2) <- df.sum2$group.lib
pheatmap::pheatmap(mat = df.sum2[,chr.seqs], #
                   annotation_colors = list( precursor=RNA_colors[c(rna2,dna)], cell=colorList[["cell"]][cellShort], lineage=colorList[["lineage"]][lineageShort] ) , # RColorBrewer::brewer.pal(6, "Set3")
                   border_color="NA", 
                   scale = "row", # row,none
                   #labels_col = 3, labels_row = 6,
                   cluster_cols = FALSE,cluster_rows = FALSE,
                   show_colnames=TRUE, show_rownames=TRUE,
                   fontsize = 10,
                   legend =  TRUE,  # Uncomment this line
                   height = 4, width = 6, # 50,30
                   color = rev( colorRampPalette( RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))( 20 ) ) , # ),  OrRd
                   # breaks = brkList, 
                   #fontsize_row = 5,
                   filename = paste0("./figure/",dst,"_peakNum_byChr_heat.png")
)
#


#opt2: by RNA
df.sum2 <- df.sum %>% 
  # dplyr::filter(RNA=="tRNA") %>% 
  dplyr::select(group.lib, RNA, gn.chr, mean.total.number) %>% 
  tidyr::pivot_wider(names_from = "RNA",id_cols = "group.lib",
                     values_from = "mean.total.number",
                     values_fn = sum
  )
for(i in c(rna2,dna)){
  if(!(i %in% colnames(df.sum2))){
    message("add ",i)
    df.sum2[[i]] <- 0
  }
}
df.sum2 <- as.data.frame(df.sum2)
rownames(df.sum2) <- df.sum2$group.lib
rownames(df.sum2) <- gsub("Tumor ","",rownames(df.sum2))
df.sum2 <- df.sum2[cellShort,]
annotation_col <- data.frame(row.names = c(rna2,dna), precursor=c(rna2,dna) )
pheatmap::pheatmap(mat = df.sum2[,c(rna2,dna)], # chr.seqs
                   annotation_col = annotation_col, #data frame that specifies the annotations shown on left side of the heatmap.
                   annotation_row = cell.table[,c("lineage","cell")],
                   annotation_colors = list( precursor=RNA_colors[c(rna2,dna)], cell=colorList[["cell"]][cellShort], lineage=colorList[["lineage"]][lineageShort] ) , # RColorBrewer::brewer.pal(6, "Set3")
                   border_color="NA",
                   scale = "row", # row,none
                   #labels_col = 3, labels_row = 6,
                   cluster_cols = FALSE,cluster_rows = FALSE,
                   # gaps_row = sum(annotation_row$trend=="up"), # invalid if cluster row/col
                   #cutree_cols = 2,cutree_rows = 3,
                   show_colnames=TRUE, show_rownames=TRUE,
                   fontsize = 10,
                   legend =  FALSE,  # Uncomment this line
                   height = 4, width = 6, # 50,30
                   color = rev( colorRampPalette( RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))( 20 ) ) , # ),  OrRd
                   # breaks = brkList, 
                   #fontsize_row = 5,
                   filename = paste0("./figure/",dst,"_peakNum_byRNA_heat.png")
)
#




#### peak stat
combinations <- outer(cellShort, c("Normal","Tumor"), FUN = function(x, y) paste0(x," ",y))
combinations <- as.vector(t(combinations))

df.sum <- df %>% 
  dplyr::group_by(RNA,group.lib) %>% 
  dplyr::summarize(mean.ratio=mean(RNA.peak.ratio,trim=0.1,na.rm=TRUE),mean.total.number=mean(total.number,na.rm=TRUE,trim=0.1)) %>%  # average biol rep within same group.lib grp
  dplyr::group_by(group.lib) %>% 
  dplyr::mutate(sum.ratio=sum(mean.ratio),norm.mean.ratio=mean.ratio/sum.ratio)
df.sum$group.lib <- factor(df.sum$group.lib, levels = combinations) 
p1 <- ggplot(df.sum,aes(x=group.lib,y=norm.mean.ratio,fill=RNA))+ 
  geom_bar(position = "stack",width = 0.8,size=0.3, stat = "identity", color="black") +
  #geom_errorbar()+
  xlab("")+
  ylab(paste0("Peak number ratio"))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  scale_fill_manual(name="Species",values = RNA_colors) + # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
  theme_bw() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.title.y = element_text(size = 24,color ="black",angle = 90),  # ,vjust = 1
        axis.text = element_text(size= 20,color = "black"),
        axis.text.x = element_text(size= 16,color = "black",angle = 90,vjust = 0.5,hjust = 1),
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = c(rep("grey50",length(other)),rep("steelblue",length(rna)),rep("firebrick",length(dna)))), #  
        # legend.position = "none",#c(.25,.6),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))
# ggsave("./PNK_group_small_domainNum_ratio.pdf",p1,width = 8,height = 10)

df.sum2 <- df.sum %>% 
  dplyr::group_by(group.lib) %>% 
  dplyr::summarize(sum.mean.total.number=sum(mean.total.number))
df.sum2$group.lib <- factor(df.sum2$group.lib, levels = combinations) 
p2 <- ggplot(df.sum2,aes(x=group.lib,y=sum.mean.total.number))+ 
  geom_bar(position = "stack",width = 0.8,size=0.3, stat = "identity", color="black") +
  #geom_errorbar()+
  xlab("")+
  ylab(paste0("Peak number"))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  #ggsci::scale_fill_d3()+
  scale_fill_manual(name="Species",values = RNA_colors) + # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
  theme_void() + 
  theme(
    plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
    axis.title = element_text(size = 24,color ="black"), 
    axis.title.y = element_text(size = 24,color ="black",angle = 90),  # ,vjust = 1
    axis.text = element_text(size= 20,color = "black",hjust = 1),
    axis.text.x = element_blank(),#element_text(size= 16,color = "black",angle = 90,vjust = 0.5,hjust = 1),
    #panel.grid=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x = element_blank(),
    # panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
    #panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    #axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = c(rep("grey50",length(other)),rep("steelblue",length(rna)),rep("firebrick",length(dna)))), #  
    legend.position = "none",#c(.25,.6),
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16),
    plot.margin = unit(x=c(10,0,-20,0),units="pt") ) # top, right, bottom, left
#cowplot::plot_grid(plotlist = list(p2,p1), nrow = 2, align = "h") # , rel_hights = c(1, 1.3),  common.legend = TRUE, legend = "right"
#ggarrange(plotlist = list(p2,p1), align = "hv",ncol = 1)
p <- ggpubr::ggarrange(p2,p1,ncol = 1,nrow = 2,common.legend = TRUE, align = "v",legend="right",heights=c(0.07,0.2))
# +theme(aspect.ratio = 1.2)
#p
#dev.off()
ggsave( paste0("figure/",dst,"_group_small_domain_num.pdf"),p,width = 11,height = 9)
#





#### tx stat
df[1:3,]
df.sum <- df %>% 
  dplyr::group_by(RNA,group.lib) %>% 
  dplyr::summarize(mean.ratio=mean(RNA.tx.ratio,trim=0.1,na.rm=TRUE)) %>%  # average biol rep within same group.lib grp
  dplyr::group_by(group.lib) %>% 
  dplyr::mutate(sum.ratio=sum(mean.ratio),norm.mean.ratio=mean.ratio/sum.ratio)
df.sum$group.lib <- factor(df.sum$group.lib, levels = combinations) 
p1 <- ggplot(df.sum,aes(x=group.lib,y=norm.mean.ratio,fill=RNA))+ 
  geom_bar(position = "stack",width = 0.8,size=0.3, stat = "identity", color="black") +
  #geom_errorbar()+
  xlab("")+
  ylab(paste0("Tx number ratio"))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  scale_fill_manual(name="Species",values = RNA_colors) +
  theme_bw() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.title.y = element_text(size = 24,color ="black",angle = 90),  # ,vjust = 1
        axis.text = element_text(size= 20,color = "black"),
        axis.text.x = element_text(size= 16,color = "black",angle = 90,vjust = 0.5,hjust = 1),
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = c(rep("grey50",length(other)),rep("steelblue",length(rna)),rep("firebrick",length(dna)))), #  
        # legend.position = "none",#c(.25,.6),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))
#p1# ggsave("./PNK_group_small_domainNum_ratio.pdf",p1,width = 8,height = 10)

df[1:3,]
df.sum2 <- df %>% 
  dplyr::group_by(group.lib) %>% 
  dplyr::summarize(sum.mean.total.number=n_distinct(chr))
df.sum2$group.lib <- factor(df.sum2$group.lib, levels = combinations) 
p2 <- ggplot(df.sum2,aes(x=group.lib,y=sum.mean.total.number))+ 
  geom_bar(position = "stack",width = 0.8,size=0.3, stat = "identity", color="black") +
  #geom_errorbar()+
  xlab("")+
  ylab(paste0("Tx number"))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  #ggsci::scale_fill_d3()+
  scale_fill_manual(name="Species",values = RNA_colors) + # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
  theme_void() + 
  theme(
    plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
    axis.title = element_text(size = 24,color ="black"), 
    axis.title.y = element_text(size = 24,color ="black",angle = 90),  # ,vjust = 1
    axis.text = element_text(size= 20,color = "black",hjust = 1),
    axis.text.x = element_blank(),#element_text(size= 16,color = "black",angle = 90,vjust = 0.5,hjust = 1),
    #panel.grid=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x = element_blank(),
    # panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
    #panel.grid.minor.y = element_blank(),
    panel.border = element_blank(),
    #axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = c(rep("grey50",length(other)),rep("steelblue",length(rna)),rep("firebrick",length(dna)))), #  
    legend.position = "none",#c(.25,.6),
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16),
    plot.margin = unit(x=c(10,0,-20,0),units="pt") ) # top, right, bottom, left
#p2
p <- ggpubr::ggarrange(p2,p1,ncol = 1,nrow = 2,common.legend = TRUE, align = "v",legend="right",heights=c(0.07,0.2))
ggsave(paste0("figure/",dst,"_group_small_tx_num.pdf"),p,width = 11,height = 9)
#


sample.table <- sample.table0[sample.table0$sample %in% colnames(logcpm0),]
rownames(sample.table) <- sample.table$sample

colnames(sample.table)
sample.table <- sample.table[,c("cell","group","tumor_stage","years_smoked","days_to_death","gender","race","age","estimate","purity")] # "LUMP","IHC","immune"
sample.table <- sample.table[order(sample.table$cell,sample.table$group,sample.table$tumor_stage,sample.table$race),]
#dim(logcpm0)
log2cpm <- logcpm0[,rownames(sample.table)]
cor.mat <- cor(log2cpm,method = "spearman")
cor.mat <- as.data.frame(cor.mat)
dim(cor.mat)
cor.mat <- cor.mat[rownames(sample.table),rownames(sample.table)]


brkList <- seq(-1, 1, by = 0.1)
colorList.tmp <- colorList
colorList.tmp[['cell']] <- colorList.tmp[['cell']][names(colorList.tmp[['cell']]) %in% unique(sample.table$cell)]
pheatmap::pheatmap(mat = cor.mat,
                   annotation_col = sample.table, #data frame that specifies the annotations shown on left side of the heatmap.
                   annotation_row = sample.table,
                   annotation_colors = colorList.tmp, # RColorBrewer::brewer.pal(6, "Set3")
                   border_color="NA",
                   scale = "none", # row,none
                   #labels_col = 3, labels_row = 6,
                   cluster_cols = FALSE,cluster_rows = FALSE,
                   # gaps_row = sum(annotation_row$trend=="up"), # invalid if cluster row/col
                   #cutree_cols = 2,cutree_rows = 3,
                   show_colnames=FALSE, show_rownames=FALSE,
                   fontsize = 10,
                   legend =  TRUE,  # Uncomment this line
                   height = 12, width = 14, # 50,30
                   color = rev( colorRampPalette( RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))(length(brkList)) ) , # ),  OrRd
                   breaks = brkList, 
                   #fontsize_row = 5,
                   filename = paste0("./figure/",dst,"_sample_cor_heat.png")
                   # filename = paste0("./output/",dst,"/allSmpP01_exPeak_smallDomain_diff_",disease.label,"vs",normal.label,".cpm.heatmap.pdf")
)
#             


# basic stat (peak): dataset peak length/num-per-tx/full-len-ratio hist -------
ref <- data.table::fread(file.path(DATA_DIR, "exSeek-dev/genome/hg38/chrom_sizes/tx_gn_length_newTxID.txt"),header = TRUE,sep = "\t", stringsAsFactors = FALSE)
ref <- ref[!is.na(ref$tx.length),]
colnames(ref)

## read consensus bed
bed <- paste0(DATA_DIR,"/output/",dst,"/call_peak_all/cfpeakCNN/b5_d50_p1.bed")

df <- as.data.frame(rtracklayer::import.bed(bed))
hist(df$width,breaks = 100,xlim=c(0,50))
df$RNA <- peak.table$gn_type[match(df$name,peak.table$name)]
table(df$RNA )
df$RNA <- factor(df$RNA,levels = c(rna2,dna))
df$tx.len <- ref$tx.length[match(df$seqnames,ref$transcript_id)]
df$tx.ratio <- df$width/df$tx.len
df[1:3,]

## across peak width
stack.all <- as.data.frame(table(df$RNA,df$width))
#hist(df$width,breaks = 100,xlim=c(0,50))
colnames(stack.all) <- c("Precursor","Width","Freq")
#str(stack.all)
stack.all$Width <- as.numeric(unfactor(stack.all$Width))
stack.all$Precursor <- factor(stack.all$Precursor,levels = c(rna2,dna))

# library(tidyverse)
lab.df <- stack.all %>% 
  dplyr::group_by(Precursor) %>% 
  dplyr::summarise(Num=sum(Freq),y.max=max(Freq))
#lab.df$label <- as.character(lab.df$Num/sum(lab.df$Num))
lab.df$label <- paste0(round(lab.df$Num/sum(lab.df$Num)*100,digits = 2),"%")
#stack.all[stack.all$Width<=10,]
p.split <- ggplot(stack.all, aes(x = Width, y = Freq, fill = Precursor)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(. ~ Precursor, ncol = 4,scales='free_y') + # as.table = FALSE
  geom_text(data = lab.df, mapping=aes(x=20, y=0.8*y.max,label = label),hjust = 0,vjust = 0, size = 6,  color="black") +
  scale_fill_manual(values = RNA_colors) + 
  scale_x_continuous(breaks = seq(10,65,20), labels = seq(10,65,20), limits = c(0,65)) +
  labs(x = "Length (nt)", y = "Peak number", title = paste0("total: ",sum(stack.all$Freq)) ) + 
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 1,
    axis.line = element_line()
    , panel.grid.minor = element_blank() 
    , panel.grid.major = element_blank()
    , panel.background = element_rect(fill = "transparent", colour = NA)
    , plot.background = element_rect(fill = "transparent", colour = NA)
    , panel.border = element_rect(fill = "transparent"),
    # axis.text = element_text(size = 16, colour="black"),
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 16,hjust = 0.5,vjust = 0.5,angle = 0), #hjust = 1,vjust = 0.5,angle = 90
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "none", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
#p.split
ggsave(plot = p.split, filename = paste0("figure/",dst,"_peak_len.png"), width = 10,height = 10)
#




## CPM (202406)
sample.table <- sample.table0[sample.table0$sample %in% colnames(logcpm0),]
rownames(sample.table) <- sample.table$sample

colnames(sample.table)
sample.table <- sample.table[,c("cell","group","tumor_stage","race")]
sample.table <- sample.table[order(sample.table$cell,sample.table$group,sample.table$tumor_stage,sample.table$race),]
#dim(logcpm0)
logcpm <- logcpm0[,rownames(sample.table)]
logcpm <- 2^logcpm-1
logcpm <- log10(logcpm+1)

rownames(logcpm) <- unlist(sapply(strsplit(x = rownames(logcpm), "|", fixed = TRUE), "[", 4))
tmp <- data.frame(name=rownames(logcpm),logCPM=rowMax(as.matrix(logcpm)))
tmp$logCPM <- round(tmp$logCPM, digits = 1)
df2 <- as.data.frame(left_join(df,tmp,by="name"))
#df2$logCPM <- NULL
stack.all <- as.data.frame(table(df2$RNA,df2$logCPM))
colnames(stack.all) <- c("Precursor","Width","Freq")
stack.all$Width <- as.numeric(unfactor(stack.all$Width))
stack.all$Precursor <- factor(stack.all$Precursor,levels = c(rna2,dna))

# library(tidyverse)
lab.df <- stack.all %>% 
  dplyr::group_by(Precursor) %>% 
  dplyr::summarise(Num=sum(Freq),y.max=max(Freq))
#lab.df$label <- as.character(lab.df$Num/sum(lab.df$Num))
lab.df$label <- paste0(round(lab.df$Num/sum(lab.df$Num)*100,digits = 2),"%")

summary(stack.all$Width)
length(unique(stack.all$Width))
p.split <- ggplot(stack.all,aes(x = Width, y = Freq, fill = Precursor)) + 
  geom_bar(stat = "identity", width = 0.1) + 
  facet_wrap(. ~ Precursor, ncol = 4,scales='free_y') + # as.table = FALSE
  geom_text(data = lab.df, mapping=aes(x=1, y=0.8*y.max, label = label),hjust = 0,vjust = 0, size = 6,  color="black") +
  scale_fill_manual(values = RNA_colors) + 
  scale_x_continuous(breaks = seq(0,4,1), labels = seq(0,4,1), limits = c(0,5)) +
  labs(x = "log10CPM", y = "Peak number", title = paste0("total: ",sum(stack.all$Freq)) ) + 
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 1,
    axis.line = element_line()
    , panel.grid.minor = element_blank() 
    , panel.grid.major = element_blank()
    , panel.background = element_rect(fill = "transparent", colour = NA)
    , plot.background = element_rect(fill = "transparent", colour = NA)
    , panel.border = element_rect(fill = "transparent"),
    # axis.text = element_text(size = 16, colour="black"),
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 16,hjust = 0.5,vjust = 0.5,angle = 0), #hjust = 1,vjust = 0.5,angle = 90
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "none", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
#p.split
ggsave(plot = p.split, filename = paste0("figure/",dst,"_peak_log10cpm.png"),width = 10,height = 10,limitsize = TRUE)
#



## across peak full-len-ratio
df$tx.ratio <- round(df$tx.ratio, digits = 2)
length(unique(df$tx.ratio))

stack.all <- as.data.frame(table(df$RNA,df$tx.ratio))
#stack.all[1:3,]
colnames(stack.all) <- c("Precursor","Width","Freq")
#str(stack.all)
stack.all$Width <- as.numeric(unfactor(stack.all$Width))
#stack.all$Precursor <- factor(stack.all$Precursor,levels = c(dna,rna))
stack.all$Precursor <- factor(stack.all$Precursor,levels = c(rna2,dna))

# library(tidyverse)
lab.df <- stack.all %>% 
  dplyr::group_by(Precursor) %>% 
  dplyr::summarise(Num=sum(Freq),y.max=max(Freq))
#lab.df$label <- as.character(lab.df$Num/sum(lab.df$Num))
lab.df$label <- paste0(round(lab.df$Num/sum(lab.df$Num)*100,digits = 2),"%")

summary(stack.all$Width)
#str(stack.all)
p.split <- ggplot(stack.all, aes(x = Width, y = Freq, fill = Precursor)) + 
  geom_bar(stat = "identity", width = 0.1) + 
  facet_wrap(. ~ Precursor, ncol = 4,scales='free_y') + # as.table = FALSE
  geom_text(data = lab.df,mapping=aes(x=0.2, y=0.8*y.max, label = label),hjust = 0,vjust = 0, size = 6,  color="black") +
  scale_fill_manual(values = RNA_colors) + 
  scale_x_continuous(breaks = seq(0,1,0.5), labels = c("0","0.5","1") ) + # , limits = c(-0.001,1.001)  #add limits bug ?
  labs(x = "Full-length ratio [0,1]", y = "Peak number", title = paste0("total: ",sum(stack.all$Freq)) ) + 
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 1,
    axis.line = element_line()
    , panel.grid.minor = element_blank() 
    , panel.grid.major = element_blank()
    , panel.background = element_rect(fill = "transparent", colour = NA)
    , plot.background = element_rect(fill = "transparent", colour = NA)
    , panel.border = element_rect(fill = "transparent"),
    # axis.text = element_text(size = 16, colour="black"),
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 16,hjust = 0.5,vjust = 0.5,angle = 0), #hjust = 1,vjust = 0.5,angle = 90
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "none", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
#p.split
ggsave(plot = p.split, filename = paste0("figure/",dst,"_peak_lenRatio.png"),width = 10,height = 10,limitsize = TRUE)
#



## tx num across per tx peak-number
stack.all <- as.data.frame(table(df$RNA,df$seqnames))
stack.all <- stack.all[stack.all$Freq!=0,]
table(duplicated(stack.all$Var2)) # Not ALL FALSE for gn assign. peaks of the same tx may assign to diff gn_type 
summary(stack.all$Freq)
stack.all <- as.data.frame(table(stack.all$Var1,stack.all$Freq))
colnames(stack.all) <- c("Precursor","Width","Freq") # changed order !
#str(stack.all)
stack.all$Width <- as.numeric(unfactor(stack.all$Width))
#stack.all$Precursor <- factor(stack.all$Precursor,levels = c(dna,rna))
stack.all$Precursor <- factor(stack.all$Precursor,levels = c(rna2,dna))
#stack.all$Freq <- as.character(stack.all$Freq)

# library(tidyverse)
lab.df <- stack.all %>% 
  dplyr::group_by(Precursor) %>% 
  dplyr::summarise(Num=n_distinct(Freq),y.max=max(Freq))
#lab.df$label <- as.character(lab.df$Num/sum(lab.df$Num))
lab.df$label <- paste0(round(lab.df$Num/sum(lab.df$Num)*100,digits = 2),"%")
#str(stack.all)
summary(stack.all$Width)
summary(stack.all$Freq)
p.split <- ggplot(stack.all, aes(x = Width, y = Freq, fill = Precursor)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(. ~ Precursor, ncol = 4,scales='free_y') + # as.table = FALSE
  geom_text(data = lab.df, mapping=aes(x=2, y=0.8*y.max,label = label),hjust = 0,vjust = 0, size = 6,  color="black") +
  scale_fill_manual(values = RNA_colors) + 
  scale_x_continuous(breaks = seq(0,10,5), labels = seq(0,10,5), limits = c(0,10)) +
  labs(x = "Peak num per tx", y = "Tx number", title = paste0("total: ",sum(stack.all$Freq) ) ) + 
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 1,
    axis.line = element_line()
    , panel.grid.minor = element_blank() 
    , panel.grid.major = element_blank()
    , panel.background = element_rect(fill = "transparent", colour = NA)
    , plot.background = element_rect(fill = "transparent", colour = NA)
    , panel.border = element_rect(fill = "transparent"),
    # axis.text = element_text(size = 16, colour="black"),
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 16,hjust = 0.5,vjust = 0.5,angle = 0), #hjust = 1,vjust = 0.5,angle = 90
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "none", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
#p.split
ggsave(plot = p.split, filename = paste0("figure/",dst,"_tx_peakNum.png"), width = 10,height = 10)
#








# cell/tissue/stage specific peak bed jaccard/similarity mat ------
# Pairwise BED Jaccard between conditions is precomputed by
# scripts/bash/01_bed_jaccard.sh.
pre <- DATA_DIR

jac <- read.table(paste0(DATA_DIR,"/exSeek-dev/test/",dst,"/jaccard/jac.txt"),header = TRUE,check.names = FALSE,sep = "\t",quote = "",fill=TRUE)
jac$from1 <- unlist(sapply(strsplit(jac$from,"/",fixed = TRUE),"[",9))
jac$to1 <- unlist(sapply(strsplit(jac$to,"/",fixed = TRUE),"[",9))
jac$from2 <- unlist(sapply(strsplit(jac$from,"/",fixed = TRUE),"[",10))
jac$to2 <- unlist(sapply(strsplit(jac$to,"/",fixed = TRUE),"[",10))
jac$from2[grepl("Tumor",jac$from2)] <- "Tumor"
jac$from2[grepl("Normal",jac$from2)] <- "Normal"
jac$to2[grepl("Tumor",jac$to2)] <- "Tumor"
jac$to2[grepl("Normal",jac$to2)] <- "Normal"
#table(jac$to2)
jac$from <- paste0(jac$from1,".",jac$from2)
jac$to <- paste0(jac$to1,".",jac$to2)
jac$from <- gsub("TCGA-|_small","",jac$from,perl = TRUE)
jac$to <- gsub("TCGA-|_small","",jac$to,perl = TRUE)
jac$from1 <- gsub("TCGA-|_small","",jac$from1,perl = TRUE)
jac$to1 <- gsub("TCGA-|_small","",jac$to1,perl = TRUE)
jac$overlap <- 0
for(i in 1:nrow(jac)){
  jac$overlap[i] <- max(jac$overlap_n_from[i],jac$overlap_n_to[i])
}
jac$jaccard2 <- jac$overlap/(jac$from_only_n+jac$to_only_n+jac$overlap) # num jac

jac$min <- 0
for(i in 1:nrow(jac)){
  jac$min[i] <- min(jac$to_n[i],jac$from_n[i])
}
jac$Simpson <- jac$overlap/(jac$min) # Overlap Coefficient (Szymkiewicz-Simpson Coefficient)


jac$from1 <- factor(jac$from1,levels = cellShort)
jac$to1 <- factor(jac$to1,levels = cellShort)
jac$from2 <- factor(jac$from2,levels = c("Normal","Tumor"))
jac$to2 <- factor(jac$to2,levels = c("Normal","Tumor"))
jac <- jac[order(jac$from1,jac$from2,jac$to1,jac$to2,decreasing = FALSE),]
cor.test(jac$jaccard, jac$jaccard2) # 0.98
cor.test(jac$jaccard, jac$Simpson) # 0.54

col1 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
tmp.df <- cell.table
tmp.df$g1 <- "Tumor"
tmp.df2 <- cell.table
tmp.df2$g1 <- "Normal"
tmp.df <- as.data.frame(rbind(tmp.df,tmp.df2))
tmp.df$cell <- factor(tmp.df$cell,levels = cellShort)
tmp.df$g1 <- factor(tmp.df$g1,levels = c("Normal","Tumor"))
tmp.df <- tmp.df[order(tmp.df$cell,tmp.df$g1,decreasing = FALSE),]
this.order <- paste0(tmp.df$cell,".",tmp.df$g1)

## plot num jaccard mat
library(corrplot)
tmp1 <- jac[,c("from","to","jaccard2")]
tmp2 <- jac[,c("to","from","jaccard2")]
tmp <- flat2CorrMat(rbind(tmp1,tmp2))
tmp <- tmp[this.order,this.order]
dir.create(paste0("output/",dst,"/figures/"),recursive = TRUE,showWarnings = FALSE)
pdf(paste0("./figure/peak_num_jac.pdf"), width = 12,height = 12)
p1 <- corrplot(tmp, type = "upper", # order = "hclust", # is.corr = FALSE,
               col= rev(col2(200)), #cl.lim = c(0, 0.9), 
               na.label.col = "grey90",na.label	= "square",
               tl.col = "black", tl.srt = 45)
# p1
dev.off()

## plot Simpson mat
library(corrplot)
tmp1 <- jac[,c("from","to","Simpson")]
tmp2 <- jac[,c("to","from","Simpson")]
tmp <- flat2CorrMat(rbind(tmp1,tmp2))
tmp <- tmp[this.order,this.order]
dir.create(paste0("output/",dst,"/figures/"),recursive = TRUE,showWarnings = FALSE)
pdf(paste0("./figure/peak_Simpson.pdf"), width = 12,height = 12)
p2 <- corrplot(tmp, type = "lower", # order = "hclust", # is.corr = FALSE,
               col= rev(col2(200)), #cl.lim = c(0, 0.9), 
               na.label.col = "grey90",na.label	= "square",
               tl.col = "black", tl.srt = 45)
# p2
dev.off()



## combine
p1[!upper.tri(p1)] <- 0
p2[!lower.tri(p2)] <- 0
p <- p1+p2
pdf(paste0("./figure/peak_jac.pdf"), width = 7,height = 7)
corrplot(p, #type = "upper",add = TRUE,   # order = "hclust", # is.corr = FALSE,
         col= rev(col2(200)), #cl.lim = c(0, 0.9),
         na.label.col = "grey90",na.label	= "square",
         tl.col = "black" #tl.srt = 45,
         # cl.pos = "n" #no legend
)
dev.off()
#









# diff TCGA_small16 summary (top/sig heatmap, volcano scatter) --------------
setwd(PROJECT_DIR)
outDir = paste0(diffpre,"/output/",dst,"/")
dir.create(outDir,recursive = TRUE,showWarnings = FALSE)
outSurfix="_rmPurity"
pval.cutoff1 <- 0.01
padj.cutoff1 <- 0.01
pval.cutoff2 <- 0.01 
padj.cutoff2 <- 1 # set to 0.1 not change much
fc.cutoff <- 0.6

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
    tmp <- data.table::fread(x, header = TRUE, sep="\t",check.names = FALSE,stringsAsFactors = FALSE, data.table = FALSE) # _extLibSize
    # tmp$feature <- tmp.merge
    # data.table::fwrite(tmp, file = x, quote = FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    rownames(tmp) <- tmp$feature
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
common_genes <- Reduce(intersect, lapply(diff.sum.list, rownames)) # merge all data sets and keep the shared genes
cat("Common RBP genes across all tumor types:", length(common_genes), "\n")
#tmp <- all_samples_data[["COAD"]]
diff.sum.df <- do.call(cbind,unname(lapply(diff.sum.list, function(x) x[common_genes, ])))
# diff.sum.df <- as.data.frame(do.call(cbind,diff.sum.list))
dim(diff.sum.df)
tmp <- diff.sum.df[,grepl("trend",colnames(diff.sum.df))]
hist(rowSums(tmp!="insig"))
tmp <- tmp[rowSums(tmp!="insig")>=1,]
#tmp
data.table::fwrite(cbind(feature=rownames(tmp),tmp), paste0("./output/",dst,"/diff",outSurfix,".diffSig"),quote = FALSE,row.names = FALSE,col.names = TRUE,sep = '\t')
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
tmp.df$cell <- unlist(sapply(strsplit(tmp.df$name,"|",fixed=TRUE),"[",1))
tmp.df$diff <- unlist(sapply(strsplit(tmp.df$name,"|",fixed=TRUE),"[",2))
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
data.table::fwrite(cbind(cell=rownames(tmp.df),tmp.df), paste0("./output/",dst,"/diff",outSurfix,".diffSig.summary"),quote = FALSE,row.names = FALSE,col.names = TRUE,sep = '\t')




## opt2: by topNum
outSurfix="_rmPurity" # _rmMonthOperator
for(cellID in cellShort){
  # cellID <- "STAD"
  message(cellID)
  
  diff <- list()
  peak.list <- list()
  seqs <- c(50,500,1000) # 100,200,500,1000
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
    #diff[[disease.label]] <- read.table(x, header = TRUE, sep="\t",check.names = FALSE,stringsAsFactors = FALSE)
    tmp <- data.table::fread(x, header = TRUE, sep="\t",check.names = FALSE,stringsAsFactors = FALSE,data.table = FALSE)
    rownames(tmp) <- tmp$feature
    tmp$feature <- NULL
    diff[[disease.label]] <- tmp
    tmp <- tmp[order(as.numeric(tmp$pvalue),decreasing = FALSE),]
    rownames(tmp) <- unlist(sapply(strsplit(rownames(tmp),"|",fixed=TRUE),"[",4))
    
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
all(rownames(diff[[disease.label]])==rownames(logcpm0)) # all TRUE
#



## heatmap 3grps
seqs <- c(50,500) # c(500,1000)
for(cellID in cellShort){
  # cellID <- "STAD"
  message(cellID)
  
  for(top in paste0("top",seqs)){
    # top <- "top50"
    print(top)
    
    tmp.id <- unlist(sapply(strsplit(rownames(logcpm0),"|",fixed=TRUE),"[",4))
    sample.table <- sample.table0[sample.table0$cellShort==cellID,] # 60
    #table(sample.table$group1)
    sample.table$group <- sample.table[["group1"]] # col.grp
    sample.table$group <- factor(sample.table$group,levels = c("Normal","Primary","Metastasis")) # "DC",
    sample.table <- sample.table[order(sample.table$group),]
    logcpm <- logcpm0[tmp.id %in% peak.ids[[top]],colnames(logcpm0) %in% sample.table$sample] # rownames(diff)[diff$pvalue<=0.01]
    #table(diff$pvalue<=0.01)
    logcpm <- logcpm[,sample.table$sample]
    
    table(sample.table$sample %in% colnames(logcpm))
    
    annotation_col <- sample.table[colnames(logcpm),]
    annotation_col <- annotation_col[,c("cell","group","tumor_stage","gender","race"),drop=FALSE] # "group1","group2","age","sex","ExprDate","Month","Operator","MonthOperator"
    
    #txid.RNA <- unlist(sapply(strsplit(rownames(logcpm),"|",fixed=TRUE),"[",2))
    #rownames(logcpm)[1]
    txid.RNA <- unlist(sapply(strsplit(rownames(logcpm),"|",fixed=TRUE),"[",8))
    #txid.RNA <- gsub("_for|_rev","",txid.RNA,perl = TRUE)
    #table(txid.RNA)
    annotation_row <- data.frame( RNA = txid.RNA, row.names = rownames(logcpm))
    annotation_row$RNA <- gsub("_for|_rev","",annotation_row$RNA,perl = TRUE)
    tmp.id2 <- unlist(sapply(strsplit(rownames(logcpm),"|",fixed=TRUE),"[",4))
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
    annotation_row$RNA <- factor(annotation_row$RNA, levels = c(rna2,dna))
    annotation_row <- annotation_row[order(annotation_row$RNA),]
    table(annotation_row$RNA)
    
    RNA_colors_tmp <- RNA_colors
    names(RNA_colors_tmp) <- tmp$lab[match(names(RNA_colors_tmp),tmp$Var1)]
    RNA_colors_tmp <- RNA_colors_tmp[!is.na(names(RNA_colors_tmp))]
    
    ann_colors <- list("Group"=groupAllShort[unique(as.character(sample.table$group1))], # "Group1"=c(group1Color,group2Color), "Group2"=c(group1Color,group2Color), 
                       "Tissue"=cellColor[unique(as.character(sample.table$cell))],
                       "Stage"=colorList[["tumor_stage"]],
                       "Race"=colorList[["race"]],
                       "Gender"=colorList[["gender"]],
                       "Regulation"=c("Both up"="salmon", "Both down"="steelblue4","up&down"="seagreen","down&up"="chocolate"), # , 
                       "Precursor"=RNA_colors_tmp) # ,"salmon",seagreen,orange,chocolate
    #table(annotation_col$Group)
    annotation_row <- annotation_row[order(annotation_row$trend,annotation_row$RNA,decreasing = FALSE),]
    annotation_row$RNA <- NULL
    logcpm <- logcpm[rownames(annotation_row),] #positive_samples
    colnames(annotation_row) <- c("Regulation","Precursor")
    colnames(annotation_col)[1:5] <- c("Tissue","Group","Stage","Gender","Race")
    
    
    dir.create(paste0("./output/",dst,"/heatmap/"),recursive = TRUE,showWarnings = FALSE)
    pheatmap::pheatmap(mat = logcpm, 
                       annotation_col = annotation_col, #data frame that specifies the annotations shown on left side of the heatmap.
                       annotation_row = annotation_row,
                       annotation_colors = ann_colors, # RColorBrewer::brewer.pal(6, "Set3")
                       border_color="NA",
                       # labels_col = coldata$age
                       scale = "row", # row,none
                       #labels_col = 3, labels_row = 6,
                       cluster_cols = FALSE,cluster_rows = FALSE,
                       # gaps_row = sum(annotation_row$trend=="up"), # invalid if cluster row/col
                       #cutree_cols = 2,cutree_rows = 3,
                       show_colnames=FALSE, show_rownames=FALSE,
                       fontsize = 12,
                       height = 14, width = 14,
                       color = rev( colorRampPalette( RColorBrewer::brewer.pal(n = 10, name = "RdBu"))(100) ),  # steelblue/dodgerblue4-firebrick2, viridis::viridis_pal()(100), #
                       #fontsize_row = 5,
                       filename = paste0("./output/",dst,"/heatmap/",cellID,"_",top,"_diff_",disease.label,"vs",normal.label,".cpm.heatmap2.png")
    )
  }
  
  
  
  ## heatmap 2grps
    
  ## plot diff volcano 
  my_theme <-   theme(aspect.ratio=1,
                      plot.title = element_text(size = 28,color="black",hjust = 0.5),
                      axis.title = element_text(size = 28,color ="black"), 
                      axis.text = element_text(size= 24,color = "black"), #,face="bold
                      panel.grid.minor.y = element_blank(),
                      panel.grid.minor.x = element_blank(),
                      axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
                      axis.text.y = element_text( hjust = 0 ), # angle = 45,
                      panel.grid=element_blank(),
                      legend.position = "right",#c(0.5,0.3),
                      legend.text = element_text(size= 24,color = "black"),
                      legend.title= element_text(size= 26,color = "black"))
  
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
    tmp <- read.table(x,sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
    
    tmp$minus.log10.pvalue <- -log10(tmp$pvalue)
    tmp$log10.baseMean <- log10(tmp$baseMean)
    
    tmp$alpha <- maxmin.normalize.vec(tmp$minus.log10.pvalue)
    
    tmp$RNA <- unlist(sapply(strsplit(tmp$feature,"|",fixed=TRUE),"[",8))
    tmp$RNA <- factor(tmp$RNA,levels=rna2)
    
    
    ## color by pvalue
    diff.plot[[disease.label]] <- ggplot(tmp, aes(x=log10.baseMean,y=log2FoldChange, fill=minus.log10.pvalue)) + 
      geom_point(aes(alpha=alpha),shape=21,size=3,color="white") + 
      labs(title = paste0(disease.label," vs. ",normal.label)) +
      ggraph::scale_fill_viridis(direction = -1) + 
      geom_hline(yintercept = 0,linetype="dashed",color="salmon") + 
      theme_bw()+
      guides(fill = guide_legend(title = "-log10P" ), alpha = "none") +
      my_theme
    ggsave(plot = diff.plot[[disease.label]], paste0("./output/",dst,"/heatmap/",cellID,"_",top,"_diff_",disease.label,"vs",normal.label,"_diff_basemean.png"), device = "png",dpi = 150, width=9, height=9) # "_",sample
    
    ## color by RNA
    diff.plot[[disease.label]] <- ggplot(tmp, aes(x=log10.baseMean,y=log2FoldChange, fill=RNA)) +  # minus.log10.pvalue
      geom_point(aes(alpha=0.95),shape=21,size=3,color="white") + 
      labs(title = paste0(disease.label," vs. ",normal.label)) +
      # ggraph::scale_fill_viridis(direction = -1) + 
      scale_fill_manual(values = RNA_colors)+
      geom_hline(yintercept = 0,linetype="dashed",color="salmon") + 
      theme_bw()+
      guides(fill = guide_legend(title = "-log10P" ), alpha = "none") +
      my_theme
    ggsave(plot = diff.plot[[disease.label]], paste0("./output/",dst,"/heatmap/",cellID,"_",top,"_diff_",disease.label,"vs",normal.label,"_diff_basemean2.png"), device = "png",dpi = 150, width=9, height=9) # "_",sample
  }
}
#






# diff TCGA_small16 signature stats (diff num, jaccar cor heatmap, cell recurrence bar) ------
diff.sum.list <- list()
pval.cutoff1 <- 0.01
padj.cutoff1 <- 0.01
pval.cutoff2 <- 0.01 
padj.cutoff2 <- 1 # set to 0.1 not change much
fc.cutoff <- 0.6
outSurfix="_rmPurity" # _rmMonthOperator

for(cellID in cellShort){
  message(cellID)
  for(i in c("tumorigenesis","metastasis")){
    # i <- "tumorigenesis"
    message(i)
    
    disease.label <- cmp.df[i,"pos.grp"]
    normal.label <- cmp.df[i,"neg.grp"]
    col.grp <- cmp.df[i,"col"]
    
    x <- paste0("./output/",dst,"/diff_",cellID,"_",disease.label,"vs",normal.label,outSurfix,".diff")
    if(!file.exists(x)){
      message(x," not exist")
      next
    }
    #tmp <- read.table(x, header = TRUE, sep="\t",check.names = FALSE,stringsAsFactors = FALSE) # _extLibSize
    tmp <- data.table::fread(x, header = TRUE, sep="\t",check.names = FALSE,stringsAsFactors = FALSE, data.table = FALSE) # _extLibSize
    rownames(tmp) <- tmp$feature
    tmp$feature <- NULL
    tmp$disease <- disease.label
    tmp$normal <- normal.label
    tmp$cell <- cellID
    tmp$trend <- "insig"
    if(i=="tumorigenesis"){
      padj.cutoff <- padj.cutoff1
      pval.cutoff <- pval.cutoff1
    }else{
      padj.cutoff <- padj.cutoff2
      pval.cutoff <- pval.cutoff2
    }
    tmp$trend[tmp$padj<=padj.cutoff & tmp$pvalue<=pval.cutoff & (tmp$log2FoldChange > fc.cutoff)] <- "up"
    tmp$trend[tmp$padj<=padj.cutoff & tmp$pvalue<=pval.cutoff & (tmp$log2FoldChange < fc.cutoff)] <- "down"
    diff.sum.list[[paste0(cellID,"|",disease.label)]] <- tmp
  }
}
#dst
dir.create(paste0("./output/",dst,"/rds/"),recursive = TRUE,showWarnings = FALSE)
diff.sum.df <- as.data.frame(do.call(cbind,diff.sum.list))
diff.sum.df <- diff.sum.df[,grepl("trend|logCPM|log2FoldChange|pvalue|padj|baseMean",colnames(diff.sum.df),perl = TRUE)]
diff.sum.df$feature <- rownames(diff.sum.df)

diff.sum.df.sig <- diff.sum.df[,grepl("feature|trend",colnames(diff.sum.df),perl = TRUE)] #|\\.padj
#diff.sum.df.sig$feature <- rownames(diff.sum.df.sig)
diff.sum.df.sig <- reshape2::melt(diff.sum.df.sig,id.var=c("feature"),variable.name="class",value.name = "trend")
diff.sum.df.sig$cell <- unlist(sapply(strsplit(as.character(diff.sum.df.sig$class), "|", fixed = TRUE), "[", 1))
diff.sum.df.sig$group <- unlist(sapply(strsplit(as.character(diff.sum.df.sig$class), "|", fixed = TRUE), "[", 2))
diff.sum.df.sig$group <- gsub(".trend","",diff.sum.df.sig$group,fixed = TRUE)
diff.sum.df.sig$class <- NULL
diff.sum.df.sig2 <- tidyr::pivot_wider(data = dplyr::as_tibble(diff.sum.df.sig), id_cols = c("feature","cell"), names_from = group, values_from = trend)
#for(cellID in cellShort){
names(diff.sum.df.sig2)

pos.lab <- "Yes"
neg.lab <- "Primary"
diff.sum.df.sig2$DEGtype <- "Both insignificant"
diff.sum.df.sig2$DEGtype[(diff.sum.df.sig2[[pos.lab]] == "insig") & (diff.sum.df.sig2[[neg.lab]] == "insig")] <- "Both insignificant"
diff.sum.df.sig2$DEGtype[(diff.sum.df.sig2[[pos.lab]] %in% c("up","down")) & (diff.sum.df.sig2[[neg.lab]] == "insig")] <- "Metastasis only"
diff.sum.df.sig2$DEGtype[(diff.sum.df.sig2[[pos.lab]] == "insig") & (diff.sum.df.sig2[[neg.lab]] %in% c("up","down"))] <- "Tumorigenesis only"
diff.sum.df.sig2$DEGtype[(diff.sum.df.sig2[[pos.lab]] == "up") & (diff.sum.df.sig2[[neg.lab]] == "up")] <- "Both significant"
diff.sum.df.sig2$DEGtype[(diff.sum.df.sig2[[pos.lab]] == "down") & (diff.sum.df.sig2[[neg.lab]] == "down")] <- "Both significant"
diff.sum.df.sig2$DEGtype[((diff.sum.df.sig2[[pos.lab]] == "down") & (diff.sum.df.sig2[[neg.lab]] == "up")) | 
                           ((diff.sum.df.sig2[[pos.lab]] == "up") & (diff.sum.df.sig2[[neg.lab]] == "down"))] <- "Both significant"
table(diff.sum.df.sig2$DEGtype)
diff.sum.df.sig2$lineage <- cell.table$lineage[match(diff.sum.df.sig2$cell,cell.table$cellShort)]
table(diff.sum.df.sig2$lineage)
#write.table(diff.sum.df.sig2,paste0(diffpre,"/output/",dst,"/diffSig.summary"), quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
data.table::fwrite(x = diff.sum.df.sig2, file = paste0(diffpre,"/output/",dst,"/diffSig.summary"), quote = FALSE,sep = "\t", row.names = FALSE, col.names = TRUE)

library(dplyr)
library(magrittr)
cell.mat <- list()
{ # cellNum
  diff.sum.df.sig3 <- diff.sum.df.sig2 %>%
    # mutate(record=1) %>% 
    tidyr::pivot_wider( id_cols = c("feature"), names_from = cell, values_from = DEGtype)
  
  tmp.ncol <- ncol(diff.sum.df.sig3)
  for(i in unique(diff.sum.df.sig2$DEGtype)){
    # i <- "Both insignificant"
    message(i)
    diff.sum.df.sig3[[paste0(i," cellNum")]] <- rowSums(diff.sum.df.sig3[,2:tmp.ncol]==i)
  }
  cell.mat[["metastasis"]] <- as.data.frame(row.names = diff.sum.df.sig3$feature, diff.sum.df.sig3[,2:tmp.ncol]=="Both significant" | diff.sum.df.sig3[,2:tmp.ncol]=="Metastasis only")
  cell.mat[["tumorigenesis"]] <- as.data.frame(row.names = diff.sum.df.sig3$feature, diff.sum.df.sig3[,2:tmp.ncol]=="Both significant" | diff.sum.df.sig3[,2:tmp.ncol]=="Tumorigenesis only")
}
diff.sum.df.sig3[[paste0("tumorigenesis.signature"," cellNum")]] <- rowSums(cell.mat[["tumorigenesis"]]) # Both significant Metastasis only    Tumorigenesis only
diff.sum.df.sig3[[paste0("metastasis.signature"," cellNum")]] <- rowSums(cell.mat[["metastasis"]]) # Both significant Metastasis only    Tumorigenesis only

lineage.mat <- list()
for(i in c("tumorigenesis","metastasis")){
  tmp <- cell.mat[[i]]
  tmp$feature <- rownames(tmp)
  tmp <- tidyr::pivot_longer(tmp,cols = 1:(ncol(tmp)-1),names_to = "cell")
  tmp$lineage <- cell.table$lineage[match(tmp$cell,cell.table$cellShort)]
  tmp <- tmp %>% 
    dplyr::group_by(feature,lineage) %>% 
    dplyr::summarise(value2=sum(value)) %>% 
    tidyr::pivot_wider(id_cols = c("feature"), names_from = lineage, values_from = value2)
  tmp.ncol <- ncol(tmp)
  lineage.mat[[i]] <- as.data.frame(row.names = tmp$feature, tmp[,2:tmp.ncol]>=1 | tmp[,2:tmp.ncol]>=1)
  lineage.mat[[i]] <- lineage.mat[[i]][diff.sum.df.sig3$feature,]
}
names(lineage.mat)
#all(1:4==NULL)
all(rownames(lineage.mat[['tumorigenesis']])==diff.sum.df.sig3$feature)
diff.sum.df.sig3[[paste0("metastasis",".signature"," lineageNum")]] <- rowSums(lineage.mat[["metastasis"]])
diff.sum.df.sig3[[paste0("tumorigenesis",".signature"," lineageNum")]] <- rowSums(lineage.mat[["tumorigenesis"]])

diff.sum.df$feature <- rownames(diff.sum.df)
all((diff.sum.df$feature)==diff.sum.df.sig3$feature)
#"feature" %in% colnames(diff.sum.df)
diff.sum.df <- dplyr::left_join(diff.sum.df, diff.sum.df.sig3, by="feature") # cbind(diff.sum.df.sig3,diff.sum.df)
diff.sum.df <- diff.sum.df[,c("feature",colnames(diff.sum.df)[colnames(diff.sum.df)!="feature"])]
# write.table(diff.sum.df,paste0(diffpre,"/output/",dst,"/diff.summary"), quote = FALSE,sep = "\t",row.names = TRUE,col.names = TRUE)
data.table::fwrite(diff.sum.df,paste0(diffpre,"/output/",dst,"/diff.summary"), quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
diff.sum.df <- data.table::fread(paste0("./output/",dst,"/diff.summary"), quote = "",sep = "\t",header = TRUE,check.names = FALSE,data.table = FALSE,stringsAsFactors = FALSE)



# fig3B
tmp <- splitPeakFeatureID(unique(diff.sum.df.sig2$feature))
table(tmp$RNA,tmp$RNA2)
tmp$RNA <- tmp$RNA2  #gsub("_for|_rev","",tmp$RNA,perl = TRUE)
diff.sum.df.sig2$RNA <- tmp$RNA[match(diff.sum.df.sig2$feature,tmp$feature)]
df <- as.data.frame(table(diff.sum.df.sig2$cell,diff.sum.df.sig2$DEGtype,diff.sum.df.sig2$RNA))
colnames(df) <- c("Cell","DEGtype","RNA","Number")
df.bak <- df
df.bak <- df.bak %>% 
  tidyr::pivot_wider(id_cols = c("Cell","DEGtype"), names_from = "RNA", values_from = "Number" )
df.bak[1:3,]
df.bak$Number <- apply(df.bak[,3:ncol(df.bak)],1,sum)
data.table::fwrite(df.bak,paste0("./output/",dst,"/cell_diff.summary2"), quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

df <- df %>% 
  dplyr::group_by(Cell,DEGtype) %>% 
  dplyr::summarise(Number=sum(Number))
df <- df[df$DEGtype != "Both insignificant",]
df$Cell <- factor(df$Cell, levels = cellShort)
df$Lineage <- cell.table$lineage[match(df$Cell,cell.table$cellShort)]
df$Lineage <- factor(df$Lineage, levels = lineageShort)
df$signature <- ""
df$signature[df$DEGtype %in% c("Both significant","Metastasis only")] <- "metastasis.signature"
df$signature[df$DEGtype %in% c("Both significant","Tumorigenesis only")] <- "tumorigenesis.signature"
library(dplyr)
table(df$signature)
# add signature.gene & jacard similarity
df2 <- df %>% 
  tidyr::pivot_wider( id_cols = c("Cell","Lineage"), names_from = DEGtype, values_from = Number) %>% 
  mutate(sumNum=(`Both significant` + `Metastasis only` + `Tumorigenesis only`),
         jacard=`Both significant`/sumNum)
df2$jacard[is.na(df2$jacard)] <- 0
jacard.cutoff <- 0.15
sumNum.cutoff <- 5
df2$cellGroup <- ""
df2$cellGroup[df2$`Metastasis only`>df2$`Tumorigenesis only` & df2$jacard<=jacard.cutoff] <- "Metastasis dominant"
df2$cellGroup[df2$`Metastasis only`<df2$`Tumorigenesis only` & df2$jacard<=jacard.cutoff] <- "Tumorigenesis dominant"
df2$cellGroup[df2$jacard>jacard.cutoff | df2$sumNum<=sumNum.cutoff] <- "Shared"
table(df2$cellGroup)
df2$cellGroup <- factor(df2$cellGroup, levels = c('Tumorigenesis dominant','Metastasis dominant','Shared'))
#str(df2)
df2 <- df2[order(df2$cellGroup,df2$sumNum),]
cellShort2 <- unique(df2$Cell)
#write.table(df2,paste0("./output/",dst,"/cell_diff.summary"), quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
data.table::fwrite(df2,paste0("./output/",dst,"/cell_diff.summary"), quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

table(df.bak$DEGtype)
df.bak3 <- df.bak %>% 
  dplyr::select(c("Cell","DEGtype","enhancer":"Y_RNA")) %>% 
  tidyr::pivot_longer(cols = "enhancer":"Y_RNA", names_to = "RNA", values_to = "RNA.Number") %>% 
  dplyr::filter(DEGtype %in% c("Both significant","Metastasis only","Tumorigenesis only")) 
df.bak3$Cell <- factor(df.bak3$Cell, levels = cellShort2) # cellShort2
df.bak3$RNA <- factor(df.bak3$RNA, levels = c(rna2,dna))
#table(df.bak3$RNA)
df.bak2 <- df.bak3 %>% 
  dplyr::group_by(Cell, RNA) %>%
  dplyr::summarise(Number=sum(RNA.Number))
df.bak2$Cell <- factor(df.bak2$Cell, levels = cellShort2) # cellShort2
df.bak2$RNA <- factor(df.bak2$RNA, levels = c(rna2,dna))

my_theme_bar <-   theme(#axis.ticks.x=element_blank(),
  strip.text = element_text(size=20),
  axis.title.x = element_text(size=20),
  axis.title.y = element_text(size=20),
  axis.text.x = element_text(size = 20,vjust = 1,hjust = 1, angle = 45), #
  axis.text.y = element_text(size = 20),
  plot.title = element_text(size=20),
  # strip.text = element_blank(),
  legend.position = "none", #c(0.9,0.8),#,#
  legend.text = element_text(size= 16),
  legend.title= element_text(size= 16))
my_theme_bar_top <-   theme(#axis.ticks.x=element_blank(),
  strip.text = element_text(size=20),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size=20),
  axis.text.x = element_blank(), #
  axis.text.y = element_text(size = 20),
  plot.title = element_text(size=20),
  # strip.text = element_blank(),
  legend.position = "none", #c(0.9,0.8),#,#
  legend.text = element_text(size= 16),
  legend.title= element_text(size= 16)
)

p0 <- 
  ggplot(df.bak2,aes(x=Cell,y=Number,fill=RNA))+ 
  geom_bar(position = "stack",width = 0.8,size=0.3, stat = "identity", color="black") +
  scale_fill_manual(name="Precursor",values = RNA_colors ) +
  labs(title="",x="", y = "", fill="Precursor")+
  # xlim(c(0,1))+   
  # ylim(c(0,1))+
  theme_minimal() +  # base_size=12
  my_theme_bar_top + 
  theme(
    plot.margin = unit(x=c(10,0,-10,0),units="pt") # top, right, bottom, left
  )
p01 <- 
  ggplot(df.bak3,aes(x=Cell,y=RNA.Number,fill=RNA))+ 
  geom_bar(position = "stack",width = 0.8,size=0.3, stat = "identity", color="black") +
  scale_fill_manual(name="Precursor",values = RNA_colors ) +
  # scale_x_discrete(label=c("miRNA","primary miRNA peaks","all RNA peaks","all DNA peaks","all peaks"))+
  labs(title="",x="", y = "", fill="Precursor")+
  # xlim(c(0,1))+   
  # ylim(c(0,1))+
  theme_minimal() +  # base_size=12
  my_theme_bar_top + 
  theme(
    plot.margin = unit(x=c(10,0,-20,0),units="pt") # top, right, bottom, left
    # plot.margin = unit(c(1,1,0,1), 'lines')) 
  )+
  facet_grid(DEGtype~.,scale="free_y")
p01
ggsave(filename = paste0("./figure/sigDiff_arrange0.png"),plot = p01,width = 8,height = 6)

df$Cell <- factor(df$Cell, levels = cellShort2) # cellShort2
p1 <- ggplot(df, aes(x=Cell, y=Number, fill=Lineage)) +  # DEGtype
  geom_bar(stat="summary",color="black"  # identity
           # position=position_fill()
  ) +
  scale_fill_manual(name="Lineage",values = lineageColor ) +
  labs(title="",x="", y = "", fill="Lineage")+
  theme_minimal() +  # base_size=12
  my_theme_bar_top + 
  theme(
    plot.margin = unit(x=c(-10,0,-20,0),units="pt") # top, right, bottom, left
    # plot.margin = unit(c(1,1,0,1), 'lines')) 
  )# facet_grid(lambda~.)

p2 <- ggplot(df, aes(x=Cell, y=Number, fill=DEGtype)) +  # 
  geom_bar(stat="summary",color="black",  # identity
           position=position_fill()
  ) +
  scale_fill_manual(name="DEGtype",values = DEGColor ) +  
  labs(title="",x="", y = "", fill="DEGtype")+
  # xlim(c(0,1))+   
  # ylim(c(0,1))+
  theme_minimal() +  # base_size=12
  my_theme_bar  

p012 <- ggpubr::ggarrange(align = "v",legend="none",common.legend = TRUE,ncol = 1,nrow = 2,heights = c(2,2),plotlist = list(p1,p2)) # p0 # legend = "none",
p012
ggsave(filename = paste0("./figure/sigDiff_arrange.png"),plot = p012,width = 8,height = 6)


# fig3C: cell DEGtype recurrence hist/bar plot
#sanity check
table(diff.sum.df$`metastasis.signature cellNum`>=diff.sum.df$`metastasis.signature lineageNum`)
table(diff.sum.df$`tumorigenesis.signature cellNum`>=diff.sum.df$`tumorigenesis.signature lineageNum`)
tmp.df <- diff.sum.df[,grepl("cellNum|feature",colnames(diff.sum.df),perl = TRUE)]
#colnames(tmp.df)
#tmp.df[1:3,]
tmp.df <- tidyr::pivot_longer(data = tmp.df, cols = 2:ncol(tmp.df),names_to = "cellClass", values_to = "cellNum")

tmp.df2 <- diff.sum.df[,grepl("lineageNum|feature",colnames(diff.sum.df),perl = TRUE)]
colnames(tmp.df2)
#length(unique(tmp.df2$feature))
summary(tmp.df2$`metastasis.signature lineageNum`)

tmp.df$`metastasis.signature lineageNum` <- tmp.df2$`metastasis.signature lineageNum`[match(tmp.df$feature,tmp.df2$feature)] # "metastasis.signature lineageNum" "tumorigenesis.signature lineageNum"
tmp.df$`tumorigenesis.signature lineageNum` <- tmp.df2$`tumorigenesis.signature lineageNum`[match(tmp.df$feature,tmp.df2$feature)] 
#hist(tmp.df$`metastasis.signature cellNum`)
colnames(tmp.df)
#tmp.df <- tidyr::pivot_longer(tmp.df,cols = "metastasis.signature lineageNum":"tumorigenesis.signature lineageNum",names_to = "lineageClass",values_to = "lineageNum")
tmp.df <- tidyr::pivot_longer(tmp.df,cols = c("tumorigenesis.signature lineageNum","metastasis.signature lineageNum"),names_to = "lineageClass",values_to = "lineageNum")
table(tmp.df$cellClass,tmp.df$lineageClass)
tmp.df <- tmp.df[(grepl("metastasis.signature",tmp.df$cellClass) & grepl("metastasis.signature",tmp.df$lineageClass)) | 
                   (grepl("tumorigenesis.signature",tmp.df$cellClass) & grepl("tumorigenesis.signature",tmp.df$lineageClass)),] # tumorigenesis.signature
tmp.df <- tmp.df[tmp.df$cellNum>0,]
colnames(tmp.df)


library(RColorBrewer)
tmp.df$lineageNum <- factor(tmp.df$lineageNum)
table(tmp.df$cellClass)
tmp.df$cellClass[tmp.df$cellClass=="metastasis.signature cellNum"] <- "Metastasis"
tmp.df$cellClass[tmp.df$cellClass=="tumorigenesis.signature cellNum"] <- "Tumorigenesis"
tmp.df$cellClass <- factor(tmp.df$cellClass,levels = c("Tumorigenesis","Metastasis"))
ggplot(tmp.df, aes(x=cellNum, fill=lineageNum)) +  # 
  geom_histogram(binwidth = 0.5, aes(y = ..count../length(unique(tmp.df$feature)))) +
  # scale_fill_nejm_adaptive(guide = guide_legend(reverse=FALSE)) + # name = "Dose",
  scale_fill_brewer(type = "seq", direction = -1, palette = "Greens") + # RColorBrewer::brewer.pal(n = 6,name = "Greens")[3:6]
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), (max(x) + 1) * 1.1))))) +
  # scale_x_discrete(breaks=c(1:6)) +
  # scale_x_discrete(breaks = c("1","2","3","4","5","6")) + 
  labs(title="",x="Number of cell type shared by differential expression", y = "Proportion of DEGs", fill="lineage")+
  # xlim(c(0,1))+   
  # ylim(c(0,1))+
  theme(axis.text.x=element_text(colour="black",size=14),
        axis.text.y=element_text(colour="black",size=14),
        axis.title.x=element_text(colour="black",size=14),
        axis.title.y=element_text(colour="black",size=14),
        plot.title=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(colour="black",size=15),
        panel.background =element_rect(fill="transparent",colour="black",size=1),
        panel.grid = element_blank(),
        legend.position="right",
        legend.title=element_text(colour="black",size=14),
        legend.text=element_text(colour="black",size=14))+
  guides(fill=guide_legend(title="lineage"))+ 
  facet_grid(.~cellClass,)
ggsave(paste0("./figure/sigDiff_fig3C.png"),width = 6,height = 4)



# Heatmap
calcJaccard <- function(x, y) {
  sum(x == TRUE & y == TRUE)/sum(x == TRUE | y == TRUE)  
}

#Split the data column-wise (rows are records) and use outer:
tmp1 <- asplit(cell.mat[["metastasis"]], 2)
jac.mat1 <- outer(tmp1, tmp1, Vectorize(calcJaccard))
#table(rownames(jac.mat1)==cellShort)
tmp2 <- asplit(cell.mat[["tumorigenesis"]], 2)
jac.mat2 <- outer(tmp2, tmp2, Vectorize(calcJaccard))


mat <- matrix(NA, length(cellShort), length(cellShort), dimnames = list(cellShort, cellShort))
#dim(mat)
mat[lower.tri(mat)] <- jac.mat2[lower.tri(jac.mat2, diag = FALSE)] # rows: tumorigenesis
mat[upper.tri(mat)] <- jac.mat1[upper.tri(jac.mat1, diag = FALSE)] # cols: metastasis
mat[is.nan(mat)] <- 0
for(i in 1:nrow(mat)){
  mat[i,i] <- 1
}

table(rownames(mat) %in% cellShort)
mat <- mat[rownames(mat) %in% cellShort,colnames(mat) %in% cellShort]
mat <- mat[cellShort,cellShort]
#all(rownames(mat)==cellShort)
table(is.na(mat))
#mat[is.na(mat)] <- 0
dim(mat)

cormat <- reorder_cormat(mat)
mat.lab <- cormat[["lab"]]
cormat <- cormat[["mat"]]

# jaccard mat plot opt1:
#library(reshape2)
#library(ggplot2)
# melted_cormat <- melt(cormat, na.rm = TRUE) # by each columns
# colnames(melted_cormat)[1:2] <- c("Tumorigenesis","Metastasis")
# #str(melted_cormat)
# #melted_cormat <- melted_cormat[order(as.integer(melted_cormat$`Tumorigenesis`),-as.integer(melted_cormat$`Metastasis`),decreasing = FALSE),]
# ggheatmap <- ggplot(melted_cormat, aes(x = `Tumorigenesis`, y = `Metastasis`, fill = value))+
#   geom_tile(color = "white")+
#   scale_y_discrete(position = "right")+


# jaccard mat plot opt2:
library(ComplexHeatmap)
library(RColorBrewer)
jac_both <- t(cormat)
max.cutoff <- 0.5 # 0,9
jac_both[jac_both>max.cutoff] <- max.cutoff
col2 <- lineageColor
col_fun = colorRamp2::colorRamp2(seq(0,max.cutoff,length=8), rev(brewer.pal(11,"RdYlBu")[1:8]))

#column_labels = structure(cell.table$lineage[match(colnames(cormat),cell.table$cellShort)],names=colnames(cormat))

labels_state2 <- cell.table[match(mat.lab,cell.table$cellShort),]
#str(labels_state2)
# labels_state2 = left_join(data.frame(subset=labels_state),celltype_corresp%>%select(subset,lineage),by="subset")
labels_state2$lineage = factor(labels_state2$lineage,levels=lineageShort) # c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil")

ha = HeatmapAnnotation(lineage=labels_state2$lineage,col=list(lineage=col2),show_annotation_name=FALSE,
                       height=unit(0.2,"cm"),simple_anno_size_adjust=TRUE,show_legend=FALSE,
                       annotation_legend_param=list(title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=12)))
ra =     rowAnnotation(lineage=labels_state2$lineage,col=list(lineage=col2),show_annotation_name=FALSE,
                       width=unit(0.2,"cm"),simple_anno_size_adjust=TRUE,show_legend=FALSE,
                       annotation_legend_param=list(title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=12)))

p=Heatmap(jac_both,cluster_rows=FALSE,cluster_columns=FALSE,,col=col_fun,width=unit(12,"cm"),height=unit(12,"cm"),
          row_names_gp=gpar(fontsize=12), #row_labels=column_labels[rownames(jac_both)],
          column_names_gp=gpar(fontsize=12), #column_labels=column_labels[colnames(jac_both)],
          bottom_annotation = ha,
          right_annotation = ra,
          heatmap_legend_param=list(title="Jaccard similarity index",at=seq(0,max.cutoff,by=0.1), # c(0,0.2,0.4,0.6),
                                    direction="horizontal",title_position="topcenter",
                                    title_gp=gpar(fontsize=12,fontface="bold"),labels_gp=gpar(fontsize=12),
                                    grid_height=unit(0.6,"cm"),legend_width=unit(4,"cm"))
)

pdf(paste0("./figure/sigDiff_fig3D.2.pdf"),h=12,w=10)
draw(p,heatmap_legend_side="bottom",annotation_legend_side="right",
     row_title_side = "right", row_title = "Metastasis",row_title_gp = gpar(fontsize = 18),
     column_title_side = "bottom", column_title = "Tumorigenesis",column_title_gp = gpar(fontsize = 18))
dev.off()
#



# diff TCGA_small16 stat: log2FC cor scatter --------------------
diffFCcor.list <- list()
outSurfix="_rmPurity"

for(cellID in cellShort){
  message(cellID)

  disease.label <- "Yes"
  normal.label <- "No"
  diff1 <- read.table(paste0("./output/",dst,"/diff_",cellID,"_",disease.label,"vs",normal.label,outSurfix,".diff"), header = TRUE, sep="\t",check.names = FALSE)
  
  disease.label <- "Primary"
  normal.label <- "Normal"
  diff2 <- read.table(paste0("./output/",dst,"/diff_",cellID,"_",disease.label,"vs",normal.label,outSurfix,".diff"), header = TRUE, sep="\t",check.names = FALSE)
  
  all(rownames(diff1)==rownames(diff2))
  diffFCcor.list[[cellID]] <- data.frame("metastasis_logFC"=diff1$log2FoldChange,"tumorigenesis_logFC"=diff2$log2FoldChange,
                                         "metastasis.padj"=diff1$padj,"tumorigenesis.padj"=diff2$padj,
                                         cell=cellID,"feature"=rownames(diff1))
}
all(rownames(diff1)==rownames(diff2))
diffFCcor.df <- as.data.frame(do.call(rbind,diffFCcor.list))
diffFCcor.df$cell <- factor(diffFCcor.df$cell, levels = cellShort)

#str(diffFCcor.df)
p <- 
  ggplot(data=diffFCcor.df,aes(x=metastasis_logFC, y=tumorigenesis_logFC))+
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
ggsave(plot = p, filename = paste0("./figure/","logFCcor_all",".png"),device = "png", width = 12,height = 10, dpi = 200 ) # , limitsize = FALSE

p <- ggplot(data=diffFCcor.df,aes(x=metastasis_logFC, y=tumorigenesis_logFC))+ # [diffFCcor.df$cell==cellID,]
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
ggsave(plot = p, filename = paste0("./figure/","logFCcor_byCell",".png"), device = "png", width = 20,height = 20, dpi = 200, limitsize = FALSE)
#


# for visualization logFC lim =5,-5
dim(diffFCcor.df)
this.cells2 <- c("STAD","COAD","LIHC","LUAD","BRCA") #
res_sum_max5 = diffFCcor.df[diffFCcor.df$cell %in% this.cells2,]
colnames(res_sum_max5)[1:2] <- c("metastasis_logFC","tumorigenesis_logFC")
log2FC.cutoff <- 10
res_sum_max5$tumorigenesis_logFC[res_sum_max5$tumorigenesis_logFC> log2FC.cutoff] = log2FC.cutoff
res_sum_max5$tumorigenesis_logFC[res_sum_max5$tumorigenesis_logFC< -log2FC.cutoff] = -log2FC.cutoff
res_sum_max5$metastasis_logFC[res_sum_max5$metastasis_logFC> log2FC.cutoff] = log2FC.cutoff
res_sum_max5$metastasis_logFC[res_sum_max5$metastasis_logFC< -log2FC.cutoff] = -log2FC.cutoff
padj.cutoff <- 0.05
res_sum_max5$DEGtype <- "NA"
res_sum_max5$DEGtype[(res_sum_max5$metastasis.padj > padj.cutoff) & (res_sum_max5$tumorigenesis.padj > padj.cutoff)] <- "Not significant"
res_sum_max5$DEGtype[(res_sum_max5$metastasis.padj <= padj.cutoff) & (res_sum_max5$tumorigenesis.padj > padj.cutoff)] <- "Metastasis only"
res_sum_max5$DEGtype[(res_sum_max5$metastasis.padj > padj.cutoff) & (res_sum_max5$tumorigenesis.padj <= padj.cutoff)] <- "Tumorigenesis only"
res_sum_max5$DEGtype[(res_sum_max5$metastasis.padj <= padj.cutoff) & (res_sum_max5$tumorigenesis.padj <= padj.cutoff)] <- "Both significant"
res_sum_max5$DEGtype <- factor(res_sum_max5$DEGtype,levels = c("Not significant","Metastasis only","Tumorigenesis only","Both significant"))
res_sum_max5$cell <- factor(res_sum_max5$cell, levels = this.cells2)
res_sum_max5$siggroup005 <- res_sum_max5$DEGtype
res_sum_max5 = res_sum_max5[order(res_sum_max5$siggroup005,decreasing = FALSE),]
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
  labs(x=paste0("Tumorigenesis (logFC)"),y=paste0("Metastasis (logFC)"))
# p
ggsave(plot = p, filename = paste0(diffpre,"/figure/","logFCcor_byCell_withClass",".png"),height=6,width=10,device = "png")







# signature peak stats boxplot among diff signature grps --------------------
## boxplot by cell: compare signature geneset matrices like RBP, phastcons
diff.sum.df.sig2 <- data.table::fread( paste0(diffpre,"/output/",dst,"/diffSig.summary"),sep = "\t", header = TRUE, check.names=FALSE)
diff.sum.df.sig2$name <- unlist(sapply(strsplit(diff.sum.df.sig2$feature,"|",fixed=TRUE),"[",4))
df <- diff.sum.df.sig2[diff.sum.df.sig2$name %in% peak.table$name,]
df <- left_join(df,peak.table,by="name")

table(df$DEGtype)
df$DEGtype <- gsub("Both insignificant","Not significant",df$DEGtype)
df$DEGtype <- factor(df$DEGtype, levels = c("Not significant","Tumorigenesis only","Metastasis only","Both significant")) # 
table(df$DEGtype)
#table(df$cell)
df$cell <- factor(df$cell, levels = cellShort) # 
df <- df %>%
  group_by(cell,DEGtype) %>%
  mutate(DEGtypeNum=n_distinct(name))
# df$site_overlap_ratio_product_AGO2



## 1. RBP, G4, EV, Alu et al.
colnames(df)
select.columns <- colnames(df)[grepl("site_overlap_ratio_product|MFE_FDR",colnames(df),perl = TRUE)]
#colnames(df)
df2 <- as_tibble(df) %>% 
  dplyr::select(c(colnames(df)[1:7],select.columns)) %>% 
  pivot_longer(cols = select.columns, names_to = "site2",values_to = "site_overlap_ratio_product")
#table(df2$site2)
nrow(df2)==nrow(df)*length(unique(df2$site2))   # sanity check
df2$site2 <- gsub("site_overlap_ratio_product_","",df2$site2)
df2$site2[df2$site2=="MFE_FDR"] <- "structure"
df2 <- df2[df2$site2!="iM",]
table(df2$site2) # NA: Alu (peak.table)
df2$site2 <- factor(df2$site2,levels = c("AGO2","otherRBPs","RBPhotspot","EV","G4","structure","Alu_IR","Alu_nonIR"))
df2$site2 <- gsub("AGO2|otherRBPs|RBPhotspot","RBP",df2$site2,perl = TRUE)
#df2$site2[grepl("AGO2|otherRBPs|RBPhotspot",as.character(df2$site),perl = TRUE)] <- "RBP"
table(df2$site2)
table(is.na(df2$site2))
df2 <- df2[!is.na(df2$site2),]
df2$site2 <- factor(df2$site2, levels = c("RBP","EV","G4","structure","Alu_IR","Alu_nonIR"))
#table(df2$RBP_structure)
df3 <- as_tibble(df2) %>% 
  dplyr::mutate(RBP_structure=dplyr::case_when( site2=="RBP" & (site_overlap_ratio_product>=0.01) ~ TRUE,   #site2
                                                # site2=="otherRBPs" & (site_overlap_ratio_product>=0.01) ~ TRUE,
                                                # site2=="RBPhotspot" & (site_overlap_ratio_product>=0.01) ~ TRUE,
                                                site2=="EV" & (site_overlap_ratio_product>=0.01) ~ TRUE,
                                                site2=="G4" & (site_overlap_ratio_product>=0.01) ~ TRUE,
                                                site2=="Alu_IR" & (site_overlap_ratio_product>=0.01) ~ TRUE,
                                                site2=="Alu_nonIR" & (site_overlap_ratio_product>=0.01) ~ TRUE,
                                                site2=="structure" & (site_overlap_ratio_product<=0.1) ~ TRUE,
                                                TRUE ~ FALSE
  )
  ) %>% 
  dplyr::arrange(desc(site_overlap_ratio_product)) %>% #as.data.frame()
  dplyr::filter(!duplicated(cbind(feature, cell,name,site2))) # rm dup of df
#df3 <- df3[!duplicated(df3),]
data.table::fwrite(df3[,7:ncol(df3)],paste0(outpre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN/b5_d50_p1.RBP_structure.summary"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

#table(df3$site2)
for(cmp.column in c("RBP","EV","G4","structure","Alu_IR","Alu_nonIR")){
  # cmp.column <- "RBP" #"phastCons100way", #c("RBP","EV","G4","structure")
  #df$phastCons100way
  message(cmp.column)
  # table(df3$site2==cmp.column)
  df4 <- df3[df3$site2==cmp.column,]
  p =  ggplot(df4,aes(x=DEGtype,y=site_overlap_ratio_product,fill=DEGtype))+
    # geom_violin()+
    geom_boxplot(outlier.shape=NA )+ # , width=0.3/length(unique(df4$DEGtype)), alpha=0.8
    #geom_point(size=0.5,color="black",position=position_jitterdodge())+
    scale_fill_manual(values=DEGColor)+
    ylab(cmp.column)+
    facet_wrap(~cell,nrow=4) +
    # geom_text(aes(x=DEGtype,y=-10,label=format(DEGtypeNum,digits=3)),size=6,color="grey",angle=0)+
    ggpubr::stat_compare_means(comparisons = list(
      c("Not significant","Tumorigenesis only"),c("Not significant","Metastasis only"), c("Not significant","Both significant"), c("Metastasis only","Tumorigenesis only")),
      label.x.npc = 0.5,
      hide.ns=FALSE,size =8, vjust = 0.5, #step.increase = 0.06,#paired = TRUE,
      aes(label = p.signif), # p.signif, group = Group,
      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1,1),symbols = c( "****","***", "**", "*", "ns")),
      method = "wilcox.test") +
    theme_classic()+
    theme(
      aspect.ratio = 1,
      strip.text.x = element_text(size=20),
      strip.text.y = element_text(size=20),
      axis.text.x=element_text(colour="black",angle=90,hjust=1,vjust=0.5,size=24),
      axis.text.y=element_text(colour="black",size=24),
      axis.title.x=element_blank(),
      axis.title.y=element_text(colour="black",size=28),
      plot.title=element_blank(),
      strip.background=element_blank(),
      strip.text=element_blank(),
      panel.background =element_rect(fill="transparent",colour="black",size=1),
      panel.grid = element_blank(),
      legend.position="none")+
    scale_y_continuous(expand = expansion(mult = c(0, 0.30))) + #     # Add 33% space on the y-axis above the box plots
    # ylim(c(0,0.05)) +
    labs(title = cmp.column)
  # p
  ggsave(plot = p,filename = paste0("./figure/",dst,"_DEGtype_",cmp.column,".pdf"), device = "pdf", height = 14,width =14)
}

#




## 2. phastCons100way, CPC  et al.
select.columns <- colnames(df)[grepl("phastCons100way|coding_probability",colnames(df),perl = TRUE)]
df3 <- as_tibble(df) %>% 
  dplyr::select(c(colnames(df)[1:7],select.columns)) %>% 
  pivot_longer(cols = select.columns, names_to = "site2",values_to = "site_overlap_ratio_product")
df3 <- df3[!is.na(df3$site2),]

#table(df3$site2)
for(cmp.column in c("phastCons100way")){ # coding_probability
  # cmp.column <- "RBP" #"phastCons100way", #c("RBP","EV","G4","structure")
  message(cmp.column)
  #df$phastCons100way
  #df <- df[!is.na(df$MFE),]
  df4 <- df3[df3$site2==cmp.column,]
  p =  ggplot(df4,aes(x=DEGtype,y=site_overlap_ratio_product,fill=DEGtype))+
    # geom_violin()+
    geom_boxplot(outlier.shape=NA)+ # , width=0.3/length(unique(df4$DEGtype)), alpha=0.8 
    #geom_point(size=0.5,color="black",position=position_jitterdodge())+
    scale_fill_manual(values=DEGColor)+
    ylab(cmp.column)+
    facet_wrap(~cell,nrow=4) +
    # geom_text(aes(x=DEGtype,y=-10,label=format(DEGtypeNum,digits=3)),size=6,color="grey",angle=0)+
    ggpubr::stat_compare_means(comparisons = list(
      c("Not significant","Tumorigenesis only"),c("Not significant","Metastasis only"), c("Not significant","Both significant"), c("Metastasis only","Tumorigenesis only")),
      # label.x.npc = 0.5,
      hide.ns=FALSE,size =8, vjust = 0.5, #step.increase = 0.06,#paired = TRUE,
      aes(label = p.signif), # p.signif, group = Group,
      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1,1),symbols = c( "****","***", "**", "*", "ns")),
      method = "wilcox.test") +
    theme_classic()+
    theme(
      aspect.ratio = 1,
      strip.text.x = element_text(size=20),
      strip.text.y = element_text(size=20),
      axis.text.x=element_text(colour="black",angle=90,hjust=1,vjust=0.5,size=24),
      axis.text.y=element_text(colour="black",size=24),
      axis.title.x=element_blank(),
      axis.title.y=element_text(colour="black",size=28),
      plot.title=element_blank(),
      strip.background=element_blank(),
      strip.text=element_blank(),
      panel.background =element_rect(fill="transparent",colour="black",size=1),
      panel.grid = element_blank(),
      legend.position="none")+
    scale_y_continuous(expand = expansion(mult = c(0, 0.30))) #     # Add 33% space on the y-axis above the box plots
  # ylim(c(0,0.05)) +
  labs(title=cmp.column)
  # p
  ggsave(plot = p,filename = paste0("./figure/",dst,"_DEGtype_",cmp.column,".pdf"), device = "pdf", height = 14,width =14)
}
#




# signature peak candidates complexHeatmap ------------------------------------------------
RBP_structure <- data.table::fread(paste0(outpre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN/b5_d50_p1.RBP_structure.summary"),quote = FALSE,sep = "\t",data.table = FALSE, header = TRUE)
RBP_structure[1:3,]
table(duplicated(paste0(RBP_structure$name,RBP_structure$site2)))
RBP_structure <- RBP_structure[order(RBP_structure$name, RBP_structure$site2, as.numeric(RBP_structure$site_overlap_ratio_product), decreasing = TRUE),]
RBP_structure <- RBP_structure[!duplicated(paste0(RBP_structure$name,RBP_structure$site2)),]
RBP_structure <- as_tibble(RBP_structure[,c("name","site2","RBP_structure")]) %>% 
  pivot_wider(id_cols = "name", names_from = "site2", values_from = "RBP_structure")
colnames(RBP_structure)[2:ncol(RBP_structure)] <- paste0("logi_",colnames(RBP_structure)[2:ncol(RBP_structure)])

# append peak tab
peak.table.tmp <- left_join(peak.table, RBP_structure, by="name")

peak.table.tmp <- left_join(peak.table.tmp, RBP_structure[,])
rownames(peak.table.tmp) <- peak.table.tmp$name
peak.table.tmp$peakWid <- peak.table.tmp$end-peak.table.tmp$start
peak.table.tmp$gene <- ""
peak.table.tmp$gene[peak.table.tmp$tx_type %in% rna] <- peak.table.tmp$tx_symbol[peak.table.tmp$tx_type %in% rna] # mRNA... (not change to rna2 !)
peak.table.tmp$gene[peak.table.tmp$tx_type %in% dna] <- peak.table.tmp$nearestCodingTx_symbol[peak.table.tmp$tx_type %in% dna] # get precursor or nearest gene symbol
#table(peak.table.tmp$Inactive)

diff.sum.df <- data.table::fread(paste0(diffpre,"/output/",dst,"/diff.summary"), sep = "\t",check.names = FALSE,header = TRUE,data.table = FALSE)

# change cell tab
rownames(cell.table) <- cell.table$cellShort

#cmp.df
activity.lab <- "Yes"
state.lab <- "Primary"
activity.baseMean.mat <- diff.sum.df[,grepl(paste0(activity.lab,"\\.baseMean|feature"),colnames(diff.sum.df),perl = TRUE)]
colnames(activity.baseMean.mat) <- gsub(paste0("|",activity.lab,".baseMean"),"",colnames(activity.baseMean.mat),fixed = TRUE)
activity.baseMean.mat$feature <- NULL
state.baseMean.mat <- diff.sum.df[,grepl(paste0(state.lab,"\\.baseMean|feature"),colnames(diff.sum.df),perl = TRUE)]
colnames(state.baseMean.mat) <- gsub(paste0("|",state.lab,".baseMean"),"",colnames(state.baseMean.mat),fixed = TRUE)
state.baseMean.mat$feature <- NULL
dim(activity.baseMean.mat)
#max(activity.baseMean.mat[,2])
state.baseMean.mat <- log2(state.baseMean.mat+1)
activity.baseMean.mat <- log2(activity.baseMean.mat+1)

diff.sum.df.sig2 <- data.table::fread( paste0(diffpre,"/output/",dst,"/diffSig.summary"),sep = "\t", header = TRUE, check.names=FALSE)
diff.sum.df.sig2$name <- unlist(sapply(strsplit(diff.sum.df.sig2$feature,"|",fixed=TRUE),"[",4))
diff.sum.df.sig2$DEGtype <- gsub("Both insignificant","Not significant",diff.sum.df.sig2$DEGtype)
#str(diff.sum.df.sig2)

#DEGtype.mat <- pivot_wider(diff.sum.sig2, id_cols = "feature", names_from = )
DEGtype.mat <- list()
{ # cellNum
  diff.sum.df.sig3 <- diff.sum.df.sig2 %>%
    # mutate(record=1) %>% 
    tidyr::pivot_wider( id_cols = c("feature"), names_from = cell, values_from = DEGtype)
  
  tmp.ncol <- ncol(diff.sum.df.sig3)
  for(i in unique(diff.sum.df.sig2$DEGtype)){
    # i <- "Both insignificant"
    message(i)
    diff.sum.df.sig3[[paste0(i," cellNum")]] <- rowSums(diff.sum.df.sig3[,2:tmp.ncol]==i)
  }
  DEGtype.mat[["Both significant"]] <- as.data.frame(row.names = diff.sum.df.sig3$feature, diff.sum.df.sig3[,2:tmp.ncol]=="Both significant")
  DEGtype.mat[["Not significant"]] <- as.data.frame(row.names = diff.sum.df.sig3$feature, diff.sum.df.sig3[,2:tmp.ncol]=="Not significant")
  DEGtype.mat[["Metastasis only"]] <- as.data.frame(row.names = diff.sum.df.sig3$feature, diff.sum.df.sig3[,2:tmp.ncol]=="Metastasis only")
  DEGtype.mat[["Tumorigenesis only"]] <- as.data.frame(row.names = diff.sum.df.sig3$feature, diff.sum.df.sig3[,2:tmp.ncol]=="Tumorigenesis only")
}
table(diff.sum.df.sig3$ESCA)
tmp <- (as.matrix(DEGtype.mat[["Both significant"]]) + as.matrix(DEGtype.mat[["Not significant"]]) + as.matrix(DEGtype.mat[["Metastasis only"]]) + as.matrix(DEGtype.mat[["Tumorigenesis only"]]))
table(tmp)
dim(tmp)

DEGtype.mat[["Both significant"]] <- DEGtype.mat[["Both significant"]]*(state.baseMean.mat+activity.baseMean.mat)*0.5 # log2 CPM mean
DEGtype.mat[["Not significant"]] <- DEGtype.mat[["Not significant"]]*(state.baseMean.mat+activity.baseMean.mat)*0.5
DEGtype.mat[["Metastasis only"]] <- DEGtype.mat[["Metastasis only"]]*(activity.baseMean.mat)
DEGtype.mat[["Tumorigenesis only"]] <- DEGtype.mat[["Tumorigenesis only"]]*(state.baseMean.mat)
DEGtype.mat <- (as.matrix(DEGtype.mat[["Both significant"]]) + as.matrix(DEGtype.mat[["Not significant"]]) + as.matrix(DEGtype.mat[["Metastasis only"]]) + as.matrix(DEGtype.mat[["Tumorigenesis only"]]))
tumor.freq <- 3
selecSig <- rowSums(diff.sum.df.sig3[,2:tmp.ncol]!="Not significant")>=tumor.freq
table(selecSig)
# FALSE  TRUE 
# 19826 34528 
#DEGtype.mat <- DEGtype.mat[selecSig,]

rownames(DEGtype.mat) <- unlist(sapply(strsplit(rownames(DEGtype.mat),"|",fixed = TRUE),"[",4))
max(DEGtype.mat[,3])
dim(DEGtype.mat)
#DEGtype.mat.peakID <- rownames(DEGtype.mat) 
#table(duplicated(DEGtype.mat.peakID)) # all False 


#op: intersect RBP,G4 (newest,2408)
{
  # selectPeak <- rownames(DEGtype.mat) %in% peak.table.tmp$name[peak.table.tmp$site_overlap_ratio_product_G4>=0.6 & (peak.table.tmp$site_overlap_ratio_product_RBPhotspot>=0.8 | peak.table.tmp$site_overlap_ratio_product_otherRBPs>=0.8 | peak.table.tmp$site_overlap_ratio_product_AGO2>=0.6)] & selecSig # peak.table.tmp$site_overlap_ratio_product_EV>=0.5 |
  selectPeak <- rownames(DEGtype.mat) %in% peak.table.tmp$name[(peak.table.tmp$site_overlap_ratio_product_G4>=0.3 | peak.table.tmp$MFE_FDR<=0.1 | peak.table.tmp$site_overlap_ratio_product_EV>=0.3  ) & 
                                                                 (peak.table.tmp$site_overlap_ratio_product_RBPhotspot>=0.5 & peak.table.tmp$site_overlap_ratio_product_otherRBPs>=0.5 & peak.table.tmp$site_overlap_ratio_product_AGO2>=0.2) | 
                                                                 (peak.table.tmp$name=="peak_17765") # white list: TERC
  ] & selecSig
  #selectPeak <- selectPeak  
  table(selectPeak) # ? & 5003 --> 23
  figLab <- "Sig_RBP-EV-G4"
  
  
  
  # bubble plot
  mat <- DEGtype.mat[selectPeak,] #  log2(DEGtype.mat+1.0)[selectPeak,] # 
  #mat <- scale(t(scale(t(DEGtype.mat))))[selectPeak,] #
  mat2 <- as.matrix(diff.sum.df.sig3[selectPeak,2:tmp.ncol])
  rownames(mat2) <- rownames(mat)
  mat2.df <- as.data.frame(mat2)
  mat2.df$peak <- rownames(mat2.df)
  mat2.df <- mat2.df %>% 
    pivot_longer(cols = 1:(ncol(mat2.df)-1), names_to = "cell", values_to = "DEGtype")
  mat2.df$id <- paste0(mat2.df$peak, mat2.df$cell)
  mat2.df$peak <- NULL
  mat2.df$cell <- NULL
  
  mat.df <- as.data.frame(mat)
  mat.df$peak <- rownames(mat.df)
  mat.df <- mat.df %>% 
    pivot_longer(cols = 1:(ncol(mat.df)-1), names_to = "cell", values_to = "value")
  mat.df$id <- paste0(mat.df$peak, mat.df$cell)
  
  
  
  ## get anno
  peakAnno <- peak.table.tmp[rownames(mat),]
  peakAnno <- peakAnno[,c("gn_type","peakWid",colnames(peakAnno)[grepl("logi_",colnames(peakAnno))],"MFE_FDR","coding_probability","phastCons100way","gene")]
  colnames(peakAnno) <- gsub("logi_","",colnames(peakAnno))
  #colnames(peakAnno)
  peakAnno <- peakAnno %>% 
    mutate(gn_type = fct_relevel(gn_type, c(rna2,dna)[c(rna2,dna) %in% unique(as.character(peakAnno$gn_type))]  ) ) %>% 
    # arrange(desc(nchar(gene))) %>% # rev order of gene symbol len
    arrange(desc(fct_rev(gn_type))) %>%  # order by precursor
    mutate(gene_fct = fct_relevel(gene, (unique(gene))) )
  table(peakAnno$gn_type)

  cellAnno <- cell.table[colnames(mat),]
  #table(cellAnno$lineage)
  
  df <- left_join(mat.df,mat2.df,by="id")
  df$gene <- peakAnno[df$peak,"gene"]
  df <- df %>% 
    mutate(DEGtype_fct = fct_relevel(DEGtype, DEGfactor),
           cell_fct = fct_relevel(cell, cellShort),
           peak_fct = fct_relevel(peak, unique(rownames(peakAnno)))
    ) %>% 
    # mutate_at(vars(x, y), factor)
    arrange(desc(fct_rev(DEGtype_fct)))# %>% # rev
  #str(df)
  summary(df$value)
  minDotCPM <- 1
  df$value[df$DEGtype=="Not significant" & df$value<=minDotCPM] <- 0
  df <- df[order(df$peak_fct,decreasing = FALSE),]
  #table(df$cell_fct)
  pdf(paste0("./figure/",figLab,"_peak_heat.pdf"),width = 12, height =9)
  ggplot(data = df, aes(x=cell_fct, y=peak_fct, size=value,color=DEGtype_fct)) +
    geom_point( alpha=0.9) + # shape=21,color="grey30",stroke=0,
    # xlab(paste("")) +
    # ylab(paste("PC1(",pca.var.per[1],"%","Var)",sep="")) +
    # scale_y_discrete() +
    scale_color_manual(values = DEGColor, name="DEGtype") + 
    scale_size(#breaks = , labels = , limits = , 
      range = c(0, 8), name="log2CPM") +
    # scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A") +
    theme_minimal()+
    theme(
      aspect.ratio = 1.2,
      strip.text.x = element_text(size=20),
      strip.text.y = element_text(size=20),
      axis.text.x=element_text(colour="black",angle=90,hjust=1,vjust=0.5,size=16),
      axis.text.y=element_text(colour="black",size=16),
      axis.title.x=element_blank(),
      axis.title.y=element_text(colour="black",size=20),
      plot.title=element_blank(),
      strip.background=element_blank(),
      strip.text=element_blank(),
      # panel.background =element_rect(fill="transparent",colour="black",size=1),
      panel.grid = element_line(),
      legend.position="left",
      legend.text = element_text(colour="black",size=12),
      legend.title = element_text(colour="black",size=16)
    ) + 
    guides(color = guide_legend(override.aes = list(size=14))
           # size = guide_legend(override.aes = list(size=16))
    ) 
  dev.off()
  #
  

  peakAnno <- peakAnno %>% 
    rename_with(~c("Precursor","Width","2nd Structure","RBP","G4","EV","nonIR Alu","IR Alu","2nd Structure FDR","Coding prob.","PhastCons","gene","Gene"), 
                c("gn_type","peakWid","structure","RBP","G4","EV","Alu_nonIR","Alu_IR","MFE_FDR","coding_probability","phastCons100way","gene","gene_fct")) %>% 
    select(-"gene")
  peakAnno <- peakAnno[unique(df$peak),]
  ht_opt$ROW_ANNO_PADDING = unit(0.5, "cm")
  maxWid <- 180
  rowAnno <- HeatmapAnnotation(
    # "PeakID" = anno_mark(at = c(1:4), labels = rownames(peakAnno)[1:4]),
    df = peakAnno[,c("Precursor","Gene","2nd Structure","RBP","G4","EV","nonIR Alu","IR Alu")],
    
    Width = anno_barplot(peakAnno$Width,border = FALSE, baseline = 0, ylim = c(10,maxWid),width = unit(14, "mm"), 
                         gp = gpar(fill="grey30",col="white",fontsize=16),
                         axis_param=list(gp=gpar(fontsize = 8))),
    "\nCoding prob." = anno_barplot(peakAnno$`Coding prob.`,border = FALSE, baseline = 0., ylim = c(0,1),width = unit(22, "mm"), 
                                    gp = gpar(fill="firebrick",col="white",fontsize=16),
                                    axis_param=list(gp=gpar(fontsize = 8))),
    "PhastCons" = anno_barplot(peakAnno$PhastCons,border = FALSE, baseline = 0., ylim = c(0,1),width = unit(22, "mm"), 
                               gp = gpar(fill="seagreen",col="white",fontsize=16),
                               axis_param=list(gp=gpar(fontsize = 8))),
    "\n2nd Structure FDR" = anno_barplot( peakAnno$`2nd Structure FDR`, border = FALSE, width = unit(22, "mm"), 
                                          gp = gpar(fill="steelblue",col="white",fontsize=16),
                                          axis_param=list(gp=gpar(fontsize = 8))),
    
    col = list("Width" = circlize::colorRamp2(c(0, maxWid), c("white", "seagreen4")),
               "Precursor" = RNA_colors,
               "2nd Structure" = c("FALSE"="white","TRUE"="grey30"),
               "2nd Structure FDR" = circlize::colorRamp2(c(0, 1), c("firebrick", "white")),
               "nonIR Alu" = c("FALSE"="white","TRUE"="grey30"),
               "IR Alu" = c("FALSE"="white","TRUE"="grey30"),
               # "Alu" = c("FALSE"="white","TRUE"="grey30"),
               "EV" = c("FALSE"="white","TRUE"="grey30"),
               "G4" = c("FALSE"="white","TRUE"="grey30"),
               "RBP" = c("FALSE"="white","TRUE"="grey30")
               # "Gene" = tmp.color 
    ),
    na_col = "white",
    annotation_name_gp = gpar(size=0),simple_anno_size = unit(0.5, "cm"),
    show_annotation_name = TRUE,show_legend = TRUE,name = "",
    which = "row"
  )
  
  #draw(columnAnno)
  ht_opt$COLUMN_ANNO_PADDING = unit(0.5, "cm")
  columnAnno <- HeatmapAnnotation(
    "Lineage" = cellAnno$lineage,
    # "Cell" = anno_simple(cellAnno$cellShort, pch = rep("*",length(cellAnno$cellShort))),
    col = list(
      "Lineage" = lineageColor
    ),
    na_col = "white",
    annotation_name_gp = gpar(size=0), simple_anno_size = unit(0.5, "cm"),annotation_name_side = "left", # annotation_name_rot
    
    show_annotation_name = TRUE,show_legend = TRUE,name = "", #
    which = "col"
  )
  #str(cellAnno$lineage)
  tmp.mat <- matrix(0, nrow = nrow(mat),ncol=ncol(mat))
  rownames(tmp.mat) <- rownames(peakAnno)
  colnames(tmp.mat) <- cellAnno$cellShort
  ht = Heatmap(tmp.mat, name = "mat", 
               heatmap_legend_param = list(
                 title = "", at = c( 0),
                 labels = c(""),
                 color_bar = "discrete"
               ),
               col = circlize::colorRamp2(c(-5, 0, 5), c("steelblue4", "white", "firebrick")),
               na_col = "white", 
               width = unit(160, "mm"),
               
               # top_annotation = columnAnno,
               bottom_annotation = columnAnno,
               right_annotation = rowAnno,
               
               use_raster = TRUE, raster_quality = 8, # small size plot
               
               # Ensure rows are not clustered and follow the order of peakAnno
               cluster_rows = FALSE,
               row_order = rev(rownames(tmp.mat)),
               column_order = (colnames(tmp.mat))
  )
  pdf(paste0("./figure/",figLab,"_peak_heatRowAnno.pdf"),width = 16, height = 9)
  draw(object = ht)
  dev.off()
  #
}



# select sncRNA pan-cancer candidates (2408) ------
library(tidyr)
outSurfix <- "_rmPurity"

setwd(diffpre)
diff0.list <- list()
pval.cutoff1 <- 0.01
padj.cutoff1 <- 0.01
pval.cutoff2 <- 0.01 
padj.cutoff2 <- 1 # set to 0.1 not change much
fc.cutoff <- 0.6
outSurfix="_rmPurity" # _rmMonthOperator

diff0.list <- list()
for(cellID in cellShort){
  message(cellID)
  for(i in c("tumorigenesis","metastasis")){
    # i <- "tumorigenesis"
    message(i)
    
    disease.label <- cmp.df[i,"pos.grp"]
    normal.label <- cmp.df[i,"neg.grp"]
    col.grp <- cmp.df[i,"col"]
    
    x <- paste0("./output/",dst,"/diff_",cellID,"_",disease.label,"vs",normal.label,outSurfix,".diff")
    if(!file.exists(x)){
      message(x," not exist")
      next
    }
    #tmp <- read.table(x, header = TRUE, sep="\t",check.names = FALSE,stringsAsFactors = FALSE) # _extLibSize
    tmp <- data.table::fread(x, header = TRUE, sep="\t",check.names = FALSE,stringsAsFactors = FALSE, data.table = FALSE) # _extLibSize
    rownames(tmp) <- tmp$feature
    tmp$feature <- NULL
    tmp$diffType <- i
    tmp$disease <- disease.label
    tmp$normal <- normal.label
    tmp$cell <- cellID
    tmp$trend <- "insig"
    if(i=="tumorigenesis"){
      padj.cutoff <- padj.cutoff1
      pval.cutoff <- pval.cutoff1
    }else{
      padj.cutoff <- padj.cutoff2
      pval.cutoff <- pval.cutoff2
    }
    tmp$trend[tmp$padj<=padj.cutoff & tmp$pvalue<=pval.cutoff & (tmp$log2FoldChange > fc.cutoff)] <- "up"
    tmp$trend[tmp$padj<=padj.cutoff & tmp$pvalue<=pval.cutoff & (tmp$log2FoldChange < fc.cutoff)] <- "down"
    diff0.list[[paste0(cellID,"|",disease.label)]] <- tmp
  }
}

# #opt: rm rRNA
# tmp.id3 <- unlist(sapply(strsplit(rownames(tcga.diff),"|",fixed=TRUE),"[",2))!="rRNA"
# table(tmp.id3)
#   diff0.list[[i]] <- diff0.list[[i]][tmp.id3,]
# }
pos1.lab <- "Primary"
pos2.lab <- "Yes"
for (i in names(diff0.list)){
  # i <- "STAD|Primary"
  message(i)
  for(topN in  c(200) ){ # 
    # topN <- 200
    # message(topN)
    dim(diff0.list[[i]])
    tmp.diff <- diff0.list[[i]]
    tmp.diff <- filterDiffTabByTopN(tmp.diff, N = topN, rankBy = "pvalue")
    # tmp.diff <- filterDiffTabByP(diff0.list[[i]],pvalue = as.numeric(paste0("0.",topN)),padj = as.numeric(paste0("0.",topN)), abs.log2FC = 0.7)
    diff0.list[[i]][[paste0("top",topN)]] <- 0
    diff0.list[[i]][[paste0("top",topN)]][rownames(diff0.list[[i]]) %in% rownames(tmp.diff)] <- 1
  }
  for(topN in  c("up","down") ){ # 
    tmp <- diff0.list[[i]]
    
    if(unlist(sapply(strsplit(i,"|",fixed=TRUE),"[",2))==pos1.lab){
      padj.cutoff <- padj.cutoff1
      pval.cutoff <- pval.cutoff1
    }else if(unlist(sapply(strsplit(i,"|",fixed=TRUE),"[",2))==pos2.lab){
      padj.cutoff <- padj.cutoff2
      pval.cutoff <- pval.cutoff2
    }
    tmp$trend <- "insig"
    if(topN == "up"){
      tmp$trend[tmp$padj<=padj.cutoff & tmp$pvalue<=pval.cutoff & (tmp$log2FoldChange > fc.cutoff)] <- topN
    }else if(topN == "down"){
      tmp$trend[tmp$padj<=padj.cutoff & tmp$pvalue<=pval.cutoff & (tmp$log2FoldChange < fc.cutoff)] <- topN
    }
    diff0.list[[i]][[paste0("sig",topN)]] <- 0
    diff0.list[[i]][[paste0("sig",topN)]][tmp$trend == topN] <- 1
  }
  # add orphan defination
  tmp.diff <- diff0.list[[i]]
  #summary(tmp.diff$negMeanCPM)
  minCpm <- max(2, quantile(tmp.diff$negMeanCPM, probs = 0.2) )
  #table((tmp.diff$pvalue <= 0.1) & (tmp.diff$log2FoldChange > 0.5) & (tmp.diff$posGT1RatioCount >= 0.1) & (tmp.diff$negCent90CPM <= minCpm))
  tmp.diff <- tmp.diff[(tmp.diff$padj <= 0.05) & (tmp.diff$log2FoldChange > 0.6) & (tmp.diff$posGT1RatioCount >= 0.1) & (tmp.diff$negCent90CPM <= minCpm),]  # this version not use ex libsize
  #rownames(tmp.diff)
  diff0.list[[i]][[paste0("DCP")]] <- 0
  diff0.list[[i]][[paste0("DCP")]][rownames(diff0.list[[i]]) %in% rownames(tmp.diff)] <-"1"
  #table(rownames(diff0.list[[i]]) %in% rownames(tmp.diff))
  message("DCP: ",sum(as.numeric(diff0.list[[i]][[paste0("DCP")]])))
  # add dat label
  # colnames(diff0.list[[i]]) <- paste0(i, "_", colnames(diff0.list[[i]])) 
}
consensus.diff <- as.data.frame(do.call(cbind,diff0.list))

# add sig sum col
for(topN in c("up","down") ){# c(20,50,200) c("10","05","01")
  # topN <- 50
  message(topN)
  tmp <- dplyr::as_tibble(consensus.diff[,grepl(paste0("trend$"),colnames(consensus.diff),perl = TRUE)]) 
  consensus.diff[[paste0("unionSig",topN)]] <- 0
  consensus.diff[[paste0("unionSig",topN)]] <- rowSums(tmp==topN)
}

# add top sum col
for(topN in c(200) ){# c(20,50,200) c("10","05","01")
  # topN <- 50
  message(topN)
  tmp <- dplyr::as_tibble(consensus.diff[,grepl(paste0("top",topN,"$"),colnames(consensus.diff),perl = TRUE)]) 
  tmp <- tmp %>% dplyr::mutate(across(, ~as.numeric(.)))
  # str(tmp) 
  # table(rowSums(tmp)>=1)
  consensus.diff[[paste0("unionTop",topN)]] <- 0
  consensus.diff[[paste0("unionTop",topN)]] <- rowSums(tmp)
}

consensus.diff[grepl("peak_17765",rownames(consensus.diff)),  grepl("union",colnames(consensus.diff),perl = TRUE)]


filCol.sigup <- "unionSigup"
covNum.sigup <- 9  # phastCon vs tumor freq cor: 2
table(consensus.diff[[filCol.sigup]]>=covNum.sigup) 
filCol.sigdown <- "unionSigdown"
covNum.sigdown <- 9 # phastCon vs tumor freq cor: 2
table(consensus.diff[[filCol.sigdown]]>=covNum.sigdown) 

select.columns <- "PeakID|log2FoldChange|pvalue|padj|baseMean|DCP|top|Top|Sig|sig"
#consensus.diff[grepl("peak_17765",rownames(consensus.diff)), grepl(paste0("union"),colnames(consensus.diff),perl = TRUE) ] # 9 sig up...
consensus.diff <- consensus.diff[consensus.diff[[filCol.sigup]]>=covNum.sigup | consensus.diff[[filCol.sigdown]]>=covNum.sigdown, 
                                 grepl(select.columns,colnames(consensus.diff),perl = TRUE)]
consensus.diff$feature <- rownames(consensus.diff)
consensus.diff <- consensus.diff[,c(ncol(consensus.diff),1:(ncol(consensus.diff)-1))]

consensus.diff <- consensus.diff[,grepl(select.columns, colnames(consensus.diff), perl = TRUE)]
consensus.diff$name <- unlist(sapply(strsplit(rownames(consensus.diff),"|",fixed = TRUE),"[",4))

# sort cols
consensus.diff <- consensus.diff[,order(sapply(colnames(consensus.diff),function(x) paste(rev(unlist(strsplit(x,NULL))),collapse="")))]

# join meta peak
consensus.diff <- dplyr::left_join(consensus.diff,peak.table,by="name")
data.table::fwrite(consensus.diff, paste0(SMB_DIR,"/",dst,"_candidate_diff_summary.txt.tmp"), sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
colnames(consensus.diff)
#



# add TCGA survival pval (test 2k smps)
consensus.diff <- data.table::fread(paste0(SMB_DIR,"/",dst,"_candidate_diff_summary.txt.tmp"), sep = '\t', header = TRUE, check.names = FALSE ,stringsAsFactors = FALSE, data.table = FALSE)
colnames(consensus.diff)
#consensus.diff$TCGA.Surv.pval <- print(getSurvPval(exprVec = as.numeric(mat[tmp.id == x,]), sample.table = sample.table.pair, time.col = "OS_time", event.col = "OS"))
# library(tidyverse)
library(janitor)
#install.packages ("tidyquant")
library(tidyquant)
library(patchwork)
library(survival)
library(survminer)

surv <- "OS"
dir.create(paste0(SMB_DIR,"/",dst,"_figures/",surv))
for(cellID in cellShort){
  # cellID <- "LUSC"
  message(cellID)
  #metaPre <- paste0(TCGA_DIR,"/data/TCGA-all_clinical/",cellID,"_clinical/survival%2F",cellID,"_survival.txt")
  metaPre <- paste0(TCGA_DIR,"/data/TCGA-all_clinical/",cellID,"_clinical/survival%2F",cellID,"_survival.txt")
  sample.table <- rio::import(metaPre)
  #sample.table$submitter_id <- gsub(".bam","",sample.table$submitter_id)
  if("xena_sample" %in% colnames(sample.table)){
    sample.table <- sample.table[,c("xena_sample",surv,paste0(surv,".time"))]
    colnames(sample.table) <- c("sample",surv,paste0(surv,"_time"))
    sample.table$patient_id <- substr(sample.table$sample,1,12) 
  }else{
    sample.table <- sample.table[,c("sample","_PATIENT",surv,paste0(surv,".time"))]
    colnames(sample.table) <- c("smp_id","patient_id",surv,paste0(surv,"_time"))
  }
  # #read mat
  mat <- logcpm0 # count[,2:ncol(count)]
  colnames(mat) <- substr(colnames(mat),1,12) 
  sample.table <- sample.table[sample.table$patient_id %in% colnames(mat),]
  sample.table.pair <- sample.table
  sample.table.pair <- sample.table.pair[!is.na(sample.table.pair[[surv]]) & !is.na(sample.table.pair[[paste0(surv,"_time")]]),]
  sample.table.pair[[paste0(surv,"_time")]] <- sample.table.pair[[paste0(surv,"_time")]]/30
  mat <- mat[,sample.table.pair$patient_id]
  tmp.id <- unlist(sapply(strsplit(rownames(mat),"|",fixed = TRUE),"[",4))
  dim(mat)
  
  consensus.diff[[paste0(cellID,".",surv,".pval")]] <- NaN
  rm(list=c("exprVec.cut","exprVec","sfit_tmp","res.suvr","sample.table.tmp","surv.pval"))
  for(pID in consensus.diff$name ){ 
    message(pID)
    tryCatch(
      {
        surv.pval <- getSurvPval(exprVec = as.numeric(mat[tmp.id == pID,]), sample.table = sample.table.pair)
        # print(surv.pval)
        consensus.diff[consensus.diff$name == pID,paste0(cellID,".",surv,".pval")] <-  getSurvPval(exprVec = as.numeric(mat[tmp.id == pID,]), sample.table = sample.table.pair)
        plotSurv(exprVec = as.numeric(mat[tmp.id == pID,]), sample.table = sample.table.pair,facet = "NULL",time.col = paste0(surv,"_time"), event.col = surv,palette="npg", outFile = paste0(SMB_DIR,"/",dst,"_figures/",surv,"/",pID,"_",surv,".pdf"),width = 5,height = 5)
        
      },
      error = function(e) {
        # stop("too many zero expr? skip")
        # surv.pval <- NaN
        # exprVec.cut=NaN
      }
    )
    # 
  }
}



#summary(consensus.diff[,paste0(cellID,".Surv.pval")] )
consensus.diff <- consensus.diff[,c("name",colnames(consensus.diff)[colnames(consensus.diff)!="name"])]
dir.create(paste0(diffpre, "/output/",dst,"/candidate/"),recursive = TRUE,showWarnings = FALSE)
data.table::fwrite(consensus.diff, paste0(diffpre, "/output/",dst,"/candidate/candidate_diff_summary.txt") , sep = "\t",row.names = FALSE,col.names = TRUE, quote = FALSE)
consensus.diff <- data.table::fread(paste0(diffpre, "/output/",dst,"/candidate/candidate_diff_summary.txt") , sep = "\t",check.names = FALSE, stringsAsFactors = FALSE,header = TRUE, data.table = FALSE)
#

## ORA
m_t <- read.table(file.path(REF_DIR, "geneset/PATH_ID_NAME_KEGGplusHallmarkReactomeImmReg.txt"),sep = "\t",header = TRUE, stringsAsFactors = FALSE) # ,row.names = 1
m_t2g <- m_t[,c("short_DESCRPTION","ensembl_gene_id")]
gene.top.up <- unique(consensus.diff$nearestCodingTx_ENSG[consensus.diff$unionSigup>=consensus.diff$unionSigdown])
gene.top.down <- unique(consensus.diff$nearestCodingTx_ENSG[consensus.diff$unionSigup<consensus.diff$unionSigdown])
gene.top.up <- gene.top.up[!is.na(gene.top.up)]
gene.top.down <- gene.top.down[!is.na(gene.top.down)]
tmp <- ora.kegg.offline(gene.top.up = gene.top.up, gene.top.down = gene.top.down, m_t = m_t2g )
tmp <- tmp[tmp$ID!="",]

df <- tmp
df$lab <- "all" #gsub(sufix,"",df$lab)
# df$lab <- factor(df$lab, levels = select.RNAs.lab) # mRNAmiRNA
table(df$trend)
df$trend <- factor(df$trend, levels = c("up","dw","all"))
# df <- df[df$trCend=="all",]

df.sig <- df[df$p.adjust<=0.1,]
df.sig$source <- m_t$source[match(df.sig$ID,m_t$short_DESCRPTION)]
table(df.sig$source)
#df.sig <- df.sig[!(df.sig$source %in% c("Reactome")),]
if(nrow(df.sig)<=1){
  message("not enough, skip")
  next
}

df.sig$nchr <- nchar(as.character(df.sig$ID))
max.nchr <- 65
df.sig$ID[df.sig$nchr>max.nchr] <- paste0(substr(df.sig$ID[df.sig$nchr>max.nchr],1,max.nchr),"...")
# top <- seqs
df.sig.tmp <- df.sig #[df.sig$top == top,]
df.sig.tmp <- df.sig.tmp[order(df.sig.tmp$trend,df.sig.tmp$lab,-df.sig.tmp$nchr,decreasing = TRUE),] # ,df.sig.tmp$nchr
df.sig.tmp$ID <- factor(df.sig.tmp$ID, levels = unique(df.sig.tmp$ID))
enrichType <- "ORA"
if(enrichType=="GSEA"){
  df.sig.tmp$Count <- df.sig.tmp$NES #setSize
  sizeLab <- "abs(NES)"
}else if(enrichType=="ORA"){
  sizeLab = "Count"
}
df.sig.tmp$minus.log10.pvalue <- -log10(df.sig.tmp$pvalue)
ggplot(df.sig.tmp, aes(x=trend, y=ID) ) + # type2, type3
  geom_point(aes(color=minus.log10.pvalue, size=Count)) + # fill="grey",color="grey",  minus.log10.padj
  scale_color_continuous(name="-log10.P",type = "viridis") + # YlOrRd, viridis RColorBrewer::brewer.pal(n=8,name="Reds")
  scale_size_continuous(name = sizeLab) +
  ylab("") + 
  xlab("") + 
  theme_minimal() +  # base_size=12 void
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    # aspect.ratio = 10,
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 28, hjust = 1,vjust = 1, angle=45), #
    axis.text.y = element_text(size = 20, hjust = 1),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    # legend.position = "none", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
p_wid <- 4+max(df.sig.tmp$nchr)*0.1+length(unique(df.sig.tmp$lab))*0.3
p_het <- 4+length(unique(df.sig.tmp$ID))*0.5
ggsave(paste0("./figure/pancancer_kegg_dot_",enrichType,".pdf"), width = p_wid, height = p_het, limitsize = FALSE) # width = 15, height = 15
#


## plot venn/upset 
# colnames(consensus.diff)
pos1.lab <- "Primary"
pos2.lab <- "Yes"
for(pos.lab in c(pos1.lab,pos2.lab)){
  # pos.lab <- "Primary"
  message(pos.lab)
  
  for(topN in c("200")){ #c(20,50,200): 20 fail, due to no overlap?  c("10","05","01")
    # topN <- 200
    message(topN)
    pdf(paste0("./figure/peak_candidates_intersect_upset_top",topN,"_",pos.lab,".pdf"),width = 6,height = 4)
    
    setLab <- paste0(paste0(cellShort,"|",c(pos.lab),".top"),topN) # GSE62809miR.top
    tmp.df <- consensus.diff[,setLab] %>% mutate(across(, ~as.numeric(.)))
    colnames(tmp.df) <- unlist(sapply(strsplit(colnames(tmp.df) ,"|",fixed=TRUE),"[",1))
    setLab <- unlist(sapply(strsplit(setLab,"|",fixed=TRUE),"[",1))
    # str(tmp.df)
    p1 <- UpSetR::upset(data = tmp.df, 
                        sets = setLab, mb.ratio = c(0.4, 0.6),
                        nsets = length(cellShort), keep.order = TRUE,
                        matrix.color = "salmon",set_size.numbers_size = TRUE,
                        sets.bar.color = cellColor # colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(5)[1:4]
    )
    print(p1)
    dev.off()
  }
  
  for(topN in c("up","down")){ #c(20,50,200): 20 fail, due to no overlap?  c("10","05","01")
    # topN <- 200
    message(topN)
    pdf(paste0("./figure/peak_candidates_intersect_upset_sig",topN,"_",pos.lab,".pdf"),width = 6,height = 4)
    setLab <- paste0(paste0(cellShort,"|",c(pos.lab),".sig"),topN) # GSE62809miR.top
    tmp.df <- consensus.diff[,setLab] %>% mutate(across(, ~as.numeric(.)))
    colnames(tmp.df) <- unlist(sapply(strsplit(colnames(tmp.df) ,"|",fixed=TRUE),"[",1))
    setLab <- unlist(sapply(strsplit(setLab,"|",fixed=TRUE),"[",1))
    p1 <- UpSetR::upset(data = tmp.df, 
                        sets = setLab, mb.ratio = c(0.4, 0.6),
                        nsets = length(cellShort), keep.order = TRUE,
                        matrix.color = "salmon",set_size.numbers_size = TRUE,
                        sets.bar.color = cellColor # colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(5)[1:4]
    )
    print(p1)
    dev.off()
  }
}

dim(consensus.diff)
dat <- as.data.frame(table(consensus.diff$gn_type))
colnames(dat) <- c("RNA","ratio")
#str(dat)
dat$RNA <- factor(dat$RNA, levels=c(rna2,dna))
ggplot(dat,aes(x=RNA,y=ratio,fill=RNA))+
  geom_bar(position = "stack",width = 0.8,size=0.3, stat = "sum",color="black")+
  xlab("")+
  ylab(paste0("Pan-cancer DES number"))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  scale_fill_manual(name="Precursor",values = c("pri_miRNA"="#5E576FFF","other"="grey50",RNA_colors)) + # c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(9, "Set1"))
  theme_bw() +
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"),
        axis.text = element_text(size= 20,color = "black"),
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), #   , color = c(rep("steelblue",length(rna)),rep("firebrick",length(dna)))
        legend.position = "right",#c(.25,.6),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))
ggsave(paste0("./figure/",dst,"_pan_","_sig_num.png"),width = 8,height = 6,units = "in",dpi = 320)
#

cand.snc <- read.table(paste0(diffpre, "/output/",dst,"/candidate/candidate_diff_summary.txt"),header = TRUE,sep = "\t")
cand.snc <- cand.snc[cand.snc$name %in% c("peak_7567","peak_7155"),]
dim(cand.snc)
select.col <- c("TCGAmiR.log2FoldChange","TCGAmiR.pvalue","TCGAmiR.baseMean","name","gn","gn_start","gn_end","gn_name","gn_strand","tx_type","tx_len","tx_symbol","nearestCodingTx_symbol","MFE_FDR","MFE","MFE_structure","MFE_fasta","peak_fasta","phastCons100way","TCGA.Surv.pval")
cand.snc <- cand.snc[,select.col]
#




# peak logCPM boxplot among grps ----------
this.peak.ids <- c("peak_6200") 
this.cells <- c("LUAD", "LUSC", "HNSC", "LIHC")
candidate.class <- "-"
sample.table <- sample.table0 # 1042
sample.table$group <- sample.table$group1
table(sample.table$group)
sample.table$group <- factor(sample.table$group,levels = group3)
logcpm <- logcpm0[,sample.table$sample]
dim(logcpm)

library(dplyr)
tmp.feature <- splitPeakFeatureID(rownames(logcpm))

gene.mat <- as.data.frame(t(logcpm[tmp.feature$PeakID %in% this.peak.ids,]))
colnames(gene.mat) <- tmp.feature$PeakID[match(colnames(gene.mat),tmp.feature$feature)]
gene.mat.t <- gene.mat
gene.mat.t$sample <- rownames(gene.mat.t)
#gene.mat.t[1:4,]
gene.mat.t <- left_join(gene.mat.t,sample.table,by="sample")
dim(gene.mat.t)
## plot expr&pathway expr
#dir.create(paste0(diffpre,"/output/",dst,"/candidate/"),showWarnings = FALSE,recursive = TRUE)
for (g in this.peak.ids){
  # g <- "peak_50581"
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
    labs(title = g)+
    #ylab("ssGSEA score")+
    ylab("log2(CPM+1)")+
    xlab("") +
    ggpubr::stat_compare_means(
      comparisons = list( c("Normal","Primary"),c("Primary","Metastasis"),c("Normal","Metastasis") ),
      label.x.npc = 0.5,hide.ns=FALSE,size =8, vjust = 0.5, step.increase = 0.06,#paired = TRUE,
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
  ggsave(plot = p,filename = paste0("./figure/",dst,"_candidate_box_",g,".",candidate.class,".pdf"), device = "pdf", height = 6,width =8)
  # ggsave(plot = p,filename = paste0("./figure/",dst,"_candidate_box_",g,".",candidate.class,".pdf"), device = "pdf", height = 12,width =3)
}
#




# plot IGV: multi sample bigwig ------------
library(IRanges)
library(GenomicRanges)
# read candidate summary table

this.peak.ids <- c("peak_6200") #c(peak.union,peak.ImmuneRBP,peak.ExprSurv,"peak_17765")
this.cells <- c("LUAD", "LUSC", "HNSC", "LIHC")
candidate.ids <- this.peak.ids #c("peak_23","peak_25704","peak_28847")
sncCandidate <- peak.table[peak.table$name %in% candidate.ids, ]

bw.res.list <- list()
for(cellID in c(this.cells)){
  # this.cells
  #cellID <- "Naive.CD8.T" #cellShort
  sample.table <- sample.table0[sample.table0$cellShort==cellID,] # 1042
  #table(sample.table0$group,sample.table0$treatmentnaive) # seems all inactive was treated
  sample.table$group <- sample.table$group1
  table(sample.table$group1)
  
  # prepare meta.list
  meta.col <- list()
  genes.list <- list()
  annotation_row_merge <- sncCandidate[,c("tx","start","end","name","tx_len")]
  colnames(annotation_row_merge)[c(1,4)] <- c("chr","peak")
  annotation_row_merge$width <- annotation_row_merge$end-annotation_row_merge$start
  expFold <- 2
  
  annotation_row_merge$exp_start <- 0
  annotation_row_merge$exp_end <- annotation_row_merge$tx_len
  annotation_row_merge$exp_width <- 0
  for(i in 1:nrow(annotation_row_merge)){
    annotation_row_merge$exp_start[i] <- max(0,annotation_row_merge$start[i]-expFold*annotation_row_merge$width[i])
    annotation_row_merge$exp_end[i] <- min(annotation_row_merge$tx_len[i],annotation_row_merge$end[i]+expFold*annotation_row_merge$width[i])
    annotation_row_merge$exp_width[i] <- annotation_row_merge$exp_end[i]-annotation_row_merge$exp_start[i]
    genes.list[[i]] <- data.frame(chr=annotation_row_merge[i,"chr"],start=annotation_row_merge[i,"exp_start"],width=annotation_row_merge[i,"exp_width"],peak=annotation_row_merge[i,"peak"],peak.start=annotation_row_merge[i,"start"],peak.width=annotation_row_merge[i,"width"]) #mir
  }
  
  
  # prepare IGV meta table
  randNum <- 8 # 16, 5
  randFrac <- 1
  library(tidyverse)
  library(tidydr)
  #conflict_prefer("select", "dplyr")
  
  # dataset: >18 cell types
  dsts.list <- list()
  
  #smp.tbl <- read.table(paste0(pre,"/exSeek-dev/data/",dst,"/sample_table_202311.txt"),header = TRUE, sep = "\t",quote = "",check.names = FALSE,stringsAsFactors = FALSE)
  smp.tbl <- sample.table
  smp.tbl <- smp.tbl[!(is.na(smp.tbl$group)),]
  #colnames(smp.tbl)
  #smp.tbl[1:3,]
  #colnames(smp.tbl)[c(1,4,5)] <- c("sample","group")
  smp.tbl <- smp.tbl %>% 
    dplyr::select(c("sample","group")) %>% 
    dplyr::group_by(group) %>% 
    dplyr::sample_n(randNum)
    # dplyr::sample_frac(randFrac)
  smp.tbl$bwPath <- paste0(outpre,"/output/",dst,"/call_peak",dedup,"/tbigwig_RNA_EM/",smp.tbl$sample,".transcriptome.bigWig")
  smp.tbl$dataset <- cellID
  #table(smp.tbl$origin) # Buccal Mucosa Floor of mouth    Oral Cavity    Oral Tongue
  dsts.list[[dst]] <- smp.tbl
  
  
  dsts.df <- do.call(rbind,dsts.list)
  #dsts.df$dataset <- factor(dsts.df$dataset,levels = cellIDshort)
  dsts.df$group <- factor(dsts.df$group,levels = group3)
  dsts.df <- dsts.df[order(dsts.df$dataset,dsts.df$group),] # bw order !
  
  
  # add color each dst group
  dsts <- unique(as.character(dsts.df$dataset))
  #cols <- c("seagreen4","steelblue","salmon","salmon3","salmon4") # groupShortColors
  cols <- colorList[["group1"]] #groupColor
  #names(cols) <- groupShort
  dsts.df$color <- "grey"#FA80721A
  #j2 <- 1
  for (i in 1:length(dsts)){
    tmp.df <- dsts.df[dsts.df$dataset==dsts[i],]
    #tmp.df <- tmp.df[order(tmp.df$group),]
    grps <- unique(as.character(tmp.df$group))
    #print(grps)
    for (j in 1:length(grps)){
      dsts.df$color[dsts.df$dataset==dsts[i] & dsts.df$group==grps[j]] <- cols[grps[j]]
      # j2 <- j2+1
    }
  }
  
  
  for (i in 1:length(genes.list)){
    # i <- 1 #"ENST00000385214_____1_30_49_+|pri_miRNA|ENST00000385214_____1|peak_6181|ENST00000385214_____1|30|49"
    chr <- genes.list[[i]]$chr
    start <- as.numeric(genes.list[[i]]$start)
    width <- as.numeric(genes.list[[i]]$width)
    peak.start <- as.numeric(genes.list[[i]]$peak.start)
    peak.width <- as.numeric(genes.list[[i]]$peak.width)
    peakID <- genes.list[[i]]$peak
    print(paste0(chr,".",start,".",width))
    
    res.list <- list()
    for (dst.tmp in dsts){ # names(meta.list)
      # dst.tmp <- "BRCA"
      smps <- dsts.df$sample[dsts.df$dataset==dst.tmp]
      grps <- dsts.df$group[dsts.df$dataset==dst.tmp]
      
      for (j in 1:length(smps)){
        #j <- 1
        smp <- smps[j]
        grp <- grps[j]
        # print(smp)
        # print(grp)
        bw.color <- dsts.df$color[dsts.df$dataset==dst.tmp & dsts.df$sample==smp]
        plotXbreak <- ifelse(dst.tmp==dsts[length(dsts)] & smp==smps[length(smps)], TRUE, FALSE)
        dst.lab <- dst.tmp
        bw.res.df <- readBW(bw=dsts.df$bwPath[dsts.df$dataset==dst.tmp & dsts.df$sample==smp],chr=chr,start=start,width=width)

        title <- "" #paste0(chr,":",start,"-",start+width)#""
        res.list[[paste0(dst.lab,"_",smp)]] <- plotBW(single_base_bw=bw.res.df,chr=chr,start=start,width=width,highlight=c(peak.start,peak.start+peak.width),ylab="",title=title,plotXbreak = plotXbreak, plotYbreak = TRUE, plotYtext = FALSE, color = bw.color, fill = bw.color, libsize = externalLibSize[smp]) # , maxY = maxY.list[[dst.lab]] , maxY=5200, , 
        bw.res.df$name <- paste(cellID,i,dst.tmp,grp,smp,sep =  "|")
        bw.res.list[[paste(cellID,i,dst.tmp,grp,smp,sep =  "|")]] <- bw.res.df
      }
    }
    # ## collect
    # tmp <- ggpubr::ggarrange(plotlist = res.list, ncol = 1,  align = "hv") + #heights = c(7,1,1,1,1,1,1,1,1,1,1), nrow = length(res.list),
    #   ggtitle(label = paste0(c(dsts,unique(grps))))
    # ggsave(plot = tmp,filename = paste0("./figure/",dst,"_candidate_",peakID,"_",cellID,"_bw","-",chr,".",start,".",width,".png"),width = 15, height = 1.3*length(res.list)) # 15
  }
}
#




# plot meta cov cpm line plot with ribbon (single-base, merge smps)
#names(bw.res.list)
df <- as.data.frame(do.call(rbind,bw.res.list))
dim(df)
tmpL <- strsplit(df$name,"|",fixed = TRUE)
df$cell <- unlist(sapply(tmpL,"[",1))
df$txid <- unlist(sapply(tmpL,"[",2))
df$dst <- unlist(sapply(tmpL,"[",3))
df$grp <- unlist(sapply(tmpL,"[",4))
df$smp <- unlist(sapply(tmpL,"[",5))
all((df$end-df$start)==1) # sanity check !
all(df$dst==df$cell) # sanity check !

df <- dplyr::as_tibble(df) %>% 
  dplyr::group_by(seqnames,start,end,strand,txid,cell,grp) %>% 
  dplyr::summarise(mean = mean(score, na.rm=TRUE), #,trim = 0.05, na.rm=TRUE
                   sd = sd(score),
                   upper = Rmisc::CI(score,ci = 0.95)[1],
                   lower = Rmisc::CI(score,ci = 0.95)[3],
                   max = max(score),
                   min = min(score)
  ) %>% as.data.frame()

df$grp <- factor(df$grp,levels = group3)
#df$mean <- 
df$start <- df$start+1
for (i in sort(as.numeric(unique(df$txid)))){
  # i <- 1
  message(i)
  chr <- genes.list[[i]]$chr
  start <- as.numeric(genes.list[[i]]$start)
  width <- as.numeric(genes.list[[i]]$width)
  peak.start <- as.numeric(genes.list[[i]]$peak.start)
  peak.width <- as.numeric(genes.list[[i]]$peak.width)
  peakID <- genes.list[[i]]$peak
  print(paste0(chr,".",start,".",width))
  
  highlight=c(peak.start,peak.start+peak.width)
  df.tmp <- df[df$txid==i,]
  df.tmp$cell <- factor(df.tmp$cell, levels = this.cells)
  ggplot(df.tmp,aes(x = start,y = mean, group=grp)) + # , fill=grp, color=grp
    # geom_vline(xintercept = aes(x=)) +
    # add 0.95 interval
    geom_ribbon(aes(ymin = lower,
                    ymax = upper,
                    fill = grp), # sample
                alpha = 0.5) +
    geom_line(aes(color = grp), size = 2) + # 
    geom_vline(xintercept = highlight, linetype="dashed", color="grey70") +
    theme_classic(base_size = 16) +
    scale_color_manual(name="Group",values = groupAllShort) +
    # scale_color_d3(name = 'Data type') +
    scale_fill_manual(name="Group",values = groupAllShort) +
    # geom_vline(xintercept = c(up_bins+1,up_bins+target_bins),color="grey50",linetype="dashed")+
    # geom_vline(xintercept = c((1+up_bins+target_bins+down_bins)*0.5),color="grey50",linetype="dashed")+ # need +1 for 1-based coord
    # x label
    scale_x_continuous() +
    # scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.0), limits = c(0,1), labels = c("0","0.25","0.50","0.75","1.00")) +
    xlab('Precuror') + ylab('CPM') +
    theme(aspect.ratio = 0.6,
          strip.background = element_rect(color = NA,fill = 'white'),
          strip.text = element_text(size=24),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=24),
          axis.text.x = element_text(size=16,angle = 0,hjust = 0.5),
          axis.text.y = element_text(size=16),
          legend.position =  "right",#c(0.9,0.8),
          legend.text = element_text(size= 24),
          legend.title= element_text(size= 24) ) + 
    facet_grid(cell ~ grp, scales = "free_y")
  ggsave(paste0("./figure/",dst,"_candidate_",peakID,"_",cellID,"_bwLine","-",chr,".",start,".",width,".pdf"),width = 12,height = 8)
}
#




# plot Survival (separate) --------
library(tidyverse)
library(janitor)
#install.packages ("tidyquant")
library(tidyquant)
library(patchwork)
library(survival)
library(survminer)
this.peaks <- c("peak_6200")
surv <- "OS"
for(cellID in cellShort){
  # cellID <- "LUAD"
  message(cellID)
  dir.create(paste0(SMB_DIR,"/",dst,"_figures/",cellID,"_",surv))
  
  #metaPre <- paste0(TCGA_DIR,"/data/TCGA-all_clinical/",cellID,"_clinical/survival%2F",cellID,"_survival.txt")
  metaPre <- paste0(TCGA_DIR,"/data/TCGA-all_clinical/",cellID,"_clinical/survival%2F",cellID,"_survival.txt")
  sample.table <- rio::import(metaPre)
  #sample.table$submitter_id <- gsub(".bam","",sample.table$submitter_id)
  if("xena_sample" %in% colnames(sample.table)){
    sample.table <- sample.table[,c("xena_sample",surv,paste0(surv,".time"))]
    colnames(sample.table) <- c("sample",surv,paste0(surv,"_time"))
    sample.table$patient_id <- substr(sample.table$sample,1,12) 
  }else{
    sample.table <- sample.table[,c("sample","_PATIENT",surv,paste0(surv,".time"))]
    colnames(sample.table) <- c("smp_id","patient_id",surv,paste0(surv,"_time"))
  }
  # #read mat
  mat <- logcpm0 # count[,2:ncol(count)]
  colnames(mat) <- substr(colnames(mat),1,12) 
  sample.table <- sample.table[sample.table$patient_id %in% colnames(mat),]
  sample.table.pair <- sample.table
  sample.table.pair <- sample.table.pair[!is.na(sample.table.pair[[surv]]) & !is.na(sample.table.pair[[paste0(surv,"_time")]]),]
  sample.table.pair[[paste0(surv,"_time")]] <- sample.table.pair[[paste0(surv,"_time")]]/30
  mat <- mat[,sample.table.pair$patient_id]
  tmp.id <- unlist(sapply(strsplit(rownames(mat),"|",fixed = TRUE),"[",4))
  dim(mat)
  
  rm(list=c("exprVec.cut","exprVec","sfit_tmp","res.suvr","sample.table.tmp","surv.pval"))
  for(pID in this.peaks ){  
    message(pID)
    tryCatch(
      {
        surv.pval <- getSurvPval(exprVec = as.numeric(mat[tmp.id == pID,]), sample.table = sample.table.pair)
        plotSurv(exprVec = as.numeric(mat[tmp.id == pID,]), sample.table = sample.table.pair,facet = "NULL",time.col = paste0(surv,"_time"), event.col = surv,palette="npg", outFile = paste0(SMB_DIR,"/",dst,"_figures/",cellID,"_",surv,"/",pID,"_",surv,".pdf"),width = 5,height = 5)
        
      },
      error = function(e) {
        # stop("too many zero expr? skip")
        # surv.pval <- NaN
        # exprVec.cut=NaN
      }
    )
    # 
  }
}
#


# PCC (with clinical. eg: immuneScore, age) ------------------------
library(parallel)


# logcpm0.raw <- data.table::fread( paste0(outpre,"/output/",dst,"/call_peak",dedup,"/count_matrix/cfpeakCNN_b5_d50_p1_CPMtmm.txt"), sep = "\t",check.names = FALSE,header = TRUE,data.table = FALSE,stringsAsFactors = FALSE )
# rownames(logcpm0.raw) <- logcpm0.raw$gene_id
# #min(logcpm0.raw$`TCGA-BH-A0BJ-01A-11R-A057-13_mirna_gdc_realn`)
# logcpm0.raw$gene_id <- NULL
# logcpm0.raw <- log2(logcpm0.raw+1)
# dim(logcpm0.raw) # 44706   376

# logcpm0 need rm batch by limma or combat
dir.create(paste0("./output/",dst,"/pcc_clinical"),recursive = TRUE,showWarnings = FALSE)

#outSurfix <- "_rmMonthOperator"
#pcc_var <- "" #"age"
method="pearson" # "spearman"
dir.create(paste0("./output/",dst,"/pcc_clinical"),recursive = TRUE,showWarnings = FALSE)
#dir.create(paste0("./output/",dst,"/clinical/"),recursive = TRUE,showWarnings = FALSE)

for(cellID in cellShort){
  # cellID <- "HNSC"
  message(cellID)
  
  # prep smpTbl/clinical
  sample.table <- sample.table0[sample.table0$cell==cellID,] #[(sample.table0$Threatment=="-") & (sample.table0$treatmentnaive %in% c("-","0")),] # 1042
  sample.table <- sample.table[,c("age","gender","race","tumor_stage","ajcc_stage","ajcc_n","cigarettes_per_day","years_smoked","days_to_death","purity")] #tumor_stage
  sample.table$ajcc_stage <- tolower(sample.table$ajcc_stage)
  dim(sample.table)
  str(sample.table)
  
  ## numeric
  x <- sample.table$days_to_death # age, days_to_death
  summary(x)
  table(is.na(x))
  table(x=="")
  clinical.num <- sample.table[!is.na(sample.table$days_to_death),c("age","days_to_death")]
  clinical.num <- t(clinical.num)
  data.table::fwrite(x = cbind(feature=rownames(clinical.num),clinical.num), paste0("./output/",dst,"/pcc_clinical/",cellID,"_clinical_num.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )
  
  ## category
  x <- sample.table$race #race #gender tumor_stage
  table(x)
  table(is.na(x))
  table(x=="")
  clinical.category <- sample.table[!is.na(sample.table$ajcc_stage) & sample.table$race!="not reported", c("race","gender","ajcc_stage")]
  clinical.category <- t(clinical.category)
  data.table::fwrite(x = cbind(feature=rownames(clinical.category),clinical.category), paste0("./output/",dst,"/pcc_clinical/",cellID,"_clinical_category.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )
  
  # prep sncRNA (done in mRNA_sncRNA_PCC)
  #...
  
  # prep sncRNA PCC
  
}



# PCC of clinical variables: the cluster jobs are submitted by
# scripts/bash/02_pcc_clinical.sh
cellID <- "HNSC"
method <- "pearson"

clinical.cor <- data.table::fread(paste0("output/",dst,"/pcc_clinical/",cellID,"_clinical_num_",method,".txt"), sep = "\t", check.names = FALSE, header = TRUE, data.table = FALSE )
#




# PCC (with mRNA_TPM or CNV eg: TP53.Expr, SEPTIN9.CNV) ------------------------
## read gene reference
ref <- read.table(file.path(REF_DIR, "geneset/encodev27_gene_v2.txt"),sep = "\t",header = TRUE)
table(duplicated(ref$ensg.full))
table(ref$type)
#ref <- ref[ref$type!="miRNA",]
ref <- ref[ref$type=="protein_coding",]
table(duplicated(ref$ensg))
ref <- ref[!duplicated(ref$ensg),]

ref.mrna.lncrna <- ref[ref$type %in% c("protein_coding","lincRNA","antisense_RNA","macro_lncRNA","non_coding"),]
table(duplicated(ref.mrna.lncrna$ensg))
ref.mrna.lncrna <- ref.mrna.lncrna[!duplicated(ref.mrna.lncrna$ensg),]

library(parallel)
# read logcpm0.raw from util

#reduce.enh$V4

dir.create(paste0("./output/",dst,"/pcc_mRNA"), recursive = TRUE, showWarnings = FALSE)

method="pearson" # "spearman"
outSurfix=""
for(cellID in cellShort){
  # cellID <- "BLCA"
  message(cellID)
  
  sample.table <- sample.table0[sample.table0$cell==cellID & (sample.table0$sample %in% testID),] #[(sample.table0$Threatment=="-") & (sample.table0$treatmentnaive %in% c("-","0")),] # 1042
  sample.table <- sample.table[order(sample.table$group),]
  #sample.table <- sample.table[!is.na(sample.table$purity),]
  dim(sample.table)
  
  ## get long mRNA mat
  matchTPM <- paste0(TCGA_DIR,"/data/",cellID,"_long/log2TPM_filter.txt") # use not rm Purity
  if(!file.exists(matchTPM)){
    message("filter ",cellID)
    
    longTPM <-   data.table::fread( file = paste0(TCGA_DIR,"/data/",cellID,"_long/TPM.txt"), sep = "\t", check.names = FALSE, header = TRUE, data.table = FALSE )
    rownames(longTPM) <- longTPM$feature
    longTPM$feature <- NULL
    table(rownames(longTPM) %in% ref$ensg)
    longTPM <- longTPM[ref$ensg,]
    longTPM <- na.omit(longTPM)
    #max(longTPM$`TCGA-CV-7177-11A`)
    #table(sample.table$tissue_id %in% colnames(longTPM))
    #colnames(longTPM)[grepl("TCGA-A6-2674-01A",(colnames(longTPM)))]
    longTPM <- longTPM[,!duplicated(colnames(longTPM))]
    longTPM <- longTPM[,colnames(longTPM) %in% sample.table$tissue_id]
    colnames(longTPM) <- sample.table$sample[match(sample.table$tissue_id, colnames(longTPM) )]
    longTPM <- longTPM[,sample.table$sample]
    logTPM <- log2(longTPM+1)
    #summary(logTPM[,3])
    table(logTPM>=0)
    data.table::fwrite(x = cbind(feature=rownames(logTPM),logTPM), matchTPM, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )
  }else{
    message("skip ",cellID)
  }
  
  ## get long mRNA+lncRNA mat
  matchTPM2 <- paste0(TCGA_DIR,"/data/",cellID,"_long/log2TPM_filter_withLnc.txt") # use not rm Purity
  if(!file.exists(matchTPM2)){
    message("filter ",cellID)
    
    longTPM <-   data.table::fread( file = paste0(TCGA_DIR,"/data/",cellID,"_long/TPM.txt"), sep = "\t", check.names = FALSE, header = TRUE, data.table = FALSE )
    rownames(longTPM) <- longTPM$feature
    longTPM$feature <- NULL
    table(rownames(longTPM) %in% ref.mrna.lncrna$ensg)
    longTPM <- longTPM[ref.mrna.lncrna$ensg,]
    longTPM <- na.omit(longTPM)
    longTPM <- longTPM[,!duplicated(colnames(longTPM))]
    longTPM <- longTPM[,colnames(longTPM) %in% sample.table$tissue_id]
    colnames(longTPM) <- sample.table$sample[match(sample.table$tissue_id, colnames(longTPM) )]
    longTPM <- longTPM[,sample.table$sample]
    logTPM <- log2(longTPM+1)
    #summary(logTPM[,3])
    table(logTPM>=0)
    data.table::fwrite(x = cbind(feature=rownames(logTPM),logTPM), matchTPM2, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )
  }else{
    message("skip ",cellID)
  }
  ## get long diff
  #data.table::fwrite(cbind(feature=rownames(res[["diffTable"]]),res[["diffTable"]]),paste0(outDir,"/diff_",cellID,"_",pos.grp,"vs",neg.grp,"",outSurfix,".diff"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
  
  
  # read peak diff (keep topDiff or topVar gene)
  # disease.label <- "Primary"
  # normal.label <- "Normal"
  # outSurfix <- ""
  # diff <- read.table(paste0("./output/",dst,"/diff_",cellID,"_",disease.label,"vs",normal.label,outSurfix,".diff"), header = TRUE, sep="\t",check.names = FALSE,stringsAsFactors = FALSE, row.names = 1) # _extLibSize
  # diff <- diff[order(diff$pvalue,-abs(diff$log2FoldChange),decreasing = FALSE),]
  # diff <- diff[1:as.integer(0.05*nrow(diff)),] # use top 10
  #
  
  # read small cpm
  logCPM <- logcpm0.raw[,sample.table$sample] # logcpm0.raw need rm batch by limma or combat
  logCPM <- logCPM[rowSums(as.matrix(logCPM))>0,]
  tmp.sd <- as.numeric(apply(as.matrix(logCPM),1,sd))
  logCPM <- logCPM[order(tmp.sd,decreasing = TRUE),]
  
  
  # ## opt1: filter top var sncRNA
  # topVar <- rownames(logCPM)[1:min(5000,as.integer(nrow(logCPM)*0.1))]
  # #logCPM <- logCPM[tmp.sd>0,]
  # # table((rownames(logCPM) %in% topVar ))
  # # table((rownames(logCPM) %in% rownames(diff)))
  # # table( (rownames(logCPM) %in% rownames(diff)) | (rownames(logCPM) %in% topVar ))
  # table(tmp.sd>0)
  # # logCPM <- logCPM[ ( (rownames(logCPM) %in% rownames(diff)) | (rownames(logCPM) %in% topVar ) ) & tmp.sd>0,]
  # logCPM <- logCPM[ ( (rownames(logCPM) %in% topVar ) ) & tmp.sd>0,]
  # # all(colnames(logCPM)==colnames(logTPM))
  # message(nrow(logCPM))
  # data.table::fwrite(x = cbind(feature=rownames(logCPM),logCPM), paste0("./output/",dst,"/",cellID,"_sncRNA_log2CPM_topVar",outSurfix,".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )
  
  
  ## opt2: filter interested sncRNA
  goi_snc <- c("peak_17765","peak_199")
  pid.split <- splitPeakFeatureID(rownames(logCPM))
  logCPM <- logCPM[ ( (pid.split$PeakID %in% goi_snc ) ) & tmp.sd>0,]
  # all(colnames(logCPM)==colnames(logTPM))
  message(nrow(logCPM))
  data.table::fwrite(x = cbind(feature=rownames(logCPM),logCPM), paste0("./output/",dst,"/",cellID,"_sncRNA_log2CPM_select",outSurfix,".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )
  
  # # read small eRNA cpm
  # logCPM <- logcpm0[,sample.table$sample]
  # logCPM <- logCPM[rowSums(as.matrix(logCPM))>0,]
  # tmp.sd <- as.numeric(apply(as.matrix(logCPM),1,sd))
  # logCPM <- logCPM[order(tmp.sd,decreasing = TRUE),]
  # 
  # table(tmp.sd>0)
  # logCPM <- logCPM[ tmp.sd>0 , ] # (rownames(logCPM) %in% topVar )  & 
  # logCPM <- logCPM[grepl("enhancer_",rownames(logCPM)),]
  # logCPM <- logCPM[(unlist(sapply(strsplit(rownames(logCPM),"|",fixed = TRUE),"[",3)) %in% unique(reduce.enh$V4)),]
  # message(nrow(logCPM))
  # data.table::fwrite(x = cbind(feature=rownames(logCPM),logCPM), paste0("./output/",dst,"/",cellID,"_sncRNA_log2CPM_enhancer",outSurfix,".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )
}
# # raw PCC
# res.list <- parallel::mclapply(X=1:nrow(logCPM),FUN = getCor,  pcc_var=pcc_var, method=method,conf.level=0.95,pcc_type="original",covar="Null", mc.cores = 10)
# res.df <- as.data.frame(do.call(rbind,res.list))
# res.df <- res.df[order(res.df$id,decreasing = FALSE),]
# res.df[["feature"]] <- rownames(logCPM)
# res.df[["cell"]] <- cellID
# data.table::fwrite( res.df, paste0("./output/",dst,"/pcc/",pcc_var,"","_",cellID,"_",method,outSurfix,".txt"), row.names = FALSE,col.names = TRUE,quote = FALSE,sep = "\t")






# PCC against mRNA expression: the cluster jobs are submitted by
# scripts/bash/03_pcc_mrna.sh
method    <- "pearson"
outSurfix <- "_rmPurity"

# Pathway GSEA on the PCC tables: the cluster jobs are submitted by
# scripts/bash/04_snc_pcc_gsea_gmt.sh
# (output names there carry no leading underscore).





# longRNA-seq Immune Cell Frac. (Cibersortx)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
conflicted::conflict_prefer("which", "base")
source(file.path(CIBERSORT_DIR, "scripts/cibersort/CIBERSORTx.R")) # CIBERSORTx.R has abs mode option

pre <- DATA_DIR
outpre <- file.path(TCGA_DIR, "data")
ImmCellNum <- 22
LM22.file <- file.path(CIBERSORT_DIR, "scripts/cibersort/LM22_ENSG.txt")

for(cellID in cellShort){ # cellShort: c("HNSC","STAD","COAD","LUAD","PRAD","BRCA","LIHC")
  message(cellID)
  dst <- paste0(cellID,"_long")
  
  cpm.file <- paste0(outpre,"/",dst,"/log2TPM_filter.txt")
  outdir <- paste0(outpre,"/",dst)
  dir.create(outdir,recursive = TRUE, showWarnings = FALSE)
  
  cibersort.res <- CIBERSORT(LM22.file, cpm.file, absolute = TRUE, abs_method = "no.sumto1", perm = 100, QN = FALSE) # 'no.sumto1':absolute space (returns scores) or 'sig.score': relative space (returns fractions)
  write.csv(cibersort.res, paste0(outdir,"/CIBERSORTx_abs.csv"))
  cibersort.res <- read.table(paste0(outdir,"/CIBERSORTx_abs.csv"),sep = ',',header = TRUE,row.names = 1,stringsAsFactors = FALSE,check.names = FALSE)
  
  cibersort.res <- read.table(paste0(TCGA_DIR,"/data/",cellID,"_long","/CIBERSORTx_abs.csv"),sep = ',',header = TRUE,row.names = 1,stringsAsFactors = FALSE,check.names = FALSE)
  cibersort.res <- t(cibersort.res[,1:ImmCellNum])
  data.table::fwrite(cbind(feature=rownames(cibersort.res),cibersort.res),paste0(TCGA_DIR,"/data/",cellID,"_long","/CIBERSORTx_abs.txt"),sep = '\t',row.names = FALSE, col.names = TRUE,quote=FALSE)
}
#


# PCC against CIBERSORTx immune cell fractions: the cluster jobs are
# submitted by scripts/bash/05_pcc_cibersort.sh
outSurfix <- "rmPurity"
#


# ImmReg: ----
## KEGG
outSurfix="rmPurity"
#method="pearson"
df.list <- list()
for(cellID in cellShort){
  # cellID <- "LUAD"
  message(cellID)
  # output/${dst}/pcc_mRNA/${outSurfix}
  df <- data.table::fread( paste0("output/",dst,"/pcc_mRNA/",outSurfix,"_",dst,"_",cellID,".gmtLong"), stringsAsFactors = FALSE, check.names = FALSE, data.table = FALSE, header = FALSE, sep = '\t' )
  colnames(df) <- c("mRNA","sncRNA","source")
  df <- df[(df$source=="ImmReg"),]
  
  df$cell <- cellID
  df.list[[cellID]] <- df
}
df <- as.data.frame(do.call(rbind,df.list))
#table(df$cell)
table(df$cell,df$mRNA)


df_count <- df %>%
  group_by(mRNA, cell) %>%
  summarise(count = n(),log10.number=log10(count+1))
ggplot(data = df_count, aes(x=cell, y=mRNA, size=log10.number, fill=log10.number)) +
  geom_point(shape=21,color="white",stroke=0, alpha=0.9) +
  xlab(paste("")) +
  ylab(paste("")) +
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9, "Blues") )+
  theme_minimal()+
  theme(
    aspect.ratio = 1.2,
    strip.text.x = element_text(size=20),
    strip.text.y = element_text(size=20),
    axis.text.x=element_text(colour="black",angle=90,hjust=1,vjust=0.5,size=16),
    axis.text.y=element_text(colour="black",size=16),
    axis.title.x=element_blank(),
    axis.title.y=element_text(colour="black",size=20),
    plot.title=element_blank(),
    strip.background=element_blank(),
    strip.text=element_blank(),
    # panel.background =element_rect(fill="transparent",colour="black",size=1),
    # panel.grid = element_blank(),
    panel.grid = element_line(),
    legend.position="right",
    legend.title = element_text(colour="black",size=16)
  ) + 
  guides(fill = guide_legend(override.aes = list(size=8))
         # size = guide_legend(override.aes = list(size=16))
  ) 
ggsave("./figure/ImmReg_pathway.pdf",width = 10, height =8)
#



## CIBERSORT
outSurfix="rmPurity"
method="pearson"
df.list <- list()
for(cellID in cellShort){
  # cellID <- "LUAD"
  message(cellID)
  df <- data.table::fread( paste0("output/",dst,"/pcc_CIBERSORT/",cellID,"_",method,"_",outSurfix,".txt"), stringsAsFactors = FALSE, check.names = FALSE, data.table = FALSE, header = TRUE, sep = '\t' )
  df <- df[!is.na(df$pvalue),]
  #table(df$padj<0.01)
  df <- df[df$padj<0.01,] # --sncRES_cutoff 0.9 --p_label p.adjust --p_cutoff 0.05
  df$cell <- cellID
  df.list[[cellID]] <- df
}
df <- as.data.frame(do.call(rbind,df.list))
#table(df$cell)
table(df$cell,df$mRNA)


df_count <- df %>%
  group_by(mRNA, cell) %>%
  summarise(count = n(),log10.number=log10(count+1))
ggplot(data = df_count, aes(x=cell, y=mRNA, size=log10.number, fill=log10.number)) +
  geom_point(shape=21,color="white",stroke=0, alpha=0.9) +
  xlab(paste("")) +
  ylab(paste("")) +
  scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9, "Blues") )+
  theme_minimal()+
  theme(
    aspect.ratio = 1.2,
    strip.text.x = element_text(size=20),
    strip.text.y = element_text(size=20),
    axis.text.x=element_text(colour="black",angle=90,hjust=1,vjust=0.5,size=16),
    axis.text.y=element_text(colour="black",size=16),
    axis.title.x=element_blank(),
    axis.title.y=element_text(colour="black",size=20),
    plot.title=element_blank(),
    strip.background=element_blank(),
    strip.text=element_blank(),
    panel.grid = element_line(),
    legend.position="right",
    legend.title = element_text(colour="black",size=16)
  ) + 
  guides(fill = guide_legend(override.aes = list(size=8))
         # size = guide_legend(override.aes = list(size=16))
  ) 
ggsave("./figure/immreg.pdf",width = 10, height =8)
#



# run GO/KEGG ORA (up+down+allMerged) + GSEA func. enrichment ( by precursor or nearest mRNA, 240805) ---------------

#still use rna+dna
outSurfix="_rmPurity" # _rmMonthOperator
#seqs <- c(1000) # opt1 by rank: 100,200,500,1000
seqs <- c(0.01) # opt2 by p,padj: 0.1, 0.05, 0.01
nearGeneDist <- 100000 # cutoff of nearest distance to gene for intron/promoter/enhancer/repeats peaks

pval.cutoff1 <- 0.01
padj.cutoff1 <- 0.01
pval.cutoff2 <- 0.01 
padj.cutoff2 <- 1 # set to 0.1 not change much
fc.cutoff <- 0.6

dir.create(paste0("./output/",dst,"/pathway/"))

##offline KEGG
m_t <- read.table(file.path(REF_DIR, "geneset/PATH_ID_NAME_KEGGplusHallmarkReactomeImmReg.txt"),sep = "\t",header = TRUE, stringsAsFactors = FALSE) # ,row.names = 1
m_t2g <- m_t[,c("short_DESCRPTION","ensembl_gene_id")]
m_t2e <- m_t[,c("short_DESCRPTION","ENTREZID")]

set.seed(1234)
for(cellID in cellShort){
  # cellID <- "ESCA"
  message(cellID)
  
  # # opt1: enrich by topN
  # peak.list <- list()
  # sigType <- "top"
  # for (top in seqs){
  #   peak.list[[paste0(sigType,top)]] <- list()
  # }
  # for(i in c("tumorigenesis","metastasis")){
  #   # i <- "tumorigenesis"
  #   # message(i)
  #   disease.label <- cmp.df[i,"pos.grp"]
  #   normal.label <- cmp.df[i,"neg.grp"]
  #   col.grp <- cmp.df[i,"col"]
  #   ## select peak by top
  #   x <- paste0("./output/",dst,"/diff_",cellID,"_",disease.label,"vs",normal.label,outSurfix,".diff")
  #   if(!file.exists(x)){
  #     message(x,"not exist, skip next")
  #     next
  #   }
  #   tmp <- data.table::fread(x, header = TRUE, sep="\t",check.names = FALSE,stringsAsFactors = FALSE,data.table = FALSE)
  #   rownames(tmp) <- tmp$feature
  #   tmp$feature <- NULL
  #   tmp <- tmp[order(as.numeric(tmp$pvalue),decreasing = FALSE),]
  #   rownames(tmp) <- unlist(sapply(strsplit(rownames(tmp),"|",fixed=TRUE),"[",4))
  #
  #   # op1: fixed top
  #   for (j in seqs){
  #     peak.id <- rownames(tmp)[1:j]
  #     # print(length(peak.id))
  #     peak.list[[paste0("top",i)]][[disease.label]] <- peak.id
  #    #tmp2 <- filterDiffTabByP(diff = tmp, pvalue = 0.01, padj = 1, abs.log2FC = 0.6)
  
  #   }
  # }
  
  # opt2: enrich by sig
  sigType <- "top"
  peak.list <- list()
  for (top in seqs){
    peak.list[[paste0(sigType,top)]] <- list()
  }
  for(i in c("tumorigenesis","metastasis")){
    # i <- "tumorigenesis"
    # message(i)
    disease.label <- cmp.df[i,"pos.grp"]
    normal.label <- cmp.df[i,"neg.grp"]
    col.grp <- cmp.df[i,"col"]
    ## select peak by top
    x <- paste0("./output/",dst,"/diff_",cellID,"_",disease.label,"vs",normal.label,outSurfix,".diff")
    if(!file.exists(x)){
      message(x,"not exist, skip next")
      next
    }
    tmp <- data.table::fread(x, header = TRUE, sep="\t",check.names = FALSE,stringsAsFactors = FALSE,data.table = FALSE)
    rownames(tmp) <- tmp$feature
    tmp$feature <- NULL
    tmp <- tmp[order(as.numeric(tmp$pvalue),decreasing = FALSE),]
    rownames(tmp) <- unlist(sapply(strsplit(rownames(tmp),"|",fixed=TRUE),"[",4))
    
    # op1: fixed cutoff
    if(i=="tumorigenesis"){
      padj.cutoff <- padj.cutoff1
      pval.cutoff <- pval.cutoff1
    }else{
      padj.cutoff <- padj.cutoff2
      pval.cutoff <- pval.cutoff2
    }
    
    for (j in seqs){
      tmp <- filterDiffTabByP(diff = tmp, pvalue = pval.cutoff, padj = padj.cutoff, abs.log2FC = fc.cutoff)
      peak.id <- rownames(tmp)
      peak.list[[paste0("top",i)]][[disease.label]] <- peak.id
    }
  }
  
  
  for(i in c("tumorigenesis","metastasis")){ # "tumorigenesis","metastasis"
    # i <- "metastasis"
    message(i)
    
    disease.label <- cmp.df[i,"pos.grp"]
    normal.label <- cmp.df[i,"neg.grp"]
    col.grp <- cmp.df[i,"col"]
    
    ### 5.2 read diff tbl
    # diff <- read.table(paste0("./output/",dst,"/diff_",cellID,"_",disease.label,"vs",normal.label,outSurfix,".diff"), header = TRUE, sep="\t",check.names = FALSE)
    diff <- data.table::fread(paste0("./output/",dst,"/diff_",cellID,"_",disease.label,"vs",normal.label,outSurfix,".diff"), header = TRUE, sep="\t",check.names = FALSE,stringsAsFactors = FALSE,data.table = FALSE)
    #cor.test(diff$log2FoldChange,diff2$log2FoldChange)
    rownames(diff) <- diff$feature
    diff$feature <- NULL
    diff <- diff[,c("log2FoldChange","pvalue","baseMean")]
    diff$peak.id <- unlist(sapply(strsplit(rownames(diff),"|",fixed=TRUE),"[",4))
    # diff$RNA <- unlist(sapply(strsplit(rownames(diff),"|",fixed=TRUE),"[",2))
    diff$RNA <- unlist(sapply(strsplit(rownames(diff),"|",fixed=TRUE),"[",2)) #8
    diff$RNA <- gsub("_for|_rev","",diff$RNA,perl = TRUE)
    diff$RNA <- factor(diff$RNA,levels = c(rna,dna)) # rna2
    #table(diff$RNA)
    diff$enst <- unlist(sapply(strsplit(rownames(diff),"|",fixed=TRUE),"[",5))  # 3
    diff$enst2 <- peak.table$tx_ENST[match(diff$enst,peak.table$tx)]
    diff$ensg2 <- peak.table$tx_ENSG[match(diff$enst,peak.table$tx)]
    
    
    ### 5.3 get KEGG ORA/GSEA tbl
    diff2 <- list()
    
    
    ## 5.3.1. mRNA precursor 
    #seem only mRNA,pri_miRNA included in KEGG annotation !!
    #gene universe list (better use all in annotation, instead of all peak enst)
    enst.lab <- "mRNApre" # miRNA
    enst.types <- c("mRNA") # pri_miRNA
    enst.col <- "ensg2"
    #table(diff$RNA %in% enst.types)
    diff2[[enst.lab]] <- diff[diff$RNA %in% enst.types,]
    
    ### ORA
    gene.top.up <- unique(diff2[[enst.lab]][[enst.col]][(diff2[[enst.lab]]$peak.id %in% peak.list[[paste0("top",i)]][[disease.label]]) & (diff2[[enst.lab]]$log2FoldChange>0)])
    gene.top.down <- unique(diff2[[enst.lab]][[enst.col]][(diff2[[enst.lab]]$peak.id %in% peak.list[[paste0("top",i)]][[disease.label]]) & (diff2[[enst.lab]]$log2FoldChange<0)])
    gene.top.up <- gene.top.up[!is.na(gene.top.up)]
    gene.top.down <- gene.top.down[!is.na(gene.top.down)]
    tmp <- ora.kegg.offline(gene.top.up = gene.top.up, gene.top.down = gene.top.down, m_t = m_t2g )
    tmp <- tmp[tmp$ID!="",]
    message(nrow(tmp))
    if(nrow(tmp)>0){
      tmp$top <- top
      tmp$lab <- enst.lab
      write.table(tmp, paste0("./output/",dst,"/pathway/",cellID,"_",disease.label,"vs",normal.label,"_",enst.lab,outSurfix,"_ORA.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
    rm(list = "tmp")
    
    ### GSEA
    tmp <- diff2[[enst.lab]]
    tmp <- tmp[order(tmp$log2FoldChange,decreasing = TRUE),]
    #tmp <- tmp
    rank.list <- tmp$log2FoldChange
    names(rank.list) <- tmp[[enst.col]]
    rank.list <- rank.list[!(is.na(names(rank.list)))]
    rank.list <- rank.list[!duplicated(names(rank.list))]
    tmp <- gsea.kegg.offline(gene.list = rank.list, m_t = m_t2g )
    tmp <- tmp[tmp$ID!="",]
    if(nrow(tmp)>0){
      tmp$top <- top
      tmp$lab <- enst.lab
      write.table(tmp, paste0("./output/",dst,"/pathway/",cellID,"_",disease.label,"vs",normal.label,"_",enst.lab,outSurfix,"_GSEA.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
    rm(list = "tmp")
    
    
    # 5.3.2. 11RNA + 4DNA nearest mRNA
    for(enst.types in c(rna,dna,"all")){
      # enst.types <- "all"
      message(enst.types)
      enst.lab <- paste0(enst.types,"Nearest")
      enst.col <- "nearest.ensg" #"ensg2"
      # enst.lab <- "allNearest"
      if(enst.types == "all"){
        enst.types <- c(rna,dna)
      }
      diff2[[enst.lab]] <- diff[diff$RNA %in% enst.types,]
      diff2[[enst.lab]][[enst.col]] <- peak.table$nearestCodingTx_ENSG[match(diff2[[enst.lab]]$peak.id,peak.table$name)]
      #tmp <- peak.table[peak.table$tx_type=="rRNA",]
      
      ### ORA
      gene.top.up <- unique(diff2[[enst.lab]][[enst.col]][(diff2[[enst.lab]]$peak.id %in% peak.list[[paste0("top",i)]][[disease.label]]) & (diff2[[enst.lab]]$log2FoldChange>0)])
      gene.top.down <- unique(diff2[[enst.lab]][[enst.col]][(diff2[[enst.lab]]$peak.id %in% peak.list[[paste0("top",i)]][[disease.label]]) & (diff2[[enst.lab]]$log2FoldChange<0)])
      gene.top.up <- gene.top.up[!is.na(gene.top.up)]
      #gene.top.up <- gene.top.up[!gene.top.up==""] # not affect
      gene.top.down <- gene.top.down[!is.na(gene.top.down)]
      #gene.top.down <- gene.top.down[!gene.top.down==""]
      tmp <- ora.kegg.offline(gene.top.up = gene.top.up, gene.top.down = gene.top.down, m_t = m_t2g )
      tmp <- tmp[tmp$ID!="",]
      message(nrow(tmp))
      if(nrow(tmp)>0){
        tmp$top <- top
        tmp$lab <- enst.lab
        #write.table(df, paste0(pre,"/tmp/",dst,"-",cellID,"_",enst.lab,"_ORA.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
        write.table(tmp, paste0("./output/",dst,"/pathway/",cellID,"_",disease.label,"vs",normal.label,"_",enst.lab,outSurfix,"_ORA.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
      }
      rm(list = "tmp")
      
      ### GSEA
      tmp <- diff2[[enst.lab]]
      tmp <- tmp[order(tmp$log2FoldChange,decreasing = TRUE),]
      #tmp <- tmp
      rank.list <- tmp$log2FoldChange
      names(rank.list) <- tmp[[enst.col]]
      rank.list <- rank.list[!(is.na(names(rank.list)))]
      rank.list <- rank.list[!duplicated(names(rank.list))]
      tmp <- gsea.kegg.offline(gene.list = rank.list, m_t = m_t2g )
      tmp <- tmp[tmp$ID!="",]
      message(nrow(tmp))
      if(nrow(tmp)>0){
        tmp$top <- top
        tmp$lab <- enst.lab
        write.table(tmp, paste0("./output/",dst,"/pathway/",cellID,"_",disease.label,"vs",normal.label,"_",enst.lab,outSurfix,"_GSEA.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE) 
      }
      rm(list = "tmp")
    }
    
    # # 5.3.3. miR target (deprecated)
    # library(multiMiR)
    
    
    # 5.4.1 ORA dot plot
    #... moved
  }
}
#RColorBrewer::display.brewer.all()








# run GO/KEGG ORA (up+down+allMerged) + GSEA func. enrichment (topN+sig target mRNA by snc_gmt from PCC, 240805)  ---------------

# The cluster jobs for both steps are submitted by
# scripts/bash/06_pathway_ora.sh       - all features together
# scripts/bash/07_pathway_ora_by_rna.sh - one RNA class at a time



# plot GO/KEGG ORA + GSEA --------------------
m_t <- read.table(file.path(REF_DIR, "geneset/PATH_ID_NAME_KEGGplusHallmarkReactomeImmReg.txt"),sep = "\t",header = TRUE, stringsAsFactors = FALSE) # ,row.names = 1
seqs <- 0.01 # 1000
dir.create(paste0("./output/",dst,"/pathway/figure/"))
outSurfix <- "_rmPurity"

## mRNA peak precursor
select.RNAs.lab <- c("mRNA")
sufix <- "pre"

## all peak nearest
# select.RNAs.lab <- c(rna,dna,"all")
# sufix <- "Nearest"

select.RNAs <- paste0(select.RNAs.lab,sufix)
for(enrichType in c("ORA","GSEA")){
  # enrichType <- "ORA"
  message(enrichType)
  
  for(cellID in cellShort){
    # cellID <- "ESCA"
    message(cellID)
    
    for(i in c("tumorigenesis","metastasis")){ # "tumorigenesis","metastasis"
      # i <- "tumorigenesis"
      message(i)
      
      disease.label <- cmp.df[i,"pos.grp"]
      normal.label <- cmp.df[i,"neg.grp"]
      col.grp <- cmp.df[i,"col"]
      
      #str(df)
      df <- list()
      
      for(enst.lab in select.RNAs){ # mRNAmiRNA
        # enst.lab <- "mRNApre"
        x <- paste0("./output/",dst,"/pathway/",cellID,"_",disease.label,"vs",normal.label,"_",enst.lab,outSurfix,"_",enrichType,".txt")
        if(!file.exists(x)){
          next
        }
        df[[enst.lab]] <- read.table(x, sep = "\t", header = TRUE)
      }
      df <- as.data.frame(do.call(rbind,df))
      if(nrow(df)<=1){
        message("not enough, skip")
        next
      }
      #df$minus.log10.pvalue <- -log10(df$pvalue)
      table(df$lab)
      df$lab <- gsub(sufix,"",df$lab)
      df$lab <- factor(df$lab, levels = select.RNAs.lab) # mRNAmiRNA
      table(df$trend)
      df$trend <- factor(df$trend, levels = c("up","dw","all"))
      df <- df[df$trend=="all",]
      
      df.sig <- df[df$pvalue<=0.01,]
      df.sig$source <- m_t$source[match(df.sig$ID,m_t$short_DESCRPTION)]
      table(df.sig$source)
      df.sig <- df.sig[!(df.sig$source %in% c("Reactome")),]
      if(nrow(df.sig)<=1){
        message("not enough, skip")
        next
      }
      
      df.sig$nchr <- nchar(as.character(df.sig$ID))
      max.nchr <- 65
      df.sig$ID[df.sig$nchr>max.nchr] <- paste0(substr(df.sig$ID[df.sig$nchr>max.nchr],1,max.nchr),"...")
      
      top <- seqs
      df.sig.tmp <- df.sig[df.sig$top == top,]
      df.sig.tmp <- df.sig.tmp[order(df.sig.tmp$trend,df.sig.tmp$lab,-df.sig.tmp$nchr,decreasing = TRUE),] # ,df.sig.tmp$nchr
      df.sig.tmp$ID <- factor(df.sig.tmp$ID, levels = unique(df.sig.tmp$ID))
      if(enrichType=="GSEA"){
        df.sig.tmp$Count <- df.sig.tmp$NES #setSize
        sizeLab <- "abs(NES)"
      }else if(enrichType=="ORA"){
        sizeLab = "Count"
      }
      df.sig.tmp$minus.log10.pvalue <- -log10(df.sig.tmp$pvalue)
      ggplot(df.sig.tmp, aes(x=lab, y=ID) ) + # type2, type3
        geom_point(aes(color=minus.log10.pvalue, size=Count)) + # fill="grey",color="grey",  minus.log10.padj
        # scale_fill_manual(values = c("grey50","#ED782F")) + #"#ED782F","#B52B1A","#3E8C9D"
        # scale_color_manual(values = c("grey50","#ED782F")) +
        scale_color_continuous(name="-log10.P",type = "viridis") + # YlOrRd, viridis RColorBrewer::brewer.pal(n=8,name="Reds")
        scale_size_continuous(name = sizeLab) +
        # scale_fill_brewer(type = "seq", palette = RColorBrewer::brewer.pal(n=8,name="YlOrRd")) +
        ylab("") + 
        xlab("") + 
        # scale_x_discrete(breaks = select.RNAs) + # , labels = c("mRNA precursor","all nearest")
        # facet_grid(.~lab) + # ,trend
        theme_minimal() +  # base_size=12 void
        theme(#axis.ticks.x=element_blank(),
          #strip.text.y = element_blank(),
          # aspect.ratio = 10,
          strip.text = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          axis.text.x = element_text(size = 28, hjust = 1,vjust = 1, angle=45), #
          axis.text.y = element_text(size = 20, hjust = 1),
          plot.title = element_text(size=20),
          # strip.text = element_blank(),
          # legend.position = "none", #c(0.9,0.8),#,#
          legend.text = element_text(size= 16),
          legend.title= element_text(size= 16))
      p_wid <- 4+max(df.sig.tmp$nchr)*0.1+length(unique(df.sig.tmp$lab))*0.3
      p_het <- 4+length(unique(df.sig.tmp$ID))*0.5
      ggsave(paste0("./output/",dst,"/pathway/figure/",cellID,"_",top,"_diff_",disease.label,"vs",normal.label,outSurfix,"_kegg_dot_",sufix,"_",enrichType,".pdf"), width = p_wid, height = p_het, limitsize = FALSE) # width = 15, height = 15
      
    }
  }
}




# dot plot of gmt (not by RNA types)
sufix <- "Pcc"
for(enrichType in c("ORA","GSEA")){ # ORA
  #enrichType <- "GSEA" # ORA, GSEA
  for(cellID in cellShort){
    # cellID <- "LUAD"
    message(cellID)
    
    for(i in c("tumorigenesis","metastasis")){ # "tumorigenesis","metastasis"
      # i <- "metastasis"
      message(i)
      disease.label <- cmp.df[i,"pos.grp"]
      normal.label <- cmp.df[i,"neg.grp"]
      col.grp <- cmp.df[i,"col"]
      
      df <- list()
      x <- paste0("./output/",dst,"/pathway/",cellID,"_",disease.label,"vs",normal.label,"/",enrichType,".KEGG.txt") # outSurfix
      if(!file.exists(x)){
        next
      }
      df <- read.table(x, sep = "\t", header = TRUE)
      if(nrow(df)<=1){
        next
      }
      # }
      # df <- as.data.frame(do.call(rbind,df))
      df$minus.log10.pvalue <- -log10(df$pvalue)
      df$minus.log10.padj <- -log10(df$p.adjust)
      df$minus.log10.qvalue <- -log10(df$qvalue)
      df$lab <- 'allPcc'
      df$lab <- gsub(sufix,"",df$lab)
      df$trend <- "all"
      
      df.sig <- df[df$pvalue<=0.05,]
      df.sig$source <- m_t$source[match(df.sig$ID,m_t$short_DESCRPTION)]
      table(df.sig$source)
      df.sig <- df.sig[!(df.sig$source %in% c("Reactome")),]
      if(nrow(df.sig)<=1){
        next
      }
      df.sig$nchr <- nchar(as.character(df.sig$ID))
      max.nchr <- 65
      df.sig$ID[df.sig$nchr>max.nchr] <- paste0(substr(df.sig$ID[df.sig$nchr>max.nchr],1,max.nchr),"...")
      
      #top <- seqs
      df.sig.tmp <- df.sig #[df.sig$top == top,]
      df.sig.tmp <- df.sig.tmp[order(df.sig.tmp$pvalue,decreasing = FALSE),]
      topN <- 30
      if(nrow(df.sig.tmp)>topN){
        df.sig.tmp <- df.sig.tmp[1:topN,] 
      }
      df.sig.tmp <- df.sig.tmp[order(df.sig.tmp$trend,df.sig.tmp$lab,-df.sig.tmp$nchr,decreasing = TRUE),] # ,df.sig.tmp$nchr
      df.sig.tmp$ID <- factor(df.sig.tmp$ID, levels = unique(df.sig.tmp$ID))
      if(enrichType=="GSEA"){
        df.sig.tmp$Count <- df.sig.tmp$NES #setSize
        sizeLab <- "abs(NES)"
      }else if(enrichType=="ORA"){
        sizeLab = "Count"
      }
      
      ggplot(df.sig.tmp, aes(x=lab, y=ID) ) + # type2, type3
        geom_point(aes(color=minus.log10.pvalue, size=Count)) + # fill="grey",color="grey",  minus.log10.padj
        scale_color_continuous(name="-log10.P",type = "viridis") + # YlOrRd, viridis RColorBrewer::brewer.pal(n=8,name="Reds")
        scale_size_continuous(name = sizeLab) +
        ylab("") + 
        xlab("") + 
        theme_minimal() +  # base_size=12 void
        theme(#axis.ticks.x=element_blank(),
          strip.text = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          axis.text.x = element_text(size = 28, hjust = 1,vjust = 1, angle=45), #
          axis.text.y = element_text(size = 20, hjust = 1),
          plot.title = element_text(size=20),
          legend.text = element_text(size= 16),
          legend.title= element_text(size= 16))
      p_wid <- 4+max(df.sig.tmp$nchr)*0.1+length(unique(df.sig.tmp$lab))*0.3
      p_het <- 4+length(unique(df.sig.tmp$ID))*0.5
      # ggsave(paste0("./figure/",cellID,"","_diff_",disease.label,"vs",normal.label,outSurfix,"_kegg_dot_",sufix,"_",enrichType,".pdf"), width = p_wid, height = p_het, limitsize = FALSE) # width = 15, height = 15
      ggsave(paste0("./output/",dst,"/pathway/figure/",cellID,"","_diff_",disease.label,"vs",normal.label,outSurfix,"_kegg_dot_",sufix,"_",enrichType,".pdf"), width = p_wid, height = p_het, limitsize = FALSE) # width = 15, height = 15
    }
  }
}


# dot plot of gmt ( by RNA types)
sufix <- "Pcc"
dst <- "TCGA_small16"
#seqs
outSurfix <- "_rmPurity"
for(enrichType in c("ORA","GSEA")){ #"ORA",
  #enrichType <- "GSEA" # ORA, GSEA
  for(cellID in cellShort){
    # cellID <- "LIHC"
    message(cellID)
    
    for(i in c("tumorigenesis","metastasis")){ # "tumorigenesis","metastasis"
      # i <- "metastasis"
      message(i)
      disease.label <- cmp.df[i,"pos.grp"]
      normal.label <- cmp.df[i,"neg.grp"]
      col.grp <- cmp.df[i,"col"]
      
      df <- list()
      for(RNA in c(rna,dna,"all")){
        # RNA <- "mRNA"
        # message(RNA)
        if(RNA == "all"){
          x <- paste0("./output/",dst,"/pathway/",cellID,"_",disease.label,"vs",normal.label,"/",enrichType,".KEGG.txt") # outSurfix
        }else{
          x <- paste0("./output/",dst,"/pathway/",cellID,"_",disease.label,"vs",normal.label,"_",RNA,"/",enrichType,".KEGG.txt") # outSurfix
        }
        if(!file.exists(x)){
          next
        }
        if(length(readLines(x))<2){
          next
        }
        df.tmp <- read.table(x, sep = "\t", header = TRUE)
        if(nrow(df.tmp)<=1){
          next
        }
        df.tmp$RNA <- RNA
        df[[RNA]] <- df.tmp
      }
      # }
      df <- as.data.frame(do.call(rbind,df))
      df$minus.log10.pvalue <- -log10(df$pvalue)
      df$minus.log10.padj <- -log10(df$p.adjust)
      df$minus.log10.qvalue <- -log10(df$qvalue)
      df$lab <- df$RNA #'allPcc'
      df$lab <- gsub(sufix,"",df$lab)
      #table(df$lab)
      df$lab <- factor(df$lab,levels = c(rna,dna,"all"))
      df$trend <- "all"
      
      df.sig <- df[df$pvalue<=0.05,]
      df.sig$source <- m_t$source[match(df.sig$ID,m_t$short_DESCRPTION)]
      table(df.sig$source)
      df.sig <- df.sig[!(df.sig$source %in% c("Reactome")),]
      if(nrow(df.sig)<=1){
        next
      }
      df.sig$nchr <- nchar(as.character(df.sig$ID))
      max.nchr <- 65
      df.sig$ID[df.sig$nchr>max.nchr] <- paste0(substr(df.sig$ID[df.sig$nchr>max.nchr],1,max.nchr),"...")
      
      #top <- seqs
      df.sig.tmp <- df.sig #[df.sig$top == top,]
      df.sig.tmp <- df.sig.tmp[order(df.sig.tmp$pvalue,decreasing = FALSE),]
      topN <- 30
      if(nrow(df.sig.tmp)>topN){
        df.sig.tmp <- df.sig.tmp[1:topN,] 
      }
      df.sig.tmp <- df.sig.tmp[order(df.sig.tmp$trend,df.sig.tmp$lab,-df.sig.tmp$nchr,decreasing = TRUE),] # ,df.sig.tmp$nchr
      df.sig.tmp$ID <- factor(df.sig.tmp$ID, levels = unique(df.sig.tmp$ID))
      if(enrichType=="GSEA"){
        df.sig.tmp$Count <- df.sig.tmp$NES #setSize
        sizeLab <- "abs(NES)"
      }else if(enrichType=="ORA"){
        sizeLab = "Count"
      }
      
      # table(df.sig$lab)
      ggplot(df.sig.tmp, aes(x=lab, y=ID) ) + # type2, type3
        geom_point(aes(color=minus.log10.pvalue, size=Count)) + # fill="grey",color="grey",  minus.log10.padj
        scale_color_continuous(name="-log10.P",type = "viridis") + # YlOrRd, viridis RColorBrewer::brewer.pal(n=8,name="Reds")
        scale_size_continuous(name = sizeLab) +
        ylab("") + 
        xlab("") + 
        theme_minimal() +  # base_size=12 void 
        theme(#axis.ticks.x=element_blank(),
          strip.text = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          axis.text.x = element_text(size = 28, hjust = 1,vjust = 1, angle=45), #
          axis.text.y = element_text(size = 20, hjust = 1),
          plot.title = element_text(size=20),
          legend.text = element_text(size= 16),
          legend.title= element_text(size= 16))
      p_wid <- 4+max(df.sig.tmp$nchr)*0.1+length(unique(df.sig.tmp$lab))*0.3
      p_het <- 4+length(unique(df.sig.tmp$ID))*0.5
      ggsave(paste0("./output/",dst,"/pathway/figure/",cellID,"_",seqs,"_diff_",disease.label,"vs",normal.label,outSurfix,"_kegg_dot_",sufix,"_",enrichType,".pdf"), width = p_wid, height = p_het, limitsize = FALSE) # width = 15, height = 15
    }
  }
}
#






## KEGG num summary plot
m_t <- read.table(file.path(REF_DIR, "geneset/PATH_ID_NAME_KEGGplusHallmarkReactomeImmReg.txt"),sep = "\t",header = TRUE, stringsAsFactors = FALSE) # ,row.names = 1
cancerPaths <- data.table::fread(file.path(REF_DIR, "geneset/cancer_hallmark/PATH_ID_NAME_hallmark.tsv"),header = TRUE,data.table = FALSE,sep = '\t')

outSurfix <- "_rmPurity"
enrichType <- "GSEA" # "ORA","GSEA"
kegg.sum.list0 <- list()
for(sufix in c("pre","Nearest")){
  # sufix <- "pre"
  message(sufix)
  if(sufix=="pre"){
    select.RNAs <- "mRNApre"
    select.RNAs.lab <- c("mRNA")
  }else if(sufix=="Nearest"){
    select.RNAs <- "allNearest"
    select.RNAs.lab <- "all"  #c(rna,dna,"all")
  }
  
  for(cellID in cellShort){
    # cellID <- "HNSC"
    message(cellID)
    
    for(i in c("tumorigenesis","metastasis")){ # "tumorigenesis","metastasis"
      # i <- "tumorigenesis"
      message(i)
      
      disease.label <- cmp.df[i,"pos.grp"]
      normal.label <- cmp.df[i,"neg.grp"]
      col.grp <- cmp.df[i,"col"]
      
      df <- list()
      for(enst.lab in select.RNAs){ # mRNAmiRNA
        # enst.lab <- "mRNApre"
        x <- paste0("./output/",dst,"/pathway/",cellID,"_",disease.label,"vs",normal.label,"_",enst.lab,outSurfix,"_",enrichType,".txt")
        if(!file.exists(x)){
          next
        }
        df[[enst.lab]] <- read.table(x, sep = "\t", header = TRUE)
      }
      df <- as.data.frame(do.call(rbind,df))
      df <- df[df$ID!="",]
      if(nrow(df)<=1){
        message("not enough, skip")
        next
      }
      #df$minus.log10.pvalue <- -log10(df$pvalue)
      table(df$lab)
      # df$lab <- gsub(sufix,"",df$lab)
      if(sufix=="pre"){
        df <- df[df$lab=="mRNApre",]
      }else if(sufix=="Nearest"){
        df <- df[df$lab=="allNearest",]
      }
      #df$lab <- factor(df$lab, levels = select.RNAs.lab) # mRNAmiRNA
      # table(df$trend)
      # df$trend <- factor(df$trend, levels = c("up","dw","all"))
      df <- df[df$trend=="all",]
      
      tmp <- df
      tmp$cell <- cellID
      tmp$disease <- disease.label
      tmp$normal <- normal.label
      tmp$enrichType <- enrichType
      table(tmp$lab)
      if(sufix=="pre"){
        tmp$targetType <- "cis-precursor"
      }else if(sufix=="Nearest"){
        tmp$targetType <- "cis-nearest"  #c(rna,dna,"all")
      }
      kegg.sum.list0[[paste(sufix,cellID,disease.label,normal.label,enrichType)]] <- tmp
    }
  }
}
df0 <- as.data.frame(do.call(rbind,kegg.sum.list0))
#df <- df[df$pvalue!="",]
df0 <- na.omit(df0)
#table(df0$targetType)

kegg.sum.list <- list()
for(cellID in cellShort){
  # cellID <- "STAD"
  for(i in c("tumorigenesis","metastasis")){ # "tumorigenesis","metastasis"
    # i <- "tumorigenesis"
    message(i)
    disease.label <- cmp.df[i,"pos.grp"]
    normal.label <- cmp.df[i,"neg.grp"]
    col.grp <- cmp.df[i,"col"]
    
    # for(enst.lab in c("mRNAmiRNA","DNA")){
    x <- paste0("./output/",dst,"/pathway/",cellID,"_",disease.label,"vs",normal.label,"/",enrichType,".KEGG.txt") # outSurfix
    if(!file.exists(x)){
      next
    }
    df <- read.table(x, sep = "\t", header = TRUE)
    df <- df[df$ID!="",]
    df$lab=="Pcc"
    if(nrow(df)<=1){
      next
    }
    
    tmp <- df
    tmp$cell <- cellID
    tmp$disease <- disease.label
    tmp$normal <- normal.label
    tmp$enrichType <- enrichType
    tmp$targetType <- "trans"
    kegg.sum.list[[paste(cellID,disease.label,normal.label,enrichType)]] <- tmp
    # }
  }
}
df2 <- as.data.frame(do.call(rbind,kegg.sum.list))
#df <- df[df$pvalue!="",]
df2 <- na.omit(df2)
#dim(df0)
#colnames(df0)

grp1pos <- "Primary"
grp2pos <- "Yes"
df0 <- as.data.frame(rbind(df0[,intersect(colnames(df0),colnames(df2))],df2[,intersect(colnames(df0),colnames(df2))]))
df0$cell <- factor(df0$cell, levels = cellShort)
#df0$lab <- factor(df0$lab, levels = c("mRNAmiRNA","DNA"))
#df0$trend <- factor(df0$trend, levels = c("up","dw"))
df0$lab <- "allRNA" #factor(df$lab, levels = c("mRNAmiRNA","DNA"))
df0$trend <- "allMerge"  #factor(df$trend, levels = c("up","dw"))
df0$minus.log10.padj <- -log10(df0$p.adjust)
df0$minus.log10.pval <- -log10(df0$pvalue)
table(df0$disease)
df0$disease[df0$disease==grp1pos] <- "Tumorigenesis"
df0$disease[df0$disease==grp2pos] <- "Metastasis"
df0$disease <- factor(df0$disease, levels = c("Tumorigenesis","Metastasis"))

#Cancer Immune
df0$source <- m_t$source[match(df0$ID,m_t$short_DESCRPTION)]
table(df0$source)

#Cancer Pathway
df0$source[grepl("notch|nf-|wnt|tnf|p53|hedgehog",tolower(df0$ID), perl = TRUE)] <- "CancerPathway"
table(df0$source)

#Cancer Hallmark
#table(cancerPaths$DESCRPTION)
#table(df0$ID %in% unique(cancerPaths$DESCRPTION))
df0$source[df0$ID %in% unique(cancerPaths$DESCRPTION)] <- "CancerHallmark"



#table(df0$top,df0$lab)
#table(df0$cell,df0$disease)
table(df0$pvalue<=0.05)
table(df0$p.adjust<=0.05)
target.color <- c("cis-nearest"="skyblue","cis-precursor"="salmon","trans"="orange")
cutoffORApadj <- 0.1
df <- df0[df0$p.adjust<=cutoffORApadj,]
#df <- as.data.frame(table(df$cell,df$disease, df$normal, df$lab,df$trend, df$targetType))
#colnames(df) <- c("Cell","Disease","Normal","RNA","Regulation","targetType","Number") # "DEGtype"
df <- as.data.frame(table(df$source,df$cell,df$disease, df$normal, df$lab,df$trend, df$targetType))
colnames(df) <- c("Source","Cell","Disease","Normal","RNA","Regulation","targetType","Number") # "DEGtype",
table(df$Source)
df$targetType <- factor(df$targetType, levels = c("cis-nearest","cis-precursor","trans"))
df$Cell <- factor(df$Cell, levels = cellShort)
#df$RNA <- factor(df$RNA, levels = c("mRNAmiRNA","DNA"))
#df$Regulation <- factor(df$Regulation, levels = c("up","dw"))
df$Disease <- factor(df$Disease, levels = c("Tumorigenesis","Metastasis"))
#table(df$Disease)
#table(df$RNA)
# for(disease.label in c("Yes","Primary")){
#disease.label <- "Yes"
table(df$Source)
#table(df)
ggplot(df[df$Source %in% c("CancerHallmark","CancerPathway","ImmReg"),], aes(x=Cell, y=Number, fill=targetType)) +  # [df$Disease==disease.label,] Regulation
  geom_bar(stat="identity",color="black",  # identity, identity
           position=position_dodge() #fill()
  ) +
  # geom_errorbar(aes(ymin=AUROC-std, ymax=AUROC+std), width=.4,
  #               position=position_dodge(.9))  +
  # geom_text(aes(x=type,y=0.18,label=format(AUROC,digits=3)),size=6,color="white",angle=90)+
  # scale_fill_nejm_adaptive(guide = guide_legend(reverse=FALSE)) + # name = "Dose",
  scale_fill_manual(name="Type",values = target.color) +
  # scale_x_discrete(label=c("miRNA","primary miRNA peaks","all RNA peaks","all DNA peaks","all peaks"))+
  labs(title="",x="", y = "Signif. pathway number", fill="Type")+
  # xlim(c(0,1))+   
  # ylim(c(0,1))+
  theme_minimal() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 0.5,
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,vjust = 0.5,hjust = 1,angle = 90), #
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "right", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16)) +
  facet_grid(Source~Disease,scales = "free_y")
# ggsave("19auroc.pdf")
ggsave(paste0("./figure/path_sigDiff_num","_",enrichType,".png"),width = 12,height = 8)
# ggsave(paste0("./output/",dst,"/figures/path_sigDiff_num_",disease.label,".pdf"),width = 12,height = 10)
# }
#



### by specific pathway
for(path.key in c("G2M checkpoint","Cell cycle") ) {
  #c("Human papillomavirus infection","NF-kappa B signaling","Wnt signaling","Notch signaling")
  #path.key <- "Human papillomavirus infection"
  df.IFN <- df0[grepl(path.key,df0$ID),]
  
  for(pathway in unique(df.IFN$ID)){
    # pathway <- "Ribosome"
    table(df.IFN$cell)
    df.IFN$cell <- factor(df.IFN$cell,levels = cellShort)
    ggplot(df.IFN[df.IFN$ID==pathway,], aes(x=cell, y=minus.log10.padj, fill=targetType)) + 
      geom_bar(stat="identity",color="black",  # identity
               position=position_dodge()
      ) +
      geom_hline(yintercept = 1, color="grey40",linetype="dashed") +
      # geom_errorbar(aes(ymin=AUROC-std, ymax=AUROC+std), width=.4,
      #               position=position_dodge(.9))  +
      # geom_text(aes(x=type,y=0.18,label=format(AUROC,digits=3)),size=6,color="white",angle=90)+
      # scale_fill_nejm_adaptive(guide = guide_legend(reverse=FALSE)) + # name = "Dose",
      scale_fill_manual(name="Type",values = target.color) +
      # scale_x_discrete(label=c("miRNA","primary miRNA peaks","all RNA peaks","all DNA peaks","all peaks"))+
      labs(title=pathway,x="", y = "-log10 P.adjust", fill="Type")+
      # xlim(c(0,1))+   
      # ylim(c(0,1))+
      theme_minimal() +  # base_size=12
      theme(#axis.ticks.x=element_blank(),
        #strip.text.y = element_blank(),
        # aspect.ratio = 0.5,
        strip.text = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size = 20,vjust = 0.5,hjust = 1,angle = 90), #
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size=20),
        # strip.text = element_blank(),
        legend.position = "right", #c(0.9,0.8),#,#
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16)) + 
      facet_grid(disease~lab) #targetType # trend~lab
    # ggsave("19auroc.pdf")
    ggsave(paste0("./figure/","path_sigDiff_padj_",pathway,"_",enrichType,".png"),width = 8,height = 10)
  }
}
#


### by specific sncRNA
## 
sufix <- "PCC"
enrichType="GSEA"
pcc_list <- list()
for(cellID in cellShort){
  message(cellID)
  pcc_list_tmp <- data.table::fread(paste0("output/TCGA_small16/pcc_mRNA/rmPurity_TCGA_small16_",cellID,".txt.gz"), sep = "\t", header = TRUE, nThread = 8, data.table = FALSE)
  pcc_list_tmp$cell <- cellID
  pcc_list_tmp <- pcc_list_tmp[pcc_list_tmp$sncRNA=="peak_199",]  
  if(nrow(pcc_list_tmp)<=1){
    next
  }
  pcc_list[[cellID]] <- pcc_list_tmp
}
df <- as.data.frame(do.call(rbind,pcc_list))
# pcc_df

# df <- as.data.frame(do.call(rbind,df))
df$minus.log10.pvalue <- -log10(df$pvalue)
df$minus.log10.padj <- -log10(df$p.adjust)
df$minus.log10.qvalue <- -log10(df$qvalue)
df$lab <- df$cell #'allPcc'
df$lab <- gsub(sufix,"",df$lab)
#table(df$lab)
df$lab <- factor(df$lab,levels = cellShort) # c(rna,dna,"all")
df$trend <- "all"

df.sig <- df[df$pvalue<=0.05,]
# df.sig$source <- m_t$source[match(df.sig$ID,m_t$short_DESCRPTION)]
# table(df.sig$source)
# df.sig <- df.sig[!(df.sig$source %in% c("Reactome")),]
if(nrow(df.sig)<=1){
  next
}
df.sig$nchr <- nchar(as.character(df.sig$ID))
max.nchr <- 65
df.sig$ID[df.sig$nchr>max.nchr] <- paste0(substr(df.sig$ID[df.sig$nchr>max.nchr],1,max.nchr),"...")


# #opt1: filter by top pvalue
# df.sig.tmp <- df.sig #[df.sig$top == top,]
# df.sig.tmp <- df.sig.tmp[order(df.sig.tmp$pvalue,decreasing = FALSE),]
# topN <- 30
# if(nrow(df.sig.tmp)>topN){
#   df.sig.tmp <- df.sig.tmp[1:topN,] 
# }

#opt2: filter by recur across tumor/cell
topN <- 10
df.sig.tmp2 <- as.data.frame(table(df.sig$ID)) # ,df.sig$cell
df.sig.tmp2 <- df.sig.tmp2[order(df.sig.tmp2$Freq,decreasing = TRUE),]
df.sig.tmp <- df.sig[df.sig$ID %in% df.sig.tmp2$Var1[1:topN],]

df.sig.tmp <- df.sig.tmp[order(df.sig.tmp$trend,df.sig.tmp$lab,-df.sig.tmp$nchr,decreasing = TRUE),] # ,df.sig.tmp$nchr
df.sig.tmp$ID <- factor(df.sig.tmp$ID, levels = unique(df.sig.tmp$ID))
if(enrichType=="GSEA"){
  df.sig.tmp$Count <- df.sig.tmp$NES #setSize
  sizeLab <- "abs(NES)"
}else if(enrichType=="ORA"){
  sizeLab = "Count"
}

# table(df.sig$lab)
ggplot(df.sig.tmp, aes(x=lab, y=ID) ) + # type2, type3
  geom_point(aes(color=minus.log10.pvalue, size=Count)) + # fill="grey",color="grey",  minus.log10.padj
  scale_color_continuous(name="-log10.P",type = "viridis") + # YlOrRd, viridis RColorBrewer::brewer.pal(n=8,name="Reds")
  scale_size_continuous(name = sizeLab) +
  ylab("") + 
  xlab("") + 
  theme_minimal() +  # base_size=12 void 
  theme(#axis.ticks.x=element_blank(),
    strip.text = element_text(size=20),
    axis.title.x = element_text(size=24),
    axis.title.y = element_text(size=24),
    axis.text.x = element_text(size = 28, hjust = 1,vjust = 1, angle=45), #
    axis.text.y = element_text(size = 24, hjust = 1),
    plot.title = element_text(size=20),
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
p_wid <- 6+max(df.sig.tmp$nchr)*0.1+length(unique(df.sig.tmp$lab))*0.3
p_het <- 4+length(unique(df.sig.tmp$ID))*0.5
ggsave(paste0("./figure/",cellID,"_",seqs,"_diff_",disease.label,"vs",normal.label,outSurfix,"_kegg_dot_",sufix,"_",enrichType,".pdf"), width = p_wid, height = p_het, limitsize = FALSE) # width = 15, height = 15
#ggsave(paste0("./output/",dst,"/pathway/figure/",cellID,"_",seqs,"_diff_",disease.label,"vs",normal.label,outSurfix,"_kegg_dot_",sufix,"_",enrichType,".pdf"), width = p_wid, height = p_het, limitsize = FALSE) # width = 15, height = 15



### by all pathway (scatter plot)
#unique(df0$ID)
#colnames(df0)
df0$pval.sig <- 0
df0$pval.sig[df0$p.adjust<=0.1] <- 1 # pvalue

pathType <- "trans" #  cis-nearest cis-precursor trans
df <- df0 %>% 
  tidyr::pivot_wider(id_cols = c("ID","cell","enrichType","targetType","lab","trend","source"),names_from = "disease",values_from = minus.log10.pval,values_fill = 0) %>% 
  dplyr::filter(targetType==pathType)
#table(df0$cell,df0$targetType)
#colnames(df)
#table(df0$cell)
df.sum <- df0 %>% 
  dplyr::group_by(ID,disease,enrichType,targetType,lab,trend,source) %>% 
  dplyr::summarise(recur=sum(pval.sig)) %>% 
  tidyr::pivot_wider(id_cols = c("ID","enrichType","targetType","lab","trend","source"),names_from = "disease",values_from = recur,values_fill = 0) %>% 
  dplyr::filter(targetType==pathType)
colnames(df.sum)[7:8] <- c("Tumorigenesis.Recur","Metastasis.Recur")
#table(df.sum$Metastasis)

df <- left_join(df,df.sum)
df$Recur <- df$Metastasis.Recur+df$Tumorigenesis.Recur # change to max(pri,Metastasis) each
#str(df)
df$Recur <- as.integer(df$Recur)
df$path.label <- df$ID
df$path.label[!df$source %in% c("CancerPathway","CancerHallmark","ImmReg")] <- ""
p <- ggplot(df, aes(x=Tumorigenesis, y=Metastasis, color=Recur)) + 
  geom_point(size=3) +
  # scale_color_continuous( type = "gradient") + # name = "Dose", viridis
  scale_color_gradientn(colors=RColorBrewer::brewer.pal(9, "Blues") )+
  ggrepel::geom_label_repel(col="black",aes(label = path.label,fill=source,alpha=0.8),size = 4, show.legend = TRUE, max.overlaps=10) + # , max.overlaps = Inf
  scale_fill_nejm_adaptive() +
  labs(title="Enrich KEGG",x="-log10 P adjust (Tumorigenesis)", y = "-log10 P adjust (Metastasis)", fill="Pathway")+
  # xlim(c(0,1))+   
  # ylim(c(0,1))+
  theme_classic() +  # base_size=12
  theme(#axis.ticks.x=element_blank(),
    #strip.text.y = element_blank(),
    aspect.ratio = 1,
    strip.background =element_rect(fill="grey90"),
    strip.text = element_text(colour = 'grey20',size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20), #vjust = 0.5,hjust = 1,angle = 90
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    # strip.text = element_blank(),
    legend.position = "right", #c(0.9,0.8),#,#
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16)) + 
  facet_wrap(.~cell,ncol = 4) #targetType # trend~lab
#ggsave("19auroc.pdf")
ggsave(plot = p,filename = paste0("./figure/","path_padj_",pathType,"_",enrichType,"_scatter",".png"),device = "png",dpi = 200,width = 14,height = 12)
#



#by tumor
for(cellID in cellShort){
  # cellID <- "HNSC"
  message(cellID)
  df.tmp <- df[df$cell==cellID,]
  table(duplicated(df.tmp$ID))
  # df.tmp$Metastasis <- scale(df.tmp$Metastasis)
  # df.tmp$Tumorigenesis <- scale(df.tmp$Tumorigenesis)
  p <- ggplot(df.tmp, aes(x=Tumorigenesis, y=Metastasis, color=Recur)) + 
    geom_point(size=3) +
    # scale_color_continuous( type = "viridis") + # name = "Dose",
    scale_color_gradientn(colors=RColorBrewer::brewer.pal(9, "Blues") )+
    ggrepel::geom_label_repel(col="black",aes(label = path.label,fill=source,alpha=0.8),size = 4, show.legend = TRUE, max.overlaps=20) + # , max.overlaps = Inf
    scale_fill_nejm_adaptive() +
    labs(title="Enrich KEGG",x="-log10 P adjust (Tumorigenesis)", y = "-log10 P adjust (Metastasis)", fill="Pathway")+
    # xlim(c(0,1))+   
    # ylim(c(0,1))+
    theme_classic() +  # base_size=12
    theme(#axis.ticks.x=element_blank(),
      #strip.text.y = element_blank(),
      aspect.ratio = 1,
      strip.background =element_rect(fill="grey90"),
      strip.text = element_text(colour = 'grey20',size=20),
      axis.title.x = element_text(size=20),
      axis.title.y = element_text(size=20),
      axis.text.x = element_text(size = 20), #vjust = 0.5,hjust = 1,angle = 90
      axis.text.y = element_text(size = 20),
      plot.title = element_text(size=20),
      # strip.text = element_blank(),
      legend.position = "right", #c(0.9,0.8),#,#
      legend.text = element_text(size= 16),
      legend.title= element_text(size= 16)) 
  # facet_wrap(.~cell,ncol = 3) #targetType # trend~lab
  #ggsave("19auroc.pdf")
  ggsave(plot = p,filename = paste0("./figure/","path_padj_",cellID,"_",pathType,"_",enrichType,"_scatter",".png"),device = "png",dpi = 200,width = 9,height = 7)
  #summary(df$Tumorigenesis)
}
#



# cor 
#
cor.test(df$Tumorigenesis.Recur,df$Tumorigenesis)
cor.test(df$Metastasis.Recur,df$Metastasis)

df.tmp <- data.frame(ID=df$ID,phastCons=df$Tumorigenesis, Pathway_freq=df$Tumorigenesis.Recur, cell=df$cell, Recur=df$Recur)
df.tmp$trend <- "Tumorigenesis"
df.tmp1 <- data.frame(ID=df$ID,phastCons=df$Metastasis,Pathway_freq=df$Metastasis.Recur, cell=df$cell, Recur=df$Recur)
df.tmp1$trend <- "Metastasis"
df.tmp <- as.data.frame(rbind(df.tmp,df.tmp1))
table(duplicated(df$ID))
colnames(df.tmp)
df.tmp <- as_tibble(df.tmp) %>% 
  dplyr::group_by(ID,Pathway_freq,Recur,trend) %>% 
  dplyr::summarise(phastCons=mean(phastCons))
table(duplicated(df.tmp$ID))
# dim(df)
# dim(df.tmp)
p <- 
  ggplot(data=df.tmp,aes(x=Pathway_freq, y=phastCons))+
  # geom_point(alpha=0.5, color="white",fill="grey30", size=2, shape=21) + # fill="grey",color="grey",
  geom_point( aes(fill=Recur),alpha=0.8 ,color="white",shape=21) + # fill="grey",color="grey", # , size=Pathway_freq/3
  # geom_bin2d(bins = 200) +
  # geom_hex(bins = 100) +
  # xlim(c(-10,10)) +
  # xlim(c(-6,6)) +
  scale_fill_viridis_c(name="",direction=1,na.value = "grey30" )+ # values=c(0.1,0.4,0.7,1) # option = "E"
  # scale_fill_manual(values = colorList[["cell"]])+
  scale_color_nejm_adaptive() +
  geom_smooth(method='lm', formula= y~x,color="grey30") +
  ggpubr::stat_cor(label.x = 1,size=8,label.y = 0.8,method = "spearman")+ #this means at 35th unit in the y axis, the r squared and p value will be shown
  # ggpubr::stat_regline_equation(label.x = 1,size=8,label.y = 0.8) +
  # scale_fill_manual(values = c("grey50","#ED782F")) + #"#ED782F","#B52B1A","#3E8C9D"
  # scale_color_manual(values = c("grey50","#ED782F")) +
  geom_vline(xintercept = 0, color="grey", alpha=0.8, linetype="dashed")+
  geom_hline(yintercept = 0, color="grey", alpha=0.8, linetype="dashed")+
  xlab("Signif. pathway frequency") + 
  ylab("-log10 P adjust") + 
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
    legend.title= element_text(size= 16)) +
  facet_grid(.~trend)
ggsave(plot = p, filename = paste0("./figure/","sig_path_cor",".png"),device = "png", width = 12,height = 10, dpi = 200 ) # , limitsize = FALSE
#
#







# unsupervised tumor subtype cluster + boxplot (240808) ---------------------
library(ConsensusClusterPlus)

tmp1 <- strsplit(rownames(logcpm0),"|",fixed=TRUE)
tRNA.table <-  data.frame( row.names = rownames(logcpm0),
                           tRNA=unlist(sapply(tmp1,"[",4)),
                           frag=unlist(sapply(tmp1,"[",2)),
                           RNA=unlist(sapply(tmp1,"[",8))
)


## add gsva
pathwayList <- list()
for(i in cellShort){
  message(i)
  gsva <- read.table(paste0(TCGA_DIR,"/data/",i,"_long/pathway_ssGSEA.txt"),sep = '\t',check.names = FALSE, stringsAsFactors = FALSE,header = TRUE,row.names = 1)
  pathwayList[[i]] <- rownames(gsva)
}
pathways <- (unique(Reduce(intersect,pathwayList)))
length(pathways)
#pathways
# this.paths <- c("Hepatocellular carcinoma","NaturalKiller_Cell_Cytotoxicity","p53 signaling pathway","TNF signaling pathway","NF-kappa B signaling pathway","Hippo signaling pathway") # LIHC
this.paths <- c("Human papillomavirus infection", "Wnt signaling pathway","p53 signaling pathway","TNF signaling pathway","NF-kappa B signaling pathway","Hippo signaling pathway") # LUAD
#

# if(this.cells=="all"){
gsvaList <- list()
for(i in cellShort){
  # i <- "HNSC"
  message(i)
  gsva <- read.table(paste0(TCGA_DIR,"/data/",i,"_long/pathway_ssGSEA.txt"),sep = '\t',check.names = FALSE, stringsAsFactors = FALSE,header = TRUE,row.names = 1)
  print(dim(gsva))
  table(substr(colnames(gsva),14,16))
  gsva <- gsva[,substr(colnames(gsva),14,16)=="01A"] # only keep tumor
  #table(duplicated(colnames(gsva)))
  tmp1 <- substr(colnames(gsva),1,12)
  tmp1[1:3]
  
  table(tmp1 %in% sample.table0$patient)
  gsva <- gsva[,tmp1 %in% sample.table0$patient]
  tmp1 <- substr(colnames(gsva),1,12)
  colnames(gsva) <- sample.table0$sample[match(tmp1,sample.table0$patient)]
  table(colnames(gsva) %in% sample.table0$sample)
  gsva <- gsva[pathways,colnames(gsva) %in% sample.table0$sample]
  gsvaList[[i]] <- gsva
}
names(gsvaList) <- NULL # avoid append name
gsva <- as.data.frame(do.call(cbind,gsvaList))
gsva.bak <- gsva





## run sncRNA cluster
tumor.columns <- c("sample","lineage","cell","group1","tumor_stage","gender","age","race","ESTIMATE","purity")
this.cells <- c("LUAD")  # all, "HNSC","LUAD","LIHC", "COAD"
subtypeLab <- "sncRNA" #"sncRNA", mRNA, tRNA, pri_miRNA
this.RNAs <- c(rna2,dna) # filter only by RNA type: c(rna2,dna),  

outSurfix <- "_rmPurity"
outDir <- paste0(diffpre,"/output/",dst,"/subtype_",paste(this.cells,collapse = "&"))
dir.create(outDir)

#featureLab
for (featureLab in c("pri_miRNA","mRNA","lncRNA","enhancer")){
  # RNA opt1: use single species
  # featureLab <- "pri_miRNA" 

  # RNA opt2: use combined species
  # DES: 
  # featureLab <- c("pri_miRNA","mRNA","lncRNA","enhancer") # c(rna2,dna) c("pri_miRNA","mRNA"),c("pri_miRNA","lncRNA"),c("pri_miRNA","enhancer"),c("pri_miRNA","mRNA","lncRNA"),c("pri_miRNA","mRNA","enhancer"),c("pri_miRNA","lncRNA","enhancer"),c("pri_miRNA","mRNA","lncRNA","enhancer"),
  # surval: 
  # featureLab <- c("pri_miRNA","mRNA","lncRNA","enhancer","promoter")
  print(featureLab)
  
  
  # ## filter by given features
  # #opt1: pan-cancer
  # {
  #   # tmp.diff <- read.table( paste0("./output/",dst,"/diff",outSurfix,".diffSig"),header = TRUE,sep = '\t',check.names = FALSE,stringsAsFactors = FALSE)
  #   tmp.diff2 <- read.table( paste0("./output/",dst,"/","diff.summary"),header = TRUE,sep = '\t',check.names = FALSE,stringsAsFactors = FALSE)
  #   colnames(tmp.diff2)
  #   tmp.diff2 <- tmp.diff2[,c("feature","tumorigenesis.signature cellNum","metastasis.signature cellNum")]
  #   hist(tmp.diff2$`tumorigenesis.signature cellNum`)
  #   hist(tmp.diff2$`metastasis.signature cellNum`)
  #   tmp.diff2 <- tmp.diff2[tmp.diff2$`tumorigenesis.signature cellNum`>=5 & tmp.diff2$`metastasis.signature cellNum`>=3,] 
  #   select.features <- tmp.diff2$feature
  #   length(select.features) #228 sncRNA
  # }
  
  # #opt2: LUAD DES
  # {
  #   tmp.diff2 <- data.table::fread( paste0("./output/",dst,"/","diff.summary"),header = TRUE,sep = '\t',check.names = FALSE,stringsAsFactors = FALSE,data.table = FALSE)
  #   colnames(tmp.diff2)
  #   tmp.diff2 <- tmp.diff2[,grepl(paste0( c("feature",this.cells), collapse = "|"),colnames(tmp.diff2),perl = TRUE)]
  #   table(tmp.diff2$LUAD)
  #   tmp.diff2 <- tmp.diff2[tmp.diff2$LUAD %in% c("Both significant"),]
  #   select.features <- tmp.diff2$feature
  #   tmp <- unlist(sapply(strsplit(select.features,"|",fixed = TRUE),"[",8))
  #   table(tmp)
  #   #old: miRNA:15, mRNA:25, lncRNA:19, enhancer:8
  #   #new: miRNA:14, mRNA:18, lncRNA:15, enhancer:17
  #   select.features <- select.features[tmp %in% c(featureLab)]
  #   message(length(select.features)) 
  #   this.RNAs <- featureLab
  # }

  #opt2: LUAD surv
  {
    tmp.diff2 <- data.table::fread( paste0("./output/",dst,"/","diff.summary"),header = TRUE,sep = '\t',check.names = FALSE,stringsAsFactors = FALSE,data.table = FALSE)
    colnames(tmp.diff2)
    tmp.diff2$name <- unlist(sapply(strsplit(tmp.diff2$feature,"|",fixed = TRUE),"[",4))
    # tmp.diff2 <- tmp.diff2[,grepl(paste0( c("feature",this.cells), collapse = "|"),colnames(tmp.diff2),perl = TRUE)]
    # table(tmp.diff2[[this.cells]])
    # tmp.diff2 <- tmp.diff2[tmp.diff2[[this.cells]] %in% c("Both significant"),]
    # select.features <- tmp.diff2$feature

    cand.snc <- data.table::fread(paste0(diffpre, "/output/",dst,"/candidate/candidate_diff_summary.txt"),header = TRUE,sep = "\t",data.table = FALSE)
    rownames(cand.snc) <- cand.snc$name
    colnames(cand.snc)
    table(cand.snc[[paste0(this.cells,".OS.pval")]]<=0.01) # old: n=11  new: n=31
    cand.snc$feature <- tmp.diff2$feature[match(cand.snc$name,tmp.diff2$name)]
    select.features <- cand.snc$feature[cand.snc[[paste0(this.cells,".OS.pval")]]<=0.01]

    tmp <- unlist(sapply(strsplit(select.features,"|",fixed = TRUE),"[",8))
    tmp.id <- unlist(sapply(strsplit(select.features,"|",fixed = TRUE),"[",4))
    table(tmp)
    # enhancer    lncRNA      mRNA pri_miRNA     snRNA
    # 2         1         3         4         1
    
    #new (all 3k samples):
    # enhancer    intron    lncRNA      mRNA pri_miRNA  promoter     snRNA      tRNA   tucpRNA 
    # 3         1         2         8         6         1         8         1         1 
    #new (test 2k samples):
    # enhancer    lncRNA      mRNA pri_miRNA  promoter     snRNA 
    # 1         1         3         4         1         1 
    select.features <- select.features[tmp %in% c(featureLab)]
    length(select.features) # miRNA:16, mRNA:25, lncRNA:19
    this.RNAs <- featureLab
    cand.snc[tmp.id,c("name",c("gn_type",paste0(this.cells,".OS.pval")))]
  }

  
  cc4ClassList <- list()
  # for(cellID in this.cells ){ # "all" 
  cellID <- this.cells #"all"
  message(cellID)
  if(cellID[1]=="all"){
    sample.table <- sample.table0[ sample.table0$group1!="Normal" & sample.table0$sample %in% colnames(logcpm0),] # 
  }else{
    sample.table <- sample.table0[ sample.table0$group1!="Normal" & sample.table0$sample %in% colnames(logcpm0)  & (sample.table0$cellShort %in% cellID) ,] # 
  }
  #rownames(sample.table) <- sample.table$sample
  
  #colnames(sample.table)
  sample.table <- sample.table[,tumor.columns]
  #str(sample.table)
  sample.table <- sample.table[order(sample.table$lineage,sample.table$cell,sample.table$group1),]
  dim(sample.table)
  
  vsd <- logcpm0[(tRNA.table$RNA %in% this.RNAs) & (rownames(logcpm0) %in% select.features), rownames(sample.table)]
  mads=apply(vsd,1,mad)
  length(mads)
  #hist(mads, breaks=nrow(vsd)*0.1)
  mad2k=vsd[rev(order(mads))[1:min(10000,nrow(vsd))],]
  
  clusterOut <- cellID
  if(ncol(mad2k)>300){
    plotType <- "png"
    # png(filename = paste0("./figure/",clusterOut,"_consensus.",plotType), width = 1200, height = 1200, res = 300)
    
  }else{
    plotType <- "pdf"
    # pdf(filename = paste0("./figure/",clusterOut,"_consensus.",plotType), width = 12, height = 12)
    
  }
  results = ConsensusClusterPlus(as.matrix(mad2k),
                                 tmyPal=viridis::viridis_pal()(100),
                                 maxK=8, # max cluster num
                                 reps=50,
                                 pItem=0.8,
                                 pFeature=1,
                                 title= clusterOut,
                                 clusterAlg="km", #"pam",
                                 distance="spearman",
                                 seed=1262118388.71279,
                                 plot=plotType)
  #dev.off()
  if(ncol(mad2k)<=200){
    file.copy(from = paste0(clusterOut,"/consensus.",paste(featureLab,collapse = "&"),".",plotType), to = paste0("./figure/",clusterOut,"_consensus.",paste(featureLab,collapse = "&"),".",plotType), overwrite = TRUE)
    file.remove(paste0(clusterOut,"/consensus.",paste(featureLab,collapse = "&"),".",plotType))
    unlink(clusterOut, recursive = TRUE)
  }
  
  #Silhouette width analysis
  cc4 = results[[4]] # use 4 cluster result
  # # calcultaing Silhouette width using the cluster package 
  # library(cluster)
  # cc4Sil = silhouette(x = cc4[[3]], # x is a numeric vector that indicates cluster assignment for each data point
  #                     dist = as.matrix(1- cc4[[4]])) # dist should be a distance matrix and NOT a similarity matrix, so I subtract the matrix from one to get that dist matrix
  # library(factoextra)
  # fviz_silhouette(cc4Sil, palette = "jco",
  #                 ggtheme = theme_classic())
  
  #omics subtype overlap
  cc4Class = data.frame(cc4$consensusClass)
  cc4Class$sample = rownames(cc4Class)
  cc4ClassList[[cellID]] <- cc4Class
  # }
  cc4ClassDf <- as.data.frame(do.call(rbind,cc4ClassList))
  dim(cc4ClassDf) # 
  
  
  if(this.cells=="all"){
    sample.table <- sample.table0[  sample.table0$group1!="Normal" & sample.table0$sample %in% colnames(logcpm0),] # 
  }else{
    sample.table <- sample.table0[  sample.table0$group1!="Normal" & sample.table0$sample %in% colnames(logcpm0)  & (sample.table0$cellShort %in% this.cells),] # 
  }
  sample.table <- sample.table[,c(tumor.columns,"patient")  ]
  sample.table$lineage <- factor(sample.table$lineage, levels = lineageShort)
  sample.table$cell <- factor(sample.table$cell, levels = cellShort)
  sample.table$group1 <- factor(sample.table$group1, levels = group3Short)
  sample.table <- sample.table[order(sample.table$lineage,sample.table$cell,sample.table$group1,decreasing = FALSE),]
  table(duplicated(sample.table$sample))
  dim(sample.table)
  
  dim(cc4ClassDf)
  sample.table <- left_join(sample.table,cc4ClassDf,by="sample")
  table(duplicated(sample.table$sample))
  #cc4Class = data.frame(cc4Class[match(sample.table$sample,cc4Class$ID),])
  #all(cc4Class$ID == (sample.table$sample))
  colnames(sample.table)[length(colnames(sample.table))] <- subtypeLab
  #sample.table <- sample.table[!duplicated(sample.table$patient),]
  #rownames(sample.table) <- sample.table$patient
  table(duplicated(sample.table$sample))
  rownames(sample.table) <- sample.table$sample
  dim(sample.table)
  #sample.table$cell <- NULL
  
  subtype.col <- subtypeLab #this.cells
  subtype.col <- subtype.col[subtype.col %in% colnames(sample.table)]
  df2 <- as.data.frame(sample.table[,c(subtype.col),drop=FALSE])
  # df2 <- df2[!is.na(df2$Integrative),,drop=FALSE]
  # df2[,subtype.col] <- lapply(df2[,subtype.col], factor)
  # df2[,subtype.col] <- lapply(df2[,subtype.col], as.numeric)
  table(df2[[subtype.col]])
  data.table::fwrite(df2,file = paste0(outDir,"/",paste(this.cells,collapse = "&"),"_",paste(this.RNAs,collapse = "&"),"_",subtypeLab,"_",paste(featureLab,collapse = "&"),"_subtype.txt"),row.names = FALSE)
  #rownames(df2) <- df2$patient
  #df2$sample <- NULL
  #df2$Integrative <- NULL
  #df2 <- data.table::fread(paste0(outDir,"/",paste(this.cells,collapse = "&"),"_",paste(this.RNAs,collapse = "&"),"_subtype.txt"),row.names = FALSE)
  
  
  gsva <- gsva.bak
  table(colnames(gsva) %in% sample.table$sample)
  gsva <- gsva[pathways,colnames(gsva) %in% sample.table$sample]
  dim(gsva)
  gsva[nrow(gsva)+1,] <- colnames(gsva)
  rownames(gsva)[nrow(gsva)] <- "sample"
  gsva <- as.data.frame(t(gsva))
  dim(gsva)
  gsva <- gsva[sample.table$sample,]
  rownames(gsva) <- sample.table$sample
  gsva$sample <- sample.table$sample
  #gsva <- na.omit(gsva)
  table(is.na(gsva))
  #gsva[is.na(gsva)] <- NaN
  gsva[,colnames(gsva)[colnames(gsva)!="sample"]] <- lapply(gsva[,colnames(gsva)[colnames(gsva)!="sample"]],as.numeric)
  
  
  #sample.table <- as.data.frame(cbind(sample.table,gsva))
  sample.table <- as.data.frame(dplyr::left_join(sample.table,gsva[,c("sample",this.paths)],by="sample"))
  rownames(sample.table) <- sample.table$sample
  sample.table[,this.paths] <- lapply(sample.table[,this.paths],as.numeric)
  if(length(this.cells)==1){
    sample.table[[subtype.col]] <- as.factor(sample.table[[subtype.col]])
  }else{
    sample.table[,subtype.col] <- lapply(sample.table[,subtype.col],as.factor)
  }
  
  
  # single omic RNA cluster subtypes
  mat <- logcpm0[(tRNA.table$RNA %in% this.RNAs) & (rownames(logcpm0) %in% select.features), rownames(sample.table)]
  table(is.na(as.matrix(mat)) | is.nan(as.matrix(mat)) )
  dim(sample.table)
  table(colnames(mat) %in% rownames(sample.table) )
  #brkList <- seq(-1, 1, by = 0.1)
  #sample.table[["sncRNA"]] <- sample.table[["mRNA"]]
  colorList.tmp <- colorList
  for(i in subtype.col){
    sample.table[[i]] <- as.factor(sample.table[[i]])
    colorList.tmp[[i]] <- pal_npg_adaptive()( length(unique(sample.table[[i]])) )
    names(colorList.tmp[[i]]) <- sort(unique(as.character(sample.table[[i]])))
  }
  #colnames(sample.table)
  if(length(this.cells)==1){
    sample.table <- sample.table[order(sample.table[[subtype.col]],sample.table[["lineage"]],sample.table[["cell"]],sample.table[["group1"]],decreasing = FALSE),]
  }else{
    sample.table <- sample.table[order(sample.table[["lineage"]],sample.table[["cell"]],sample.table[["group1"]],decreasing = FALSE),]
  }
  
  #dim(mat)
  table(rownames(mat) %in% rownames(tRNA.table))
  #str(sample.table)
  mat <- mat[,sample.table$sample]
  colorList.tmp[['RNA']] <- colorList.tmp[['precursor']][names(colorList.tmp[['precursor']]) %in% c(rna2,dna) ]
  colorList.tmp[['frag']] <- colorList.tmp[['precursor']][names(colorList.tmp[['precursor']]) %in% c(rna,dna) ]
  colorList.tmp[c(tumor.columns[tumor.columns!="sample"])]
  #colnames(sample.table)
  #table(is.na(sample.table$CNA))
  tmpL <- table(sample.table[[subtype.col]])
  dim(mat)
  pheatmap::pheatmap(mat = mat,
                     annotation_col = sample.table[colnames(mat),c(tumor.columns[!tumor.columns %in% c("sample","patient")], this.paths, subtype.col)], #data frame that specifies the annotations shown on left side of the heatmap.
                     annotation_row = tRNA.table[rownames(mat),c("frag","RNA")],
                     annotation_colors = colorList.tmp, # RColorBrewer::brewer.pal(6, "Set3")
                     border_color="NA",
                     scale = "row", # row,none
                     #labels_col = 3, labels_row = 6,
                     cluster_cols = FALSE,cluster_rows = TRUE,
                     # gaps_row = sum(annotation_row$trend=="up"), # invalid if cluster row/col
                     gaps_col = as.integer(cumsum(tmpL)),
                     #cutree_cols = 2,cutree_rows = 3,
                     show_colnames=FALSE, show_rownames=FALSE,
                     fontsize = 10,
                     legend =  FALSE,
                     height = 6+0.2*length(c(this.cells,this.paths)), width = 18,
                     color = rev( colorRampPalette( RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))(100) ) , # ),  RdYlBu  viridis::viridis_pal()(length(brkList))
                     # breaks = brkList,
                     #fontsize_row = 5,
                     filename = paste0("./figure/",dst,"_",paste(this.cells,collapse = "&"),"_",paste(this.RNAs,collapse = "&"),"_",subtypeLab,"_subtype_heat","_",paste(featureLab,collapse = "&"),".png")
  )
  sample.table.bak <- sample.table
  #
  
  
  
  
  
  
  
  # plot meta/ssGSEA mat boxplot among groups/subtypes 
  #add treatment/drug (xena TCGA not GDA, LIHC et al. no response data)
  treatDf <- data.table::fread(paste0(TCGA_DIR,"/data/TCGA-all_clinical/",this.cells,"_clinical/TCGA.",this.cells,".sampleMap%2F",this.cells,"_clinicalMatrix"),sep='\t',data.table = FALSE,header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
  treatDf <- treatDf[,grepl("sample|therap|treat|response",colnames(treatDf),ignore.case = TRUE,perl = TRUE)]
  #treatDf[1:4,]
  treatDf$patient <- substr(treatDf$sampleID,1,12)
  treatDf <- treatDf[!duplicated(treatDf$patient),]
  treatDf$sample <- sample.table0$sample[match(treatDf$patient,sample.table0$patient)]
  table(duplicated(treatDf$sample))
  treatDf <- treatDf[!duplicated(treatDf$sample),]
  colnames(treatDf) # primary_therapy_outcome_success: category
  table(treatDf$primary_therapy_outcome_success)
  treatDf$response <- "nonresponse"
  treatDf$response[treatDf$primary_therapy_outcome_success %in% c("Complete Remission/Response","Partial Remission/Response")] <- "response"
  table(treatDf$response)
  
  sample.table <- sample.table.bak
  sample.table <- as.data.frame(dplyr::left_join(sample.table,treatDf[,c("sample","response")],by="sample")) # 
  #summary(gsva$`TCGA-4P-AA8J-01A`)
  rownames(sample.table) <- sample.table$sample

  
  for(varType in c("num", "category")){
    # varType <- "num"
    message(varType)
    if(varType=="num"){
      this.vars <- c( # "age",
        "ESTIMATE", # num clinical
        # "Human papillomavirus infection",
        "Wnt signaling pathway",
        # "Hepatocellular carcinoma",
        # "NaturalKiller_Cell_Cytotoxicity",
        "p53 signaling pathway"
        # "TNF signaling pathway",
        # "NF-kappa B signaling pathway",
        # "Hippo signaling pathway" # path
      )
    }else if(varType=="category"){
      this.vars <- c(
        # "tumor_stage", # num clinical
        "group1", #
        # "subtype",
        "gender",
        # "race",
        "response" # treatment/response
      )    
    }

    
    subGrp.list <- list()
    subGrps <- c("sncRNA") # ,"protein""sncRNA","miRNA",
    for(subGrp in subGrps){
      # subGrp <- "sncRNA"
      sample.table$group <- sample.table[[subGrp]]
      logcpm <- as.data.frame(t(sample.table[,c(this.vars),drop=FALSE]))
      
      #library(dplyr)
      # tmp.feature <- splitPeakFeatureID(rownames(logcpm))
      gene.mat <- as.data.frame(t(logcpm[rownames(logcpm) %in% this.vars,]))
      # colnames(gene.mat) <- tmp.feature$PeakID[match(colnames(gene.mat),tmp.feature$feature)]
      gene.mat.t <- gene.mat
      gene.mat.t$sample <- rownames(gene.mat.t)
      gene.mat.t[1:4,]
      gene.mat.t <- left_join(gene.mat.t,sample.table[,c("sample","group","cell")],by="sample")
      dim(gene.mat.t)
      gene.mat.t$subGrp <- subGrp
      subGrp.list[[subGrp]] <- gene.mat.t
    }
    subGrp.df <- as.data.frame(do.call(rbind,subGrp.list))
    
    ## plot expr&pathway expr
    df <- reshape2::melt(subGrp.df,id.var=c("sample","group","cell","subGrp"),measure.var=this.vars)
    df$group <- factor(df$group)
    df <- df[df$cell %in% this.cells,]
    set.seed(1)
    b <- runif(nrow(df), -0.2, 0.2)
    library(ggplot2)
    df$variable <- gsub(" signaling pathway|Human|","",df$variable)
    conflict_prefer("fisher.test", "stats")
    if(varType=="num"){
      p0 <- ggplot(df,aes(x=group,y=value,fill=group))+  # ,label=label
        geom_boxplot(aes(x=(group),y=value,fill=group),color="black",alpha=0.6,outlier.size = 0) +
        geom_point(aes(x=(as.numeric(group)+b),y=value,fill=group),
                   shape=21,size=2,color="black",stroke=0.5) +
        ggpubr::stat_compare_means(
          ref.group="1",
          label.x.npc = 0.5,hide.ns=FALSE,size =8, vjust = 0.5, step.increase = 0.04,#paired = TRUE,
          # aes(label = p.signif), # p.signif, group = Group,
          label="p.signif",
          symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.1,1),symbols = c( "****","***", "**", "*", "ns")),
          method = "wilcox.test") +
        facet_grid(variable~subGrp, scales = "free")
      
      tmpL <- pal_npg_adaptive()( length(unique(1:8)) )
      names(tmpL) <- unique(1:8)
    }else if(varType == "category"){
      df <- df %>% 
        dplyr::mutate(n=1) %>%
        dplyr::group_by(group,cell,subGrp,variable,value) %>%  #  sample
        dplyr::summarize(num=sum(n))
      # colnames(df)
      res.List <- list()
      for(i in unique(df$subGrp)){
        for(j in unique(df$variable)){
          #table(subGrp.df$group)
          contingency_table <- subGrp.df[subGrp.df$subGrp==i,c("group",j)]
          contingency_table <- as.data.frame.matrix(table(contingency_table[[j]],contingency_table[["group"]]))
          ## opt1: chi test
          # chi_squared_test <- chisq.test(contingency_table)
          ## opt2: fisher
          res <- fisherExactFromCountTable(contingency_table)
          # res$min.padj <- min(res$adjusted_p_value)
          # res$min.pval <- min(res$p_value)
          res$subGrp <- i
          res$variable <- j
          res.List[[paste(i,j)]] <- res
        }
      }
      res.List <- as.data.frame(do.call(rbind,res.List))
      res.List <- as_tibble(res.List) %>% 
        dplyr::group_by(subGrp,variable,grp1,grp2) %>% 
        dplyr::summarise(min.padj=min(adjusted_p_value),
                         min.pval=min(p_value)
        )
      res.List <- res.List[order(res.List$subGrp,res.List$variable,decreasing = FALSE),]
      p0 <- ggplot(df,aes(x=group,y=num,fill=value))+  # ,label=label
        geom_bar(position = "fill",width = 0.8,size=0.3, stat = "identity", color="black") +  # identity,count
        # ylim(c(0,1.25)) +
        # geom_text(position = position_fill(vjust = 0.5), aes(label = ), size = 7)
        # geom_text(vjust = -0.2, hjust = 0.5, size = 6) +  # Add sig labels
        geom_text(aes(label = num),stat = "identity", size=6,position = position_fill(.5)) + # stat = "prop",
        facet_grid(subGrp~variable, scales = "free_y")
      # p0
      tmpL <- colorList[this.vars]
      names(tmpL) <- NULL
    }
    
    p <- p0 +
      # scale_fill_nejm_adaptive() +
      scale_fill_manual(values = unlist(tmpL) ) + 
      # scale_x_discrete(label=c("CT","CN","PT","PN","RT","RN"))+
      # labs(title = g)+
      ylab("")+
      xlab("") +
      theme_bw(base_size=12) +  # minimal
      theme(#axis.ticks.x=element_blank(),
        # panel.border = element_rect(colour="white", fill = NA), #for panel.border # 
        strip.text.y = element_text(family="sans",size=20), # serif
        strip.text.x = element_text(family="sans",size=20), # serif
        axis.title.x = element_text(family="sans",size=20),
        axis.title.y = element_text(family="sans",size=20),
        axis.text.x = element_text(size = 20,family="sans"), # ,angle = 45,hjust=1,vjust=1
        axis.text.y = element_text(size = 20,family="sans"),
        plot.title = element_text(size=20),
        legend.position = "right",
        legend.title = element_text(family="sans",size=24),
        legend.text = element_text(family="sans",size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
    # facet_grid(variable~subGrp, scales = "free")
    # facet_grid(subGrp~variable, scales = "free_y")
    print(p)
    if(varType=="category"){
      ggsave(plot = p,filename = paste0("./figure/subtype_",cellID,"_",varType,"_",paste(featureLab,collapse = "&"),".png"), device = "png", height = 4+0.5*length(subGrps),width =7+1*length(unique(df$variable)) ) # category
      #res.List$min.padj <- pvalue2plabel(res.List$min.padj)
      print(res.List)
    }else if(varType == "num"){
      ggsave(plot = p,filename = paste0("./figure/subtype_",cellID,"_",varType,"_",paste(featureLab,collapse = "&"),".png"), device = "png", height = 8+length(unique(df$variable)),width =6+0.5*length(subGrps) ) # num
    }
  }
  # 
  
  
  
  
  
  # plot meta/ssGSEA mat KM-surv among groups/subtypes
  #add survival (xena TCGA)
  survDf <- data.table::fread(paste0(TCGA_DIR,"/data/TCGA-all_clinical/",this.cells,"_clinical/","survival%2F",this.cells,"_survival.txt"),sep='\t',data.table = FALSE,header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
  survDf <- survDf[,grepl("sample|DFI|PFI|OS",colnames(survDf),ignore.case = TRUE,perl = TRUE)]
  #survDf[1:4,]
  survDf$patient <- substr(survDf$sample,1,12)
  survDf <- survDf[!duplicated(survDf$patient),]
  survDf$OS.time <- round(survDf$OS.time/30,digits = 0)
  survDf$DFI.time <- round(survDf$DFI.time/30,digits = 0)
  survDf$PFI.time <- round(survDf$PFI.time/30,digits = 0)
  
  sample.table <- sample.table.bak
  sample.table$patient <- sample.table0$patient[match(sample.table$sample,sample.table0$sample)]
  sample.table <- as.data.frame(dplyr::left_join(sample.table,survDf,by="patient")) #
  
  # # uncomment this two lines for late stage patient
  # table(sample.table$tumor_stage)
  # sample.table <- sample.table[sample.table$tumor_stage %in% c("stage iii","stage iv"),]
  
  ## plot KM survial
  # R 3.6 compatibility patch
  if(getRversion() < "4.0.0") {
    deparse1 <- function(expr, collapse = " ", width.cutoff = 500L, ...) {
      paste(deparse(expr, width.cutoff, ...), collapse = collapse)
    }
    assignInNamespace("deparse1", deparse1, ns = "base")
    assign("deparse1", deparse1, envir = .GlobalEnv)
  }
  # library(tidyverse)
  library(janitor)
  #install.packages ("tidyquant")
  library(tidyquant)
  library(patchwork)
  library(survival)
  library(survminer)
  
  xs <- c("sncRNA","tumor_stage") #  ,"miRNA","mRNA" ,"subtype","tumor_stage" (tumor_stage actually only need run once)
  for (x in xs){
    # x <- "sncRNA"
    message(x)
    for(surv in c("OS")){ # ,"PFI","DFI"
      # surv <- "OS"
      message(surv)
      
      rm(list=c("exprVec.cut","exprVec","sfit_tmp","res.suvr","sample.table.tmp","surv.pval"))
      sample.table.tmp <- sample.table[!is.na(sample.table[[x]]) & !is.nan(sample.table[[x]]),]
      tryCatch(
        {
          plotSurv(exprVec = as.factor(sample.table.tmp[[x]]), sample.table = sample.table.tmp,facet = "NULL",time.col = paste0(surv,".time"), event.col = surv,palette="npg", outFile = paste0("./figure/",dst,"_",cellID,"_",x,"_surv_",surv,"_",paste(featureLab,collapse = "&"),".pdf"),width = 5,height = 5)
        },
        error = function(e) {
          # stop("too many zero expr? skip")
          # surv.pval <- NaN
          # exprVec.cut=NaN
        }
      )
    }
  }
}

#



# supervised tumor subtyping ----------------

this.cells <- "LUAD"
sample.table <- sample.table0[sample.table0$cell %in% this.cells,] # COAD   BRCA
treatDf <- data.table::fread(paste0(TCGA_DIR,"/data/TCGA-all_clinical/",this.cells,"_clinical/TCGA.",this.cells,".sampleMap%2F",this.cells,"_clinicalMatrix"),sep='\t',data.table = FALSE,header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
treatDf <- treatDf[,grepl("sample|therap|treat|response",colnames(treatDf),ignore.case = TRUE,perl = TRUE)]
#treatDf[1:4,]
treatDf$patient <- substr(treatDf$sampleID,1,12)
treatDf <- treatDf[!duplicated(treatDf$patient),]
treatDf$sample <- sample.table0$sample[match(treatDf$patient,sample.table0$patient)]
table(duplicated(treatDf$sample))
treatDf <- treatDf[!duplicated(treatDf$sample),]
colnames(treatDf) # primary_therapy_outcome_success: category
table(treatDf$primary_therapy_outcome_success)
treatDf$response <- "nonresponse"
treatDf$response[treatDf$primary_therapy_outcome_success %in% c("Complete Remission/Response","Partial Remission/Response")] <- "response"
table(treatDf$response)
sample.table <- left_join(sample.table,treatDf,by="sample")
sample.table <- sample.table[,c("sample","response")]
colnames(sample.table) <- c("sample","group")
table(is.na(sample.table$group))
sample.table <- sample.table[!is.na(sample.table$group),]
write.table(sample.table[,c("sample","group")],paste0(DATA_DIR,"/exSeek-dev/data/",dst,"/sample_table_",this.cells,"_subtype.txt"),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
#


# train & test in TCGA LUAD
# The job is submitted by: MODE=subtype CELL=LUAD bash scripts/bash/08_random_forest.sh


# supervised normal tissue-of-origin ----------------
# prep mat
## all normal
sample.table <- sample.table0[sample.table0$group=="Normal",]
#unique(sample.table$cell[!sample.table$cell %in% names(tumor2tissue)])
sample.table$group <- sample.table$cell
sample.table$group <- gsub("THCA","thyroid",sample.table$group,perl = TRUE)
sample.table$group <- gsub("ESCA","esophagus",sample.table$group,perl = TRUE)
sample.table$group <- gsub("PRAD","prostate",sample.table$group,perl = TRUE)
sample.table$group <- gsub("LUSC|LUAD","lung",sample.table$group,perl = TRUE)
sample.table$group <- gsub("STAD","stomach",sample.table$group,perl = TRUE)
sample.table$group <- gsub("LIHC","liver",sample.table$group,perl = TRUE)
sample.table$group <- gsub("KIRC|KIRP|KICH","kidney",sample.table$group,perl = TRUE)
sample.table$group <- gsub("BLCA","bladder",sample.table$group,perl = TRUE)
sample.table$group <- gsub("COAD","bowel-colon",sample.table$group,perl = TRUE)
sample.table$group <- gsub("HNSC","head-neck",sample.table$group,perl = TRUE)
sample.table$group <- gsub("BRCA","breast",sample.table$group,perl = TRUE)
sample.table$group <- gsub("UCEC","uterine",sample.table$group,perl = TRUE)
sample.table$group <- gsub("CHOL","bile-duct",sample.table$group,perl = TRUE)
table(sample.table$group)
write.table(sample.table[,c("sample","group")],paste0(DATA_DIR,"/exSeek-dev/data/",dst,"/sample_table_NC_TOO.txt"),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
tumor2tissue <- c("THCA"="thyroid","ESCA"="esophagus", "PRAD"="prostate","LUSC,LUAD"="lung","STAD"="stomach","LIHC"="liver","KIRC,KIRP,KICH"="kidney","BLAD"="bladder","COAD"="bowel-colon")
#tissue2tumor <- setNames(names(tumor2tissue), tumor2tissue)
write.table(sample.table[sample.table$group %in% tumor2tissue,c("sample","group")],paste0(DATA_DIR,"/exSeek-dev/data/",dst,"/sample_table_NC_TOO_GSE163534.txt"),quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)
#table(sample.table$group %in% tumor2tissue)



# train & test in TCGA normal
# The job is submitted by: MODE=too bash scripts/bash/08_random_forest.sh


# re-plot confusion mat in TCGA use all normal
average_conf_matrix <- read.csv(paste0(PROJECT_DIR,"/output/",dst,"/TOO/average_confusion_matrix.csv"),row.names = 1,check.names = FALSE,header = TRUE,stringsAsFactors = FALSE)
conf_matrix_df <- as.data.frame((as.matrix(average_conf_matrix))) # as.table
colnames(conf_matrix_df) <- c("Prediction", "Reference", "Frequency")
conf_matrix_df$Frequency <- as.numeric(conf_matrix_df$Frequency)
tissueShort <- c("head-neck","thyroid","lung","bowel-colon","stomach","esophagus","liver","bile-duct","bladder","kidney","prostate","uterine","breast")
#unique(conf_matrix_df$Prediction)
conf_matrix_df$Prediction <- factor(conf_matrix_df$Prediction, levels = tissueShort)
conf_matrix_df$Reference <- factor(conf_matrix_df$Reference, levels = tissueShort)
# Plot confusion matrix with ggplot2
ggplot(conf_matrix_df, aes(x = Prediction, y = Reference, fill = Frequency)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(title = "", x = "Prediction", y = "True") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0( "./figure/TCGA_average_confusion_matrix_plot.png"), width = 5, height = 4)
#






# evaluate peak/domain mature miR + AGO2 + DASHRdb boundary hist and overlap --------------------------------
## get pri-miR & miR hg38 coordinates
setwd(DATA_DIR)
suppressPackageStartupMessages(library(bedtoolsr))
library(tidyverse)
library(dplyr)

## read overlap bed file
read.mir.intersect <- function(x){
  #x <- "output/AGO2_IP/call_domain/domains/5/05_gn_interectMIR.bed"
  # x <- "output/GSE133684/call_peak_dedupByPos/domains/b5_p01_gn_interectMIR.bed"
  print(x)
  rbp <- data.table::fread(x,data.table = FALSE, header = FALSE,sep = "\t",stringsAsFactors = FALSE)
  idx <- colnames(rbp)[ncol(rbp)]
  
  rbp$intersect.ratio <- as.numeric(rbp[[idx]]/(rbp$V3-rbp$V2)) * as.numeric(rbp[[idx]]/(rbp$V15-rbp$V14))
  
  rbp <- rbp[order(rbp$V4,rbp$intersect.ratio,rbp[[idx]],decreasing = TRUE),] # rm dup
  rbp <- rbp[!duplicated(rbp$V4),]
  rbp <- na.omit(rbp) # na exist in intersect.ratio
  # table(rbp[[idx]]==0)
  rbp <- rbp[(rbp[[idx]]>=10 | rbp$intersect.ratio>=0.01) & (rbp$V15-rbp$V14)<=65,]
  
  rbp$delta.5 <- rbp$V2 - rbp$V14 # start - start
  #hist(rbp$delta.5,breaks = 10,col="grey50",xlim = c(-10,10))
  rbp$delta.3 <- rbp$V3 - rbp$V15 # start - start
  #hist(rbp$delta.3,breaks = 10,col="grey50",xlim = c(-10,10))
  rbp <- rbp[,c(1:6,ncol(rbp)-4,ncol(rbp)-3,ncol(rbp)-2,ncol(rbp)-1,ncol(rbp))]
  #plot(x=rbp$delta.5,y=rbp$delta.3)
  rbp$path <- x
  #rbp$p <- as.numeric(paste0("0.",gsub("_gn_interectMIR.bed","",basename(x),fixed = TRUE)))
  return(rbp)
}

#need define best decay and p cutoff statistics/param
#AGO2_IP/call_peak_dedup_bk, GSE71008/call_peak_all
dst <- "TCGA_small16" # "TCGA-COAD_small_NC_NCpool" # GSE71008_NCpool, GSE110381_NCpool, WSQ_SMARTer_NEB_NCpool (NCpool_NEBNext)
smp <- "b5_d50_p1"  #"NCpool_NEBNext"
pre <- paste0(dst,"/call_peak_all") # "AGO2_IP/call_peak_dedup_bk", "GSE71008/call_peak_all", "GSE110381/call_peak_all" 
res1 <- lapply(c(
  paste0("output/",pre,"/cfpeakCNN/",smp,"_gn_interectAGO2.bed"),
  paste0("output/",pre,"/cfpeakCNN/",smp,"_gn_interectDASHR.bed"),
  paste0("output/",pre,"/cfpeakCNN/",smp,"_gn_interectMIR.bed") 
),
read.mir.intersect)
#not apply to long seq, no intersect records?

res1.df <- as_tibble(do.call("rbind",res1))
table(res1.df$path)
res1.df$method <- ""
res1.df$method[grepl("AGO2",res1.df$path)] <- "AGO2"
res1.df$method[grepl("DASHR",res1.df$path)] <- "DASHRdb"
# res1.df$method[grepl("localmax",res1.df$path)] <- "LocalMax"
res1.df$method[grepl("MIR",res1.df$path)] <- "miR"
table(res1.df$method)

x <- "05"
res1.df$p <- as.numeric(paste0("0.",x))
table(res1.df$p )
# 0 0.05  0.1  0.2  0.5  0.9 
# 351  405  405  407  410  414
res1.df2 <- res1.df %>% 
  dplyr::filter(p==0.05) %>% # 0.01
  tidyr::pivot_longer(cols = delta.3:delta.5, names_to = "direction",values_to = "distance")
res1.df2$direction <- ifelse(res1.df2$direction=="delta.3","3'","5'")
res1.df2$direction <- factor(res1.df2$direction,levels = c("5'","3'"))
table(res1.df2$direction)


## plot histplot: boundary postion 
library(plyr); library(dplyr)
mu <- ddply(res1.df2, "method", plyr::summarise, grp.mean=mean(distance))
library(dplyr)
mu <- res1.df2 %>% 
  dplyr::group_by(direction,method) %>% 
  dplyr::summarise(grp.mean=mean(distance), grp.sum=dplyr::n(), grp.zero=sum(distance==0), ratio=grp.zero/grp.sum )
#mu$method <- factor(mu$method,levels = c("Piranha","CLIPper","CLAM","exPeak"))
summary(mu$ratio)

#head(res1.df2)
library(ggplot2)
ggplot(res1.df2, aes(x=distance, color=method, fill=method)) +
  geom_histogram( position="identity", alpha=0.9, color="black")+
  geom_vline(xintercept = 0,color="grey30",linetype="dashed")+
  scale_fill_manual(values = c("grey30","salmon","orange")) + # pal_nejm_adaptive()(4)
  # scale_color_manual(values = pal_nejm_adaptive()(4))+
  labs(title="",x="MIR end relative position (5'/3')", y = "Predicted position frequency")+
  facet_grid(method~direction,scales = "free")+
  # ylim(c(0,1))+
  xlim(c(-50,50))+
  theme_classic(base_size=12) + 
  theme(#axis.ticks.x=element_blank(),  
    #strip.text.y = element_blank(),
    #strip.text.x = element_text(face="bold",family="arial",size=20),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size = 20,hjust = 0.5,vjust = 0.5), # ,angle = 90
    axis.text.y = element_text(size = 20),
    plot.title = element_text(size=20),
    strip.text = element_text(size = 20),
    legend.position = "none",# 
    legend.text = element_text(size= 16),
    legend.title= element_text(size= 16))
ggsave("./figure/mir_boundary.pdf",width = 9,height = 6) # ,height = 14
mu






# split by RNA
res1.df2$RNA <- peak.table$gn_type[match(res1.df2$V4,peak.table$name)]
table(res1.df2$RNA)

for( i in c("AGO2","DASHRdb","miR")){
  # i <- "AGO2"
  message(i)
  
  ggplot(res1.df2[res1.df2$method==i,], aes(x=distance, color=method, fill=method)) +
    geom_histogram( position="identity", alpha=0.9, color="black")+
    geom_vline(xintercept = 0,color="grey30",linetype="dashed")+
    scale_fill_manual(values = c("grey30","salmon","orange")) + # pal_nejm_adaptive()(4)
    # scale_color_manual(values = pal_nejm_adaptive()(4))+
    labs(title="",x= paste0(i," end relative position (5'/3')"), y = "Predicted position frequency")+
    facet_grid(RNA~direction,scales = "free")+
    # ylim(c(0,1))+
    xlim(c(-50,50))+
    theme_classic(base_size=12) + 
    theme(#axis.ticks.x=element_blank(),  
      #strip.text.y = element_blank(),
      #strip.text.x = element_text(face="bold",family="arial",size=20),
      axis.title.x = element_text(size=20),
      axis.title.y = element_text(size=20),
      axis.text.x = element_text(size = 20,hjust = 0.5,vjust = 0.5), # ,angle = 90
      axis.text.y = element_text(size = 20),
      plot.title = element_text(size=20),
      strip.text = element_text(size = 20),
      legend.position = "none",# 
      legend.text = element_text(size= 16),
      legend.title= element_text(size= 16))
  ggsave( paste0("./figure/",i,"_boundary.pdf"), width = 9, height = 26 ) # ,height = 14
}
#


getIntersectPeak <- function(x){
  # x <- paste0(outpre,"/output/",dst,"/call_peak_all","/cfpeakCNN/",smp,"_gn_interectDASHR.bed")
  rbp <- data.table::fread( x, data.table = FALSE, header = FALSE,sep = "\t",stringsAsFactors = FALSE )
  idx <- colnames(rbp)[ncol(rbp)]
  
  rbp$intersect.ratio <- as.numeric(rbp[[idx]]/(rbp$V3-rbp$V2)) * as.numeric(rbp[[idx]]/(rbp$V15-rbp$V14))
  
  rbp <- rbp[order(rbp$V4,rbp$intersect.ratio,rbp[[idx]],decreasing = TRUE),] # rm dup
  rbp <- rbp[!duplicated(rbp$V4),]
  totalPeakNum <- nrow(rbp)
  
  rbp <- na.omit(rbp) # na exist in intersect.ratio
  # table(rbp[[idx]]==0)
  #rbp <- rbp[(rbp[[idx]]>=10 | rbp$intersect.ratio>=0.01) & (rbp$V15-rbp$V14)<=65,]
  rbp <- rbp[(rbp[[idx]]>=5 | rbp$intersect.ratio>=0.01) & (rbp$V15-rbp$V14)<=200,]
  return(rbp$V4)
}

res2 <- lapply(c(
  paste0(outpre,"/output/",dst,"/call_peak_all","/cfpeakCNN/","b5_d50_p1","_gn_interectAGO2.bed"),
  paste0(outpre,"/output/",dst,"/call_peak_all","/cfpeakCNN/","b5_d50_p1","_gn_interectDASHR.bed"),
  paste0(outpre,"/output/",dst,"/call_peak_all","/cfpeakCNN/","b5_d50_p1","_gn_interectMIR.bed")
),
getIntersectPeak)
#res2[[1]] #

VennDiagram::venn.diagram(x = list(res2[[1]], res2[[2]], res2[[3]]),category.names = c("AGO2","DASHR","MIR"),
                          fill=ggplot2::alpha(RColorBrewer::brewer.pal(3, "Set3"),0.96),
                          force.unique = TRUE,print.mode = "raw",
                          filename = paste0("./figure/","intersect_venn.png"),imagetype="png",height = 4000, width = 4000,
                          cex = 2.5,
                          cat.cex = 2
                          #cat.pos = c(-10,+10,0),cat.dist = c(0.09,0.02)
)
#







# unitas: diff TCGA_small16 tsRNA (by cell type, grp1+2 edgeR, 202408) ----------------
#diff count & map strategy, need seperate rm batch
count0.tRNA <- data.table::fread(paste0(outpre,"/output/",dst,"/tsRNA_unitas_matrix/tsRNA_unitas_count.txt"), header = TRUE, sep="\t",check.names = FALSE,data.table=FALSE) # in line with long-seq
rownames(count0.tRNA) <- count0.tRNA$feature
count0.tRNA$feature <- NULL
dim(count0.tRNA)

## run diff
outDir = paste0(diffpre,"/output/",dst,"/")
dir.create(outDir,recursive = TRUE,showWarnings = FALSE)

outSurfix="_rmPurity"
#cellShort <- cellShort[cellShort!="EV"]
for(cellID in cellShort){
  # cellID <- "STAD"
  message(cellID)
  
  for(i in c("tumorigenesis","metastasis")){ # "tumorigenesis","metastasis"
    # i <- "metastasis"
    message(i)
    
    pos.grp <- cmp.df[i,"pos.grp"]
    neg.grp <- cmp.df[i,"neg.grp"]
    col.grp <- cmp.df[i,"col"]
    
    method="edger_glmlrt"
    norm_method="TMM"
    
    sample.table <- sample.table0[ (sample.table0$cellShort==cellID),] # 1042
    sample.table[["group"]] <- sample.table[[col.grp]] 
    
    positive_samples <- sample.table[sample.table$group==pos.grp,"sample"]  # 
    negative_samples <- sample.table[sample.table$group==neg.grp,"sample"] #
    if(length(positive_samples)==0 | length(negative_samples)==0 ) {
      message(cellID,i,"skip for zero smp")
      next
    }
    samples <- c(positive_samples, negative_samples)
    group <- c(rep("positive",length(positive_samples)),rep("negative",length(negative_samples)))
    sample.table <- sample.table[samples,]
    sample.table$group <- factor(sample.table$group,levels = c(neg.grp,pos.grp))
    # sample.table$Month <- factor(sample.table$Month)
    
    mat <- count0.tRNA[,samples]
    useExternalLibSize <- externalLibSize[samples]
    #str(sample.table)
    
    form <- ~ group + purity # + Month + Operator
    res <- diff.v3.dcb(mat = mat, form = form,sample.table = sample.table,samples = samples, group = group,
                       method = method, norm_method = norm_method, filterType = "NULL", featureType = "domain", 
                       getCountGt1SmpFrac = TRUE, getPercExpByGrp = TRUE, useExternalLibSize = useExternalLibSize)
    data.table::fwrite(cbind(feature=rownames(res[["diffTable"]]), res[["diffTable"]]),paste0(outDir,"/unitas_",cellID,"_",pos.grp,"vs",neg.grp,"",outSurfix,".diff"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE) # 
  }
}
#



# unitas: tRFs length,ratio ------------------------
## unitas tRFs reads meta distribution
#unitasDir <- paste0(outpre,"/output/TCGA-",dst,"_small/tsRNA_unitas/",sid)
#unitasList <- list()
readUnitasTbl <- function(unitasDir){
  x <- paste0(unitasDir,"/unitas.tRF-table.txt")
  if(file.exists(x)){
    tmp <- read.table(x,skip = 3, check.names = FALSE, sep = "\t",header = TRUE,row.names = 1)
    tmp <- tmp[,!grepl("absolute",colnames(tmp))] # 2020,PNAS,Quake used fractionated
    #colnames(tmp) <- gsub(" (fractionated)","",colnames(tmp),perl = FALSE)
    colnames(tmp) <- unlist(sapply(strsplit(colnames(tmp)," "),"[",1))
    # unique(rownames(tmp))
    tmp$id <- rownames(tmp)
    tmp$tx.type <- unlist(sapply(strsplit(tmp$id,"-"),"[",2))
    # tmp$tx.type2 <- unlist(sapply(strsplit(tmp$id,"-"),"[",3))
    # table(tmp$tx.type)
    tmp$path <- unitasDir
    #unitasList[[unitasDir]] <- tmp
    return(tmp)
  }
}

outpre <- DATA_DIR
# cellID <- "STAD"
# sample.table <- sample.table0[sample.table0$cellShort==cellID,]
# table(sample.table$group)
sid <- sample.table0$sample
#sid <- sid[1]
frag.types <- c("tRNA-leader","5p-tRFs","5p-tR-halves","misc-tRFs","3p-tR-halves","3p-tRFs","3p-CCA-tRFs","tRF-1") #sort(unique(unitasDF.tsR2$frag.type))


unitasList <- mclapply( paste0(outpre,"/output/",dst,"/tsRNA_unitas/",sid), FUN = readUnitasTbl, mc.cores = 10)
#lapply( paste0(outpre,"/output/",dst,"/tsRNA_unitas/",sid), FUN = readUnitasTbl)
unitasDF <- as.data.frame(do.call(rbind,unitasList))
unitasDF[1:3,]


dim(unitasDF)
unitasDF$tx.type2 <- unlist(sapply(strsplit(unitasDF$id,"-"),"[",3))
unitasDF$tx.type2[is.na(unitasDF$tx.type2)] <- ""
unitasDF$tx.type3 <- paste0(unitasDF$tx.type, "-",unitasDF$tx.type2)
unitasDF$tx.type3[(unitasDF$tx.type2)==""] <- unitasDF$id[(unitasDF$tx.type2)==""]
#table(unitasDF$tx.type3)
unitasDF$dst <- splitVec(unitasDF$path,delim = "/",fixed = TRUE,nth = 8)
unitasDF$smp <- splitVec(unitasDF$path,delim = "/",fixed = TRUE,nth = 10)
unitasDF$cell <- sample.table0$cellShort[match(unitasDF$smp,sample.table0$sample)]
unitasDF$group <- sample.table0$group[match(unitasDF$smp,sample.table0$sample)]
unitasDF$tx.source <- "Nuc"
unitasDF$tx.source[unitasDF$tx.type %in% c(paste0("T",toupper(letters)),paste0("T",toupper(letters),1))] <- "Mt"
table(unitasDF$tx.source)
#table(unitasDF$cell)
#table(unitasDF$tx.type)
smpN <- length(unique(unitasDF$smp))

unitasDF$libsize <- externalLibSize[unitasDF$smp]

unitasDF.long <- as_tibble(unitasDF) %>% 
  tidyr::pivot_longer(cols = "5p-tR-halves":"misc-tRFs",names_to = "frag.type",values_to = "count") %>% 
  dplyr::mutate(cpm=count/libsize*1e6)
# unitasDF.wid <- unitasDF.long %>% 
#   tidyr::pivot_wider(names_from = c(smp), values_from = "cpm") # err...
#data.table::fwrite(unitasDF.wid,"sample_matrix")  

#summary(unitasDF$libsize)
#unitasDF$norm.count <- unitasDF$


## get tsRNA peak info
#colnames(unitasDF)
unitasDF.tsR <- unitasDF.long %>% 
  dplyr::select("id","tx.type","tx.type2","tx.type3","tx.source","frag.type") %>% 
  dplyr::distinct_all() %>% 
  as.data.frame()
#rownames(unitasDF.tsR)
#colnames(unitasDF.tsR)
unitasDF.tsR2 <- unitasDF.tsR[!duplicated(paste0(unitasDF.tsR$tx.type3,"|",unitasDF.tsR$frag.type)),c("tx.type", "tx.type2","tx.type3","tx.source","frag.type")]
rownames(unitasDF.tsR2) <- paste0(unitasDF.tsR2$tx.type3,"|",unitasDF.tsR2$frag.type)
unitasDF.tsR3 <- unitasDF.tsR[!duplicated(paste0(unitasDF.tsR$tx.type3)),c("tx.type", "tx.type2","tx.type3","tx.source")]
rownames(unitasDF.tsR3) <- paste0(unitasDF.tsR3$tx.type3)

## aver by Glu-AGG
df1 <- as_tibble(unitasDF.long) %>% 
  # tidyr::pivot_longer(cols = "5p-tR-halves":"misc-tRFs",names_to = "frag.type",values_to = "count") %>% 
  dplyr::group_by(tx.type3,frag.type,cell,group) %>% 
  dplyr::summarize(mean.cpm=mean(cpm))
summary(df1$mean.cpm)
table(is.na(df1$mean.cpm))
df1.wid <- df1 %>%
  tidyr::pivot_wider(names_from = c("cell","group"), values_from = "mean.cpm") %>% 
  as.data.frame()
rownames(df1.wid) <- paste0(df1.wid$tx.type3,"|",df1.wid$frag.type)
df1.wid[,c("tx.type3","frag.type")] <- NULL
df1.wid[(is.na(df1.wid))] <- 0
df1.wid <- df1.wid[apply(df1.wid,1,sd)>0,]

# hc3 <- hclust(stats::dist(tmp.mat2))
# pdf("tRNAfrag_hclust.pdf",width = 8,height = 4)
# plot(hc3, hang = -1)
# dev.off()
#brkList <- seq(0, 1, by = 0.01)
# brkList <- seq(-1, 1, by = 0.1)
#min(cor.mat)
colorList[["cell"]] <- colorList[["cellShort"]]

tmp <- pal_npg_adaptive()(8)
names(tmp) <- frag.types
# 3p-CCA-tRFs 3p-tR-halves      3p-tRFs 5p-tR-halves      5p-tRFs    misc-tRFs        tRF-1  tRNA-leader 
# "#E64B35FF"  "#4DBBD5FF"  "#00A087FF"  "#3C5488FF"  "#F39B7FFF"  "#8491B4FF"  "#91D1C2FF"  "#DC0000FF" 
colorList[["frag.type"]] <- tmp
colorList[["tx.source"]] <- c("Nuc"="grey30","Mt"="orange")
#dim(df1.wid)
annotation_col <- data.frame(id=colnames(df1.wid))
annotation_col$cell <- unlist(sapply(strsplit(annotation_col$id,"_"),"[",1))
annotation_col$cell <- factor(annotation_col$cell,levels = cellShort)
annotation_col$group <- unlist(sapply(strsplit(annotation_col$id,"_"),"[",2))
annotation_col$group <- factor(annotation_col$group,levels = groupShort)
rownames(annotation_col) <- annotation_col$id
annotation_col$id <- NULL
pheatmap::pheatmap(mat = df1.wid,
                   annotation_col = annotation_col[colnames(df1.wid),], #data frame that specifies the annotations shown on left side of the heatmap.
                   annotation_row = unitasDF.tsR2[rownames(df1.wid),c("tx.source","frag.type")],
                   annotation_colors = colorList, # RColorBrewer::brewer.pal(6, "Set3")
                   border_color="NA",
                   scale = "row", # row,none
                   #labels_col = 3, labels_row = 6,
                   cluster_cols = TRUE,cluster_rows = TRUE,
                   # gaps_row = sum(annotation_row$trend=="up"), # invalid if cluster row/col
                   #cutree_cols = 2,cutree_rows = 3,
                   show_colnames=FALSE, show_rownames=FALSE,
                   fontsize = 10,
                   legend =  TRUE, 
                   height = 12, width = 14,
                   color = viridis::viridis_pal()(100), #rev( colorRampPalette( RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))( 100 ) ) , # ),  OrRd
                   # breaks = brkList,
                   #fontsize_row = 5,
                   filename = paste0("./figure/",dst,"_tsR_heat.png")
                   # filename = paste0("./output/",dst,"/allSmpP01_exPeak_smallDomain_diff_",disease.label,"vs",normal.label,".cpm.heatmap.pdf")
)
df1.wid <- df1.wid[,grepl("_Normal",colnames(df1.wid))]
table(unitasDF.tsR2$tx.source[match(rownames(df1.wid),rownames(unitasDF.tsR2))])
df1.wid <- df1.wid[unitasDF.tsR2$tx.source[match(rownames(df1.wid),rownames(unitasDF.tsR2))]=="Nuc",]
df1.wid <- df1.wid[apply(df1.wid,1,sd)>0,]
dim(df1.wid)
pheatmap::pheatmap(mat = df1.wid,
                   annotation_col = annotation_col[colnames(df1.wid),], #data frame that specifies the annotations shown on left side of the heatmap.
                   annotation_row = unitasDF.tsR2[rownames(df1.wid),c("tx.source","frag.type")],
                   annotation_colors = colorList, # RColorBrewer::brewer.pal(6, "Set3")
                   border_color="NA",
                   scale = "row", # row,none
                   #labels_col = 3, labels_row = 6,
                   cluster_cols = TRUE,cluster_rows = TRUE,
                   # gaps_row = sum(annotation_row$trend=="up"), # invalid if cluster row/col
                   #cutree_cols = 2,cutree_rows = 3,
                   show_colnames=FALSE, show_rownames=FALSE,
                   fontsize = 10,
                   legend =  TRUE, 
                   height = 12, width = 14,
                   color = viridis::viridis_pal()(100), #rev( colorRampPalette( RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))( 100 ) ) , # ),  OrRd
                   # breaks = brkList,
                   #fontsize_row = 5,
                   filename = paste0("./figure/",dst,"_tsR_heat_NC1.png")
                   # filename = paste0("./output/",dst,"/allSmpP01_exPeak_smallDomain_diff_",disease.label,"vs",normal.label,".cpm.heatmap.pdf")
)
#


## zscore across cell
df1.wid2 <- df1 %>%
  dplyr::group_by(tx.type3,frag.type,group) %>% 
  dplyr::mutate(scale.cpm=scale(mean.cpm)) %>%  # zscore across cell
  dplyr::select(-mean.cpm) #%>% 
# tidyr::pivot_wider(names_from = c("cell","group","frag.type"), values_from = "scale.cpm") %>% 
# as.data.frame()
#table(duplicated(df1.wid2$tx.type3))
#rownames(df1.wid2) <- paste0(df1.wid2$tx.type3)
#df1.wid2[,c("tx.type3")] <- NULL
table((is.na(df1.wid2)))
#table(is.nan(df1.wid2))
#tmp <- unlist(sapply(strsplit(colnames(df1.wid2),"[",fixed = TRUE),"[",1))
#table(tmp)
table(duplicated(colnames(df1.wid2)))
#df1.wid2[(is.na(df1.wid2))] <- 0
#df1.wid2 <- df1.wid2[apply(df1.wid2,1,sd)>0,]
#df1.wid2 <- df1.wid2[,grepl("_HC",colnames(df1.wid2))]

#facet need ggplot
#df1.wid2$leftPos <- factor(df1.wid2$leftPos)
df1.wid2$cell <- factor(df1.wid2$cell,levels = cellShort)
unitasDF.tsR3 <- unitasDF.tsR3[order(unitasDF.tsR3$tx.source),]
df1.wid2$tx.type3 <- factor(df1.wid2$tx.type3, levels = unique(rownames(unitasDF.tsR3)))
tmp <- df1.wid2[df1.wid2$group=="Normal",]
library(ggh4x)
# 3p-CCA-tRFs 3p-tR-halves      3p-tRFs 5p-tR-halves      5p-tRFs    misc-tRFs        tRF-1  tRNA-leader
# "#E64B35FF"  "#4DBBD5FF"  "#00A087FF"  "#3C5488FF"  "#F39B7FFF"  "#8491B4FF"  "#91D1C2FF"  "#DC0000FF"
#table(tmp$tx.type3,tmp$frag.type)
tmp$frag.type <- factor(tmp$frag.type,levels = frag.types)
#colorList[["frag.type"]] 
ggplot(tmp, aes(cell, tx.type3)) + 
  geom_tile(aes(fill = scale.cpm)) + 
  # scale_fill_gradient(low = "white", high = "red") + 
  # facet_wrap(~frag.type,nrow = 1)+ 
  ylab("") + xlab("") +
  scale_fill_continuous(type = "viridis") +
  facet_wrap2(.~frag.type, nrow = 1, #scales = 'free',
              strip = strip_themed(
                background_x = list(element_rect(fill = "#E64B35FF"),
                                    element_rect(fill = "#4DBBD5FF"),
                                    element_rect(fill = "#00A087FF"),
                                    element_rect(fill = "#3C5488FF"),
                                    element_rect(fill = "#F39B7FFF"),
                                    element_rect(fill = "#8491B4FF"),
                                    element_rect(fill = "#91D1C2FF"),
                                    element_rect(fill = "#DC0000FF"))
              )
  ) +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.text = element_text(size= 16,color = "black"),
        axis.text.x = element_text(size= 16,color = "black",angle = 45,hjust = 1), # angle = 45,,hjust = 1
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        strip.text = element_text(size= 28,color = "black"),
        panel.border = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = c(rep("grey50",length(other)),rep("steelblue",length(rna)),rep("firebrick",length(dna)))), #  
        legend.position = "right",#c(.25,.6),
        legend.text = element_text(size= 20),
        legend.title= element_text(size= 24))
ggsave(paste0("./figure/",dst,"_tsR_heat_NC2.png"),width = 40, height = 14,limitsize = FALSE)


## zscore across frag.type
df1.wid2 <- df1 %>%
  dplyr::group_by(tx.type3,cell,group) %>% 
  dplyr::mutate(scale.cpm=scale(mean.cpm)) %>%  # zscore across cell
  dplyr::select(-mean.cpm) #%>% 
table((is.na(df1.wid2)))
table(duplicated(colnames(df1.wid2)))
df1.wid2$cell <- factor(df1.wid2$cell,levels = cellShort)
unitasDF.tsR3 <- unitasDF.tsR3[order(unitasDF.tsR3$tx.source),]
df1.wid2$tx.type3 <- factor(df1.wid2$tx.type3, levels = unique(rownames(unitasDF.tsR3)))
tmp <- df1.wid2[df1.wid2$group=="Normal",]
library(ggh4x)
ggplot(tmp, aes(frag.type, tx.type3)) + 
  geom_tile(aes(fill = scale.cpm)) + 
  # scale_fill_gradient(low = "white", high = "red") + 
  # facet_wrap(~frag.type,nrow = 1)+ 
  ylab("") + xlab("") +
  scale_fill_continuous(type = "viridis") +
  facet_wrap2(.~cell, nrow = 1) +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.text = element_text(size= 16,color = "black"),
        axis.text.x = element_text(size= 16,color = "black",angle = 45,hjust = 1), # angle = 45,,hjust = 1
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        strip.text = element_text(size= 28,color = "black"),
        strip.background=element_blank(),
        panel.border = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = c(rep("grey50",length(other)),rep("steelblue",length(rna)),rep("firebrick",length(dna)))), #  
        legend.position = "right",#c(.25,.6),
        legend.text = element_text(size= 20),
        legend.title= element_text(size= 24))
ggsave(paste0("./figure/",dst,"_tsR_heat_NC3.png"),width = 40, height = 14,limitsize = FALSE)


#opt1: 
unitasDF.long$group1 <-  sample.table0$group1[match(unitasDF.long$smp,sample.table0$sample)]
#opt2: 
# unitasDF.long$group1 <-  sample.table0$cell[match(unitasDF.long$smp,sample.table0$sample)]

#table(unitasDF.long$group)
df2 <- as_tibble(unitasDF.long) %>%
  dplyr::filter(group %in% c("Tumor","Normal")) %>% # opt1
  # dplyr::filter(group=="Normal") %>% # opt2
  dplyr::group_by(tx.source,tx.type,frag.type,cell,group1) %>%
  dplyr::summarize(cpm=mean(cpm,na.rm=TRUE))
# dplyr::summarize(cpm=mean(count/libsize*1e6)/smpN)
summary(df2$cpm)
table(df2$group1)

#df2$tx.type <- factor(df2$tx.type,levels = hc3$labels[hc3$order])
df2$frag.type <- factor(df2$frag.type,levels = frag.types)
#str(df2)
#unique(df2$frag.type)
table(df2$tx.type)
df2$log2cpm <- log2(df2$cpm+1)
df2$group <- df2$group1
table(df2$cell,df2$group)
df2$cell <- factor(df2$cell,cellShort)
#df2$group <- factor(df2$group,groupShort)

df2.cell <- df2[as.character(df2$cell) %in% c("HNSC"),] # opt1:
# df2.cell <- df2[as.character(df2$cell) %in% c("HNSC","LUAD","LIHC"),] # opt2:
p1 <- ggplot(df2.cell,aes(x=log2cpm,y=tx.type,fill=frag.type))+ 
  # facet_grid(tx.source~.,space = "free_y",scale="free_y") +
  facet_grid(tx.source~group,space = "free_y",scale="free_y") +
  geom_bar(position = "fill",width = 0.8,linewidth=0.3, stat = "identity", color="black") +
  #geom_errorbar()+
  scale_x_continuous(name = "Ratio", breaks = seq(0,1,0.5), labels = seq(0,1,0.5) ) +
  xlab("")+
  ylab(paste0(""))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  #ggsci::scale_fill_d3()+
  scale_fill_npg_adaptive(name="Frag Type")+
  # scale_fill_manual(values = c("grey30","steelblue","firebrick","grey90","steelblue2","firebrick2","steelblue1","salmon"))+
  theme_void() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.text = element_text(size= 20,color = "black"),
        axis.text.x = element_text(size= 20,color = "black",vjust = 0.5), # angle = 90,,hjust = 1
        #panel.grid=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        strip.text = element_text(size= 28,color = "black"),
        panel.border = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = c(rep("grey50",length(other)),rep("steelblue",length(rna)),rep("firebrick",length(dna)))), #  
        legend.position = "right",#c(.25,.6),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))
#ggsave("./unitas_frag_reads_ratio.pdf",width = 8,height = 12)
p2 <- ggplot(df2.cell,aes(x=log2cpm,y=tx.type))+ 
  # facet_grid(tx.source~.,space = "free_y",scale="free_y") +
  facet_grid(tx.source~group,space = "free_y",scale="free_y") +
  geom_bar(width = 0.8,linewidth=0.3, stat = "identity", fill="grey50", color="grey50") + # position = "fill",
  #geom_errorbar()+
  scale_x_continuous(name = "log2CPM" ) +
  xlab("")+
  ylab(paste0(""))+
  geom_vline(xintercept = 0,linetype="dashed",color="grey10")+
  #ggsci::scale_fill_d3()+
  # scale_x_continuous(trans = log10) + 
  scale_fill_npg_adaptive(name="Frag Type")+
  # scale_fill_manual(values = c("grey30","steelblue","firebrick","grey90","steelblue2","firebrick2","steelblue1","salmon"))+
  theme_void() + 
  theme(plot.title = element_text(size = 24,color="black",hjust = 0.5,face="bold"),
        axis.title = element_text(size = 24,color ="black"), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size= 20,color = "black",vjust = 0.5), # angle = 90,,hjust = 1
        #panel.grid=element_blank(),
        panel.grid.major=element_blank(),
        # panel.grid.minor.x = element_blank(),
        # panel.grid.major.y = element_line(color = "grey50",linetype = "dashed"), #size= 1,
        #panel.grid.minor.y = element_blank(),
        strip.text = element_text(size= 28,color = "black"),
        panel.border = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = c(rep("grey50",length(other)),rep("steelblue",length(rna)),rep("firebrick",length(dna)))), #  
        legend.position = "right",#c(.25,.6),
        legend.text = element_text(size= 16),
        legend.title= element_text(size= 16))
p <- ggpubr::ggarrange(plotlist = list(p1,p2),common.legend = TRUE,align = "h", widths=rep(1,5), heights= c(0.16,0.1),legend = "none") # ncol = 5,nrow=2,

ggsave(plot = p, filename = "./figure/unitas_frag_reads.png",width = 24,height = 12)



