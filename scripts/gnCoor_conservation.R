setwd("/lulabdata3/baopengfei/projects/motif-RBP-EDA/")
options(scipen = 4,stringsAsFactors = F,digits = 5,stringsAsFactors = F)
options(bedtools.path = "/lulabdata3/baopengfei/biosoft")
library(bedtoolsr)
library(GenomicRanges)

source("/lulabdata3/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/util.R")
#



dst <- "TCGA_small16"
dedup <- "_all"
pre <- "/lulabdata3/baopengfei/projects/WCHSU-FTC"
# dst <- "SLE"
# dedup <- "_dedup"
# pre <- "/data2/lulab1/bpf/projects/WCHSU-FTC"
outSurfix <- "_rmPurity"
species1 <- "hg38"
extLen <- 75 #0
fromCenter <- T #F
gnDir <- "/lulabdata3/baopengfei/shared_reference/genomes"
outDir <- paste0("output/",dst,"/conservation")
dir.create(outDir)
liftOverVec <- c("panTro5","gorGor5","macFas5","mm10","rn6","susScr11","bosTau8","canFam3","felCat8","monDom5")



# prep genomes in batch ---------
# Define a list of species combinations with their corresponding chain file URLs
# species_combinations <- list(
#   c("hg38", "mm10", paste0("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38To",tocapital("mm10"),".over.chain.gz")),
#   c("hg38", "rn6", "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/mm10ToHg19.over.chain.gz"),
#   # Add more species combinations here as needed
# )

# Create a temporary directory to store the chain files
#temp_dir <- tempdir()
chain_dir <- "/lulabdata3/baopengfei/shared_reference/ucsc/liftover"
# Loop through the species combinations and download the chain files
liftOverVec <- c("panTro5","gorGor5","macFas5","mm10","rn6","susScr11","bosTau8","canFam3","felCat8","monDom5")

for (species2 in liftOverVec) {
  species1 <- "hg38" #spe[1]
  # species2 <- "rn6" # spe[2]
  message(species2)
  chain_url <- paste0("https://hgdownload.soe.ucsc.edu/goldenPath/",species1,"/liftOver/",species1,"To",tocapital(species2),".over.chain.gz")

  chain_file <- file.path(chain_dir,basename(chain_url))
  if(!file.exists(chain_file)){
    message("skip ",chain_file)
  }else{
    download.file(chain_url, chain_file)
  }

}
#


# run liftOver on given bed --------------------
setwd("/lulabdata3/baopengfei/projects/motif-RBP-EDA/")
# Set the source and target genome assemblies
liftOverBed12ToBed6 <- function(bed12){
  ## 1 to 1 # take wider one each record
  # non-dup allowed in V4
  # only keep paste0("chr",c(1:22,"X","Y","M"))
  
  # Create a GRanges object with the coordinates you want to convert
  # bed12 <- sncGnBed12
  keep.chr <- paste0("chr",c(1:22,"X","Y","M"))
  #coords <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(start = 10000, end = 20000))
  #coords <- rtracklayer::import.bed(bed12)
  coords <- bedtoolsr::bt.bed12tobed6(i = bed12 )
  #coords[1:3,]
  #dim(coords)
  #summary(coords$V3-coords$V2)
  coords <- coords[order(coords$V3-coords$V2,decreasing = T),]
  coords <- coords[coords$V1 %in% keep.chr,] # only kept common chr
  #table(coords$V1)
  #summary(coords$V3-coords$V2)
  coords <- coords[!duplicated(coords$V4),] # split bed12 and take longer one each id
  #head(coords$blocks)
  #coords <- coords$blocks
  colnames(coords) <- c("seqnames","start","end","name","score","strand")
  rownames(coords) <- coords$name
  return(coords) # bed6 df (0-base)
}
getLiftOverGnGr <- function(bed6,species1,species2,chain_dir="/lulabdata3/baopengfei/shared_reference/ucsc/liftover"){
  #species1 <- "hg38"
  #species2 <- "rn6"
  # bed6 <- gnBed6
  # chain_dir <- "/lulabdata3/baopengfei/shared_reference/ucsc/liftover"
  chain <- rtracklayer::import.chain(file.path(chain_dir,paste0(species1,"To",tocapital(species2),".over.chain"))) # fail for .gz
  colnames(bed6) <- c("seqnames","start","end","name","score","strand")
  rownames(bed6) <- bed6$name
  #coords$start <- bed6$start+1 # 0- to 1-based
  bed6 <- GenomicRanges::makeGRangesFromDataFrame(bed6, starts.in.df.are.0based=TRUE)
  lifted_coords <- rtracklayerø::liftOver(bed6, chain = chain) # no much param like shell cmd line... (eg. -bedPlus=N)
  return(lifted_coords) # 1-base GR
}
getExtBed6 <- function(bed6,chrSizePath,extLen,fromCenter=T){
  # bed6 <- tmp
  # chrSizePath <- gnChrSize2
  # extLen <- 50
  colnames(bed6) <- c("seqnames","start","end","name","score","strand")
  if(fromCenter){
    bed6$start <- as.integer(0.5*(bed6$start+bed6$end))-extLen
    bed6$start[bed6$start<0] <- 0
    bed6$end <- bed6$start+1
    bed6 <- bedtoolsr::bt.slop(i = bed6, g = chrSizePath, l=0, r = 2*extLen-1)
  }else{
    bed6 <- bedtoolsr::bt.slop(i = bed6, g = chrSizePath, l=extLen, r = extLen-1)
  }

  # bed6[bed6$V4=="peak_26353",]
  #colnames(bed6) <- colnames(gnBed6)
  #rownames(tmp.ext2) <- tmp.ext2$name
  colnames(bed6) <- c("seqnames","start","end","name","score","strand")
  rownames(bed6) <- bed6$name
  return(bed6)
}



# prep hg38 bed6 & fa
sncGnBed12 <- paste0(pre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN/b5_d50_p1_gn.bed")
gnChrSize1 <- paste0(gnDir,"/",species1,"/",species1,".chrom.sizes")
gnFa1 <- paste0(gnDir,"/",species1,"/",species1,".fa")
gnBed6 <- liftOverBed12ToBed6(bed12 = sncGnBed12)
#gnBed6[1:3,]
gnBed6 <- getExtBed6(gnBed6,chrSizePath = gnChrSize1, extLen = extLen, fromCenter=fromCenter)
bedtoolsr::bt.getfasta(fi = gnFa1, fo = paste0(outDir,"/",species1,".fa"), nameOnly = T,bed = gnBed6,split = T,s=T) # need pre-index fa

# prep liftOver bed6 & fa
liftCoordList <- list()
for(species2 in liftOverVec){
  # species2 <- "rn6"
  message(species2)
  # Perform the liftover
  liftCoordList[[paste(species1,species2)]] <- getLiftOverGnGr(bed6 = gnBed6, species1 = species1, species2=species2, chain_dir = "/lulabdata3/baopengfei/shared_reference/ucsc/liftover")
}

for(species2 in liftOverVec){
  # species2 <- "rn6"
  message(species2)
  gnChrSize2 <- paste0(gnDir,"/",species2,"/",species2,".chrom.sizes")
  gnFa2 <- paste0(gnDir,"/",species2,"/",species2,".fa")
  
  tmp <- as.data.frame(liftCoordList[[paste(species1,species2)]]) # 1-based
  tmp <- tmp[order(tmp$width,decreasing = T),]
  tmp <- tmp[!duplicated(tmp$group_name),]
  tmp$score <- 1
  tmp <- tmp[,c("seqnames","start","end","group_name","score","strand")]
  tmp$start <- tmp$start-1 # 0-based
  colnames(tmp)[4] <- "name"
  tmp <- na.omit(tmp)
  #table(duplicated(tmp$name))
  tmp.ext <- getExtBed6(bed6=tmp,chrSizePath=gnChrSize2,extLen=extLen,fromCenter=fromCenter)
  #table(tmp.ext$name=="peak_26353")
  rownames(tmp.ext) <- tmp.ext$name
  table(gnBed6$name %in% tmp.ext$name )
  #table(tmp.ext$name=="peak_26353")
  # all(tmp.ext$name==gnBed6$name) # sanity check 
  
  ## get sequence (compare identity)
  bedtoolsr::bt.getfasta(fi = gnFa2, fo = paste0(outDir,"/",species2,".fa"), nameOnly = T, bed = tmp.ext, split = T, s=T) # need pre-index fa
}
#




# run msa ------------------------------------------------
#library(adegenet)
#library(ape)
#library(ggtree)
library(ggplot2)
library(stats)
#library(ips)
library(msa)
library(seqinr)
library(Biostrings)

# Function to calculate pairwise percentage identity
pairwise_pid <- function(seq1, seq2) {
  alignment <- pairwiseAlignment(seq1, seq2)
  pid(alignment, type = "PID1")
}

# Main function to calculate the identity matrix
calculate_identity_matrix <- function(sequences) {
  # Perform multiple sequence alignment
  aligned <- msa(sequences)
  aligned_sequences <- DNAStringSet(aligned)
  
  # Initialize the identity matrix
  n <- length(aligned_sequences)
  identity_matrix <- matrix(NA, nrow = n, ncol = n)
  
  # Get all unique pairs of indices
  pairs <- combn(seq_along(aligned_sequences), 2, simplify = TRUE)
  
  # Calculate percentage identity for each unique pair using apply
  pairwise_pids <- apply(pairs, 2, function(pair) {
    i <- pair[1]
    j <- pair[2]
    pid_value <- pairwise_pid(aligned_sequences[[i]], aligned_sequences[[j]])
    list(i = i, j = j, pid_value = pid_value)
  })
  
  # Fill the identity matrix with the calculated pid values
  for (result in pairwise_pids) {
    identity_matrix[result$i, result$j] <- result$pid_value
    identity_matrix[result$j, result$i] <- result$pid_value
  }
  
  # Set the diagonal to 100% identity (a sequence is identical to itself)
  diag(identity_matrix) <- 100
  
  # Set row and column names
  colnames(identity_matrix) <- names(aligned_sequences)
  rownames(identity_matrix) <- names(aligned_sequences)
  
  return(identity_matrix)
}

#mySequenceFile <- system.file("examples", "exampleAA.fasta", package="msa")
mySequences <- list()
for(i in c("hg38",liftOverVec)){
  message(i)
  # tmp1 <- (seqinr::read.fasta(paste0(outDir,"/",i,".fa")))
  # tmp <- as.data.frame(readDNAStringSet(paste0(outDir,"/",i,".fa"),format = "fasta",use.names=T))
  
  tmp = Biostrings::readDNAStringSet(paste0(outDir,"/",i,".fa"))
  name = names(tmp)
  seq = paste(tmp)
  tmp <- data.frame(seq,name)
  
  # colnames(tmp) <- "seq"
  tmp$sp <- i
  tmp$path <- paste0(outDir,"/",i,".fa")
  # tmp$name <- rownames(tmp)
  mySequences[[i]] <- tmp
}
#all(rownames(mySequences[["hg38"]])==rownames(mySequences[["monDom5"]]))
#mySequences

#myFirstAlignmentSeq <- unlist(sapply(mySequences,"[",1)) #c(mySequences[[species1]][1],mySequences[[species2]][1])
myFirstAlignmentSeq <- as.data.frame(do.call(rbind,mySequences))
length(unique(myFirstAlignmentSeq$name))
myFirstAlignmentSeq$name <- gsub("\\(\\+\\)|\\(-\\)","",myFirstAlignmentSeq$name,perl = T)
length(unique(myFirstAlignmentSeq$name))
table(duplicated(myFirstAlignmentSeq)) # should not dup

# library(tidyverse)
sp.freq <- myFirstAlignmentSeq %>% 
  # dplyr::group_by(name) %>%
  dplyr::mutate(value=1) %>%
  dplyr::select(sp,name,value) %>% 
  tidyr::pivot_wider(id_cols = name,names_from = sp, values_from = value) %>% 
  as.data.frame()
#table(duplicated(sp.freq$name))
rownames(sp.freq) <- sp.freq$name
sp.freq$name <- NULL
table(is.na(sp.freq))
sp.freq[is.na(sp.freq)] <- 0
sp.freq$freq <- rowSums(as.matrix(sp.freq))
hist(sp.freq$freq) # most exist in 11 sp.
sp.freq <- sp.freq[order(sp.freq$freq,decreasing = T),]

## select candidates
#opt1: pan
this.cells <- cellShort  # all, "HNSC","LUAD","LIHC"
featureLab <- c(rna2,dna)
{
    tmp.diff <- read.table( paste0("./output/",dst,"/diff",outSurfix,".diffSig"),header = T,sep = '\t',check.names = F,stringsAsFactors = F)
  dim(tmp.diff)
  tmp <- matrix(as.matrix(tmp.diff[,2:ncol(tmp.diff)]) %in% c("up","down"),nrow = nrow(tmp.diff))
  dim(tmp)
  tmp <- rowSums(tmp)
  table(tmp >= 16)
  tmp.diff <- tmp.diff[tmp >= 16,]
  select.features <- tmp.diff$feature
  
  # tmp.diff2 <- data.table::fread( paste0("./output/",dst,"/","diff.summary"),header = T,sep = '\t',check.names = F,stringsAsFactors = F,data.table = F)
  # #colnames(tmp.diff2)
  # tmp.diff2 <- tmp.diff2[,grepl(paste0( c("feature",this.cells), collapse = "|"),colnames(tmp.diff2),perl = T)]
  # table(tmp.diff2$LUAD)
  # tmp.diff2 <- tmp.diff2[tmp.diff2$LUAD %in% c("Both significant","Metastasis only","Tumorigenesis only"),]
  # select.features <- tmp.diff2$feature
  
  tmp <- unlist(sapply(strsplit(select.features,"|",fixed = T),"[",8))
  table(tmp)
  #old:
  # lncRNA      mRNA     piRNA pri_miRNA      tRNA 
  # 1         3         1         4         6 
  #new:
  # enhancer    lncRNA      mRNA     piRNA pri_miRNA      Rfam      tRNA 
  # 1         1         1         1         4         1         6 
  select.features <- select.features[tmp %in% c(featureLab)]
  length(select.features) # recur.diff>-16:  15
  select.features <-  unlist(sapply(strsplit(select.features,"|",fixed = T),"[",4))

}
#peak_17765: chr3:169764610-169764655  45nt (lncRNA TERC)
#peak_14288: chr14:20343149-20343273 24nt (lncRNA AL355075)
#this.peak <- "peak_124" # peak_16597
#peak.table.tmp <- peak.table[peak.table$gn=="chr3" & abs(peak.table$gn_start-169764610)<=50 & abs(peak.table$gn_end-169764655)<=50,]
#peak.table.tmp <- peak.table[peak.table$gn=="chr14" & abs(peak.table$gn_start-20343115)<=50 & abs(peak.table$gn_end-20343273)<=50,]

#peak_24038
# peak_48574,hsa_piR_020345,chr19:47907760-47907781
# peak_28885,LTR16E2__chr3___84299204____84299581_neg,chr3:84299383-84299399 (-)

#rownames(peak.table) <- peak.table$name
peak.table.tmp <- peak.table[select.features,]
this.peaks <- c("peak_17765","peak_14288",select.features )
for(this.peak in this.peaks){
  # this.peak <- "peak_39963"
  message(this.peak)
  if(!(this.peak %in% myFirstAlignmentSeq$name)){
    message("skip ",this.peak)
    next
  }

  sp.freq[this.peak,"freq"]
  {
    #TERC: chr3:169764450-169764876 peak_16598+peak_16597
    seqs <- myFirstAlignmentSeq$seq[myFirstAlignmentSeq$name==this.peak]
    names(seqs) <- myFirstAlignmentSeq$sp[myFirstAlignmentSeq$name==this.peak]

    ## predefined order
    myFirstAlignment <- msa::msa(seqs, 
                                 type = "dna", #order = "input",
                                 method = "ClustalW"
                                 # gapOpening = 12,gapExtension = 3
                                 )   # , "ClustalOmega", "Muscle"
    ## use default substitution matrix
    #myFirstAlignment
    #print(myFirstAlignment, show="complete")
    
    ## op1: plot multi-align
    # msaPrettyPrint(myFirstAlignment, output="pdf", showNames="none",
    #                showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
    # msaPrettyPrint(myFirstAlignment, output="pdf", y=c(164, 213),
    #                subset=c(1:6), showNames="none", showLogo="top",
    #                logoColors="rasmol", shadingMode="similar",
    #                showLegend=FALSE, askForOverwrite=FALSE)
    alignFa <- paste0(outDir,"/",this.peak,".fasta")
    outFig <- paste0(outDir,"/",this.peak,".seq.pdf")
    tmpFig <- paste0("./",this.peak,".seq.pdf")
    msaPrettyPrint(myFirstAlignment, output=c("pdf"), file=outFig, alFile=alignFa, # save .fasta not .fa, otherwise got error.
                   # y=c(164, 213),
                   # paperWidth=11, paperHeight=8.5,
                   paperWidth=2+max(width(unmasked(myFirstAlignment)))*0.1,paperHeight=0.1+length(myFirstAlignment@unmasked@ranges*0.1),
                   showNames="left", shadingMode="functional",
                   shadingModeArg="structure",
                   askForOverwrite=FALSE)
    file.copy(from = tmpFig, to   = outFig, overwrite = T)
    file.remove(tmpFig)

    #opt1: identity score
    tmp <- calculate_identity_matrix(Biostrings::DNAStringSet(myFirstAlignment))

    #opt2: identity score (seem not correct)
    # myFirstAlignment2 <- msaConvert(myFirstAlignment, type="seqinr::alignment")
    # d <- dist.alignment(myFirstAlignment2, matrix="identity",gap=F) # for protein: c("identity", "similarity"), for DNA/RNA: identity (and seem not correct)
    # tmp <- as.matrix(d)  # distance matrix
    
    #hsa-miR-148b
    #dist.mat <- distance(birds.t, method = "jaccard")
    #true.dist.mat <- as.dist(dist.mat)
    #clust.res <- hclust(true.dist.mat, method = "complete")
    
    
    #d <- colnames(d)
    #library(ape)
    #hemoTree <- ape::nj(d)  # Neighbor-Joining Tree Estimation, d may be an object of class “dist”
    # hemoTree2 <- stats::hclust(d,method = "complete")  # method: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"
    hemoTree2 <- stats::hclust(as.dist(tmp),method = "complete")  # method: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"
    
    #plot with ggdendro
    library("ggplot2")
    library("ggdendro")
    ggdendrogram(as.dendrogram(hemoTree2), rotate = TRUE, theme_dendro = FALSE)+
      ggsci::scale_color_nejm()+
      theme_classic() +
      xlab(this.peak) + ylab("") +
      theme_minimal() + 
      theme(plot.title = element_text(size = 28,color="black",hjust = 0.5,face="bold"),
            axis.title = element_text(size = 28,color ="black",face="bold"),
            axis.text = element_text(size= 24,color = "black"), #,face="bold
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.x = element_text( hjust = 0.5 ), # angle = 45,
            axis.text.y = element_text( hjust = 1 ), # angle = 45,
            panel.grid=element_blank(),
            #legend.position = "top",
            legend.text = element_text(size= 26,color = "black",face="bold"),
            legend.title= element_text(size= 26,color = "black",face="bold"))
    width <- 7 #2+max(width((unmasked(myFirstAlignment))))*0.1
    height <- 2+length(unmasked(myFirstAlignment))*0.5
    ggsave(filename = paste0(outDir,"/",this.peak,".tree.pdf"),
           width = width,height = height,limitsize = F
    )
    #
    
    
    ## seq.identity.cluster order
    myFirstAlignment <- msa::msa(seqs[rev(hemoTree2$labels[hemoTree2$order])], 
                                 type = "dna", 
                                 method = "ClustalW", order = "input",
                                 # gapOpening = 12,gapExtension = 3
    )   # , "ClustalOmega", "Muscle"
    alignFa <- paste0(outDir,"/",this.peak,".fasta")
    outFig <- paste0(outDir,"/",this.peak,".seq.reorder.pdf")
    tmpFig <- paste0("./",this.peak,".seq.reorder.pdf")
    msaPrettyPrint(myFirstAlignment, output=c("pdf"), file=outFig, alFile=alignFa, # save .fasta not .fa, otherwise got error.
                   # y=c(164, 213),
                   # paperWidth=11, paperHeight=8.5,
                   paperWidth=2+max(width(unmasked(myFirstAlignment)))*0.1,paperHeight=0.1+length(myFirstAlignment@unmasked@ranges*0.1),
                   showNames="left", shadingMode="functional",
                   shadingModeArg="structure",
                   askForOverwrite=FALSE)
    file.copy(from = tmpFig, to   = outFig, overwrite = T)
    file.remove(tmpFig)
    
    # ## op2: ggmsa need R>4.1
    # library(ggmsa)
    # extLen <- 50
    # #protein_sequences <- system.file("extdata", "sample.fasta", package = "ggmsa")
    # myFirstAlignment2 <- readDNAMultipleAlignment(alignFa)
    # p <- ggmsa( myFirstAlignment2,
    #       font = "helvetical",
    #       color = "Chemistry_NT",custom_color = NULL,by_conservation=T,
    #       # position_highlight = c(40,60),ignore_gaps = FALSE,
    #       # start = 0, end = extLen*2,
    #       # consensus_views	=T,use_dot	=F,
    #       char_width = 0.5, seq_name = T, none_bg=F
    #       ) + 
    #   geom_seqlogo() +geom_msaBar()
    # width <- 2+max(width(unmasked(myFirstAlignment2)))*0.1
    # height = 0.1+length(myFirstAlignment2@unmasked@ranges*0.1)
    # ggsave(plot = p, filename = paste0(outDir,"/",this.peak,".align.pdf"),
    #        width = width,height = height,limitsize = F
    #        )
  }
}
#cat ./output/TCGA_small7/diff_LUAD_YesvsNo_rmPurity.diff | cut -f1 | grep "peak_13512"





# plot conservation matrix heatmap ----------------------------------------
# get all peak identity mat identity & jaccard heatmap
getMSAidentityMat <- function(idx,seqs.list,type="dna",method="ClustalW"){
  #seqs # named vec of seq strings
  # idx <- 24
  library(seqinr)
  library(msa)
  message(idx)
  seqs <- seqs.list[[idx]]
  myFirstAlignment <- msa::msa(seqs, 
                               type = type, #order = "input",
                               method = method
                               # gapOpening = 12,gapExtension = 3
  )  
  
  ## prepare tree
  #hemoSeq <- readAAStringSet(system.file("examples/HemoglobinAA.fasta",package="msa"))
  #hemoAln <- msa(hemoSeq)
  #hemoAln2 <- msaConvert(hemoAln, type="seqinr::alignment")
  
  ##opt1:
  # myFirstAlignment2 <- msaConvert(myFirstAlignment, type="seqinr::alignment")
  # d <- dist.alignment(myFirstAlignment2, matrix = "identity") #  matrix = c("identity", "similarity")  no difference in real ?
  # tmp <- as.data.frame(as.matrix(d))  # distance matrix
  
  ##opt2:
  tmp <- calculate_identity_matrix(Biostrings::DNAStringSet(myFirstAlignment))
  tmp <- as.data.frame(tmp)  # distance matrix
  
  # similarity is the likeness (resemblance) between two sequences in comparison 
  #while identity is the number of characters that match exactly between two different sequences.
  # tmp$id <- names(seqs.list)[[idx]]
  return(tmp)
}
flattern_triMat <- function(m) {
  ut <- upper.tri(m,diag = T)
  data.frame(from = rownames(m)[row(m)[ut]],
             to = rownames(m)[col(m)[ut]],
             value=t(m)[ut]
             # p=t(m)[ut]
             )
}
duplicated(c(1:4,2:5))
table(duplicated(myFirstAlignmentSeq$name)) # rm uniq rows
dup.name <- unique(myFirstAlignmentSeq$name[duplicated(myFirstAlignmentSeq$name)])
myFirstAlignmentSeq <- myFirstAlignmentSeq[myFirstAlignmentSeq$name %in% dup.name,]
#myFirstAlignmentSeq <- myFirstAlignmentSeq[,]


## select candidates
#opt1: pan
this.cells <- cellShort  # all, "HNSC","LUAD","LIHC"
featureLab <- c(rna2,dna)
{
  tmp.diff <- read.table( paste0("./output/",dst,"/diff",outSurfix,".diffSig"),header = T,sep = '\t',check.names = F,stringsAsFactors = F)
  dim(tmp.diff)
  tmp <- matrix(as.matrix(tmp.diff[,2:ncol(tmp.diff)]) %in% c("up","down"),nrow = nrow(tmp.diff))
  dim(tmp)
  tmp <- rowSums(tmp)
  table(tmp >= 6)
  tmp.diff <- tmp.diff[tmp >= 6,]
  select.features <- tmp.diff$feature
  
  # tmp.diff2 <- data.table::fread( paste0("./output/",dst,"/","diff.summary"),header = T,sep = '\t',check.names = F,stringsAsFactors = F,data.table = F)
  # #colnames(tmp.diff2)
  # tmp.diff2 <- tmp.diff2[,grepl(paste0( c("feature",this.cells), collapse = "|"),colnames(tmp.diff2),perl = T)]
  # table(tmp.diff2$LUAD)
  # tmp.diff2 <- tmp.diff2[tmp.diff2$LUAD %in% c("Both significant","Metastasis only","Tumorigenesis only"),]
  # select.features <- tmp.diff2$feature
  
  tmp <- unlist(sapply(strsplit(select.features,"|",fixed = T),"[",8))
  table(tmp)
  select.features <- select.features[tmp %in% c(featureLab)]
  tmp <- unlist(sapply(strsplit(select.features,"|",fixed = T),"[",8))
  table(tmp)
  length(select.features) # recur.diff>-16:  15
  select.features <-  unlist(sapply(strsplit(select.features,"|",fixed = T),"[",4))
  #table(is.na(select.features))
  select.features <- select.features[select.features %in% unique(c(myFirstAlignmentSeq$name,peak.table$name))]
  # lncRNA      mRNA     piRNA pri_miRNA      tRNA 
  # 1         3         1         4         6 
  select.features.list <- list()
  for(i in unique(tmp)){
    select.features.list[[i]]  <- select.features[tmp == i]
  }
  
}
#peak_17765: chr3:169764610-169764655  45nt (lncRNA TERC
#peak_14288: chr14:20343149-20343273 24nt (lncRNA AL355075)
#this.peak <- "peak_124" # peak_16597
#peak.table.tmp <- peak.table[peak.table$gn=="chr3" & abs(peak.table$gn_start-169764610)<=50 & abs(peak.table$gn_end-169764655)<=50,]
#peak.table.tmp <- peak.table[peak.table$gn=="chr14" & abs(peak.table$gn_start-20343115)<=50 & abs(peak.table$gn_end-20343273)<=50,]



## by each RNA topN conserved & diff
# outSurfix <- "_rmMonthOperator"
# diffSig <- read.table(paste0("./output/",dst,"/cfPeakCNN_smallDomain_diff",outSurfix,".diffSig"),header = T,sep = '\t')(
#dim(diffSig)
c(rna,dna)
outpre="/lulabdata3/baopengfei/projects/WCHSU-FTC"
#outpre="/data2/lulab1/bpf/projects/WCHSU-FTC"
for(this.RNA in c("mRNA","lncRNA","tRNA","enhancer","promoter","ktRNA","pri_miRNA","piRNA","repeats")){ # rRNA: error
# this.RNA <- "pri_miRNA"
  message(this.RNA)
  
  # peak.table <- data.table::fread(paste0(outpre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN/b5_d50_p1.summary"),header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",data.table = F)
  peak.table.tmp <- peak.table[order(peak.table$phastCons100way,decreasing = T),]
  peak.table.tmp <- peak.table.tmp[peak.table.tmp$gn_type==this.RNA,]
  peak.table.tmp <- peak.table.tmp[select.features.list[[this.RNA]],c("name","phastCons100way","MFE_FDR", "coding_probability")] # 1 tx_type
  peak.table.tmp <- peak.table.tmp[1:min(50,nrow(peak.table.tmp)),]
  table(peak.table.tmp$tx_type)
  this.peaks <- unique(myFirstAlignmentSeq$name)
  this.peaks <- this.peaks[this.peaks %in% peak.table.tmp$name]
  #this.peaks <- this.peaks[sample(1:length(this.peaks),250,replace = F)]
  #this.peaks
  seqs.list <- list()
  for(this.peak in this.peaks){
    # message(this.peak)
    # sp.freq[this.peak,"freq"]
    seqs <- myFirstAlignmentSeq$seq[myFirstAlignmentSeq$name==this.peak]
    names(seqs) <- myFirstAlignmentSeq$sp[myFirstAlignmentSeq$name==this.peak]
    seqs.list[[this.peak]] <- seqs
  }
  
  library(parallel)
  # tmp.list <- parallel::mclapply(1:length(seqs.list), getMSAidentityMat, seqs.list=seqs.list,mc.cores = 6) # both opt: error...
  tmp.list <- lapply(1:length(seqs.list), FUN = getMSAidentityMat, seqs.list=seqs.list)
  
  names(tmp.list) <- names(seqs.list)
  #names(tmp)
  tmp2 <- lapply(tmp.list,flattern_triMat)
  tmp2 <- as.data.frame(do.call(rbind,tmp2))
  tmp2$id <- splitVec(rownames(tmp2),delim = ".",fixed = T,nth = 1)
  table(tmp2$from==tmp2$to)
  if(max(tmp2$value)>1 & max(tmp2$value)<=100 ){
    tmp2$value <- tmp2$value/100
  }
  tmp2$value[tmp2$from==tmp2$to] <- 1
  table(tmp2$from)
  table(tmp2$to)
  for(i in 1:nrow(tmp2)){
    # i <- 2
    # tmp2[i,]
    tmp2[i,c("from","to")] <- sort(tmp2[i,c("from","to")])
  }
  
  
  # seq ident heatmap: focus on hg38
  select.sp <- "hg38"
  tmp3 <- tmp2[tmp2$from==select.sp | tmp2$to==select.sp,]
  hist(tmp3$value)
  tmp3$from2 <- tmp3$from
  tmp3$to2 <- tmp3$to
  tmp3$from2[tmp3$from!=select.sp] <- select.sp
  tmp3$to2[tmp3$from!=select.sp] <- tmp3$from[tmp3$from!=select.sp]
  table(tmp3$from2) # only hg38
  table(tmp3$to2)
  
  cons.tbl <- dplyr::as_tibble(tmp3) %>% 
    dplyr::select(-c("from","to")) %>% 
    tidyr::pivot_wider(id_cols = id,names_from = to2,values_from = value) %>% 
    as.data.frame()
  cons.tbl <- cons.tbl[,c("id",select.sp,liftOverVec)]
  data.table::fwrite(x = cons.tbl,file = paste0(outDir,"/",this.RNA,"_seq_identity_mat.txt"),quote = F,row.names = F,col.names = T,sep = '\t')
  
  rownames(cons.tbl) <- cons.tbl$id
  cons.tbl$id <- NULL
  cons.tbl <- cons.tbl[,c(select.sp,liftOverVec)]
  
  brkList <- seq(0, 1, by = 0.01)
  # brkList <- seq(-1, 1, by = 0.1)
  #min(cor.mat)
  i <- "sp"
  colorList <- list()
  colorList[[i]] <- pal_npg_adaptive()( length( unique( c(select.sp, liftOverVec) ) ) )
  names(colorList[[i]]) <- unique( c(select.sp, liftOverVec) )
  colorList[["RNA"]] <- RNA_colors[c(rna2,dna)]

  #peak.table.tmp <- data.table::fread(paste0(outpre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN/b5_d50_p1.summary"),header = T,sep = "\t",stringsAsFactors = F,check.names = F,quote = "",data.table = F)
  #peak.table.tmp[1:3,1:4]
  rownames(peak.table.tmp) <- peak.table.tmp$name
  peak.table.tmp <- peak.table.tmp[rownames(cons.tbl),]
  dim(peak.table.tmp)
  #colnames(peak.table.tmp)
  peak.table.tmp <- peak.table.tmp[,c("MFE_FDR","coding_probability","phastCons100way")] # 
  # colnames(peak.table.tmp)[1] <- "RNA"
  # peak.table.tmp$RNA <- NULL
  pheatmap::pheatmap(mat = cons.tbl,
                     # annotation_col = sample.table, #data frame that specifies the annotations shown on left side of the heatmap.
                     annotation_row = peak.table.tmp,
                     annotation_colors = colorList, # RColorBrewer::brewer.pal(6, "Set3")
                     border_color="white",na_col = "white",
                     scale = "none", # row,none
                     #labels_col = 3, labels_row = 6,
                     cluster_cols = F,cluster_rows = T,
                     # gaps_row = sum(annotation_row$trend=="up"), # invalid if cluster row/col
                     #cutree_cols = 2,cutree_rows = 3,
                     show_colnames=T, show_rownames=F,
                     fontsize = 10,
                     legend =  T, 
                     height = 7, width = 6,
                     color = ( colorRampPalette( RColorBrewer::brewer.pal(n = 9, name = "OrRd"))(length(brkList)) ) , # ),  OrRd
                     breaks = brkList,
                     #fontsize_row = 5,
                     filename = paste0(outDir,"/",this.RNA,"_seq_identity_mat.png")
  )
  #
  
  
  
  # all sp. jaccard heatmap (with filter) 
  #(2024coPARSE: BLAST v2.12.0 bl2seq (E value < 10−4, hit length >50 nt, overall sequence identity >50%))
  tmp4 <- tmp2
  ident.cutOff <- 0.9 #not ext: 0.9 #mean(tmp2$value,na.rm=T)
  tmp4[1:3,]
  tmp4$value[tmp4$value>=ident.cutOff] <- 1
  tmp4$value[tmp4$value<ident.cutOff] <- 0
  table(duplicated(paste0(tmp4$from,tmp4$to,tmp4$id)))
  
  # tmp4$from <- factor(tmp4$from,levels = c(select.sp,liftOverVec))
  # tmp4$to <- factor(tmp4$to,levels = c(select.sp,liftOverVec))
  cons.tbl2 <- dplyr::as_tibble(tmp4) %>% 
    # dplyr::
    dplyr::group_by(from,to) %>%
    dplyr::summarise(frac=sum(value)/dplyr::n()) %>% 
    tidyr::pivot_wider(names_from = to,values_from = frac) %>%
    as.data.frame()
  #summary(cons.tbl2$frac)
  data.table::fwrite(x = cons.tbl2,file = paste0(outDir,"/",this.RNA,"_seq_jaccard_mat.txt"),quote = F,row.names = T,col.names = T,sep = '\t')
  
  
  #get full-matrix from upper half tri mat
  rownames(cons.tbl2) <- cons.tbl2$from
  cons.tbl2$from <- NULL
  liftOverVecAll <- c(select.sp, liftOverVec)
  for (i in liftOverVecAll){
    for(j in liftOverVecAll[liftOverVecAll!=i]){
      if(is.na(cons.tbl2[j,i])){
        cons.tbl2[j,i] <- cons.tbl2[i,j] 
      }
    }
  }
  cons.tbl2 <- cons.tbl2[liftOverVecAll,liftOverVecAll]
  
  brkList <- seq(0, 1, by = 0.01)
  # brkList <- seq(-1, 1, by = 0.1)
  #min(cor.mat)
  i <- "sp"
  colorList <- list()
  colorList[[i]] <- pal_npg_adaptive()( length( unique( c(select.sp, liftOverVec) ) ) )
  names(colorList[[i]]) <- unique( c(select.sp, liftOverVec) )
  pheatmap::pheatmap(mat = cons.tbl2,
                     # annotation_col = sample.table, #data frame that specifies the annotations shown on left side of the heatmap.
                     # annotation_row = sample.table,
                     annotation_colors = colorList, # RColorBrewer::brewer.pal(6, "Set3")
                     border_color="white",na_col = "white",
                     scale = "none", # row,none
                     #labels_col = 3, labels_row = 6,
                     cluster_cols = F,cluster_rows = F,
                     # gaps_row = sum(annotation_row$trend=="up"), # invalid if cluster row/col
                     #cutree_cols = 2,cutree_rows = 3,
                     show_colnames=T, show_rownames=T,
                     fontsize = 10,
                     legend =  T, 
                     height = 2.6, width = 3,
                     color =  ( colorRampPalette( RColorBrewer::brewer.pal(n = 9, name = "OrRd"))(length(brkList)) ) , # ),  OrRd, rev(RdYlBu)
                     breaks = brkList,
                     #fontsize_row = 5,
                     filename = paste0(outDir,"/",this.RNA,"_seq_jac_mat.png")
  )
}
# peak.table <- peak.table[!grepl(",",peak.table$gn_blockSizes),]
# peak.table <- peak.table[,7:12]
# colnames(peak.table) <- c("seqnames","start","end","name","score","strand")
# 
# signal.gr <- rtracklayer::import.bw(con=as.character("/lulabdata3/baopengfei/shared_reference/ucsc/conservation/hg38.phastCons100way.bw"),which=GenomicRanges::makeGRangesFromDataFrame(peak.table)) # import(.bed) got error if input bed not proper formated

#



# optinal: plot conservation deltaPhastCons (peak-flank) score --------------------
#should be reasonable if merged into snc_cluster (<20? bp)

dst <- "TCGA_small16"
dedup <- "_all"
pre <- "/lulabdata3/baopengfei/projects/WCHSU-FTC"
#peak.table <- data.table::fread(paste0(pre,"/output/",dst,"/call_peak",dedup,"/cfpeakCNN/b5_d50_p1.summary"),sep='\t',header = T,data.table = F,check.names = F,stringsAsFactors = F)
#rownames(peak.table) <- peak.table$name
peak.table$deltaCons <- peak.table$phastCons100way-peak.table$phastCons100way_flank
summary(peak.table$deltaCons,breaks = 100)

table((peak.table$phastCons100way==0) & (peak.table$phastCons100way_flank==0))
peak.table.tmp <- peak.table[(peak.table$phastCons100way!=0) & !is.na(peak.table$phastCons100way),]
peak.table.tmp <- peak.table.tmp[order(peak.table.tmp$deltaCons,decreasing = T),]
peak.table.tmp$rank <- 1:nrow(peak.table.tmp)
peak.table.tmp <- peak.table.tmp[order(peak.table.tmp$deltaCons,decreasing = T),c("name","rank","deltaCons","gn_type","tx")]
peak.table.tmp$gn_type <- factor(peak.table.tmp$gn_type, levels = c(rna2,dna))
table(peak.table.tmp$gn_type)
p1 <- ggplot(peak.table.tmp,aes(x=gn_type,y=deltaCons,fill=gn_type))+  # ,label=label
  geom_boxplot(aes(x=(gn_type),y=deltaCons,fill=gn_type),color="black",alpha=0.6,outlier.size = 0) +
  # geom_point(
  #   # aes(x=(as.numeric(group)+b),y=value,fill=group),
  #            shape=21,
  #            size=2,stroke=0.2,alpha=0.9,color="grey20"
  #   ) +
  scale_fill_manual(values = RNA_colors) +
  # scale_x_discrete(label=c("CT","CN","PT","PN","RT","RN"))+
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "deltaCons")+
  #ylab("ssGSEA score")+
  ylab("deltaCons")+
  xlab("") +
  # ggpubr::stat_compare_means(
  #   comparisons = list( c("HC","LDA"),c("LDA","HDA"),c("HC","HDA") ),
  #   label.x.npc = 0.5,hide.ns=F,size =8, vjust = 0.5, step.increase = 0.04,#paired = TRUE,
  #   aes(label = p.signif), # p.signif, group = Group,
  #   symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1,1),symbols = c( "****","***", "**", "*", "ns")),
  #   method = "wilcox.test") +
  theme_minimal(base_size=12) +
  theme(#axis.ticks.x=element_blank(),
    aspect.ratio = 1,
    # panel.border = element_rect(colour="white", fill = NA), #for panel.border # 
    strip.text.y = element_blank(),
    strip.text.x = element_text(family="sans",size=20), # serif
    axis.title.x = element_text(family="sans",size=20),
    axis.title.y = element_text(family="sans",size=20),
    axis.text.x = element_text(size = 20,family="sans",angle = 45,hjust=1,vjust=1),
    axis.text.y = element_text(size = 20,family="sans"),
    plot.title = element_text(size=20),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
  # facet_wrap(facets = "tx_type", scales = "free_y",ncol = 4 )
ggsave(plot = p1, "./figure/rank_delta_gn.pdf",height = 8, width = 10)


#peak.table.tmp2[1:3,]
colnames(peak.table)
table((peak.table$phastCons100way==0) & (peak.table$phastCons100way_flank==0))
table(is.na(peak.table$phastCons100way))
peak.table.tmp2 <- peak.table[(peak.table$phastCons100way!=0) & !is.na(peak.table$phastCons100way),]
peak.table.tmp2 <- dplyr::as_tibble(peak.table.tmp2) %>% 
  dplyr::select(name,phastCons100way,phastCons100way_flank,gn_type) %>% 
  tidyr::pivot_longer(cols = c("phastCons100way","phastCons100way_flank"),names_to = "region",values_to = "phastCons") %>% 
  as.data.frame()
# peak.table.tmp2$tx_type <- factor(peak.table.tmp2$tx_type, levels = c(dna,rna))
peak.table.tmp2$gn_type <- factor(peak.table.tmp2$gn_type, levels = c(rna2,dna))
table(peak.table.tmp2$gn_type)
library(ggpubr)
p2 <- ggplot(peak.table.tmp2,aes(x=gn_type,y=phastCons,fill=region))+  # ,label=label
  geom_boxplot(color="black",alpha=0.8,outlier.size = 0,outlier.colour = "white") +
  # geom_point(
  #   # aes(x=(as.numeric(group)+b),y=value,fill=group),
  #   shape=21,
  #   size=2,stroke=0.2,alpha=0.9,color="grey20"
  # ) +
  scale_fill_manual(values = c("phastCons100way"="red","phastCons100way_flank"="salmon")) + #RNA_colors
  # scale_x_discrete(label=c("CT","CN","PT","PN","RT","RN"))+
  ylim(c(0,1.2)) +
  labs(title = "phastCons")+
  #ylab("ssGSEA score")+
  ylab("phastCons")+
  xlab("") +
  ggpubr::stat_compare_means(
    # ref.group = "phastCons100way",
    # comparisons = list( c("HC","LDA"),c("LDA","HDA"),c("HC","HDA") ),
    label.x.npc = 0.5,hide.ns=F,size =8, vjust = 0.5, step.increase = 0.04,#paired = TRUE,
    # aes(label = "p.signif"), # p.signif, group = Group,
    label="p.signif",
    symnum.args=list(cutpoints = c(0, 0.0001,0.001, 0.01, 0.05,1),symbols = c( "****","***", "**", "*", "ns")),
    method = "wilcox.test",
    method.args = list(alternative = "less")  #
  ) +
  theme_minimal(base_size=12) +
  theme(#axis.ticks.x=element_blank(),
    aspect.ratio = 0.5,
    # panel.border = element_rect(colour="white", fill = NA), #for panel.border # 
    strip.text.y = element_blank(),
    strip.text.x = element_text(family="sans",size=20), # serif
    axis.title.x = element_text(family="sans",size=20),
    axis.title.y = element_text(family="sans",size=20),
    axis.text.x = element_text(size = 20,family="sans",angle = 45,hjust=1,vjust=1),
    axis.text.y = element_text(size = 20,family="sans"),
    plot.title = element_text(size=20),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
  # facet_wrap(facets = "tx_type", scales = "free_y",ncol = 4 )
ggsave(plot=p2, "./figure/phastCons_flank_box_gn.pdf",height = 8, width = 12)
#




# merge plot
library(cowplot)
my_theme <- theme_minimal(base_size = 12) +
  theme(
    aspect.ratio = NULL,
    strip.text.y = element_blank(),
    strip.text.x = element_text(family = "sans", size = 20),
    axis.title.x = element_text(family = "sans", size = 20),
    axis.title.y = element_text(family = "sans", size = 20),
    axis.text.x  = element_text(size = 20, family = "sans", angle = 45, hjust = 1, vjust = 1),
    axis.text.y  = element_text(size = 20, family = "sans"),
    plot.title   = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# 3) p1
p1 <- ggplot(peak.table.tmp, aes(x = gn_type, y = deltaCons, fill = gn_type)) +
  geom_boxplot(color = "black", alpha = 0.95, outlier.size = 0) +
  scale_fill_manual(values = RNA_colors) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  labs(y = "Delta phastCons", x = NULL, fill = "Precursor") +
  my_theme +
  theme(
    legend.position = "right",
    plot.margin  = margin(0, 5.5, 0, 5.5)   # 底部 margin 设 0
  )

# 4) p2
p2 <- ggplot(peak.table.tmp2, aes(x = gn_type, y = phastCons, fill = region)) +
  geom_boxplot(color = "black", alpha = 0.95, outlier.size = 0, outlier.colour = "white") +
  scale_fill_manual(values = c("phastCons100way" = "red",
                               "phastCons100way_flank" = "salmon")) +
  ylim(c(0, 1.2)) +
  geom_hline(yintercept = 0.5, color = "grey", linetype = "dashed", linewidth = 1) +
  labs(y = "PhastCons", x = "") +
  ggpubr::stat_compare_means(
    label.x.npc = 0.5, hide.ns = FALSE, size = 8, vjust = 0.5, step.increase = 0.04,
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                       symbols = c("****", "***", "**", "*", "ns")),
    method = "wilcox.test",
    method.args = list(alternative = "less")
  ) +
  my_theme +
  theme(
    legend.position = "none",
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x  = element_blank(),
    plot.margin  = margin(5.5, 5.5, 0, 5.5)   #
  )

# 下面的 p1：保留 x 轴
p12 <- plot_grid(
  p2, p1,
  ncol        = 1,
  align       = "v",       
  axis        = "tblr",   
  rel_heights = c(1, 1.2),
  labels      = NULL
) 
ggsave(plot=p12, "./figure/phastCons_flank_box_gn_merge.png",height = 7, width = 12)
#



# cor: DEG Freq. vs PhastCons ----------------
#consensus.diff <- data.table::fread(paste0("/lulabdata3/baopengfei/smb/",dst,"_candidate_diff_summary.txt.tmp"), sep = "\t",header = T, stringsAsFactors = F, check.names = F, data.table = F)
consensus.diff <- data.table::fread(paste0("/lulabdata3/baopengfei/smb/",dst,"_candidate_diff_summary.txt.tmp2"), sep = "\t",header = T, stringsAsFactors = F, check.names = F, data.table = F)
#rownames(peak.table)
#dim(consensus.diff)
#consensus.diff$phastCons100way
tmp1 <- peak.table[consensus.diff$name,c("name","phastCons100way","phastCons100way_flank")]
#tmp1 <- tmp1
cor.test(tmp1$phastCons100way,consensus.diff$unionSigup)
cor.test(tmp1$phastCons100way,consensus.diff$unionSigdown)
df.tmp <- data.frame(name=tmp1$name,phastCons=tmp1$phastCons100way,DES_Freq=consensus.diff$unionSigup)
df.tmp$trend <- "up"
df.tmp1 <- data.frame(name=tmp1$name,phastCons=tmp1$phastCons100way,DES_Freq=consensus.diff$unionSigdown)
df.tmp1$trend <- "down"
df.tmp <- as.data.frame(rbind(df.tmp,df.tmp1))
table(duplicated(df.tmp$name))
ggplot(df.tmp, aes(x=DES_Freq, y=phastCons) ) + # type2, type3
  geom_point() + # fill="grey",color="grey",  minus.log10.padj
  # scale_color_continuous(name="-log10.P",type = "viridis") + # YlOrRd, viridis RColorBrewer::brewer.pal(n=8,name="Reds")
  # scale_size_continuous(name = sizeLab) +
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
    legend.title= element_text(size= 16)) + 
  facet_grid(trend~.)
p <- 
  ggplot(data=df.tmp,aes(x=DES_Freq, y=phastCons))+
  # geom_point(alpha=0.5, color="white",fill="grey30", size=2, shape=21) + # fill="grey",color="grey",
  geom_point( aes(fill=phastCons,alpha=DES_Freq, size=2) ,color="white",shape=21) + # fill="grey",color="grey",
  # geom_bin2d(bins = 200) +
  # geom_hex(bins = 100) +
  # xlim(c(-10,10)) +
  # xlim(c(-6,6)) +
  scale_fill_viridis_c(direction=1,na.value = "grey30" )+ # values=c(0.1,0.4,0.7,1) # option = "E"
  scale_color_nejm_adaptive() +
  geom_smooth(method='lm', formula= y~x,color="grey30") +
  ggpubr::stat_cor(label.x = 1,size=8,label.y = 0.8,method = "spearman")+ #this means at 35th unit in the y axis, the r squared and p value will be shown
  # ggpubr::stat_regline_equation(label.x = 1,size=8,label.y = 0.8) +
  # scale_fill_manual(values = c("grey50","#ED782F")) + #"#ED782F","#B52B1A","#3E8C9D"
  # scale_color_manual(values = c("grey50","#ED782F")) +
  geom_vline(xintercept = 0, color="grey", alpha=0.8, linetype="dashed")+
  geom_hline(yintercept = 0, color="grey", alpha=0.8, linetype="dashed")+
  xlab("DES Frequency") + ylab("phastCons100way") +
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
ggsave(plot = p, filename = paste0("./figure/","phastCon_DES_cor",".png"),device = "png", width = 12,height = 10, dpi = 200 ) # , limitsize = FALSE
#
#