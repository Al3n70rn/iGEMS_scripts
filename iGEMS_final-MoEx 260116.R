rm(list=ls(all=T))

##-------------Install libraries required in R for iGEMS--------------##

source("https://bioconductor.org/biocLite.R")
biocLite("limma")
biocLite("fdrtool")
biocLite("statmod")
biocLite("stringr")
biocLite("parallel")
biocLite("affxparser")
biocLite("affy")
biocLite("aroma.affymetrix")
biocLite("biomaRt")
biocLite("dplyr")
biocLite("data.table")

## FIRMAGene is not available through BioConductor or CRAN, to install it please follow these instructions
# 1. Download FIRMAGene R package from the link https://r-forge.r-project.org/R/?group_id=424
# 2. Open the terminal and change to the directory where the package is downloaded eg: cd Downloads
# 3. Use the command R CMD INSTALL FIRMAGene_0.9.8.tar to install the package

#*Please ignore R install command listed on FIRMAgene website at it does not work and gives error.

##-------------end of installation instructions-------------------##


##-------sessionInfo() for the run that gives the version of all the packages used--------#

# R version 3.2.3 (2015-12-10)
# Platform: x86_64-apple-darwin13.4.0 (64-bit)
# Running under: OS X 10.11.1 (El Capitan)
# 
# locale:
#   [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] data.table_1.9.6        dplyr_0.4.3             stringr_1.0.0           statmod_1.4.22          biomaRt_2.26.1         
# [6] limma_3.26.3            fdrtool_1.2.15          FIRMAGene_0.9.8         aroma.light_3.0.0       aroma.affymetrix_2.14.0
# [11] aroma.core_2.14.0       R.devices_2.13.2        R.filesets_2.9.0        R.utils_2.2.0           R.oo_1.19.0            
# [16] R.methodsS3_1.7.0       affxparser_1.42.0      
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_0.12.2          base64enc_0.1-3      bitops_1.0-6         tools_3.2.3          digest_0.6.8         aroma.apd_0.6.0     
# [7] RSQLite_1.0.0        R.cache_0.12.0       DBI_0.3.1            parallel_3.2.3       R.rsp_0.21.0         S4Vectors_0.8.5     
# [13] PSCBS_0.60.0         globals_0.6.0        IRanges_2.4.6        stats4_3.2.3         Biobase_2.30.0       R6_2.1.1            
# [19] listenv_0.6.0        DNAcopy_1.44.0       AnnotationDbi_1.32.3 XML_3.98-1.3         magrittr_1.5         codetools_0.2-14    
# [25] matrixStats_0.50.1   BiocGenerics_0.16.1  assertthat_0.1       future_0.9.0         R.huge_0.9.0         stringi_1.0-1       
# [31] RCurl_1.95-4.7       chron_2.3-47    


#----------******iGEMS******----------------------#

# /*Load aroma.affymetrix and other libraries required in R*/

library(affxparser)
library(aroma.affymetrix)
library(affy)
library(FIRMAGene)
library(fdrtool)
library(limma)
library(biomaRt)
library(statmod)
library(stringr)
library(dplyr)
library(data.table)

####################----USER INPUT

# Please run the script form line 1 to 307  before making analysis setting/threshold choices.

# user must provide the name of folder in rawData with CEL files
name <- 'name of your study'

# the dataset to be used for annotation using biomart:for human replace "mmusculus_gene_ensembl" with "hsapiens_gene_ensembl"
dataset_mart<-"mmusculus_gene_ensembl" 

# add the reference genome/genome assembly for the cdf files used in analysis (info available on brainarray)
genome_assembly<-"GrCh38"

# add the version of brain array cdf used in the analysis (info available on brainarray)
cdf_version<-"BarrayV17"

# user needs to provide the path of working directory
path <- setwd('~/Splicing/iGEMS/')

# user specific setting dependent on hardware
setOption(aromaSettings, "memory/ram", 100)

# user needs to provide a design file specific to the study
# Design file (tab delimited) should comprise two columns, with cel file names in column 1 and designated group in column 2
# Read in the experiment design file and define your groups as "A" and "B" 

design<-read.table(paste(name,".txt",sep=""),sep="\t",header=FALSE)

# Load gene mapping coordinates: this file has seven columns,first column is the ID to which probes are mapped and correspond to gene ID 
ENSG_mapping <- read.delim("annotationData/chipTypes/MoEx-1_0-st-v1_ENSG_17/MoEx-1_0-st-v1_ENSG_17_mapping.txt")
colnames(ENSG_mapping)[1] <- 'gene'

#list of genes from ENSG mapping
genes_unique <- subset(ENSG_mapping, !duplicated(gene))

# Load exon mapping coordinates: this file has seven columns,first column is the ID to which probes are mapped and correspond to exon ID 
ENSE_Mapping <- read.delim("annotationData/chipTypes/MoEx-1_0-st-v1_ENSE_17/MoEx-1_0-st-v1_ENSE_17_mapping.txt")
colnames(ENSE_Mapping)[1] <- 'exon'

#Load annotation file containing list of exons and corresponding genes 
ENSE_desc <- read.delim("annotationData/chipTypes/MoEx-1_0-st-v1_ENSE_17/MoEx-1_0-st-v1_ENSE_17_desc.txt")

gene_exon_temp <- data.frame(ENSE_desc$Probe.Set.Name, str_extract(ENSE_desc$Description, '\\ENSMUSG...........'))

colnames(gene_exon_temp)[2]<- 'gene'
gene_exon_temp$at <- "at"
gene_exon_temp$gene <- paste(gene_exon_temp$gene,gene_exon_temp$at,sep='_')
gene_exon_at <- data.frame(gene_exon_temp$gene, gene_exon_temp$ENSE_desc.Probe.Set.Name)
colnames(gene_exon_at)<- c('gene','exon')
# exons per gene
ge_ex_DT <- data.table(gene_exon_at, key='gene')
grouped_by_gene <- group_by(ge_ex_DT, gene)
exons_per_gene <- summarise(grouped_by_gene, count = n_distinct(exon))

# upLimit may have to be defined depending on computer resources. e.g. macPro upLimit=500.

# gene_probenumber <- within(ENSG_mapping,{num <- as.numeric(ave(as.character(gene), as.character(gene), FUN=seq_along)) })
# set_upLimit <- max(gene_probenumber$num)+1

#OR
set_upLimit <- 500

#rawData folder should contain folder with the name provided above and the named subfolder should contain folder with the name of chipType

#Name of the cdf file  - gene level  - for example for HTA brain array gene level: hta20_Hs_ENSG_19
chipType_gene <- 'MoEx-1_0-st-v1_ENSG_17'
#Name of the cdf file  - exon level - for example for HTA brain array gene level: hta20_Hs_ENSE_19
chipType_exon <- 'MoEx-1_0-st-v1_ENSE_17'

### Start - CONVERT CDF
### convert the ASCII cdf to binary format:much faster and more efficient-need to RUN the following lines ONLY ONCE to create the binary cdfs ###
convertCdf("annotationData/chipTypes/MoEx-1_0-st-v1_ENSG_17/MoEx-1_0-st-v1_ENSG_17.cdf","annotationData/chipTypes/MoEx-1_0-st-v1_ENSG_17/MoEx-1_0-st-v1_ENSG_17,binary.cdf",verbose=TRUE)  
convertCdf("annotationData/chipTypes/MoEx-1_0-st-v1_ENSE_17/MoEx-1_0-st-v1_ENSE_17.cdf","annotationData/chipTypes/MoEx-1_0-st-v1_ENSE_17/MoEx-1_0-st-v1_ENSE_17,binary.cdf",verbose=TRUE)  

### END - CONVERT CDF ###

##################################################################################################################

#Load gene probes for determining number of probes per gene (FDR modelling step)
verbose <- Arguments$getVerbose(-30); timestampOn(verbose)

#####----- Change aroma setting according to cell definition file(ENSG_mapping)
setOption(aromaSettings, "medianPolishThreshold", c(set_upLimit,length(design$V1)))
setOption(aromaSettings, "models/RmaPlm/medianPolishThreshold", c(set_upLimit,length(design$V1)))
getOption(aromaSettings)

###---------------Gene level preprocessing-----------------###

# Read all the cel files that are there in you folder
cdf <- AffymetrixCdfFile$byChipType(chipType_gene, tags="binary")
cs <- AffymetrixCelSet$byName(name,cdf=cdf,verbose=verbose)

#Read a file with cel files of experiment design - it deals with >CEL than currently being used
mt<-match(design[,1],getNames(cs))
ds <- extract(cs, mt)

setCdf(ds,cdf)

#Background correction
bc <- RmaBackgroundCorrection(ds)
csBC <- process(bc,verbose=verbose)

#Quantile normalizationd
qn <- QuantileNormalization(csBC, typesToUpdate="pm")
csN <- process(qn, verbose=verbose)

#convert to unique cell
cdfu <- getUniqueCdf(cdf,verbose=-80)
csNU <- convertToUnique(csN,verbose=-20)

# probe level model
plmTr <- ExonRmaPlm(csNU, mergeGroups=TRUE)
fit(plmTr, verbose=verbose)

cesTr <- getChipEffectSet(plmTr)
trFit <- extractDataFrame(cesTr, addNames=TRUE)

# ,1 refers to column with probeset ID and ,6 onwards refers to expression values
trFit <- data.frame(trFit$unitName, trFit[,6:ncol(trFit)],check.names=FALSE)
colnames(trFit)[1] <- 'gene'

###-----------------Exon level preprocessing-------------------###

cdf_exon <- AffymetrixCdfFile$byChipType(chipType_exon, tags="binary")
cs_exon<-AffymetrixCelSet$byName(name,cdf=cdf_exon, verbose=verbose)

mt_exon<-match(design[,1],getNames(cs_exon))
ds_exon <- extract(cs_exon, mt_exon)

setCdf(ds_exon,cdf_exon)

#Background correction for exons
bc_exon <- RmaBackgroundCorrection(ds_exon)
csBC_exon <- process(bc_exon,verbose=verbose)

#Quantile normalization for exons
qn_exon <- QuantileNormalization(csBC_exon, typesToUpdate="pm")
csN_exon <- process(qn_exon, verbose=verbose)

#convert to unique cell for exons
cdfu_exon <- getUniqueCdf(cdf_exon,verbose=-80)
csNU_exon <- convertToUnique(csN_exon,verbose=-20)

#probe level model for exons
plmTr_exon <- ExonRmaPlm(csNU_exon, mergeGroups=TRUE)
fit(plmTr_exon, verbose=verbose)
cesTr_exon <- getChipEffectSet(plmTr_exon)
trFit_exon <- extractDataFrame(cesTr_exon, addNames=TRUE)
trFit_exon <- data.frame(trFit_exon$unitName,trFit_exon[,6:ncol(trFit_exon)],check.names=FALSE)
colnames(trFit_exon)[1] <- 'exon'

#Calculate Median of the replicates for each group/class
ind <- match(design[,1], getNames(ds))
cls <- design[,2][ind]

#mark which columns belong to which condition group
col_grA <- grep("A", cls)+1
col_grB <- grep("B", cls)+1

#calculate median expression values for each condition group per gene
trFit$grA_gene_med <- apply(trFit[,col_grA],1,median)
trFit$grB_gene_med <- apply(trFit[,col_grB],1,median)

#calculate mean expression and SD values  per exon
mean_exon <- apply(trFit_exon[,-1], 1, mean)
trFit_exon_sd <- data.frame(exon=trFit_exon$exon) 

#calculate median expression values for each condition group per exon
trFit_exon_sd$grA_exon_med <- apply(trFit_exon[,col_grA],1,median)
trFit_exon_sd$grB_exon_med <- apply(trFit_exon[,col_grB],1,median)

#calculate mean expression values for each condition group per exon as used for plotting
trFit_exon_sd$grA_exon_mean <- apply(trFit_exon[,col_grA],1,mean)
trFit_exon_sd$grB_exon_mean <- apply(trFit_exon[,col_grB],1,mean)

trFit_exon_sd$sd_exon <- apply(trFit_exon[,-1],1,sd)
trFit_exon_sd$sd_exon_grA <- apply(trFit_exon[,col_grA],1,sd)
trFit_exon_sd$sd_exon_grB <- apply(trFit_exon[,col_grB],1,sd) 

trFit_exon_sd$cov_exon_grA <- (trFit_exon_sd$sd_exon_grA/trFit_exon_sd$grA_exon_mean)*100
trFit_exon_sd$cov_exon_grB <- (trFit_exon_sd$sd_exon_grB/trFit_exon_sd$grB_exon_mean)*100

# group level sds
grA_ex_sd <- data.table(exon=as.character(trFit_exon_sd[,1]), grA_sd=trFit_exon_sd$sd_exon_grA)
grA_ex_sd$gene <- with(gene_exon_at, gene_exon_at[,1][match(gene_exon_at[,2], exon)]) # add genes
setkey(grA_ex_sd, key='gene')

grB_ex_sd <- data.table(exon=as.character(trFit_exon_sd[,1]), grB_sd=trFit_exon_sd$sd_exon_grB)
grB_ex_sd$gene <- with(gene_exon_at, gene_exon_at[,1][match(gene_exon_at[,2], exon)]) # add genes
setkey(grB_ex_sd, key='gene')

grA_grB_ex_SD<-merge(grA_ex_sd,grB_ex_sd,by="exon")
grA_grB_ex_SD<-subset(grA_grB_ex_SD,select=c(1:4))
colnames(grA_grB_ex_SD)[3]<-"gene"

# Plot the SD distribution to choose the SD threshold for your study (note: PDF version of plots is saved as well)
#options(scipen=10)
#par(mfrow=c(1,1))

# Distribution of SD values between range  0 to 30 in Group A
#zoom_in_A<-subset(grA_grB_ex_SD,grA_grB_ex_SD$grA_sd<=30)
#h_zoom_A <- hist(zoom_in_A$grA_sd, breaks=750, plot=FALSE)
#plot(h_zoom_A, main=('Distribution of exon level Standard Deviation in group A(SD<=30)'),xlab='Standard Deviation', col="black",axes=F,cex.main=0.75,cex.lab=0.75)
#axis(side=1,at=c(seq(0,30,2)),cex.axis=0.6)
#axis(side=2,at=c(seq(0,max(h_zoom_A$counts),500)),cex.axis=0.6)

# Distribution of SD values between range  0 to 30 in Group B
#zoom_in_B<-subset(grA_grB_ex_SD,grA_grB_ex_SD$grB_sd<=30)
#h_zoom_B <- hist(zoom_in_B$grB_sd, breaks=750, plot=FALSE)
#plot(h_zoom_B, main=('Distribution of exon level Standard Deviation in group B(SD<=30)'),xlab='Standard Deviation',border="gray",col="lightgray",axes =F,cex.main=0.75,cex.lab=0.75)
#axis(side=1,at=c(seq(0,30,2)),cex.axis=0.6)
#axis(side=2,at=c(seq(0,max(h_zoom_A$counts),500)),cex.axis=0.6)

# Plot PDF: SD distribution--this is saved in your working sirectory
pdf(paste(name,"Distribution of Standard deviation.pdf"),height=8,width=12)
options(scipen=10)
par(mfrow=c(1,1))

# Distribution of SD values between range  0 to 30 in Group A
zoom_in_A<-subset(grA_grB_ex_SD,grA_grB_ex_SD$grA_sd<=30)
h_zoom_A <- hist(zoom_in_A$grA_sd, breaks=750, plot=FALSE)
plot(h_zoom_A, main=('Distribution of exon level Standard Deviation in group A(SD<=30)'),xlab='Standard Deviation', col="black",axes=F,cex.main=1,cex.lab=0.75)
axis(side=1,at=c(seq(0,30,1)),cex.axis=1)
axis(side=2,at=c(seq(0,max(h_zoom_A$counts),500)),cex.axis=1)

# Distribution of SD values between range  0 to 30 in Group B
zoom_in_B<-subset(grA_grB_ex_SD,grA_grB_ex_SD$grB_sd<=30)
h_zoom_B <- hist(zoom_in_B$grB_sd, breaks=750, plot=FALSE)
plot(h_zoom_B, main=('Distribution of exon level Standard Deviation in group B(SD<=30)'),xlab='Standard Deviation',border="gray",col="gray",axes =F,cex.main=1,cex.lab=0.75)
axis(side=1,at=c(seq(0,30,1)),cex.axis=1)
axis(side=2,at=c(seq(0,max(h_zoom_A$counts),500)),cex.axis=1)

dev.off()
#---end of plot-----#

# setup threshold values for the different flters to be used in the downstream analysis- can be changed by the user

sd_threshold <- 0 #threshold for SD: step 1 filter to remove invariant genes, user can change this based on distribution plot of SD

nSamples <- 10000 #establishing FDR for MUF scores - number of permutations per probe_number

fdr_threshold_value <- 0.01 #fdr threshold value: step 2 filter for MUF scores

si_threshold_value <- 0.1 # si threshold: step 3 filter for exons with SI at extreme ends of distribution

limma_fdr_threshold <- 0.01 # limma FDR threshold: step 4 false positive filtering


###-----------------------------START SD FILTER--------------------###

# TRUE for genes where AT LEAST all-1 exons have sd < sd_threshold
# also drops single exon genes
grA_sd_filter <- grA_ex_sd[, length(which(grA_sd < sd_threshold)) >= length(exon)-1, by='gene'] 
grB_sd_filter <- grB_ex_sd[, length(which(grB_sd < sd_threshold)) >= length(exon)-1, by='gene'] 

# filter genes for removal
grA_out_genes <- grA_sd_filter[V1==TRUE] 
grB_out_genes <- grB_sd_filter[V1==TRUE] 

# RETAIN SINGLE EXON GENES FOR BG FILE
# Extract single exon genes
exons_per_gene_single <- subset(exons_per_gene, exons_per_gene$count==1)

grA_ex_sd_single<-subset(grA_ex_sd,gene %in% exons_per_gene_single$gene)
grB_ex_sd_single<-subset(grB_ex_sd,gene %in% exons_per_gene_single$gene)

# Filter genes for removal
grA_sd_filter_out_genes_single <- grA_ex_sd_single[grA_sd< sd_threshold] 
grB_sd_filter_out_genes_single <- grB_ex_sd_single[grB_sd< sd_threshold] 

# Extract expressed genes single
grA_expressed_genes_single <- grA_ex_sd_single[grA_sd > sd_threshold]
grB_expressed_genes_single <- grB_ex_sd_single[grB_sd > sd_threshold]

grA_all_out_genes <- data.frame(gene = setdiff(grA_out_genes$gene,grA_expressed_genes_single$gene))
grB_all_out_genes <- data.frame(gene = setdiff(grB_out_genes$gene,grB_expressed_genes_single$gene))
all_out_genes <- data.frame(genes_removed = union(as.character(grA_all_out_genes[,"gene"]), as.character(grB_all_out_genes[,"gene"])))
names(all_out_genes) <- 'gene'

grA_bg_genes <- data.frame(gene= setdiff(grA_ex_sd$gene,grA_all_out_genes$gene))
grB_bg_genes <- data.frame(gene= setdiff(grB_ex_sd$gene,grB_all_out_genes$gene))

# add corresponding SD values 
all_out_genes_sd<-merge(all_out_genes,grA_grB_ex_SD,by="gene")

grA_all_out_genes_sd<-subset((merge(grA_all_out_genes,grA_grB_ex_SD,by="gene")),select=c(1,2,3))
grB_all_out_genes_sd<-subset((merge(grB_all_out_genes,grA_grB_ex_SD,by="gene")),select=c(1,2,4))

grA_bg_genes_sd <- subset((merge(grA_bg_genes,grA_grB_ex_SD,by="gene")),select=c(1,2,3))
grB_bg_genes_sd <- subset((merge(grB_bg_genes,grA_grB_ex_SD,by="gene")),select=c(1,2,4))

## END sd FILTER ##

###Calculate Gene level Residuals

rs <- calculateResidualSet(plmTr, verbose=verbose)
rsu <- readUnits(rs, verbose = verbose)

#Collapse residuals from samples to groups
new_rsu <- list()

for( i in 1:length(rsu))
{
  new_rsu[[i]] <- rsu[[i]][[1]][[1]]
  if (length(rsu[[i]]) > 1)
  {
    for (j in 2:length(rsu[[i]]))
    {
      new_rsu[[i]] <- rbind(new_rsu[[i]],rsu[[i]][[j]][[1]])
    }
  }
}
names(new_rsu) <- names(rsu)

# remove genes that fail the sd filter BEFORE we calculate MUF scores
# new_rsu is list, remove based on name match
outs <- which(names(new_rsu)%in%all_out_genes[,1])

# if sd_threshold >0 then run the following line and if  sd_threhold=0 then uncomment and run second line
new_rsu <- new_rsu[-outs]
#new_rsu <- new_rsu

# low sd genes removed from new_rsu

#Calculate Median of the replicates for each group/class
ind <- match(design[,1], getNames(ds))
cls <- design[,2][ind]

#Define a function to calculate median
calcMeans_median <- function(x, cls)
{
  ucls <- unique(cls)
  v <- matrix(, nrow = nrow(x), ncol = length(ucls))
  for (i in 1:length(ucls)) {
    ind <- which(cls == ucls[i])
    if (length(ind) > 1) {
      xx <- as.matrix(x[, ind])
      v[, i] <- apply(xx,1,median)
    }
    else {
      v[, i] <- x[, ind]
    }
  }
  colnames(v) <- ucls
  v
}
# now applying to gene filtered new_rsu
rsu1 <- lapply(new_rsu, FUN = function(u, cls) {
  r <-log2(u)
  calcMeans_median(r, cls)
}, cls = cls)

#Calculate number of probes per Gene
nProbe <- sapply(new_rsu, FUN = function(u) nrow(u))
minProbes <- 3
w <- (nProbe >= minProbes) & (nProbe <= set_upLimit)
nProbe <- nProbe[w]

gc()
cat("Calculating MUF score (observed data).\n")
#Calculate muf scores
rsu1_diff <- lapply(rsu1,FUN=function(u) matrix(data = u[,2]-u[,1],ncol=1))
mufScores <- lapply(rsu1_diff[w], FUN = function(u) c(mufColumns(u)))
mufScores <- as.data.frame(unlist(mufScores))
colnames(mufScores) <- c("mufScores")
mufScores$abs_mufscore <- abs(mufScores$mufScores)

# Extracting residual matrix for permutation (empirical distribution)
cat("Extracting residual matrix.\n")
resMatrix <- matrix(unlist(rsu1, use.names = FALSE), byrow = FALSE, ncol = length(unique(cls)))
resMatrix <- resMatrix[rowSums(is.na(resMatrix) | is.infinite(resMatrix)) ==0, ]
gc()
cat("Calculating MUF score (residual permutation).\n")
v <- unique(sort(nProbe))

#Generate empirical distributions
set.seed(12345)
z <- lapply(v, FUN = function(u) {
  if (verbose)
    cat(u, " ", sep = "")
  abs((mufColumns(matrix(sample(resMatrix[,2], size = u * nSamples, replace=TRUE),nc = nSamples)-matrix(sample(resMatrix[,1], size = u * nSamples, replace=TRUE),nc = nSamples))))
})

mufScores$nProbes <- nProbe
distributions <- lapply(z, FUN=function(anon)
{
  ecdf(anon)
}
)
names(distributions) <- v

# calculate FDR
fdr <- matrix(0,nrow=nrow(mufScores),ncol=2)
colnames(fdr)<- c('mufscore','fdr')
rownames(fdr) <- rownames(mufScores)

for (gene in 1:nrow(mufScores)){
  print(gene)
  fdr[gene,'fdr'] <- 1-distributions[[which(names(distributions)==mufScores$nProbes[gene])]](abs(mufScores$mufScores[gene]))
}

fdr_df <- as.data.frame(fdr)
#add mufScores and nProbes
fdr_df$mufscore <- mufScores$mufScores
fdr_df$nProbes <- mufScores$nProbes
fdr_df$abs_mufscore <- mufScores$abs_mufscore

#######---Exon level analysis - calculations for  third filter(SI)########

# function to calculate fold change at exon level
fc <- function(x,y){
  if (x>y) {
    fc <- x/y
  }
  else
    fc <- -1/(x/y)
}
trFit_exon$fc <- mapply(fc,trFit_exon_sd[,'grB_exon_med'],trFit_exon_sd[,'grA_exon_med'])

exon_expr <- trFit_exon[,2:(length(col_grA)+length(col_grB)+1)]

trFit_exon_all <- data.frame(trFit_exon_sd$exon, trFit_exon_sd$grA_exon_med, trFit_exon_sd$grB_exon_med,trFit_exon_sd$grA_exon_mean, trFit_exon_sd$grB_exon_mean,trFit_exon$fc,trFit_exon_sd$sd_exon_grA,trFit_exon_sd$sd_exon_grB)
colnames(trFit_exon_all) <- c('exon','grA_exon_med', 'grB_exon_med','grA_exon_mean', 'grB_exon_mean','fc','sd_exon_grA','sd_exon_grB')

# merge gene id with exon id
ge <- merge(gene_exon_at, trFit, by='gene')

# extract only columns of interest
ge_med <- ge[,c('gene','exon','grA_gene_med','grB_gene_med')]

# merge with exon data
ge_ex <- merge(ge_med, trFit_exon_all, by='exon')

# calculate Splicing Index 
ge_ex$NI_grA <- ge_ex[,'grA_exon_med']/ge_ex[,'grA_gene_med']
ge_ex$NI_grB <- ge_ex[,'grB_exon_med']/ge_ex[,'grB_gene_med']
ge_ex$si <- log2(ge_ex[,'NI_grB']/ge_ex[,'NI_grA'])
ge_ex$abs_si <- abs(ge_ex$si)

# add consecutive exon number
exon_map <- as.data.frame(ENSE_Mapping)
colnames(exon_map)[1]<-"exon"
colnames(exon_map)<-gsub(" ","",colnames(exon_map))# remove empty spaces from column names
exon_map_neg<- subset(exon_map, exon_map$Chr.Strand == '-')
exon_map_pos<- subset(exon_map, exon_map$Chr.Strand == '+')
exon_map_neg<- exon_map_neg[,1:4][with(exon_map_neg[,1:4], order(exon_map_neg[,1:4]$exon, -exon_map_neg$Chr.From)),]
exon_map_pos <- exon_map_pos[,1:4][with(exon_map_pos[,1:4], order(exon_map_pos[,1:4]$exon, exon_map_pos$Chr.From)),]

exon_map <- rbind(exon_map_neg, exon_map_pos)
exon_map_u <- subset(exon_map, !duplicated(exon))
exon_df <- merge(ge_ex, exon_map_u, by='exon')

exon_df_neg<- subset(exon_df, exon_df$Chr.Strand == '-')
exon_df_pos<- subset(exon_df, exon_df$Chr.Strand == '+')
exon_df_neg<- exon_df_neg[with(exon_df_neg,order(exon_df_neg$gene,-exon_df_neg$Chr.From)),]
exon_df_seq_neg <- within(exon_df_neg,{num <- as.numeric(ave(as.character(gene), as.character(gene), FUN=seq_along)) })

exon_df_pos <- exon_df_pos[with(exon_df_pos,order(exon_df_pos$gene,exon_df_pos$Chr.From)),]
exon_df_seq_pos <- within(exon_df_pos,{num <- as.numeric(ave(as.character(gene), as.character(gene), FUN=seq_along)) })
exon_df_seq <- rbind(exon_df_seq_neg, exon_df_seq_pos)

# add data from previous step. i.e. nProbes, fdr and mufScores
fdr_df$gene <- row.names(fdr_df) 

##Extract information about, associated to ensemble gene id, gene names
#connect to biomart
mart<- useMart('ENSEMBL_MART_ENSEMBL',dataset=dataset_mart, host="www.ensembl.org")
#create vector with ensembl geneIds from object gene_exon
ensembl_genes <- gsub("_at","",fdr_df$gene)
#Extract associated gene names from biomart using lastest version of ensembl annotation (GRCH38.v3)
gene_names <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "external_gene_name"),values= ensembl_genes,mart= mart)
#Rename columns in object gene_names in order to merge them with gene_exon object
colnames(gene_names) <- c('gene', 'gene_name')

#Create a dataframe with missing genes in database and assign NA as gene names
mis_genes <- data.frame(gene=setdiff(unique(ensembl_genes), gene_names$gene),gene_name=NA)
#Bind ensembl gene ids with associated gene name with ensembl gene ids without associated gene name
gene_names_all <- rbind(gene_names, mis_genes)

## add back the _at suffix to the gene Ids to make it consistent with other data frames
gene_names_all$gene<-paste(gene_names_all$gene,"_at",sep="")
fdr_df_final<-merge(fdr_df,gene_names_all,by="gene")

final_all <- subset((merge(exon_df_seq, fdr_df_final, by='gene')),select=c(1:2,24,3:23))

final_all$gene_med_ratio <- final_all$grA_gene_med/final_all$grB_gene_med
# list of all genes that were taken for consideration in the analysis for GO and for mufScore vs fdr plot (only gene level information of interest).
final_uniqueGene <- subset(final_all, !duplicated(gene))

#-----Selection of significant MUF scores based on FDR----#

final_fdr001 <- subset(final_all, final_all$fdr < fdr_threshold_value)
final_fdr001_u <- subset(final_fdr001, !duplicated(gene)) # unique fdr001 genes

#------Selection of exons with SI  based on thresholds------#

si_threshold_sel_pos <- quantile(final_fdr001$si[final_fdr001$si>0], 1-si_threshold_value)
si_threshold_sel_neg <- quantile(final_fdr001$si[final_fdr001$si<0], si_threshold_value)


#top % SI 
final_fdr001_si <- subset(final_fdr001, final_fdr001$si>si_threshold_sel_pos | final_fdr001$si < si_threshold_sel_neg)
final_fdr001_si_u <- subset(final_fdr001_si, !duplicated(gene))


####---Final filter:differential expression on exon using limma filter----#####

data.t_fil <- trFit_exon[,1:(length(cls)+1)]

all.levels <-c(rep("0",length(cls)))
grA<-grep("A",cls)
grB<-grep("B",cls)

all.levels[grA] <- c("grA")
all.levels[grB] <- c("grB")

all.levels = factor(all.levels, levels =c("grA","grB"))
all.levels
design <- model.matrix(~0 + all.levels)
colnames(design) <- levels(all.levels)
contrast.matrix <- makeContrasts(grB-grA, levels=design)
contrast.matrix
data_fit <- lmFit(log2(data.t_fil[2:(length(grA)+length(grB)+1)]), design)
data_fit2 <- contrasts.fit(data_fit, contrast.matrix)
data_fit2 <- eBayes(data_fit2)
sel.diff<-as.matrix(p.adjust(data_fit2 $F.p.value,method="fdr"))
sel.diff_per<- as.matrix(apply(sel.diff,1,FUN=function(x) x*100))

result_limma <- data.frame(cbind(as.character(data.t_fil[,1]),as.numeric(sel.diff[,1]),sel.diff_per[,1]))
colnames(result_limma) <- c("exon","FDR_LIMMA","FDR_LIMMA(in%)")

# get the exon  level FDR values for exons in our splicing list 
final_fdr001_si_fdr_ex <- merge(final_fdr001_si, result_limma, by='exon')
final_fdr001_si_fdr_ex$is_exon_expr_lessthen_gene_expr <- (final_fdr001_si_fdr_ex$grA_exon_med<final_fdr001_si_fdr_ex$grA_gene_med) & (final_fdr001_si_fdr_ex$grA_exon_med<final_fdr001_si_fdr_ex$grB_gene_med) & (final_fdr001_si_fdr_ex$grB_exon_med<final_fdr001_si_fdr_ex$grA_gene_med) & (final_fdr001_si_fdr_ex$grB_exon_med<final_fdr001_si_fdr_ex$grB_gene_med) 


#------Selection of final list of events of AEU------#

# top 10% lower expre A and B than gene A expr and gene B expr
final_fdr001_si_fp <- subset(final_fdr001_si_fdr_ex,final_fdr001_si_fdr_ex$grA_exon_med<final_fdr001_si_fdr_ex$grA_gene_med & final_fdr001_si_fdr_ex$grA_exon_med<final_fdr001_si_fdr_ex$grB_gene_med & final_fdr001_si_fdr_ex$grB_exon_med<final_fdr001_si_fdr_ex$grA_gene_med & final_fdr001_si_fdr_ex$grB_exon_med<final_fdr001_si_fdr_ex$grB_gene_med &as.numeric(as.character(final_fdr001_si_fdr_ex$FDR_LIMMA))>=limma_fdr_threshold)

mt <- match(final_fdr001_si_fp$exon,final_fdr001_si_fdr_ex$exon)
final_fdr001_si_tp <- final_fdr001_si_fdr_ex[-mt,]
final_fdr001_si_tp_u <- subset(final_fdr001_si_tp, !duplicated(gene))
rownames(final_fdr001_si_tp_u) <- NULL

summary_final_fdr001_si_tp<-subset(final_fdr001_si_tp, select=c(1,2,3,20,23,4,5,6,7,11,12,24,22,15,27,28))

# dataframe of metrics and statistics for all the exons in the analysis after SD filter
final_all_metrics<-merge(final_all,result_limma,by="exon")

#dataframe including all genes for expression plots and visualisation:
final_all_plot <- final_all[,c('gene','exon','grA_exon_mean','grB_exon_mean','num','sd_exon_grA','sd_exon_grB')]
final_all_plot <- final_all_plot[with(final_all_plot,order(final_all_plot$gene,final_all_plot$num)),]


########################################################################
#------------------------------PLOTS-----------------------------------#
########################################################################

# Plot 1: mufScore vs fdr
par(mfrow=c(1,1))
pdf(paste(name,"Distribution of mufScore vs FDR.pdf"),height=8,width=12)
final_uniqueGene_plot <- data.frame(final_uniqueGene$mufscore, final_uniqueGene$fdr)
colnames(final_uniqueGene_plot) <- c('mufscore','fdr')
final_uniqueGene_plot$fdr_perc <- final_uniqueGene_plot$fdr*100
plot(abs(final_uniqueGene_plot$mufscore) ~ (final_uniqueGene_plot$fdr_perc), main='MUF scores vs FDR', xlab='FDR (%)', ylab='MUF scores')
abline(v=1, col=2)
abline(v=5, col=5)
text(1,10,'1%',cex=0.9,col=2)
text(5,12,'5%',cex=0.9,col=5)
dev.off()
#---end of plot 1

# Plot 2: SI distribution with coloured thresholds
pdf(paste(name,"Distribution of splicing index.pdf"),height=8,width=12)
options(scipen=10)
h <- hist(final_all$si, breaks=100, plot=FALSE)
cuts <- cut(h$breaks, c(-Inf,si_threshold_sel_neg,si_threshold_sel_pos, Inf))
plot(h, main=('Distribution of Splicing Index'),xlab='Splicing Index', col=c("red","grey","red")[cuts],axes=F)
axis(side=1,at=c(seq(min(h$breaks),max(h$breaks),1)))
axis(side=2,at=c(seq(min(h$counts),max(h$counts),5000)))
dev.off()
#---end of plot 2

# Plot 3:SI boxplot
pdf(paste(name,"Boxplot for splicing index.pdf"),height=8,width=12)
boxplot(final_fdr001_si_fdr_ex$si, main='Splicing Index distribution')
dev.off()
#---end of plot 3

###########----Table outputs------#########

# list of genes available on the particular platform and cdf
write.table(genes_unique[,1], file=paste(name,'all genes_list',genome_assembly,cdf_version,'.txt',sep='_'), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

write.table(trFit, file=paste(name,'genelevel_expression',genome_assembly,cdf_version,'.txt',sep='_'), row.names=FALSE, quote=FALSE, sep='\t')
write.table(trFit_exon, file=paste(name,'exonlevel_expression',genome_assembly,cdf_version,'.txt',sep='_'), row.names=FALSE, quote=FALSE, sep='\t')

write.table(all_out_genes_sd, paste('sd_limit', sd_threshold,'genes_removed.txt', sep='_'), quote=FALSE, sep='\t', row.names=FALSE)
write.table(grA_all_out_genes_sd, file=paste('sd_limit',sd_threshold,'grA_genes_removed.txt', sep='_'), quote=FALSE, sep='\t', row.names=FALSE)
write.table(grB_all_out_genes_sd, file=paste('sd_limit', sd_threshold, 'grB_genes_removed.txt', sep='_'), quote=FALSE, sep='\t', row.names=FALSE)

write.table(grA_bg_genes_sd, paste('sd_limit', sd_threshold,'grA_expressed_genes.txt', sep='_'), quote=FALSE, sep='\t', row.names=FALSE)
write.table(grB_bg_genes_sd, paste('sd_limit', sd_threshold,'grB_expressed_genes.txt', sep='_'), quote=FALSE, sep='\t', row.names=FALSE)

write.table(final_all_metrics, file=paste(name,"all_data_statistics_metrics.txt",sep='_'), row.names=FALSE, quote=FALSE, sep='\t')

write.table(final_fdr001, file=paste(name,sd_threshold,"SD_filtered_FDR_filtered",fdr_threshold_value,'MUF_scores.txt',sep='_'), row.names=FALSE, quote=FALSE, sep='\t')

write.table(result_limma,file=paste(name,"FDR_limma.csv",sep='_'),sep=",",quote=FALSE)

write.table(final_fdr001_si, file=paste(name,sd_threshold,'SDfiltered_FDR_filtered',fdr_threshold_value,'MUF_scores',si_threshold_value,'SI_filtered.txt',sep='_'), row.names=FALSE, quote=FALSE, sep='\t')

write.table(final_fdr001_si_fp, file=paste(name,sd_threshold,'SDfiltered_FDR_filtered',fdr_threshold_value,'MUF_scores',si_threshold_value,'SI_filtered_falsepositives.txt',sep='_'), row.names=FALSE, quote=FALSE, sep='\t')

write.table(final_fdr001_si_tp, file=paste(name,sd_threshold,'SDfiltered_FDR_filtered',fdr_threshold_value,'MUF_scores',si_threshold_value,'SI_filtered_truepositives.txt',sep='_'), row.names=FALSE, quote=FALSE, sep='\t')

write.table(summary_final_fdr001_si_tp, file=paste('Summary',name,sd_threshold,'SDfiltered_FDR_filtered',fdr_threshold_value,'MUF_scores',si_threshold_value,'SI_filtered_truepositives.txt',sep='_'), row.names=FALSE, quote=FALSE, sep='\t')

write.table(final_all_plot, file=paste(name,'plotting_data_mean_error_bars',genome_assembly,cdf_version,'CDF.txt',sep='_'), row.names=FALSE, quote=FALSE, sep='\t')

####---To plot Intensity plots use following-----#####
# without log transformation
par(mar=c(8.0, 4.0, 1.0, 1.0))
gene <- 'ENSG00000005471_at' # pick a gene

ylim <- c(0, max(c(max(final_all_plot[final_all_plot$gene==gene,'grA_exon_mean'])+max(final_all_plot[final_all_plot$gene==gene,'sd_exon_grA']),max(final_all_plot[final_all_plot$gene==gene,'grB_exon_mean'])+max(final_all_plot[final_all_plot$gene==gene,'sd_exon_grB']))))

pdf(paste(gene,"intensity plot.pdf"),height=10,width=12)
plot(final_all_plot[final_all_plot$gene==gene,'grA_exon_mean'], ylab='Linear intensities', xlab='', ylim=ylim ,pch=20, xaxt='n', main=(gsub('_at', '', gene)))

x0 <- final_all_plot[final_all_plot$gene==gene,'num']
y0 <- final_all_plot[final_all_plot$gene==gene,'grA_exon_mean']
x1 <- final_all_plot[final_all_plot$gene==gene,'num']
y1 <- final_all_plot[final_all_plot$gene==gene,'grA_exon_mean']-final_all_plot[final_all_plot$gene==gene,'sd_exon_grA']

x2 <- final_all_plot[final_all_plot$gene==gene,'num']
y2 <- final_all_plot[final_all_plot$gene==gene,'grA_exon_mean']+final_all_plot[final_all_plot$gene==gene,'sd_exon_grA']

arrows(x0,y0,x1,y1, code=2, angle=90, length=0.05, col='black')
arrows(x0,y0,x2,y2,code=2, angle=90, length=0.05, col='black')
points(final_all_plot[final_all_plot$gene==gene,'grA_exon_mean'], ylab='Linear intensities', xlab='', ylim=ylim ,pch=20, xaxt='n')

axis(1,at=1:nrow(final_all_plot[final_all_plot$gene==gene,]),labels=substring(final_all_plot$exon[final_all_plot$gene==gene],1,18),cex.axis=0.5,las=2)
xb0 <- final_all_plot[final_all_plot$gene==gene,'num']
yb0 <- final_all_plot[final_all_plot$gene==gene,'grB_exon_mean']
xb1 <- final_all_plot[final_all_plot$gene==gene,'num']
yb1 <- final_all_plot[final_all_plot$gene==gene,'grB_exon_mean']-final_all_plot[final_all_plot$gene==gene,'sd_exon_grB']
xb2 <- final_all_plot[final_all_plot$gene==gene,'num']
yb2 <- final_all_plot[final_all_plot$gene==gene,'grB_exon_mean']+final_all_plot[final_all_plot$gene==gene,'sd_exon_grB']
arrows(xb0,yb0,xb1,yb1, code=2, angle=90, length=0.05, col=2)
arrows(xb0,yb0,xb2,yb2,code=2, angle=90, length=0.05, col=2)
points(final_all_plot[final_all_plot$gene==gene,'grB_exon_mean',],col=2, pch=18)

legend(2,(max(ylim)-20),c("Group A","Group B"),col=c("black","red"),pch=16)
dev.off()
