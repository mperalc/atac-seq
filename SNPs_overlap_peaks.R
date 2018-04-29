# overlap ATAC peaks and T2D GWAS SNPs

### load libraries
library(dplyr)
library(biomaRt)
#library(readr)
#library(vcfR)
library(GenomicRanges)
library(BSgenome) # for SNP by overlaps function

#### input files
counts <- read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/counts_all_samples_peaks_narrowpeaks_conservative_set_highQualityBams_renamed.txt",
                     header = T) # bed file with start and end of peaks
peaks = counts[1:3]
rm(counts)

credset = read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_T2D/DIAGRAM_T2D_1000G_99_Credsets/credible_sets_diagram_forMAGENTA.txt",
                     header = F) # 99 credible sets T2D 2017 DIAGRAM
credset=unique(credset[1:2]) # SNPs with unique positions (some might be shared between loci?)
colnames(credset)=c("Chr","Pos")
nrow(peaks)
any(duplicated(peaks)) # there are no duplicated peaks
#[1] FALSE
nrow(credset)

# save the step of assigning SNP id for after intersecting peaks
# start_time <- Sys.time()
# snps = read.vcfR( "/Users/Marta/Documents/WTCHG/DPhil/Data/SNPs_b150_GRCh37p13/common_all_20170710.vcf.gz", verbose = FALSE ) # common SNPs build 150 GRCh37p13
# end_time - start_time

# remove sex chromosomes
r1=rownames(peaks[grep("chrY", rownames(peaks), value=F),])
r2=rownames(peaks[grep("chrX", rownames(peaks), value=F),])

peaks=peaks[!rownames(peaks) %in% c(r1,r2),]
rm(r1,r2)
nrow(peaks)

########### 
# 
# # check SNP names
# # load for SNPs 
# ensembl_snp = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_snp")
# 
# 
# getSNPsInDf = function(credset){
#   # biomart query forward strand
#   filterlist <- list(paste(credset[,"Chr"],credset[,"Pos"],credset[,"Pos"],1,sep = ":"))
#   
#   results1 <- getBM(attributes = c("refsnp_id",'chr_name', "chrom_start"),
#                    filters = c("chr_name"),values = list("1"), mart = ensembl_snp)
#   # biomart query reverse strand
#    filterlist <- list(paste(credset[n,"Chr"],credset[n,"Pos"],credset[n,"Pos"],-1,sep = ":"))
#   
#   results2 <- getBM(attributes = c("refsnp_id",'chr_name', "chrom_start", "chrom_end"),
#                     filters = c("chromosomal_region"),values = filterlist, mart = ensembl_snp)
#   res=rbind(results1,results2)
#   
#   return(res)
# }
# 
# counts$Chr=gsub("\\chr*","",counts$Chr) #remove "chr" in every element in location
# 


colnames(credset)=c("chr","start")
credset$end=credset$start

credset.gr = credset %>% 
  mutate(chr=paste0('chr', chr)) %>%   
  makeGRangesFromDataFrame

peaks.gr = makeGRangesFromDataFrame(peaks) 

# Get all credible set SNPs overlapping ATAC peaks:
overlapSNP = subsetByOverlaps(credset.gr,peaks.gr)


# Get all peaks overlapping credible set SNPs :
overlapPeak = subsetByOverlaps(peaks.gr,credset.gr)
#overlapPeak=as.data.frame(overlapPeak)

# more SNPs overlaping than peaks, so there must be several SNPs within each peak in some cases

# % of SNPs from credset in peaks
(length(overlapSNP)/nrow(credset))*100

# % of peaks with intersecting SNPs from credset
(length(overlapPeak)/nrow(peaks))*100

#### Annotate by promoter, coding, etc
# 
# 

source('~/WTCHG/R scripts/atac-seq/annotate_GRanges_TxDb.R') # loading whoile genome annotation function

# ATAC peaks and GWAS credible set SNPs
allPeaksAnnotated = annotate_GRanges_TxDb(peaks.gr) 
allPeaksAnnotated$anno = rownames(allPeaksAnnotated) # shaping for ggplot
allPeaksAnnotated$type= "All ATAC-seq peaks"
peaksOverGWASAnnotated = annotate_GRanges_TxDb(overlapPeak) 
peaksOverGWASAnnotated$anno = rownames(peaksOverGWASAnnotated)
peaksOverGWASAnnotated$type= "ATAC-seq peaks overlapping credible sets"
GWASAnnotated = annotate_GRanges_TxDb(credset.gr)
GWASAnnotated$anno = rownames(GWASAnnotated)
GWASAnnotated$type= "All SNPs in credible sets"

source('~/WTCHG/R scripts/annotate_whole_genome_TxDb.R') # loading whole genome annotation function

whole = annotate_whole_genome_TxDb() # whole genome annotation
whole$anno = rownames(whole)
whole$type = "Whole genome"
  
combined = rbind(allPeaksAnnotated,peaksOverGWASAnnotated,GWASAnnotated, whole)

combined$anno = factor(combined$anno,levels = unique(combined$anno),ordered = T,
                       labels = c("Promoter","Exon","Intron","Intergenic")) # ordering annotations
combined$type = factor(combined$type,levels = unique(combined$type),ordered = T) # ordering labels of data type

# colour scale for plot

pal=c("#96165a","#bf3b5c","#3D8DC1","#95C8EB")   # red and blue scale

source('~/WTCHG/R scripts/plotting/ggplot_stacked_barplot_percent.R') # loading barplot function

p = ggplot_stacked_barplot_percent(combined,pal,0.75 )

png("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/TxDb_annotation_peaks_credible_set_T2D_DIAGRAM.png", type="cairo",
    width=14,height=8,units="in",res=500,pointsize = 12)
p 
dev.off()


## save peaks that overlap GWAS, and get nearest gene

library(ChIPpeakAnno)
data(TSS.human.GRCh37) # load TSS info from CRCh37

# range data and annotate position with respect to closest TSS
# This way I can see if it's in a promoter, inside a gene or in an intergenic region

range_data = function(data){
  rangedpeak = RangedData(IRanges(start=data$Start,end=data$End), # passing peak tables as ranged data
                          names=rownames(data),space=data$Chr)
  return(rangedpeak)
}

rangedPeaks = range_data(peaks)

# Two methods of getting the nearest gene, both could be valid in principle

annotatedPeakNearestLoc = annotatePeakInBatch(overlapPeak, AnnotationData=TSS.human.GRCh37) # annotate with hg19. 
annotatedPeakShortestDistance = annotatePeakInBatch(overlapPeak, AnnotationData=TSS.human.GRCh37,output = "shortestDistance")
# "shortestDistance" will output nearest features to peaks
# Displays nearest features to TSS

annotatedPeakNearestLoc = as.data.frame(annotatedPeakNearestLoc )
annotatedPeakShortestDistance = as.data.frame(annotatedPeakShortestDistance)

write.table(annotatedPeakNearestLoc,file="/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/Peaks_over_SNPs_AllPeakToGeneAnnotated.txt",
            sep="\t", quote = F, row.names = F, col.names = T)

write.table(annotatedPeakShortestDistance,file="/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/Peaks_over_SNPs_annotatedPeakShortestDistance.txt",
            sep="\t", quote = F, row.names = F, col.names = T)

# annotate with gene names



write.table(as.data.frame(overlapPeak),file="/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/Peaks_over_SNPs.txt",
            sep="\t", quote = F, row.names = F, col.names = T)

## save SNPs that overlap peaks, and get nearest gene

overlapSNPdf = as.data.frame(overlapSNP)
colnames(overlapSNPdf) = c("chr","start","end","width","strand")
overlapSNPdf$chr=gsub("\\chr*","",overlapSNPdf$chr) #remove "chr" in every element in location

# load ensembl
ensembl_SNP = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_snp")
# # get regions to query in correct format

chr.regionF=paste(overlapSNPdf$chr,overlapSNPdf$start,overlapSNPdf$end, "1", sep = ":") # forward strand
#chr.regionR=paste(overlapSNPdf$chr,overlapSNPdf$start,overlapSNPdf$end, "-1", sep = ":") # reverse strand


# biomart in R is a pain of timeouts
# saving regions to query with browser
write.table(data.frame(chr.regionF),file="/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/overlapSNPDF.txt",
            +             sep="\t", quote = F, row.names = F, col.names = F)
# start_time <- Sys.time()
# 
# # batch queries of regions 
# # I had to reduce the regions queried at one time because some take too long and cause a connection timeout
# # at the same time, calling getBM multiple times seems to cause longer wait 
# # I'm pushing their servers a bit with so many calls because I don't see another way around it
# 
# resultsF = list()
# 
# get_SNP_ids_at_SNP_df_intervals <- function(overlapSNPdf,ensembl_SNP, i, interval) {
#   # getting from i to f at interval from 1 to end of dataframe
#   for(f in seq(from=i + interval, to=nrow(overlapSNPdf), by=interval)){
#     if( nrow(overlapSNPdf) - f > interval){
#       resultsF[[i]] <- getBM(attributes = c('refsnp_id','chr_name','chrom_start',"associated_gene","associated_variant_risk_allele"),
#                              filters = c('chromosomal_region'),
#                              values = list(as.factor(chr.regionF[i:f])), mart = ensembl_SNP) # chromosomal region needs to be a list of factors for this to work
#       
#       
#     }
#     else{ # last batch goes through this "else"
#       resultsF[[i]] <- getBM(attributes = c('refsnp_id','chr_name','chrom_start','chrom_strand'),
#                              filters = c('chromosomal_region'),
#                              values = list(as.factor(chr.regionF[i:nrow(overlapSNPdf)])), mart = ensembl_SNP) 
#      
#     }
#     i = f
#     assign("i", i, envir = .GlobalEnv) # to pass the value of i to outside the functions
#     message(paste("f is ", f,sep = ""))
#   }
#   return(resultsF)
# }
# 
# attempt <- 1
# i = 1
# while( length(resultsF)==0 && attempt <= 200 ) {
#   try(
#     if(attempt==1 | i==1){
#     resultsF<- get_SNP_ids_at_SNP_df_intervals(overlapSNPdf,ensembl_SNP, i = 1, interval = 100)
#     }
#     else{ # here i has been received from within function in previous attempts
#       resultsF <- get_SNP_ids_at_SNP_df_intervals(overlapSNPdf,ensembl_SNP, i = i, interval = 100)
#       
#     }
#   )
#   attempt <- attempt + 1
#   
# } # if you want to stop this function while running you need to terminate R
# 
# end_time <- Sys.time()
# end_time - start_time
# 
# results = do.call("rbind",results)
# 
# write.table(results,file="/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/SNPs_in_ATAC_peaks.txt",
#             sep="\t", quote = F, row.names = F, col.names = T)
