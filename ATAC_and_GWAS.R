# Investigating correlations and enrichment of ATAC-seq peaks and GWAS hits

# How many ATAC-seq peaks are present inside GWAS credible regions? 

# T2D DIAGRAM - 150,000 indiv imputed 1000Genomes Scott et al 2017

T2D_regions= read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_T2D/filenames_credible_regions_Feb17.txt",header = F)

# modifying region table
T2D_regions$V1=gsub("\\chr*","",T2D_regions$V1) #remove "chr" in every element in location
region_names=sub(".*_","",T2D_regions$V1) # remove everything before last "_"
region_names=sub("\\..*","",region_names) # remove everything after dot
T2D_regions$V1=gsub("\\_[^_]*$","",T2D_regions$V1) # remove everything after last "_"
T2D_regions$V1=gsub("_", ":", T2D_regions$V1)  # replace _ for :
T2D_regions$V1=gsub("-", ":", T2D_regions$V1)  # replace - for :
T2D_regions=cbind(T2D_regions,region_names)
colnames(T2D_regions)[1]=c("location")
rm(region_names)

# read in ATAC-seq peaks file
counts <- read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/counts_all_samples_peaks_renamed.txt",header = T) # bed file with start and end of peaks
counts$Chr=as.character(counts$Chr)
