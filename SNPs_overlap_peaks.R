# overlap ATAC peaks and T2D GWAS SNPs

### load libraries
#library(biomaRt)
library(dplyr)
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)  # SNP data 
snps <- SNPlocs.Hsapiens.dbSNP142.GRCh37

#### input files
counts <- read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/counts_all_samples_peaks_narrowpeaks_conservative_set_highQualityBams_renamed.txt",
                     header = T) # bed file with start and end of peaks
credset = read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/Feb_17_credible_set_T2D/DIAGRAM_T2D_1000G_99_Credsets/credible_sets_diagram_forMAGENTA.txt",
                     header = F) # 99 credible sets T2D 2017 DIAGRAM
nrow(counts)
# remove sex chromosomes
r1=rownames(counts[grep("chrY", rownames(counts), value=F),])
r2=rownames(counts[grep("chrX", rownames(counts), value=F),])

counts=counts[!rownames(counts) %in% c(r1,r2),]
rm(r1,r2)
nrow(counts)

########### script

counts$Chr=gsub("\\chr*","",counts$Chr) #remove "chr" in every element in location

# load ensembl
# ensembl = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_snp")
# # get regions to query in correct format
chr.region=paste(counts$Chr,counts$Start,counts$End,sep = ":")
filterlist <- list(as.character(chr.region))
# # biomart query
# results <- getBM(attributes = c("refsnp_id",'chr_name', "chrom_start", "chrom_end"),
#                  filters = c("chromosomal_region"),values = filterlist, mart = ensembl)

# biomart queries this big take too long


mycoords.gr = lapply(filterlist, function (x) {res=strsplit(x, ':')}) %>%
  unlist %>%
  as.numeric %>%
  matrix(ncol=3, byrow=T) %>%
  as.data.frame %>%
  select(chrom=V1, start=V2, end=V3) %>%
  mutate(chrom=paste0('ch', chrom)) %>%   # ch and not chr so that seq levels are equal in the snp object and this one
  makeGRangesFromDataFrame

## Get all SNPs overlapping some genomic region of interest:
test=snpsByOverlaps(snps,  mycoords.gr)


rm(test)
