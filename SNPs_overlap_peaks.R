# overlap ATAC peaks and T2D GWAS SNPs

### load libraries
library(dplyr)
library(biomaRt)
#library(readr)
#library(vcfR)
library(GenomicRanges)
library(BSgenome) # for SNP by overlaps function

#### input files
counts <-
  read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/counts_all_samples_peaks_narrowpeaks_conservative_set_highQualityBams_renamed.txt",
             header = T) # bed file with start and end of peaks
peaks = counts[1:3]
rm(counts)

### update credset!! ##########
credset = read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/HRC.credible.set.bed",
                     header = F) # 99 credible sets T2D 2017 DIAGRAM
credset = unique(credset[1:2]) # SNPs with unique positions (some are shared between loci, because they belong to different independent signals and have different PPAs)
#credset = read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/RREB1_rs9505097_rs35742417_credible_set.txt",
#          header = T)
colnames(credset) = c("Chr", "Pos")
credset = credset[, c("Chr", "Pos")]
nrow(peaks)
any(duplicated(peaks)) # there are no duplicated peaks
#[1] FALSE
nrow(credset)

# save the step of assigning SNP id for after intersecting peaks
# start_time <- Sys.time()
# snps = read.vcfR( "/Users/Marta/Documents/WTCHG/DPhil/Data/SNPs_b150_GRCh37p13/common_all_20170710.vcf.gz", verbose = FALSE ) # common SNPs build 150 GRCh37p13
# end_time - start_time

# remove sex chromosomes
r1 = rownames(peaks[grep("chrY", rownames(peaks), value = F),])
r2 = rownames(peaks[grep("chrX", rownames(peaks), value = F),])

peaks = peaks[!rownames(peaks) %in% c(r1, r2),]
rm(r1, r2)
nrow(peaks)


colnames(credset) = c("chr", "start")
credset$end = credset$start

credset.gr = credset %>%
 # mutate(chr = paste0('chr', chr)) %>%
  makeGRangesFromDataFrame

peaks.gr = makeGRangesFromDataFrame(peaks)

# Get all credible set SNPs overlapping ATAC peaks:
overlapSNP = subsetByOverlaps(credset.gr, peaks.gr)


# Get all peaks overlapping credible set SNPs :
overlapPeak = subsetByOverlaps(peaks.gr, credset.gr)
#overlapPeak=as.data.frame(overlapPeak)


write.table(
  as.data.frame(overlapPeak),
  file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/Peaks_over_SNPs.txt",
  sep = "\t",
  quote = F,
  row.names = T,
  col.names = T
)

# subset peaks over 10CPM
load(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/session_objects/dge_atac-seq_10CPM_trim.xz", verbose =
       TRUE)  #loading again the dge object to do the selected QC

peaksToSubset = as.data.frame(overlapPeak)
peaksToSubset = peaksToSubset[which(rownames(peaksToSubset) %in% dge$genes$Geneid), ]


write.table(
  peaksToSubset,
  file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/over10CPM_peaks_over_SNPs.txt",
  sep = "\t",
  quote = F,
  row.names = T,
  col.names = T
)


write.table(
  as.data.frame(overlapSNP),
  file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/SNPs_over_peaks.txt",
  sep = "\t",
  quote = F,
  row.names = T,
  col.names = T
)

# more SNPs overlaping than peaks, so there must be several SNPs within each peak in some cases

# % of SNPs from credset in peaks
(length(overlapSNP) / nrow(credset)) * 100

# % of peaks with intersecting SNPs from credset
(length(overlapPeak) / nrow(peaks)) * 100

#### Annotate by promoter, coding, etc
#
#

source('~/WTCHG/R scripts/atac-seq/annotate_GRanges_TxDb.R') # loading whole genome annotation function

# ATAC peaks and GWAS credible set SNPs
allPeaksAnnotated = annotate_GRanges_TxDb(peaks.gr)
allPeaksAnnotated$anno = rownames(allPeaksAnnotated) # shaping for ggplot
allPeaksAnnotated$type = "All ATAC-seq peaks"
peaksOverGWASAnnotated = annotate_GRanges_TxDb(overlapPeak)
peaksOverGWASAnnotated$anno = rownames(peaksOverGWASAnnotated)
peaksOverGWASAnnotated$type = "ATAC-seq peaks overlapping credible sets"
GWASAnnotated = annotate_GRanges_TxDb(credset.gr)
GWASAnnotated$anno = rownames(GWASAnnotated)
GWASAnnotated$type = "All SNPs in credible sets"

source('~/WTCHG/R scripts/annotate_whole_genome_TxDb.R') # loading whole genome annotation function

whole = annotate_whole_genome_TxDb() # whole genome annotation
whole$anno = rownames(whole)
whole$type = "Whole genome"

combined = rbind(allPeaksAnnotated,
                 peaksOverGWASAnnotated,
                 GWASAnnotated,
                 whole)

combined$anno = factor(
  combined$anno,
  levels = unique(combined$anno),
  ordered = T,
  labels = c("Promoter", "Exon", "Intron", "Intergenic")
) # ordering annotations
combined$type = factor(combined$type,
                       levels = unique(combined$type),
                       ordered = T) # ordering labels of data type

# colour scale for plot

pal = c("#96165a", "#bf3b5c", "#3D8DC1", "#95C8EB")   # red and blue scale

source('~/WTCHG/R scripts/plotting/ggplot_stacked_barplot_percent.R') # loading barplot function

p = ggplot_stacked_barplot_percent(combined, pal, 0.75)

png(
  "/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/TxDb_annotation_peaks_credible_set_T2D_DIAGRAM.png",
  type = "cairo",
  width = 14,
  height = 8,
  units = "in",
  res = 500,
  pointsize = 12
)
p
dev.off()


## save peaks that overlap GWAS, and get nearest gene

library(ChIPpeakAnno)
data(TSS.human.GRCh37) # load TSS info from CRCh37

# range data and annotate position with respect to closest TSS
# This way I can see if it's in a promoter, inside a gene or in an intergenic region

range_data = function(data) {
  rangedpeak = RangedData(
    IRanges(start = data$Start, end = data$End),
    # passing peak tables as ranged data
    names = rownames(data),
    space = data$Chr
  )
  return(rangedpeak)
}

rangedPeaks = range_data(peaks)

# Two methods of getting the nearest gene, both could be valid in principle

annotatedPeakNearestLoc = annotatePeakInBatch(overlapPeak, AnnotationData =
                                                TSS.human.GRCh37) # annotate with hg19.
annotatedPeakShortestDistance = annotatePeakInBatch(overlapPeak, AnnotationData =
                                                      TSS.human.GRCh37, output = "shortestDistance")
# "shortestDistance" will output nearest features to peaks
# Displays nearest features to TSS

annotatedPeakNearestLoc = as.data.frame(annotatedPeakNearestLoc)
annotatedPeakShortestDistance = as.data.frame(annotatedPeakShortestDistance)

write.table(
  annotatedPeakNearestLoc,
  file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/Peaks_over_SNPs_annotatedPeakNearestLoc.txt",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

write.table(
  annotatedPeakShortestDistance,
  file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/Peaks_over_SNPs_annotatedPeakShortestDistance.txt",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

# annotate with gene names
# load ensembl
ensembl_genes = useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  host = "grch37.ensembl.org",
  path = "/biomart/martservice",
  dataset = "hsapiens_gene_ensembl"
)

annotated_peaks_closeGenes = rbind(annotatedPeakNearestLoc, annotatedPeakShortestDistance)
ensemblID = annotated_peaks_closeGenes$feature

results = getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters = "ensembl_gene_id",
  values = ensemblID,
  mart = ensembl_genes
)
colnames(annotated_peaks_closeGenes)[7] = "ensembl_gene_id"
results = merge(annotated_peaks_closeGenes, results, by = "ensembl_gene_id")
results = results[c(
  "peak" ,
  "seqnames",
  "start",
  "end",
  "width",
  "ensembl_gene_id",
  "external_gene_name",
  "start_position",
  "end_position",
  "feature_strand",
  "insideFeature",
  "distancetoFeature",
  "shortestDistance",
  "fromOverlappingOrNearest"
)]

write.table(
  results,
  file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/Peaks_over_SNPs_annotatedPeakOverlappingandNearest.txt",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

#### format overlapping peaks for HOMER
overlapPeak_homer = as.data.frame(overlapPeak)
overlapPeak_homer$peakID = rownames(overlapPeak_homer)
overlapPeak_homer = overlapPeak_homer[c(6, 1, 2, 3, 5)]
overlapPeak_homer$strand = rep(".", nrow(overlapPeak_homer))
write.table(
  overlapPeak_homer,
  file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/Peaks_over_SNPs_forHOMER.txt",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = F
)

overlapPeak_homer = as.data.frame(peaksToSubset)
overlapPeak_homer$peakID = rownames(overlapPeak_homer)
overlapPeak_homer = overlapPeak_homer[c(6, 1, 2, 3, 5)]
overlapPeak_homer$strand = rep(".", nrow(overlapPeak_homer))

write.table(
  overlapPeak_homer,
  file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/over10CPM_peaks_over_SNPs_forHOMER.txt",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = F
)

## save SNPs that overlap peaks, and get nearest gene

overlapSNPdf = as.data.frame(overlapSNP)
colnames(overlapSNPdf) = c("chr", "start", "end", "width", "strand")
overlapSNPdf$chr = gsub("\\chr*", "", overlapSNPdf$chr) #remove "chr" in every element in location

# load ensembl
ensembl_SNP = useMart(
  biomart = "ENSEMBL_MART_SNP",
  host = "grch37.ensembl.org",
  path = "/biomart/martservice" ,
  dataset = "hsapiens_snp"
)
# # get regions to query in correct format

chr.regionF = paste(overlapSNPdf$chr,
                    overlapSNPdf$start,
                    overlapSNPdf$end,
                    "1",
                    sep = ":") # forward strand
#chr.regionR=paste(overlapSNPdf$chr,overlapSNPdf$start,overlapSNPdf$end, "-1", sep = ":") # reverse strand


# biomart in R is a pain of timeouts
# saving regions to query with browser
write.table(
  data.frame(chr.regionF),
  file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/overlapSNPDF.txt",
  +sep = "\t",
  quote = F,
  row.names = F,
  col.names = F
)
