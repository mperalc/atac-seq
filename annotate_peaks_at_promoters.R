# annotate peaks at promoters
# all peaks and remaining after 10CPM filter
# previous step to promoter_plot_longitudinal_variation_peak_counts.R

# Load necessary libraries
library(readr)  # to read big tables fast
library(ggplot2)
library(reshape2)  # to modify dataframes for ggplot2
library(GenomicRanges)
library(ChIPpeakAnno)
data(TSS.human.GRCh37) # load TSS info from CRCh37
library(biomaRt)
currentDate <- Sys.Date() # to save date in name of output files

#####################################
### no CPM filter

# all atac-seq data
load(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/session_objects/dge_atac-seq.xz"
)
#if no CPM filter
  peaks = as.data.frame(cpm(dge$counts))

# getting location info for peaks
peaks$chr = dge$genes[which(rownames(dge$counts) %in% rownames(peaks)),"Chr"]
peaks$start = dge$genes[which(rownames(dge$counts) %in% rownames(peaks)),"Start"]
peaks$end = dge$genes[which(rownames(dge$counts) %in% rownames(peaks)),"End"]

# make into GRanges
peaks.gr = makeGRangesFromDataFrame(peaks,keep.extra.columns = T)

annotatedPeak = annotatePeakInBatch(peaks.gr, AnnotationData = TSS.human.GRCh37, output = "nearestLocation")
annotatedPeak = as.data.frame(annotatedPeak)

# select peaks at promoters
peaks_at_promoter = rbind(annotatedPeak[which(annotatedPeak$insideFeature ==
                                                "overlapStart"), ],
                          annotatedPeak[which(annotatedPeak$insideFeature ==
                                                "upstream"  & annotatedPeak$shortestDistance < 2001), ],
                          annotatedPeak[which(annotatedPeak$insideFeature ==
                                                "inside" & annotatedPeak$shortestDistance < 501), ])

# search on biomart ensembl gene names
# ensembl37 = useMart(
#   biomart = "ENSEMBL_MART_ENSEMBL",
#   host = "grch37.ensembl.org",
#   path = "/biomart/martservice" ,
#   dataset = "hsapiens_gene_ensembl"
# )
ensembl = useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  host = "www.ensembl.org",
  path = "/biomart/martservice" ,
  dataset = "hsapiens_gene_ensembl"
)

#look for gene names from ids:

all_ensembl_info_ensemblFound <-
  getBM(
    attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'),
    filters = 'ensembl_gene_id',
    values = peaks_at_promoter$feature,
    mart = ensembl
  )

# change column name for emsebl gene ids to merge with new df
colnames(peaks_at_promoter)[names(peaks_at_promoter) == "feature"] = "ensembl_gene_id"
peaks_at_promoter = merge( all_ensembl_info_ensemblFound, peaks_at_promoter,by =
                             "ensembl_gene_id") # merge old and new info


peaks_at_promoter = peaks_at_promoter[which(peaks_at_promoter$gene_biotype %in% c("protein_coding", "lincRNA")), ]


write.table(
  peaks_at_promoter,
  file = paste(
    "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/",
    currentDate,
    "annotated_peaks_at_promoters_allCPMs.txt",
    sep = ""
  ),
  sep = "\t",
  col.names = T,
  quote = F,
  row.names = F
)

#####################################
### 10CPM filter

  peaks = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/CPMs_atac-seq_10CPM_trim.txt")

# getting location info for peaks
peaks$chr = dge$genes[which(rownames(dge$counts) %in% rownames(peaks)),"Chr"]
peaks$start = dge$genes[which(rownames(dge$counts) %in% rownames(peaks)),"Start"]
peaks$end = dge$genes[which(rownames(dge$counts) %in% rownames(peaks)),"End"]

# make into GRanges
peaks.gr = makeGRangesFromDataFrame(peaks,keep.extra.columns = T)

annotatedPeak = annotatePeakInBatch(peaks.gr, AnnotationData = TSS.human.GRCh37, output = "nearestLocation")
annotatedPeak = as.data.frame(annotatedPeak)

# select peaks at promoters
peaks_at_promoter = rbind(annotatedPeak[which(annotatedPeak$insideFeature ==
                                                "overlapStart"), ],
                          annotatedPeak[which(annotatedPeak$insideFeature ==
                                                "upstream"  & annotatedPeak$shortestDistance < 2001), ],
                          annotatedPeak[which(annotatedPeak$insideFeature ==
                                                "inside" & annotatedPeak$shortestDistance < 501), ])

# search on biomart ensembl gene names
# ensembl37 = useMart(
#   biomart = "ENSEMBL_MART_ENSEMBL",
#   host = "grch37.ensembl.org",
#   path = "/biomart/martservice" ,
#   dataset = "hsapiens_gene_ensembl"
# )
ensembl = useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  host = "www.ensembl.org",
  path = "/biomart/martservice" ,
  dataset = "hsapiens_gene_ensembl"
)

#look for gene names from ids:

all_ensembl_info_ensemblFound <-
  getBM(
    attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'),
    filters = 'ensembl_gene_id',
    values = peaks_at_promoter$feature,
    mart = ensembl
  )

# change column name for emsebl gene ids to merge with new df
colnames(peaks_at_promoter)[names(peaks_at_promoter) == "feature"] = "ensembl_gene_id"
peaks_at_promoter = merge( all_ensembl_info_ensemblFound, peaks_at_promoter,by =
                             "ensembl_gene_id") # merge old and new info


peaks_at_promoter = peaks_at_promoter[which(peaks_at_promoter$gene_biotype %in% c("protein_coding", "lincRNA")), ]


write.table(
  peaks_at_promoter,
  file = paste(
    "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/",
    currentDate,
    "annotated_peaks_at_promoters_over10CPMs.txt",
    sep = ""
  ),
  sep = "\t",
  col.names = T,
  quote = F,
  row.names = F
)


