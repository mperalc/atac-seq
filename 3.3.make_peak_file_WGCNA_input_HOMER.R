# make peak file for homer from peak matrix for WGCNA
# rename 3.4
# after running fGWAS

library(readr)

dir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/WGCNA/WGCNA_P12_S120_deepSplit2_signedHybrid_rmOutliersT/"
peaks = read_delim(file = paste(dir,"ATAC_pval_selected_MM.txt", sep = ""),
                   delim = "\t",
                   col_names = T)
# read in expression/atac data
# ATAC
load(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/session_objects/dge_atac-seq_10CPM_trim_conservativeCounts_allQualitySamples.xz"
)
colnames(dge$genes) = c("gene","Chr",  "Start",    "End", "Strand", "Length")
merged = merge(peaks,dge$genes, by = "gene")
rownames(merged) = merged$gene
merged = merged[,c(3,4,5,2,1)]
modules_selected = c("magenta", "pink", "black", "turquoise", "grey") # change every time according to significantly enriched modules

############################################################

colnames(merged) = c("Chr", "Start", "End", "Module","Name")

p = list()
for (s in modules_selected) {
  p[[s]] = merged[which(merged$Module == s), c("Name", "Chr", "Start", "End")]
  p[[s]]["Strand"] = rep(".", nrow(p[[s]]))
  write.table(
    p[[s]],
    file = paste(
      dir,"/HOMER/HOMER_input/",
      s,
      "_peaks_WGCNA_forHOMER.txt",
      sep = ""
    ),
    quote = F,
    row.names = F,
    col.names = F,
    sep = "\t"
  )
}

### annotate per module


library(ChIPpeakAnno)
library(biomaRt)
data(TSS.human.GRCh37) # load TSS info from CRCh37

# range data and annotate position with respect to closest TSS
# This way I can see if it's in a promoter, inside a gene or in an intergenic region

range_data = function(data) {
  rangedpeak = RangedData(
    IRanges(start = data$Start, end = data$End),
    # passing peak tables as ranged data
    names = data$Name,
    space = data$Chr
  )
  return(rangedpeak)
}

for (s in names(p)) {
  rangedPeaks = range_data(p[[s]])
  
  # Two methods of getting the nearest gene, both could be valid in principle
  
  annotatedPeakNearestLoc = annotatePeakInBatch(rangedPeaks, AnnotationData =
                                                  TSS.human.GRCh37) # annotate with hg19.
  annotatedPeakShortestDistance = annotatePeakInBatch(rangedPeaks, AnnotationData =
                                                        TSS.human.GRCh37, output = "shortestDistance")
  # "shortestDistance" will output nearest features to peaks
  # Displays nearest features to TSS
  
  annotatedPeakNearestLoc = as.data.frame(annotatedPeakNearestLoc)
  annotatedPeakShortestDistance = as.data.frame(annotatedPeakShortestDistance)
  
  write.table(
    annotatedPeakNearestLoc,
    file = paste(
      dir,
      s,
      "_annotatedPeakNearestLoc.txt",
      sep = ""
    ),
    sep = "\t",
    quote = F,
    row.names = F,
    col.names = T
  )
  
  write.table(
    annotatedPeakShortestDistance,
    file = paste(
      dir,
      s,
      "_annotatedPeakShortestDistance.txt",
      sep = ""
    ),
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
    file = paste(
      dir,
      s,
      "_annotatedPeakOverlappingandNearest.txt",
      sep = ""
    ),
    sep = "\t",
    quote = F,
    row.names = F,
    col.names = T
  )
  
}