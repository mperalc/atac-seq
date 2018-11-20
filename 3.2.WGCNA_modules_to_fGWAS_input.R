# WGCNA modules to fGWAS input file
# rename 3.3

# load peak matrix file (with needed peak annotation)
peak = read.table(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/CPMs_atac-seq_10CPM_trim_conservativeCounts_allQualitySamples.txt",
  header = T
)

# with minimum module membership
Results = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/WGCNA/WGCNA_P12_S120_deepSplit2_signedHybrid_rmOutliersT/"
WGCNA = read.table(paste(Results, "ATAC_pval_selected_MM.txt", sep = ""), header = T) # less peaks than above because of MM filter

# how many modules are there
length(unique(WGCNA$module))

# match peaks in WGCNA file with their chr, start and end info.
matched = cbind(peak[match(WGCNA$gene, rownames(peak)), 1:3], WGCNA)
matched = matched[c("Chr", "Start", "End", "module")] # select needed columns
matched = matched[which(matched$Chr %in% paste("chr", 1:22, sep = "")),]
dim(matched)


write.table(
  matched,
  paste(Results,"fGWAS_WGCNA_ATAC_pval_selected_MM.bed", sep = ""),
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)
