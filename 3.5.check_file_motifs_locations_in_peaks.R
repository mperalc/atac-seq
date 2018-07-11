# 3.5
# convert MOTIF location file into something easier to browse
library(reshape2)

# genes to examine from tables
genes = c("CDKAL1", "KCNQ1", "INS")
# WGCNA modules significant after fgwas
modules = c("blue", "green", "lightyellow", "pink", "purple", "royalblue")
output = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/WGCNA/10CPM_P12_S120_deepSplit2_signedHybrid_noOutliers/HOMER/HOMER_output/"
for (module in modules) {
  input = read.table(
    file = paste(
      output,
      module,
      "/homerResults/allMotifs_annotatePeaks.txt",
      sep = ""
    ),
    header = T,
    sep = "\t",
    fill = T,
    na.strings = "",
    stringsAsFactors = F,
    quote = "\""
  )
  
  colnames(input)[1] = "PeakID"
  
  input = input[, c(1:4, 16, 22:ncol(input))]
  # subset by genes
  input = input[which(input$Gene.Name %in% genes),]
  # put in long format
  melted = melt(
    input,
    id.vars = c("PeakID", "Chr", "Start", "End", "Gene.Name"),
    na.rm = T
  )
  
  colnames(melted)[c(6, 7)] = c("Motif HOMER", "Sequence found")
  write.table(
    melted,
    file = paste(
      output,
      module,
      "/homerResults/allMotifs_annotatePeaks_selectedGenes.txt",
      sep = ""
    ),
    quote = F,
    row.names = F,
    col.names = T,
    sep = "\t"
  )
}
