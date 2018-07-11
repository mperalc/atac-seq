# peaks that overlap SNPs from credible sets
# from WGCNA HOMER TF annotation results
library(GenomicRanges)
library(dplyr)
library(reshape2)

modules = c("blue", "green", "lightyellow", "pink", "purple", "royalblue")

output = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/WGCNA/10CPM_P12_S120_deepSplit2_signedHybrid_noOutliers/HOMER/HOMER_output/"

credset = read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/HRC.credible.set.bed",
                               header = F) # 99% credible sets T2D 2017 HRC (Mahajan et al 2018)


####################

anyDuplicated(credset) # duplicated SNPs from different loci?
credset = unique(credset)
colnames(credset) = c("chr","start","end")
credset$chr = as.character(credset$chr)
rownames(credset) = paste("snp",credset$chr,credset$start, credset$end,sep = "_")
nrow(credset)

credset.gr = credset %>%
  makeGRangesFromDataFrame

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
  input.gr = input[,c("Chr","Start","End")]
  
  input.gr  = makeGRangesFromDataFrame(input.gr )
  
  # Get all credible set SNPs overlapping ATAC peaks in module:
  overlapSNP = subsetByOverlaps(credset.gr, input.gr)
  
  # Get all peaks in module overlapping credible set SNPs :
    overlapPeak = subsetByOverlaps(input.gr, credset.gr)
    hits <- findOverlaps( input.gr, credset.gr) # to add SNP metadata to overlapPeak table
    rsid <- CharacterList(split(names(credset.gr)[subjectHits(hits)], queryHits(hits))) 
    rsid = paste(rsid,collapse = ",")
    
  mcols(overlapPeak) <- DataFrame(mcols(overlapPeak), rsid) 
  
  # re-shaping input file
  colnames(input)[1] = "PeakID"
  
  input = input[, c(1:4,8,10,15,16, 22:ncol(input))]
  
  overlapPeak = as.data.frame(overlapPeak)
  overlapPeak = overlapPeak[,c(1:3,6)]
  colnames(overlapPeak) = c("Chr", "Start", "End","SNPpos")
  overlapPeak$PeakID = paste("peak",overlapPeak$Chr,overlapPeak$Start,sep = "_")
  
  input = input[which(input$PeakID %in% overlapPeak$PeakID),]
  input = merge(overlapPeak[c("PeakID","SNPpos")],input,by = "PeakID")
  # put in long format
  melted = melt(
    input,
    id.vars = c("SNPpos","PeakID", "Chr", "Start", "End","Annotation","Distance.to.TSS","Nearest.Ensembl", "Gene.Name"),
    na.rm = T
  )
  
  colnames(melted)[c(6, 7)] = c("Motif HOMER", "Sequence found")
  write.table(
    melted,
    file = paste(
      output,
      module,
      "/homerResults/allMotifs_annotatePeaks_peaksOverCredsets.txt",
      sep = ""
    ),
    quote = F,
    row.names = F,
    col.names = T,
    sep = "\t"
  )
}
