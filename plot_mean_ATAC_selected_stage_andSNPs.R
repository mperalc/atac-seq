# several ATAC plots + SNP info

# a variation of plot_meth_ATAC

# mainly to check regulatory data from various sources

############## libraries

library(ggplot2)
library(readr)
library(reshape2)
library(grid)
library(gridExtra)
library(gtable)
library(ChAMPdata)
library(digest)
library(Homo.sapiens)
library(dplyr)
library(biomaRt)
library(ggbio)
library(Gviz)
library(GenomicFeatures)
library(stringr)
source("/Users/Marta/Documents/WTCHG/R scripts/locString2GRanges.R")
currentDate <- Sys.Date() # to save date in name of output files

############### variable ##########
#region <- c(15,53744228,54075228)  # Chr, position start, position finish
region <- c(20,21453433,21496433)  # Chr, position start, position finish
#highlight = c(181427200,181428200)
top_snp_file = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/top_credset_SNPs_in_peaks_PPA.txt")
top_snp_file = unique(top_snp_file)
stage = "PE"

############ get regions
# for (t in top_snp_file$V1) {
#   region = as.numeric(unlist(str_split(t, fixed("_"))))
#   highlight = c(region[2] - 1000, region[2] + 1000)
#   region[3] = region[2] + 500000
#   region[2] = region[2] - 500000




stages = c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

#output_name = paste("SLC2A10_GT","chr",sep = "_")
output_name = paste(region[1],region[2],region[3],sep = "_")


########### plot x axis

region.gr = locString2GRanges(paste("chr", region[1],":", region[2], "-", region[3], sep = ""))

gtrack = GenomeAxisTrack(range =  region.gr, min.height = 1000 ) # plot selected region


##################### 1. prepare ATAC data
prepareATAC = function(){
# QC'd ATAC data
load(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/session_objects/dge_atac-seq_10CPM_trim.xz"
)
  CPMs = read.table(file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/CPMs_atac-seq_10CPM_trim.txt")
  
# Normalize ATAC data
datATAC = cbind(dge$genes[1:4], CPMs)

# get peaks in region
peaks.gr = makeGRangesFromDataFrame(datATAC[1:4])

# input must be in char string format
region.gr = locString2GRanges(paste("chr", region[1],":", region[2], "-", region[3], sep = ""))

peak_subset = subsetByOverlaps(peaks.gr, region.gr)

# get peak counts of subset
datATAC = datATAC[which(datATAC$Geneid %in% rownames(as.data.frame(peak_subset))),]

if(nrow(datATAC) == 0){
  message("No overlapping peaks in the region. Aborting!")
  break()
}
nc = datATAC[, grepl(stage, names(datATAC))]  #takes columns whose name matches the selected stage

datATAC = cbind(datATAC[1:4],rowMeans(nc))
colnames(datATAC)[5] = "mean_CPM"

return(datATAC)
}

peak_subset = prepareATAC()

############# 2. plot ATAC data

# give fake floor of 0 to the peaks to plot nicely
fakeFloor = function(){
  
  # ifelse is correcting unnecesary zero points

  p_m = peak_subset$Start - 1 # start of peaks -1
  p_p = peak_subset$End + 1  # end of peaks +1
  
  floor = cbind(c(p_p, region[2],p_m, region[3]),
                 rep(0,length(c(p_p, region[2],p_m, region[3]))))
  colnames(floor) = c("pos","mean_CPM")
  
  # adding peak data
  floor = rbind(floor,
                cbind(peak_subset$Start,peak_subset$mean_CPM),
                cbind(peak_subset$End,peak_subset$mean_CPM))
  floor = as.data.frame(floor)
  floor = floor[order(floor$pos),] # order by location
  return(floor)
}

pointsToPlot = fakeFloor()

 dat <- pointsToPlot$mean_CPM
 coords = pointsToPlot$pos

 ################ bedformat
 # 
 # bamfile = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/bam/PE_filtered_bams_merged.sorted.bam"
 # 
 # 
 # alTrack <- AlignmentsTrack(range = bamfile, 
 #                            isPaired = TRUE,
 #                            chromosome = paste("chr",region[1],sep=""),
 #                            from = region[2],to = region[3] )
 # 
 # plotTracks(
 #   list(alTrack, gtrack),
 #   chromosome = region[1],
 #   from = region[2] ,
 #   to = region[3] ,
 #   name = ""
 # )
 # 
 # 
 ################## ATAC track #################
 
 ATACtrack <- DataTrack(
   data = dat[-1],
   start = coords[-length(coords)],
   end = coords[-1],
   chromosome = paste("chr", region[1], sep = ""),
   genome = "hg19",
   name = "mean open chromatin (CPM)",
   type = "histogram",
   col.histogram = "royalblue4",
   fill.histogram = "royalblue4")
 
 ##################### gene track #######


 library(TxDb.Hsapiens.UCSC.hg19.knownGene)
 txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

 txTr <-
   GeneRegionTrack(txdb,
                   chromosome = paste("chr", region[1], sep = ""),
                   start = region[2],
                   end = region[3],
                   name = "genes",
                   transcriptAnnotation = "symbol",
                   collapseTranscripts = "meta")  # collapses all transcripts, showing a "metatranscript" with all exons
 
 # to substitute UCSC gene symbols for ensembl symbols
 symbols <- unlist(mapIds(org.Hs.eg.db, gene(txTr), "SYMBOL", "ENTREZID", multiVals = "first"))
 symbol(txTr) <- symbols[gene(txTr)]
 
 
 # biomart alternative
 # 
 # 
 # biomTrack <- BiomartGeneRegionTrack(
 #   genome = "hg19",
 #   chromosome = paste("chr", region[1], sep = ""),
 #   start = region[2],
 #   end = region[3],
 #   name = "ENSEMBL",
 #   filter = list(biotype = "protein_coding"),
 #   collapseTranscripts = TRUE
 #   #stacking = "squish"
 # )
 # 
 # 
 # plotTracks(
 #   biomTrack,
 #   col.line = NULL,
 #   col = NULL,
 #   stackHeight = 0.3,
 #   transcriptAnnotation = "symbol"
 # )
 # 
 
 ######################### GWAS credible set track ###############
 
 # if the region overlaps any credible set

 prepareCredSet = function(){
   
   credset <-
     read.table(
       "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/HRC.credible.set.bed",
       header = F,
       check.names = F,
       skip = 1,
       sep = "\t"
     )
   colnames(credset) = c("chr","start","end")
   
   credset.gr = makeGRangesFromDataFrame(credset)
   region.gr = locString2GRanges(paste("chr", region[1],":", region[2], "-", region[3], sep = ""))
   
   credset_subset = subsetByOverlaps(credset.gr, region.gr)
   
   # if (nrow(credset_subset) == 0) {
   #   message("No overlapping SNPs in the region. Aborting!")
   #   break()
   # }
  
   # get SNPs overlapped by peaks to give them specific colors
   load(
     "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/session_objects/dge_atac-seq_10CPM_trim.xz"
   )
   
   datATAC = dge$genes[1:4]
   peaks.gr = makeGRangesFromDataFrame(datATAC)
   overlaps = findOverlaps(credset_subset, peaks.gr)
   
   # initialize dataframe with colors
   cols = data.frame(credset_subset)
   cols$cols = rep("royalblue4", nrow(cols))
   values(credset_subset) = cols$cols
   values(credset_subset[overlaps@from]) = "red" 
   
   return(credset_subset)
 }
 
credset_subset = prepareCredSet()
# credset_subset = as.data.frame(credset_subset)
cols = as.data.frame(credset_subset)$X

credsetTrack  <-
  AnnotationTrack(
    credset_subset,
    name = "credible set SNPs",
    stacking = "dense",
    genome = "hg19",
    fill = cols
  )
 

################## Wang chromHMM data #########
if (stage %in% c("DE","GT","PF","PE")){
prepareChromHMM = function(){
  
  if (stage == "DE") {
    cHMM <-
      read.table(
        "/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/ChromHMM_islet_diff_stages_Wang_2015/DE_14_dense.bed",
        header = F,
        check.names = F,
        skip = 1,
        sep = "\t"
      )
  }
  if (stage == "GT") {
    cHMM <-
      read.table(
        "/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/ChromHMM_islet_diff_stages_Wang_2015/GT_14_dense.bed",
        header = F,
        check.names = F,
        skip = 1,
        sep = "\t"
      )
    
  }
  if (stage == "PF") {
    cHMM <-
      read.table(
        "/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/ChromHMM_islet_diff_stages_Wang_2015/PF_14_dense.bed",
        header = F,
        check.names = F,
        skip = 1,
        sep = "\t"
      )
    
  }
  if (stage == "PE") {
    cHMM <-
      read.table(
        "/Users/Marta/Documents/WTCHG/DPhil/Data/Reference/ChromHMM_islet_diff_stages_Wang_2015/PE_14_dense.bed",
        header = F,
        check.names = F,
        skip = 1,
        sep = "\t"
      )
  }
    colnames(cHMM) = c("chr","start","finish","state","x","y","start2","finish2","rgb")
    
    cHMM.gr = makeGRangesFromDataFrame(df = cHMM, keep.extra.columns = T,
                                          seqnames.field = "chr" ,start.field = "start", end.field = "finish")
    region.gr = locString2GRanges(paste("chr", region[1],":", region[2], "-", region[3], sep = ""))
    
    cHMM_subset = subsetByOverlaps(cHMM.gr, region.gr)
    cHMM_subset$rgb = as.factor(cHMM_subset$rgb)
    cHMM_subset$rgb = factor(cHMM_subset$rgb, levels = c(levels(cHMM_subset$rgb),
                                                         c("seagreen2" ,"mediumblue","steelblue3","seagreen4","springgreen3","red4",
                                                           "midnightblue","magenta3","olivedrab4","purple3","khaki1","plum2","yellow2",
                                                           "lightgoldenrod2")))
# 
#     # change colors
#     values(cHMM_subset[cHMM_subset$rgb == "51,255,153","rgb"]) = factor("seagreen2")
#     values(cHMM_subset[cHMM_subset$rgb == "0,0,255","rgb"]) = factor("mediumblue" )
#     values(cHMM_subset[cHMM_subset$rgb == "0,153,204","rgb"]) = factor("steelblue3")
#     values(cHMM_subset[cHMM_subset$rgb == "0,102,0","rgb"]) = factor("seagreen4" )
#     values(cHMM_subset[cHMM_subset$rgb == "0,204,51","rgb"]) = factor("springgreen3" )
#     values(cHMM_subset[cHMM_subset$rgb == "204,0,51","rgb"]) = factor("red4")
# 
#     values(cHMM_subset[cHMM_subset$rgb == "0,0,51","rgb"]) = factor("midnightblue")
#     values(cHMM_subset[cHMM_subset$rgb == "255,0,204","rgb"]) = factor("magenta3" )
#     values(cHMM_subset[cHMM_subset$rgb == "102,153,51","rgb"]) = factor("olivedrab4")
#     values(cHMM_subset[cHMM_subset$rgb == "102,51,153","rgb"]) = factor("purple3" )
#     values(cHMM_subset[cHMM_subset$rgb == "255,255,204","rgb"]) = factor("khaki1")
#     values(cHMM_subset[cHMM_subset$rgb == "204,153,255","rgb"]) = factor("plum2")
#     values(cHMM_subset[cHMM_subset$rgb == "255,255,0","rgb"]) = factor("yellow2" )
#     values(cHMM_subset[cHMM_subset$rgb == "204,204,102","rgb"]) = factor("lightgoldenrod2")


    cHMM_subset = as.data.frame(cHMM_subset)
    cHMM_subset$rgb = as.character(cHMM_subset$rgb)

    # change colors
   cHMM_subset[cHMM_subset$rgb == "51,255,153","rgb"] =  "seagreen2"
     cHMM_subset[cHMM_subset$rgb == "0,0,255","rgb"] =  "mediumblue" 
     cHMM_subset[cHMM_subset$rgb == "0,153,204","rgb"] =  "steelblue3"
     cHMM_subset[cHMM_subset$rgb == "0,102,0","rgb"] =  "seagreen4" 
     cHMM_subset[cHMM_subset$rgb == "0,204,51","rgb"] =  "springgreen3" 
     cHMM_subset[cHMM_subset$rgb == "204,0,51","rgb"] =  "red4"

     cHMM_subset[cHMM_subset$rgb == "0,0,51","rgb"] =  "midnightblue"
     cHMM_subset[cHMM_subset$rgb == "255,0,204","rgb"] =  "magenta3" 
     cHMM_subset[cHMM_subset$rgb == "102,153,51","rgb"] =  "olivedrab4"
     cHMM_subset[cHMM_subset$rgb == "102,51,153","rgb"] =  "purple3" 
     cHMM_subset[cHMM_subset$rgb == "255,255,204","rgb"] =  "khaki1"
     cHMM_subset[cHMM_subset$rgb == "204,153,255","rgb"] =  "plum2"
     cHMM_subset[cHMM_subset$rgb == "255,255,0","rgb"] =  "yellow2" 
     cHMM_subset[cHMM_subset$rgb == "204,204,102","rgb"] =  "lightgoldenrod2"

  return(cHMM_subset)
}

cHMM_subset = prepareChromHMM()


# Broad1 <-
#   UcscTrack(
#     track = paste(stage,'ChromHMM',sep = " "),
#     table = "wgEncodeBroadHmmGm12878HMM",
#     trackType = "AnnotationTrack",
#     genome = 'hg18',
#     chromosome = 'chr18',
#     name = '12878',
#     from = 44675486,
#     to = 44679944,
#     start = "chromStart",
#     end = "chromEnd",
#     feature = "itemRgb",
#     id = "name",
#     collapse = FALSE,
#     stacking = "dense"
#   )

#### annotations and colors are wrong!!!1!!!


#feat <- unique(feature(Broad1)) featCol <- setNames(as.list(rgb(t(sapply(strsplit(feat, ","), as.numeric)), maxColorValue=255)), feat)

# cHMMTrack = AnnotationTrack(cHMM_subset,
#                             name = 'ChromHMM',
#                             stacking = "dense",
#                             collapse = F,
#                             genome = "hg19",
#                             fill = cols,
#                             featureAnnotation = "state",
#                             fontcolor.feature = "darkblue")

cHMM_subset = as.data.frame(cHMM_subset)

state = factor(cHMM_subset$state, levels = unique(cHMM_subset$state) , ordered = T)
cols = factor(cHMM_subset$rgb, levels = unique(cHMM_subset$rgb) , ordered = T)

cHMMTrack = AnnotationTrack(start = cHMM_subset$start,
                            end = cHMM_subset$end,
chromosome = paste("chr", region[1], sep = ""),
genome = "hg19",
fill = cols,
feature = state,
featureAnnotation = "feature",
fontcolor.feature = "darkblue")

}
############ highlight specific region #################
# 
# if (stage %in% c("DE","GT","PF","PE")) {
#   ht <- HighlightTrack(trackList = list(ATACtrack, txTr, credsetTrack, cHMMTrack), start = highlight[1], end = highlight[2],
#                      chromosome = region[1])
# } else{
#   ht <- HighlightTrack(trackList = list(ATACtrack, txTr, credsetTrack), start = highlight[1], end = highlight[2],
#                        chromosome = region[1])
# }


tiff(
  paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/peak_profiles/",
        stage,"_",
        output_name,
        ".tiff",
        sep = ""),
  type = "cairo",
  width = 4.5,
  height = 6,
  units = "in",
  res = 1000,
  pointsize = 12,
  compression = "none"
)
# plotTracks(list(ht, gtrack),
#            from = region[2] ,
#            to = region[3] ,
#            name = "")

plotTracks(list(ATACtrack, txTr, credsetTrack, cHMMTrack, gtrack),
           from = region[2] ,
           to = region[3] ,
           name = "")
dev.off()


 
