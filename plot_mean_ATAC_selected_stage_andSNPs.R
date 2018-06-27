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
source("/Users/Marta/Documents/WTCHG/R scripts/locString2GRanges.R")
currentDate <- Sys.Date() # to save date in name of output files

############### variable ##########

region <- c(2,60750630,60780633)  # Chr, position start, position finish

stage = "DE"


stages = c("iPSC","DE","GT","PF","PE","EP","EN","BLC")

#gene_of_interest = "RREB1"



########### plot x axis

region.gr = locString2GRanges(paste("chr", region[1],":", region[2], "-", region[3], sep = ""))

gtrack = GenomeAxisTrack(range =  region.gr ) # plot selected region


##################### 1. prepare ATAC data
prepareATAC = function(){
# QC'd ATAC data
load(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/session_objects/dge_atac-seq_10CPM_trim.xz"
)

# Normalize ATAC data
datATAC = cbind(dge$genes[1:4], dge$counts)
datATAC[5:ncol(datATAC)]  = cpm(datATAC[5:ncol(datATAC)])

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
# GWAS <- read.table( "/Users/Marta/Documents/WTCHG/DPhil/Data/GWAS_list/RREB1_rs9505097_rs35742417_credible_set.txt", 
#                     header=T,check.names=F,sep="\t")

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

 ######## create data track
 
 dtrack <- DataTrack(
   data = dat[-1],
   start = coords[-length(coords)],
   end = coords[-1],
   chromosome = paste("chr", region[1], sep = ""),
   genome = "hg19",
   name = "mean open chromatin (CPM)" )
 
 #### plot all tracks  
 plotTracks(list(gtrack, dtrack),
            from = region[2] ,
            to = region[3] , type = "histogram")

