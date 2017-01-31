# check if sequences fall within atac-seq peaks

# read in queries

exon4=read.delim(file="/Users/Marta/Documents/WTCHG/DPhil/Katia_012017/2017-01-25\ RREB1\ ex4\ off\ targets.txt",sep="\t",header = T,fileEncoding="UTF-16LE")
exon10=read.delim(file="/Users/Marta/Documents/WTCHG/DPhil/Katia_012017/2017-01-25\ RREB1\ ex10\ off\ targets.txt",sep="\t",header = T,fileEncoding="UTF-16LE")
exon4$chromosome=as.character(exon4$chromosome)
exon10$chromosome=as.character(exon10$chromosome)

# read in peaks file
counts <- read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/counts_all_samples_peaks_renamed.txt",header = T) # bed file with start and end of peaks
counts$Chr=as.character(counts$Chr)


exon4_peaks <- list()
exon10_peaks <- list()
for(n in 1:nrow(exon4)){
  x=exon4[n,]
  exon4_peaks[[n]]=subset(counts, x$chromosome == Chr & (x$chrstart <= End & Start <= x$chrend))
  
}

exon4_peaks <- do.call("rbind", exon4_peaks)

for(n in 1:nrow(exon10)){
  x=exon10[n,]
  exon10_peaks[[n]]=subset(counts, x$chromosome == Chr & (x$chrstart <= End & Start <= x$chrend))
  
}
exon10_peaks <- do.call("rbind", exon10_peaks)

write.table(exon4_peaks,file="/Users/Marta/Documents/WTCHG/DPhil/Katia_012017/exon4_peaks.txt",sep="\t",row.names = F)
write.table(exon10_peaks,file="/Users/Marta/Documents/WTCHG/DPhil/Katia_012017/exon10_peaks.txt",sep="\t",row.names = F)

# two regions intersect if there is at least 1 bp  where x1 <= bp <= x2 and y1 <= bp <= y2, therefore  x1 <= y2 && y1 <= x2



# the peaks are narrow-ish (300 bp- 1000 bp), so I sould add 200 bp around them and check again

counts_extended=counts
counts_extended$Start=counts_extended$Start - 200
counts_extended$End=counts_extended$End + 200


exon4_peaks_extended <- list()
exon10_peaks_extended <- list()
for(n in 1:nrow(exon4)){
  x=exon4[n,]
  exon4_peaks_extended[[n]]=subset(counts_extended, x$chromosome == Chr & (x$chrstart <= End & Start <= x$chrend))
  
}

exon4_peaks_extended <- do.call("rbind", exon4_peaks_extended)


for(n in 1:nrow(exon10)){
  x=exon10[n,]
  exon10_peaks_extended[[n]]=subset(counts_extended, x$chromosome == Chr & (x$chrstart <= End & Start <= x$chrend))
  
}

exon10_peaks_extended <- do.call("rbind", exon10_peaks_extended)



allpeaks=read.delim(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/merged_narrowpeaks.sorted.bed",sep="\t",header=F)


exon4_peaks_chrX <- list()
exon10_peaks_chrX <- list()

for(n in 1:nrow(exon4)){
  exon4_peaks_chrX[[n]]=allpeaks[which(allpeaks$V1=="chrX" & ((exon4[n,2]<=(allpeaks$V3)) & ((allpeaks$V2)<=exon4[n,3]))),]
}
exon4_peaks_chrX <- do.call("rbind", exon4_peaks_chrX)


for(n in 1:nrow(exon10)){
  exon10_peaks_chrX[[n]]=allpeaks[which(allpeaks$V1=="chrX" & ((exon10[n,2]<=(allpeaks$V3)) & ((allpeaks$V2)<=exon10[n,3]))),]
}

exon10_peaks_chrX <- do.call("rbind", exon10_peaks_chrX)

# WARNING
# it's picking up most of the peaks correctly, but also includes peaks within coordinates in other chromosomes!!



