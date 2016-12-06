# length of regions in bed peaks file

setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq")


currentDate <- Sys.Date() # to save date in name of output files


table <- read.table(file="merged_narrowpeaks.sorted.bed") # bed file with start and end of peaks

red_table <- table[c(2,3)]

red_table$lengths=red_table$V3-red_table$V2

total_coverage=sum(red_table$lengths)

# make SAF file for featureCounts

table$peaks=paste("peak",table$V1,table$V2,sep = "_")  # names of peaks = chr + start position
colnames(table)=c("Chr","Start","End","Samples","GeneID") # need to change peak for GeneID to be recognised  by featureCounts

table=table[c(5,1,2,3)] # reorder for SAF format
table$Strand=rep(".",length(table$GeneID))  # . because we don't care about the strand

table=table[which(table$Chr!="chrM" & table$Chr!="chrX" & table$Chr!="chrY"),]  # exclude those chromosomes
write.table(table,"all_peaks.saf",sep="\t",quote = F, col.names = T, row.names = F)
