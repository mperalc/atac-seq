# read in ATAC peak count matrix from featureCounts
#  step 1 after counting the ATAC-seq reads with featureCounts
#  unfiltered (no QC after coming out of peak calling)
# renames the columns and saves the file


counts <- read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/peak_calling/IDR_peaks_from_bam_bds_pipeline/conservative_set/counts_all_samples_peaks_narrowpeaks_conservative_set_highQualityBams.txt",
                     header = T,
                     check.names=F,
                     skip=1) # file with feature counts & info of start and end of peaks (all peaks present in at least one sample as detected by MACS2)

colnames(counts)[7:ncol(counts)]=paste(rep(c("A","B","C","D"),each=6),
                                       rep(as.character(c(1:6)),4),
                                       sep="")  # rename columns for samples and stages

sample=factor(c("sbad2.1","sbad3.1","neo1.1"),levels=c("sbad2.1","sbad3.1","neo1.1"))   # data comes from three donors

stage=factor(c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC"),levels=c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC")) # 8 stages

rownames(counts)=counts[,1]  # peaks are row names
counts=counts[,-1]
names <- read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/names_samples.txt",header = T)

names$stage_sample <- paste(names$stage,names$sample,sep="-")

colnames(counts)[c(6:ncol(counts))]=names$stage_sample


# sort by stage and sample
source("/Users/Marta/Documents/WTCHG/R scripts/sort_by_2_strings.R")
counts=cbind(counts[1:5],sort_by_2_strings(counts[c(6:ncol(counts))],stage,sample))

write.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/peak_calling/IDR_peaks_from_bam_bds_pipeline/conservative_set/counts_all_samples_peaks_narrowpeaks_conservative_set_highQualityBams_renamed.txt",
            counts,col.names = T,row.names = T,sep="\t")
