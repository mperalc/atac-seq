# from merged peaks to fGWAS input
source("/Users/Marta/Documents/WTCHG/R scripts/atac-seq/peaks_merged_to_binary_df.R")

peaks=read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/peak_calling/IDR_peaks_from_bam_bds_pipeline/conservative_set/merged_narrowpeaks_conservative_set_highQualityBams.bed",
                 header = F,
                 check.names=F) # file with info of start and end of peaks, and samples that have it 
colnames(peaks)=c("Chr","Start","End","Peaks_present")

binary=peaks_merged_to_binary_df(peaks)

stage=c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC")

peak_per_stage=list()

for(s in stage){
  peak_per_stage[[s]]=binary[which(binary[s]==1),c("Chr","Start","End")]
  peak_per_stage[[s]]$Stage=rep(s,nrow(peak_per_stage[[s]]))
}
fgwas_annotation=do.call("rbind",peak_per_stage)

fgwas_annotation=fgwas_annotation[order(fgwas_annotation[,1],fgwas_annotation[,2]),]

write.table(fgwas_annotation,
            "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/peak_calling/IDR_peaks_from_bam_bds_pipeline/conservative_set/fGWAS_merged_narrowpeaks_conservative_set_highQualityBams.bed",
            row.names = F,col.names = F,quote = F)
