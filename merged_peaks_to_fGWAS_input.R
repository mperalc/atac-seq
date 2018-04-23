# from merged peaks to fGWAS input
source("/Users/Marta/Documents/WTCHG/R scripts/atac-seq/peaks_merged_to_binary_df.R")
stages=c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC")


#################### from peaks directly####################################
####################### read in the file with all peaks called (for all samples), and stages that have them#############
peaks=read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/peak_calling/IDR_peaks_from_bam_bds_pipeline/conservative_set/merged_narrowpeaks_conservative_set_highQualityBams.bed",
                 header = F,
                 check.names=F) # file with info of start and end of peaks, and samples that have it 
colnames(peaks)=c("Chr","Start","End","Peaks_present")  # add header
peaks$Name=paste(peaks$Chr,peaks$Start,sep="_") # add name column needed for the binary function
peaks=peaks[c(5,1:4)] # reorder


binary=peaks_merged_to_binary_df(peaks) # this function needs columns "Name","Chr","Start","End","Peaks_present", in that order


#### to save fGWAS file with all peaks (minus those in sex chr, no filter by counts)#############################
peak_per_stage=list()

for(s in stages){ # get name of stage per peak
  peak_per_stage[[s]]=binary[which(binary[s]==1),c("Chr","Start","End")]
  peak_per_stage[[s]]$Stage=rep(s,nrow(peak_per_stage[[s]]))
}
fgwas_annotation=do.call("rbind",peak_per_stage) # bind all peaks

fgwas_annotation=fgwas_annotation[order(fgwas_annotation[,1],fgwas_annotation[,2]),]
fgwas_annotation=fgwas_annotation[!(fgwas_annotation$Chr=="chrY"|fgwas_annotation$Chr=="chrX"),]  # remove sex chromosomes
write.table(fgwas_annotation,
            "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/peak_calling/IDR_peaks_from_bam_bds_pipeline/conservative_set/fGWAS_merged_narrowpeaks_conservative_set_highQualityBams.bed",
            sep="\t",row.names = F,col.names = F,quote = F)


#### to save fGWAs file with peaks that are present in only one stage (minus those in sex chr, no filter by counts) ###########################
binary=binary[which(rowSums(binary[5:ncol(binary)])==1),] # get peaks present in only one stage

peak_per_stage=list()

for(s in stages){  # get name of stage per peak
  peak_per_stage[[s]]=binary[which(binary[s]==1),c("Chr","Start","End")]
  peak_per_stage[[s]]$Stage=rep(s,nrow(peak_per_stage[[s]]))
}
fgwas_annotation=do.call("rbind",peak_per_stage) # bind all peaks

fgwas_annotation=fgwas_annotation[order(fgwas_annotation[,1],fgwas_annotation[,2]),]
fgwas_annotation=fgwas_annotation[!(fgwas_annotation$Chr=="chrY"|fgwas_annotation$Chr=="chrX"),]

write.table(fgwas_annotation,
            "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/peak_calling/IDR_peaks_from_bam_bds_pipeline/conservative_set/fGWAS_merged_narrowpeaks_uniquepeaks_conservative_set_highQualityBams.bed",
            sep="\t",row.names = F,col.names = F,quote = F)

#################### from diff chromatin accessibility data ####################################

###### to save fGWAS files with peaks assigned to stage with maximum diff chromatin accessibility ###############
# stage-specific (one stage vs all others)

setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/stage-specific/")

# loading Diff open peaks per stage (across-stages method)
BDS_DEA_tables=list()
for(s in stages){
  BDS_DEA_tables[[s]]=read.csv(paste("peak2017-12-04_sig_maxvals_",s,"_diff_peak_analysis_stage-specific.csv",sep=""))
  BDS_DEA_tables[[s]]=BDS_DEA_tables[[s]][,2:ncol(BDS_DEA_tables[[s]])]
  BDS_DEA_tables[[s]]=BDS_DEA_tables[[s]][,c(2:ncol(BDS_DEA_tables[[s]]),1)]
  BDS_DEA_tables[[s]]=BDS_DEA_tables[[s]][,1:3]
  colnames(BDS_DEA_tables[[s]])=c("Chr","Start","End")
  BDS_DEA_tables[[s]]$Stage=rep(s,nrow(BDS_DEA_tables[[s]]))
}

fgwas_annotation=do.call("rbind",BDS_DEA_tables) # bind all peaks

write.table(fgwas_annotation,
            "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/stage-specific/fGWAS_merged_narrowpeaks_diffATAC_stage-specific_conservative_set_highQualityBams.bed",
            sep="\t",row.names = F,col.names = F,quote = F)

###### to save fGWAS files with peaks assigned to stage with maximum diff chromatin accessibility ###################
# across-stages (each stage vs iPSC baseline)


setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/across-stages/")


# loading Diff open peaks per stage (across-stages method)
BDS_DEA_tables=list()
for(s in stages){
  BDS_DEA_tables[[s]]=read.csv(paste("peak2017-12-04_sig_maxvals",s,"_diff_peak_analysis_across-stages.csv",sep=""))
  BDS_DEA_tables[[s]]=BDS_DEA_tables[[s]][,2:ncol(BDS_DEA_tables[[s]])]
  BDS_DEA_tables[[s]]=BDS_DEA_tables[[s]][,c(2:ncol(BDS_DEA_tables[[s]]),1)]
  BDS_DEA_tables[[s]]=BDS_DEA_tables[[s]][,1:3]
  BDS_DEA_tables[[s]]$Stage=rep(s,nrow(BDS_DEA_tables[[s]]))
  
}
fgwas_annotation=do.call("rbind",BDS_DEA_tables) # bind all peaks

write.table(fgwas_annotation,
            "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/across-stages/fGWAS_merged_narrowpeaks_diffATAC_across-stages_conservative_set_highQualityBams.bed",
            sep="\t",row.names = F,col.names = F,quote = F)