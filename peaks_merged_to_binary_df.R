# transforming merged narrowpeaks file into presence/absence matrix per peak position (or feature)

# Transforms the counts into a presence/absence matrix
# The input is the merged narrowpeak file generated with bedtools from the combined narrowpeaks
# It must have the following columns: "Chr","Start","End","Peaks_present", "Name"
# AND only those!!
# "Peaks_present" is the list of samples that share the peak "Name"
# It is a long character vector that is splitted inside the function
# The output is a bed file with feature name, presence/absence per sample & info of start and end of peaks 
# (for all peaks present in at least one sample)



peaks_merged_to_binary_df = function(peaks){
  if(ncol(peaks)>5){
    stop("Data frame with wrong number of columns")
  }
  # to get sample names for columns
  samples=paste(peaks$Peaks_present,collapse=",") 
  samples=unique(unlist(strsplit(samples, ",")))
  
  df=data.frame(matrix(ncol = length(samples)+4,nrow = nrow(peaks))) # make matrix to fill with n columns for samples 
                                                  # + 4 more (Name of feature, chr, start, end)
  colnames(df)=c(c("Name","Chr","Start","End"),samples) # filling in new df
  df$Name=peaks$Name
  df$Chr=peaks$Chr
  df$Start=peaks$Start
  df$End=peaks$End
  df[,-1:-4] = 0  # all count columns to 0 (absence of peak)
  
  for (c in colnames(df[,-1:-4])){
    
    rows=peaks[grep(c,peaks$Peaks_present),"Name"] # peaks present in sample "c"
    df[which(df$Name %in% rows),which(colnames(df)==c)] = 1 # put 1 in column for "c" if peak is present in "c"
  }
  

  return(df)
}