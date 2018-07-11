# 3.5
# combine motif files in single .motif file
# e.g. to find location of specific motifs in my data
# 
# directory = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/WGCNA/10CPM_P12_S120_deepSplit2_signedHybrid_noOutliers/HOMER/HOMER_output/blue/homerResults"
# 
# #####################
# 
# setwd(directory)
# 
# file_list <- list.files()
# # grep files that contain .motif
# file_list = file_list[grep("\\b.motif\\b", file_list)]
# 
# for (file in file_list) {
#   # if the merged dataset doesn't exist, create it
#   if (!exists("dataset")) {
#     dataset <- readLines(file)
#   }
#   
#   # if the merged dataset does exist, append to it
#   if (exists("dataset")) {
#     temp_dataset <- read.table(file, header = F, fill=TRUE)
#     dataset <- rbind(dataset, temp_dataset)
#     rm(temp_dataset)
#   }
#   
# }

# much easier in bash