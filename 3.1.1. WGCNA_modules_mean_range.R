# WGCNA modules stats


# locations
ATAC_dir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/WGCNA/WGCNA_P12_S120_deepSplit2_signedHybrid_rmOutliersT/"

# Read in module files for RNA and ATAC-seq
ATAC = read.table(file = paste(ATAC_dir, "WGCNA.gene2module.30M.txt", sep = ""),
                  header = T)

peaks_per_module = list()
for (m in unique(ATAC$module)){
  peaks_per_module[[m]] = ATAC[which(ATAC$module == m),"gene"]
}

 module_peaknumber = lapply(peaks_per_module,length)
 module_peaknumber = as.data.frame(do.call("rbind",module_peaknumber))
 max(module_peaknumber$V1) - min(module_peaknumber$V1) # range
 median(module_peaknumber$V1)# median
 