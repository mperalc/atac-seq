# plot variation of peak counts accross stages for peak nearby (non-promoter) of selected gene
# also plot RNA-seq data for said gene?
# plotting longitudinal variation of peak counts


# Load necessary libraries
library(readr)  # to read big tables fast
library(ggplot2)
library(reshape2)  # to modify dataframes for ggplot2
library(ggpubr)
library(GenomicRanges)
library(ChIPpeakAnno)
data(TSS.human.GRCh37) # load TSS info from CRCh37
library(biomaRt)
################## variable #################
gene = "SOX4"
output_name = gene
output_type = "multiple_tiff"  # can be:  "multiple_png", "multiple_tiff"
output_directory = "/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/peak_profiles/"
ncol = 4
TenCPM_filter = F  # peaks that have over 10CPM? (T or F) 

## aesthetics
# colour palette
diaPalette <-
  c("#C15858", "#6DA567", "#7883BA")  # Diabetologia palette

origin = c("Ad2.1", "Ad3.1", "Neo1.1")  # original samples, here called origINS

stages = c("iPSC", "DE", "PGT", "PFG", "PE", "EP", "EN", "BLC") # 8 differentiation stages


#####################################
# select dataset to work on
if (TenCPM_filter == F ) {
  peaks = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/2018-07-09annotated_peaks_nearby_non_promoters_allCPMs.txt",
                     header = T)
} else {
  peaks = read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/2018-07-09annotated_peaks_nearby_non_promoters_over10CPMs.txt",
                     header = T)
}

QC_and_reshaping = function() {
  
  plot_long = peaks[match(gene, peaks$external_gene_name), ]
  plot_long = na.omit(plot_long)              # remove NAs, in case there was not a match for a gene
  
  gene_number = nrow(plot_long)   
  
  # melt data for ggplot2
  
  long = melt(plot_long, measure.vars = c(9:32))
  head(long)
  
  # rename stages and samples
  stage_2 = c("iPSC", "DE", "GT", "PF", "PE", "EP", "EN", "BLC")  #shortening EN names EN7= BLC (beta-like cells)
  # stage_2=old_stages
  samples <- c(rep(origin, 8))
  # samples <- c(rep(c("Sample 1","Sample 2","Sample 3"),8)) # for paper
  
  # samples <- c(rep(origin,6))
  long$variable = rep(stage_2, each = 3 * gene_number)                    # sample size times number of genes
  # long$variable=rep(stage_2,each=2*gene_number)
  
  colnames(long)[which(names(long) == "variable")] <- "stage"
  long$Sample = rep(samples, each = gene_number)
  long$stage <- factor(long$stage, levels = stage_2)
  long$Sample = as.factor(long$Sample)
  long$peak = factor(long$peak, levels = unique(long$peak))
  head(long)
  return(long)
  
} 




###################### plot ########################


if (output_type == "multiple_png" | output_type == "multiple_tiff") {
  ################ 1st part as QC script
  long = QC_and_reshaping()
  
  ###################### plot ########################
  
  p <-
    ggplot(data = long, aes(x = stage, y = value, group = Sample)) +
    ylab ("Accessibility counts [CPM]") +
    expand_limits(y = 0) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      size = 1,
      col = "#DCDCDC"
    ) +
    geom_line(aes(linetype = Sample, col = Sample), size = 1) +
    geom_point(size = 3, aes(shape = Sample, col = Sample)) +
    scale_color_manual(values = diaPalette) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_rect(size = 2),
      axis.text.x = element_text(
        size = 16,
        face = "bold",
        angle = 45,
        vjust = 0.55
      ),
      axis.text.y = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 16, face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(size = 16, face = "bold.italic")
    ) +  # strip controls subfigure titles
    
    facet_wrap( ~ peak, scales = "free", ncol = ncol)
  
  if (output_type == "multiple_png") {
    png(
      paste(
        output_directory,
        output_name,"_ATAC_peaks_nearby_non_promoter.png",
        sep = ""
      ),
      type = "cairo",
      antialias = "default",
      width = 7,
      height = 7,
      units = "in",
      res = 600,
      pointsize = 13
    )
    
    print(p)
    dev.off()
  }
  if (output_type == "multiple_tiff") {
    tiff(
      paste(
        output_directory,
        output_name,"_ATAC_peaks_nearby_non_promoter.tiff",
        sep = ""
      ),
      type = "cairo",
      compression = "lzw",
      antialias = "default",
      width = 14,
      height = 18,
      units = "in",
      res = 600,
      pointsize = 13
    )
    
    print(p)
    dev.off()
  }
}
