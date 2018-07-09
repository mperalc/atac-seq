# plotting longitudinal variation of peak counts


# Load necessary libraries
library(readr)  # to read big tables fast
library(ggplot2)
library(reshape2)  # to modify dataframes for ggplot2
library(ggpubr)

################## variable #################

list = unlist(read.table("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/SNP_overlap_peaks/peaks_over_SNPs_names_for_r.txt"))
list = unique(list)
list = list[80:length(list)]

output_name = "peak_chr17_7788128"
output_type = "tiff"  # can be: "png", "tiff", "multiple_png", "multiple_tiff", "pdf"
output_directory = "/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/peak_profiles/"
ncol = 4
# colour palette
diaPalette <-
  c("#C15858", "#6DA567", "#7883BA")  # Diabetologia palette

origin = c("Ad2.1", "Ad3.1", "Neo1.1")  # original samples, here called origINS

stages = c("iPSC", "DE", "PGT", "PFG", "PE", "EP", "EN", "BLC") # 8 differentiation stages

#stage= old_stages=c("iPSC","DE","GT","PF","PE","EN")
# origin= old_origin=c("Sbad2.1","Sbad3.4")

load(
  "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/session_objects/dge_atac-seq.xz"
)

cpm = as.data.frame(cpm(dge$counts))


QC_and_reshaping = function() {
  
  plot_long = cpm[match(list, rownames(cpm)), ]
  plot_long$Name = rownames(plot_long)
  plot_long = na.omit(plot_long)              # remove NAs, in case there was not a match for a gene
  
  gene_number = nrow(plot_long)   # how many genes to plot
  
  # melt data for ggplot2
  
  long = melt(plot_long, measure.vars = c(1:(ncol(plot_long) - 1)))
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
  long$Name = factor(long$Name, levels = list)
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
    
    facet_wrap( ~ Name, scales = "free", ncol = ncol)
  
  if (output_type == "multiple_png") {
    png(
      paste(
        output_directory,
        output_name,"_ATAC.png",,
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
        output_name,"_ATAC.tiff",
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

# 
# 
# else{
#   if (output_type == "pdf") {
#     long = QC_and_reshaping()  # reshape for plotting
#     head(long)
#     
#     #create plots first
#     plot_list = list()
#     
#     for (i in unique(long$GeneName))
#     {
#       long2 = long[long$GeneName == i, ] # subset each gene
#       p <-
#         ggplot(data = long2, aes(x = stage, y = value, group = group)) +
#         ggtitle(unique(long2$GeneName)) +
#         xlab ("Differentiation stages") +
#         ylab ("Expression [TPM]") +
#         expand_limits(y = 0) +
#         geom_hline(
#           yintercept = 0,
#           linetype = "dashed",
#           size = 1,
#           col = "#DCDCDC"
#         ) +
#         geom_line(aes(linetype = genotype, col = crispr), size = 1) +
#         scale_colour_manual(values = pal) +  # pallete
#         geom_point(size = 3, aes(shape = genotype, col = crispr)) +
#         #scale_colour_manual(values="#000000") +  # for black and white, otherwise map lines and point colours to samples
#         theme_bw() +
#         theme(
#           panel.grid.minor = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.border = element_rect(size = 2),
#           axis.text = element_text(size = 12, face = "bold"),
#           axis.title = element_text(size = 14, face = "bold"),
#           plot.title = element_text(size = 16, face = "bold.italic"),
#           legend.text = element_text(size = 11, face = "bold"),
#           legend.key.height = unit(1, "cm"),
#           legend.title = element_text(size = 13, face = "bold")
#         )
#       plot_list[[i]] = p
#       print(i)
#     }
#     # create progress bar
#     library(tcltk)
#     total = length(unique(long$GeneName))
#     pb <- tkProgressBar(
#       title = "progress bar",
#       min = 0,
#       max = total,
#       width = 300
#     )
#     
#     
#     # I have to open pdf connection before looping so that it saves one gene on each page
#     somePDFPath = paste(output_directory,
#                         output_name,"_CRISPR_RREB1",
#                         ".pdf",
#                         sep = "")
#     pdf(file = somePDFPath)
#     j = 1
#     for (i in unique(long$GeneName)) {
#       print(plot_list[[i]])
#       Sys.sleep(0.1)
#       setTkProgressBar(pb, j, label = paste(round(j / total * 100, 0),
#                                             "% done"))
#       j = j + 1
#     }
#     close(pb)
#     dev.off()
#     
#   } else{
#     for (genes in genes) {
#       # individually gene by gene
#       
#       long = QC_and_reshaping()  # reshape for plotting
#       head(long)
#       p <-
#         ggplot(data = long, aes(x = stage, y = value, group = group)) +
#         ggtitle(unique(long$GeneName)) +
#         xlab ("Differentiation stages") +
#         ylab ("Expression [TPM]") +
#         expand_limits(y = 0) +
#         geom_hline(
#           yintercept = 0,
#           linetype = "dashed",
#           size = 1,
#           col = "#DCDCDC"
#         ) +
#         geom_line(aes(linetype = genotype, col = crispr), size = 1) +
#         scale_colour_manual(values = pal) +  # pallete
#         geom_point(size = 3, aes(shape = genotype, col = crispr)) +
#         #scale_colour_manual(values="#000000") +  # for black and white, otherwise map lines and point colours to samples
#         theme_bw() +
#         theme(
#           panel.grid.minor = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.border = element_rect(size = 2),
#           axis.text = element_text(size = 12, face = "bold"),
#           axis.title = element_text(size = 14, face = "bold"),
#           plot.title = element_text(size = 16, face = "bold.italic"),
#           legend.text = element_text(size = 11, face = "bold"),
#           legend.key.height = unit(1, "cm"),
#           
#           legend.title = element_text(size = 13, face = "bold")
#         )
#       
#       if (output_type == "png") {
#         png(
#           paste(
#             output_directory,
#             output_name,"_CRISPR_RREB1",
#             ".png",
#             sep = ""
#           ),
#           type = "cairo",
#           width = 8,
#           height = 5,
#           units = "in",
#           res = 300,
#           pointsize = 12
#         )
#         print(p)
#         dev.off()
#         
#       }
#       if (output_type == "tiff") {
#         tiff(
#           paste(
#             output_directory,
#             output_name,"_CRISPR_RREB1",
#             ".tiff",
#             sep = ""
#           ),
#           type = "cairo",
#           compression = "lzw",
#           antialias = "default",
#           width = 8,
#           height = 5,
#           units = "in",
#           res = 1000,
#           pointsize = 13
#         )
#         print(p)
#         dev.off()
#       }
#     }
#   }
#   
# for (n in unique(long$Name)) {
#   long2 = long[which(long$Name == n), ]
#   diaPalette <-
#     c("#C15858", "#6DA567", "#7883BA")  # Diabetologia palette
#   p <- ggplot(data = long2, aes(x = stage, y = value, group = Sample)) +
#     ggtitle(unique(long2$Name)) +
#     xlab ("Differentiation stages") +
#     ylab ("Accessibility counts [CPM]") +
#     expand_limits(y = 0) +
#     geom_hline(
#       yintercept = 0,
#       linetype = "dashed",
#       size = 1,
#       col = "#DCDCDC"
#     ) +
#     geom_line(aes(linetype = Sample, col = Sample), size = 1) +
#     scale_colour_manual(values = diaPalette) +  # diabetologia pallete
#     geom_point(size = 3, aes(shape = Sample, col = Sample)) +
#     #scale_colour_manual(values="#000000") +  # for black and white, otherwise map lines and point colours to samples
#     theme_bw() +
#     theme(
#       panel.grid.minor = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.border = element_rect(size = 2),
#       axis.text = element_text(size = 12, face = "bold"),
#       axis.title = element_text(size = 14, face = "bold"),
#       plot.title = element_text(size = 16, face = "bold.italic"),
#       legend.text = element_text(size = 11, face = "bold"),
#       legend.title = element_text(size = 13, face = "bold")
#     )
#   
#   png(
#     paste(
#       "/Users/Marta/Documents/WTCHG/DPhil/Plots/",
#       n,
#       "_ATAC.png",
#       sep = ""
#     ),
#     type = "cairo",
#     width = 8,
#     height = 5,
#     units = "in",
#     res = 300,
#     pointsize = 12
#   )
#   print(p)
#   dev.off()
#   
# }
