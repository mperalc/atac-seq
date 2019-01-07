# WGCNA ATAC-seq only
# ATAC = raw counts filytered 10 CPM

Power = 10
Size = 120
deepsplit = 2
method = "signed hybrid"
remove_outliers = T  # remove outliers? T / F

## load libraries

library(edgeR)
library("DESeq2")
library(RDAVIDWebService)
library(ggplot2)
library("WGCNA")
allowWGCNAThreads()

# directories
 input = "../data/"
 output = "../results/"
#input = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/session_objects/"
#output = "/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/WGCNA/10CPM_P12_S120_deepSplit2_signedHybrid_noOutliers/"

## load data
load(file = paste(input, "dge_atac-seq_10CPM_trim.xz", sep = ""),
     verbose = TRUE)

if (remove_outliers) {
  keep = c(
    "iPSC-sbad2.1",
    "iPSC-sbad3.1",
    "iPSC-neo1.1",
    "DE-sbad2.1",
    "DE-sbad3.1",
    "DE-neo1.1",
    "GT-sbad2.1",
    "GT-sbad3.1",
    "GT-neo1.1",
    "PF-sbad2.1",
    "PF-sbad3.1",
    "PF-neo1.1",
    "PE-sbad2.1",
    "PE-sbad3.1",
    "PE-neo1.1",
    "EP-sbad3.1",
    "EP-neo1.1",
    "EN-sbad3.1",
    "EN-neo1.1",
    "BLC-sbad3.1",
    "BLC-neo1.1"
  )
  dge <-
    dge[, keep, keep.lib.sizes = TRUE]
}

dim(dge)
#36,917    24

## renaming
# ATAC
counts_ATAC = dge$counts


rownames(counts_ATAC)[which(duplicated(rownames(counts_ATAC)))]
# no duplicates


subject = sapply(strsplit(colnames(counts_ATAC), split = "-"), function(x)
  x[2]) # takes second element of column names after splitting by "_"
design = data.frame(
  row.names = colnames(counts_ATAC),
  stage = dge$samples$group,
  subject = subject
) # design is the same for both
dds <- DESeqDataSetFromMatrix(
  countData = counts_ATAC,
  colData = design,
  design = ~ subject + stage
)

design$stage = factor(design$stage, levels(design$stage))

############ Run Vst normalization on the whole dataset (with DESeq2) ############

vsd <- varianceStabilizingTransformation(dds)
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

vstMat = assay(vsd)
write.table(
  vstMat,
  file = paste(output, "vstMat.txt", sep = ""),
  sep = "\t",
  quote = F
)

#### blind=TRUE should be used for comparing samples in an manner unbiased by prior  ###
#### information on samples, for example to perform sample QA (quality assurance).   ###
#### blind=FALSE should be used for transforming data for downstream analysis, where ###
#### the full use of the design information should be made							 ###

# Fix this to do the filtering after the transformation?

vsd_b <- varianceStabilizingTransformation(dds, blind = F)
vstMat_b = assay(vsd)
write.table(
  vstMat_b,
  file = paste(output, "vstMat.blind_F.txt", sep = ""),
  sep = "\t",
  quote = F
)

### MDS plot
pdf(paste(output, "mds.pdf", sep = ""), width = 10)
par(mfrow = c(1, 2))
plotMDS(vstMat, col = as.numeric(dds$stage))
plotMDS(vstMat_b, col = as.numeric(dds$stage))
dev.off()
## no difference in vst transform here


##################### run WGCNA  #####################
#if (!file.exists("net.xz")) {
  degData = t(vstMat)
  gsg = goodSamplesGenes(degData, verbose = 3)
  
  gsg$allOK
  # TRUE
  
  ### pick power
  powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
  sft = pickSoftThreshold(degData, powerVector = powers, verbose = 5)
  pdf(
    file = paste(output, "WGCNA_softThreshold.pdf", sep = ""),
    width = 9,
    height = 5
  )
  
  par(mfrow = c(1, 2))
  
  cex1 = 0.9
  
  plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit,signed R^2",
    type = "n",
    main = paste("Scale independence")
  )
  
  text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  
  # this line corresponds to using an R^2 cut-off of h
  abline(h = 0.90, col = "red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    type = "n",
    main = paste("Mean connectivity")
  )
  text(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  dev.off()
  
  # looks like power of 10 is fine
  ##### run WGCNA module detection
  # the blocks need to be in the same place where the script is run
  
  net = blockwiseModules(
    degData,
    power = Power,
   # maxBlockSize = 15500,
    # 1 block only, takes longer but is more precise
    maxBlockSize = 40000,
    TOMType = "signed",
    networkType = method,
    minModuleSize = Size,
    mergeCutHeight = 0.15,
    deepSplit = deepsplit,
    numericLabels = TRUE,
    saveTOMs = T,
    saveTOMFileBase = paste(output, "blockwiseTOM",sep = ""),
    #loadTOM = T,
    verbose = 3
  )
  
  save(net,
       file = paste(output, "net.xz", sep = ""),
       compress = "xz")
  table(net$colors)
  
  
  pdf(
    paste(output, "WGCNA_dendrogram_block1.30M.pdf", sep = ""),
    width = 12,
    height = 9
  )
  moduleColors = labels2colors(net$colors)
  # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(
    net$dendrograms[[1]],
    moduleColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
  dev.off()
  
  # 
  # pdf(
  #   paste(output, "WGCNA_dendrogram_block2.30M.pdf", sep = ""),
  #   width = 12,
  #   height = 9
  # )
  # moduleColors = labels2colors(net$colors)
  # # Plot the dendrogram and the module colors underneath
  # plotDendroAndColors(
  #   net$dendrograms[[2]],
  #   moduleColors[net$blockGenes[[2]]],
  #   "Module colors",
  #   dendroLabels = FALSE,
  #   hang = 0.03,
  #   addGuide = TRUE,
  #   guideHang = 0.05
  # )
  # dev.off()
  # 
  # 
  # pdf(
  #   paste(output, "WGCNA_dendrogram_block3.30M.pdf", sep = ""),
  #   width = 12,
  #   height = 9
  # )
  # moduleColors = labels2colors(net$colors)
  # # Plot the dendrogram and the module colors underneath
  # plotDendroAndColors(
  #   net$dendrograms[[3]],
  #   moduleColors[net$blockGenes[[3]]],
  #   "Module colors",
  #   dendroLabels = FALSE,
  #   hang = 0.03,
  #   addGuide = TRUE,
  #   guideHang = 0.05
  # )
  # dev.off()
  # 
  # 
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs
  

  
  save(
    MEs,
    moduleLabels,
    moduleColors,
    geneTree,
    file = paste(output, "ATAC_only_network.RData", sep = "")
  )
  
  # Plotting dendrogram of eigenvectors
  
  datME = moduleEigengenes(degData, moduleColors)$eigengenes
  ###### This is Agata's, but it is not used to calculate the plot below
  #dissimME = (1 - t(WGCNA::cor(datME, method = "p"))) / 2
  #hclustdatME = hclust(as.dist(dissimME), method = "average")
  #########################
  
  #pheno = data.frame(stage = as.numeric(design$stage))
  
  #MET = orderMEs(cbind(datME, pheno))
  pdf(
    paste(output, "WGCNA_cluster_eigenvectors.30M.pdf", sep = ""),
    width = 8,
    height = 12
  )
  par(cex = 0.9)
  plotEigengeneNetworks(
    orderMEs(datME),
    #MET,
    "",
    marDendro = c(0, 4, 1, 2),
    marHeatmap = c(3, 4, 1, 2),
    cex.lab = 0.8,
    xLabelsAngle = 90
  )
 
  dev.off()

  # Plot eigengenes with dissimilarity threshold line
  # If I choose Agata's method of calculating dissimilarity I get same separation of modules 
  # but different height
  #dissimME = (1 - t(WGCNA::cor(datME, method = "p"))) / 2
  #hclustdatME = hclust(as.dist(dissimME), method = "average")
  
  dissimME = (1 - WGCNA::cor(datME, method = "p")) # This is the dissimilarity value used in the 
  # function plotEigengeneNetworks
  
  hclustdatME = hclust(as.dist(dissimME), method = "average")
  
  pdf(
    paste(output, "WGCNA_cluster_eigenvectors.dissimline.pdf", sep = ""),
    width = 8,
    height = 12
  )
  par(mfrow = c(1,1))
  plot(hclustdatME,main = "Clustering of module eigengenes",xlab = "", sub = "")
  
  # Plot chosen dissimilarity value
  
  MEDissThres = 0.07
  abline(h = MEDissThres,col = "red")
  dev.off()
  
  ## Plotting barplots
  
  pdf(
    paste(output, "Module_eigengenes.30M.barplots.pdf", sep = ""),
    width = 10,
    height = 3
  )
  for (i in 1:length(unique(moduleColors))) {
    which.module = unique(moduleColors)[i]
    ME = datME[, paste("ME", which.module, sep = "")]
    barplot(
      ME,
      col = as.numeric(dds$stage),
      ylab = "eigengene expression",
      xlab = "array sample",
      main = which.module
    )
  }
  dev.off()
  
  gene2module = data.frame(gene = colnames(degData), module = moduleColors)
  write.table(
    gene2module,
    file = paste(output, "WGCNA.gene2module.30M.txt", sep = ""),
    sep = "\t",
    row.names = F,
    quote = F
  )
  
  #### make ribbon plots for each module  #####
  pdf(paste(output, "Modules.ribbon_plots.pdf", sep = ""))
  for (j in 1:dim(datME)[2]) {
    stages = levels(design$stage)
    min = rep(0, length(stages))
    max = rep(0, length(stages))
    mean = rep(0, length(stages))
    for (i in 1:length(stages)) {
      min[i] = min(datME[which(design$stage == stages[i]), j])
      max[i] = max(datME[which(design$stage == stages[i]), j])
      mean[i] = mean(datME[which(design$stage == stages[i]), j])
    }
    plot_data = data.frame(
      stages = c(1:length(stages)),
      min = min,
      max = max,
      mean = mean
    )
    
    p = ggplot(plot_data, aes(stages, mean)) + geom_ribbon(aes(ymin = min, ymax =
                                                                 max),
                                                           colour = "lightgrey",
                                                           fill = "lightgrey") + geom_line(color = "steelblue4", lwd =
                                                                                             1) + theme_bw() + scale_x_continuous(breaks = c(1:8)) + ylab("module eigengene") +
      ggtitle(colnames(datME)[j])
    
    p = ggplot(plot_data, aes(stages, mean)) + geom_ribbon(aes(ymin = min, ymax =
                                                                 max),
                                                           colour = "lightgrey",
                                                           fill = "lightgrey") + geom_line(color = "steelblue4", lwd =
                                                                                             1) + theme_bw() + scale_x_continuous(breaks = c(1:8), labels = stages) +
      ylab("module eigengene") + ggtitle(colnames(datME)[j])
    
    print(p)
  }
  dev.off()
  
  tiff(
    paste(
      output,
      "module.ribbon_plots_all.tiff",
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
  

#     for (j in 1:dim(datME)[2]) {
#     stages = levels(design$stage)
#     min = rep(0, length(stages))
#     max = rep(0, length(stages))
#     mean = rep(0, length(stages))
#     for (i in 1:length(stages)) {
#       min[i] = min(datME[which(design$stage == stages[i]), j])
#       max[i] = max(datME[which(design$stage == stages[i]), j])
#       mean[i] = mean(datME[which(design$stage == stages[i]), j])
#     }
#     plot_data = data.frame(
#       stages = c(1:length(stages)),
#       min = min,
#       max = max,
#       mean = mean
#     )
#     
#     p = ggplot(plot_data, aes(stages, mean)) + geom_ribbon(aes(ymin = min, ymax =
#                                                                  max),
#                                                            colour = "lightgrey",
#                                                            fill = "lightgrey") + geom_line(color = "steelblue4", lwd =
#                                                                                              1) + theme_bw() + scale_x_continuous(breaks = c(1:8)) + ylab("module eigengene") +
#       ggtitle(colnames(datME)[j])
#     
#     p = ggplot(plot_data, aes(stages, mean)) + geom_ribbon(aes(ymin = min, ymax =
#                                                                  max),
#                                                            colour = "lightgrey",
#                                                            fill = "lightgrey") + geom_line(color = "steelblue4", lwd =
#                                                                                              1) + theme_bw() + scale_x_continuous(breaks = c(1:8), labels = stages) +
#       ylab("module eigengene") + ggtitle(colnames(datME)[j]) + 
#       facet_wrap( ~ colnames(datME), scales = "free", ncol = ncol)
#     
#     
#     print(p)
#   dev.off()
# }

if (file.exists(paste(output, "vstMat.txt", sep = ""))) {
  vstMat = read.table(paste(output, "vstMat.txt", sep = ""))
  degData = t(vstMat)
  load(paste(output, "net.xz", sep = ""))
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  ### get connectivity values for each gene within its module
  ### get connectivity values:
  ADJ1 = abs(WGCNA::cor(degData, use = "p")) ^ Power
  Alldegrees1 = intramodularConnectivity(ADJ1, moduleColors)
  Alldegrees1$Module = moduleColors
  
  write.table(
    Alldegrees1,
    file = paste(output, "WGCNA.connectivity.gene2module.txt", sep = ""),
    sep = "\t",
    quote = F
  )
  
  # split table based on modules
  
  Alldegrees1_list <- split(Alldegrees1 , f = Alldegrees1$Module)
  for (i in 1:length(Alldegrees1_list)) {
    write.table(
      Alldegrees1_list[i],
      file = paste0(
        output,
        "WGCNA.connectivity.gene2module.",
        names(Alldegrees1_list)[i],
        ".txt"
      ),
      sep = "\t",
      quote = F
    )
  }
  
}