# correlate WGCNA modules

corWGCNAmodules = function(modulesDataSet1,modulesDataSet2, dataSet1, dataSet2) {
  # dataSet1 and 2 contain genomic data in dge format
  # modulesDataSet1,modulesDataSet2 contain WGCNA gene2module results
  
  # Number of ATAC-seq peaks, genes and samples
  nPeaks = nrow(modulesDataSet1)
  nGenes = nrow(modulesDataSet2)
  nSamples = 24
  
  datATAC = t(dataSet1$counts) # ATAC data needs to be in format row-sample and col-peak
  datExpr = t(dataSet2$counts)
  # The eigengene is the summary profile of each module we are going to correlate
  
  MEs0 = moduleEigengenes(datATAC, modulesDataSet1$module)$eigengenes
  MEs = orderMEs(MEs0)
  
  #same for expression data
  MEs0_Expr = moduleEigengenes(datExpr, modulesDataSet2$module)$eigengenes
  MEs_Expr = orderMEs(MEs0_Expr)
  
  # correlation of eigenvectors of RNA and ATAC-seq results of WGCNA
  # columns are ATAC modules
  moduleTraitCor = cor(MEs_Expr, MEs, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  # df to examine more easily
  corDF = as.data.frame(moduleTraitCor)
  
  ########### graphical representation of the correlation
  sizeGrWindow(10, 6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2),
                     "\n(",
                     signif(moduleTraitPvalue, 1),
                     ")",
                     sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  
  par(mar = c(6, 8.5, 3, 3))
  
  tiff(
    file = "Heatmap_cor_modules_WGCNA.tiff",
    width = 25,
    height = 25,
    units = "in",
    compression = "lzw",
    res = 300
  )
  # Display the correlation values within a heatmap plot
  labeledHeatmap(
    Matrix = moduleTraitCor,
    xLabels = names(MEs),
    yLabels = names(MEs_Expr),
    ySymbols = names(MEs_Expr),
    colorLabels = FALSE,
    colors = greenWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.5,
    zlim = c(-1, 1),
    main = paste("Module-trait relationships")
  )
  
  dev.off()
}