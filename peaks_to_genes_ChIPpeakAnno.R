# Annotating peaks with ChIPpeakAnno
# from differential atac-seq analysis

  library(ChIPpeakAnno)
  data(TSS.human.GRCh37) # load TSS info from CRCh37
  library(biomaRt)
  library(xlsx) # export as excel workbook
  
  currentDate <- Sys.Date() # to save date in name of output files
  
  jgc <- function()
  {
    .jcall("java/lang/System", method = "gc")
  }   # to manage garbage for xlsx writing
  options(java.parameters = "-Xmx8000m") # giving more memory to java
  
  ###### variable
  date="2017-12-04"  # date for read in files
  stages=c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC")
  
  
  ##### differential open chromatin data
  
  # initializing lists
  stage_specific=list()
  across_stages=list()
  
  for(s in stages){
      stage_specific[[s]]=read.csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/stage-specific/peak",
                                      date,"_sig_maxvals_",s,
                                      "_diff_peak_analysis_stage-specific.csv",sep=""),
                                     header = T,stringsAsFactors = F)
      stage_specific[[s]]=stage_specific[[s]][,c(2:ncol(stage_specific[[s]]))] # dirty fix
      colnames(stage_specific[[s]])=c("Peak_names","Chr","Start","End","Length","log2FC","adj.P.Val")
      across_stages[[s]]=read.csv(file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/across-stages/peak",
                                      date,"_sig_maxvals",s,
                                      "_diff_peak_analysis_across-stages.csv",sep=""),
                                      header=T,stringsAsFactors = F)
      across_stages[[s]]=across_stages[[s]][,c(2:ncol(across_stages[[s]]))] # dirty fix
      colnames(across_stages[[s]])=c("Peak_names","Chr","Start","End","Length","log2FC","adj.P.Val")

    }

  
  ######################### stage-specific
  
  for(s in stages){
    rangedpeak=RangedData(IRanges(start=stage_specific[[s]]$Start,end=stage_specific[[s]]$End), # passing peak tables as ranged data
                          names=stage_specific[[s]]$Peak_names,space=stage_specific[[s]]$Chr)
    # annotatedPeak = annotatePeakInBatch(rangedpeak, AnnotationData=TSS.human.GRCh37) # annotate with hg19. Displays distance to nearest TSS for EVERY PEAK (default: output=nearestLocation)
    annotatedPeak = annotatePeakInBatch(rangedpeak, AnnotationData=TSS.human.GRCh37,output = "shortestDistance")
    # "shortestDistance" will output nearest features to peaks
    # Displays nearest features to TSS
    
    annotatedPeak = as.data.frame(annotatedPeak)
    #write.table(annotatedPeak, outfilename, sep="\t", row.names=F,col.names = T, quote=F)
    
    # plot distribution of distances to TSS
    # hist(annotatedPeak$shortestDistance)
    # hist(annotatedPeak[which(annotatedPeak$shortestDistance<10000),"shortestDistance"])
    # hist(annotatedPeak[which(annotatedPeak$insideFeature=="overlapStart"),"shortestDistance"]) # overlap the start of the gene
    # hist(annotatedPeak[which(annotatedPeak$insideFeature=="upstream"  & annotatedPeak$shortestDistance<1001),"shortestDistance"]) 
    # hist(annotatedPeak[which(annotatedPeak$insideFeature=="downstream" & annotatedPeak$shortestDistance<401),"shortestDistance"]) 
    # hist(annotatedPeak[which(annotatedPeak$insideFeature=="inside" & annotatedPeak$shortestDistance<401),"shortestDistance"]) 
    # 
    # select those with the following insideFeature(s): overlapStart, inside 400bp (or downstream?) & upstream 1000bp
    # combine dfs
    
    peaks_at_promoter=rbind(annotatedPeak[which(annotatedPeak$insideFeature=="overlapStart"),],
                            annotatedPeak[which(annotatedPeak$insideFeature=="upstream"  & annotatedPeak$shortestDistance<2001),],
                            annotatedPeak[which(annotatedPeak$insideFeature=="inside" & annotatedPeak$shortestDistance<501),])
    # search on biomart ensembl gene names
    
    ensembl37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") 
    
    # look for gene names from ids:
    
    all_ensembl_info_ensemblFound37<- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'),
                                            filters = 'ensembl_gene_id', values = peaks_at_promoter$feature, mart = ensembl37)

    # change column name for emsebl gene ids to merge with new df
    colnames(peaks_at_promoter)[7]= "ensembl_gene_id"
    peaks_at_promoter=merge(peaks_at_promoter,all_ensembl_info_ensemblFound37,by="ensembl_gene_id") # merge old and new info
    peaks_at_promoter$Peak_names=paste("peak",peaks_at_promoter$seqnames,peaks_at_promoter$start,sep="_")
    peaks_at_promoter=merge(peaks_at_promoter,stage_specific[[s]],by="Peak_names")
    
    # reorder and rename columns
    peaks_at_promoter=peaks_at_promoter[c("Peak_names","seqnames","start","end","width","log2FC","adj.P.Val","distancetoFeature","shortestDistance",       
                                          "insideFeature","fromOverlappingOrNearest", "external_gene_name","ensembl_gene_id", "gene_biotype",       
                                          "start_position","end_position","feature_strand")]
    colnames(peaks_at_promoter)=c("Peak_names","chr","peak_start","peak_end","peak_width","log2FC","adj.P.Val","distancetoFeature","shortestDistance",       
                                  "insideFeature","fromOverlappingOrNearest", "external_gene_name","ensembl_gene_id", "gene_biotype",       
                                  "start_position","end_position","feature_strand")
    
    write.table(peaks_at_promoter,
                file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/stage-specific/all_biotypes/DEA",
                           currentDate,"_peaks_at_promoters_",s,
                           "_stage-specific.txt",sep = ""),
                sep="\t",col.names = T,quote = F,row.names = F)
    
    gc()   # garbage collector for java
    jgc() # garbage collector for java
    
      write.xlsx(peaks_at_promoter, 
             file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/stage-specific/all_biotypes/DEA",
                        currentDate,"_peaks_at_promoters_stage-specific.xlsx",sep = ""),
             append=TRUE,sheetName=s, row.names=F) # needs garbage collecting so it doesn't crash
      
      ########## get only protein-coding and lincRNAs
      peaks_at_promoter=peaks_at_promoter[which(peaks_at_promoter$gene_biotype %in% c("protein_coding","lincRNA")),]
      
      write.table(peaks_at_promoter,
                  file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/stage-specific/DEA",
                             currentDate,"_peaks_at_promoters_",s,
                             "_stage-specific.txt",sep = ""),
                  sep="\t",col.names = T,quote = F,row.names = F)
      
      gc()   # garbage collector for java
      jgc() # garbage collector for java
      
      write.xlsx(peaks_at_promoter, 
                 file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/stage-specific/DEA",
                            currentDate,"_peaks_at_promoters_stage-specific.xlsx",sep = ""),
                 append=TRUE,sheetName=s, row.names=F) # needs garbage collecting so it doesn't crash
      
      ############ after this, correlate with genes in TADs
  }
              
              
              
  #################### across-stages
  
  for(s in stages){
  rangedpeak=RangedData(IRanges(start=across_stages[[s]]$Start,end=across_stages[[s]]$End), # passing peak tables as ranged data
                        names=across_stages[[s]]$Peak_names,space=across_stages[[s]]$Chr)
  # annotatedPeak = annotatePeakInBatch(rangedpeak, AnnotationData=TSS.human.GRCh37) # annotate with hg19. Displays distance to nearest TSS for EVERY PEAK (default: output=nearestLocation)
  annotatedPeak = annotatePeakInBatch(rangedpeak, AnnotationData=TSS.human.GRCh37,output = "shortestDistance")
  # "shortestDistance" will output nearest features to peaks
  # Displays nearest features to TSS
  
  annotatedPeak = as.data.frame(annotatedPeak)
  #write.table(annotatedPeak, outfilename, sep="\t", row.names=F,col.names = T, quote=F)
  
  # plot distribution of distances to TSS
  # hist(annotatedPeak$shortestDistance)
  # hist(annotatedPeak[which(annotatedPeak$shortestDistance<10000),"shortestDistance"])
  # hist(annotatedPeak[which(annotatedPeak$insideFeature=="overlapStart"),"shortestDistance"]) # overlap the start of the gene
  # hist(annotatedPeak[which(annotatedPeak$insideFeature=="upstream"  & annotatedPeak$shortestDistance<1001),"shortestDistance"]) 
  # hist(annotatedPeak[which(annotatedPeak$insideFeature=="downstream" & annotatedPeak$shortestDistance<401),"shortestDistance"]) 
  # hist(annotatedPeak[which(annotatedPeak$insideFeature=="inside" & annotatedPeak$shortestDistance<401),"shortestDistance"]) 
  # 
  
  
  peaks_at_promoter=rbind(annotatedPeak[which(annotatedPeak$insideFeature=="overlapStart"),],
                          annotatedPeak[which(annotatedPeak$insideFeature=="upstream"  & annotatedPeak$shortestDistance<2001),],
                          annotatedPeak[which(annotatedPeak$insideFeature=="inside" & annotatedPeak$shortestDistance<501),])
  # search on biomart ensembl gene names
  
  ensembl37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") 
  
  # look for gene names from ids:
  
  all_ensembl_info_ensemblFound37<- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'),
                                          filters = 'ensembl_gene_id', values = peaks_at_promoter$feature, mart = ensembl37)
  
  # change column name for emsebl gene ids to merge with new df
  colnames(peaks_at_promoter)[7]= "ensembl_gene_id"
  peaks_at_promoter=merge(peaks_at_promoter,all_ensembl_info_ensemblFound37,by="ensembl_gene_id") # merge old and new info
  peaks_at_promoter$Peak_names=paste("peak",peaks_at_promoter$seqnames,peaks_at_promoter$start,sep="_")
  peaks_at_promoter=merge(peaks_at_promoter,across_stages[[s]],by="Peak_names")
  
  # reorder and rename columns
  peaks_at_promoter=peaks_at_promoter[c("Peak_names","seqnames","start","end","width","log2FC","adj.P.Val","distancetoFeature","shortestDistance",       
                                        "insideFeature","fromOverlappingOrNearest", "external_gene_name","ensembl_gene_id", "gene_biotype",       
                                        "start_position","end_position","feature_strand")]
  colnames(peaks_at_promoter)=c("Peak_names","chr","peak_start","peak_end","peak_width","log2FC","adj.P.Val","distancetoFeature","shortestDistance",       
                                "insideFeature","fromOverlappingOrNearest", "external_gene_name","ensembl_gene_id", "gene_biotype",       
                                "start_position","end_position","feature_strand")
  
  write.table(peaks_at_promoter,
              file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/across-stages/all_biotypes/DEA",
                         currentDate,"_peaks_at_promoters_",s,
                         "_across-stages.txt",sep = ""),
              sep="\t",col.names = T,quote = F,row.names = F)
  
  gc()   # garbage collector for java
  jgc() # garbage collector for java
  
  write.xlsx(peaks_at_promoter, 
             file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/across-stages/all_biotypes/DEA",
                        currentDate,"_peaks_at_promoters_across-stages.xlsx",sep = ""),
             append=TRUE,sheetName=s, row.names=F) # needs garbage collecting so it doesn't crash
  
  ########## get only protein-coding and lincRNAs
  peaks_at_promoter=peaks_at_promoter[which(peaks_at_promoter$gene_biotype %in% c("protein_coding","lincRNA")),]
  
  write.table(peaks_at_promoter,
              file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/across-stages/DEA",
                         currentDate,"_peaks_at_promoters_",s,
                         "_across-stages.txt",sep = ""),
              sep="\t",col.names = T,quote = F,row.names = F)
  
  gc()   # garbage collector for java
  jgc() # garbage collector for java
  
  write.xlsx(peaks_at_promoter, 
             file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/across-stages/DEA",
                        currentDate,"_peaks_at_promoters_across-stages.xlsx",sep = ""),
             append=TRUE,sheetName=s, row.names=F) # needs garbage collecting so it doesn't crash
  
  ############ after this, correlate with genes in TADs
  }
  
  # plot diff expr numbers per stage vs diff peaks
  # plot long counts variation for peaks (cpm or voom normalized?)
  # % peaks at promoters assigned to same stage, or to stage before, for both methods
  