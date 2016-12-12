# Filtering of atac-seq data and differential "expression"  analysis of peaks


# with edgeR?
# with limma+voom

# 1st: trimming & normalising

#3rd DE

library(edgeR)
library(knitr)   # for kable tables
library(ggplot2)

currentDate <- Sys.Date() # to save date in name of output files

load(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/session_objects/dge_atac-seq.xz",verbose=TRUE)  #loading the dge object for conservative counts


sample=factor(c("sbad2.1","sbad3.1","neo1.1"),levels=c("sbad2.1","sbad3.1","neo1.1"))   # data comes from three donors

stage=factor(c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC"),levels=c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC")) # 8 stages


dge <- calcNormFactors(dge) # Calculating normalization factors

# plot PCA before filtering  ######################

samples <- rep(sample,8)   

stages <-rep(stage,each=3) 


##plot PCA for un/normalized counts
# par(mfrow=c(1,1))
# 
# # looks better with scaling scaling
# plot_pca=function(x,s=samples,st=stages){
#   pca1<-prcomp(t(x), retx=TRUE, scale. = F)
# 
#   plot(pca1, type = "l") #variance vs first 10 components
#   summary(pca1)          #importance of each component (important line is "proportion of variance")
# 
#   percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
#   percentVar <- round(100 * percentVar)
#   pcs <- as.data.frame(pca1$x)
#   pcs <- cbind(pcs,sample=samples,stage=stages)
#   pcs$stage <- ordered(pcs$stage, levels = stage)
# 
#   p <- ggplot(pcs, aes(PC1, PC2, colour=stage, shape=samples)) +
#     geom_point(size = 3) + xlab (paste0( "PC1:" ,percentVar[ 1 ],"% variance")) +
#     ylab (paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
#   p <- p + theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                     panel.background = element_blank(),
#                     panel.border = element_rect(fill = NA, colour = "black"),
#                     legend.key = element_blank(),# legend.position = c(0.5,0.5),
#                     axis.title.y = element_text(face="bold", angle=90, size=12, vjust=0.2),
#                     axis.title.x = element_text(face="bold", size=12, vjust=0),
#                     axis.text.x = element_text(face="bold", colour = "black", angle=90, size=12, vjust=0.2, hjust =1 ),
#                     axis.text.y = element_text(face="bold", colour = "black"),
#                     axis.ticks = element_line(colour = "black"),
#                     axis.line = element_line(colour = "black"))
#   # ggsave("/Users/Marta/Documents/WTCHG/DPhil/Plots/conservative_counts/Diff_v2_PCA.tiff",p,compression="lzw")
#   return(p)
# }
# 
# p=plot_pca(dge$counts)
# p
# ggsave(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/PCA_unfiltered_counts_atac-seq",currentDate,".jpg",sep=""),p,width=6,height=5,units="in",dpi=300)
# 
# CPMs=cpm(dge, normalized.lib.sizes=TRUE, log=FALSE)  #calculate counts per million using the normalized library sizes
# 
# CPMs=as.data.frame(CPMs)
# p=plot_pca(CPMs)
# p
# ggsave(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/PCA_unfiltered_CPMs_atac-seq",currentDate,".jpg",sep=""),p,width=6,height=5,units="in",dpi=300)




###############

old_peaks=dim(dge)[1]# save how many peaks there are initially

trim_by_stage_and_cpm=function(dge,CPM=1,old_peaks=dim(dge)[1]){
  keep <- rowSums(cpm(dge,normalized.lib.sizes = TRUE)>CPM) >= 3   # getting peaks with at least n cpm in at least one stage (3 rep per stage).
  dge <- dge[keep, , keep.lib.sizes=FALSE]  # adjusting the library size
  
  dge <- calcNormFactors(dge)# recalculating normalization factors after altering library sizes
  # if counts are distributed roughly equally across all samples, it shouldn0t change much
  
  # go through the table 3 by 3
  # if there are more than 2 values in three columns above 0, go to next row.
  # if not, go to next 3
  # if gets to end of loop and not found the values, get that row. After loop, subset matrix excluding those rows
  
  cpm_data=cpm(dge, normalized.lib.sizes=TRUE)
  
  trim_by_stage=function(dge){
    data=cpm(dge, normalized.lib.sizes=TRUE)
    data[data<CPM]=0
    data[data>CPM]=1 # convert to table where 0=0 and 1 equals any value above n
    
    c=c() # initiate vector of numbers, to save rows to take out of matrix
    
    for (r in 1:nrow(data)){
      for (i in seq(1,24,3)){   # loop through stages 3 by 3 (to jump over samples of same stage)
        
        if (sum(data[r,c(i,i+1,i+2)]) > 2) {
          
          break
        }
        if (i==22){
          
          c=c(c,r)  # if it gets to the end and there are no valid values, append to vector
        }
        
        
      }
      
    }
    
    rm(data)
    return(c)
    
  }
  
  c=trim_by_stage(dge)  # calling previous function
  
  fail=100*(length(c)/nrow(dge$counts))
  genes_fail=length(c)
  ditch=dge$counts[c,] # table to throw away
  
  keep=!(rownames(dge$counts) %in% rownames(ditch))
  
  dge <- dge[keep, , keep.lib.sizes=FALSE] # counts now has 17892 genes
  
  rm(c,keep,ditch)
  
  
  
  return(dge)
}

# call function for CPM=1

dge <- trim_by_stage_and_cpm(dge,CPM=1,old_peaks=dim(dge)[1])
dge <- calcNormFactors(dge)# recalculating normalization factors after altering library sizes

nrow(dge)/old_peaks #proportion of peaks that remain

# save(dge, file = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/session_objects/dge_atac-seq_trimmed.xz" , compress="xz")   # saving the dge object
# 
# # plot PCA for un/normalized counts, after filtering
# par(mfrow=c(1,1))
# 
# # looks better with scaling scaling
# plot_pca=function(x,s=samples,st=stages){
#   pca1<-prcomp(t(x), retx=TRUE, scale. = TRUE)
# 
#   plot(pca1, type = "l") #variance vs first 10 components
#   summary(pca1)          #importance of each component (important line is "proportion of variance")
# 
#   percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
#   percentVar <- round(100 * percentVar)
#   pcs <- as.data.frame(pca1$x)
#   pcs <- cbind(pcs,sample=samples,stage=stages)
#   pcs$stage <- ordered(pcs$stage, levels = stage)
# 
#   p <- ggplot(pcs, aes(PC1, PC2, colour=stage, shape=samples)) +
#     geom_point(size = 3) + xlab (paste0( "PC1:" ,percentVar[ 1 ],"% variance")) +
#     ylab (paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
#   p <- p + theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                     panel.background = element_blank(),
#                     panel.border = element_rect(fill = NA, colour = "black"),
#                     legend.key = element_blank(),# legend.position = c(0.5,0.5),
#                     axis.title.y = element_text(face="bold", angle=90, size=12, vjust=0.2),
#                     axis.title.x = element_text(face="bold", size=12, vjust=0),
#                     axis.text.x = element_text(face="bold", colour = "black", angle=90, size=12, vjust=0.2, hjust =1 ),
#                     axis.text.y = element_text(face="bold", colour = "black"),
#                     axis.ticks = element_line(colour = "black"),
#                     axis.line = element_line(colour = "black"))
#   # ggsave("/Users/Marta/Documents/WTCHG/DPhil/Plots/conservative_counts/Diff_v2_PCA.tiff",p,compression="lzw")
#   return(p)
# }
# 
# p=plot_pca(dge$counts)
# p
# ggsave(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/PCA_filtered_counts_atac-seq_withscaling",currentDate,".jpg",sep=""),p,width=6,height=5,units="in",dpi=300)
# 
# CPMs=cpm(dge, normalized.lib.sizes=TRUE, log=FALSE)  #calculate counts per million using the normalized library sizes
# 
# CPMs=as.data.frame(CPMs)
# p=plot_pca(CPMs)
# p
# ggsave(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/PCA_filtered_CPMs_atac-seq_withscaling",currentDate,".jpg",sep=""),p,width=6,height=5,units="in",dpi=300)


#################################### DEA  

library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)
library(limma)
library(ggbiplot)
library(ggrepel) # provides geoms for ggplot2 to repel overlapping text labels.

library(ggfortify)
library(reshape2)  # to modify dataframes for ggplot2
library(knitr)   #for tables in rmd
library(pander)  #more tables


################################################ peak stages ###################################
#create the design matrix
design1 <- model.matrix(~stages + samples)                    
kable(design1)


# This converts the counts to log-counts per million with associated precision weights. After this, the RNA-seq data can be analyzed as if it was microarray data.


v2=voom(dge,design=design1,plot=TRUE)   


par(mfrow=c(1,2))
plot(log2(dge$counts + 1)[,1:2],
     pch=16, cex=0.3, main="log2")
plot(v2$E[,1:2],
     pch=16, cex=0.3, main=" voom norm counts")

# For first 50 genes, see the distribution of counts:
random_sample=sample(nrow(dge),50)
par(mfrow=c(1,3))
plot(density(as.numeric(unlist(dge$counts[random_sample,]))), main="counts", cex.main=2)
plot(density(as.numeric(unlist(log2(dge$counts+1)[random_sample,]))), main="log2", cex.main=2)
plot(density(as.numeric(unlist(v2$E[random_sample,]))), main="voom norm counts", cex.main=2)


#function wrap so that plotMDS doesn't produce a plot when called:
plotMDS.invisible <- function(...){
  ff <- tempfile()
  png(filename=ff)
  res <- plotMDS(...)
  dev.off()
  unlink(ff)
  res
}

pretty_mds=function(v2){
  mds_p=plotMDS.invisible(v2$E,gene.selection = "pairwise")    # pairwise method (default)
  mds_c=plotMDS.invisible(v2$E,gene.selection = "common")      #common method 
  
  
  # Rearrange data for ggplot
  
  # method: pairwise
  m_p=as.data.frame(mds_p$cmdscale.out)
  m_p <- cbind(m_p,sample=samples,stage=stages)
  colnames(m_p)=c("Dimension 1", "Dimension 2", "Samples","Stages")
  
  # method: common
  m_c=as.data.frame(mds_c$cmdscale.out)
  m_c <- cbind(m_c,sample=samples,stage=stages)
  colnames(m_c)=c("Dimension 1", "Dimension 2", "Samples","Stages")
  
  
  # plot pairwise
  
  # if(qc_plots==TRUE){
  #tiff("/Users/Marta/Documents/WTCHG/DPhil/Plots/conservative_counts/Diff_v2_MDS_pairwise_only_ours.tiff", type="cairo",
  # width=10,height=10,units="in",res=300,pointsize = 13,compression="lzw")
  mp=ggplot(m_p) +
    geom_point(aes(`Dimension 1` ,`Dimension 2`), size = 2, color = 'grey') +
    geom_label_repel(
      aes(`Dimension 1` , `Dimension 2`, fill = factor(Samples), label = Stages),
      fontface = 'bold', color = 'white',
      box.padding = unit(0.25, "lines"),
      point.padding = unit(0.25, "lines")
    ) +
    coord_fixed(ratio = 1.2) +
    
    theme_bw(base_size=20) + 
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank() )+
    theme(panel.border= element_rect())+
    
    theme(legend.position = "bottom")  +
    ggtitle("MDS plot:pairwise method")+
    labs(fill = "Sample")
  # dev.off()
  
  # plot common
  # tiff("/Users/Marta/Documents/WTCHG/DPhil/Plots/conservative_counts/Diff_v2_MDS_common_only_ours.tiff", type="cairo",
  #      width=10,height=10,units="in",res=300,pointsize = 13,compression="lzw")
  mc= ggplot(m_c) +
    geom_point(aes(`Dimension 1` ,`Dimension 2`), size = 2, color = 'grey') +
    geom_label_repel(
      aes(`Dimension 1` , `Dimension 2`, fill = factor(Samples), label = Stages),
      fontface = 'bold', color = 'white',
      box.padding = unit(0.25, "lines"),
      point.padding = unit(0.25, "lines")
    ) +
    coord_fixed(ratio = 2) +  #fix x-y ratio
    
    theme_bw(base_size=20) + 
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank() )+
    theme(panel.border= element_rect())+
    
    theme(legend.position = "bottom")  +
    ggtitle("MDS plot: common method")+
    labs(fill = "Sample")
  # dev.off()
  
  #  }
  return(list(mp,mc))   #return both plots as a list
}

mds_plots=pretty_mds(v2)

plot(mds_plots[[1]])
plot(mds_plots[[2]])


# get sample distance cluster  #


plot_sdc=function(x){
  sampleDists <- dist(t(x))   
  #This function computes and returns the distance matrix computed by using the specified distance measure to compute the distances between the rows of a data matrix.
  # by default, "euclidean"
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(stages, samples, sep=" - ")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  #png("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/distance_matrix_voom_peaks.png", type="cairo",
  #    width=7,height=5,units="in",res=200,pointsize = 13)
  
  sdc= pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
  
  # dev.off()
  
  return(sdc)
}



sdc_voom=plot_sdc(v2$E) # I used the design v2, but it gives the same result with any other (depends on the expression values, v$E doesn't vary with design matrix)


print(sdc_voom)



# tried for different designs: exactly the same
par(mfrow=c(1,1))
plot_pca=function(x,s=samples,st=stages){
  pca1<-prcomp(t(x), retx=TRUE,scale. = T)
  
  plot(pca1, type = "l") #variance vs first 10 components
  summary(pca1)          #importance of each component (important line is "proportion of variance")
  
  percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs,sample=samples,stage=stages)
  pcs$stage <- ordered(pcs$stage, levels = stage)
  
  p <- ggplot(pcs, aes(PC1, PC2, colour=stage, shape=samples)) + 
    geom_point(size = 3) + xlab (paste0( "PC1:" ,percentVar[ 1 ],"% variance")) + 
    ylab (paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
  p <- p + theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    panel.border = element_rect(fill = NA, colour = "black"), 
                    legend.key = element_blank(),# legend.position = c(0.5,0.5),
                    axis.title.y = element_text(face="bold", angle=90, size=12, vjust=0.2),
                    axis.title.x = element_text(face="bold", size=12, vjust=0),
                    axis.text.x = element_text(face="bold", colour = "black", angle=90, size=12, vjust=0.2, hjust =1 ),
                    axis.text.y = element_text(face="bold", colour = "black"),
                    axis.ticks = element_line(colour = "black"),
                    axis.line = element_line(colour = "black"))
  # ggsave("/Users/Marta/Documents/WTCHG/DPhil/Plots/conservative_counts/Diff_v2_PCA.tiff",p,compression="lzw")
  return(p)
}

p=plot_pca(v2$E)
p

ggsave(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/PCA_voom_atac-seq_withscaling_",currentDate,".jpg",sep=""),p,width=6,height=5,units="in",dpi=300)



###DEA per se


d=as.data.frame(design1)

stagesiPSC=as.factor(c(1,1,1, rep(0,21)))  # add iPSC to design matrix, because the function to make it ignores this (is the baseline)
d=cbind(d,stagesiPSC)
d=d[c(1,11,2:10)] # re-order to loop through it without testing intercept or samples

#initialize list of data frames to save them inside the loop

DE_list<-list()



for (i in colnames(d[2:9])) {
  
  test=as.factor(d[,which(colnames(d)==i)]) 
  design2=model.matrix(~samples+test)
  v2=voom(dge,design=design2,save.plot=TRUE) # voom normalize before fitting the linear model
  
  
  name <- gsub("stages","",i)
  
  
  
  # tiff(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/conservative_counts/Diff_v2_MVtrend",currentDate,name,"vs_others.tiff"), type="cairo",
  #  width=10,height=10,units="in",res=300,pointsize = 13,compression="lzw") 
  plot(v2$voom.xy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", 
       pch = 16, cex = 0.25)
  title(paste("voom: Mean-variance trend -",name,"vs other stages",sep=" "))
  lines(v2$voom.line, col = "red")
  # dev.off()
  # 
  
  
  #   #  After this, the usual limma pipelines for differential expression can be applied, for example:
  
  
  fit2 <- lmFit(v2,design2)
  fit2 <- eBayes(fit2)
  
  diff_exp=topTable(fit2,coef=ncol(design2),sort.by = "none",number=nrow(fit2$coefficients))
  
  
  #save here table for later
  
  DE_list[[i]]<-diff_exp
  
  diff_exp_sig=diff_exp[which(diff_exp$adj.P.Val<0.01),]  #pval<0.01
  
  
  
  
  plotMD(fit2,coef=ncol(design2),status=rownames(fit2) %in% rownames(diff_exp_sig),main=paste("Mean-difference plot:",name,"vs other stages"),ylab = "log-fold change")
  
  # dev.off()
  
  print("In red, peaks that are significantly (p<0.01) differentially *expressed*")
  
  
  # Q-Q plot of moderated t-statistics. Points off the line may be differentially expressed
  
  
  
  datap=as.data.frame(rownames(fit2))
  datap$Colour="black"
  # Set new column values to appropriate colours
  datap$Colour[rownames(fit2) %in% rownames(diff_exp_sig)]="red"
  
  
  # tiff(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/conservative_counts/Diff_v2_QQplot_mod_t_stats",currentDate,name,"vs_others.tiff"), type="cairo",
  #   width=10,height=8,units="in",res=300,pointsize = 13,compression="lzw") 
  qqt(fit2$t[,4],df=fit2$df.total,col=datap$Colour,main=paste("Q-Q plot of moderated t-statistics:",name,"vs other stages"))  # careful to select the roght
  # column in data
  abline(0,1)
  #dev.off()
  
  print("In red, peaks that are significantly (p<0.01) differentially *expressed*")
  
  
  
  # Volcano plot
  # tiff(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/conservative_counts/Diff_v2_volcano_",currentDate,name,"vs_others.tiff"), type="cairo",
  #  width=10,height=10,units="in",res=300,pointsize = 13,compression="lzw") 
  volcanoplot(fit2,coef=ncol(design2),highlight=10,names=rownames(fit2),main=paste("Volcano plot:",name,"vs other stages"))
  # dev.off()
  
  
  
  
  
  
}   



# DE_list_sig <- lapply(DE_list, subset, adj.P.Val<0.01)   # list of tables with all significant results

#take logFC and adj P values of each, and combine in single dataframe

combined_df <- lapply(DE_list, "[", c(1,5))   # subsetting list of dataframes with columns I want

combined_df <- do.call("cbind", combined_df)   # merging into one dataframe



#get max conditions per gene
maxVals <- apply(combined_df,1,function(x) which.max(x[c(1,3,5,7,9,11,13,15)]))

#find significant genes per stage
#rowsums checks that there's at least one positive logFC 
# maxvals selects which one to take as max value for each stage
DE_stages<-list()  

sig_iPSC_stage <- combined_df[combined_df$stagesiPSC.adj.P.Val < 0.01 & maxVals == 1,c(1,2)]
sig_iPSC_stage=sig_iPSC_stage[order(sig_iPSC_stage[2]),] #order by adj p values
sig_iPSC_stage$peak_names=rownames(sig_iPSC_stage)
rownames(sig_iPSC_stage)=NULL # take out row names
sig_iPSC_stage=sig_iPSC_stage[c(3,1,2)]
DE_stages [["iPSC"]]<-sig_iPSC_stage

sig_DE_stage <- combined_df[combined_df$stagesDE.adj.P.Val < 0.01 &  maxVals == 2, c(3,4) ]
sig_DE_stage=sig_DE_stage[order(sig_DE_stage[2]),]
sig_DE_stage$peak_names=rownames(sig_DE_stage)
rownames(sig_DE_stage)=NULL # take out row names
sig_DE_stage=sig_DE_stage[c(3,1,2)]
DE_stages [["DE"]]<-sig_DE_stage

sig_PGT_stage <- combined_df[ combined_df$stagesPGT.adj.P.Val < 0.01 & maxVals == 3,c(5,6)  ]
sig_PGT_stage=sig_PGT_stage[order(sig_PGT_stage[2]),]
sig_PGT_stage$peak_names=rownames(sig_PGT_stage)
rownames(sig_PGT_stage)=NULL # take out row names
sig_PGT_stage=sig_PGT_stage[c(3,1,2)]
DE_stages [["PGT"]]<-sig_PGT_stage

sig_PFG_stage <- combined_df[ combined_df$stagesPFG.adj.P.Val < 0.01 & maxVals == 4, c(7,8) ]
sig_PFG_stage=sig_PFG_stage[order(sig_PFG_stage[2]),]
sig_PFG_stage$peak_names=rownames(sig_PFG_stage)
rownames(sig_PFG_stage)=NULL # take out row names
sig_PFG_stage=sig_PFG_stage[c(3,1,2)]
DE_stages [["PFG"]]<-sig_PFG_stage

sig_PE_stage <- combined_df[ combined_df$stagesPE.adj.P.Val < 0.01 & maxVals == 5,c(9,10)  ]
sig_PE_stage=sig_PE_stage[order(sig_PE_stage[2]),]
sig_PE_stage$peak_names=rownames(sig_PE_stage)
rownames(sig_PE_stage)=NULL # take out row names
sig_PE_stage=sig_PE_stage[c(3,1,2)]
DE_stages [["PE"]]<-sig_PE_stage

sig_EP_stage <- combined_df[ combined_df$stagesEP.adj.P.Val < 0.01 & maxVals == 6, c(11,12)  ]
sig_EP_stage=sig_EP_stage[order(sig_EP_stage[2]),]
sig_EP_stage$peak_names=rownames(sig_EP_stage)
rownames(sig_EP_stage)=NULL # take out row names
sig_EP_stage=sig_EP_stage[c(3,1,2)]
DE_stages [["EP"]]<-sig_EP_stage

sig_EN_stage <- combined_df[ combined_df$stagesEN.adj.P.Val < 0.01 & maxVals == 7, c(13,14) ]
sig_EN_stage=sig_EN_stage[order(sig_EN_stage[2]),]
sig_EN_stage$peak_names=rownames(sig_EN_stage)
rownames(sig_EN_stage)=NULL # take out row names
sig_EN_stage=sig_EN_stage[c(3,1,2)]
DE_stages [["EN"]]<-sig_EN_stage

sig_BLC_stage <- combined_df[ combined_df$stagesBLC.adj.P.Val < 0.01 & maxVals == 8, c(15,16)  ]
sig_BLC_stage=sig_BLC_stage[order(sig_BLC_stage[2]),]
sig_BLC_stage$peak_names=rownames(sig_BLC_stage)
rownames(sig_BLC_stage)=NULL # take out row names
sig_BLC_stage=sig_BLC_stage[c(3,1,2)]
DE_stages [["BLC"]]<-sig_BLC_stage


# way too many peaks?
# saving results


DE_stages_logFC1<-list()  # just interested in peaks with logFC>1

for(i in names(DE_stages)){
  
  DE_stages_logFC1[[i]] = DE_stages[[i]][which(DE_stages[[i]][[2]] > 1),]
}


for (i in names(DE_stages_logFC1)){
  write.csv(DE_stages_logFC1[i],quote=F,row.names=F,file=paste("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/DEA/peak",currentDate,"_sig_maxvals_",i,"_diff_peak_analysis.csv",sep=""))
  
}





# across stages (contrasts)

