# PCA from normalized atac-seq data in matrix format

library(ggplot2)
library(grid)
library(gridExtra)

# 1st create matrix with all samples and stages (A-D, 1-6)
setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/normalized/from_peaks/reduced_genomesize")

currentDate <- Sys.Date() # to save date in name of output files
  
samples=c("A","B","C","D")
samples_2=rep(samples,each=6)
stage=as.character(c(1:6))
stage_2=rep(stage,4)

table <- read.table(file="normalized_matrix_peaks_multiBigwigSummary.tab") # first 3 columns: chr, start and end of peak

colnames(table)[4:ncol(table)]=paste(samples_2,stage_2,sep="")  # give names
colnames(table)[1:3] = c("chr","start","end")

counts=table[4:ncol(table)]
#read in table of names

names <- read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/names_samples.txt",header = T)


# then PCA  (NO SCALING)

#PC1 vs 2

plot_pca=function(x){
  pca1<-prcomp(t(x), retx=TRUE,scale. = F)   # transpose so that samples are in rows, and genes/peaks, etc in columns
  
  # plot(pca1, type = "l") #variance vs first 10 components
  summary(pca1)          #importance of each component (important line is "proportion of variance")
  
  percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs,sample=names$sample,stage=names$stage)
  pcs$sample=factor(pcs$sample,levels=c("sbad2.1","sbad3.1","neo1.1"))
  pcs$stage=factor(pcs$stage,levels=c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC"))
  # pcs$stage <- ordered(pcs$stage, levels = stage)
  
  p <- ggplot(pcs, aes(PC1, PC2, colour=stage, shape=sample)) + 
    geom_point(size = 3) + xlab (paste0( "PC1:" ,percentVar[ 1 ],"% variance")) + 
    ylab (paste0( "PC2: ",percentVar[ 2 ],"% variance" ))
  p <- p + theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                  panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, colour = "black"), 
                  legend.key = element_blank(),# legend.position = c(0.5,0.5),
                  axis.title.y = element_text(face="bold", angle=90, size=12, vjust=0.2),
                  axis.title.x = element_text(face="bold", size=12, vjust=0),
                  axis.text.x = element_text(face="bold", colour = "black", angle=90, size=12, vjust=0.2, hjust =1 ),
                  axis.text.y = element_text(face="bold", colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  axis.line = element_line(colour = "black"),
                  legend.position = "bottom")
  
  return(p)
}

p1=plot_pca(counts)
ggsave(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/PCA1",currentDate,".jpg",sep=""),p1,width=8,height=8,units="in",dpi=300)



#PC2 vs 3
plot_pca=function(x){
  pca1<-prcomp(t(x), retx=TRUE,scale. = F)   # transpose so that samples are in rows, and genes/peaks, etc in columns
  
  # plot(pca1, type = "l") #variance vs first 10 components
  summary(pca1)          #importance of each component (important line is "proportion of variance")
  
  percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs,sample=names$sample,stage=names$stage)
  pcs$sample=factor(pcs$sample,levels=c("sbad2.1","sbad3.1","neo1.1"))
  pcs$stage=factor(pcs$stage,levels=c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC"))
  # pcs$stage <- ordered(pcs$stage, levels = stage)
  
  p <- ggplot(pcs, aes(PC2, PC3, colour=stage, shape=sample)) + 
    geom_point(size = 3) + xlab (paste0( "PC2:" ,percentVar[ 2 ],"% variance")) + 
    ylab (paste0( "PC3: ",percentVar[ 3 ],"% variance" ))
  p <- p + theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                  panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, colour = "black"), 
                  legend.key = element_blank(),# legend.position = c(0.5,0.5),
                  axis.title.y = element_text(face="bold", angle=90, size=12, vjust=0.2),
                  axis.title.x = element_text(face="bold", size=12, vjust=0),
                  axis.text.x = element_text(face="bold", colour = "black", angle=90, size=12, vjust=0.2, hjust =1 ),
                  axis.text.y = element_text(face="bold", colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  axis.line = element_line(colour = "black"),
                  legend.position="none")
  
  return(p)
}

p2=plot_pca(counts)


#PC3 vs 4
plot_pca=function(x){
  pca1<-prcomp(t(x), retx=TRUE,scale. = F)   # transpose so that samples are in rows, and genes/peaks, etc in columns
  
   #plot(pca1, type = "l") #variance vs first 10 components
  summary(pca1)          #importance of each component (important line is "proportion of variance")
  
  percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs,sample=names$sample,stage=names$stage)
  pcs$sample=factor(pcs$sample,levels=c("sbad2.1","sbad3.1","neo1.1"))
  pcs$stage=factor(pcs$stage,levels=c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC"))
  # pcs$stage <- ordered(pcs$stage, levels = stage)
  
  p <- ggplot(pcs, aes(PC1, PC3, colour=stage, shape=sample)) + 
    geom_point(size = 3) + xlab (paste0( "PC1:" ,percentVar[ 1 ],"% variance")) + 
    ylab (paste0( "PC3: ",percentVar[ 3 ],"% variance" ))
  p <- p + theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                  panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, colour = "black"), 
                  legend.key = element_blank(),# legend.position = c(0.5,0.5),
                  axis.title.y = element_text(face="bold", angle=90, size=12, vjust=0.2),
                  axis.title.x = element_text(face="bold", size=12, vjust=0),
                  axis.text.x = element_text(face="bold", colour = "black", angle=90, size=12, vjust=0.2, hjust =1 ),
                  axis.text.y = element_text(face="bold", colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  axis.line = element_line(colour = "black"),
                  legend.position="none")
  
  return(p)
}

p3=plot_pca(counts)


plot_pca=function(x){
  pca1<-prcomp(t(x), retx=TRUE,scale. = F)   # transpose so that samples are in rows, and genes/peaks, etc in columns
  
  #plot(pca1, type = "l") #variance vs first 10 components
  summary(pca1)          #importance of each component (important line is "proportion of variance")
  
  percentVar <- (pca1$sdev)^2 / sum(pca1$sdev^2)
  percentVar <- round(100 * percentVar)
  pcs <- as.data.frame(pca1$x)
  pcs <- cbind(pcs,sample=names$sample,stage=names$stage)
  pcs$sample=factor(pcs$sample,levels=c("sbad2.1","sbad3.1","neo1.1"))
  pcs$stage=factor(pcs$stage,levels=c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC"))
  # pcs$stage <- ordered(pcs$stage, levels = stage)
  
  p <- ggplot(pcs, aes(PC3, PC4, colour=stage, shape=sample)) + 
    geom_point(size = 3) + xlab (paste0( "PC3:" ,percentVar[ 3 ],"% variance")) + 
    ylab (paste0( "PC4: ",percentVar[ 4 ],"% variance" ))
  p <- p + theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                  panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, colour = "black"), 
                  legend.key = element_blank(),# legend.position = c(0.5,0.5),
                  axis.title.y = element_text(face="bold", angle=90, size=12, vjust=0.2),
                  axis.title.x = element_text(face="bold", size=12, vjust=0),
                  axis.text.x = element_text(face="bold", colour = "black", angle=90, size=12, vjust=0.2, hjust =1 ),
                  axis.text.y = element_text(face="bold", colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  axis.line = element_line(colour = "black"),
                  legend.position="none")
  
  return(p)
}

p4=plot_pca(counts)


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(p1)
p1 <- p1 + theme(legend.position="none")

pca_no_scaling=grid.arrange(p1, p2, p3, p4, legend, ncol=4, nrow = 3, layout_matrix = rbind(c(1,1,2,2), c(3,3,4,4),c(NA,5,5,NA)),
             widths = c(2,1,1, 2), heights = c(3,3,0.5))
ggsave(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/PCA_all_samples",currentDate,".jpg",sep=""),pca_no_scaling,width=8,height=8,units="in",dpi=300)

# MDS version

library(limma)
library(ggrepel) 
#function wrap so that plotMDS doesn't produce a plot when called:
plotMDS.invisible <- function(...){
  ff <- tempfile()
  png(filename=ff)
  res <- plotMDS(...)
  dev.off()
  unlink(ff)
  res
}

pretty_mds=function(table){
  mds_p=plotMDS.invisible(counts,gene.selection = "pairwise")    # pairwise method (default)
  mds_c=plotMDS.invisible(counts,gene.selection = "common")      #common method 
  
  
  # Rearrange data for ggplot
  
  # method: pairwise
  m_p=as.data.frame(mds_p$cmdscale.out)
  m_p <- cbind(m_p,sample=names$sample,stage=names$stage)
  colnames(m_p)=c("Dimension 1", "Dimension 2", "Samples","Stages")
  
  # method: common
  m_c=as.data.frame(mds_c$cmdscale.out)
  m_c <- cbind(m_c,sample=names$sample,stage=names$stage)
  colnames(m_c)=c("Dimension 1", "Dimension 2", "Samples","Stages")
  
  
  # plot pairwise
  

  mp=ggplot(m_p) +
    geom_point(aes(`Dimension 1` ,`Dimension 2`), size = 2, color = 'grey') +
    geom_label_repel(
      aes(`Dimension 1` , `Dimension 2`, fill = Stages, label = Samples),
      fontface = 'bold', color = 'white',
      box.padding = unit(0.25, "lines"),
      point.padding = unit(0.25, "lines")
    ) +
    coord_fixed(ratio = 4) +
    
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
      aes(`Dimension 1` , `Dimension 2`, fill = Stages, label = Samples),
      fontface = 'bold', color = 'white',
      box.padding = unit(0.25, "lines"),
      point.padding = unit(0.25, "lines")
    ) +
    coord_fixed(ratio = 4) +  #fix x-y ratio
    
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

mds_plots=pretty_mds(counts)

plot(mds_plots[[1]])
plot(mds_plots[[2]])



####### investigate distribution of counts

mean_peak_length=mean(table$end-table$start)  # calculate mean peak length
median_peak_length=median(table$end-table$start)  # calculate mean peak length

each_sample_counts_mean=colMeans(counts[1:ncol(counts)]) # calculate mean counts per sample
appended_columns=unlist(counts[1:ncol(counts)],use.names = F)
all_sample_counts_mean=mean(appended_columns)   # calculate mean counts for all samples
# also, for each peak:
each_peak_counts_mean=rowMeans(counts[1:ncol(counts)])
mean(each_peak_counts_mean)  # gives same result as above

all_sample_counts_variance=var(appended_columns) # calculate variance for counts for all samples

# mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

mode= getmode(appended_columns)

####plot distribution of reads per peak

#### empty rows when I aggregate, not recognised as empty or NA
# what is going on??

library(ggplot2)
library(grid)
library(gridExtra)

unique_values_counts=aggregate(data.frame(count = appended_columns), list(value = appended_columns), length)

colnames(unique_values_counts)=c("counts_in_peak","number")

sum(unique_values_counts$number)  # total times a counts_in_peak values is repeated 
# 5570832
length(unique_values_counts$counts_in_peak)
# frequency for each unique value of counts per peak:

unique_values_counts$frequency=unique_values_counts$number/sum(unique_values_counts$number)

sum(unique_values_counts$frequency) # sum of all frequencies gives 1 

p1 <- ggplot(unique_values_counts, aes(x=counts_in_peak, y=frequency))
p1 <- p1 + geom_point() + scale_x_continuous(limits = c(0, 100))
plot(p1)

# another way

library(reshape2)
library(plyr)
library(geneplotter)
mdata <- melt(counts,na.rm = T) 

mdata_aggregated=count(mdata, c("variable", "value"))
# sum(mdata_aggregated[which(mdata_aggregated$value==0.0000),3])  # sum of all 0s in all samples is equal to the total sum of 0s in unique_values_counts

colnames(mdata_aggregated)=c("samples","normalized_counts_per_peak","number")
# fix 
#mdata_aggregated$frequency_per_sample=mdata_aggregated$number/sum(mdata_aggregated$number[which(mdata_aggregated)])


p1 <- ggplot(mdata_aggregated,aes(x=normalized_counts_per_peak, y=number, colour = samples)) +
  geom_point() +  scale_x_continuous(limits=c(0,10),breaks=c(1:10))
ggsave(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/normalised_counts_from_deeptools_distribution_",currentDate,".jpg",sep=""),p1,width=6,height=5,units="in",dpi=300)


p2 <- ggplot(mdata_aggregated,aes(x=normalized_counts_per_peak, y=number, colour = samples)) +
  geom_point() +  scale_x_continuous(limits=c(0,1000))
ggsave(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/normalised_counts_from_deeptools_distribution_zoomout",currentDate,".jpg",sep=""),p2,width=6,height=5,units="in",dpi=300)

