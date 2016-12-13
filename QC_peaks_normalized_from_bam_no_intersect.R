# PCA from normalized atac-seq data in matrix format
# normalized with bamCoverage from Normal bams
# then bigwigsummary, not intersecting bam files with bedtools

library(ggplot2)
library(grid)
library(gridExtra)

# 1st create matrix with all samples and stages (A-D, 1-6)
setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/Results/ATAC-seq/normalized/from_bams")

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


# then PCA  

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
ggsave(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/PCA_deeptools_normalized_normal_bams_not_intersected",currentDate,".jpg",sep=""),p1,width=8,height=8,units="in",dpi=300)


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
ggsave(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/normalised_counts_from_deeptools_distribution_from_normal_bams_no_intersect",currentDate,".jpg",sep=""),p1,width=6,height=5,units="in",dpi=300)


p2 <- ggplot(mdata_aggregated,aes(x=normalized_counts_per_peak, y=number, colour = samples)) +
  geom_point() +  scale_x_continuous(limits=c(0,1000))
ggsave(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/normalised_counts_from_deeptools_distribution_zoomoutfrom_normal_bams_no_intersect",currentDate,".jpg",sep=""),p2,width=6,height=5,units="in",dpi=300)



#####################filter
# how many peaks do I lose with different filter threshold? ###################################

###### fix this


sample=factor(c("sbad2.1","sbad3.1","neo1.1"),levels=c("sbad2.1","sbad3.1","neo1.1"))   # data comes from three donors

stage=factor(c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC"),levels=c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC")) # 8 stages


colnames(counts)=paste(names$stage,names$sample,sep="-")

#re-order before trimming
order_by_stages=function(counts,stage.=stage){
  nc=counts[ , grepl( "iPSC" , names( counts) ) ]  #takes columns whose name matches x
  
  #Do the same for the other stages
  
  for (s in stage.[-1])  {         
    
    i= counts[ , grepl( s , names( counts ) ) ]
    nc=cbind(nc,i)
  }
  
  return(nc)
}


counts=order_by_stages(counts)


old_peaks=dim(counts)[1]# save how many peaks there are initially

trim_by_stage_and_threshold=function(counts,threshold=1,old_peaks=dim(dge)[1]){
  keep <- rowSums(counts>threshold) >= 3   # getting peaks with at least n cpm in at least one stage (3 rep per stage).
  counts <- counts[keep, ]  # not losing anything
  
  # go through the table 3 by 3
  # if there are more than 2 values in three columns above 0, go to next row.
  # if not, go to next 3
  # if gets to end of loop and not found the values, get that row. After loop, subset matrix excluding those rows
  
  data=counts
  
  trim_by_stage=function(data){
    data2=data
    data2[data2<threshold]=0
    data2[data2>threshold]=1 # convert to table where 0=0 and 1 equals any value above n
    
    c=c() # initiate vector of numbers, to save rows to take out of matrix
    
    for (r in 1:nrow(data2)){
      for (i in seq(1,24,3)){   # loop through stages 3 by 3 (to jump over samples of same stage)
        
        if (sum(data2[r,c(i,i+1,i+2)]) > 2) {
          
          break
        }
        if (i==22){
          
          c=c(c,r)  # if it gets to the end and there are no valid values, append to vector
        }
        
        
      }
      
    }
    
    rm(data2)
    return(c)
    
  }
  
  c=trim_by_stage(data)  # calling previous function
  
  fail=100*(length(c)/nrow(counts))
  genes_fail=length(c)
  ditch=counts[c,] # table to throw away
  
  keep=!(rownames(counts) %in% rownames(ditch))
  
  counts <- counts[keep,] # subset for keep
  
  rm(c,keep,ditch)
  
  peaks_remain=((nrow(counts)/old_peaks)*100)
  
  return(peaks_remain)
}

# call function for all CPMs
peaks_remain=c()
for(i in 1:7){
  peaks_remain[i]=trim_by_stage_and_threshold(counts,threshold=i,old_peaks=dim(counts)[1])
}


df=data.frame(threshold = c(1:7),peaks_remain=peaks_remain)

# change function to return counts
trim_by_stage_and_threshold=function(counts,threshold=1,old_peaks=dim(dge)[1]){
  keep <- rowSums(counts>threshold) >= 3   # getting peaks with at least n cpm in at least one stage (3 rep per stage).
  counts <- counts[keep, ]  # not losing anything
  
  # go through the table 3 by 3
  # if there are more than 2 values in three columns above 0, go to next row.
  # if not, go to next 3
  # if gets to end of loop and not found the values, get that row. After loop, subset matrix excluding those rows
  
  data=counts
  
  trim_by_stage=function(data){
    data2=data
    data2[data2<threshold]=0
    data2[data2>threshold]=1 # convert to table where 0=0 and 1 equals any value above n
    
    c=c() # initiate vector of numbers, to save rows to take out of matrix
    
    for (r in 1:nrow(data2)){
      for (i in seq(1,24,3)){   # loop through stages 3 by 3 (to jump over samples of same stage)
        
        if (sum(data2[r,c(i,i+1,i+2)]) > 2) {
          
          break
        }
        if (i==22){
          
          c=c(c,r)  # if it gets to the end and there are no valid values, append to vector
        }
        
        
      }
      
    }
    
    rm(data2)
    return(c)
    
  }
  
  c=trim_by_stage(data)  # calling previous function
  
  fail=100*(length(c)/nrow(counts))
  genes_fail=length(c)
  ditch=counts[c,] # table to throw away
  
  keep=!(rownames(counts) %in% rownames(ditch))
  
  counts <- counts[keep,] # subset for keep
  
  rm(c,keep,ditch)
  
  peaks_remain=((nrow(counts)/old_peaks)*100)
  
  return(counts)
}

counts=trim_by_stage_and_threshold(counts,threshold=1,old_peaks=dim(counts)[1])

# pca plot from filtered counts

### not getting samples right

#########fix filter and PCA!

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
ggsave(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/PCA_deeptools_normalized_normal_bams_not_intersected_filtered",currentDate,".jpg",sep=""),p1,width=8,height=8,units="in",dpi=300)


