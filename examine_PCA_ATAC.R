# playing with PCA atac-seq
library(ggplot2)

load(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/session_objects/dge_atac-seq.xz",verbose=TRUE)  #loading again the dge object to do the selected QC


dge <- calcNormFactors(dge) # Calculating normalization factors


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
  peaks_fail=length(c)
  ditch=dge$counts[c,] # table to throw away
  
  keep=!(rownames(dge$counts) %in% rownames(ditch))
  
  dge <- dge[keep, , keep.lib.sizes=FALSE] # counts now has 17892 peaks
  
  rm(c,keep,ditch)
  
  
  
  return(dge)
}

# call function for CPM=1

# call trimming function for CPM=1

dge <- trim_by_stage_and_cpm(dge,CPM=1,old_peaks=dim(dge)[1])
dge <- calcNormFactors(dge)# recalculating normalization factors after altering library sizes

nrow(dge)/old_peaks #proportion of peaks that remain

sample=factor(c("sbad2.1","sbad3.1","neo1.1"),levels=c("sbad2.1","sbad3.1","neo1.1"))   # data comes from three donors

stage=factor(c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC"),levels=c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC")) # 8 stages

samples <- rep(sample,8)   

stages <-rep(stage,each=3) 

#create the design matrix
design1 <- model.matrix(~stages + samples)                    
kable(design1)


# This converts the counts to log-counts per million with associated precision weights. After this, the RNA-seq data can be analyzed as if it was microarray data.


v2=voom(dge,design=design1,plot=TRUE)   


par(mfrow=c(1,1))

apply(X = v2$E,2,var)
# variances are quite different, so R book recommends scaling
# scale = T?
pca1<-prcomp(t(v2$E), retx=TRUE,scale. = F)
plot(pca1, type = "l") #variance vs first 10 components
summary(pca1)          #importance of each component (important line is "proportion of variance")
# biplot(pca1) # too many peaks to plot

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



p
