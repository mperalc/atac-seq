# total mapped counts vs peak number

setwd("/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/")

currentDate <- Sys.Date() # to save date in name of output files


reads <- read.table(file="all_nics_iPSC_differentiated_A1_D6_samples.October_2016.txt", header=T) # tab separated file with coverage from 1 million regions sampled randomly


samples=c("A","B","C","D")
samples_2=rep(samples,each=6)
stage=as.character(c(1:6))
stage_2=rep(stage,4)
colnames(table)[4:ncol(table)]=paste(samples_2,stage_2,sep="")   # samples created


mapped_reads_no_dupl=data.frame(matrix(nrow=24,ncol=4))
mapped_reads_no_dupl$X1=paste(samples_2,stage_2,sep="")

mapped_reads_no_dupl$X2=reads[which(reads$Sample %in% mapped_reads_no_dupl$X1 & reads$Type=="Mapped_no_dupl"),1]
mapped_reads_no_dupl$X3=reads[which(reads$Sample %in% mapped_reads_no_dupl$X1 & reads$Type=="Mapped_reads"),1]
mapped_reads_no_dupl$X4=reads[which(reads$Sample %in% mapped_reads_no_dupl$X1 & reads$Type=="Total_reads"),1]

colnames(mapped_reads_no_dupl)=c("samples","mapped_reads_no_duplicates","mapped_reads","total_reads")


#read in table of names

names <- read.table(file="/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/names_samples.txt",header = T)

# read in peaks file

peaks <- read.table(file="merged_narrowpeaks.sorted.bed")  

# count for every sample in how many rows (peaks) it appears
sums=data.frame(matrix(nrow=24,ncol=2)) 
sums$X1=paste(samples_2,stage_2,sep="")
colnames(sums)=c("samples","total_peaks")

for(s in sums$samples){
  
  print(s)
  
  sums[which(sums$samples==s),2]=sum(grepl(s,peaks[,4]))
  
  
}

# plot dotplot
# x axis: peak number
# y axis: total counts

merged=merge(sums,mapped_reads_no_dupl,by="samples") # merge counts and peaks
merged$samples=names$sample      # rename, names and merged samples are in the same order
merged=cbind(names$stage,merged)  # add stages
colnames(merged)[1]="stage"
merged$samples=factor(merged$samples,levels=c("sbad2.1","sbad3.1","neo1.1"))
merged$stage=factor(merged$stage,levels=c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC"))


#plot

library(ggplot2)
library(grid)
library(gridExtra)



p1 <- ggplot(merged, aes(total_peaks, mapped_reads_no_duplicates, colour=stage, shape=samples)) + 
  geom_point(size = 3) + xlab ("Total peaks") + 
  ylab ("Total mapped reads without duplicates")
p1 <- p1 + theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.background = element_blank(),
                panel.border = element_rect(fill = NA, colour = "black"), 
                legend.key = element_blank(),
                axis.title.y = element_text(face="bold", angle=90, size=10, vjust=0.2),
                axis.title.x = element_text(face="bold", size=10, vjust=0),
                axis.text.x = element_text(face="bold", colour = "black", angle=90, size=10, vjust=0.2, hjust =1 ),
                axis.text.y = element_text(face="bold", colour = "black"),
                axis.ticks = element_line(colour = "black"),
                axis.line = element_line(colour = "black"),
                legend.position = "bottom")
  p1<- p1 + expand_limits(y=0,x=0)  # so that axes include 0
         

  plot(p1)
  ggsave(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/peaks_vs_reads_mapped_wo_duplicates",currentDate,".jpg",sep=""),p1,width=8,height=8,units="in",dpi=300)
  
  
  p2 <- ggplot(merged, aes(total_peaks, mapped_reads, colour=stage, shape=samples)) + 
    geom_point(size = 3) + xlab ("Total peaks") + 
    ylab ("Total mapped reads")
  p2 <- p2 + theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    panel.border = element_rect(fill = NA, colour = "black"), 
                    legend.key = element_blank(),
                    axis.title.y = element_text(face="bold", angle=90, size=10, vjust=0.2),
                    axis.title.x = element_text(face="bold", size=10, vjust=0),
                    axis.text.x = element_text(face="bold", colour = "black", angle=90, size=10, vjust=0.2, hjust =1 ),
                    axis.text.y = element_text(face="bold", colour = "black"),
                    axis.ticks = element_line(colour = "black"),
                    axis.line = element_line(colour = "black"),
                    legend.position = "none")
  p2<- p2 + expand_limits(y=0,x=0)  # so that axes include 0
  
  
  plot(p2)
  
  p3 <- ggplot(merged, aes(total_peaks, total_reads, colour=stage, shape=samples)) + 
    geom_point(size = 3) + xlab ("Total peaks") + 
    ylab ("Total reads")
  p3 <- p3 + theme(	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    panel.border = element_rect(fill = NA, colour = "black"), 
                    legend.key = element_blank(),
                    axis.title.y = element_text(face="bold", angle=90, size=10, vjust=0.2),
                    axis.title.x = element_text(face="bold", size=10, vjust=0),
                    axis.text.x = element_text(face="bold", colour = "black", angle=90, size=10, vjust=0.2, hjust =1 ),
                    axis.text.y = element_text(face="bold", colour = "black"),
                    axis.ticks = element_line(colour = "black"),
                    axis.line = element_line(colour = "black"),
                    legend.position = "none")
  p3<- p3 + expand_limits(y=0,x=0)  # so that axes include 0
  
  
  plot(p3)
  
  
  get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  legend <- get_legend(p1)
  p1 <- p1 + theme(legend.position="none")
  
  plots_grid=grid.arrange(p1, p2, p3, legend, ncol=4, nrow = 3, layout_matrix = rbind(c(1,1,2,2), c(3,3,NA,NA),c(NA,5,5,NA)),
                              widths = c(2,1,1, 2), heights = c(3,3,0.5))
  ggsave(paste("/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/peaks_vs_reads",currentDate,".jpg",sep=""),plots_grid,width=8,height=8,units="in",dpi=300)
  

