# log-transformed histogram of read frequency (normalized) vs fragment size
# as in fig2 Buenrostro et al. 2013

dir = "/Users/Marta/Documents/WTCHG/DPhil/Data/Regulation/atac-seq/fragment_length_counts_filtered_bams/"


samples = paste(rep(c("A", "B", "C", "D"), each = 6), 
                rep(as.character(c(1:6)), 4), sep = "")   # samples created

merged=merge(sums,mapped_reads_no_dupl,by="samples") # merge counts and peaks
merged$samples=names$sample      # rename, names and merged samples are in the same order
merged=cbind(names$stage,merged)  # add stages
colnames(merged)[1]="stage"
merged$samples=factor(merged$samples,levels=c("sbad2.1","sbad3.1","neo1.1"))
merged$stage=factor(merged$stage,levels=c("iPSC","DE","PGT","PFG","PE","EP","EN","BLC"))


samples_list = list()
read_number = list()
for (i in samples) {
  samples_list[[i]] = read.table(paste(dir,i,"_fragment_length_count.txt", sep = ""))
  colnames(samples_list[[i]]) = c("count","fragment")
  samples_list[[i]]["sample"] = rep(i, nrow(samples_list[[i]]))
  read_number[[i]] = sum(samples_list[[i]]["count"])
  
  # normalize: count / total counts
  samples_list[[i]]["count"] = samples_list[[i]]["count"] / sum(samples_list[[i]]["count"])
}

samples_list = do.call("rbind", samples_list)
read_number = do.call("rbind", read_number)

####### plot
# transform 1000 more freq for x axis clarity

toPlot = samples_list
toPlot$count = toPlot$count * 1000

p <-
  ggplot(data = toPlot, aes(x = fragment, y = count, group = sample)) +
  #ggtitle(unique(long$GeneName)) +
  # xlab("Differentiation stages") +
  # ylab("Normalized read density x ") +
  labs(x = "Fragment length (bp)", y = expression(paste("Normalized read density"," ", x, " ", 10^{-3}))) +
  expand_limits(y = 0) +
  scale_x_continuous(limits = c(0, 1000)) +
  # geom_hline(
  #   yintercept = 0,
  #   linetype = "dashed",
  #   size = 1,
  #   col = "#DCDCDC"
  # ) +
 # geom_density() +
   geom_line(
  #   #aes(linetype = sample, col = sample), 
     size = 1) +
  #scale_colour_manual(values = diaPalette) +  # diabetologia pallete
  #scale_colour_manual(values="#000000") +  # for black and white, otherwise map lines and point colours to samples
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(size = 2),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold.italic"),
    legend.text = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 13, face = "bold")
  )
plot(p)

# in log scale

toPlotlog = samples_list
#toPlotlog$count = log(toPlotlog$count)
pal2 = c(rep("#002147",4),"red", "#002147", # A
        rep("#002147",6), # B
        "yellow",rep("#002147",2), "pink",rep("#002147",2), # C
        rep("#002147",3),"red",rep("#002147",2)) # D  
 p <-
  ggplot(data = toPlotlog, aes(x = fragment, y = count, group = sample, color = sample)) +
  #ggtitle(unique(long$GeneName)) +
  # xlab("Differentiation stages") +
  # ylab("Normalized read density x ") +
  labs(x = "Fragment length (bp)", y = expression("Normalized read density")) +
  expand_limits(y = 0) +
  scale_x_continuous(limits = c(0, 1000)) +
  scale_color_manual(values = pal2) +
  
 
  scale_y_log10(breaks = trans_breaks("log10",  function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(scaled = T, sides = "l") +
  # geom_hline(
  #   yintercept = 0,
  #   linetype = "dashed",
  #   size = 1,
  #   col = "#DCDCDC"
  # ) +
  # geom_density() +
  geom_line(
    #   #aes(linetype = sample, col = sample), 
    size = 1) +
  #scale_colour_manual(values = diaPalette) +  # diabetologia pallete
  #scale_colour_manual(values="#000000") +  # for black and white, otherwise map lines and point colours to samples
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(size = 2),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold.italic"),
    legend.text = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 13, face = "bold")
  )
plot(p) 
#

png(
  paste(
    "/Users/Marta/Documents/WTCHG/DPhil/Plots/atac-seq/",
    "log10_norm_read_density_from_filtered_bams",
    ".png",
    sep = ""
  ),
  type = "cairo",
  width = 8,
  height = 5,
  units = "in",
  res = 300,
  pointsize = 12
)
print(p)
dev.off()
