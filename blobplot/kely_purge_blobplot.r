# this script is used to plot the blobplot for Keiferia genome
# please modify this script for the use on the other cases
library(ggplot2)


path <- getwd()
contig_dat <- read.table(paste0(path, "/Kely_purge_summary.tsv"), header=T, sep="\t")

x11()
ggplot(contig_dat, aes(x=gc, y=purging_aligned_cov, size=length, col=bestsumorder_phylum)) + 
  scale_color_manual(values=c("blue", "darkgreen")) +
  geom_point(alpha=0.2) +
  scale_size(range = c(2, 10), name="length") +
  scale_x_continuous(name="Speed of cars", limits=c(0, 1)) +
  ylab("Contig coverage") +
  xlab("GC proportion") +
  guides(size = "none") +
  guides(color=guide_legend(override.aes = list(size=3), title="best blast phylum")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

