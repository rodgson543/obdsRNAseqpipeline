#Day 2: plotting with ggplot2
library(ggplot2)
#Import coding gene bed and change col names
coding_bed <- read.table("/Users/rhodgson/OBDStestdata/Week5/coding_gene_region.bed",  header=FALSE, sep="\t", stringsAsFactors = FALSE)
colnames(coding_bed)
cols <- c("Chr", "Start", "Stop", "Name","Score", "Strand") 
colnames(coding_bed) <- cols
colnames(coding_bed)
#Add column containing Length - will be stop minus start
coding_bed$Length <- coding_bed$Stop - coding_bed$Start
head(coding_bed)


#Plot a histogram of the lengths usign ggplot2
#Add a plot title
#Change the size and angle of the x tick labels
#Change the x and y axis titles and sizes
#Change the x axis upper limit to 500,000
#Change the number of bins or the bin width
#Change the fill and border colour of the bars
Plot1 <- ggplot(coding_bed, aes(x=Length))+
  geom_histogram(bins = 30, color="black", fill="lightblue") +
  xlim(0,500000)+
  labs(title = "Length of Genomic intervals", x= "Length (bp)", y= "Count")+
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)) 
#could also have added scale_x_continuous(format(, scientific =FALSE))

#Using the diamonds dataset, plot a scatter plot of diamond length by
#price, coloured by the diamond colour and sized according to diamond width:
  #Use one of the ggplot built-in themes to alter the plot appearance
#Change the x axis upper limit to 12 and the intervals to 1.5
#Add x and y axis titles and change their size
#Plot the two plots that you have just made side-by-side using a ggplot2 function

#Load data
data(diamonds)
install.packages("wesanderson")
library(wesanderson)
library(RColorBrewer)
plot2 <- ggplot(diamonds, aes(x=x, y=price, colour=color, size=y))+geom_point(shape = 18)+
  scale_x_continuous(breaks = c(0,1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12), limits=c(0,12))+
  labs(title = "Diamond price by size", x= "diamond length (mm)", y= "price ($)")+
  theme_classic()+
  scale_color_manual(values=brewer.pal(n=7, name="Set3"))+
  theme(axis.title = element_text(size = 10))

library(cowplot)
#Plot
dev.off()
plot_grid(Plot1,plot2)  

#Plot as a grid of price by diamond length and color/size
plot2 + facet_grid(cut ~ color)
plot2+facet_wrap(~color, ncol=1)


#Tidyverse
library(tidyverse)
