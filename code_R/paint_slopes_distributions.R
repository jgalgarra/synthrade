# Cumulative degree-cumulative strength distributions 
#  
# Author: Javier Garcia Algarra
#
# Invocation: Rscript paint_slopes_distributions
#                   
#
#  Plot at /figures/linksstrength/ALL_FILTERED_Slopes.png

library(grid)
library(gridExtra)
library(ggplot2)
source("parse_command_line_args.R")

slopes <- read.delim("../results/Slopes.txt")
slopes_emp <- slopes[slopes$Experiment==1,]
dataslope <- data.frame("Year"=c(),"slope"=c(),"Node"=c())
for (i in 1:nrow(slopes_emp)){
  datatemp <- data.frame("Year"=slopes_emp$Year[i],"slope"=slopes_emp$ExpSlopeEmp[i],"Node"="Exporter")
  dataslope <- rbind(dataslope,datatemp)
  datatemp <- data.frame("Year"=slopes_emp$Year[i],"slope"=slopes_emp$ImpSlopeEmp[i],"Node"="Importer")
  dataslope <- rbind(dataslope,datatemp)
}

p <- ggplot(data=dataslope,aes(Node,slope)) + 
     geom_boxplot(aes(group=Node,color = Node,fill=Node),alpha=0.3) +
     ylab("Slope")+xlab("")+ggtitle("Cumulative Normalized Strength vs Degree")+
     theme_bw()+
     theme(panel.border = element_blank(),
        legend.key = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(linetype = 2, color="ivory3"),
        panel.grid.major.x = element_blank(), 
        legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=12, face="bold"),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(lineheight=.8, size=12, face="bold",hjust = 0.5),
        axis.text = element_text(face="bold", size=13),
        axis.title.x = element_text(face="bold", size=13),
        axis.title.y  = element_text(face="bold", size=13) )

dir.create("../figures/linksstrength/", showWarnings = FALSE)
fsal <- paste0("../figures/linksstrength/ALL_FILTERED_Slopes.png")
ppi <- 600
png(fsal, width=6*ppi, height=6*ppi, res=ppi)
print(p)
dev.off()