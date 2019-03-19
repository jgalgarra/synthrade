library(grid)
library(gridExtra)
library(ggplot2)

PaintKS <- function(series,fillcol,title)
{
  r <- ggplot(data=series,aes(Year,p.value)) +
    geom_boxplot(aes(group=Year),fill=fillcol, alpha = 0.5)+
    geom_hline(aes(yintercept=0.1), colour="green", alpha = 0.8, size = 1)+ylab("KS test p.value\n")+xlab("")+
    ggtitle(title)+
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
          axis.text.y = element_text(face="bold", size=13),
          axis.text.x = element_text(face="bold", size=13, angle = 45, hjust = 1),
          axis.title.x = element_text(face="bold", size=13),
          axis.title.y  = element_text(face="bold", size=13) )
  return(r)
}

KSdata <- read.table("../results/KSTEST.txt",header=TRUE);
Exdata <- data.frame("Year" = KSdata$Year, "p.value" = KSdata$KSexport)
KSexp <- PaintKS(Exdata,"blue","Exporters Strength")
Impdata <- data.frame("Year" = KSdata$Year, "p.value" = KSdata$KSimport)
KSImp <- PaintKS(Impdata,"red","Importers Strength")
dir.create("../figures/", showWarnings = FALSE)
dir.create("../figures/tests/", showWarnings = FALSE)
fsal <- paste0("../figures/tests/KSplots.png")
ppi <- 600
png(fsal, width=10*ppi, height=8*ppi, res=ppi)
grid.arrange(KSexp, KSImp, ncol=1, nrow=2,top="" )
dev.off()