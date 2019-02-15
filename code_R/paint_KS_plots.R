library(grid)
library(gridExtra)
library(ggplot2)

PaintKS <- function(series,fillcol,title)
{
  r <- ggplot(data=series,aes(Year,p.value)) +
    #geom_point( color=fillcol,size=1.5,shape=21)+
    geom_boxplot(aes(group=Year),fill=fillcol)+
    geom_hline(aes(yintercept=0.1), colour="green", linetype="dashed", alpha = 0.8, size = 1)+
    ggtitle(title)+
    theme_bw()
  print(r)
}

KSdata <- read.table("../results/KSTEST.txt",header=TRUE);
Exdata <- data.frame("Year" = KSdata$Year, "p.value" = KSdata$KSexport)
PaintKS(Exdata,"blue","Exporter")
Impdata <- data.frame("Year" = KSdata$Year, "p.value" = KSdata$KSimport)
PaintKS(Impdata,"orange","Importer")