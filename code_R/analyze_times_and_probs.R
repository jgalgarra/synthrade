# Analysis of numlinks log files
#
# numlinks fields:
#         Simulation Step
#         Number of inks
#         Number of tokens
#         Probability matrix Tukey five number
#         Probability of last existing link
#         Cumulative probability of existing links

library("ggplot2")
library(scales)
source("read_filter_condition.R")
source("parse_command_line_args.R")
# Read simulation log
symlog <- read.csv("../results/symlog.txt", header=FALSE, sep=";")
datalog <- symlog
names(datalog) <- c("Instant", "Year", "NumExper", "SimStep", "TokenCount", "LinkCount", "TotalLinks")

experiment <- 1

years <- seq(ini_seq,end_seq)
for (year in years){
  logfilt <- datalog[datalog$Year==year & datalog$Instant=="FT" & datalog$NumExper == experiment,]
  ftdata <- logfilt[1,]
  file_name <- paste0("numlinks_",year,"_FILT_W_",experiment,".txt")
  simdata <- read.csv(paste0("../results/numlinks/",file_name), header=FALSE, sep=";")
  names(simdata) <- c("SimStep","Links","Tokens","T5Min","T5Q25","T5Median","T5Q75","T5Max",
                      "LastLink","EmptyCells","meanprob","varsigma")
  levplot <- c("Min","Median","Max","LastLink","EmptyCells")
  rowsEvo <- nrow(simdata)*length(levplot)
  EvoDataProb <- data.frame("Step"=rep(0,rowsEvo),"Prob"=rep(0,rowsEvo),"Type"=rep("Min",rowsEvo))
  levels(EvoDataProb$Type) <- levplot
  j=1;
  for (i in seq(1:nrow(simdata))){
    EvoDataProb$Step[j] <- simdata$SimStep[i]
    EvoDataProb$Prob[j] <- simdata$T5Min[i]
    EvoDataProb$Type[j] <- "Min"
    j = j + 1
    EvoDataProb$Step[j] <- simdata$SimStep[i]
    EvoDataProb$Prob[j] <- simdata$T5Median[i]
    EvoDataProb$Type[j] <- "Median"
    j = j + 1
    EvoDataProb$Step[j] <- simdata$SimStep[i]
    EvoDataProb$Prob[j] <- simdata$T5Max[i]
    EvoDataProb$Type[j] <- "Max"
    j = j + 1
    EvoDataProb$Step[j] <- simdata$SimStep[i]
    EvoDataProb$Prob[j] <- simdata$LastLink[i]
    EvoDataProb$Type[j] <- "LastLink"
    j = j + 1
    EvoDataProb$Step[j] <- simdata$SimStep[i]
    EvoDataProb$Prob[j] <- simdata$EmptyCells[i]
    EvoDataProb$Type[j] <- "EmptyCells"
    j = j + 1
  }
  
  if (max(EvoDataProb$Step)<=80000){
    timebreak <- 20000
  } else
    timebreak <- 50000
  
  p <- ggplot(data=EvoDataProb,aes(x=Step, y=10**Prob, color=Type))+geom_line(size=0.7)+
       scale_x_continuous(trans = sqrt_trans(), 
                          breaks = c(50,500,ftdata$SimStep,10000,seq(0,(1+(max(EvoDataProb$Step)%/%timebreak))*timebreak,by=timebreak))) +
       scale_y_continuous(trans = log10_trans(),
                                         breaks = trans_breaks("log10", function(x) 10^x),
                                         labels = trans_format("log10", math_format(10^.x)))+
       ggtitle(year)+xlab("Simulation Step")+ylab("Probability")+
       geom_vline(data=ftdata, aes(xintercept=SimStep),
               linetype="dashed", size=0.7, colour="lightblue")+
       theme_bw() +  
       theme(legend.title = element_blank(), 
          panel.border = element_blank(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.text = element_text(size=9, face="bold"),
          plot.title = element_text(size=12,lineheight=.5, face="bold",hjust = 0.5),
          axis.text = element_text(size=10),
          axis.title.x = element_text(face="bold", size=10),
          axis.title.y  = element_text(face="bold", size=11) )
  
  dir.create("../figures/probevolution/", showWarnings = FALSE)
  ppi <- 300
  png(paste0("../figures/probevolution/",year,"_probevo.png"), width=(8*ppi), height=4*ppi, res=ppi)
  print(p)
  dev.off()
}