
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

for (year in years) {
  logfilt <- datalog[datalog$Year==year & datalog$Instant=="FT" & datalog$NumExper == experiment,]
  ftdata <- logfilt[1,]
  file_name <- paste0("numlinks_",year,"_FILT_W_",experiment,".txt")
  simdata <- read.csv(paste0("../results/numlinks/",file_name), header=FALSE, sep=";")
  names(simdata) <- c("SimStep","Links","Tokens","Exporters","Importers","T5Min","T5Q25","T5Median","T5Q75","T5Max",
                      "LastLink","EmptyCells","meanprob","varsigma")
  
  # Relación lineal entre tokens y 1/p
  plot(1/exp(simdata$EmptyCells),simdata$Tokens, main="Tokens vs prob empty")
  
  # Relación cuadrática entre links y 1/p
  plot(1/exp(simdata$EmptyCells),simdata$Links^2, main="Links vs prob empty")
}