library(grid)
library(gridExtra)
library(ggplot2)
library(kcorebip)

file_name <- "RedAdyCom????.txt"
data_files <- Sys.glob(paste0("../data/",file_name))
nldf <- data.frame("file"=c(),"year"=c(),"links"=c(), "connectance" = c())
for (i in data_files)
{
  raw_data <- read.table(i)
  max_possible_links <- nrow(raw_data)*ncol(raw_data)
  nlinks <- sum(raw_data>0)
  nfile <- strsplit(strsplit(i,"data/")[[1]][2],".txt")[[1]][1]
  year <- strsplit(strsplit(i,"Com")[[1]][2],".txt")[[1]][1]
  nldf <- rbind(nldf,data.frame("file"=nfile,"year"=year,"links"=nlinks, "connectance" = nlinks/max_possible_links))
}

write.table(nldf,"../results/NUMLINKS.txt",row.names = FALSE)