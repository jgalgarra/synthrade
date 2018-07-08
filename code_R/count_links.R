library(grid)
library(gridExtra)
library(ggplot2)
library(kcorebip)

file_name <- "RedAdyCom????_FILT.txt"
data_files <- Sys.glob(paste0("../data/",file_name))
nldf <- data.frame("file"=c(),"year"=c(),"nexp"=c(),
                  "nimp"=c(),"links"=c(), "connectance" = c())
for (i in data_files)
{
  raw_data <- read.table(i)
  max_possible_links <- sum(rowSums(raw_data)>0)*sum(colSums(raw_data)>0)
  nlinks <- sum(raw_data>0)
  nfile <- strsplit(strsplit(i,"data/")[[1]][2],".txt")[[1]][1]
  year <- gsub("_FILT","",strsplit(strsplit(i,"Com")[[1]][2],".txt")[[1]][1])
  nldf <- rbind(nldf,data.frame("file"=nfile,"year"=year,"links"=nlinks, 
                                "nexp"=sum(rowSums(raw_data)>0),"nimp"=sum(colSums(raw_data)>0),
                                "connectance" = nlinks/max_possible_links))
}
plot(nldf$year,nldf$connectance)
plot(nldf$year,nldf$links)
print(paste("Average connectnce",mean(nldf$connectance)))

write.table(nldf,"../results/NUMLINKS.txt",row.names = FALSE)