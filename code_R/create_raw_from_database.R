library(sqldf)
year <- 1983
dataf <- read.delim(paste0("../data/raw_data/",year,".txt"), header=FALSE)
names(dataf) <- c("year","origin","dest","prod","export_val","import_val")
u <- sqldf("select orig,dest,money from dataf group by origin,dest,prod")

