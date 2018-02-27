nfiles <- seq(2011,2011)
for (file in nfiles){
  ctdata <- read.csv(paste0("../data/raw_data/TodosYR",file,".csv"), sep=";")
  outdata <- data.frame("or"=ctdata$origin,"dest"=ctdata$dest,"imp"=ctdata$import_val)
  write.table(outdata,paste0("../data/raw_data/RedRaw",file,".txt"),sep=" ",row.names = FALSE, col.names = FALSE)
}
