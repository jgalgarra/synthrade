
files <- seq(1962,2014)
for (nfile in files)
{
  print(nfile)
  raw_data <- read.table(paste0("../data/RedCom",nfile,".txt"))
  names(raw_data) <- c("or","dest","imp")
  clean_data <- raw_data[raw_data$or != 'wld',]
  clean_data <- clean_data[clean_data$dest != 'wld',]
  countries1 <- unique(clean_data$or)
  countries2 <- unique(clean_data$dest)
  if (length(countries1) > length(countries2))
    countries = countries1
  else
    countries = countries2
  dim = length(countries)
  B <- matrix( rep(0,dim*dim), nrow=dim, ncol=dim)
  dfsal = as.data.frame(B)
  rownames(dfsal) <- countries
  colnames(dfsal) <- countries
  for (j in 1:nrow(clean_data))
    dfsal[as.character(clean_data$or[j]),as.character(clean_data$dest[j])] <- as.numeric(clean_data$imp[j])
  write.csv(dfsal,paste0("../data/RedAdyNames",nfile,".txt"),row.names = TRUE)
  write.table(dfsal,paste0("../data/RedAdyCom",nfile,".txt"),row.names = FALSE, col.names = FALSE, sep ="\t")
}