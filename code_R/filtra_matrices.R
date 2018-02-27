get_data <- function(file_name,filter=FALSE)
{
r_df <- read.table(paste0("../data/",file_name,".txt"),sep="\t")
r_matrix <- as.matrix(r_df);
r_matrix <- r_matrix[r_matrix>0]
r_aux <- r_matrix[order(r_matrix)]
sum_tot <- sum(r_matrix)
cs <- cumsum(r_aux)
partial <- r_aux[cs < 0.001*sum_tot]
min_allowed <- partial[length(partial)]
maximo <- max(r_matrix)
if (filter)
  r_matrix <- r_matrix[r_matrix >= min_allowed]
maximo <- max(r_matrix)
minimo <- min(r_matrix)
ratio <- maximo/minimo
numlinks <- sum(r_matrix > 0)
print(file_name)
print(paste0("Sum: ",sum(r_matrix)," links: ",numlinks," max:",maximo," min: ", minimo," ratio: ",ratio))
for (i in 1:nrow(r_df))
  for (j in 1:ncol(r_df))
    if (r_df[i,j] <= min_allowed)
      r_df[i,j] = 0
print(paste0("Porcentaje eliminado: ", 100*(sum_tot-sum(r_matrix))/sum_tot))
return(r_df)
}

#files = c("RedAdyCom2000","RedAdyCom1984","RedAdyCom1970")
files = paste0("RedAdyCom",seq(2005,2011))
for (nf in files){
  print("Unfiltered")
  get_data(nf)
  print("Filtered")
  r <- get_data(nf, filter = TRUE)
  write.table(r,paste0("../data/",nf,"_FILT.txt"),row.names = FALSE,col.names = FALSE,sep="\t")
}
