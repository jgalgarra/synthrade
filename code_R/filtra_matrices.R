get_data <- function(file_name,filter_factor = 0)
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
#min_allowed <-  maximo*(10^-filter_factor)
if (filter_factor != 0)
  r_matrix <- r_matrix[r_matrix >= min_allowed]
maximo <- max(r_matrix)
minimo <- min(r_matrix)
ratio <- maximo/minimo
numlinks <- sum(r_matrix > 0)
print(file_name)
print(paste0("Sum: ",sum(r_matrix)," links: ",numlinks," max:",maximo," min: ", minimo," ratio: ",ratio))
if (filter_factor != 0){
  for (i in 1:nrow(r_df))
    for (j in 1:ncol(r_df))
      if (r_df[i,j] <= min_allowed)
        r_df[i,j] = 0
  print(paste0("Porcentaje eliminado: ", 100*(sum_tot-sum(r_matrix))/sum_tot))
  return(r_df)
  }
}

#files = c("RedAdyCom2000","RedAdyCom1984","RedAdyCom1970")
files = paste0("RedAdyCom",seq(1962,2014))
for (nf in files){
  nivel_filtrado = 1
  print("Sin filtro")
  get_data(nf)
  print("Filtrado")
  r <- get_data(nf, filter_factor = nivel_filtrado)
  write.table(r,paste0("../data/",nf,"_ff_",nivel_filtrado,".txt"),row.names = FALSE,col.names = FALSE,sep="\t")
}
