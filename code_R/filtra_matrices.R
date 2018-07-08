get_data <- function(file_name,filter=FALSE, connectance = 0)
{
r_df <- read.table(paste0("../data/",file_name,".txt"),sep="\t")
r_df <- r_df[rowSums(r_df)>0,]
r_df <- r_df[,colSums(r_df)>0]
r_matrix <- as.matrix(r_df);
r_matrix <- r_matrix[r_matrix>0]
r_aux <- r_matrix[order(r_matrix)]
sum_tot <- sum(r_matrix)
cs <- cumsum(r_aux)
partial <- r_aux[cs < 0.01*sum_tot]
min_allowed <- partial[length(partial)]
if (filter)
  r_matrix <- r_matrix[r_matrix > min_allowed]


maximo <- max(r_matrix)
minimo <- min(r_matrix)
ratio <- maximo/minimo
numlinks <- sum(r_matrix > 0)
print(file_name)
print(paste0("Sum: ",sum(r_matrix)," orig_links ",sum(r_df>0)," filt links: ",numlinks," max:",maximo," min: ", minimo," ratio: ",ratio))
for (i in 1:nrow(r_df))
  for (j in 1:ncol(r_df))
    if (r_df[i,j] < min_allowed)
      r_df[i,j] = 0
print(paste0("Porcentaje eliminado: ", 100*(sum_tot-sum(r_matrix))/sum_tot))
return(r_df)
}

files = paste0("RedAdyCom",seq(1965,2014))
for (nf in files){
  # print("Unfiltered")
  # get_data(nf)
  print("Filtered")
  r <- get_data(nf, filter = TRUE)
  r <- r[,colSums(r)>0]
  r <- r[rowSums(r)>0,]
  print(paste("connectance",sum(r>0)/(nrow(r)*ncol(r)) ))
  write.table(r,paste0("../data/",nf,"_FILT.txt"),row.names = FALSE,col.names = FALSE,sep="\t")
}
