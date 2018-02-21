get_data <- function(file_name,filter_factor = 0)
{
r_df <- read.table(paste0("../data/",file_name,".txt"),sep="\t")
r_matrix <- as.matrix(r_df);
r_matrix <- r_matrix[r_matrix>0]
maximo <- max(r_matrix)
min_allowed <-  maximo*(10^-filter_factor)
if (filter_factor != 0)
  r_matrix <- r_matrix[r_matrix > min_allowed]
maximo <- max(r_matrix)
minimo <- min(r_matrix)
ratio <- maximo/minimo
numlinks <- sum(r_matrix > 0)
print(file_name)
print(paste0("Sum: ",sum(r_matrix)," links: ",numlinks," max:",maximo," min: ", minimo," ratio: ",ratio))
for (i in 1:nrow(r_df))
  for (j in 1:ncol(r_df))
    if (r_df[i,j] < min_allowed)
      r_df[i,j] = 0
return(r_df)
}

print("Sin filtro")
get_data("RedAdyCom2012")
get_data("RedAdyCom2013")
get_data("RedAdyCom2014")

print("Filtrado")
r <- get_data("RedAdyCom2012", filter_factor = 6)
get_data("RedAdyCom2013", filter_factor = 6)
get_data("RedAdyCom2014", filter_factor = 6)