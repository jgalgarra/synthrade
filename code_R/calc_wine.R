library(grid)
library(gridExtra)
library(ggplot2)
library(kcorebip)

MPack <- function(matrix,normalize = TRUE)
{
  sum_row <- rep(0,nrow(matrix))
  sum_col <- rep(0,ncol(matrix))
  if (normalize)
    matrix = matrix/max(matrix)
  for (i in 1:nrow(matrix))
    sum_row[i] <- sum(matrix[i,])
  for (i in 1:ncol(matrix))
    sum_col[i] <- sum(matrix[,i])
  ord_matrix <- matrix[rev(order(sum_row)),rev(order(sum_col))]
  return(t(ord_matrix))       # Transpose because of order of python-written matrix
}


NREPS <- 1000
orig_file <- "RedAdyCom2014"
file_name <- paste0(orig_file,"_ff_1")
#file_name <- "kaka1"
experiment_files <- Sys.glob(paste0("../results/",file_name,"_W_*.txt"))

# zero_matrix <- read.table(experiment_files[1],sep="\t")
# for (l in 1:nrow(zero_matrix))
#   for(m in 1:ncol(zero_matrix))
#     zero_matrix[l,m]<-0

numexper <- length(experiment_files)


or_matrix <- read.table(paste0("../data/",orig_file,".txt"),sep="\t")
sum_row <- rep(0,nrow(or_matrix))
sum_col <- rep(0,ncol(or_matrix))
ind_matrix_p <- MPack(or_matrix)
# Remove all zeroes columns and rows
dfint <- as.data.frame(ind_matrix_p)
w <- wine(dfint,nreps=NREPS)
obswine <- w$wine
print(paste0(orig_file," wine ",obswine))

emp_matrix <- read.table(paste0("../data/",file_name,".txt"),sep="\t")
sum_row <- rep(0,nrow(emp_matrix))
sum_col <- rep(0,ncol(emp_matrix))
dfanid <- data.frame("wine"=c(),"exper"=c())
ind_matrix_p <- MPack(emp_matrix)
# Remove all zeroes columns and rows
dfint <- as.data.frame(ind_matrix_p)
w <- wine(dfint,nreps=NREPS)
obswine <- w$wine
print(paste0(file_name," wine ",obswine))
dfanid <- rbind(dfanid,data.frame("wine"=obswine,"exper"=0))


for (i in 1:numexper){
  ind_matrix <- read.table(experiment_files[i],sep="\t")
  ind_matrix_p <- MPack(ind_matrix)
  dfint <- as.data.frame(ind_matrix_p)
  w <- wine(dfint,nreps=NREPS)
  obswine <- w$wine
  print(paste0(experiment_files[i]," wine ",obswine))
  dfanid <- rbind(dfanid,data.frame("wine"=obswine,"exper"=i))
}

write.table(dfanid,paste0("../nestedness/",file_name,"_nestvalues.csv"),row.names = FALSE,
            sep=";")
