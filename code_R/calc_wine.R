library(grid)
library(gridExtra)
library(ggplot2)
library(kcorebip)

MPack <- function(matrix)
{
  sum_row <- rep(0,nrow(matrix))
  sum_col <- rep(0,ncol(matrix))
  for (i in 1:nrow(matrix))
    sum_row[i] <- sum(matrix[i,])
  for (i in 1:ncol(matrix))
    sum_col[i] <- sum(matrix[,i])
  ord_matrix <- matrix[rev(order(sum_row)),rev(order(sum_col))]
  return(t(ord_matrix))       # Transpose because of order of python-written matrix
}

file_name <- "RedAdyCom2014"
experiment_files <- Sys.glob(paste0("../results/",file_name,"_W_*.txt"))

zero_matrix <- read.table(experiment_files[1],sep="\t")
for (l in 1:nrow(zero_matrix))
  for(m in 1:ncol(zero_matrix))
    zero_matrix[l,m]<-0

numexper <- length(experiment_files)



emp_matrix <- read.table(paste0("../data/",file_name,".txt"),sep="\t")
dfanid <- data.frame("NDOF"=c(),"WNODF"=c(),"wine"=c(),"exper"=c())
ind_matrix_p <- MPack(emp_matrix)
# Remove all zeroes columns and rows
dfint <- as.data.frame(ind_matrix_p)
write.csv(dfint,"dfint.csv")
result_analysis <- analyze_network("dfint.csv", directory = "", guild_a = "Exp", guild_b = "Imp", plot_graphs = FALSE, only_NODF = FALSE)
obsnodf <- result_analysis$nested_values["NODF"]
obswnodf <- result_analysis$nested_values["weighted NODF"]
obswine <- result_analysis$nested_values["wine"]
print(paste0(file_name," WNODF ",obswnodf))
print(paste0(file_name," wine ",obswine))
dfanid <- rbind(dfanid,data.frame("NDOF"=obsnodf,"WNODF"=obswnodf,
                                  "wine"=obswine,"exper"=0))

for (i in 1:numexper){
  ind_matrix <- read.table(experiment_files[i],sep="\t")
  ind_matrix_p <- MPack(ind_matrix)
  dfint <- as.data.frame(ind_matrix_p)
  write.csv(dfint,"dfint.csv")
  result_analysis <- analyze_network("dfint.csv", directory = "", guild_a = "Exp", guild_b = "Imp", plot_graphs = FALSE, only_NODF = FALSE)
  obsnodf <- result_analysis$nested_values["NODF"]
  obswnodf <- result_analysis$nested_values["weighted NODF"]
  obswine <- result_analysis$nested_values["wine"]
  print(paste0(experiment_files[i]," WNODF ",obswnodf))
  print(paste0(experiment_files[i]," wine ",obswine))
  dfanid <- rbind(dfanid,data.frame("NDOF"=obsnodf,"WNODF"=obswnodf,
                                    "wine"=obswine,"exper"=i))
}

write.table(dfanid,paste0("../anidamientos/",file_name,"_anid.csv"),row.names = FALSE,
            sep=";")

mean_matrix <- zero_matrix/numexper



# 
# hm_emp <- crea_lista_heatmap(MPack(emp_matrix))
# hm_mean <- crea_lista_heatmap(MPack(mean_matrix))
# 
# maxleg <- (1+round(max(max(emp_matrix),max(mean_matrix)))%/%100)*100


# paint_int_matrix <- function(mq,titulo="",maximo=100)
# {
#   zp1 <- ggplot(mq,
#                 aes(x = X, y = rev(Y)))
#   b <- c(0,1,maximo/2)
#   zp1 <- zp1 + geom_tile(aes(fill=cuenta+0.00001)) + scale_fill_gradientn(colours=c("grey95",
#                                                                                     "blue","red"),
#                                                    trans = "log",  
#                                                    breaks =b, 
#                                                    labels=b)
#   zp1 <- zp1 + coord_equal() + ggtitle(titulo) + xlab("E-exporter") + ylab("I-importer")
#   
#   zp1 <- zp1 + theme_bw() +theme(panel.grid.major = element_blank(),
#                                  panel.grid.minor = element_blank(),
#                                  legend.position = "right",
#                                  axis.text = element_blank(),
#                                  legend.title = element_blank(),
#                                  axis.ticks = element_blank(),
#                                  panel.border = element_blank(),
#                                  plot.title = element_text(hjust = 0.5))
#   return(zp1)
# }

# 
# dir.create("../figuras", showWarnings = FALSE)
# fsal <- paste0("../figuras/",file_name,"_nexper_",numexper,"_IntMatrix.png")
# 
# ppi <- 600
# png(fsal, width=10*ppi, height=4*ppi, res=ppi)
# grid.arrange(m_emp, m_mean, ncol=2)
# dev.off()
# 

