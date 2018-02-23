library(grid)
library(gridExtra)
library(ggplot2)

file_name <- "RedAdyCom2006_FILT"
experiment_files <- Sys.glob(paste0("../results/",file_name,"_W_*.txt"))

zero_matrix <- read.table(experiment_files[1],sep="\t")
for (l in 1:nrow(zero_matrix))
  for(m in 1:ncol(zero_matrix))
    zero_matrix[l,m]<-0

numexper <- length(experiment_files)

numexper <- 1

crea_lista_heatmap <- function(matriz)
{
df <- data.frame("X"=c(),"Y"=c(),"cuenta"=c())
#im <- as.data.frame(matrix(rep(0,3*nrow(matriz)*ncol(matriz)),nrow(matriz)*ncol(matriz),3))
#names(im) <- c("X","Y","cuenta")
# im$Y <- rep(seq(1,ncol(matriz)),nrow(matriz))
# for (i in 1:nrow(im)){
#   rowindex <- 1 + ((i-1) %/% nrow(matriz))
#   im$X[i] <- rowindex
#   im$cuenta[i] <- matriz[rowindex,im$Y[i]]
# }
for (l in 1:nrow(matriz))
  for(m in 1:ncol(matriz)){
    dfaux <- data.frame("X"=l,"Y"=m,"cuenta"=matriz[l,m])
    df <- rbind(df,dfaux)
  }
return(df)
}

for (i in 1:numexper){
  ind_matrix <- read.table(experiment_files[i],sep="\t")
  zero_matrix <- zero_matrix + ind_matrix
}

mean_matrix <- zero_matrix/numexper
emp_matrix <- read.table(paste0("../data/",file_name,".txt"),sep="\t")
# eliminar filas y columnas a cero
sum_row <- rep(0,nrow(emp_matrix))
sum_col <- rep(0,ncol(emp_matrix))
for (i in 1:nrow(emp_matrix))
  sum_row[i] <- sum(emp_matrix[i,])
for (i in 1:ncol(emp_matrix))
  sum_col[i] <- sum(emp_matrix[,i])
emp_matrix <- emp_matrix[sum_row>0,sum_col>0]

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

hm_emp <- crea_lista_heatmap(MPack(emp_matrix))
hm_mean <- crea_lista_heatmap(MPack(mean_matrix))

maxleg <- (1+round(max(max(emp_matrix),max(mean_matrix)))%/%100)*100


paint_int_matrix <- function(mq,titulo="",maximo=100)
{
  zp1 <- ggplot(mq,
                aes(x = X, y = rev(Y)))
  b <- c(0,1,maximo/2)
  zp1 <- zp1 + geom_tile(aes(fill=cuenta+0.00000001)) + scale_fill_gradientn(colours=c("grey95",
                                                                                    "blue","red"),
                                                   trans = "log",  
                                                   breaks =b, 
                                                   labels=b)
  zp1 <- zp1 + coord_equal() + ggtitle(titulo) + xlab("E-exporter") + ylab("I-importer")
  
  zp1 <- zp1 + theme_bw() +theme(panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 legend.position = "right",
                                 axis.text = element_blank(),
                                 legend.title = element_blank(),
                                 axis.ticks = element_blank(),
                                 panel.border = element_blank(),
                                 plot.title = element_text(hjust = 0.5))
  return(zp1)
}

hist(log10(hm_emp$cuenta),breaks=10)
hist(log10(hm_mean$cuenta),breaks=10)

m_emp <- paint_int_matrix(hm_emp,titulo=paste(file_name,"Empirical Matrix"))
m_mean <- paint_int_matrix(hm_mean,titulo=paste0("Simulated Matrix. #Experiments: ",numexper ))

dir.create("../figures", showWarnings = FALSE)
fsal <- paste0("../figures/",file_name,"_nexper_",numexper,"_IntMatrix.png")

ppi <- 600
png(fsal, width=10*ppi, height=4*ppi, res=ppi)
grid.arrange(m_emp, m_mean, ncol=2)
dev.off()


