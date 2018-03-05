library(grid)
library(gridExtra)
library(ggplot2)

crea_lista_heatmap <- function(matriz)
{
  df <- data.frame("X"=c(),"Y"=c(),"cuenta"=c())
  for (l in 1:nrow(matriz))
    for(m in 1:ncol(matriz))
      if (matriz[l,m] >0 )
      {
        dfaux <- data.frame("X"=l,"Y"=m,"cuenta"=matriz[l,m])
        df <- rbind(df,dfaux)
      }
  return(df)
}

MSimp <- function(matrix,normalize = TRUE)
{
  sum_row <- rep(0,nrow(matrix))
  sum_col <- rep(0,ncol(matrix))
  if (normalize)
    matrix = matrix/max(matrix)
  for (i in 1:nrow(matrix))
    sum_row[i] <- sum(matrix[i,])
  for (i in 1:ncol(matrix))
    sum_col[i] <- sum(matrix[,i])
  matrix <- matrix[sum_row>0,sum_col>0]
  return(t(matrix))       # Transpose because of order of python-written matrix
}


anyos <- seq(1998,2014)
for (year in anyos){
  
  file_name <- paste0("RedAdyCom",year,"_FILT")
  file_orig <- paste0("RedAdyCom",year)
  experiment_files <- Sys.glob(paste0("../results/",file_name,"_W_1.txt"))
  numexper <- 1
  filt_matrix <- read.table(paste0("../data/",file_name,".txt"),sep="\t")
  orig_matrix <- read.table(paste0("../data/",file_orig,".txt"),sep="\t")
  sim_matrix <- read.table(experiment_files[1],sep="\t")
  
  hm_filt <- crea_lista_heatmap(MSimp(filt_matrix,normalize = TRUE))
  hm_sim <- crea_lista_heatmap(MSimp(sim_matrix,normalize = TRUE))
  hm_orig <- crea_lista_heatmap(MSimp(orig_matrix,normalize = TRUE))
  
  
  hist(log10(hm_filt$cuenta),breaks=20)
  hist(log10(hm_sim$cuenta),breaks=20)
  hist(log10(hm_orig$cuenta),breaks=20)
  
  hm_filt$collection <- "Filtered"
  hm_sim$collection <- "Synthetic"
  hm_orig$collection <- "Original"
  
  hm_all <- rbind(hm_filt,hm_sim,hm_orig)
  
  
  
  p <- ggplot() + geom_histogram(aes(x= cuenta, color = collection, fill = collection),  alpha = .1,
                               data=hm_all, position = "identity",bins = 50)+ xlab("Normalized strength")+ylab("Count\n")+
    ggtitle(year)+ scale_x_log10()+scale_fill_manual(values=c("blue","white","red"))+
    scale_color_manual(values=c("blue","grey","red"))+
    theme_bw() +
    theme(panel.border = element_blank(),
          legend.key = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(linetype = 2, color="ivory3"),
          panel.grid.major.x = element_blank(), 
          legend.title = element_text(size=14, face="bold"),
          legend.text = element_text(size=8, face="bold"),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(lineheight=.8, face="bold"),
          axis.text = element_text(face="bold", size=13),
          axis.title.x = element_text(face="bold", size=13),
          axis.title.y  = element_text(face="bold", size=13) )
  
  
  

  dir.create("../figures/densities/", showWarnings = FALSE)
  fsal <- paste0("../figures/densities/Density_",year,".png")

  ppi <- 600
  png(fsal, width=8*ppi, height=4*ppi, res=ppi)
  print(p)
  dev.off()
}

