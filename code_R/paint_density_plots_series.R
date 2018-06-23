library(grid)
library(gridExtra)
library(ggplot2)
source("aux_functions_matrix.R")
source("read_filter_condition.R")

MSimp <- function(matrix,normalize = TRUE)
{
  sum_row <- rep(0,nrow(matrix))
  sum_col <- rep(0,ncol(matrix))
  if (normalize)
    matrix = matrix/sum(matrix)
  sum_row <- rowSums(matrix)
  sum_col <- colSums(matrix)
  matrix <- matrix[sum_row>0,]
  matrix <- matrix[,sum_col>0]
  return(t(matrix))       # Transpose because of order of python-written matrix
}


PaintDensPlot <- function(datos,titletext,xlabel)
{
  p <- ggplot() + geom_density(aes(x= cuenta, color = collection, fill = collection),  alpha = .1,
                               data=datos, position = "identity", adjust=1)+ 
    xlab(xlabel)+ylab("Count\n")+
    ggtitle(titletext)+ scale_x_log10()+
    scale_fill_manual(values=c("blue","white","red"))+
    scale_color_manual(values=c("blue","grey","red"))+
    theme_bw() +
    theme(panel.border = element_blank(),
          legend.key = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(linetype = 2, color="ivory3"),
          panel.grid.major.x = element_blank(), 
          legend.title = element_blank(),
          legend.text = element_text(size=12, face="bold"),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(lineheight=.8, face="bold"),
          axis.text = element_text(face="bold", size=13),
          axis.title.x = element_text(face="bold", size=13),
          axis.title.y  = element_text(face="bold", size=13) )
  
  return(p)
}

PaintBoxPlot <- function(datos,titletext,xlabel)
{
  p <- ggplot() + geom_boxplot(aes(x=as.factor(collection),y= cuenta, color = collection, fill = collection),  alpha = .1,
                               data=datos, position = "identity")+ 
    xlab(xlabel)+ylab("Count\n")+
    ggtitle(titletext)+ scale_y_log10()+
    scale_fill_manual(values=c("blue","white","red"))+
    scale_color_manual(values=c("blue","grey","red"))+
    theme_bw() +
    theme(panel.border = element_blank(),
          legend.key = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(linetype = 2, color="ivory3"),
          panel.grid.major.x = element_blank(), 
          legend.title = element_blank(),
          legend.text = element_text(size=12, face="bold"),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(lineheight=.8, face="bold"),
          axis.text = element_text(face="bold", size=13),
          axis.title.x = element_text(face="bold", size=13),
          axis.title.y  = element_text(face="bold", size=13) )
  
  return(p)
}

source("parse_command_line_args.R")

anyos <- seq(ini_seq,end_seq)

sbestlillies <- FALSE       # If set to TRUE searches the best GOF in BestLillies.txt
                             # else chooses experiment number 1

if (nchar(filtered_string)>1){
  if (sbestlillies)
    bestlillies <- read.table("../results/BestLillies.txt",header=TRUE)
  fstring <- "FILT"
} else {
  if (sbestlillies)
    bestlillies <- read.table("../results/BestLilliesUnfiltered.txt",header=TRUE)
  fstring <- "UNFILT"
}

for (year in anyos){
  if (sbestlillies)
    posbest <- bestlillies[bestlillies$Year==year,]$Experiment
  else
    posbest <- 1
  file_name <- paste0("RedAdyCom",year,filtered_string)
  file_orig <- paste0("RedAdyCom",year)
  experiment_files <- Sys.glob(paste0("../results/",file_name,"_W_",posbest,".txt"))
  numexper <- 1
  filt_matrix <- read.table(paste0("../data/",file_name,".txt"),sep="\t")
  orig_matrix <- read.table(paste0("../data/",file_orig,".txt"),sep="\t")
  sim_matrix <- read.table(experiment_files[1],sep="\t")
  hm_filt <- crea_lista_heatmap(MSimp(filt_matrix,normalize = FALSE),justcount = TRUE)
  hm_sim <- crea_lista_heatmap(MSimp(sim_matrix,normalize = FALSE),justcount = TRUE)
  hm_orig <- crea_lista_heatmap(MSimp(orig_matrix,normalize = FALSE),justcount = TRUE)
  hm_filt$collection <- "Filtered"
  hm_sim$collection <- "Synthetic"
  hm_orig$collection <- "Original"
  hm_all_deg <- rbind(hm_filt,hm_sim,hm_orig)
  hm_all_importers_deg <- hm_all_deg[hm_all_deg$type=="IMP",]
  hm_filt <- crea_lista_heatmap(MSimp(filt_matrix,normalize = TRUE))
  hm_sim <- crea_lista_heatmap(MSimp(sim_matrix,normalize = TRUE))
  hm_orig <- crea_lista_heatmap(MSimp(orig_matrix,normalize = TRUE))
  hm_filt$collection <- "Filtered"
  hm_sim$collection <- "Synthetic"
  hm_orig$collection <- "Original"
  hm_all_weight <- rbind(hm_filt,hm_sim,hm_orig)
  hm_all_importers_weight <- hm_all_weight[hm_all_weight$type=="IMP",]
  
  hm_all_exporters_deg <- hm_all_deg[hm_all_deg$type=="EXP",]
  hm_filt <- crea_lista_heatmap(MSimp(filt_matrix,normalize = TRUE))
  hm_sim <- crea_lista_heatmap(MSimp(sim_matrix,normalize = TRUE))
  hm_orig <- crea_lista_heatmap(MSimp(orig_matrix,normalize = TRUE))
  hm_filt$collection <- "Filtered"
  hm_sim$collection <- "Synthetic"
  hm_orig$collection <- "Original"
  hm_all_weight <- rbind(hm_filt,hm_sim,hm_orig)
  hm_all_exporters_weight <- hm_all_weight[hm_all_weight$type=="EXP",]
 
  q <- PaintDensPlot(hm_all_importers_deg,"Importers","Degree")
  r <- PaintDensPlot(hm_all_importers_weight,"Importers","Normalized strength")
  s <- PaintDensPlot(hm_all_exporters_deg,"Exporters","Degree")
  t <- PaintDensPlot(hm_all_exporters_weight,"Exporters","Normalized strength")
  
  bq <- PaintBoxPlot(hm_all_importers_deg,"Importers","Degree")
  br <- PaintBoxPlot(hm_all_importers_weight,"Importers","Normalized strength")
  bs <- PaintBoxPlot(hm_all_exporters_deg,"Exporters","Degree")
  bt <- PaintBoxPlot(hm_all_exporters_weight,"Exporters","Normalized strength")
  dir.create("../figures/densities/", showWarnings = FALSE)
  fsal <- paste0("../figures/densities/Density_DegStr_",year,"_",fstring,".png")
  ppi <- 600
  png(fsal, width=12*ppi, height=6*ppi, res=ppi)
  grid.arrange(q,r,s,t, ncol=2, nrow=2,top=year )
  dev.off()
  fsal2 <- paste0("../figures/densities/Boxplot_DegStr_",year,"_",fstring,".png")
  ppi <- 600
  png(fsal2, width=12*ppi, height=6*ppi, res=ppi)
  grid.arrange(bq,br,bs,bt, ncol=2, nrow=2,top=year )
  dev.off()

  
}

