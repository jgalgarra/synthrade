library(nortest)

crea_lista_heatmap <- function(matriz, justcount = FALSE)
{
  df <- data.frame("N"=c(),"cuenta"=c(),"type"=c())
  # Only sum 1 per filled cell to return degree instead of weight
  if (justcount)
    matriz[matriz>0] = 1
  for (l in 1:nrow(matriz))
  {
    dfaux <- data.frame("N"=l,"cuenta"=sum(matriz[l,]),"type"="EXP")
    df <- rbind(df,dfaux)
  }
  for(m in 1:ncol(matriz))
  {
    dfaux <- data.frame("N"=m,"cuenta"=sum(matriz[,m]),"type"="IMP")
    df <- rbind(df,dfaux)
  }
  return(df)
}

MSimp <- function(matrix,normalize = TRUE)
{
  sum_row <- rep(0,nrow(matrix))
  sum_col <- rep(0,ncol(matrix))
  if (normalize)
    matrix = matrix/mean(as.matrix(matrix))
  sum_row <- rowSums(matrix)
  sum_col <- colSums(matrix)
  return(t(matrix))       # Transpose because of order of python-written matrix
}

source("parse_command_line_args.R")

anyos <- seq(ini_seq,end_seq)

anyos <- seq(1962,2014)

dfbestlillies <- data.frame("Year"=c(),"Experiment"=c())

for (year in anyos){
  
  file_name <- paste0("RedAdyCom",year,"_FILT")
  file_orig <- paste0("RedAdyCom",year)
  filt_matrix <- read.table(paste0("../data/",file_name,".txt"),sep="\t")
  orig_matrix <- read.table(paste0("../data/",file_orig,".txt"),sep="\t")
  
  experiment_files <- Sys.glob(paste0("../results/",file_name,"_W_*.txt"))
  sim_matrix <- read.table(experiment_files[1],sep="\t")
  hm_filt <- crea_lista_heatmap(MSimp(filt_matrix,normalize = TRUE))
  hm_sim <- crea_lista_heatmap(MSimp(sim_matrix,normalize = TRUE))
  maxlilliescore <- 0
  poslilliescore <- 0
  for (i in seq(1,length(experiment_files))){
    other_sim_matrix <- read.table(experiment_files[i],sep="\t")
    other_hm_sim <- crea_lista_heatmap(MSimp(other_sim_matrix,normalize = FALSE))
    hm_sim_exp <- other_hm_sim[other_hm_sim$type == "EXP",]$cuenta
    hm_sim_imp <- other_hm_sim[other_hm_sim$type == "IMP",]$cuenta
    plillieexp <- lillie.test(log(hm_sim_exp))$p.value
    plillieimp <- lillie.test(log(hm_sim_imp))$p.value
    plillie <- sqrt(plillieexp*plillieimp)
    if (plillie > maxlilliescore){
      maxlilliescore <- plillie
      poslilliescore <- i
    }
  }
  hm_filt <- hm_filt[hm_filt$cuenta>0,]
  ll_filt_exp <- lillie.test(log(hm_filt[hm_filt$type == "EXP",]$cuenta))$p.value
  ll_filt_imp <- lillie.test(log(hm_filt[hm_filt$type == "IMP",]$cuenta))$p.value
  print(paste("Year",year,"Experiment",poslilliescore))
  dfbestlillies <- rbind(dfbestlillies,data.frame("Year"=year,"Experiment"=poslilliescore,"Geom_mean"=sqrt(plillie),
                                                  "Synthetic_exporter" = plillieexp,
                                                  "Synthetic_importer" = plillieimp,
                                                  "Empirical_exporter"=ll_filt_exp,
                                                  "Empirical_importer"=ll_filt_imp))
  write.table(dfbestlillies,"../results/BestLillies.txt",sep="\t",row.names = FALSE)

}

