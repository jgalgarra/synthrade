library(nortest)
source("read_filter_condition.R")
if (nchar(filtered_string)>1) {
  fcond <- "YES"
} else
  fcond <- "NO"

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


dfbestlillies <- data.frame("Year"=c(),"Experiment"=c())

for (year in anyos){
  
  unfilt_name <- paste0("RedAdyCom",year)
  file_name <- paste0("RedAdyCom",year,filtered_string)
  file_orig <- paste0("RedAdyCom",year)
  experiment_files <- Sys.glob(paste0("../results/",file_name,"_W_*.txt"))
  if (length(experiment_files)>0){
    
    filt_matrix <- read.table(paste0("../data/",file_name,".txt"),sep="\t")
    unfilt_matrix <- read.table(paste0("../data/",unfilt_name,".txt"),sep="\t")
    orig_matrix <- read.table(paste0("../data/",file_orig,".txt"),sep="\t")
    
    
    sim_matrix <- read.table(experiment_files[1],sep="\t")
    hm_unfilt <- crea_lista_heatmap(MSimp(unfilt_matrix,normalize = TRUE))
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
    
    hm_unfilt <- hm_unfilt[hm_unfilt$cuenta>0,]
    ll_unfilt_exp <- lillie.test(log(hm_unfilt[hm_unfilt$type == "EXP",]$cuenta))$p.value
    ll_unfilt_imp <- lillie.test(log(hm_unfilt[hm_unfilt$type == "IMP",]$cuenta))$p.value
    print(paste("Year",year,"Experiment",poslilliescore))
    dfbestlillies <- rbind(dfbestlillies,data.frame("Year"=year,"Experiment"=poslilliescore,"Geom_mean"=sqrt(plillie),
                                                    "Synthetic_exporter" = plillieexp,
                                                    "Synthetic_importer" = plillieimp,
                                                    "Empirical_exporter_filtered"=ll_filt_exp,
                                                    "Empirical_importer_filtered"=ll_filt_imp,
                                                    "Empirical_exporter_unfiltered"=ll_unfilt_exp,
                                                    "Empirical_importer_unfiltered"=ll_unfilt_imp))
    if (fcond == "YES")
      write.table(dfbestlillies,"../results/BestLillies.txt",sep="\t",row.names = FALSE)
    else
      write.table(dfbestlillies,"../results/BestLilliesUnfiltered.txt",sep="\t",row.names = FALSE)

  }
}

