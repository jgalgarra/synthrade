# Synthetic trade network creation

ReadMatrix <- function(year){
  original_file <- read.delim(paste0("../data/RedAdyCom",year,"_FILT.txt"), header=FALSE)
  #original_file <- read.delim(paste0("../data/RedAdyCom",year,".txt"), header=FALSE)
  or_matrix <- as.matrix(original_file)
  # Clean rows and cols full of zeroes
  clean_matrix <- or_matrix[,colSums(or_matrix) > 0]
  clean_matrix <- clean_matrix[rowSums(clean_matrix) > 0,]
  return(clean_matrix)
}

PrefAttachment <- function(vecprob,lvec)
{
  listanodes = c()
  for (i in 1:lvec){
    if (vecprob[i] != 0)
      if (rbinom(1,1,vecprob[i])>0)
        listanodes <- append(listanodes,i)
  }
  #node = sample(seq(1,lvec),1,vecprob[1:lvec],replace=FALSE)
  return(listanodes)
}

UpdatableLinks <- function(matrixprob){  
  links_found = FALSE
  listaedges = c()
  tam = c(nrow(matrixprob),ncol(matrixprob))
  # mpos <- matrix(seq(1,tam[1]*tam[2]),nrow=tam[1],ncol=tam[2])
  # positions <- which(sample(mpos,1,prob=matrixprob)==mpos,arr.ind=T)
  newlinks <- 0
  while (sum(newlinks) == 0){
    randunif <- runif(tam[1]*tam[2],0,1)
    randunif = matrix(randunif,nrow=tam[1],ncol=tam[2])
    newlinks = randunif < matrixprob
    positions <- which(newlinks !=0, arr.ind = T)
  }
  return(positions)
}

SynthMatrix <- function(matrixemp, year){
  n_imp <- ncol(matrixemp)
  n_exp <- nrow(matrixemp)
  numlinks <- sum(matrixemp > 0)
  totweight <- sum(matrixemp)
  print(paste("exporters",n_exp,"importers",n_imp,"numlinks",numlinks))
  
  # Create a synthetic matrix full of zeroes
  msynth <- matrix(rep(0.0,n_imp*n_exp), nrow = n_exp, byrow = TRUE)
  exp_max <- 3
  imp_max <- 3
  min_token <- 1
  msynth[1,2] <- min_token 
  msynth[1,3] <- min_token
  msynth[2,1] <- min_token
  msynth[2,3] <- min_token
  msynth[3,1] <- min_token
  msynth[3,2] <- min_token
  lambda_imp = (n_imp^2-n_imp)/(2*numlinks)
  lambda_exp = (n_exp^2-n_exp)/(2*numlinks)
  
  cuenta_links <- sum(msynth > 0)
  min_links <- cuenta_links
  print(paste("lambda imp:",lambda_imp,"lambda_exp",lambda_exp))
  Pr_E <- rowSums(msynth)/sum(msynth)
  Pr_I <- colSums(msynth)/sum(msynth)
  void_prob <- matrix(rep(0.0,exp_max*imp_max), nrow = imp_max, byrow = TRUE)
  prob_new_links <- void_prob
  morenewnodes <- TRUE
  cuenta_antciclo <- 0
  
  while ((morenewnodes)|| (cuenta_links < numlinks))
  {
    new_node <- FALSE
    if (cuenta_antciclo != cuenta_links){
      cuenta_antciclo <- cuenta_links
      if (cuenta_links %% 1000 == 0) 
        print(paste(cuenta_links,"links out of",numlinks,
                    "exporters",exp_max,"out of",n_exp,"importers",
                    imp_max,"out of",n_imp))
    }
    if (exp_max < n_exp)
      if (rbinom(1,1,min(1,lambda_exp/exp_max))>0)
      {
        linkstoI <- c()
        while (length(linkstoI) == 0)
          linkstoI <- PrefAttachment(Pr_I,imp_max)
        for (i in linkstoI)
          if ((cuenta_links < numlinks) && (exp_max < n_exp)){
            exp_max <- exp_max + 1
            msynth[exp_max,i] <- 1#/length(linkstoI)
            cuenta_links <-  cuenta_links + 1
            new_node <- TRUE
          }
      }
    
    if (imp_max < n_imp)
      if (rbinom(1,1,min(1,lambda_imp/imp_max))>0)
      {
        linkstoE <- c()
        while (length(linkstoE) == 0)
          linkstoE <- PrefAttachment(Pr_E,exp_max)
        for (i in linkstoE)
          if ((cuenta_links < numlinks) && (imp_max < n_imp)){
            imp_max <- imp_max + 1
            msynth[i,imp_max] <- 1#/length(linkstoE)
            cuenta_links <-  cuenta_links + 1
            new_node <- TRUE
          }
      }
    
    if (new_node){
      smat <- sum(msynth)
      Pr_E <- rowSums(msynth)/smat
      Pr_I <- colSums(msynth)/smat
      prob_new_links <- t(Pr_I[1:imp_max] %o% Pr_E[1:exp_max])
    }
    else if (cuenta_links > min_links){
      update_links <- UpdatableLinks(prob_new_links)
      lupdate <- (length(update_links)/2)
      for (m in 1:lupdate)
        if (cuenta_links < numlinks) {
          rowl <- update_links[m,1]
          coll <- update_links[m,2]
          msynth[rowl,coll] <- msynth[rowl,coll] + 1#/lupdate
        }
    }
    
    cuenta_links <- sum(msynth > 0)
    smat <- sum(msynth)
    Pr_E <- rowSums(msynth)/smat
    Pr_I <- colSums(msynth)/smat
    prob_new_links <- t(Pr_I[1:imp_max] %o% Pr_E[1:exp_max])
    morenewnodes <- (exp_max < n_exp) || (imp_max < n_imp)
  }
  return(msynth)
}

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  ini_seq <- 1962
  end_seq <- 1962
  maxexper <- 1
} else{
  ini_seq <- as.numeric(args[1])
  end_seq <- as.numeric(args[2])
  maxexper <- as.numeric(args[3]) 
}

years <- seq(ini_seq,end_seq)
for (lyear in years)
  for (nexper in seq(1,maxexper)){
    print(paste(lyear,"Experiment",nexper))
    matrix_emp <- ReadMatrix(lyear)
    nlinks <- sum(matrix_emp>0)
    matrix_experiment <- SynthMatrix(matrix_emp,lyear)
    dir.create("../results", showWarnings = FALSE)
    write.table(matrix_experiment,paste0("../results/RedAdyCom",lyear,"_FILT_W_",nexper,".txt"),
                row.names = FALSE, col.names = FALSE, sep = "\t")
  }