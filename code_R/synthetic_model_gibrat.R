# Synthetic trade network creation

lyear <- 1950
nexper <- 1

ReadMatrix <- function(year){
  original_file <- read.delim(paste0("../data/RedAdyCom",year,"_FILT.txt"), header=FALSE)
  original_matrix <- as.matrix(original_file)
  # Clean rows and cols full of zeroes
  clean_matrix <- original_matrix[,colSums(original_matrix) > 0]
  clean_matrix <- clean_matrix[rowSums(clean_matrix) > 0,]
  return(clean_matrix)
}

PrefAttachment <- function(vecprob,lvec,newnode = FALSE)
{
  if (newnode)
    listanodes = sample(seq(1,lvec),1,vecprob[1:lvec],replace=TRUE)
  else{
  listanodes = c()
    for (i in 1:lvec)
    {
        if (vecprob[i] != 0)
          if (rbinom(1,1,vecprob[i])>0)
            listanodes <- append(listanodes,i)
    }
  }
return(listanodes)
}

UpdatableLinks <- function(matrixprob){  
  
  links_found = FALSE
  tam = c(nrow(matrixprob),ncol(matrixprob))
  while (!links_found){
    randunif <- runif(tam[1]*tam[2],0,1)
    newlinks <- c()
    while (length(newlinks) == 0){
      randunif = matrix(randunif,nrow=tam[1],ncol=tam[2])
      newlinks = randunif < matrixprob  
      positions <- which(newlinks !=0, arr.ind = T)
      links_found = length(positions>0)
    }
  }
  return(positions)
}

NewLinks <- function(m,nrows,ncols)
{
  raw_mat <- m[1:nrows,1:ncols]
  links_found <- FALSE
  nvoids <- sum(raw_mat==0)
  pr_mat <- raw_mat
  pr_mat[raw_mat > 0] <- 0
  pr_mat[raw_mat == 0] <- 1/nvoids
  pos_voids <- which(pr_mat==sample(pr_mat[pr_mat>0],1),arr.ind = T)
  one_at_random <- sample(seq(1:nrow(pos_voids)),1)
  positions <- pos_voids[one_at_random,]
  return(positions)
}

IncLinks <- function(m,nrows,ncols)
{
  raw_mat <- m[1:nrows,1:ncols]
  links_found <- FALSE
  nincs <- sum(raw_mat>0)
  pr_mat <- raw_mat
  pr_mat[raw_mat == 0] <- 0
  pr_mat<- pr_mat/nincs
  pos_voids <- which(pr_mat==sample(pr_mat[pr_mat>0],1),arr.ind = T)
  one_at_random <- sample(seq(1:nrow(pos_voids)),1)
  positions <- pos_voids[one_at_random,]
  return(positions)
}

SynthMatrix <- function(matrixemp, year){
  n_imp <- ncol(matrixemp)
  n_exp <- nrow(matrixemp)
  print(paste("exporters",n_exp,"importers",n_imp))
  numlinks <- sum(matrixemp > 0)
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
  cuenta_links <- sum(msynth >0)
  print(paste("lambda imp:",lambda_imp,"lambda_exp",lambda_exp))
  Pr_E <- rowSums(msynth)/sum(msynth)
  Pr_I <- colSums(msynth)/sum(msynth)
  void_prob <- matrix(rep(0.0,exp_max*imp_max), nrow = imp_max, byrow = TRUE)
  prob_new_links <- void_prob
  morenewnodes <- TRUE
  cuenta_antciclo <- 0
  while ((morenewnodes) || (cuenta_links < numlinks))
  {
    new_node <- FALSE
    if (cuenta_antciclo != cuenta_links){
      cuenta_antciclo <- cuenta_links
      if (cuenta_links %% 100 == 0) 
        print(paste(cuenta_links,"links out of",numlinks))
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
            msynth[exp_max,i] <- 1
            cuenta_links <-  cuenta_links + 1
            new_node <- TRUE
            print(paste("cuenta_links",cuenta_links," out of",numlinks,"exp_max",exp_max))
          }
          # else
          #   break
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
            msynth[i,imp_max] <- 1
            cuenta_links <-  cuenta_links + 1
            new_node <- TRUE
            print(paste("cuenta_links",cuenta_links," out of",numlinks,"imp_max",imp_max))
          }
          # else
          #   break
      }
    if (new_node){
        Pr_E <- rowSums(msynth)/sum(msynth)
        Pr_I <- colSums(msynth)/sum(msynth)
        prob_new_links <- matrix(apply(expand.grid(Pr_I[1:imp_max],Pr_E[1:exp_max]), 1, 
                                       FUN = function(x) {x[1] * x[2]}),
                                       nrow = exp_max, ncol = imp_max, byrow = TRUE)
    }

    
    
    nlink <- NewLinks(msynth,exp_max,imp_max)
    msynth[nlink["row"],nlink["col"]] <- 1
    
    inclink <- IncLinks(msynth,exp_max,imp_max)
    msynth[inclink["row"],inclink["col"]] <- 1
    # update_links <- UpdatableLinks(prob_new_links)
    # for (m in 1:(length(update_links)/2)){
    #   if (msynth[update_links[m,1],update_links[m,2]] == 0)
    #     if ((cuenta_links < numlinks) && (length(update_links) > 0)){
    #       msynth[update_links[m,1],update_links[m,2]] <- 1 + msynth[update_links[m,1],update_links[m,2]] 
    #       b<- 1-prob_new_links[update_links[m,1],update_links[m,2]]
    #       print("antes")
    #       print(update_links)
    #     }
    # }
    cuenta_links <- sum(msynth > 0)
    
    Pr_E <- rowSums(msynth)/sum(msynth)
    Pr_I <- colSums(msynth)/sum(msynth)
    morenewnodes <- (exp_max < n_exp) || (imp_max < n_imp)
  }
  return(msynth)
}

matrix_emp <- ReadMatrix(lyear)
matrix_experiment <- SynthMatrix(matrix_emp,lyear)
write.table(t(matrix_experiment),paste0("../results/RedAdyCom",lyear,"_FILT_W_",nexper,".txt"),
            row.names = FALSE, col.names = FALSE, sep = "\t")