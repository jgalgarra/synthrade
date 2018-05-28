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

MPack <- function(matrix,normalize = TRUE)
{
  sum_row <- rep(0,nrow(matrix))
  sum_col <- rep(0,ncol(matrix))
  if (normalize)
    matrix = matrix/sum(matrix)
  for (i in 1:nrow(matrix))
    sum_row[i] <- sum(matrix[i,])
  for (i in 1:ncol(matrix))
    sum_col[i] <- sum(matrix[,i])
  ord_matrix <- matrix[rev(order(sum_row)),rev(order(sum_col))]
  return(t(ord_matrix))       # Transpose because of order of python-written matrix
}

