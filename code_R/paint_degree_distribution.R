library(grid)
library(gridExtra)
library(igraph)
library(ggplot2)


lread_network <- function(namenetwork, guild_astr = "pl", guild_bstr = "pol", directory="")
{
  # Reading species names
  namesred <- read.csv(paste0(directory,namenetwork),header=FALSE,stringsAsFactors=FALSE)
  names_guild_a <- namesred[1,2:ncol(namesred)]
  names_guild_b <- namesred[2:nrow(namesred),1]
  
  #Reading matrix data
  m <- read.csv(paste0(directory,namenetwork),header=TRUE,row.names=1)
  
  # Calc number of species of each guild
  num_guild_a <- ncol(m)
  num_guild_b <- nrow(m)
  # Create an graph object
  g <- graph.empty()
  # Add one node for each species and name it
  for (i in 1:num_guild_a){
    g <- g + vertices(paste0(guild_astr,i),color="white",guild_id="a",name_species=names_guild_a[i],id=i)
  }
  for (i in 1:num_guild_b){
    g <- g + vertices(paste0(guild_bstr,i),color="red",guild_id="b",name_species=names_guild_b[i],id=i)
  }
  
  # Adding links to the graph object
  mm <- matrix(unlist(list(m)),nrow=num_guild_b,ncol=num_guild_a)
  listedgesn <- which(mm!=0, arr.ind = T)
  listedgesn <- listedgesn[order(listedgesn[,1],listedgesn[,2]),]
  listedgesn[,1] <- paste0(guild_bstr,listedgesn[,1])
  listedgesn[,2] <- paste0(guild_astr,listedgesn[,2])
  g <- g + graph.edgelist(listedgesn)
  # Return values
  calc_values <- list("graph" = g, "matrix" = m, "num_guild_b" = num_guild_b, "num_guild_a" = num_guild_a,
                      "names_guild_a" = names_guild_a, "names_guild_b"=names_guild_b)
  return(calc_values)
  
}

gen_deg_distribution <- function(red,series, colors, seq_breaks = c(1,5,10,20,30), TFS = "")
{
  
  gen_deg_data_frame <- function(input_matrix,tipo,tamanyo,nalpha,serie, TFS = TFS)
  {
  
  # Remove all zeroes columns and rows
  dfint <- as.data.frame(input_matrix)
  write.csv(dfint,"dfcab.csv")
  dr <- lread_network("dfcab.csv", guild_astr = "Exporter", guild_bstr = "Importer", directory="")
  grafo <- as.undirected(dr[["graph"]])
  ddegree <- igraph::degree(grafo,mode = c("out"), loops = TRUE, normalized = FALSE)
  dfdeg <- data.frame("degree" = as.numeric(ddegree))
  dfdeg$type <- "Exporter"
  dfdeg$tamanyo <- tamanyo
  dfdeg$nalpha <- nalpha
  dfdeg$nodename <- names(ddegree)
  dfdeg[grepl("Importer",dfdeg$nodename),]$type <- as.character("Importer")
  dfdeg <- dfdeg[order(dfdeg$degree),]
  for (k in 1:nrow(dfdeg)){
    if (dfdeg$type[k] == "Importer")
    {
      indice <- as.numeric(strsplit(dfdeg$nodename[k],"Importer")[[1]][2])
      dfdeg$weight[k] <- sum(dr$m[indice,])
    }
    else
    {
      indice <- as.numeric(strsplit(dfdeg$nodename[k],"Exporter")[[1]][2])
      dfdeg$weight[k] <- sum(dr$m[,indice])
    }
  }
  dfdeg$weight <- dfdeg$weight/sum(dfdeg$weight)
  ddeg_exporter <- dfdeg[dfdeg$type == "Exporter",]
  occur <- ddeg_exporter$degree

  woccur <- ddeg_exporter$weight
  woccur <- woccur[order(woccur)]

  alpha_level = 0.25
  p = occur/sum(occur)
  dy = rev(cumsum(rev(p)))
  dx = occur
  wp = woccur/sum(woccur)
  wdy = rev(cumsum(rev(wp)))
  wdx = woccur
  type = ddeg_exporter$type
  tamanyo = ddeg_exporter$tamanyo
  nalpha = ddeg_exporter$nalpha
  auxdf_exporter <- data.frame(dx,dy,type,tamanyo,nalpha)
  auxdfw_exporter <- data.frame(wdx,wdy,type,tamanyo,nalpha)

  ddeg_importer <- dfdeg[dfdeg$type == "Importer",]
  occur <- ddeg_importer$degree

  woccur <- ddeg_importer$weight

  woccur <- woccur[order(woccur)]
  alpha_level = 0.5
  p = occur/sum(occur)
  dy = rev(cumsum(rev(p)))
  dx = occur
  wp = woccur/sum(woccur)
  wdy = rev(cumsum(rev(wp)))
  wdx = woccur
  type = ddeg_importer$type
  tamanyo = ddeg_importer$tamanyo
  nalpha = ddeg_importer$nalpha
  auxdf_importer <- data.frame(dx,dy,type,tamanyo,nalpha)
  auxdfw_importer <- data.frame(wdx,wdy,type,tamanyo,nalpha)
  if (serie == "Exporter"){
    auxdf <- auxdf_exporter
    auxdfw <- auxdfw_exporter
  }
  if (serie == "Importer"){
    auxdf <- auxdf_importer
    auxdfw <- auxdfw_importer
  }
  if (serie == "Both"){
    auxdf <- rbind(auxdf_exporter, auxdf_importer)
    auxdfw <- rbind(auxdfw_exporter, auxdfw_importer)
  }
  auxdf$method <- tipo
  auxdfw$method <- tipo
  calc_values <- list("auxdf" = auxdf, "auxdfw" = auxdfw)
  return(calc_values)
  }
  dred <- gsub(TFstring,"",red)
  emp_matrix <- read.table(paste0("../data/",dred,".txt"),sep="\t")
  auxdf_emp <- gen_deg_data_frame(emp_matrix,"Empirical",1,0.2,series)
  auxdf <- auxdf_emp[["auxdf"]]
  auxdfw <- auxdf_emp[["auxdfw"]]
  if (TFS == "")
    subdir <- ""
  else
    subdir <- "TFMatrix"
  ficheros <- Sys.glob(paste0("../results/",subdir,red,"_W_*",".txt"))
  for (j in ficheros){
  #for (j in ficheros[2]){           # Solo para pruebas
    sim_matrix <- read.table(j,sep="\t")
    auxdf_sim <- gen_deg_data_frame(sim_matrix,"Simulated",0.5,0.02,series)
    auxdf <- rbind(auxdf,auxdf_sim[["auxdf"]])
    auxdfw <- rbind(auxdfw,auxdf_sim[["auxdfw"]])
  }

  auxdf <- auxdf[auxdf$dx > 0,]
  dist_deg <- ggplot(data = auxdf, aes(x = dx, y = dy)) + 
    geom_point(aes(alpha = nalpha,shape=method,size=tamanyo,stroke=tamanyo),color=colors) +
    scale_x_log10(breaks = seq_breaks) + scale_y_log10(breaks=c(0.1,0.2,0.5,1.0)) + xlab("Degree") + 
    ylab(cumulativetxt) + ggtitle("") +  scale_shape_manual(values=c(21, 15)) +
    scale_alpha(guide = 'none') +  scale_size_identity() + ggtitle(series) +
    theme_bw() +
    theme(
      axis.title.x = element_text(color="grey30", size=15),
      axis.title.y = element_text(color="grey30", size=15),
      legend.title=element_blank(),
      legend.text=element_text(size=14),
      axis.text.x = element_text(face="bold", color="grey30", size=13),
      axis.text.y = element_text(face="bold", color="grey30", size=13)
    )
  
  auxdfw <- auxdfw[auxdfw$wdx > 0,]
  dist_wdeg <- ggplot(data = auxdfw, aes(x = wdx, y = wdy)) + 
    geom_point(aes(alpha = nalpha,shape=method,size=tamanyo,stroke=tamanyo),color=colors) +
    scale_x_log10(breaks=c(0.001,0.01,0.1,1)) + scale_y_log10(breaks=c(0.1,0.2,0.5,1.0)) + xlab("Normalized strength") + 
    ylab(cumulativetxt) + ggtitle("") +  scale_shape_manual(values=c(21, 15)) +
    scale_alpha(guide = 'none') +  scale_size_identity() + ggtitle(series) +
    theme_bw() +
    theme(
      axis.title.x = element_text(color="grey30", size=15),
      axis.title.y = element_text(color="grey30", size=15),
      legend.title=element_blank(),
      legend.text=element_text(size=14),
      axis.text.x = element_text(face="bold", color="grey30", size=13),
      axis.text.y = element_text(face="bold", color="grey30", size=13)
    )
  calc_values <- list("dist_deg" = dist_deg, "dist_wdeg" = dist_wdeg)
  return(calc_values)
}


languageEl <- "EN"
if (languageEl == "EN"){
  cumulativetxt = "Cumulative distribution probability"
  xscale = "degree scale"
} else {
  cumulativetxt = "Distribuci?n acumulada de probabilidad"
  xscale = "escala degree"
}
source("parse_command_line_args.R")

# Third command argument allows to plot densities at build up time
TFstring <- as.character(args[3])
if (is.na(TFstring)){
  TFstring <- ""
} else
  TFstring <- "TF_"

# ini_seq <- 2000
# end_seq <- 2000

files <- paste0(TFstring,"RedAdyCom",seq(ini_seq,end_seq))
for (orig_file in files)
{
  red <- paste0(orig_file,"_FILT")
  series = "Exporter"
  grafs <- gen_deg_distribution(paste0(red),series,"blue")
  e_degree <- grafs$dist_deg
  e_weight <- grafs$dist_wdeg
  series = "Importer"
  grafs <- gen_deg_distribution(paste0(red),series,"red",TFS = TFstring)
  i_degree <- grafs$dist_deg
  i_weight <- grafs$dist_wdeg
  titulo=strsplit(red,"RedAdyCom")[[1]][-1]
  title1=textGrob(paste0(titulo,"\n"), gp=gpar(fontface="bold",fontsize=30))
  ppi <- 300
  dir.create("../figures/degdistributions/", showWarnings = FALSE)
  png(paste0("../figures/degdistributions/ALLdist_",red,"_",languageEl,".png"), width=(16*ppi), height=16*ppi, res=ppi)
  grid.arrange(e_degree,e_weight,i_degree,i_weight, ncol=2, nrow=2,top=title1 )
  dev.off()
}