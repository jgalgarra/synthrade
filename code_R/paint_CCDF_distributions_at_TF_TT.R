library("grid")
library("gridExtra")
library("igraph")
library("ggplot2")
source("aux_functions_matrix.R")
source("parse_command_line_args.R")


gen_links_strength_distribution <- function(red,series, colors, seq_breaks = c(1,5,10,20,50,100), empirical = FALSE)
{
  
  
  gen_ls_data_frame <- function(input_matrix,tipo,tamanyo,nalpha,serie,titlestr)
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
    dfdeg$weight <- 0
    for (k in 1:nrow(dfdeg)){
      if (dfdeg$type[k] == "Importer")
      {
        indice <- as.numeric(strsplit(dfdeg$nodename[k],"Importer")[[1]][2])
        dfdeg$weight[k] <- sum(as.numeric(dr$m[indice,]))
      }
      else
      {
        indice <- as.numeric(strsplit(dfdeg$nodename[k],"Exporter")[[1]][2])
        dfdeg$weight[k] <- sum(as.numeric(dr$m[,indice]))
      }
    }
    dfdeg$weight <- dfdeg$weight /max(as.numeric(dfdeg$weight))
    
    ddeg_exporter <- dfdeg[dfdeg$type == "Exporter",]
    degree <- ddeg_exporter$degree
    weight <- ddeg_exporter$weight
    datosplot <- data.frame("degree" = degree, "strength" = weight)
    dpexp <- datosplot
    
    mod <- lm(datosplot$strength ~ datosplot$degree)
    etmodel <- sprintf("log10 s = %.4f log10 d %.4f     Adj. R^2 = %0.3f",as.numeric(mod[[1]][2]),as.numeric(mod[[1]][1]),summary(mod)$adj.r.squared)
    exptf <- ggplot(datosplot,aes(x=degree,y=strength))+geom_point(color="blue",alpha=0.5)+scale_x_log10()+scale_y_log10()+
      ggtitle(paste0("Exporters at ",titlestr))+
      geom_smooth(method = "lm", se = FALSE, show.legend = TRUE,color="grey50",linetype = "dashed")+
      geom_text(x=1, y=0,label=etmodel, size = 5)+xlab("Degree")+ylab("Normalized strength")+
      theme_bw() +  theme(plot.title = element_text(hjust = 0.5, size = 18),
                          axis.title.x = element_text(color="grey30", size = 15, face="bold"),
                          axis.title.y = element_text(color="grey30", size= 15, face="bold"),
                          legend.title=element_blank(),
                          legend.position = "top",
                          legend.text=element_text(size=10),
                          panel.grid.minor = element_blank(),
                          axis.text.x = element_text(face="bold", color="grey30", size=14),
                          axis.text.y = element_text(face="bold", color="grey30", size=14)
      )


    ddeg_importer <- dfdeg[dfdeg$type == "Importer",]
    degree <- ddeg_importer$degree
    weight <- ddeg_importer$weight
    datosplot <- data.frame("degree" = degree, "strength" = weight)
    dpimp <- datosplot
    mod <- lm(datosplot$strength ~ datosplot$degree)
    etmodel <- sprintf("log10 s = %.4f log10 d %.4f     Adj. R^2 = %0.3f",as.numeric(mod[[1]][2]),as.numeric(mod[[1]][1]),summary(mod)$adj.r.squared)

    imptf <- ggplot(datosplot,aes(x=degree,y=strength))+geom_point(colour="red",alpha=0.5)+scale_x_log10()+scale_y_log10()+
         ggtitle(paste0("Importers at ",titlestr))+
         geom_smooth(method = "lm", se = FALSE, show.legend = TRUE,color="grey50",linetype = "dashed")+
         geom_text(x=1, y=max(log10(datosplot$strength)),label=etmodel, size = 5)+xlab("Degree")+ylab("Normalized strength")+
      theme_bw() +  theme(plot.title = element_text(hjust = 0.5, size = 18),
                          axis.title.x = element_text(color="grey30", size = 15, face="bold"),
                          axis.title.y = element_text(color="grey30", size= 15, face="bold"),
                          legend.title=element_blank(),
                          legend.position = "top",
                          legend.text=element_text(size=10),
                          panel.grid.minor = element_blank(),
                          axis.text.x = element_text(face="bold", color="grey30", size=14),
                          axis.text.y = element_text(face="bold", color="grey30", size=14)
      )


    calc_values <- list("imptf" = imptf, "exptf" = exptf, "data_exp" = dpexp, "data_imp" = dpimp)
    return(calc_values)
  }
  
  experiment <- 1
  
  if (!empirical)
  {
    dred <- gsub(TFstring,"",red)
    subdir <- "TFMatrix/"
    ficheros <- Sys.glob(paste0("../results/",subdir,red,"_W_",experiment,".txt"))
    for (j in ficheros){
      sim_matrix <- read.table(j,sep="\t")
      plots_TF <- gen_ls_data_frame(sim_matrix,"Simulated",0.5,0.02,series,"TF")
    }
    subdir <- ""
    ficheros <- Sys.glob(gsub("TF_","",paste0("../results/",subdir,red,"_W_1",".txt")))
    for (j in ficheros){
      sim_matrix <- read.table(j,sep="\t")
      plots_final <- gen_ls_data_frame(sim_matrix,"Simulated",0.5,0.02,series,"TT")
    }
  }
  
  else
  {
    dred <- gsub(TFstring,"",red)
    subdir <- "data/"
    ficheros <- Sys.glob(paste0("../",subdir,dred,".txt"))
    for (j in ficheros){
      sim_matrix <- read.table(j,sep="\t")
      plots_final <- gen_ls_data_frame(sim_matrix,"Empirical",0.5,0.02,series,"TT")
      plots_TF <- plots_final
    }

  }
  
  calc_values <- list("plots_TF" = plots_TF, "plots_final" = plots_final)
  return(calc_values)

}

plot_ccdf_fit <- function(datosinput,titlestr="",dcol="red")
{
  
  datosinput$strength <- datosinput$strength/sum(as.numeric(datosinput$strength))
  occur <- datosinput$degree
  woccur <- datosinput$strength
  woccur <- woccur[order(woccur)]
  alpha_level = 0.8
  p = occur/sum(occur)
  dy = rev(cumsum(rev(p)))
  dx = occur
  wp = woccur/sum(woccur)
  wdy = rev(cumsum(rev(wp)))
  wdx = woccur
  auxdf <- data.frame(dx,dy)
  auxdfw <- data.frame(wdx,wdy)
  auxdf <- auxdf[auxdf$dx > 0,]
  distd_plot <- ggplot(data = auxdf) +
    geom_point(aes(x = sqrt(dx), y = dy), color=dcol)+
    #scale_x_log10() + scale_y_log10() 
  scale_x_sqrt()
  distw_plot <- ggplot(data = auxdfw) +
    geom_point(aes(x = wdx, y = wdy), color=dcol) +
    scale_x_sqrt()
    #scale_x_log10() + scale_y_log10() 
    # xlab(paste(year,series,"Degree")) + 
    # ylab("CCDF") + #scale_shape_manual(values=c(1,16)) +
    # scale_alpha(guide = 'none') + scale_size_identity()  +
    # theme_bw() +
    # theme(
    #   axis.title.x = element_text(color="grey30", size = 11),
    #   axis.title.y = element_text(color="grey30", size= 11),
    #   legend.title=element_blank(),
    #   legend.position = "top",
    #   legend.text=element_text(size=10),
    #   axis.text.x = element_text(face="bold", color="grey30", size=10),
    #   axis.text.y = element_text(face="bold", color="grey30", size=10)
    # )
  
  # datatrf <- datosplot
  # datatrf$log10_degree <- log10(datosplot$degree)^2
  # datatrf$log10_strength <- log10(datosplot$strength)
  # mod <- lm(datatrf$log10_strength ~ datatrf$log10_degree^2)
  # minx <- min(sqrt(datatrf$log10_degree))
  # maxx <- round(max(sqrt(datatrf$log10_degree)))
  # 
  # etmodel <- sprintf("log10 s = %.4f (log10 d)^2 %.4f     Adj. R^2 = %0.3f",as.numeric(mod[[1]][2]),as.numeric(mod[[1]][1]),summary(mod)$adj.r.squared)
  # imptf <- ggplot(datatrf,aes(x=log10_degree,y=log10_strength))+geom_point(color=dcol,alpha=0.5)+
  #   ggtitle(titlestr)+xlab("Degree")+ylab("Normalized strength")+
  #   scale_x_continuous(breaks=c(0,1,4),labels=c(1,10,100))+
  #   scale_y_continuous(breaks=c(0,-2,-4),labels=c("1","1e-02","1e-04"))+
  #   geom_smooth(method = "lm", se = FALSE, show.legend = TRUE,color="grey50",linetype = "dashed")+
  #   geom_text(x=2, y=max(datatrf$log10_strength),label=etmodel, size = 5)+
  #   theme_bw() +  theme(plot.title = element_text(hjust = 0.5, size = 18),
  #                       axis.title.x = element_text(color="grey30", size = 15, face="bold"),
  #                       axis.title.y = element_text(color="grey30", size= 15, face="bold"),
  #                       legend.title=element_blank(),
  #                       legend.position = "top",
  #                       legend.text=element_text(size=10),
  #                       panel.grid.minor = element_blank(),
  #                       axis.text.x = element_text(face="bold", color="grey30", size=14),
  #                       axis.text.y = element_text(face="bold", color="grey30", size=14)
  #   )
  # return(imptf)
  print(dist_plot)
  return()
}


plot_sq_fit <- function(datosplot,titlestr="",dcol="red")
{

  datatrf <- datosplot
  datatrf$log10_degree <- log10(datosplot$degree)^2
  datatrf$log10_strength <- log10(datosplot$strength)
  mod <- lm(datatrf$log10_strength ~ datatrf$log10_degree^2)
  minx <- min(sqrt(datatrf$log10_degree))
  maxx <- round(max(sqrt(datatrf$log10_degree)))

  etmodel <- sprintf("log10 s = %.4f (log10 d)^2 %.4f     Adj. R^2 = %0.3f",as.numeric(mod[[1]][2]),as.numeric(mod[[1]][1]),summary(mod)$adj.r.squared)
  imptf <- ggplot(datatrf,aes(x=log10_degree,y=log10_strength))+geom_point(color=dcol,alpha=0.5)+
    ggtitle(titlestr)+xlab("Degree")+ylab("Normalized strength")+
    scale_x_continuous(breaks=c(0,1,4),labels=c(1,10,100))+
    scale_y_continuous(breaks=c(0,-2,-4),labels=c("1","1e-02","1e-04"))+
    geom_smooth(method = "lm", se = FALSE, show.legend = TRUE,color="grey50",linetype = "dashed")+
    geom_text(x=2, y=max(datatrf$log10_strength),label=etmodel, size = 5)+
    theme_bw() +  theme(plot.title = element_text(hjust = 0.5, size = 18),
                        axis.title.x = element_text(color="grey30", size = 15, face="bold"),
                        axis.title.y = element_text(color="grey30", size= 15, face="bold"),
                        legend.title=element_blank(),
                        legend.position = "top",
                        legend.text=element_text(size=10),
                        panel.grid.minor = element_blank(),
                        axis.text.x = element_text(face="bold", color="grey30", size=14),
                        axis.text.y = element_text(face="bold", color="grey30", size=14)
                        )
  return(imptf)
}

TFstring = "TF_"
files <- paste0(TFstring,"RedAdyCom",seq(ini_seq,end_seq))

for (orig_file in files)
{
  red <- paste0(orig_file,"_FILT")
  redorig <- gsub(TFstring,"",red)                #Empirical data
  series = "Exporter"
  year=gsub("_FILT","",strsplit(red,"RedAdyCom")[[1]][-1])
  grafs <- gen_links_strength_distribution(red,series,"blue",empirical = FALSE)
  data_e <- grafs$plots_final$data_exp
  data_i <- grafs$plots_final$data_imp
  
  sqe <- plot_sq_fit(data_e, titlestr = "Synthetic Exporters at TT", dcol="blue")
  sqi <- plot_sq_fit(data_i, titlestr = "Synthetic Importers at TT", dcol="red")
  
  grafsemp <-  gen_links_strength_distribution(red,series,"blue",empirical = TRUE)
  data_e_emp <- grafsemp$plots_final$data_exp
  data_i_emp <- grafsemp$plots_final$data_imp
  sqe_emp <- plot_sq_fit(data_e_emp, titlestr = "Empirical Exporters", dcol="blue")
  sqi_emp <- plot_sq_fit(data_i_emp, titlestr = "Empricial Importers", dcol="red")
  
  plot_ccdf_fit(data_e_emp, titlestr = "Empirical Exporters", dcol="blue")

  dir.create("../figures/linksstrength/", showWarnings = FALSE)
  ppi <- 300
  png(paste0("../figures/linksstrength/LS_SYNTH_",red,".png"), width=(16*ppi), height=12*ppi, res=ppi)
  # grid.arrange(grafs$plots_TF$imptf,grafs$plots_final$imptf, sqi, grafs$plots_TF$exptf,
  #             grafs$plots_final$exptf,sqe,ncol=3, nrow=2)
  grid.arrange(grafs$plots_TF$imptf,sqi, grafs$plots_TF$exptf,sqe,ncol=2, nrow=2)
  dev.off()
  
  ppi <- 300
  png(paste0("../figures/linksstrength/LS_EMP_",red,".png"), width=(16*ppi), height=12*ppi, res=ppi)
  grid.arrange(sqi,sqi_emp,sqe,sqe_emp,ncol=2, nrow=2)
  dev.off()
  
  dir.create("../figures/linksstrength/", showWarnings = FALSE)
  ppi <- 300
  png(paste0("../figures/linksstrength/LS_SYNTH_",red,".png"), width=(22*ppi), height=12*ppi, res=ppi)
  grid.arrange(grafs$plots_TF$imptf, sqi, sqi_emp, grafs$plots_TF$exptf,
             sqe, sqe_emp, ncol=3, nrow=2)
  dev.off()
}