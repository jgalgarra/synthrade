library("goft")
library(MASS)

experiment_files <- Sys.glob("../results/RedAdyCom1970_FILT_W_*.txt")
first_file = TRUE
for (file in experiment_files)
{
  exp_data <- read.delim(file, header=FALSE)
  smat <- sum(exp_data)
  exp_data <- exp_data/smat
  if (first_file){
    tot_data <- exp_data
    first_file <- FALSE
  }else{
    tot_data <- tot_data+exp_data
  }
}

#exp_data <- tot_data

Pr_E <- rowSums(exp_data)/smat
Pr_I <- colSums(exp_data)/smat
dfPr_Tot <- Pr_E %o% Pr_I
Pr_Tot <- as.vector(dfPr_Tot)/max(dfPr_Tot)
lPr_Tot <- log(Pr_Tot)
lPr_Tot <- lPr_Tot - min(lPr_Tot) + 0.00000000000000001

mydata <- lPr_Tot
mydata <- exp_data[exp_data>0]

fit_ln <- fitdist(mydata, "lnorm")
summary(fit_ln)
fit_w  <- fitdist(mydata, "weibull")
fit_g  <- fitdist(mydata, "gamma")
par(mfrow=c(2,2))
plot.legend <- c("Weibull", "lognormal", "gamma")
denscomp(list(fit_w, fit_g, fit_ln), legendtext = plot.legend)
cdfcomp (list(fit_w, fit_g, fit_ln), legendtext = plot.legend)
qqcomp  (list(fit_w, fit_g, fit_ln), legendtext = plot.legend)
ppcomp  (list(fit_w, fit_g, fit_ln), legendtext = plot.legend)

gofstat(list(fit_ln, fit_w, fit_g), fitnames = c("lognormal", "Weibull", "gamma"))
hist(log(Pr_I),breaks=10)
hist(log(Pr_E),breaks=10)
hist(log(Pr_Tot),breaks=10)


lt_Tot <- lillie.test(log(Pr_Tot))

exp_data <- read.delim("../results/RedAdyCom1000_FILT_W_1.txt", header=FALSE)
smat <- sum(exp_data)
Pr_E <- rowSums(exp_data)/smat
Pr_I <- colSums(exp_data)/smat
Pr_Tot <- Pr_E %o% Pr_I
hist(log(Pr_I),breaks=20,main="PR_I tf")
hist(log(Pr_E),breaks=20,main="PR_E tf")
hist(log(Pr_Tot),breaks=20,main="P_tot tf")