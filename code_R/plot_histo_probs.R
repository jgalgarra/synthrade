exp_data <- read.delim("../results/RedAdyCom1970_FILT_W_1.txt", header=FALSE)
smat <- sum(exp_data)
Pr_E <- rowSums(exp_data)/smat
Pr_I <- colSums(exp_data)/smat
Pr_Tot <- Pr_E %o% Pr_I
hist(log(Pr_I),breaks=10)
hist(log(Pr_E),breaks=10)
hist(log(Pr_Tot),breaks=10)

exp_data <- read.delim("../tfresults/RedAdyCom1970_FILT_W_1.txt", header=FALSE)
smat <- sum(exp_data)
Pr_E <- rowSums(exp_data)/smat
Pr_I <- colSums(exp_data)/smat
Pr_Tot <- Pr_E %o% Pr_I
hist(log(Pr_I),breaks=20,main="PR_I tf")
hist(log(Pr_E),breaks=20,main="PR_E tf")
hist(log(Pr_Tot),breaks=20,main="P_tot tf")