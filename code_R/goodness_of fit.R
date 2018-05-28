library("goft")
library("nortest")

datosn <- rnorm(1000,0,1)

st <- shapiro.test(datosn)
llt <- lillie.test(datosn)

datosln <- exp(datosn)

lnt <- lnorm_test(datosln)