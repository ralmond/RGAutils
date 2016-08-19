## Generated data used for various tests.

library(mvtnorm)

twoskill.cormat <- diag(2)
twoskill.cormat[1,2] <- twoskill.cormat[2,1] <- .75

theta <- rmvnorm(100,sigma=cormat)
colnames(theta) <- paste("Skill",1:2,sep="")

twoskill.B <- matrix(c(10,0,
                  0,10,
                  7,8,
                  9,4), 4,2,byrow=TRUE,
                dimnames=list(obs=c("Test1","Test2","IntTest1","IntTest2"),
                              skill=colnames(theta))
                )
twoskill.b0 <- c(Test1=50,Test2=100,IntTest1=200,IntTest2=75)
twoskill.Q <- twoskill.B !=0

obs <- sweep(theta%*%t(twoskill.B),2,twoskill.b0,"+")
twoskill.resid.sd <- diag(c(2,4,6,8))
obs1 <- obs + rmvnorm(100,sigma=twoskill.resid.sd)


twoskill.data0 <- data.frame(theta,obs)
twoskill.data1 <- data.frame(theta,obs1)
save(twoskill.data0,twoskill.data1,twoskill.B,twoskill.b0,twoskill.Q,
     twoskill.resid.sd,
     file="RGAutils/data/twoskill.RData")
