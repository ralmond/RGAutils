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



####
## BQMiss Model

bqm.n <- 1000

bqm.Sigma <- diag(3)
bqm.Sigma[1,2] <- bqm.Sigma[2,1] <- .6
bqm.Sigma[1,3] <- bqm.Sigma[3,1] <- -.5
bqm.Sigma[2,3] <- bqm.Sigma[3,2] <- -.4

bqm.Theta <- rmvnorm(bqm.n,sigma=bqm.Sigma)
colnames(bqm.Theta) <- paste("Skill",1:ncol(bqm.Theta),sep="")

bqm.B <- matrix(c(1/3,1/3,1/3,
                  1/2,1/4,1/4,
                  .8,.2,.2,
                  .2,.8,.2,
                  .2,.2,.8,
                  2/3,0,1/3,
                  1/2,1/2,0,
                  1,0,0,
                  .8,0,0,
                  .4,.6,0
                  ),10,3,byrow=TRUE,
                dimnames=list(obs=paste("Test",1:10,sep=""),
                              skill=colnames(bqm.Theta))
                )
bqm.Q <- bqm.B > 0

bqm.B0 <- c(-1,-1,1,1,1,-1,0,0,0,0)

bqm.obs0 <- sweep(bqm.Theta%*%t(bqm.B),2,bqm.B0,"+")
bqm.resid.sd <- diag(rep(.4,10))
bqm.obs1 <- bqm.obs0 + rmvnorm(bqm.n,sigma=bqm.resid.sd)

## Randomly make some data missing.
bqm.obs2 <- bqm.obs1
bqm.m <- matrix(FALSE,bqm.n,10)
for (k in 1:5) {
  bqm.m[sample.int(bqm.n,round(.1*k*bqm.n)),k] <- TRUE
  is.na(bqm.obs2[bqm.m[,k],k]) <- TRUE
}

bqm.data0 <- data.frame(bqm.Theta,bqm.obs0)
bqm.data1 <- data.frame(bqm.Theta,bqm.obs1)
bqm.data2 <- data.frame(bqm.Theta,bqm.obs2)

save(bqm.data0,bqm.data1,bqm.data2,
     bqm.Theta,bqm.obs0,bqm.obs1,bqm.obs2,
     bqm.B,bqm.B0,bqm.Q,
     bqm.resid.sd,bqm.Sigma,
     file="~/Projects/RGAutils/data/bqm.RData")
