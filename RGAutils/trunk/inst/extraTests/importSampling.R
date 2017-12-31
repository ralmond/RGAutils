### Test using the BQ Regression to do importance sampling.

library(mvtnorm)
library(RGAutils)
#data(bqm)
load("~/Projects/RGAutils/data/bqm.RData")

dat <- bqm.obs0

Npart <- 5L
Nburnin <- 100L
Niter <- 100L
Nresamp <- 5L
## Number of M-steps to take during various phases of the algorithm.
Mburnin <- 5L
Miter <- 25L

## Determine base parameters from Q matrix and items.
Ntheta <- ncol(bqm.Q)
Nitem <- nrow(bqm.Q)
NN <- nrow(dat)

## Replicate data
repdat <- matrix(rep(dat,each=Npart),Npart*NN,Nitem)

## Initial PM parameters
Sigma <- diag(Ntheta)
## Start based on unit loading and 0 intercept for all measures.
mu <- scale(dat%*%bqm.Q)


##Set up initial draws for particles (theta realizations)
thetas <- matrix(rep(mu,each=Npart),Npart*NN,Ntheta) +
  rmvnorm(Npart*NN,sigma=Sigma)
weight <- rep(1/Npart,Npart*NN)

### Burn-in/warm-up loop

for (icyc in 1L:Nburnin) {
  cat("Burn in Cycle",icyc,"\n")

  bqrout <- BQreg(thetas,repdat,bqm.Q,weight,Mburnin)

  preddat <- sweep(thetas%*%t(bqrout$B),2,bqrout$b0,"+")
  sig <- bqrout$Syy.x
  if (!is.pdm(sig)) {
    sig <- diag(diag(sig))
  }

  ww <- dmvnorm(repdat-preddat,sigma=sig)
  w0 <- apply(matrix(ww,Npart),2,sum)
  weight <- ww/rep(w0,each=Npart)
}



## Guts of dmvnorm, which was mysteriously retruning 0's when sigma
## was not positive definite.
## dec <- chol(bqrout$Syy.x)
## tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
## rss <- colSums(tmp^2)
## logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 *
##           pi) - 0.5 * rss

