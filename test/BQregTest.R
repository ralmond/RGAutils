#library(mvnmle)

## Little & Rubin (2002, p 154)
#data(missvals)
#mvMux <- missvals[,c(3,5,1,2,4)]
data(mvMux)                             # locally cached version.

mvX <- mvMux[,1:2]
mvY <- mvMux[,3:5]

## For comparision purposes, start with numbers based on complete data only.
mv.cccov <- cov.wt(na.omit(mvMux),method="ML")
mv.ccm <- apply(mvMux,2,mean,na.rm=TRUE)
mv.ccT <- rbind(c(-1/nrow(mvMux),mv.cccov$center),
                 cbind(mv.cccov$center,mv.cccov$cov))

mv.ccTswp12 <- matSweep(mv.ccT,2:3)
mv.BB <- mv.ccTswp12[4:6,2:3]
mv.B0 <- mv.ccTswp12[4:6,1]
mv.B <- cbind(mv.B0,mv.BB)

mv.Syy.x <- mv.ccTswp12[4:6,4:6]

## This is actually another test of the sweep operator.
## Mainly used to make sure I'm getting the matrix multiplications
## right.

## mvXc <- sweep(as.matrix(mvX),2,mv.ccm[1:2])
## mvYc <- sweep(as.matrix(mvY),2,mv.ccm[3:5])

## alt.BB <- cov(mvYc[1:6,],mvXc[1:6,])%*%solve(var(mvXc[1:6,]))
## alt.B0 <- mv.ccm[3:5] - mv.BB%*%(mv.ccm[1:2])
## stopifnot(all(abs(alt.BB-mv.BB) < .00001),
##           all(abs(alt.B0-mv.B0) < .00001))

## ## Now checking how to run the prediciton equations.
## mv.Xplus <- as.matrix(cbind(1,mvX))
## t(mv.B) %*% mv.Xplus
## alt.Xpred <- sweep(as.matrix(mvX)%*%t(alt.BB),2,alt.B0,"+")
## mv.Xpred <- tcrossprod(mv.Xplus,mv.B)

## stopifnot(all(abs(alt.Xpred-mv.Xpred)<.00001))


## Reproducing internal joint matrix contstruction steps
exsT0 <- matrix(NA,6,6)
exsT0[1:3,1:3] <- matSweep(SSX(mvX),1)
exsT0[2:3,2:3] <- exsT0[2:3,2:3]/nrow(mvX)
exsT0[1:3,1:3] <- matSweep(exsT0[1:3,1:3],2:3)
exsT0[1,4:6] <- mv.B0
exsT0[2:3,4:6] <- t(mv.BB)
exsT0[4:6,1] <- mv.B0
exsT0[4:6,2:3] <- mv.BB
exsT0[4:6,4:6] <- mv.Syy.x

## Do a series of regression with the complete data part to get the
## expected values.
## Complete Data part
msig35c <-cov.wt(mvMux[1:6,1:2],method="ML")

## (x1, x2) ~ x3 + x5
lm1.35c <- lm(x1 ~ x3 + x5,data=mvMux[1:6,])
lm2.35c <- lm(x2 ~ x3 + x5, data=mvMux[1:6,])
sig12.35c <-cov.wt(cbind(x1=resid(lm1.35c),x2=resid(lm2.35c)),method="ML")$cov

## x4 ~ x1+ x2+ x3 + x5
## Order is important here!
lm4.1235c <- lm(x4 ~ x3 + x5 + x1 +x2, data=mvMux[1:6,])
sig4.1235c <- cov.wt(matrix(resid(lm4.1235c),ncol=1),method="ML")$cov
lm4.35c <- lm(x4 ~ x3 + x5, data=mvMux[1:6,])
sig4.35c <- cov.wt(matrix(resid(lm4.35c),ncol=1),method="ML")$cov

## Set up to test expected value of imputed matrix:
exXY1 <- cbind(1,mvMux)
p1 <- 7:9                               #Rows with only X4 missing
p2 <- 10:13                             #Rows with X1, X2, X4 missing

exXY1[p1,6] <- predict(lm4.1235c,newdata=exXY1[p1,])
exXY1[p2,6] <- predict(lm4.35c,newdata=exXY1[p2,])
exXY1[p2,5] <- predict(lm2.35c,newdata=exXY1[p2,])
exXY1[p2,4] <- predict(lm1.35c,newdata=exXY1[p2,])

## Code to track down small discrepency between results due to
## starting with different means.
## XY1mean <- apply(mvMux,2,mean,na.rm=TRUE)
## XYp1mean <- apply(mvMux[1:6,],2,mean,na.rm=TRUE)
## XY1means <- rbind(p0=XY1mean,p1=XYp1mean)

## p1x4 <- predict(lm4.1235c,newdata=as.data.frame(XY1means))
## p2x4 <- predict(lm4.35c,newdata=as.data.frame(XY1means))

## XY1meanss <- rbind(p0=XY1mean,p1=XYp1mean,
##                    p12=apply(mvMux[1:9,],2,mean,na.rm=TRUE))
## p2x2 <- predict(lm2.35c,newdata=as.data.frame(XY1meanss))
## p2x1 <- predict(lm1.35c,newdata=as.data.frame(XY1meanss))




## This is the test using the complete data only matrix as starting
## point.
Tstep <- EbuildT(mvX, mvY, mv.BB, mv.B0, mv.Syy.x, w=1)

stopifnot(all.equal(Tstep$XY1,as.matrix(exXY1),check.attributes=FALSE))

## Problem:  This is basically 1 step from an EM algorithm, so
## although it has done the imputations correctly, there is no easy
## way to build the expected values.

## Test agains ML estimatate from mvmle package.
## mle <- mlest(missvals)
## Hmm, this number doesn't match Little & Rubin
## I'll abandon this as a standard.


## This is basically recreating the calculations in Little and Rubin
## (2002, p 155).

## Complete Data part
msig35 <-cov.wt(mvMux[,1:2],method="ML")

## (x1, x2) ~ x3 + x5
lm1.35 <- lm(x1 ~ x3 + x5, na.action=na.omit,data=mvMux)
lm2.35 <- lm(x2 ~ x3 + x5, na.action=na.omit,data=mvMux)
sig12.35 <-cov.wt(cbind(x1=resid(lm1.35),x2=resid(lm2.35)),method="ML")$cov

## x4 ~ x1+ x2+ x3 + x5
## Order is important here!
lm4.1235 <- lm(x4 ~ x3 + x5 + x1 +x2, na.action=na.omit,data=mvMux)
sig4.1235 <- cov.wt(matrix(resid(lm4.1235),ncol=1),method="ML")$cov
lm4.35 <- lm(x4 ~ x3 + x5, na.action=na.omit,data=mvMux)
sig4.35 <- cov.wt(matrix(resid(lm4.35),ncol=1),method="ML")$cov

lrT1 <- rbind(x0=c(x0=-1,msig35$center),cbind(msig35$center,msig35$cov))
lrT2 <- rbind(cbind(matSweep(lrT1,2:3),x1=coef(lm1.35),x2=coef(lm2.35)),
              cbind(rbind(x1=coef(lm1.35),x2=coef(lm2.35)),
                    sig12.35))
lrT3 <- rbind(cbind(matSweep(lrT2,4:5),x4=coef(lm4.1235)),
              x4=c(coef(lm4.1235),sig4.1235))
lrT <- revSweep(lrT3,2:5)

lrTumS <-lrT[c(4,5,2,6,3),c(4,5,2,6,3)]
lrTumm <- lrT[1,c(4,5,2,6,3)]

## Actuall results from L & R
lrmuhat <- c( x1= 6.655,  x2=49.965, x3=11.769, x4=27.047, x5=95.423)
lrSigmahat <- matrix(c( 21.826,   20.864, -24.900,  -11.473,   46.953,
                        20.864,  238.012, -15.817, -252.072,  195.604,
                       -24.900,  -15.817,  37.870,   -9.599,  -47.566,
                       -11.473, -252.072,  -9.599,  294.183, -190.599,
                        46.953,  195.604, -47.566, -190.599,  208.905),
                     5,5,dimnames=list(paste("x",1:5,sep=""),
                                    paste("x",1:5,sep="")))

stopifnot(all.equal(lrmuhat,lrTumm,tolerance=.001))
all.equal(lrTumS,lrSigmahat,tolerance=.001)

## This is basically the input matrix.  As this is the MLE, we should
## reproduce it.
lrT3a <- revSweep(lrT3,4:5)


Tstep1 <- EbuildT(mvX, mvY, lrT3a[4:6,2:3], lrT3a[1,4:6], lrT3a[4:6,4:6], w=1)

Tsigmahat <- Tstep1$T[c(4,5,2,6,3),c(4,5,2,6,3)]
Tmuhat <- Tstep1$T[1,c(4,5,2,6,3)]

## Check the outputs
stopifnot(all.equal(lrTumm,Tmuhat))
stopifnot(all.equal(lrTumS,Tsigmahat))

data(twoskill)

## data1 has a bit of residual standard deviation.
That1 <- matSweep(SSX(as.matrix(twoskill.data1),1),1)
## Normalize covariance estimates
That1[-1,-1] <- That1[-1,-1]/nrow(twoskill.data1)

Best1 <- TQB(That1,twoskill.Q)

lmt1 <- lm(Test1~Skill1,data=twoskill.data1)
lmt2 <- lm(Test2~Skill2,data=twoskill.data1)
lmI1 <- lm(IntTest1~Skill1+Skill2,data=twoskill.data1)
lmI2 <- lm(IntTest2~Skill1+Skill2,data=twoskill.data1)

lm.B <- rbind(c(coef(lmt1),0),c(coef(lmt2)[1],0,coef(lmt2)[2]),
              coef(lmI1),coef(lmI2))
Best1.B <- cbind(Best1$b0,Best1$B)
stopifnot(all.equal(lm.B,Best1.B,check.attributes=FALSE))

lm.rsde <- c(var(resid(lmt1)),var(resid(lmt2)),
             var(resid(lmI1)),var(resid(lmI2)))
stopifnot(all.equal(lm.rsde*99/100,diag(Best1$Syy.x)))

## Once more with weights
## data1 has a bit of residual standard deviation.
wei <- runif(nrow(twoskill.data1),.1,.3)
That1w <- matSweep(SSX(as.matrix(twoskill.data1),wei),1)
## Normalize covariance estimates
That1w[-1,-1] <- That1w[-1,-1]/sum(wei)

Best1w <- RGAutils:::TQB(That1w,twoskill.Q)

lmt1w <- lm(Test1~Skill1,data=twoskill.data1,weights=wei)
lmt2w <- lm(Test2~Skill2,data=twoskill.data1,weights=wei)
lmI1w <- lm(IntTest1~Skill1+Skill2,data=twoskill.data1,weights=wei)
lmI2w <- lm(IntTest2~Skill1+Skill2,data=twoskill.data1,weights=wei)

lmw.B <- rbind(c(coef(lmt1w),0),c(coef(lmt2w)[1],0,coef(lmt2w)[2]),
              coef(lmI1w),coef(lmI2w))
Best1w.B <- cbind(Best1w$b0,Best1w$B)
stopifnot(all.equal(lmw.B,Best1w.B,check.attributes=FALSE))

lmw.resids <- data.frame(t1=resid(lmt1w),t2=resid(lmt2w),
                         I1=resid(lmI1w),I2=resid(lmI2w))

lmw.rsde <- cov.wt(lmw.resids,wei,method="ML")$cov

stopifnot(
    all.equal(diag(lmw.rsde),diag(Best1w$Syy.x),check.attributes=FALSE)
    )



####################################
###
## BQ Regression tests

Est0 <- BQreg(twoskill.data0[,1:2],twoskill.data0[,3:6],
              twoskill.Q)

## data0 has no noise, so should recover parameters exactly
stopifnot(all.equal(Est0$B,twoskill.B,check.attributes=FALSE),
          all.equal(Est0$b0,twoskill.b0,check.attributes=FALSE),
          Est0$converged==TRUE,Est0$iterations==0)


Est1 <- BQreg(twoskill.data1[,1:2],twoskill.data1[,3:6],
              twoskill.Q)
Est1.B <- cbind(Est1$b0,Est1$B)
stopifnot(all.equal(lm.B,Est1.B,check.attributes=FALSE))

lm.resids <- data.frame(t1=resid(lmt1),t2=resid(lmt2),
                        I1=resid(lmI1),I2=resid(lmI2))
lm.rcov <- cov.wt(lm.resids,center=FALSE,method="ML")$cov

stopifnot(
    all.equal(lm.rcov,Est1$Syy.x, check.attributes=FALSE)
    )

####
## Missing Data test
data(mvMux)

mvX <- mvMux[,1:2] # Complete Data
mvY <- mvMux[,3:5] # Data with missing values.

mvEst <- BQreg(mvX,mvY,maxstep=50)

## This is basically the input matrix.  As this is the MLE, we should
## reproduce it.

lrT3a <- revSweep(lrT3,4:5)

stopifnot(
  all.equal(mvEst$B,lrT3a[4:6,2:3],check.attributes=FALSE),
  all.equal(mvEst$b0,lrT3a[1,4:6],check.attributes=FALSE),
  all.equal(mvEst$Syy.x,lrT3a[4:6,4:6],check.attributes=FALSE),
  mvEst$converged==TRUE)








