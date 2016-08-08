library(mvnmle)

## Little & Rubin (2002, p 154)
data(missvals)

mvMux <- missvals[,c(3,5,1,2,4)]

mvX <- mvMux[,1:2]
mvY <- mvMux[,3:5]

mv.cccov <- cov(mvMux,use="complete")
mv.ccm <- apply(mvMux,2,mean,na.rm=TRUE)
mv.ccT <- rbind(c(-1/nrow(mvMux),mv.ccm),
                cbind(mv.ccm,mv.cccov))

mv.ccTswp12 <- matSweep(mv.ccT,2:3)
mv.BB <- mv.ccTswp12[4:6,2:3]
mv.B0 <- mv.ccTswp12[4:6,1]
mv.B <- cbind(mv.B0,mv.BB)

mv.Syy.x <- mv.ccTswp12[4:6,4:6]

mvXc <- sweep(as.matrix(mvX),2,mv.ccm[1:2])
mvYc <- sweep(as.matrix(mvY),2,mv.ccm[3:5])

## This is actually another test of the sweep operator.
alt.BB <- cov(mvYc[1:6,],mvXc[1:6,])%*%solve(var(mvXc[1:6,]))
alt.B0 <- mv.ccm[3:5] - mv.BB%*%(mv.ccm[1:2])
stopifnot(all(abs(alt.BB-mv.BB) < .00001),
          all(abs(alt.B0-mv.B0) < .00001))

## Now checking how to run the prediciton equations.
mv.Xplus <- as.matrix(cbind(1,mvX))
t(mv.B) %*% mv.Xplus
alt.Xpred <- sweep(as.matrix(mvX)%*%t(alt.BB),2,alt.B0,"+")
mv.Xpred <- tcrossprod(mv.Xplus,mv.B)

stopifnot(all(abs(alt.Xpred-mv.Xpred)<.00001))



Tstep <- MbuildT(mvX, mvY, mv.BB, mv.B0, mv.Syy.x, w=1)
Ts1 <- matSweep(Tstep$T,1)
Ts1[c(4,5,2,6,3),c(4,5,2,6,3)]/13

mle <- mlest(missvals)
## Hmm, this number doesn't match Little & Rubin
## Complete Data part
msig35 <-cov.wt(missvals[,c(3,5)],method="ML")

## (x1, x2) ~ x3 + x5
lm1.35 <- lm(x1 ~ x3 + x5, na.action=na.omit,data=missvals)
lm2.35 <- lm(x2 ~ x3 + x5, na.action=na.omit,data=missvals)
sig12.35 <-cov.wt(cbind(x1=resid(lm1.35),x2=resid(lm2.35)),method="ML")$cov

## x4 ~ x1+ x2+ x3 + x5
## Order is important here!
lm4.1235 <- lm(x4 ~ x3 + x5 + x1 +x2, na.action=na.omit,data=missvals)
sig4.1235 <- cov.wt(matrix(resid(lm4.1235),ncol=1),method="ML")$cov

lrT1 <- rbind(x0=c(x0=-1,msig35$center),cbind(msig35$center,msig35$cov))
lrT2 <- rbind(cbind(matSweep(lrT1,2:3),x1=coef(lm1.35),x2=coef(lm2.35)),
              cbind(rbind(x1=coef(lm1.35),x2=coef(lm2.35)),
                    sig12.35))
lrT3 <- rbind(cbind(matSweep(lrT2,4:5),x4=coef(lm4.1235)),
              x4=c(coef(lm4.1235),sig4.1235))
lrT <- revSweep(lrT3,2:6)

