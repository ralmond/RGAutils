
## Lets set up an artificial data set


tqr.H <- matrix(c(1/3,1/3,1/3,
              1/2,1/4,1/4,
              2/3,0,1/3,
              1/2,1/2,0),4,3,byrow=TRUE)
tqr.Q <- tqr.H > 0

tqr.n <- 100

tqr.XX <- matrix(rnorm(tqr.n*ncol(tqr.Q)),tqr.n)
tqr.YY <- tqr.XX%*%t(tqr.H) + matrix(rnorm(tqr.n*nrow(tqr.Q)),tqr.n)

tqr.Hout <- Qreg(tqr.YY,tqr.XX,tqr.Q)

colnames(tqr.XX) <- paste("X",1:ncol(tqr.XX),sep="")
colnames(tqr.YY) <- paste("Y",1:ncol(tqr.YY),sep="")
tqr.df <- as.data.frame(cbind(tqr.XX,tqr.YY))

tqr.Hexp <- matrix(NA,4,3)
tqr.Hexp[1,] <- lm(Y1~X1+X2+X3-1,data=tqr.df)$coef
tqr.Hexp[2,] <- lm(Y2~X1+X2+X3-1,data=tqr.df)$coef
tqr.Hexp[3,c(1,3)] <- lm(Y3~X1+X3-1,data=tqr.df)$coef
tqr.Hexp[4,1:2] <- lm(Y4~X1+X2-1,data=tqr.df)$coef

stopifnot(all(is.na(tqr.Hexp)==is.na(tqr.Hout)),
          max(abs(tqr.Hexp-tqr.Hout),na.rm=TRUE)< 1e-8)

