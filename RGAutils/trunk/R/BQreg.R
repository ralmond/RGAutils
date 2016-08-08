BQreg <- function(Y,X,Q=matrix(TRUE,ncol(Y),ncol(X)),weights=1) {
  if (nrow(Y) != nrow(X)) {
    stop("X and Y not same length")
  }
  if (ncol(Y) != nrow(Q) || ncol(X) != ncol(Q)) {
    stop("Q matrix must be ",ncol(Y),"by",ncol(X))
  }
  ## Apply weights to X
  wX <- sweep(X,1,weights,"*")
  XtwX <- t(wX)%*%X
  ## Force Q to be logical
  Q <- array(as.logical(Q),dim(Q))
  ## Find unique patterns
  Qpat <- apply(sweep(Q,2,2^(0:(ncol(Q)-1)),"*"),1,sum)

  ## Loop unique patterns.
  H <- array(NA_real_,dim(Q))
  for (qp in unique(Qpat)) {
    whichvars <- which(qp==Qpat)
    mask <- Q[whichvars[1],]            #Same for all matching vars.
    ## Solve normal equations for appropriate subset
    beta <- solve(XtwX[mask,mask])%*%t(wX[,mask])%*%Y[,whichvars]
    H[whichvars,mask] <- t(beta)
  }
  H
}


## Builds sufficient statistics for BQ-regression from current
## paramters and data.

## Assume that X is completely observed, but that there may be missing
## data in Y.

MbuildT <- function(X,Y,B,b0,Syy.x,w=1) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  K <- ncol(X)
  J <- ncol(Y)


  ## Build prior augmented covariance matrix.
  sxT0 <- matrix(NA,K+J+1,K+J+1)

  ## Upper left corner is swept covariance matrix.
  totweight <- ifelse(length(w)==1,nrow(X)*w,sum(w))
  Txx <- matSweep(SSX(X,w),1)
  ## Change to variance esimate
  Txx[2:(K+1),2:(K+1)] <- Txx[2:(K+1),2:(K+1)]/totweight
  ## And sweep on X
  sxT0[1:(K+1),1:(K+1)] <- matSweep(Txx,2:(K+1))

  ## Off diagonal blocks are regression coefficients,
  ## Constant is in row/column 1.
  sxT0[1,(K+2):(K+J+1)] <- b0
  sxT0[2:(K+1),(K+2):(K+J+1)] <- t(B)
  sxT0[(K+2):(K+J+1),1] <- b0
  sxT0[(K+2):(K+J+1),2:(K+1)] <- B
  ## Lower right corner is residual vairance matrix.
  sxT0[(K+2):(K+J+1),(K+2):(K+J+1)] <- Syy.x
  ## Unsweep to get augmented covariance matrix
  T0 <- revSweep(sxT0,2:(K+1))
  # T0 is the old covariance matrix.

  ## Now build the final sum of cross product matrix.
  T <- matrix(0,K+J+1,K+J+1)
  XY1 <- cbind(1,X,Y)                   #will impute in this matrix
  M <- is.na(Y)

  ## Find unique patterns
  Mpat <- apply(sweep(M,2,2^(0:(ncol(M)-1)),"*"),1,sum)

  ## Loop unique patterns.
  for (mp in unique(Mpat)) {
    if (as.integer(mp)==0) next         #Complete data case
    whichobs <- which(mp==Mpat)
    whichvars <- 1+K+which(M[whichobs[1],])
    ## Produce regression
    swpT0 <- matSweep(T0,setdiff(1:nrow(T0),whichvars))
    ## Mean imputation into missing cells
    XY1[whichobs,whichvars] <-
      XY1[whichobs,-whichvars,drop=FALSE] %*%
      swpT0[-whichvars,whichvars,drop=FALSE]
    sumwt <- ifelse(length(w)==1,w*length(whichobs),sum(w[whichobs]))
    T[whichvars,whichvars] <- T[whichvars,whichvars] +
      sumwt*swpT0[whichvars,whichvars]
  }
  T <- crossprod(sweep(XY1,1,w,"*"),XY1) + T
  list(T=T,XY1=XY1,M=M)
}
