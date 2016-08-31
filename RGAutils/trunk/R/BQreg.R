BQreg <- function(X,Y,Q=matrix(TRUE,ncol(Y),ncol(X)),weights=1,
                  maxstep=25L,tolerance=10*sqrt(.Machine$double.eps),
                  trace=FALSE) {
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  if (is.data.frame(Y)) {
    Y <- as.matrix(Y)
  }
  if (nrow(Y) != nrow(X)) {
    stop("X and Y not same length")
  }
  J <- ncol(Y)
  K <- ncol(X)
  if (nrow(Q) != J || ncol(Q) !=K) {
    stop("Q matrix must be ",J,"by",K)
  }
  ## Calculate sufficient statistic matrix, T
  ## First, do a first pass estimate using complete cases.
  complete <- !apply(is.na(Y),1,any)
  wc <- weights

  T0 <- SSX(cbind(X[complete,],Y[complete,]),
            ifelse(length(weights)>1,weights[complete],weights))
  sumwc <- T0[1,1]                      #Sum of weights
  T0 <- matSweep(T0,1)
  T0[-1,-1] <- T0[-1,-1]/sumwc
  if (all(complete)) {
    converged <- TRUE
    iterations <- 0
    BQ <- TQB(T0,Q)
    B <- BQ$B
    b0 <- BQ$b0
    Syy.x0 <- BQ$Syy.x
  } else {
    ## Missing data, EM time:
    converged <- FALSE
    for (iterations in 1L:maxstep) {
      BSest <- TQB(T0,Q)
      T1 <- EbuildT(X,Y,BSest$B,BSest$b0,BSest$Syy.x,weights)$T
      diff <- max(abs(T1-T0))
      T0 <- T1
      if (trace)
        cat("Iteration,",iterations,"max difference:",diff,"\n")
      if (diff < tolerance) {
        converged <- TRUE
        break
      }
    }                                   #End of loop
    BQ <- TQB(T0,Q)
    B <- BQ$B
    b0 <- BQ$b0
    Syy.x0 <- BQ$Syy.x
  }                                     #End of EM
  ## Calculate residual covariance matrix.
  resid <- Y - sweep(tcrossprod(X,B),2,b0,"+")
  na.resid <- is.na(resid)
  resid[na.resid] <- 0                  #Zero residual if NA
  if (length(weights)==1) {
    weights <- rep(weights,nrow(X))
  }
  Syy.x <- crossprod(sweep(resid,1,weights,"*"),resid)
  if (any(na.resid)) {
    ## Adjust elements of covariance matrix as needed.
    for (j1 in 1:J) {
      whichmiss <- na.resid[,j1]
      if (!any(whichmiss)) next         #Skip complete data
      ## Diagonal adjustment
      Syy.x[j1,j1] <- Syy.x[j1,j1] + Syy.x0[j1,j1]*sum(weights[whichmiss])
      if (j1==J) break                  #Cross product terms finished
      for (j2 in (j1+1):J) {            #Cross product terms
        whichmiss <- na.resid[,j1] & na.resid[,j2]
        if (!any(whichmiss)) next       #Skip complete data where both
                                        #are not missing.
        Syy.x[j1,j2] <- Syy.x[j1,j2] + Syy.x0[j1,j2]*sum(weights[whichmiss])
        Syy.x[j2,j1] <- Syy.x[j1,j2]    #Symmetry
      }                                 #next j2
    }                                   #next j1
  Syy.x <- Syy.x/sum(weights)
  }                                     #End of missing data processing

  list(B=B,b0=b0,Syy.x=Syy.x,
       converged=converged,iterations=iterations)
}

TQB <- function (T,Q) {

  ## Force Q to be logical
  Q <- array(as.logical(Q),dim(Q))
  K <- ncol(Q)
  J <- nrow(Q)

  ## Coefficients
  B <- array(0,dim(Q))
  b0 <- rep(0,nrow(Q))
  ## Residual covariance matrix.
  Syy.x <- matrix(0,J,J)

  ## Find unique patterns
  Qpat <- apply(sweep(Q,2,2^(0:(ncol(Q)-1)),"*"),1,sum)

  ## Loop unique patterns.

  for (qp in unique(Qpat)) {
    whichY <- which(qp==Qpat)
    whichX <- which(Q[whichY[1],])   #Same for all matching vars.
    sT <- matSweep(T,1+whichX)
    B[whichY,whichX] <- sT[1+K+whichY,1+whichX]
    b0[whichY] <- sT[1+K+whichY,1]
    Syy.x[whichY,whichY] <- sT[1+K+whichY,1+K+whichY]
  }
  list(B=B,b0=b0,Syy.x=Syy.x)
}


## Builds sufficient statistics for BQ-regression from current
## paramters and data.

## Assume that X is completely observed, but that there may be missing
## data in Y.

EbuildT <- function(X,Y,B,b0,Syy.x,w=1) {
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
    ## Already swept on row 1, so start there.
    swpT0 <- matSweep(T0,setdiff(2:nrow(T0),whichvars))
    ## Mean imputation into missing cells
    XY1[whichobs,whichvars] <-
      XY1[whichobs,-whichvars,drop=FALSE] %*%
      swpT0[-whichvars,whichvars,drop=FALSE]
    sumwt <- ifelse(length(w)==1,w*length(whichobs),sum(w[whichobs]))
    T[whichvars,whichvars] <- T[whichvars,whichvars] +
      sumwt*swpT0[whichvars,whichvars]
  }
  T <- crossprod(sweep(XY1,1,w,"*"),XY1) + T
  T <- matSweep(T,1)
  sumwt <- ifelse(length(w)==1,nrow(XY1),sum(w))
  T[-1,-1] <- T[-1,-1]/sumwt
  list(T=T,XY1=XY1,M=M)
}
