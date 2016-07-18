Qreg <- function(Y,X,Q=matrix(TRUE,ncol(Y),ncol(X)),weights=1) {
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
  ## Find unique patters
  Qpat <- apply(sweep(Q,2,2^(0:(ncol(Q)-1)),"*"),1,sum)

  ## Loop unique patters.
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



