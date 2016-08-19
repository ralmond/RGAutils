### This is the matrix sweep operator:  Beaton (1964), Dempster (1969),
### Goodnight(1979).

### Goodnight's version seems to be slightly different than Dempster's

### Version 0.1  RGA 2016-03-16 -- Based on Goodnight
### Version 0.2  RGA 2016-05-07 -- Corrected to Dempster's

### Function matSweep
### Args:
###  A -- matrix to be swept
### ploc -- sequence of pivot locations.
### tol -- Check for zero pivot values (to avoid singular inversions)
### Output:  Matrix A with rows indicated by ploc swept out.
matSweep <- function (A, ploc= 1L:nrow(A), tol=.Machine$double.eps^0.5) {
  if (!is.matrix(A) || nrow(A) != ncol(A))
    stop("The first argument to matSweep must be a symmetric matrix")
  if (any(ploc < 1L) || any(ploc > nrow(A)) ||
      any( (ploc - round(ploc)) > tol))
    stop("Pivot locations must be integers referencing array rows.")
  for (k in ploc) { #Sweep on k
    D <- A[k,k]                           #Step 1
    if (abs(D) < tol) {
      stop("Singular pivot at location",k)
    }
    A[k,] <- A[k,]/D                      #Step 2
    for (i in 1L:nrow(A)) {               #Step 3
      if (i != k ) {
        B <- A[i,k]
        A[i,] <- A[i,] - B*A[k,]
        A[i,k] <- B/D
                                        # Goodnight: A[i,k] <- - B/D
      }
    }
    A[k,k] <- -1/D                         #Step 4
                                        # Goodnight:  A[k,k] <- 1/D
  }
  A
}

### Function revSweep
### Args:
###  A -- matrix to be unswept
### ploc -- sequence of pivot locations.
### tol -- Check for zero pivot values (to avoid singular inversions)
### Output: Matrix of the same size as A with the indicated rows unswept
revSweep <- function (A, ploc= 1L:nrow(A), tol=.Machine$double.eps^0.5) {
  if (!is.matrix(A) || nrow(A) != ncol(A))
    stop("The first argument to revSweep must be a symmetric matrix")
  if (any(ploc < 1L) || any(ploc > nrow(A)) ||
      any( (ploc - round(ploc)) > tol))
    stop("Pivot locations must be integers referencing array rows.")
  for (k in ploc) { # RevSweep(k)
    D <- A[k,k]                           #Step 1
    if (abs(D) < tol) {
      stop("Singular pivot at location",k)
    }
    A[k,] <- -A[k,]/D                      #Step 2
    for (i in 1L:nrow(A)) {               #Step 3
      if (i != k ) {
        B <- A[i,k]
        A[i,] <- A[i,] + B*A[k,]
        A[i,k] <- - B/D
      }
    }
    A[k,k] <- -1/D                         #Step 4
  }
  A
}

###
## I'm having some issues here.  I'm taking the formula from Little &
## Rubin (2002), p 151, and adapting it to fit into the Algorithm
## version of Goodnight.
## Now fixed using Dempster (1969), p 67.


#### Orthoganalization using Sweep opperator.
##
orthoInvert <- function (Q,tol=.Machine$double.eps^0.5) {
  p <- nrow(Q)
  T <- diag(NA,p)
  T[1,1] <- Q[1,1]
  B <- diag(1,p)
  A <- diag(1,p)
  for (i in 1:p) {
    Q <- matSweep(Q,i,tol)
    if (i <p) {
      T[i+1,i+1] <- Q[i+1,i+1]
      #print("T"); print(T,digits=5)
      B[(i+1):p,i] <- Q[(i+1):p,i]
      #print("B"); print(B,digits=5)
      A[i+1,1:i] <- -Q[i+1,1:i]
      #print("A"); print(B,digits=5)
    }
  }
  list(Qinv=-Q,T=T,B=B,A=A)
}


### The assimilate operator is defined by Dempster (1969, p 68).

matASM <- function (swpQ,Q123,Q33,ploc=1L:nrow(swpQ)) {
  ploc <- as.integer(ploc)
  if (any(ploc < 1L) || any(ploc>nrow(swpQ))) {
    stop("Swept indexes out of bounds.")
  }
  if (nrow(swpQ) != ncol(swpQ)) {
    stop("swpQ must be a square matrix.")
  }
  if (nrow(swpQ) != nrow(Q123)) {
    stop("swpQ and Q123 must have same number of rows.")
  }
  if (nrow(Q33) != nrow(Q33)) {
    stop("Q33 must be a square matrix.")
  }
  if (ncol(Q33) != ncol(Q123)) {
    stop("Q33 and Q123 must have same number of columns.")
  }

  rr <- (nrow(swpQ)+1):(nrow(swpQ)+nrow(Q33))
  KK <- matrix(NA,nrow(swpQ)+nrow(Q33),ncol(swpQ)+ncol(Q33))
  KK[1L:nrow(swpQ),1L:ncol(swpQ)] <- swpQ
  K11 <- swpQ[ploc,ploc,drop=FALSE]
  Q13 <- Q123[ploc,,drop=FALSE]
  H13 <- -K11%*%Q13
  KK[ploc,rr] <- H13
  KK[rr,ploc] <- t(H13)
  KK[rr,rr] <- Q33 - t(H13)%*%Q13
  ## Check for empty Q22 case
  if (length(ploc) < nrow(swpQ)) {
    H21 <- swpQ[-ploc,ploc,drop=FALSE]
    Q23 <- Q123[-ploc,,drop=FALSE]
    KK[-c(ploc,rr),rr] <- Q23 - H21%*%Q13
    KK[rr,-c(ploc,rr)] <- t(KK[-c(ploc,rr),rr])
  }
  KK
}


matMST <- function (Q, C=diag(nrow(Q)), do=1L:nrow(Q), done=integer()) {
  for (r1 in do) {
    undone <-setdiff(setdiff(1:nrow(Q),done),r1)
    ## Row of Q**21 on which we are operating.
    if (length(done)>0L) {
      qr1done <- Q[r1,done]
      C[r1,] <- C[r1,] - qr1done%*%C[done,]
    } else {
      qr1done <- numeric()
    }
    Ur1.norm <- sqrt(Q[r1,r1] - sum(qr1done^2))
    C[r1,] <- C[r1,]/Ur1.norm
    if (length(undone) > 0L) {
      if (length(done) > 0L) {
        qdud <- Q[done,undone]
        Q[undone,r1] <- Q[undone,r1] - qr1done%*%qdud
      }
      Q[undone,r1] <- Q[undone,r1]/Ur1.norm
      Q[r1,undone] <- Q[undone,r1]
    }
    Q[r1,r1] <- 1
    if (length(done)>0L) {
      Q[done,r1] <- 0
      Q[r1,done] <- 0
    }
    done <- c(done,r1)
  }
  return (list(Q=Q, C=C, done=done))

}

SSX <- function (X,weights=1,side="left") {
  X <- as.matrix(X)
  sideint <- pmatch(tolower(side),c("left","right"))
  if (is.na(sideint)) {
    stop("Side must be 'left' or 'right'.")
  }
  if (sideint==1) {
    X1 <- cbind(1,X)
  } else {
    X1 <- cbind(X,1)
  }
  if (!all(weights==1)) {
    X1w <- sweep(X1,1,weights,"*")
    return( crossprod(X1w,X1))
  } else {
    return( crossprod(X1))
  }
}
