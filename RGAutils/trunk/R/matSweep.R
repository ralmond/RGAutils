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
  ### Uses tail recursion to loop through
  matSweepAux(A,ploc,tol)
}

matSweepAux <- function (A, ploc, tol) {
  ## Base Case
  if (length(ploc) == 0L) return (A)
  ## Pivot first index
  k <- ploc[1]                          #Goodnight SWEEP(k) operator
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
  ## Recurse
  matSweepAux(A,ploc[-1],tol)
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
  ### Uses tail recursion to loop through
  revSweepAux(A,ploc,tol)
}

###
## I'm having some issues here.  I'm taking the formula from Little &
## Rubin (2002), p 151, and adapting it to fit into the Algorithm
## version of Goodnight.
## Now fixed using Dempster (1969), p 67.


revSweepAux <- function (A, ploc, tol) {
  ## Base Case
  if (length(ploc) == 0L) return (A)
  ## Pivot first index
  k <- ploc[1]                          #Goodnight SWEEP(k) operator
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
  ## Recurse
  revSweepAux(A,ploc[-1],tol)
}


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



