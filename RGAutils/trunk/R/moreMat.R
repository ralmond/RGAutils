### More Matrix functions.

## This is loosely based on is.positive.definite from the package
## matrixcalc.  However, the version I was using has a bug, and
## returns an error is situtations for which I want it to quietly
## return false.
is.pdm <- function (m, tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(m) || !is.numeric(m) || any(is.na(m))) return(FALSE)
  if (ncol(m) != nrow(m)) return(FALSE)
  ## Test for symmetry
  if (any(abs(m-t(m)) > tol)) return (FALSE)
  lambda <- eigen(m,symmetric=TRUE,only.values=TRUE)$values
  return (all(lambda >= tol))
}

dmvnmar <- function(x,mean=rep(0,p), sigma=diag(p), log=FALSE) {
  ## Initial checks borrowed from dmvnorm
  if (is.vector(x))
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  if (!missing(mean)) {
    if (!is.null(dim(mean)))
      dim(mean) <- NULL
    if (length(mean) != p)
      stop("mean and sigma have non-conforming size")
  }
  if (!missing(sigma)) {
    if (p != ncol(sigma))
      stop("x and sigma have non-conforming size")
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                     check.attributes = FALSE))
      stop("sigma must be a symmetric matrix")
  }
  retval <- rep(NA,nrow(x))
  M <- is.na(x)
  Mpat <- apply(sweep(M,2,2^(0:(ncol(M)-1)),"*"),1,sum)
  for (m in unique(Mpat)) {
    pat <-M[match(m,Mpat),]
    if (all(pat)) next
    xpat <- x[Mpat==m,!pat,drop=FALSE]
    mpat <- mean[!pat]
    Spat <- sigma[!pat,!pat]
    retval[Mpat==m] <- dmvnorm(xpat,mpat,Spat,log)
  }
  retval
}
