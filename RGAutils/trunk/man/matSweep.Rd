\name{matSweep}
\alias{matSweep}
\alias{revSweep}
\title{Beaton (Dempster) Sweep operator for Matrixes}
\description{
  
  The sweep operator takes a covariance matrix and replaces the
  elements corresponding to the swept out rows and columns with the
  inverse correlation matrix.  The remaining rows and columns give the
  regression coefficients for predicting the unswept variables from the
  swept variables and the residual covariance matrix.  The reverse sweep
  operator is its inverse.

  This version of the sweep operator follows Dempster(1969) which is
  slightly different from the one given in Goodnight (1979).

}
\usage{
matSweep(A, ploc = 1L:nrow(A), tol = .Machine$double.eps^0.5)
revSweep(A, ploc = 1L:nrow(A), tol = .Machine$double.eps^0.5)
}
\arguments{
  \item{A}{
    A symmetric matrix representing a correlation matrix or an augmented
    correlation matrix.
  }
  \item{ploc}{
    A vector of integers giving the rows and columns to be swept.  
  }
  \item{tol}{
    A machine tolerance value to check for 0 pivots (singular
    matrixes). 
  }
}
\details{

  The sweep operator was introduced by Beaton (1964) as a way of
  inverting matrixes in place (important with the limited memory
  computers of the time).  Dempster (1969) uses Beaton's sweep operator
  and gives it an interpretation in terms of a regression in which the
  unswept variables are predicted from the swept variables.  The
  original purpose for inverting matrixes has been superseded by more
  numerically stable algorithms (see \code{\link[base]{solve}}), but the
  interpretation in terms of various partial regression models is
  interesting.  In particular, Little and Rubin (2002) use it in
  building imputation models for missing data.

  Let \eqn{\bold{Q}} by an \eqn{s} by \eqn{s} matrix, and let
  \eqn{\bold{Q}_11} be the first \eqn{k} rows and columns,
  \eqn{\bold{Q}_22} be the last \eqn{s-k} rows and columns, and
  \eqn{\bold{Q}_{12}} and \eqn{\bold{Q}_{21}} the remaining two parts of
  the original matrix.  Let \eqn{SWP[1:k]\bold{Q}}
  be the result of sweeping out the first \eqn{k} rows and columns.
  (The sweeping does not need to be done in order, but it makes the
  description cleaner if it is).  Then

  \deqn{SWP[1:k]\bold{Q} = \left [ \begin{array}{cc}
    -\bold{Q}_{11}^{-1} & \bold{H}_{12} \\
    \bold{H}_{21} & \bold{Q}_{22.1}
    \end{array} \right ] }

  where \eqn{\bold{H}_{12} = \bold{Q}_{11}^{-1}\bold{Q}_{12}} and
  \eqn{\bold{Q}_{22.1} = \bold{Q}_{22} -
  \bold{Q}_{21}\bold{Q}_{11}^{-1}\bold{Q}_{12}}.  Note that if all rows
  and columns of a matrix are swept out, then this produces the negative
  of the inverse of the matrix.  (Several built in R routines for
  inverting a matrix, e.g., \code{\link[base]{solve}}, are available.
  The built-in R routines are probably more numberically stable and
  certainly have been better tested.)

  The sweep operator is associative, so that \eqn{SWP[1]SWP[2]\bold{Q} =
  SWP[1:2]\bold{Q}$ and communtative so that \eqn{SWP[2:1]\bold{Q} =
  SWP[1:2]\bold{Q}$.  The reverse sweep operator, \code{revSweep()} or
  \eqn{RSW[k]\bold{Q}} undoes the sweep operation.

  The properties of the sweep operator are more interesting when the
  matrix is augmented.  Let \eqn{\bold{X}} be an \eqn{n} by \eqn{p}
  matrix of observed data and let \eqn{\bold{X}_{(+)}} be \eqn{\bold{X}} 
  augmented by setting the \eqn{p+1}st column to 1.  Let
  \eqn{\bold{Q}_{(+)} = \bold{X}_{(+)}^T\bold{X}_{(+)}}.   Note that the
  last row and column of \eqn{\bold{Q}_{(+)}} is the column sums of
  \eqn{\bold{X}_{(+)}} with $q_{(+)p+1,p+1}=n$.  Sweeping on row and
  column \eqn{p+1} produces:

  \deqn{SWP[p+1]\bold{Q}_{(+)} = \left [ \begin{array}{cc}
    -\bold{T} & \bar \bold{X} \\
    \bar \bold{X} & -1/n
    \end{array} \right ] }

  where \eqn{\bold{T}} is \eqn{n-1} times the covariance matrix for
  \eqn{\bold{X}}. In general, sweeping out all rows columns except for
  \eqn{k} produces the regression of \eqn{X_k} on the remaining
  variables, with the values in the \eqn{k}th row and column the
  regression coefficients and the
  \eqn{SWP[c(1:(k-1),(k+1):(p+1))]\bold{Q}_{(+)}[k,k]} equal to the 
  residual variance.

  Let \eqn{\bold{Q}} = \bold{X}}^T\bold{X}}} and \eqn{\bold{X}_1} is the
  first \eqn{k} columns and \eqn{\bold{X}_2} the remaining columns, and
  let \eqn{SWP[1:k]\bold{Q}} be the partitioned matrix shown above.
  Then the least squares linear predictor of \eqn{\bold{X}_2} from
  \eqn{\bold{X}_1} is \eqn{\bold{H}_{12}\bold{X}_1}, and
  \eqn{\bold{Q}_{22.1}} is the residual covariance matrix.

  

}
\value{
  A square matrix of the same size as \code{A}.  If \code{A} is positive
  definete, then the output will be as well.
  
}
\references{

  Beaton, A. E. (1964).  The Use of Special Matrix Operators in
  Statistical Calculus.  Educational Testing Service Research Report
  RB-64-51. 

  Dempsters, A. P. (1969).  \emph{Elements of Continuous Multivariate
  Analysis.}  Addison-Wesley.
  
  Goodnight, J. H. (1979).  A Tutorial on the SWEEP Operator.
  \emph{The American Statistician.}  \bold{33} (3).  140--158.

  Little, R. J. A. and Rubin, D. B. (2002).  \emph{Statistical Analysis
    with Missing Data, Second Edition.}  Wiley.
  
}
\author{Russell ALmond}
\note{
  The version of the sweep operator described in Goodnight (1979) and
  the one in Dempster (1969) (also Little and Rubin, 2002) are slightly
  different.  This implementation follows the one in Dempster.  In fact
  the example calculations found in Dempster, p. 64-65 and 151-152 are
  available in the test scripts.  (There are some minor discrepancies
  between the values calculated here and given in Dempster, probably due
  to typographical or rounding errors.)
}
\seealso{
  \code{\link{matASM}}
}
\examples{

## Following example from Dempster (1969, p 151-2).
data(eggs)
eggsPlus <- cbind(eggs,1)
eggsQplus <- t(eggsPlus)\%*\%eggsPlus

## Sweep out constant row and column.
eggsTplus <- matSweep(eggsQplus,4)
stopifnot(all(abs(eggsTplus[1:3,1:3]-(nrow(eggs)-1)*cov(eggs)) < .0001),
          all(abs(eggsTplus[4,1:3]-apply(eggs,2,mean)) < .0001),
          abs(eggsTplus[4,4]+1/nrow(eggs)) < .0001)

## Predict log(Volume) from other measures.
eggsQ124 <- matSweep(eggsQplus,c(1:2,4))
lm.eggs <- lm(logVol~logLength+logWidth,data=as.data.frame(eggs))
stopifnot(all(abs(eggsQ124[3,c(4,1:2)] - coef(lm.eggs)) < .0001),
          abs(eggsQ124[3,3] - deviance(lm.eggs)) < .0001)


## Inversion using sweep
eggsQinv <- matSweep(eggsQplus,1:4)
stopifnot(all(abs(eggsQinv+solve(eggsQplus))<.0001))


## Commutative and associative properties.
eggsQ421 <- matSweep(matSweep(eggsQplus,4),2:1)
stopifnot(all(abs(eggsQ421-eggsQ124) <.0001))

## Reverse sweep operator
eggsQinv3 <- revSweep(eggsQinv,3)
stopifnot(all(abs(eggsQinv3-eggsQ124) <.0001))

}
\keyword{ algebra }% use one of  RShowDoc("KEYWORDS")
\keyword{ array }% use one of  RShowDoc("KEYWORDS")
