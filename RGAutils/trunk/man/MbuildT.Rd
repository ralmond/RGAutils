\name{EbuildT}
\alias{EbuildT}
\title{Builds expected values of sufficient statistics in multivariate
  normal regression}
\description{

  This does the E-step of an EM algorithm for multivariate regression.
  It is assumed that \code{X} is fully observed and \code{Y} may have
  missing data.  This returns an augmented matrix which contains the
  expected value of the means and the covariance matrix for the current
  set of parameters, \code{B}, \code{b0}, and \code{Syy.x}.


}
\usage{
MbuildT(X, Y, B, b0, Syy.x, w = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A data matrix with rows representing observations and columns
    variables.  These are independent variables in a regression, and
    must be fully observed.
  }
  \item{Y}{A data matrix with rows representing observations and columns
    variables.  These are dependent variables in a regression, and
    may contain missing values.  \code{X} and \code{Y} should have the
    same number of rows.
  }
  \item{B}{A \eqn{J \times K} matrix containing the non-constant
    coefficients of a multivariate regression of \code{Y} on \code{X}.
    Here \eqn{J} is the number of columns of \code{Y} and \eqn{K} is the
    number of columns of \code{X}.
  }
  \item{b0}{A vector of the constant terms from the multivariate
    regression of \code{Y} on \code{X}. This should have length
    \eqn{J}. 
  }
  \item{Syy.x}{A \eqn{J \times J} symmetric matrix giving the residual
    covariance of the multivariate regression of \code{Y} on \code{X}.
  }
  \item{w}{If supplied this should be a vector of weights equal to the
    number of rows of \code{X}.  The total sample size is considered to
    be the sum of the weights.
  }
}
\details{

  Following Little and Rubin (2002), the sufficient statistics for a
  multivariate normal distribution can be represented by a matrix:
  \deqn{ \bold{T} = \left [ \begin{array}{ccc}
    -1 & \mu_X & \mu_Y \\
    \mu_X & \Sigma_{XX} & \Sigma_{XY} \\
    \mu_Y & \Sigma_{YX} & \Sigma_{YY}
    \end{array} \right ] }.

  Sweeping (see \code{\link{matSweep}}) the rows and columns
  corresponding to the \eqn{\bold{X}} variables produces the
  multivariate regression of \eqn{\bold{Y}} on \eqn{\bold{X}}.
  
  \deqn{ SWP[X]\bold{T} = \left [ \begin{array}{ccc}
    * & * & b_{0.X} \\
    * & -\Sigma_{XX}^{-1} & B^T_{XY.X} \\
    b_{0.X} & B_{XY.X} & \Sigma_{YY.X}
    \end{array} \right ] },

  Here \eqn{B_{XY.X}} is the matrix of regression coefficients,
  \eqn{b_{0.X}} is the constants and \eqn{\Sigma_{YY.X}} is the residual
  covariances.

  If \code{Y} was fully observed, then \eqn{\bold{T}} could be easily
  calculated by first forming a matrix \code{XY1=cbind(1,X,Y)}, and then
  calculating \eqn{SWP[1](XY1)^T W (XY1)}, where \eqn{W} is a diagonal matrix
  with the weights, \code{w}, on the diagonal, and SWP[1] is the sweep
  operator applied to the constant row/column (assumed to be the first).
  Note that to get the parameter of the normal distribution, the part
  corresponding to the covariance matrix must be scaled by the sum of
  the weights.  

  If there are missing values, then calculating the expected value of
  \eqn{\bold{T}} requires two steps.  First, the missing values for
  \code{Y} must be imputed.  Second, an adjustment needs to be made for
  any term in the sum which involves the product of two imputed
  values (or the square of a single imputed value).

  In order to perform both the imputation and the adjustment, a 
  parameter estimate \eqn{bold{T}^{(i)}} is required.  In this
  function, this is supplied through the 


}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, Y, B, b0, Syy.x, w = 1) 
{
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    K <- ncol(X)
    J <- ncol(Y)
    sxT0 <- matrix(NA, K + J + 1, K + J + 1)
    totweight <- ifelse(length(w) == 1, nrow(X) * w, sum(w))
    Txx <- matSweep(SSX(X, w), 1)
    Txx[2:(K + 1), 2:(K + 1)] <- Txx[2:(K + 1), 2:(K + 1)]/totweight
    sxT0[1:(K + 1), 1:(K + 1)] <- matSweep(Txx, 2:(K + 1))
    sxT0[1, (K + 2):(K + J + 1)] <- b0
    sxT0[2:(K + 1), (K + 2):(K + J + 1)] <- t(B)
    sxT0[(K + 2):(K + J + 1), 1] <- b0
    sxT0[(K + 2):(K + J + 1), 2:(K + 1)] <- B
    sxT0[(K + 2):(K + J + 1), (K + 2):(K + J + 1)] <- Syy.x
    T0 <- revSweep(sxT0, 2:(K + 1))
    T <- matrix(0, K + J + 1, K + J + 1)
    XY1 <- cbind(1, X, Y)
    M <- is.na(Y)
    Mpat <- apply(sweep(M, 2, 2^(0:(ncol(M) - 1)), "*"), 1, sum)
    for (mp in unique(Mpat)) {
        if (as.integer(mp) == 0) 
            next
        whichobs <- which(mp == Mpat)
        whichvars <- 1 + K + which(M[whichobs[1], ])
        swpT0 <- matSweep(T0, setdiff(2:nrow(T0), whichvars))
        XY1[whichobs, whichvars] <- XY1[whichobs, -whichvars, 
            drop = FALSE] \%*\% swpT0[-whichvars, whichvars, drop = FALSE]
        sumwt <- ifelse(length(w) == 1, w * length(whichobs), 
            sum(w[whichobs]))
        T[whichvars, whichvars] <- T[whichvars, whichvars] + 
            sumwt * swpT0[whichvars, whichvars]
    }
    T <- crossprod(sweep(XY1, 1, w, "*"), XY1) + T
    T <- matSweep(T, 1)
    sumwt <- ifelse(length(w) == 1, nrow(XY1), sum(w))
    T[-1, -1] <- T[-1, -1]/sumwt
    list(T = T, XY1 = XY1, M = M)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
