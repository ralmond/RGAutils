\name{matSweep}
\alias{matSweep}
\alias{revSweep}
\title{Beaton (Dempster) Sweep operator for Matrixes}
}
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
    \bold{Q}_{21}\bold{Q}_{11}^{-1}\bold{Q}_{12}}.

  
  

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
function (A, ploc = 1L:nrow(A), tol = .Machine$double.eps^0.5) 
{
    if (!is.matrix(A) || nrow(A) != ncol(A)) 
        stop("The first argument to matSweep must be a symmetric matrix")
    if (any(ploc < 1L) || any(ploc > nrow(A)) || any((ploc - 
        round(ploc)) > tol)) 
        stop("Pivot locations must be integers referencing array rows.")
    matSweepAux(A, ploc, tol)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
