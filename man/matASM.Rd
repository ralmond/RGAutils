\name{matASM}
\alias{matASM}
\title{Asimilates more rows into an already swept matrix.}
\description{

A swept matrix (see \code{\link{matSweep}}) is a correlation matrix that
has certain variables regressed out.  In particular, it represents a
regression of the variables corresponding to the unswept rows on
variables corresponding to the swept rows.  The \code{matASM} add
additional rows and columns corresponding to new variables, sweeping
those new rows and columns on the already swept variables.
}
\usage{
matASM(swpQ, Q123, Q33, ploc = 1L:nrow(swpQ))
}
\arguments{
  \item{swpQ}{
    A symetric matrix which is the result of a sweep operation.  In
    particular, it is assumed that this is the output of
    \code{matSweep(Q,ploc)}. 
  }
  \item{Q123}{
    A matrix of covariances where rows represent the existing variables,
    and columns the new variables.
  }
  \item{Q33}{
    A symmetric covariance matrix giving the covariance of the new variables.
  }
  \item{ploc}{
    A list of the rows on columns on which the first argument has been swept.
  }
}
\details{

  The asimilate operator is described by Dempster (1969, pp 67--70) as a
  complement to the sweep operator (\code{\link{matSweep}}).  Basically,
  it adds additional rows and columns corresponding to an already swept
  matrix corresponding to new variables added to a problem.
  
  Assume that the data matrix \eqn{\bold{X}} is divided into three sets
  of columns (corresponding to sets of variables):
  \eqn{\bold{X}_1}, \eqn{\bold{X}_2} and \eqn{\bold{X}_3}.  Let
  \eqn{\bold{Q} = \bold{X}^{T}\bold{X}} be partitioned as follows:

  \deqn{\bold{Q} = \left ( \begin{array}{ccc}
    \bold{Q}_{11} & \bold{Q}_{12} & \bold{Q}_{13} \\
    \bold{Q}_{21} & \bold{Q}_{22} & \bold{Q}_{23} \\
    \bold{Q}_{31} & \bold{Q}_{32} & \bold{Q}_{33}
    \end{array} \right)}

  Suppose that the original focus of the analysis was on the variables 
  \eqn{\bold{X}_1} and \eqn{\bold{X}_2}, and from the submatrix of
  \eqn{\bold{Q}} corresponding to those variables, the variables
  corresponding to \eqn{\bold{X}_1} have been swept out, leaving the
  matrix:
  
  \deqn{swp\bold{Q} = SWP[1]\left ( \begin{array}{cc}
    \bold{Q}_{11} & \bold{Q}_{12} \\
    \bold{Q}_{21} & \bold{Q}_{22} 
    \end{array} \right) =
    \left ( \begin{array}{cc}
    -\bold{Q}_{11}^{-1} & \bold{H}_{12} \\
    \bold{H}_{21} & \bold{Q}_{22.1}
    \end{array} \right)}

  The assimilate operator adds the variables corresponding to
  \eqn{\bold{X}_3} to \eqn{swp\bold{Q}}, sweeping out \eqn{\bold{X}_1}
  so that the result is \eqn{SWP[1]\bold{Q}}.

  The arguments of \code{matASM} are \code{swpQ} the result of sweeping
  out the variables in \code{ploc} from the original matrix.  \code{Q33}
  is the covariance matrix for the variables in \eqn{\bold{X}_3}, and
  \code{Q123} is the covariance matrix for the variables in
  \eqn{\bold{X}_1} and \eqn{\bold{X}_2} with \eqn{\bold{X}_3} (that is
  \eqn{\left(\atop{\bold{Q}_{13}}{\bold{Q}_{13}}\right)}).


}
\value{
  The value is a symmetric matric as described above.
  
}
\references{
  Dempsters, A. P. (1969).  \emph{Elements of Continuous Multivariate
  Analysis.}  Addison-Wesley.
}
\author{Russell G. Almond}
\note{
  Example from  Dempster(1969, p 69) is recreated in the test suite, but
  seems to differ in the least significant digits from what is printed
  in the book. 
}
\seealso{
  \code{\link{matSweep}}
}
\examples{
##### Example taken from Dempster 1969, p  64.

DQ <- matrix(c(19.1434,  9.0356,  9.7634, 3.2394,
                9.0356, 11.8658,  4.6232, 2.4746,
                9.7634,  4.6232, 12.2978, 3.8794,
                3.2394,  2.4746,  3.8794, 2.4604), 4,4)
## Sweep out row/column 1 from 2x2 submatrix
DQ2.1 <- matSweep(DQ[1:2,1:2],1)
## Assimilate row/columns 3 and 4
DQ2.1a34 <- matASM(DQ2.1,DQ[1:2,3:4],DQ[3:4,3:4],1)

## Sweep out row/column 1 from full matrix.
DQ234.1 <- matSweep(DQ,1)

## Should have gotten to the same place
stopifnot(all (abs(DQ234.1-DQ2.1a34) < .00001))

}
\keyword{ array }% use one of  RShowDoc("KEYWORDS")
\keyword{ algebra }
