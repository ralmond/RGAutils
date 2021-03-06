\name{numericPart}
\alias{numericPart}
\alias{factorPart}
\title{Splits a mixed data frame into a numeric matrix and a factor part.}
\description{
  The function \code{numericPart()} converts a \code{data.frame} to a
  \code{matrix}, by dropping columns which contain non-numeric data.
  The function \code{factorPart} grabs the state information by
  selecting only columns which are factors.
}
\usage{
numericPart(table)
}
\arguments{
  \item{table}{A \code{data.frame} object.}
}
\details{

  The original purpose of this funciton was to split a conditional
  probability distribution in data frame format (with a set of factor
  rows identifying the states of the parent variables; see CPTtools package) to
  a table of just the numbers, and a data frame of just the factors so
  that they can be tackled separately.

}
\value{
  A matrix containing just the numeric columns of the data frame.
}
\author{Russell Almond}
\seealso{ \code{\link[base]{data.frame}}, \code{\link[base]{matrix}},
      \code{\link[base]{data.matrix}}}
\examples{
name <-c("Shahrazad", "Marguerite")
height <- c(34, 36)
weight <- c(28, 26)
twins <- data.frame(name=I(name),height=height, weight=weight)
numericPart(twins)
}
\keyword{utilities}
