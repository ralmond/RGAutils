\name{dgetFromString}
\alias{dgetFromString}
\alias{dputToString}
\title{ Serializes an R object to a string}
\description{
  The function \code{dputToString} converts an R object to a string
  which can then be turned back into an R object using
  \code{dgetFromString}. 
}
\usage{
dgetFromString(str)
dputToString(obj)
}
\arguments{
  \item{str}{A string containing a serialized object}
  \item{obj}{An object to be serialized}
}
\details{
  These functions call the base R functions \code{\link[base]{dget}} and
  \code{\link[base]{dput}} using a string buffer as the connection.
  Thus, they serialize the R object and return a string value which can
  be stored in a database, or foreign object.    Note that the object
  must be self-contained, as links to other objects will not be properly
  recreated. 

  These are used in the RNetica package to create fields which are
  essentially R object (see \code{\link[RNetica]{NodeUserObj}})
  and \code{\link[RNetica]{NetworkUserObj}}).

}
\value{
  The function \code{dputToString} returns a character scalar containing
  the serialized object.  Note:  Sometimes R \dQuote{helpfully} adds
  line breaks, returning a vector of strings.  This can be fixed by
  using \code{paste(dputToString(obj),collapse=" ")}.
  

  The function \code{dgetFromString} returns an arbitrary R object
  depending on what was stored in \code{str}.

}
\author{Russell Almond}
\seealso{
  \code{\link[base]{dget}}, \code{\link[base]{dput}},
  Uses in the RNetica package:
  \code{\link[RNetica]{NodeUserObj}}), \code{\link[RNetica]{NetworkUserObj}}
}
\examples{
x <- sample(1L:10L)

x1 <- dgetFromString(dputToString(x))

stopifnot(all(x==x1))

}
\keyword{ manip }

