\name{Categorical}
\alias{dcat}
\alias{rcat}
\alias{docat}
\alias{pocat}
\alias{qocat}
\alias{rocat}
\title{Categorical and Ordered Categorical Distributions.
}
\description{

  A categorical distribution is one where the random variable can take
  one of a finite number of values.  The random variable can take on a
  nominal or ordinal scale.  The random variable can be represented by a
  (ordered) factor or an integer.  For unordered categories, only
  the \code{dcat} and \code{rcat} operations are supported.  For ordered
  categories, \code{pocat} and \code{qocat} are supported as well.

}
\usage{
dcat(x, prob, log = FALSE)
docat(x, prob, log = FALSE)
pocat(q, prob, log = FALSE)
qocat(p, prob, factor = TRUE, labels = NULL)
rcat(n, prob, factor = TRUE, ordered = FALSE, labels = NULL)
rocat(n, prob, factor = TRUE, labels = NULL)
}
\arguments{
  \item{n}{An integer scalar giving the number of random values to generate.}
  \item{p}{A numeric vector of probabilities.}
  \item{q}{An integer or ordered factor giving quantiles of the distribution.}
  \item{x}{An integer, factor or ordered factor giving values of the
    random variable.}
  \item{prob}{A vector, matrix or data frame of probability values.  All
    rows must add to one.  The number of rows should match \code{n} or the
    length of \code{p}, \code{q} or \code{x} (if those are note one).}
  \item{factor}{A logical scalar, indicates whether or not the output
    should be converted into a factor.}
  \item{ordered}{A logical scalar, indicates whether or not the output
    should be converted into an ordered factor.}
  \item{labels}{A character vector.  If the output is an ordered factor,
    these are the names of the levels.  The default is the names or
    column names of \code{prob}.}
  \item{log}{A logical scalar.  If true, log probabilities are returned
    instead of probabilities.}
}
\details{

  A categorical distribution is descrbed by a vector of probabilities,
  \eqn{p_1,\ldots,p_k} which sum to one.  The possible values are the
  integers \eqn{1,\ldots,k}.  The parameter \code{prob} is represented
  by a numeric vector whose values sum to one.  The random value
  \code{x} can be represented either by an integer, or a factor value.
  If it is a factor value, then the level names of the factor should
  match the names of \code{prob}.  The value of \code{prob} can be a
  matrix or a data frame.  In which case each row is treated as a
  parameter.  The column names of the matrix (or names of the data
  frame) should match the level names of the factor variable.

  An ordered categorical distribution is a categorical distribution
  whose values are considered to be ordered.  The while all categorical
  variables are on at least a nominal scale, ordered categorial
  variables are also considered to be on an interval scale.  Note that
  quantiles are defined for ordered categorical variables, but not for
  unordered categorical variables, thus \code{pocat(q,prob)} and
  \code{qocat(p,prob)} functions are available, but not \code{pcat(q,prob)} and
  \code{qcat(p,prob)}

  The functions \code{dcat(x,prob)} and \code{docat(x,prob)} calculate the
  probability of the random values \code{x}.  The two functions are
  identical, as ordered and unordered categorical variables behave the
  same for this operation.  The function is vectorized; if \code{x} is a
  vector or \code{prob} is a matrix or data frame, then the result is a
  vector of probabilities.  

  The functions \code{rcat(n,prob)} and \code{rocat(n,prob)} both
  generate a random vector of categorical values.  The only difference
  between the functions is that \code{rocat} produces an ordered factor
  by default. If \code{prob} is a matrix or data frame, then the number
  of rows should be \code{n} and a different set of probabilities is
  used for each value.

  The function \code{pocat(q,prob)} calculate the cumulative probability
  of the quantile \code{q}.  This operation is only defined for ordered
  categorical variables.  The function is vectorized; if \code{q} is a
  vector or \code{prob} is a matrix or data frame, then the result is a
  vector of probabilities.

  The function \code{pocat(p,prob)} calculate the quantile corresponding
  to a given probability \code{p}.  This operation is only defined for ordered
  categorical variables.  The function is vectorized; if \code{p} is a
  vector or \code{prob} is a matrix or data frame, then the result is a
  vector of values.

}
\value{
  For \code{rcat}, \code{rocat} and \code{qocat} the value is a vector
  of values.  These may be integers or (ordered) factors depending on
  the value of the \code{factor} argument.

  For \code{dcat}, \code{docat} and \code{pocat} the value is a vector
  of probabilities.

}
\author{Russell G Almond}
\seealso{
  \code{\link{median.ordered}}
}
\examples{

### Random Number Test
set.seed(123456)
pp <- c(L=.5,M=.3,H=.2)
x <- rcat(1000,pp)
stopifnot(is.factor(x),!is.ordered(x),
          all.equal(levels(x),names(pp)))
exp <- 1000*pp
chisq <- sum((table(x)-exp)^2/exp)
stopifnot (chisq < qchisq(.95,2))

## Ordered categorical
x1 <- rocat(10,pp)
stopifnot(is.factor(x1),is.ordered(x1),
          all.equal(levels(x1),names(pp))

### Discrete Markov Chain Generation
N <- 100
Tmax <- 10
P0 <- c(L=.25,M=.5,H=.25)
P1 <- rbind(L=c(L=.6,M=.3,H=.1),
            M=c(L=.2,M=.6,H=.2),
            H=c(L=.1,M=.3,H=.6))
chain <- matrix(ordered(NA,levels=1:3,labels=c("L","M","H")),N,Tmax)
chain[,1] <- rcat(N,P0)
for (t in 2:Tmax) {
  chain[,t] <- rcat(N,P1[as.integer(chain[,t-1L]),])
}


dd <- dcat(3:1,pp)
stopifnot(all.equal(dd,rev(pp)))
dd <- dcat(factor(3:1,levels=1:3,labels=names(pp)),pp)
stopifnot(all.equal(dd,rev(pp)))

ddd <- dcat(1:3,P1)
stopifnot(all.equal(ddd,diag(P1),check.attributes=FALSE))
ddd <- dcat(factor(1:3,levels=1:3,labels=colnames(P1)),P1)
stopifnot(all.equal(ddd,diag(P1),check.attributes=FALSE))

pp1 <- pocat(1:3,pp)
stopifnot(all.equal(pp1,cumsum(pp),check.attributes=FALSE))
pp1 <- pocat(ordered(3:1,levels=1:3,labels=names(pp)),pp)
stopifnot(all.equal(pp1,rev(cumsum(pp)),check.attributes=FALSE))

pp1 <- pocat(1:3,P1)
stopifnot(all.equal(pp1,diag(t(apply(P1,1,cumsum))),
          check.attributes=FALSE))


qq1 <- qocat(seq(.1,1,.1),pp)
stopifnot(all.equal(qq1,ordered(rep(1:3,times=10*pp),
                                levels=1:3,labels=names(pp))))
qq2 <- qocat(.5,P1)
stopifnot(all.equal(qq2,ordered(1:3,levels=1:3,labels=names(pp)),
          check.attributes=FALSE))


}



\keyword{ distribution }

