### Methods for generating random categorical variables.

rcat <- function(n,prob,factor=TRUE,ordered=FALSE, labels=NULL) {
  if (missing(labels)) {
    if (is.matrix(prob)) {
      labels <- colnames(prob)
    } else {
      labels <- names(prob)
    }
  }
  if (is.data.frame(prob))
    prob <- as.matrix(prob)
  if (is.matrix(prob)) {
    kmax <- ncol(prob)
    cp <- t(apply(prob,1L,cumsum))
    if (any(abs(cp[,kmax]-1.0)>sqrt(.Machine$double.eps))) {
      stop("Probabilities don't sum to one.")
    }
  } else {
    kmax <- length(prob)
    cp <- cumsum(prob)
    if (abs(cp[kmax]-1)>sqrt(.Machine$double.eps)) {
      stop("Probabilities don't sum to one.")
    }
  }
  u <- runif(n)
  if (is.matrix(cp)) {
    x <- apply(sweep(cp,1L,u,"<"),1L,sum) + 1L
  } else {  # < and > because order of args is reversed
    x <- apply(outer(u,cp,">"),1L,sum) +1L
  }
  if (factor) {
    levels <- 1L:kmax
    if (is.null(labels)) {
      labels <- levels
    }
    x <- factor(x,levels,labels,ordered=ordered)
  }
  x
}

rocat <- function(n,prob,factor=TRUE,ordered=TRUE, labels=NULL) {
  if (missing(labels)) {
    if (is.matrix(prob)) {
      labels <- colnames(prob)
    } else {
      labels <- names(prob)
    }
  }
  rcat(n,prob,factor,ordered,labels)
}


## Also should have a dcat function
dcat <- function (x, prob, log=FALSE) {
  if (is.data.frame(prob))
    prob <- as.matrix(prob)
  if (is.matrix(prob)) {
    if (length(x) == 1L) {
      x <- rep(x,nrow(prob))
    }
    kmax <- ncol(prob)
    cp <- apply(prob,1L,sum)
    if (any(abs(cp-1.0)>sqrt(.Machine$double.eps))) {
      stop("Probabilities don't sum to one.")
    }
  } else {
    kmax <- length(prob)
    if (abs(sum(prob)-1.0) > sqrt(.Machine$double.eps)) {
      stop("Probabilities don't sum to one.")
    }
  }
  if (is.factor(x) && !is.null(colnames(prob))) {
    if (!all.equal(levels(x),colnames(prob))) {
      stop("Probability names do not match levels of factors.")
    }
  } else {
    if (!is.factor(x) && !is.integer(x) &&
        any(abs(x-round(x)) > sqrt(.Machine$double.eps), na.rm=TRUE)) {
      stop("x must be an integer or factor variable.")
    }
  }
  x <- as.integer(x)
  if (any (x < 1,na.rm=TRUE) || any(x>kmax,na.rm=TRUE)) {
    stop("x must be a value between 1 and ",kmax)
  }
  result <- numeric(length(x))
  if (is.matrix(prob)) {
    if (length(x) != nrow(prob)) {
      stop("Length of x does not match number of rows of prob")
    }
    for (i in 1L:length(x)) {
      result[i] <- prob[i,x[i]]
    }
  } else {
    result <- prob[x]
  }
  if (log) {
    result <- log(result)
  }
  result
}
docat <- function (x, prob, log=FALSE) {
    dcat(x,prob,log)
}

pocat <- function (q, prob, log=FALSE) {
  if (is.data.frame(prob))
    prob <- as.matrix(prob)
  if (is.matrix(prob)) {
    if (length(q) == 1L) {
      q <- rep(q,nrow(prob))
    }
    kmax <- ncol(prob)
    cp <- t(apply(prob,1L,cumsum))
    if (any(abs(cp[,kmax]-1.0)>sqrt(.Machine$double.eps))) {
      stop("Probabilities don't sum to one.")
    }
  } else {
    kmax <- length(prob)
    cp <- cumsum(prob)
    if (abs(cp[kmax]-1.0)>sqrt(.Machine$double.eps)) {
      stop("Probabilities don't sum to one.")
    }
  }
  if (is.factor(q) && !is.ordered(q)) {
    stop("Only well defined for ordered categorical variables.")
  }
  if (is.factor(prob) && !is.null(colnames(prob))) {
    if (!all.equal(levels(q),colnames(prob))) {
      stop("Probability names do not match levels of factors.")
    }
  } else {
    if (!is.factor(x) && !is.integer(x) &&
        any(abs(x-round(x)) > sqrt(.Machine$double.eps), na.rm=TRUE)) {
      stop("x must be an integer or factor variable.")
    }
  }
  q <- as.integer(q)
  if (any (q < 1,na.rm=TRUE) || any(q>kmax,na.rm=TRUE)) {
    stop("q must be a value between 1 and ",kmax)
  }
  result <- numeric(length(q))
  if (is.matrix(prob)) {
    if (length(q) != nrow(prob)) {
      stop("Length of q does not match number of rows of prob")
    }
    for (i in 1L:length(q)) {
      result[i] <- cp[i,q[i]]
    }
  } else {
    result <- cp[q]
  }
  if (log) {
    result <- log(result)
  }
  result
}

qocat <- function(p,prob,factor=TRUE,labels=NULL) {
  if (missing(labels)) {
    if (is.matrix(prob)) {
      labels <- colnames(prob)
    } else {
      labels <- names(prob)
    }
  }
  if (is.data.frame(prob))
    prob <- as.matrix(prob)
  if (is.matrix(prob)) {
    if (length(p) == 1L) {
      p <- rep(p,nrow(prob))
    }
    kmax <- ncol(prob)
    cp <- t(apply(prob,1L,cumsum))
    if (any(abs(cp[,kmax]-1.0)>sqrt(.Machine$double.eps))) {
      stop("Probabilities don't sum to one.")
    }
  } else {
    kmax <- length(prob)
    cp <- cumsum(prob)
    if (abs(cp[kmax]-1)>sqrt(.Machine$double.eps)) {
      stop("Probabilities don't sum to one.")
    }
  }
  if (!is.numeric(p) || any (p<0.0) || any(p>1.0)) {
    stop("p must be a numeric value between 0 and 1.")
  }
  if (is.matrix(cp)) {
    x <- apply(sweep(cp,1L,p,"<"),1L,sum) + 1L
  } else {# < and > because order of args is reversed
    x <- apply(outer(p,cp,">"),1L,sum) +1L
  }
  if (factor) {
    levels <- 1L:kmax
    if (is.null(labels)) {
      labels <- levels
    }
    x <- ordered(x,levels,labels)
  }
  x
}

## Median method for ordered factors.
## Follows the same syntax as median.
## Adding this to your workspace will cause median() to work on ordered
## variables.
median.ordered <- function (x, na.rm = FALSE) {
  result <- round(median(as.numeric(x),na.rm=na.rm))
  ordered(levels(x)[result],levels=levels(x))
}
