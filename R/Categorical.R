### Methods for generating random categorical variables.

rcat <- function(n,p,factor=TRUE,ordered=FALSE,
                 labels=ifelse(is.matrix(p),colnames(p),names(p))) {
  if (is.data.frame(p))
    p <- as.matrix(p)
  if (is.matrix(p)) {
    kmax <- ncol(p)
    cp <- t(apply(p,1,cumsum))
    if (any(cp[,kmax]>1)) {
      stop("Probabilities sum to greater than one.")
    }
  } else {
    kmax <- length(p)
    cp <- cumsum(p)
    if (cp[kmax]>1) {
      stop("Probabilities sum to greater than one.")
    }
  }
  u <- runif(n)
  if (is.matrix(cp)) {
    x <- apply(sweep(cp,1,u,">"),1,sum) + 1
  } else {
    x <- apply(outer(u,cp,"<"),1,sum) +1
  }
  if (factor) {
    levels <- 1:max(x)
    if (is.null(labels)) {
      labels <- levels
    }
    x <- factor(x,levels,labels,ordered=ordered)
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
