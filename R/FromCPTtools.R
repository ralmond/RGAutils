"numericPart" <-
function(table) {
  which <- sapply(table,is.numeric)
  as.matrix(table[,which])
}

"factorPart" <-
function(table) {
  which <- sapply(table,is.factor)
  table[,which]
}


colorspread <- function(col,steps,maxsat=FALSE,rampval=FALSE) {
  hsvmat <- rgb2hsv(col2rgb(rep(col,steps)))
  if (maxsat) {
    hsvmat["s",] <- 1
  }
  hsvmat["s",] <- hsvmat["s",]*(1:steps)/steps
  if (rampval) {
    hsvmat["v",] <- hsvmat["v",]*(steps:1)/steps
  }
  hsv(hsvmat["h",],hsvmat["s",],hsvmat["v",])
}


## These methods support general serialization to string buffers.
dputToString <- function (obj) {
  con <- textConnection(NULL,open="w")
  tryCatch({dput(obj,con);
           textConnectionValue(con)},
           finally=close(con))
}

dgetFromString <- function (str) {
  con <- textConnection(str,open="r")
  tryCatch(dget(con), finally=close(con))
}


