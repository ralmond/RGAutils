#### This implements a factory paradigm for Models and Nodes.
## Warehouses are initialized with a default class and a manifest:
## a list of possible objects along with their metadata.  The factories
## maintain a cache of objects.  If the object is not found,
## metadata from the manifest is used to create it.

#########################################################
## Abstract Warehouse

Warehouse <- setRefClass("Warehouse",
                       fields=c(type="character",
                                manifest="data.frame",
                                inventory="list",
                                key="character",
                                packsep="character"),
                       prototype=list(type="UNDEFINED",
                                      manifest=data.frame(),
                                      inventory=list(),
                                      key="Name",
                                      packsep="/")
                       )
Warehouse$methods(
              ## Prototype doesn't work with reference classes, check for defaults.
              initialize = function(...,type="UNDEFINED",key="Name",
                                    packsep="/") {
                callSuper(...,type=type,key=key,packsep=packsep)
                })



Warehouse$methods(
            show = function () {
              cat("<Warehouse for ",type," with ",length(inventory),
                  " items in inventory>\n")
            },
            clear = function () {
              for (name in names(inventory)) {
                .self$free(name)
              }
              invisible(inventory)
            },
            manifestData = function (name) {
              name <- as.character(name)
              if (ispacked(name)) name <- unpackname(name)
              if (length(name) != length(key))
                stop("Expected name to contain elements",key)
              whch = rep(TRUE,nrow(manifest))
              for (i in 1:length(key)) {
                whch <- whch & manifest[[key[i]]] == name[i]
              }
              manifest[whch,,drop=FALSE]
            },
            supply = function (name) {
              val <- fetch(name)
              if (is.null(val))
                val <- make(name)
              val
            },
            packname = function (name) {
              ## This makes sure that the name is a single string.
              paste(name,collapse=packsep)
            },
            unpackname = function (name) {
              ## This undoes the previous operation
              strsplit(name,packsep,fixed=TRUE)[[1]]
            },
            ispacked = function (name) {
              length(name)==1 && grepl(packsep,name,fixed=TRUE)
            },
            isAvailable = function(name) {
              nrow(manifestData(name)) > 0L
            },
            fetch = function(name) {
              inventory[[packname(name)]]
            },
            make = function(name) {
              name <- as.character(name)
              pname <- packname(name)
              if (!is.null(inventory[[pname]]))
                free(name)
              dat = manifestData(name)
              val <- do.call(paste("Make",type,sep="."),
                             list(name,dat))
              inventory[[pname]] <<- val
              val
            },
            free = function(name) {
              pname <- packname(name)
              if (!is.null(inventory[[pname]]))
                Free(inventory[[pname]])
              inventory[[pname]] <<- NULL
              invisible(NULL)
            }
)

Make <- function(name,data) {
  stop("Use Warehouse$make instead.")
}

## Set this up as an S3 generic.
Free <- function (obj) {
  UseMethod("Free")
}
#setGeneric("Free") #Register it as an S4 generic, too.

## Default method is no special action required.
Free.default <- function (obj) {invisible(NULL)}


