
###########################################################
## Pnet Factory

### Create table of Model meta-information
BuildNetManifest <- function (Pnetlist) {
  ## Node Level fields
  Name <- character()
  Hub <- character()
  Title <- character()
  Pathname <- character()
  Description <- character()
  ## Main Loop
  for (net in Pnetlist) {
    if (!is.Pnet(net)) {
      stop("Expected a list of Pnets, got ",net)
    }
    Name <-c(Name,PnetName(net))
    Title <- c(Title,PnetTitle(net))
    Hub <- c(Hub,PnetHub(net))
    Pathname <-c(Pathname,PnetPathname(net))
    Description <- c(Description,PnetDescription(net))
  }
  result <- data.frame(Name,Title,Hub,Pathname,Description,
                       stringsAsFactors=FALSE)
  rownames(result) <- Name
  ## Seems there is a bug in the class checking mechanism,
  ## Easiest to just leave it as a data.frame
  ##class(result) <- c("NetManifest",class(result))
  result
}

PnetWarehouse <- setRefClass("PnetWarehouse",contains="Warehouse")


PnetWarehouse$methods(
            show = function () {
              cat("<Pnet Warehouse for ",type," with ",length(inventory),
                  " items in inventory\n>")
            },
            save = function (name,pathname) {
              if (missing(pathname)) {
                pathname <- manifestData(name)$Pathname
              }
              Save(inventory[[packname(name)]],pathname)
            },
            delete = function (name) {
              pname <- packname(name)
              if (!is.null(inventory[[pname]])) {
                Delete(inventory[[packname(name)]])
              }
              inventory[[pname]] <<- NULL
            },
            reload = function (name,pahtname,fromScratch=FALSE) {
              if (fromScratch) {
                new <- make(name)
              } else {
                if (missing(pathname)) {
                  pathname <- manifestData(name)$Pathname
                }
                free(name)
                new <- Reload(inventory[[packname(name)]],pathname)
                inventory[[packname(name)]] <<- new
              }
              invisible(new)
            })

Save<-function(net,pathname) {
  UseMethod("Save")
}
#setGeneric("Save")


Reload<-function(net,pathname)
    UseMethod("Reload")
#setGeneric("Reload")

Delete<-function(obj)
  UseMethod("Delete")
#setGeneric("Delete")
