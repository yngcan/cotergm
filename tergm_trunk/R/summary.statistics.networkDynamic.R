summary.statistics.networkDynamic <- function(object, at,..., basis=NULL){
  if(!is.null(basis)) object <- ergm.update.formula(object, basis~.)
  duration.dependent <- is.durational(object)
  t(rbind(sapply(at,
              function(t){
                nw <- network.extract.with.lasttoggle(ergm.getnetwork(object), t, duration.dependent)
                f <- ergm.update.formula(object, nw~., from.new="nw")
                summary(f,...)
              }
          )
      )
  )
}
