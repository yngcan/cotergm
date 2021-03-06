# Summary function for STERGM fits.

summary.stergm <- function (object, ...){
  out<-list(formation=summary(object$formation.fit,...),
            dissolution=summary(object$dissolution.fit,...))
  class(out)<-"summary.stergm"
  out
}

print.summary.stergm <- function(x, ...){
  cat("\n==============================\n")
    cat("Summary of formation model fit \n")
    cat("==============================\n\n")

  cat("Formula:   ")
  f <- x$formation$formula
  if(length(f)==3) f<-f[c(1,3)]
  print(f, showEnv=FALSE)
  cat("\n")

  print.summary.ergm(x$formation, ..., print.header=FALSE, print.formula=FALSE, print.degeneracy=FALSE, print.drop=FALSE, print.deviances=x$formation$estimate!="EGMME")

  cat("\n================================\n")
    cat("Summary of dissolution model fit\n")
    cat("================================\n\n")
  
  cat("Formula:   ")
  f <- x$dissolution$formula
  if(length(f)==3) f<-f[c(1,3)]
  print(f, showEnv=FALSE)
  cat("\n")

  print.summary.ergm(x$dissolution, ..., print.header=FALSE, print.formula=FALSE, print.degeneracy=FALSE, print.drop=FALSE, print.deviances=x$formation$estimate!="EGMME")

}
