\name{coef.stergm}
\alias{coef.stergm}
\alias{coefficients.stergm}
\title{Extract Model Coefficients}
\description{
  \code{coef} is a Method which extracts model coefficients
  from objects returned by the \code{\link{stergm}} function.  
  \code{coefficients} is an \emph{alias} for it.
}
\usage{
\method{coef}{stergm}(object, \dots)
\method{coefficients}{stergm}(object, \dots)
}
\arguments{
  \item{object}{A \code{\link{stergm}} fit.}
  \item{\dots}{other arguments.}
}
\value{
  Coefficients extracted from \code{object} in the form of a list with
  two elements: \code{formation}, a vector of formation
  coefficients and \code{dissolution}, a vector of dissolution
  coefficients.
}
\seealso{
  \code{\link{stergm}} and  \code{\link{ergm}}
}
\keyword{regression}
\keyword{models}
