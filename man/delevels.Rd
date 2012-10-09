\name{delevels}
\alias{delevels}
\title{
Reduce (delete) or replace levels from a factor variable (useful for preprocessing datasets).
}
\description{
Reduce (delete) or replace levels from a factor variable (useful for preprocessing datasets).
}
\usage{
delevels(x, levels, label = NULL)
}
\arguments{
  \item{x}{factor with several levels}
  \item{levels}{vector with the levels that will be replaced}
  \item{label}{the new label used for all \code{levels} examples (if \code{NULL} then \code{"_OTHER"} is assumed).}
}
\value{
Returns a factor with less levels.
}
\references{
See \code{\link{fit}}.
}
\author{
Paulo Cortez \url{http://www3.dsi.uminho.pt/pcortez}
}
\seealso{
\code{\link{fit}} and \code{\link{imputation}}.
}
\examples{
f=factor(c("A","A","B","B","C","D","E"))
print(table(f))
f1=delevels(f,"A","a")
print(table(f1))
f2=delevels(f,c("C","D","E"),"CDE")
print(table(f2))
f3=delevels(f,c("B","C","D","E"))
print(table(f3))
}
\keyword{ manip }
