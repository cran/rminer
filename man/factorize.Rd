\name{factorize}
\alias{factorize}
\title{
Converts numeric object into a factor (levelling).
}
\description{
Converts numeric object into a factor (levelling).
}
\usage{
factorize(x, limits, labels)
}
\arguments{
  \item{x}{numeric vector.}
  \item{limits}{numeric vector with the limits [min, max[ for each label of the factor (note: \code{length(limits)} needs to be equal to 
               \code{length(labels)+1}).}
  \item{labels}{the new levels of the factor.}
}
\value{
Returns a factor.
}
\references{
See \code{\link{fit}}.
}
\author{
Paulo Cortez \url{http://www3.dsi.uminho.pt/pcortez}
}
\seealso{
\code{\link{fit}}, \code{\link{delevels}}, \code{\link{imputation}}.
}
\examples{
x=seq(-1,1,length.out=5)
f=factorize(x,c(-Inf,0,1e-99,Inf),c("neg","zero","pos"))
print(x)
print(f)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
