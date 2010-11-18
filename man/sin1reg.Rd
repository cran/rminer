\name{sin1reg}
\alias{sin1reg}
\docType{data}
\title{sin1 regression dataset}
\description{
Simple synthetic dataset with 1000 points, where y=0.7*sin(pi*x1/2000)+0.3*sin(pi*x2/2000)
}
\usage{data(sin1reg)}
\format{
  The format is:
 chr "sin1reg"
}
\details{Simple synthetic dataset with 1000 points, where y=0.7*sin(pi*x1/2000)+0.3*sin(pi*x2/2000)}
\source{See references}
\references{
\itemize{
\item To cite the Importance function or sensitivity analysis method, please use:\cr
P. Cortez and M. Embrechts.\cr
Opening Black Box Data Mining Models Using Sensitivity Analysis.\cr
In Proceedings of the 2011 IEEE Symposium on Computational Intelligence and Data Mining (CIDM),
pp. 341-348, Paris, France, April, 2011.\cr
\url{http://www3.dsi.uminho.pt/pcortez/}\cr
}
}
\examples{
data(sin1reg)
print(summary(sin1reg))
}
\keyword{datasets}
