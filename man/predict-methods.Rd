\name{predict.fit}
\alias{predict-methods}
\alias{predict.fit}
\alias{predict,model-method}
\title{predict method for fit objects (rminer)}
\description{predict method for fit objects (rminer)}

\section{Methods}{
\describe{

\item{\code{signature(object = "model")}}{ describe this method here }
}}


\arguments{

  \item{object}{a model object created by \code{\link{fit}}}
  \item{newdata}{a data frame or matrix containing new data}

}

\value{
   If \code{task} is \code{prob} returns a matrix, where each column is the class probability.\cr
   If \code{task} is \code{class} returns a factor.\cr
   If \code{task} is \code{reg} returns a numeric vector.\cr
   }

\seealso{
\code{\link{fit}}, \code{\link{mining}}, \code{\link{mgraph}}, \code{\link{mmetric}}, \code{\link{savemining}}, \code{\link{CasesSeries}}, \code{\link{lforecast}} and \code{\link{Importance}}.
}

\references{
\itemize{
\item To check for more details about rminer and for citation purposes:\cr
P. Cortez.\cr
Data Mining with Neural Networks and Support Vector Machines Using the R/rminer Tool.\cr
In P. Perner (Ed.), Advances in Data Mining - Applications and Theoretical Aspects 10th Industrial Conference on Data Mining (ICDM 2010), Lecture Notes in Artificial Intelligence 6171, pp. 572-583, Berlin, Germany, July, 2010. Springer. ISBN: 978-3-642-14399-1.\cr
@Springer: \url{http://www.springerlink.com/content/e7u36014r04h0334}\cr
\url{http://www3.dsi.uminho.pt/pcortez/2010-rminer.pdf}\cr
}
}

\examples{
## simple classification example with logistic regression
data(iris)
M=fit(Species~.,iris,model="lr")
P=predict(M,iris)
print(mmetric(iris$Species,P,"CONF")) # confusion matrix

## check also fit for more examples
}

\keyword{methods}
\keyword{classif}
\keyword{regression}
\keyword{neural}
\keyword{nonlinear}
