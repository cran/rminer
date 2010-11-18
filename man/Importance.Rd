\name{Importance}
\alias{Importance}
\title{
Measure input importance given a supervised data mining model.
}
\description{
Measure input importance given a supervised data mining model.
}
\usage{
Importance(M, data, RealL = 6, method = "sens", measure = "gradient", 
           sampling = "regular", baseline = "mean", responses = TRUE, 
           outindex = NULL, task = "default", PRED = NULL, 
           interactions = NULL)
}
\arguments{
  \item{M}{fitted model, typically is the object returned by \code{\link{fit}}. Can also be any fitted model (i.e. not from rminer), provided that PRED is defined.}
  \item{data}{training data (the same data.frame that was used to fit the model).}
  \item{RealL}{for numeric inputs, the number of sensitivity analysis levels (e.g. 6).}
  \item{method}{input importance method. Options are:
                \itemize{
                         \item sens -- sensitivity analysis
                         \item sensv -- equal to \code{sens} but sets \code{metric="variance"}.
                         \item sensg -- equal to \code{sens} but sets \code{metric="gradient"}.
                         \item sensr -- equal to \code{sens} but sets \code{metric="range"}.
                         \item sensi -- equal to \code{sens} but sets full interactions mode (Ith-D) sensitivity analysis (much slower than \code{sens}).
                         \item data -- uses all data for the sensitivity, much slower than \code{sens}.
                         \item randomforest -- uses method of Leo Breiman (type=1), only makes sense when M is a randomforest.
                        }
}
  \item{measure}{sensitivity analysis measure (used to measure input importance). Options are:
                \itemize{
                         \item gradient -- average absolute gradient (y_i+1-y_i) of the responses.
                         \item variance -- variance of the responses. 
                         \item range -- maximum - minimum of the responses. 
                        }
}
  \item{sampling}{for numeric inputs, the sampling scan function. Options are:
                \itemize{
                         \item  regular -- regular sequence (uniform distribution).
                         \item quantile -- sample values from the input that are more closer to the variable distribution in \code{data}. 
                        }
}
  \item{baseline}{baseline vector used during the sensitivity analysis. Options are:
                  \itemize{
                         \item mean -- uses a vector with the mean values of each attribute from \code{data}.
                         \item median -- uses a vector with the median values of each attribute from \code{data}.
                         \item a data.frame with the baseline example (should have the same attribute names as \code{data}). 
                        }
                 }
  \item{responses}{if \code{TRUE} then all sensitivity analysis responses are stored and returned.}
  \item{outindex}{the output index (column) of \code{data} if \code{M} is not a model object (returned by fit).}
  \item{task}{the \code{task} as defined in \code{\link{fit}} if \code{M} is not a model object (returned by fit).}
  \item{PRED}{the prediction function of \code{M}, if \code{M} is not a model object (returned by fit).}
  \item{interactions}{numeric vector with the attributes (columns) used by Ith-D sensitivity analysis (2-D or higher).}
}
\details{
This function provides several algorithms for measuring input importance of supervised data mining models.
A particular emphasis is given on sensitivity analysis (SA), which is a simple method that measures the effects on the output of a given model
when the inputs are varied through their range of values. Check the references for more details. 
In future documentation, the \code{sensi} and \code{interations} options will be explained, please ignore them in current use.
}
\value{
A \code{list} with the components:
\itemize{
         \item $value -- numeric vector with the computed sensitivity analysis measure for each attribute.
         \item $imp -- numeric vector with the relative importance for each attribute.
         \item $sresponses -- vector list as described in the Value documentation of \code{\link{mining}}.
        }
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
\author{
Paulo Cortez \url{http://www3.dsi.uminho.pt/pcortez}
}
\note{
See also \url{http://www3.dsi.uminho.pt/pcortez/rminer.html}
}

\seealso{
\code{\link{fit}}, \code{\link{mining}}, \code{\link{mgraph}}, \code{\link{mmetric}}, \code{\link{savemining}}.
}
\examples{
# typical use under rminer:
data(sin1reg) # x1 should account for 70%, x2 for 30% and x3 for 0%.
M=fit(y~.,sin1reg,model="svm")
I=Importance(M,sin1reg,method="sens",measure="gradient")
print(I)
L=list(runs=1,sen=t(I$imp),sresponses=I$sresponses)
mgraph(L,graph="IMP",leg=names(sin1reg),col="gray",Grid=10)
mgraph(L,graph="VEC",xval=1,Grid=10)

# theoretical sin1reg prediction:
mypred=function(M,data)
{ return (M[1]*sin(pi*data[,1]/M[3])+M[2]*sin(pi*data[,2]/M[3])) }
M=c(0.7,0.3,2000)
I=Importance(M,sin1reg,method="sens",measure="gradient",PRED=mypred,outindex=4)
print(I$imp) # x1=72.3% and x2=27.7%
L=list(runs=1,sen=t(I$imp),sresponses=I$sresponses)
mgraph(L,graph="IMP",leg=names(sin1reg),col="gray",Grid=10)
mgraph(L,graph="VEC",xval=1,Grid=10)
}
\keyword{ classif }
\keyword{ neural }
