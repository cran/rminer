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
  \item{M}{fitted model, typically is the object returned by \code{\link{fit}}. 
           Can also be any fitted model (i.e. not from rminer), provided that the predict function PRED is defined (see examples for details).}
  \item{data}{training data (the same data.frame that was used to fit the model, currently only used to add data histogram to VEC curve).}
  \item{RealL}{for numeric inputs, the number of sensitivity analysis levels (e.g. 6). Note: you need to use \code{RealL}>=2.}
  \item{method}{input importance method. Options are:
                \itemize{
                         \item sens -- sensitivity analysis
                         \item sensv -- equal to \code{sens} but sets \code{measure="variance"}.
                         \item sensg -- equal to \code{sens} but sets \code{measure="gradient"}.
                         \item sensr -- equal to \code{sens} but sets \code{measure="range"}.
                         \item randomforest -- uses method of Leo Breiman (type=1), only makes sense when M is a randomRorest.
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
  \item{PRED}{the prediction function of \code{M}, if \code{M} is not a model object (returned by fit). 
              Note: this function should behave like the rminer \code{\link{predict-methods}}, i.e. return a numeric vector in case of regression;
a matrix of examples (rows) vs probabilities (columns) (\code{task="prob"}) or a factor (\code{task="class"}) in case of classification.
    }
  \item{interactions}{numeric vector with the attributes (columns) used by Ith-D sensitivity analysis (2-D or higher):
\itemize{
                      \item if \code{NULL} then only a 1-D sensitivity analysis is performed.
                      \item if \code{length(interactions)==1}?then a "special" 2-D sensitivity analysis is performed using the index of interactions versus
all remaining inputs. Note: the $sresponses[[interactions]] will be empty (in \code{\link{vecplot}} do not use \code{xval=interactions}).
                      \item if \code{length(interactions)>1}?then a full Ith-D sensitivity analysis is performed, where I=length(interactions).
Note: Computational effort can highly increase if I is too large, i.e. O(RealL^I). Also, you need to preprocess the 
returned list (e.g. using \code{avg_imp}) to use the \code{\link{vecplot}} function (see the examples).
}
                     }
}
\details{
This function provides several algorithms for measuring input importance of supervised data mining models.
A particular emphasis is given on sensitivity analysis (SA), which is a simple method that measures the effects on the output of a given model
when the inputs are varied through their range of values. Check the reference for more details.
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
\item To cite the Importance function or sensitivity analysis method, please use:\cr
P. Cortez and M. Embrechts.\cr
Opening Black Box Data Mining Models Using Sensitivity Analysis.\cr
In Proceedings of the 2011 IEEE Symposium on Computational Intelligence and Data Mining (CIDM),
pp. 341-348, Paris, France, April, 2011.\cr
\url{http://www3.dsi.uminho.pt/pcortez/}\cr
}
}
\author{
Paulo Cortez \url{http://www3.dsi.uminho.pt/pcortez}
}
\note{
See also \url{http://www3.dsi.uminho.pt/pcortez/rminer.html}
}

\seealso{
\code{\link{vecplot}}, \code{\link{fit}}, \code{\link{mining}}, \code{\link{mgraph}}, \code{\link{mmetric}}, \code{\link{savemining}}.
}
\examples{
### Typical use under rminer:
# 1st example, regression, 1-D sensitivity analysis
data(sin1reg) # x1 should account for 70%, x2 for 30% and x3 for 0%.
M=fit(y~.,sin1reg,model="svm")
I=Importance(M,sin1reg,method="sens",measure="gradient") # 1-D SA
print(I)
L=list(runs=1,sen=t(I$imp),sresponses=I$sresponses)
mgraph(L,graph="IMP",leg=names(sin1reg),col="gray",Grid=10)
mgraph(L,graph="VEC",xval=1,Grid=10,data=sin1reg) # or:
vecplot(I,xval=1,Grid=10,data=sin1reg,datacol="gray") # the same graph
vecplot(I,xval=c(1,2,3),pch=c(1,2,3),Grid=10,
leg=list(pos="bottomright",leg=c("x1","x2","x3"))) # all x1, x2 and x3 VEC curves

# 2nd example, regression, 2-D sensitivity analysis with 
# the most relevant input (x1, index 1):
I2=Importance(M,sin1reg,method="sensg",interactions=which.max(I$imp))
print(I2)
# influence of x1 and x2 over y
vecplot(I2,graph="VEC3",xval=2) # VEC surface
vecplot(I2,graph="VECC",xval=2) # VEC contour
# influence of x1 and x3 over y (influence of x3 is small random noise, 0%):
vecplot(I2,graph="VEC3",xval=3)
vecplot(I2,graph="VECC",xval=3)

# 3rd example, regression, full 3-D sensitivity analysis
I3=Importance(M,sin1reg,method="sensg",interactions=c(1,2,3))
print(I3)
I3_1d=avg_imp(I3,c(1)) # 1-D averaging under x1
vecplot(I3_1d,graph="VEC",xval=1,Grid=10)
I3_2d=avg_imp(I3,c(1,2)) # 2-D averaging under the pair x1,x2
vecplot(I3_2d,graph="VEC3")

### If you want to use Importance over your own model:
# 1st example, regression, uses the theoretical sin1reg function
mypred=function(M,data)
{ return (M[1]*sin(pi*data[,1]/M[3])+M[2]*sin(pi*data[,2]/M[3])) }
M=c(0.7,0.3,2000)
# 4 is the column index of y
I=Importance(M,sin1reg,method="sens",measure="gradient",PRED=mypred,outindex=4) 
print(I$imp) # x1=72.3% and x2=27.7%
L=list(runs=1,sen=t(I$imp),sresponses=I$sresponses)
mgraph(L,graph="IMP",leg=names(sin1reg),col="gray",Grid=10)
mgraph(L,graph="VEC",xval=1,Grid=10) # equal to:
vecplot(I,graph="VEC",xval=1,Grid=10)

# 2nd example, 3-class classification for iris and lda model:
data(iris)
library(MASS)
predlda=function(M,data) # the PRED function
{ return (predict(M,data)$posterior) }
LDA=lda(Species ~ .,iris, prior = c(1,1,1)/3)
# 4 is the column index of Species
I=Importance(LDA,iris,method="sensg",PRED=predlda,outindex=4)
vecplot(I,graph="VEC",xval=1,Grid=10,TC=1,
main="1-D VEC for Sepal.Lenght (x-axis) influence in setosa (prob.)")

# 3rd example, binary classification for setosa iris and lda model:
iris2=iris;iris2$Species=factor(iris$Species=="setosa")
predlda2=function(M,data) # the PRED function
{ return (predict(M,data)$class) }
LDA2=lda(Species ~ .,iris2)
I=Importance(LDA2,iris2,method="sensg",PRED=predlda2,outindex=4)
vecplot(I,graph="VEC",xval=1,
main="1-D VEC for Sepal.Lenght (x-axis) influence in setosa (class)",Grid=10)
}

\keyword{ classif }
\keyword{ neural }
