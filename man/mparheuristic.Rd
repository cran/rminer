\name{mparheuristic}
\alias{mparheuristic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that returns a list of searching (hyper)parameters for a particular classification or regression model
}
\description{
Function that returns a list of searching (hyper)parameters for a particular classification or regression model.
The result is to be put in a \code{search} argument, used by \code{\link{fit}} or \code{\link{mining}} functions. Something
like: \code{search=list(search=mparheuristic(...),...)}.
}
\usage{
mparheuristic(model, n = NA, lower = NA, upper = NA, by = NA, kernel = "rbfdot")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{model type name. See \code{\link{fit}} for details.
} 
  \item{n}{number of searches (either \code{n} or \code{by} should be used, \code{n} has prevalence over \code{by}).
}
  \item{lower}{lower bound for the (hyper)parameter (if \code{NA} a default value is assumed).
}
  \item{upper}{upper bound for the (hyper)parameter (if \code{NA} a default value is assumed).
}
  \item{by}{increment in the sequence (if \code{NA} a default value is assumed depending on \code{n}).
%%     ~~Describe \code{by} here~~
}
  \item{kernel}{optional kernel type, only used when \code{model="ksvm"}. Currently mapped kernels are \code{"rbfdot"}, \code{"polydot"} and
\code{"vanilladot"}; see \code{\link{ksvm}} for kernel details.
%%     ~~Describe \code{kernel} here~~
}
}
\details{
This function facilitates the definition of the \code{search} argument used by \code{\link{fit}} or \code{\link{mining}} functions.
Using simple heuristics, reasonable (hyper)parameter search values are suggested for several rminer models. For models not
mapped in this function, the function returns \code{NULL}, which means that no hyperparameter search is executed (often,
this implies using rminer or R function default values).

The heuristic assumes lower and upper bounds for a (hyper)parameter. If \code{n=1}, then rminer or R defaults are assumed.
Else, a search is created using \code{seq(lower,upper,by)}, where \code{by} was set by the used or computed from \code{n}.
For \code{model="ksvm"}, \code{2^seq(...)} is used for sigma and C, \code{(1/10)^seq(...)} is used for scale.
}
\value{
A list with one ore more (hyper)parameter values to be searched.
}
\references{
\itemize{
\item To check for more details about rminer and for citation purposes:\cr
P. Cortez.\cr
Data Mining with Neural Networks and Support Vector Machines Using the R/rminer Tool.\cr
In P. Perner (Ed.), Advances in Data Mining - Applications and Theoretical Aspects 10th Industrial Conference on Data Mining (ICDM 2010), Lecture Notes in Artificial Intelligence 6171, pp. 572-583, Berlin, Germany, July, 2010. Springer. ISBN: 978-3-642-14399-1.\cr
@Springer: \url{http://www.springerlink.com/content/e7u36014r04h0334}\cr
\url{http://www3.dsi.uminho.pt/pcortez/2010-rminer.pdf}\cr

\item This tutorial shows additional code examples:\cr
P. Cortez.\cr
A tutorial on using the rminer R package for data mining tasks.\cr
Teaching Report, Department of Information Systems, ALGORITMI Research Centre, Engineering School, University of Minho, Guimaraes, 
Portugal, July 2015.\cr
\url{http://hdl.handle.net/1822/36210}\cr

\item Some bounds were retrieved from:\cr
M. Fernandez-Delgado, E. Cernadas, S. Barro and D. Amorim.
Do we need hundreds of classifiers to solve real world classification problems?. 
In The Journal of Machine Learning Research, 15(1), 3133-3181, 2014.
}

}
\author{
Paulo Cortez \url{http://www3.dsi.uminho.pt/pcortez}
}
\note{
See also \url{http://hdl.handle.net/1822/36210} and \url{http://www3.dsi.uminho.pt/pcortez/rminer.html}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fit}} and \code{\link{mining}}.
}
\examples{
## "kknn"
s=mparheuristic("kknn",n=1)
print(s)
s=mparheuristic("kknn",n=10)
print(s)
s=mparheuristic("kknn",lower=5,upper=15,by=2)
print(s)

## "mlpe"
s=mparheuristic("mlpe")
print(s) # "NA" means set size with inputs/2 in fit
s=mparheuristic("mlpe",n=10)
print(s) 

## "randomForest"
s=mparheuristic("randomForest",n=10)
print(s) 

## "ksvm"
s=mparheuristic("ksvm",n=10)
print(s) 
s=mparheuristic("ksvm",n=10,kernel="vanilladot")
print(s) 
s=mparheuristic("ksvm",n=10,kernel="polydot")
print(s) 

## "rpart" and "ctree" are special cases (see help(fit,package=rminer) examples):
s=mparheuristic("rpart",n=3)
print(s) 
s=mparheuristic("ctree",n=3)
print(s) 

### examples with fit
\dontrun{
### classification
data(iris)
s=mparheuristic("ksvm",n=3,kernel="vanilladot")
print(s)
search=list(search=s,method=c("holdout",2/3,123))
M=fit(Species~.,data=iris,model="ksvm",search=search,fdebug=TRUE)
print(M@mpar)

### regression
data(sa_ssin)
s=mparheuristic("ksvm",n=3,kernel="polydot")
print(s)
search=list(search=s,metric="MAE",method=c("holdout",2/3,123))
M=fit(y~.,data=sa_ssin,model="ksvm",search=search,fdebug=TRUE)
print(M@mpar)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif}
\keyword{regression}
\keyword{neural}
\keyword{nonlinear}
