\name{crossvaldata}
\alias{crossvaldata}
\title{
Computes k-fold cross validation for rminer models.
}
\description{
Computes k-fold cross validation for rminer models.
}
\usage{
crossvaldata(x, data, theta.fit, theta.predict, ngroup = n, 
             order = FALSE, model, task, feature = "none",
             ...)
}
\arguments{
  \item{x}{See \code{\link{fit}} for details.}
  \item{data}{See \code{\link{fit}} for details.}
  \item{theta.fit}{fitting function}
  \item{theta.predict}{prediction function}
  \item{ngroup}{number of folds}
  \item{order}{if \code{TRUE} then a static ordered sampling cross-validation is adopted (e.g. useful for time series data), else the normal random sampling is adopted.}
  \item{model}{See \code{\link{fit}} for details.}
  \item{task}{See \code{\link{fit}} for details.}
  \item{feature}{See \code{\link{fit}} for details.}
  \item{\dots}{Additional parameters sent to \code{theta.fit} (e.g. \code{search}, \code{mpar}, \code{scale}, \code{transform})}
}
\details{
Standard k-fold cross-validation but adopted for rminer models.
For classification tasks ("class" or "prob") a stratified sampling is used (the class distributions are identical for each fold).
}
\note{A better control (e.g. use of several Runs) is achieved using the simpler \code{\link{mining}} function.}

\value{
Returns a list with:
      \itemize{
               \item $cv.fit -- all predictions (factor if \code{task="class"}, matrix if \code{task="prob"} or numeric if \code{task="reg"});
               \item $mpar -- matrix with the mpar for each fold;
               \item $attributes -- the selected attributes for each fold if a feature selection algorithm was adopted;
               \item $ngroup -- the number of folds;
               \item $leave.out -- the computed size for each fold (=\code{nrow(data)/ngroup});
               \item $groups -- vector list with the indexes of each group;
               \item $call -- the call of this function;
   
              }
}
\references{
Check the \code{crossval} function of the bootstrap library.
}

\author{This function was adapted by Paulo Cortez from the \code{crossval} function of the bootstrap library (S original by R. Tibshirani and R port by F. Leisch).}

\seealso{
\code{\link{holdout}}, \code{\link{fit}}, \code{\link{mining}} and \code{\link{predict.fit}}.}
\examples{
 data(iris)
 M=crossvaldata(Species~.,iris,fit,predict,ngroup=3,model="mlpe",
                task="prob",search=4,mpar=c(3,100,"holdout",2/3,"AUC"))
}
\keyword{ models }
