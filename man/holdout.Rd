\name{holdout}
\alias{holdout}
\title{
Computes indexes for holdout data split into training and test sets.
}
\description{
Computes indexes for holdout data split into training and test sets.
}
\usage{holdout(y, ratio = 2/3, internalsplit = FALSE, mode = "random", iter = 1)
}
\arguments{
  \item{y}{desired target: numeric vector; or factor -- then a stratified holdout is applied (i.e. the proportions of the classes are the same for each set).}
  \item{ratio}{split ratio (in percentage -- sets the training set size; or in total number of examples -- sets the test set size).}
  \item{internalsplit}{if \code{TRUE} then the training data is further split into training and validation sets. The same \code{ratio} parameter is used for the internal split.}
  \item{mode}{sampling mode. Options are:
              \itemize{
                       \item \code{random} -- standard randomized holdout;
                       \item \code{order} -- static mode, where the first examples are used for training and the later ones for testing (useful for time series data);
                       \item \code{incremental} -- incremental retraining mode (e.g. useful for spam detection), similar to \code{order} except that \code{ratio}=\emph{batch} size
                                                   and \code{iter} is used. In each iteration, the training set size grows, while the test set size is equal to \code{ratio} except
                                                   for the last batch (where it may be smaller).
                      }
}
  \item{iter}{iteration of the incremental retraining mode (only used when \code{mode==}"incremental", typically \code{iter} is set within a cycle, see the example below).}
}
\details{
Computes indexes for holdout data split into training and test sets (if \code{y} is a factor then a stratified holdout is applied).
}
\value{
A list with the components:
\itemize{
         \item $tr -- numeric vector with the training examples indexes;
         \item $ts -- numeric vector with the test examples indexes;
         \item $itr -- numeric vector with the internal training examples indexes;
         \item $val -- numeric vector with the internal validation examples indexes;
        }
}
\references{
See \code{\link{fit}}.
}
\author{
Paulo Cortez \url{http://www3.dsi.uminho.pt/pcortez}
}

\seealso{
\code{\link{fit}}, \code{\link{predict.fit}}, \code{\link{mining}}, \code{\link{mgraph}}, \code{\link{mmetric}}, \code{\link{savemining}}, \code{\link{Importance}}.
}
\examples{
### simple examples:
H=holdout(1:10,ratio=2,internal=TRUE,mode="order")
print(H)
H=holdout(1:10,ratio=2/3,internal=TRUE,mode="order")
print(H)
H=holdout(1:10,ratio=2/3,internal=TRUE,mode="random")
print(H)
H=holdout(1:10,ratio=2/3,internal=TRUE,mode="random")
print(H)

### classification example
data(iris)
# random stratified holdout
H=holdout(iris$Species,ratio=2/3,internal=TRUE) 
print(summary(iris[H$itr,]))
print(summary(iris[H$val,]))
print(summary(iris[H$tr,]))
print(summary(iris[H$ts,]))
M=fit(Species~.,iris[H$tr,],model="dt") # training data only
P=predict(M,iris[H$ts,]) # test data
print(mmetric(iris$Species[H$ts],P,"CONF"))

### regression example with incremental training
ts=c(1,4,7,2,5,8,3,6,9,4,7,10,5,8,11,6,9)
d=CasesSeries(ts,c(1,2,3))
for(b in 1:3) # iterations
{
 H=holdout(d$y,ratio=4,mode="incremental",iter=b)
 print(H)
 M=fit(y~.,d[H$tr,],model="mlpe",search=2)
 P=predict(M,d[H$ts,])
 cat("batch :",b,"TR size:",length(H$tr),"TS size:",
     length(H$ts),"mae:",mmetric(d$y[H$ts],P,"MAE"),"\n")
}
}
\keyword{ manip }
