\name{mining}
\alias{mining}
\alias{medianminingpar}
\title{
Powerful function that trains and tests a particular fit model under several runs and a given validation method
}
\description{
Powerful function that trains and tests a particular fit model under several runs and a given validation method.
Since there can be a huge number of models, the fitted models are not stored. Yet, several useful statistics (e.g. predictions) are returned.
}
\usage{
mining(x, data = NULL, Runs = 1, method = NULL, model = "default", 
       task = "default", search = "heuristic", mpar = NULL,
       feature="none", scale = "default", transform = "none", 
       debug = FALSE, ...)
}
\arguments{
  \item{x}{a symbolic description (formula) of the model to be fit. If \code{x} contains the data, then \code{data=NULL} (similar to x in \code{\link{ksvm}}, kernlab package).}

  \item{data}{an optional data frame (columns denote attributes, rows show examples) containing the training data, when using a formula.}

  \item{Runs}{number of runs used (e.g. 1, 5, 10, 20, 30)}

  \item{method}{a vector with c(\emph{vmethod},\emph{vpar}), where \emph{vmethod} is:
                \itemize{
                         \item \code{all} -- all \emph{NROW} examples are used as both training and test sets (no \emph{vpar} is needed).
                         \item \code{holdout} -- standard holdout method. If \emph{vpar}<1 then \emph{NROW}*vpar random samples are used for training and the remaining rows are used for testing. Else, then \emph{NROW}*vpar random samples are used for testing and the remaining are used for training. For classification tasks (\code{prob} or \code{class}) a stratified sampling is used.
                         \item \code{holdoutorder} -- similar to \code{holdout} except that instead of a random sampling, the first rows (until the split) are used for training and the remaining ones for testing (equal to \code{mode="order"} in \code{\link{holdout}}).
                         \item \code{holdoutinc} -- incremental holdout retraining (e.g. used for spam data). Here, \emph{vpar} is the batch size.
                         \item \code{kfold} -- K-fold cross-validation method, where \emph{vpar} is the number of folds.
                         \item \code{kfoldo} -- similar to \code{kfold} except that instead of a random sampling, the order of the rows is used to build the folds.
                        }
               }
  \item{model}{See \code{\link{fit}} for details.}

  \item{task}{See \code{\link{fit}} for details.}

  \item{search}{See \code{\link{fit}} for details.}

  \item{mpar}{See \code{\link{fit}} for details.\cr
             }

  \item{feature}{See \code{\link{fit}} for more details about \code{feature="none"}, \code{"sabs"} or \code{"sbs"} options.\cr
                 For the \code{mining} function, additional options are \code{feature=}\emph{fmethod}, where \emph{fmethod} can be one of:
                 \itemize{ 
                           \item \code{sens} or \code{sensg} -- compute the 1-D sensitivity analysis input importances (\code{$sen}), gradient measure. 
                           \item \code{sensv} -- compute the 1-D sensitivity analysis input importances (\code{$sen}), variance measure.
                           \item \code{sensr} -- compute the 1-D sensitivity analysis input importances (\code{$sen}), range measure.
                           \item \code{simp}, \code{simpg} or \code{s} -- equal to \code{sensg} but also computes the 1-D sensitivity responses (\code{$sresponses}, 
                                 useful for \code{graph="VEC"}).
                           \item \code{simpv} -- equal to \code{sensv} but also computes the 1-D sensitivity responses (useful for \code{graph="VEC"}).
                           \item \code{simpr} -- equal to \code{sensr} but also computes the 1-D sensitivity responses (useful for \code{graph="VEC"}).
                      }
                }

  \item{scale}{See \code{\link{fit}} for details.}

  \item{transform}{See \code{\link{fit}} for details.}

  \item{debug}{If TRUE shows some information about each run.}

  \item{\dots}{See \code{\link{fit}} for details.}

}

\details{
Powerful function that trains and tests a particular fit model under several runs and a given validation method
(see [Cortez, 2010] for more details).\cr
Several \code{Runs} are performed. In each run, the same validation method is adopted (e.g. \code{holdout}) and
several relevant statistics are stored. Warning: be patient, this function can require some computational effort, specially if
a high number of \code{Runs} is used.
}

\value{
A \code{list} with the components:
\itemize{
\item $time -- vector with time elapsed for each run.
\item $test -- vector list, where each element contains the test (target) results for each run.
\item $pred -- vector list, where each element contains the predicted results for each test set and each run.
\item $error -- vector with an error \code{metric} for each run (the error depends on the \code{metric} parameter of \code{mpar}, valid options are explained in \code{\link{mmetric}}).
\item $mpar -- data.frame with each fit model mpar parameters, the sequence repeats \code{Runs} (times \emph{vpar} if \code{kfold} is used).
\item $model -- the \code{model}.
\item $task -- the \code{task}.
\item $method -- the external validation \code{method}.
\item $sen -- a matrix with the 1-D sensitivity analysis input importances. The number of rows is \code{Runs} times \emph{vpar}, if \code{kfold}, else is \code{Runs}.
\item $sresponses -- a vector list with a size equal to the number of attributes (useful for \code{graph="VEC"}). 
Each element contains a list with the 1-D sensitivity analysis input responses
(\code{n} -- name of the attribute; \code{l} -- number of levels; \code{x} -- attribute values; \code{y} -- 1-D sensitivity responses.\cr
                  Important note: sresponses (and "VEC" graphs) are only available if \code{feature="sabs"} or \code{"simp"} related (see \code{feature}).
\item $runs -- the \code{Runs}.
\item $attributes -- vector list with all attributes (features) selected in each run (and fold if \code{kfold}) if a feature selection algorithm is used.
\item $feature -- the \code{feature}.
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
\code{\link{fit}}, \code{\link{predict.fit}}, \code{\link{mgraph}}, \code{\link{mmetric}}, \code{\link{savemining}}, \code{\link{holdout}} and \code{\link{Importance}}.
}

\examples{
### simple regression example
x1=rnorm(200,100,20); x2=rnorm(200,100,20)
y=0.7*sin(x1/(25*pi))+0.3*sin(x2/(25*pi))
M=mining(y~x1+x2,Runs=2,model="mlpe",search=2)
print(M)
print(mmetric(M,metric="MAE"))

### classification example (task="prob")
data(iris)
M=mining(Species~.,iris,Runs=10,method=c("kfold",3),model="dt")
print(mmetric(M,metric="CONF"))
print(mmetric(M,metric="AUC"))
print(meanint(mmetric(M,metric="AUC")))
mgraph(M,graph="ROC",TC=2,baseline=TRUE,Grid=10,leg="Versicolor",
       main="versicolor ROC")
mgraph(M,graph="LIFT",TC=2,baseline=TRUE,Grid=10,leg="Versicolor",
       main="Versicolor ROC")
M2=mining(Species~.,iris,Runs=10,method=c("kfold",3),model="svm")
L=vector("list",2)
L[[1]]=M;L[[2]]=M2
mgraph(L,graph="ROC",TC=2,baseline=TRUE,Grid=10,leg=c("DT","SVM"),main="ROC")

### regression example
data(sin1reg)
M=mining(y~.,data=sin1reg,Runs=3,method=c("holdout",2/3),model="mlpe",
         search="heuristic5",mpar=c(50,3,"kfold",3,"MAE"),feature="sabs")
print(mmetric(M,metric="MAE"))
print(M$mpar)
cat("median H nodes:",medianminingpar(M)[1],"\n")
print(M$attributes)
mgraph(M,graph="RSC",Grid=10,main="sin1 MLPE scatter plot")
mgraph(M,graph="REP",Grid=10,main="sin1 MLPE scatter plot",sort=FALSE)
mgraph(M,graph="REC",Grid=10,main="sin1 MLPE REC")
mgraph(M,graph="IMP",Grid=10,main="input importances",xval=0.1,leg=names(sin1reg))
mgraph(M,graph="VEC",Grid=10,main="x1 VEC curve",xval=1,leg=names(sin1reg)[1])

### another classification example
data(iris)
M=mining(Species~.,data=iris,Runs=2,method=c("kfold",2),model="svm",search="heuristic",mpar=c(NA,NA,"kfold",3,"AUC"),feature="s")
print(mmetric(M,metric="AUC",TC=2))
mgraph(M,graph="ROC",TC=2,baseline=TRUE,Grid=10,leg="SVM",main="ROC",intbar=FALSE)
mgraph(M,graph="IMP",TC=2,Grid=10,main="input importances",xval=0.1,leg=names(iris),axis=1)
mgraph(M,graph="VEC",TC=2,Grid=10,main="Petal.Width VEC curve",data=iris,xval=4)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif}
\keyword{regression}
\keyword{neural}
\keyword{nonlinear}
