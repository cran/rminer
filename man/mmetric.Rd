\name{mmetric}
\alias{mmetric}
\alias{metrics}
\title{
Compute classification or regression error metrics.
}
\description{
Compute classification or regression error metrics.
}
\usage{mmetric(y, x = NULL, metric, D = 0.5, TC = -1, val = NULL)}

\arguments{
  \item{y}{if there are predictions (\code{!is.null(x)}), \code{y} should be a numeric vector or factor with the target desired responses (or output values).\cr
           Else, \code{y} should be a list returned by the  \code{\link{mining}} function.
          }
  \item{x}{the predictions (should be a numeric vector if \code{task="reg"}, matrix if \code{task="prob"} or factor if \code{task="class"} (use if \code{y} is not a list).}
  \item{metric}{a R function or if a character valid options are:
              \itemize{
                       \item \code{AUC} -- overall area under the curve (of ROC curve, classification).
                       \item \code{NAUC} -- normalized AUC (given a fixed \code{val=}FPR, classification).
                       \item \code{TPRATFPR} -- the TPR (given a fixed \code{val=}FPR, classification).
                       \item \code{ALIFT} -- area of the accumulative percent of responses captured (LIFT accumulative curve, classification).
                       \item \code{NALIFT} -- normalized ALIFT (given a fixed \code{val=}percentage of examples, classification).
                       \item \code{ALIFTATPERC} -- ALIFT value (given a fixed \code{val=}percentage of examples, classification).
                       \item \code{ACC} -- overall classification accuracy rate (classification).
                       \item \code{BRIER} -- overall Brier score (for probabilistic classification).
                       \item \code{CONF} -- confusion matrix (classification).
                       \item \code{TPR} -- true positive rate (classification).
                       \item \code{TNR} -- true negative rate (classification).
                       \item \code{KAPPA} -- kappa index (classification).
                       \item \code{SAD} -- sum absolute error (regression).
                       \item \code{MAD} or \code{MAE} -- mean absolute error (regression).
                       \item \code{MdAE} -- median absolute error (regression).
                       \item \code{GMAE} or \code{GMAD} -- geometric mean absolute error (regression).
                       \item \code{RMAD} or \code{RAE} -- relative absolute error (regression).
                       \item \code{SSE} -- sum squared error (regression).
                       \item \code{MSE} -- mean squared error (regression).
                       \item \code{RSE} -- relative squared error (regression).
                       \item \code{RMSE} -- root mean squared error (regression).
                       \item \code{ME} -- mean error (regression).
                       \item \code{COR} -- correlation (regression).
                       \item \code{R2} -- R^2 (regression).
                       \item \code{NAREC} -- normalized REC area (given a fixed \code{val=}tolerance, regression).
                       \item \code{TOLERANCE} -- the tolerance (y-axis value) of a REC curve (given a fixed \code{val=}tolerance, regression).
                       \item \code{MdAPE} -- Mean Absolute Percentage mmetric forecasting metric (regression).
                       \item \code{RMSPE} -- Root Mean Square Percentage mmetric forecasting metric (regression).
                       \item \code{RMdSPE} -- Root Median Square Percentage mmetric forecasting metric (regression).
                       \item \code{MAPE} --  Mean Absolute Percentage mmetric forecasting metric (regression).
                       \item \code{SMAPE} --  Symmetric Mean Absolute Percentage mmetric forecasting metric (regression).
                       \item \code{SMdAPE} -- Symmetric Median Absolute Percentage mmetric forecasting metric  (regression).
                       \item \code{MRAE} -- Mean Relative Absolute mmetric forecasting metric (\code{val} includes the time series in-samples, or training data, regression).
                       \item \code{MdRAE} -- Median Relative Absolute mmetric forecasting metric (\code{val} includes the time series in-samples, or training data, regression).
                       \item \code{GMRAE} -- Geometric Mean Relative Absoluate mmetric forecasting metric  (\code{val} includes the time series in-samples, or training data, regression).
                       \item \code{THEILSU2} -- Theils'U2 forecasting metric (\code{val} includes the time series in-samples, or training data, regression).
                       \item \code{MASE} -- MASE forecasting metric  (\code{val} includes the time series in-samples, or training data, regression).
                      }

}
  \item{D}{decision threshold (for \code{task="prob"}, probabilistic classification) within [0,1]. The class is TRUE if \emph{prob>D}.}
  \item{TC}{target class (for multi-class classification class) within {1,...,\emph{Nc}}, where \emph{Nc} is the number of classes.
            Notes: if \code{TC==-1} (the default value), then it is assumed:
            \itemize{
            \item if \code{metric} is "AUC", "ACC", "CONF", "KAPPA" or "BRIER" -- the global metric (for all classes) is computed.
            \item if \code{metric} is "ALIFTATPERC", "NALIFT", "LIFT", "NAUC", "TPRATFPR", "TPR" or "TNR" -- \code{TC} is set to the index of the last class.
                    }
          }
            
  \item{val}{auxiliary value, check the \code{metric} argument for details.}
}
\details{
Compute classification or regression error metrics:
\itemize{
         \item \code{mmetric} -- compute one error metric given y and x or given a mining list.
         \item \code{metrics} -- computes several classification or regression metrics. Additional arguments are:
         \itemize{
  \item \code{AUC} -- if TRUE compute the AUC metrics.
  \item \code{BRIER} -- if TRUE compute the Brier score.
  \item \code{task} --see \code{\link{fit}}.
                 }
         In \code{metrics}, \code{tauc} is the global AUC and \code{acc class} is the each class accuracy.\cr
         The \code{metrics} usage is: \code{metrics(y, x, D = 0.5, TC = -1, AUC = TRUE, BRIER = FALSE, Run=1, task = "default")}\cr
         where:
         \itemize{
         \item \code{Run} -- is the Run of mining execution if \code{y} is a mining object.
         \item if \code{AUC=TRUE} or \code{BRIER=TRUE} then AUC/BRIER statistics are also computed.
                 }
}
}
\value{
Returns the computed error metric(s).
}
\references{
\itemize{
\item To check for more details about rminer and for citation purposes:\cr
P. Cortez.\cr
Data Mining with Neural Networks and Support Vector Machines Using the R/rminer Tool.\cr
In P. Perner (Ed.), Advances in Data Mining - Applications and Theoretical Aspects 10th Industrial Conference on Data Mining (ICDM 2010), Lecture Notes in Artificial Intelligence 6171, pp. 572-583, Berlin, Germany, July, 2010. Springer. ISBN: 978-3-642-14399-1.\cr
@Springer: \url{http://www.springerlink.com/content/e7u36014r04h0334}\cr
\url{http://www3.dsi.uminho.pt/pcortez/2010-rminer.pdf}\cr

\item About the Brier and Global AUC scores:\cr
A. Silva, P. Cortez, M.F. Santos, L. Gomes and J. Neves.\cr
Rating Organ Failure via Adverse Events using Data Mining in the Intensive Care Unit.\cr
In Artificial Intelligence in Medicine, Elsevier, 43 (3): 179-193, 2008.\cr
\url{http://dx.doi.org/10.1016/j.artmed.2008.03.010}\cr

\item About the classification and regression metrics:\cr
I. Witten and E. Frank.\cr
Data Mining: Practical machine learning tools and techniques.\cr
Morgan Kaufmann, 2005.

\item About the forecasting metrics:\cr
R. Hyndman and A. Koehler\cr
Another look at measures of forecast accuracy.\cr
In International Journal of Forecasting, 22(4):679-688, 2006.\cr
}
}
\author{
Paulo Cortez \url{http://www3.dsi.uminho.pt/pcortez}
}
\note{
See also \url{http://www3.dsi.uminho.pt/pcortez/rminer.html}
}

\seealso{
\code{\link{fit}}, \code{\link{predict.fit}}, \code{\link{mining}}, \code{\link{mgraph}}, \code{\link{savemining}} and \code{\link{Importance}}.
}
\examples{
### regression
y=c(1:94,95.01,96.1,97.2,97.8,99.3,99.7);x=rnorm(100,5,0.3)+y
print(mmetric(y,x,"MAE"))
print(mmetric(y,x,"RMSE"))
print(mmetric(y,x,"TOLERANCE",val=5))
print(mmetric(y[95:100],x[95:100],"THEILSU2",val=y[1:94])) # b = in-samples
print(mmetric(y[95:100],x[95:100],"MASE",val=y[1:94])) # b = in-samples
print(metrics(y,x)) # several regression metrics
# user defined error function example:
# myerror = number of samples with absolute error above 10% of y: 
myerror=function(y,x){return (sum(abs(y-x)>(0.1*y)))}
print(mmetric(y,x,metric=myerror))

### binary classification 
y=factor(c("a","a","b","b"))
x=factor(c("a","b","b","b"))
print(mmetric(y,x,"CONF"))
print(mmetric(y,x,"ACC"))
print(metrics(y,x))
px=matrix(nrow=4,ncol=2)
px[1,]=c(0.7,0.3)
px[2,]=c(0.4,0.6)
px[3,]=c(0.7,0.3)
px[4,]=c(0.3,0.7)
print(px)
print(mmetric(y,px,"CONF"))
print(mmetric(y,px,"ACC"))
print(mmetric(y,px,"CONF",D=0.7,TC=2))
print(metrics(y,px,D=0.7,TC=2))
px2=px[,2]
print(px2)
print(mmetric(y,px,"CONF"))
print(mmetric(y,px2,"CONF",D=0.7,TC=2))
print(mmetric(y,px,"AUC"))
print(mmetric(y,px2,"AUC"))
print(mmetric(y,px2,"AUC",TC=2))

### multi-class classification 
data(iris)
M=mining(Species~.,iris,model="dt",Runs=2)
print(mmetric(M,metric="ACC",TC=2,D=0.7))
print(mmetric(M,metric="CONF",TC=2,D=0.7))
print(mmetric(M,metric="AUC",TC=3))
print(mmetric(M,metric="AUC",TC=1))
print(mmetric(M,metric="TPR",TC=1))
print(mmetric(M,metric="TPRATFPR",TC=3,val=0.05))
print(mmetric(M,metric="NAUC",TC=3,val=0.05))
print(mmetric(M,metric="ALIFT",TC=3))
print(mmetric(M,metric="ALIFTATPERC",TC=3,val=0.1))
print(mmetric(M,metric="NALIFT",TC=3,val=0.1))
print(metrics(M,BRIER=TRUE,Run=1)) # several Run 1 classification metrics
print(metrics(M,BRIER=TRUE,Run=1,TC=3)) # several Run 1 TC=3 classification metrics
}
\keyword{classif}
\keyword{regression}
