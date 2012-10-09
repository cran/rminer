\name{mmetric}
\alias{mmetric}
\alias{metrics}
\title{
Compute classification or regression error metrics.
}
\description{
Compute classification or regression error metrics.
}
\usage{mmetric(y, x = NULL, metric, D = 0.5, TC = -1, val = NULL, aggregate = "no")}

\arguments{
  \item{y}{if there are predictions (\code{!is.null(x)}), \code{y} should be a numeric vector or factor with the target desired responses (or output values).\cr
           Else, \code{y} should be a list returned by the  \code{\link{mining}} function.
          }
  \item{x}{the predictions (should be a numeric vector if \code{task="reg"}, matrix if \code{task="prob"} or factor if \code{task="class"} (used if \code{y} is not a list).}
  \item{metric}{a R function or a character with valid options (">" means "better" if higher value; "<" means "better" if lower value):
              \itemize{
                       \item \code{ALL} -- returns all classification or regression metrics (context dependent, multi-metric). 
                       \item if vector -- returns all metrics included in the vector, vector elements can be any of the options below (multi-metric).

                       \item \code{CONF} -- confusion matrix (classification, matrix). 

                       \item \code{ACC} -- classification accuracy rate (classification, ">", [0-\%100]). 
                       \item \code{CE} -- classification error or misclassification error rate (classification, "<", [0-\%100]). 
                       \item \code{MAEO} -- mean absolute error for ordinal classification (classification, "<", [0-Inf[). 
                       \item \code{MSEO} -- mean squared error for ordinal classification (classification, "<", [0-Inf[). 
                       \item \code{KENDALL} -- Kendalls's coefficient  for ordinal classification or (mean if) ranking (classification, ">", [-1;-1]). Note: if ranking, \code{y} is a matrix and mean metric is computed.
                       \item \code{SPEARMAN} -- Mean Spearman's rho coefficient  for ranking (classification, ">", [-1;-1]). Note: if ranking, \code{y} is a matrix and mean metric is computed.
                       \item \code{MSE} -- classification error or misclassification error rate (classification, "<", [0-\%100]). 
                       \item \code{BER} -- balanced error rate (classification, "<", [0-\%100]). 
                       \item \code{KAPPA} -- kappa index (classification, "<", [0-\%100]). 
                       \item \code{CRAMERV} -- Cramer's V  (classification, ">", [0,1.0]). 

                       \item \code{ACCLASS} -- classification accuracy rate per class (classification, ">", [0-\%100]). 
                       \item \code{TPR} -- true positive rate, sensitivity or recall (classification, ">", [0-\%100]). 
                       \item \code{TNR} -- true negative rate or specificity (classification, ">", [0-\%100]). 
                       \item \code{PRECISION} -- precision (classification, ">", [0-\%100]). 
                       \item \code{F1} -- F1 score (classification, ">", [0-\%100]). 
                       \item \code{MCC} -- Matthews correlation coefficient (classification, ">", [-1,1]). 
                       
                       \item \code{BRIER} -- overall Brier score (classification "prob", "<", [0,1.0]).
                       \item \code{BRIERCLASS} -- Brier score per class (classification "prob", "<", [0,1.0]).
                       \item \code{ROC} --  Receiver Operating Characteristic curve (classification "prob", list with several components).
                       \item \code{AUC} -- overall area under the curve (of ROC curve, classification "prob", ">", domain values: [0,1.0]). 
                       \item \code{AUCCLASS} -- area under the curve per class (of ROC curve, classification "prob", ">", domain values: [0,1.0]). 
                       \item \code{NAUC} -- normalized AUC (given a fixed \code{val=}FPR, classification "prob", ">", [0,1.0]). 
                       \item \code{TPRATFPR} -- the TPR (given a fixed \code{val=}FPR, classification "prob", ">", [0,1.0]).
                       \item \code{LIFT} -- accumulative percent of responses captured (LIFT accumulative curve, classification "prob", list with several components).
                       \item \code{ALIFT} -- area of the accumulative percent of responses captured (LIFT accumulative curve, classification "prob", ">", [0,1.0]).
                       \item \code{NALIFT} -- normalized ALIFT (given a fixed \code{val=}percentage of examples, classification "prob", ">", [0,1.0]).
                       \item \code{ALIFTATPERC} -- ALIFT value (given a fixed \code{val=}percentage of examples, classification "prob", ">", [0,1.0]).

                       \item \code{SAE} -- sum absolute error/deviation (regression, "<", [0,Inf[).
                       \item \code{MAE} -- mean absolute error (regression, "<", [0,Inf[).
                       \item \code{MdAE} -- median absolute error (regression, "<", [0,Inf[).
                       \item \code{GMAE} -- geometric mean absolute error (regression, "<", [0,Inf[).
                       \item \code{MaxAE} -- maximum absolute error (regression, "<", [0,Inf[).
                       \item \code{RAE} -- relative absolute error (regression, "<", [0\%,Inf[).
                       \item \code{SSE} -- sum squared error (regression, "<", [0,Inf[).
                       \item \code{MSE} -- mean squared error (regression, "<", [0,Inf[).
                       \item \code{MdSE} -- median squared error (regression, "<", [0,Inf[).
                       \item \code{RMSE} -- root mean squared error (regression, "<", [0,Inf[).
                       \item \code{GMSE} -- geometric mean squared error (regression, "<", [0,Inf[).
                       \item \code{HRMSE} -- Heteroscedasticity consistent root mean squared error (regression, "<", [0,Inf[).
                       \item \code{RSE} -- relative squared error (regression, "<", [0\%,Inf[).
                       \item \code{RRSE} -- root relative squared error (regression, "<", [0\%,Inf[).
                       \item \code{ME} -- mean error (regression, "<", [0,Inf[).
                       \item \code{SMinkowski3} -- sum of Minkowski loss function (q=3, heavier penalty for large errors when compared with SSE, regression, "<", [0\%,Inf[).
                       \item \code{MMinkowski3} -- mean of Minkowski loss function (q=3, heavier penalty for large errors when compared with SSE, regression, "<", [0\%,Inf[).
                       \item \code{MdMinkowski3} -- median of Minkowski loss function (q=3, heavier penalty for large errors when compared with SSE, regression, "<", [0\%,Inf[).
                       \item \code{COR} -- correlation (regression, ">", [-1,1]).
                       \item \code{q2} -- =1-correlation^2 test error metric, as used by M.J. Embrechts  (regression, "<", [0,1.0]).
                       \item \code{R2} -- "press" R^2 (regression, ">", ]-Inf,1]).
                       \item \code{Q2} -- "press"/SD test error metric, as used by M.J. Embrechts (regression, "<", [0,Inf[).
                       \item \code{REC} --  Regression Error Characteristic curve (regression, list with several components).
                       \item \code{NAREC} -- normalized REC area (given a fixed \code{val=}tolerance, regression, ">", [0,1.0]).
                       \item \code{TOLERANCE} -- the tolerance (y-axis value) of a REC curve (given a fixed \code{val=}tolerance, regression, ">", [0,1.0]).
                       \item \code{MAPE} --  Mean Absolute Percentage mmetric forecasting metric (regression, "<", [0\%,Inf[).
                       \item \code{MdAPE} -- Mean Absolute Percentage mmetric forecasting metric (regression, "<"), [0\%,Inf[).
                       \item \code{RMSPE} -- Root Mean Square Percentage mmetric forecasting metric (regression, "<", [0\%,Inf[).
                       \item \code{RMdSPE} -- Root Median Square Percentage mmetric forecasting metric (regression, "<", [0\%,Inf[).
                       \item \code{SMAPE} --  Symmetric Mean Absolute Percentage mmetric forecasting metric (regression, "<", [0\%,200\%]).
                       \item \code{SMdAPE} -- Symmetric Median Absolute Percentage mmetric forecasting metric  (regression, "<", [0\%,200\%]).
                       \item \code{MRAE} -- Mean Relative Absolute mmetric forecasting metric (\code{val} should contain the last in-sample/training data value (for random walk) or full benchmark time series, regression, "<", [0,Inf[).
                       \item \code{MdRAE} -- Median Relative Absolute mmetric forecasting metric (\code{val} should contain the last in-sample/training data value (for random walk) or full benchmark time series, regression, "<", [0,Inf[).
                       \item \code{GMRAE} -- Geometric Mean Relative Absoluate mmetric forecasting metric  (\code{val} should contain the last in-sample/training data value (for random walk) or full benchmark time series, regression, "<", [0,Inf[).
                       \item \code{THEILSU2} -- Theils'U2 forecasting metric (\code{val} should contain the last in-sample/training data value (for random walk) or full benchmark time series, regression, "<", [0,Inf[).
                       \item \code{MASE} -- MASE forecasting metric  (\code{val} should contain the time series in-samples or training data, regression, "<", [0,Inf[).
                      }

}
  \item{D}{decision threshold (for \code{task="prob"}, probabilistic classification) within [0,1]. The class is TRUE if \emph{prob>D}.}
  \item{TC}{target class index or vector of indexes (for multi-class classification class) within {1,...,\emph{Nc}}, where \emph{Nc} is the number of classes:<cr>
            \itemize{
            \item if \code{TC==-1} (the default value), then it is assumed:
            \itemize{
            \item if \code{metric} is "CONF" -- \code{D} is ignored and highest probability class is assumed (if \code{TC>0}, the metric is computed for positive \code{TC} class and \code{D} is used).
            \item if \code{metric} is "ACC", "CE", "BER", "KAPPA", "CRAMERV", "BRIER", or "AUC" -- the global metric (for all classes) is computed (if \code{TC>0}, the metric is computed for positive \code{TC} class).
            \item if \code{metric} is "ACCLASS", "TPR", "TNR", "Precision", "F1", "MCC", "ROC", "BRIERCLASS", "AUCCLASS" -- it returns one result per class (if \code{TC>0}, it returns negative (e.g. "TPR1") and positive (TC, e.g. "TPR2") result).
            \item if \code{metric} is "NAUC", "TPRATFPR", "LIFT", "ALIFT", "NALIFT" or "ALIFTATPERC" -- TC is set to the index of the last class.
                    }
                    }
          }
            
  \item{val}{auxiliary value:
            \itemize{
            \item when two or more metrics need different \code{val} values, then \code{val} should be a vector list, see example.
            \item if numeric or vector -- check the \code{metric} argument for specific details of each metric val meaning.
            }
            }

  \item{aggregate}{character with type of aggregation performed when y is a \code{\link{mining}} list. Valid options are:
            \itemize{
            \item no -- returns all metrics for all \code{\link{mining}} runs. If \code{metric} includes "CONF", "ROC", "LIFT" or "REC", it returns a vector list, else if \code{metric} includes a single metric, it returns a vector; else it retuns a data.frame (runs x metrics).
            \item sum -- sums all run results. 
            \item mean -- averages all run results.
            \item note: both "sum" and "mean" only work if only \code{metric=="CONF"} is used or if \code{metric} does not contain "ROC", "LIFT" or "REC".
            }
                  }
}
\details{
Compute classification or regression error metrics:
\itemize{
         \item \code{mmetric} -- compute one or more classification/regression metrics given y and x OR a mining list.
         \item \code{metrics} -- deprecated function, same as \code{mmetric(x,y,metric="ALL")}, included here just for compatability purposes but will be removed from the package.
        }
 }
\value{
Returns the computed error metric(s): 
\itemize{
\item one value if only one \code{metric} is requested (and \code{y} is not a mining list);
\item named vector if 2 or more elements are requested in \code{metric} (and \code{y} is not a mining list); 
\item list if there is a "CONF", "ROC", "LIFT" or "REC" request on \code{metric} (other metrics are stored in field \code{$res}, and \code{y} is not a mining list).
\item if \code{y} is a mining list then there can be several runs, thus:
 \itemize{
  \item a vector list of size \code{y$runs} is returned if \code{metric} includes "CONF", "ROC", "LIFT" or "REC" and \code{aggregate="no"};
  \item a data.frame is returned if \code{aggregate="no"} and  \code{metric} does not include "CONF", "ROC", "LIFT" or "REC";
  \item a table is returned if \code{aggregate="sum" or "mean"} and \code{metric="CONF"};
  \item a vector or numeric value is returned if  \code{aggregate="sum" or "mean"}  and \code{metric} is not "CONF".
         }
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

\item About the ordinal classification metrics:\cr
J.S. Cardoso and R. Sousa.\cr
Measuring the Performance of Ordinal Classification.\cr
In International Journal of Pattern Recognition and Artificial Intelligence, 25(8):1173-1195, 2011.\cr
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
### regression examples:
y=c(95.01,96.1,97.2,97.8,99.3,99.7);x=95:100
print(mmetric(y,x,"ALL"))
print(mmetric(y,x,"MAE"))
m=mmetric(y,x,c("MAE","RMSE","RAE","RSE"))
print(m)
cat(names(m)[3],"=",round(m[3],digit=2),"\n",sep="")
print(mmetric(y,x,c("COR","R2","Q2")))
print(mmetric(y,x,c("TOLERANCE","NAREC"),val=0.1))
print(mmetric(y,x,"THEILSU2",val=94)) # val = 1-ahead random walk
print(mmetric(y,x,"THEILSU2",val=94:99)) # val = 1-ahead random walk
print(mmetric(y,x,"MASE",val=1:94)) # val = in-samples
val=vector("list",length=4)
val[[2]]=0.5;val[[3]]=94;val[[4]]=1:94
print(mmetric(y,x,c("MAE","NAREC","THEILSU2","MASE"),val=val))
# user defined error function example:
# myerror = number of samples with absolute error above 10% of y: 
myerror=function(y,x){return (sum(abs(y-x)>(0.1*y)))}
print(mmetric(y,x,metric=myerror))
# example that returns a list since "REC" is included:
print(mmetric(y,x,c("MAE","REC","TOLERANCE"),val=1))

### pure binary classification 
y=factor(c("a","a","a","a","b","b","b","b"))
x=factor(c("a","a","b","a","b","a","b","a"))
print(mmetric(y,x,"CONF")$conf)
print(mmetric(y,x,"ALL"))
print(mmetric(y,x,metric=c("ACC","TPR","ACCLASS")))

### probabilities binary classification 
y=factor(c("a","a","a","a","b","b","b","b"))
px=matrix(nrow=8,ncol=2)
px[,1]=c(1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3)
px[,2]=1-px[,1]
print(px)
print(mmetric(y,px,"CONF")$conf)
print(mmetric(y,px,"CONF",D=0.5,TC=2)$conf)
print(mmetric(y,px,"CONF",D=0.3,TC=2)$conf)
print(mmetric(y,px,metric="ALL",D=0.3,TC=2))
print(mmetric(y,px,metric=c("ACC","AUC","AUCCLASS","BRIER","BRIERCLASS","CE"),D=0.3,TC=2))

### pure multi-class classification 
y=factor(c("a","a","b","b","c","c"))
x=factor(c("a","a","b","c","b","c"))
print(mmetric(y,x,metric="CONF")$conf)
print(mmetric(y,x,metric="CONF",TC=-1)$conf)
print(mmetric(y,x,metric="CONF",TC=1)$conf)
print(mmetric(y,x,metric="ALL"))
print(mmetric(y,x,metric=c("ACC","ACCLASS","KAPPA")))
print(mmetric(y,x,metric=c("ACC","ACCLASS","KAPPA"),TC=1))

### probabilities multi-class 
y=factor(c("a","a","b","b","c","c"))
px=matrix(nrow=6,ncol=3)
px[,1]=c(1.0,0.7,0.5,0.3,0.1,0.7)
px[,2]=c(0.0,0.2,0.4,0.7,0.3,0.2)
px[,3]=1-px[,1]-px[,2]
print(px)
print(mmetric(y,px,metric=c("AUC","AUCCLASS","NAUC"),TC=-1,val=0.1))
print(mmetric(y,px,metric=c("AUC","NAUC"),TC=3,val=0.1))
print(mmetric(y,px,metric=c("ACC","ACCLASS"),TC=-1))
print(mmetric(y,px,metric=c("CONF"),TC=3,D=0.5)$conf)
print(mmetric(y,px,metric=c("ACCLASS"),TC=3,D=0.5))
print(mmetric(y,px,metric=c("CONF"),TC=3,D=0.7)$conf)
print(mmetric(y,px,metric=c("ACCLASS"),TC=3,D=0.7))

### ordinal multi-class (example in Ricardo Sousa PhD thesis 2012)
y=ordered(c(rep("a",4),rep("b",6),rep("d",3)),levels=c("a","b","c","d"))
x=ordered(c(rep("c",(4+6)),rep("d",3)),levels=c("a","b","c","d"))
print(mmetric(y,x,metric="CONF")$conf)
print(mmetric(y,x,metric=c("CE","MAEO","MSEO","KENDALL")))
# note: only y needs to be ordered
x=factor(c(rep("b",4),rep("a",6),rep("d",3)),levels=c("a","b","c","d"))
print(mmetric(y,x,metric="CONF")$conf)
print(mmetric(y,x,metric=c("CE","MAEO","MSEO","KENDALL")))

### ranking (multi-class) 
y=matrix(nrow=1,ncol=12);x=y
# http://www.youtube.com/watch?v=D56dvoVrBBE
y[1,]=1:12
x[1,]=c(2,1,4,3,6,5,8,7,10,9,12,11)
print(mmetric(y,x,metric="KENDALL"))
print(mmetric(y,x,metric="ALL"))

y=matrix(nrow=2,ncol=7);x=y
y[1,]=c(2,6,5,4,3,7,1)
y[2,]=7:1
x[1,]=1:7
x[2,]=1:7
print(mmetric(y,x,metric="ALL"))

### IMPORTANT NOTE ###
# The execution of the examples below requires more than a few seconds.
# Thus, to fulfill CRAN policies, I have added 2 lines: IF(FALSE){ ... }
# Such lines are marked with: # DELETE-OR-IGNORE ...
# If you want to test the code, please comment or ignore these 2 lines
### END OF NOTE    ###
if(FALSE){ # DELETE-OR-IGNORE-THIS AND LAST EXAMPLE LINE #
### mining, several runs, prob multi-class
data(iris)
M=mining(Species~.,iris,model="dt",Runs=2)
R=mmetric(M,metric="CONF",aggregate="no")
print(R[[1]]$conf)
print(R[[2]]$conf)
print(mmetric(M,metric="CONF",aggregate="mean"))
print(mmetric(M,metric="CONF",aggregate="sum"))
print(mmetric(M,metric=c("ACC","ACCLASS"),aggregate="no"))
print(mmetric(M,metric=c("ACC","ACCLASS"),aggregate="mean"))
print(mmetric(M,metric="ALL",aggregate="no"))
print(mmetric(M,metric="ALL",aggregate="mean"))

### mining, several runs, regression
data(sin1reg)
S=sample(1:nrow(sin1reg),40)
M=mining(y~.,data=sin1reg[S,],model="svm",search=2^3,Runs=10)
R=mmetric(M,metric="MAE")
print(mmetric(M,metric="MAE",aggregate="mean"))
miR=meanint(R) # mean and t-student confidence intervals
cat("MAE=",round(miR$mean,digits=2),"+-",round(miR$int,digits=2),"\n")
print(mmetric(M,metric=c("MAE","RMSE")))
print(mmetric(M,metric=c("MAE","RMSE"),aggregate="mean"))
R=mmetric(M,metric="REC",aggregate="no")
print(R[[1]]$rec)
print(mmetric(M,metric=c("TOLERANCE","NAREC"),val=0.2))
print(mmetric(M,metric=c("TOLERANCE","NAREC"),val=0.2,aggregate="mean"))
} # DELETE-OR-IGNORE-THIS LINE #

}
\keyword{classif}
\keyword{regression}
