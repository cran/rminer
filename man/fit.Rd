\name{fit}
\alias{fit} %- Also NEED an '\alias' for EACH other topic documented here.
\alias{model-class} %- Also NEED an '\alias' for EACH other topic documented here.
\alias{svmgrid} %- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a supervised data mining model (classification or regression) model
}
\description{
Fit a supervised data mining model (classification or regression) model. Kind of a wrapper function that allows to fit distinct data mining methods under
the same coherent function structure. Also, it tunes the hyperparameters of some models (e.g. \code{knn}, \code{mlp}, \code{mlpe} and \code{svm}) and performs
some feature selection methods.
}
\usage{
fit(x, data = NULL, model = "default", task = "default", 
    search = "heuristic", mpar = NULL, feature = "none", 
    scale = "default", transform = "none", 
    created = NULL, ...)
}

\arguments{
  \item{x}{a symbolic description (formula) of the model to be fit. If \code{x} contains the data, then \code{data=NULL} (similar to x in \code{\link{ksvm}}, kernlab package).}
  \item{data}{an optional data frame (columns denote attributes, rows show examples) containing the training data, when using a formula.}
  \item{model}{a character with the model type name (data mining method). Valid options are: 
              \itemize{
              \item \code{naive} most common class (classification) or mean output value (regression)
              \item \code{lr} (or \code{logistic}) -- logistic regression (classification, uses \code{\link{multinom}})
              \item \code{naivebayes} -- naive bayes (classification, very similar to \code{naiveBayes})
              \item \code{lda} -- linear discriminant analysis (classification, uses \code{\link{lda}}, requires \code{MASS} package, i.e. \code{library(MASS)})
              \item \code{qda} -- quadratic discriminant analysis (classification, uses \code{\link{qda}}, requires \code{MASS} package)
              \item \code{dt} -- decision tree (classification and regression, uses \code{\link{rpart}})
              \item \code{mr} -- multiple regression (regression, equivalent to \code{\link{lm}} but uses \code{\link{nnet}} with zero hidden nodes and linear output function)
              \item \code{bruto} -- additive spline model (regression, uses \code{bruto}, requires \code{mda} package, i.e. \code{library(mda)})
              \item \code{mars} -- multivariate adaptive regression splines (regression, uses \code{mars}, requires \code{mda} package, i.e. \code{library(mda)})
              \item \code{knn} -- k-nearest neighbor (classification and regression, uses \code{\link{kknn}})
              \item \code{mlp} -- multilayer perceptron with one hidden layer (classification and regression, uses \code{\link{nnet}})
              \item \code{mlpe} -- multilayer perceptron ensemble (classification and regression, uses \code{\link{nnet}})
              \item \code{svm} -- support vector machine (classification and regression, uses \code{\link{ksvm}})
              \item \code{randomforest} -- random forest algorithm (classification and regression, uses and requires \code{randomForest} package, i.e. 
              \code{library(randomForest)})
                      }
              }
  \item{task}{data mining task. Valid options are:
              \itemize{
                       \item \code{prob} (or \code{p}) -- classification with output probabilities (i.e. the sum of all outputs equals 1).
                       \item \code{class} (or \code{c}) -- classification with discrete outputs (\code{\link{factor}}) 
                       \item \code{reg} (or \code{r}) -- regression (numeric output) 
                       \item \code{default} tries to guess the best task (\code{prob} or \code{reg}) given the \code{model} and output variable type (if factor then \code{prob} else \code{reg})
                      }
             }
  \item{search}{used to tune the hyperparameter(s) of the model (only for: \code{knn} -- number of neighbors (k); 
							\code{mlp} or \code{mlpe} -- number of hidden nodes (\emph{H}) or \emph{decay};
							\code{svm} -- gaussian kernel parameter (\emph{sigma});
							\code{randomforest} -- \code{mtry} parameter). Valid options are:
              \itemize{
                       \item \code{heuristic} -- simple heuristic, one search parameter (e.g. \emph{k}=3 for \code{knn}, \emph{H}=inputs/2 for \code{mlp}, \emph{sigma}=2^-7 for \code{svm})
                       \item \code{heuristic5} -- heuristic with a 5 range grid-search (e.g. \code{seq(1,9,2)} for \code{knn}, \code{seq(0,8,2)} for \code{mlp}, 
     										 	     \code{2^seq(-15,3,4)} for \code{svm}, \code{1:5} for \code{randomforest})
                       \item \code{heuristic10} -- heuristic with a 10 range grid-search (e.g. \code{seq(1,10,1)} for \code{knn}, \code{seq(0,9,1)} for \code{mlp}, 
     										 	     \code{2^seq(-15,3,2)} for \code{svm}, \code{1:10} for \code{randomforest})
                       \item \code{UD}, \code{UD1} or \code{UD2} -- uniform design 2-Level with 13 (\code{UD} or \code{UD2}) or 21 (\code{UD1}) searches (for \code{svm}).
                       \item a-\code{\link{vector}} -- numeric vector with all hyperparameter values that will be searched within an internal grid-search (the number of searches 
                       is \code{length(search)} when \code{convex=0})

                       \item a-\code{\link{matrix}} -- numeric matrix with all hyperparameter values (one different parameter per column) that will be searched within an internal grid-search (the number of searches is \code{nrow(search)} when \code{convex=0})
                       \item a-\code{\link{list}} -- a list with: 
                                   \itemize{ 
                                    \item \code{$search} -- a-\code{\link{vector}} with all hyperparameter values or uniform design ranges
                                    \item \code{$convex} -- number that defines how many searches are performed after a local minimum/maximum is found (if >0, the search can be stopped without testing all grid-search values) 
                                    \item \code{$smethod} -- search method. Valid options are: 
                                          \itemize{
                                                   \item \code{normal} -- default value for normal grid search.
                                                   \item \code{2L} - nested 2-Level grid search. First level range is set by \code{$search} and then the 2nd level performs a fine tuning, with \code{length($search)} searches around best value in first level.
                                                   \item \code{UD}, \code{UD1} or \code{UD2} -- uniform design 2-Level with 13 (\code{UD} or \code{UD2}) or 21 (\code{UD1}) searches (for \code{svm}). Under this option,
\code{$search} should contain the first level ranges, such as \code{c(-15,3,-5,15)} for classification (\emph{gamma} min and max, \emph{C} min and max, 
after which a \code{2^} transform is applied) or \code{c(-8,0,-1,6,-8,-1)} for regression (last two values are \emph{epsilon} min and max, after which a \code{2^} transform is applied).
                                                  }
                                           }
                      }
              }

  \item{mpar}{vector with extra model parameters (used for modeling, search and feature selection) with:	
                                    \itemize{ 
                                               \item c(\emph{vmethod},\emph{vpar},\code{metric}) -- if \code{model}= \code{knn} or \code{randomforest};
                                               \item c(\emph{C},\emph{epsilon},\emph{vmethod},\emph{vpar},\code{metric}) -- if \code{svm};
					       \item c(\emph{Nr},\emph{Me},\emph{vmethod},\emph{vpar},\code{metric}) -- if \code{mlp} or \code{mlpe} and \code{search} for \emph{H} is used; or
					       \item c(\emph{Nr},\emph{Me},\emph{vmethod},\emph{vpar},\code{metric},\emph{H}) -- if \code{mlp} or \code{mlpe} and
                                               \code{search} for \emph{decay} is used.
                                            }
              \emph{C} and \emph{epsilon} are default values for \code{svm} (if any of these is \code{=NA} then heuristics are used to set the value).\cr
              \emph{Nr} is the number of \code{mlp} runs or \code{mlpe} individual models, while \emph{Me} is the maximum number of epochs (if any of these is \code{=NA} then heuristics are used to set the value).\cr
              In the \code{fit} function, \emph{vmethod} can only be set to: "all", "holdout", "holdoutorder", "kfold" or "kfoldo". For help on \emph{vmethod} and \emph{vpar} see \code{\link{mining}}.\cr
              \code{metric} is the internal error function (e.g. used by search to select the best model), valid options are explained in \code{\link{mmetric}}. 
              When \code{mpar=NULL} then default values are used. If there are \code{NA} values (e.g. \code{mpar=c(NA,NA)}) then default values are used.
             }

  \item{feature}{feature selection and sensitivity analysis control. Valid \code{fit} function options are:
                  \itemize{ 
                           \item \code{none} -- no feature selection;
                           \item a-\code{\link{vector}} -- vector with c(\emph{fmethod},\emph{deletions},\emph{Runs},\emph{vmethod},\emph{vpar},\emph{defaultsearch})
                           \item a-\code{\link{vector}} -- vector with c(\emph{fmethod},\emph{deletions},\emph{Runs},\emph{vmethod},\emph{vpar})
                          }
                 \emph{fmethod} sets the type. Valid options are: 
                  \itemize{ 
                           \item \code{sbs} -- standard backward selection;
                           \item \code{sabs} -- sensitivity analysis backward selection (faster);
                           \item \code{sabsv} -- equal to \code{sabs} but uses variance for sensitivity importance measure;
                           \item \code{sabsr} -- equal to \code{sabs} but uses range for sensitivity importance measure;
                           \item \code{sabsg} -- equal to \code{sabs} (uses gradient for sensitivity importance measure);
                          }
                 \emph{deletions} is the maximum number of feature deletions (if -1 not used).\cr
                 \emph{Runs} is the number of runs for each feature set evaluation (e.g. 1).\cr
                 For help on \emph{vmethod} and \emph{vpar} see \code{\link{mining}}.\cr
                 \emph{defaultsearch} is one hyperparameter used during the feature selection search, after selecting the best feature set then \code{search} is used (faster). 
                 If not defined, then \code{search} is used during feature selection (may be slow).\cr
                 When feature is a vector then default values are used to fill missing values or \code{NA} values.
                }

  \item{scale}{if data needs to be scaled (i.e. for \code{mlp} or \code{mlpe}). Valid options are:
                  \itemize{ 
                           \item \code{default} -- uses scaling when needed (i.e. for \code{mlp} or \code{mlpe})
                           \item \code{none} -- no scaling; 
                           \item \code{inputs} -- standardizes (0 mean, 1 st. deviation) input attributes; 
                           \item \code{all} -- standardizes (0 mean, 1 st. deviation) input and output attributes;
                          }
               If needed, the \code{predict} function of rminer performs the inverse scaling.
              }

  \item{transform}{if the output data needs to be transformed (e.g. \code{log} transform). Valid options are:
                  \itemize{ 
                           \item \code{none} -- no transform; 
                           \item \code{log} -- y=(log(y+1)) (the inverse function is applied in the \code{predict} function);
                           \item \code{positive} -- all predictions are positive (negative values are turned into zero);
                           \item \code{logpositive} -- both \code{log} and \code{logpositive}; 
                          }
                  }

  \item{created}{time stamp for the model. By default, the system time is used. Else, you can specify another time.}

  \item{\dots}{ additional and specific parameters send to each fit function model (e.g. \code{dt}, \code{randomforest}). For example, the \code{rpart} function is used for \code{dt}, 
                thus you can add: \code{control=rpart.control(cp=.05)}.
              }
}

\details{
Fits a classification or regression model given a data.frame (see [Cortez, 2010] for more details):
\itemize{
\item Neural Network: \code{mlp} trains \emph{Nr} multilayer perceptrons (with \emph{Me} epochs, \emph{H} hidden nodes 
and \emph{decay} value according to the \code{\link{nnet}} function) and selects the best network according to minimum penalized error (\code{$value}). \code{mlpe} uses an ensemble
of \emph{Nr} networks and the final prediction is given by the average of all outputs. To tune \code{mlp} or \code{mlpe} you can use the \code{search} parameter, which performs a grid
search for \emph{H} or \emph{decay}.

\item Support Vector Machine: \code{svm} adopts the gaussian kernel. For classification tasks, you can use \code{search} to tune \emph{sigma} (gaussian kernel parameter) 
and \emph{C} (complexity parameter). For regression, the epsilon insensitive function is adopted and there is an additional hyperparameter \emph{epsilon}.

\item Other methods: Random Forest -- if needed, you can tune the \code{mtry} parameter using \code{search}; k-nearest neighbor -- use \code{search} to tune \emph{k}.
}
}

\value{
Returns a model object. You can check all model elements with \code{str(M)}, where \code{M} is a model object. The slots are:
                  \itemize{ 
                           \item \code{@formula} -- the \code{x};
                           \item \code{@model} -- the \code{model};
                           \item \code{@task} -- the \code{task};
                           \item \code{@mpar} -- data.frame with the best model parameters (interpretation depends on \code{model});
                           \item \code{@attributes} -- the attributes used by the model;
                           \item \code{@scale} -- the \code{scale};
                           \item \code{@transform} -- the \code{transform};
                           \item \code{@created} -- the date when the model was created;
                           \item \code{@time} -- computation effort to fit the model;
                           \item \code{@object} -- the R object model (e.g. \code{rpart}, \code{nnet}, ...);
                           \item \code{@outindex} -- the output index (of @attributes);
                           \item \code{@levels} -- if \code{task=="prob"||task=="class"} stores the output levels; 
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

\item For the sabs feature selection:\cr
P. Cortez, A. Cerdeira, F. Almeida, T. Matos and J. Reis.\cr
Modeling wine preferences by data mining from physicochemical properties.\cr
In Decision Support Systems, Elsevier, 47(4):547-553, 2009.\cr
\url{http://dx.doi.org/10.1016/j.dss.2009.05.016}\cr

\item For the uniform design details:\cr
C.M. Huang, Y.J. Lee, D.K.J. Lin and S.Y. Huang.\cr
Model selection for support vector machines via uniform design,\cr
In Computational Statistics & Data Analysis, 52(1):335-346, 2007.\cr
}
}

\author{
Paulo Cortez \url{http://www3.dsi.uminho.pt/pcortez}
}
\note{
See also \url{http://www3.dsi.uminho.pt/pcortez/rminer.html}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mining}}, \code{\link{predict.fit}}, \code{\link{mgraph}}, \code{\link{mmetric}}, \code{\link{savemining}}, \code{\link{CasesSeries}}, \code{\link{lforecast}},
\code{\link{holdout}} and \code{\link{Importance}}. Check all rminer functions using: \code{help(package=rminer)}.
}

\examples{
### IMPORTANT NOTE ###
# The execution of the examples below require more few seconds.
# Thus, to fulfill CRAN policies, I have added 2 lines: IF(FALSE){ ... }
# Such lines are marked with: # DELETE-OR-IGNORE ...
# If you want to test the code, please comment or ignore these 2 lines
### END OF NOTE    ###

if(FALSE){ # DELETE-OR-IGNORE-THIS AND LAST EXAMPLE LINE #

### simple regression (with a formula) example.
x1=rnorm(200,100,20); x2=rnorm(200,100,20)
y=0.7*sin(x1/(25*pi))+0.3*sin(x2/(25*pi))
M=fit(y~x1+x2,model="mlpe",search=2)
new1=rnorm(100,100,20); new2=rnorm(100,100,20)
ynew=0.7*sin(new1/(25*pi))+0.3*sin(new2/(25*pi))
P=predict(M,data.frame(x1=new1,x2=new2,y=rep(NA,100)))
print(mmetric(ynew,P,"MAE"))

### simple classification example.
data(iris)
M=fit(Species~.,iris,model="dt")
P=predict(M,iris)
print(mmetric(iris$Species,P,"CONF"))
print(mmetric(iris$Species,P,"ACC"))
print(mmetric(iris$Species,P,"AUC"))
print(metrics(iris$Species,P))
mgraph(iris$Species,P,graph="ROC",TC=2,main="versicolor ROC",
baseline=TRUE,leg="Versicolor",Grid=10)

### classification example with discrete classes, probabilities and holdout
H=holdout(iris$Species,ratio=2/3)
M=fit(Species~.,iris[H$tr,],model="svm",task="class")
M2=fit(Species~.,iris[H$tr,],model="svm",task="prob")
P=predict(M,iris[H$ts,])
P2=predict(M2,iris[H$ts,])
print(mmetric(iris$Species[H$ts],P,"CONF"))
print(mmetric(iris$Species[H$ts],P2,"CONF"))
print(mmetric(iris$Species[H$ts],P,"CONF",TC=1))
print(mmetric(iris$Species[H$ts],P2,"CONF",TC=1))
print(mmetric(iris$Species[H$ts],P2,"AUC"))

### classification example with hyperparameter selection
# SVM 
M=fit(Species~.,iris,model="svm",search=2^-3,mpar=c(3)) # C=3, gamma=2^-3
print(M@mpar) # gamma, C, epsilon (not used here)
M=fit(Species~.,iris,model="svm",search="heuristic10") # 10 grid search for gamma
print(M@mpar) # gamma, C, epsilon (not used here)
M=fit(Species~.,iris,model="svm",search="heuristic10") # 10 grid search for gamma
print(M@mpar) # gamma, C, epsilon (not used here)
M=fit(Species~.,iris,model="svm",search=2^seq(-15,3,2),
      mpar=c(NA,NA,"holdout",2/3,"AUC")) # same 0 grid search for gamma
print(M@mpar) # gamma, C, epsilon (not used here)
search=svmgrid(task="prob") # grid search as suggested by the libsvm authors
M=fit(Species~.,iris,model="svm",search=search) # 
print(M@mpar) # gamma, C, epsilon (not used here)
M=fit(Species~.,iris,model="svm",search="UD") # 2 level 13 point uniform-design
print(M@mpar) # gamma, C, epsilon (not used here)
# MLPE
M=fit(Species~.,iris,model="mlpe",search="heuristic5") # 5 grid search for H
print(M@mpar)
M=fit(Species~.,iris,model="mlpe",search="heuristic5",
      mpar=c(3,100,"kfold",3,"AUC",2)) # 5 grid search for decay, inner 3-fold
print(M@mpar)
# faster grid search 
M=fit(Species~.,iris,model="mlpe",search=list(smethod="normal",convex=1,search=0:9)) 
print(M@mpar)
# 2 level grid with total of 5 searches
M=fit(Species~.,iris,model="mlpe",search=list(smethod="2L",search=c(4,8,12))) 
print(M@mpar)
# 2 level grid for decay
search=list(smethod="2L",search=c(0,0.1,0.2));mpar=c(3,100,"holdout",3,"AUC",2) 
M=fit(Species~.,iris,model="mlpe",search=search,mpar=mpar)
print(M@mpar)

### regression example
data(sin1reg)
M=fit(y~.,data=sin1reg,model="svm",search="heuristic")
P=predict(M,sin1reg)
print(mmetric(sin1reg$y,P,"MAE"))
mgraph(sin1reg$y,P,graph="REC",Grid=10)
# uniform design
M=fit(y~.,data=sin1reg,model="svm",search="UD")
print(M@mpar)
# sensitivity analysis feature selection
M=fit(y~.,data=sin1reg,model="svm",search="heuristic5",feature="sabs") 
print(M@mpar)
print(M@attributes) # selected attributes (1 and 2 are the relevant inputs)
P=predict(M,sin1reg); print(mmetric(sin1reg$y,P,"MAE"))
# sensitivity analysis feature selection
M=fit(y~.,data=sin1reg,model="mlp",search=2,feature=c("sabs",-1,1,"kfold",3)) 
print(M@mpar)
print(M@attributes)

M=fit(y~.,data=sin1reg,model="svm",search="heuristic")
P=predict(M,data.frame(x1=-1000,x2=0,x3=0,y=NA)) # P should be negative...
print(P)
M=fit(y~.,data=sin1reg,model="svm",search="heuristic",transform="positive")
P=predict(M,data.frame(x1=-1000,x2=0,x3=0,y=NA)) # P is not negative...
print(P)

} # DELETE-OR-IGNORE-THIS LINE #

}

\keyword{classif}
\keyword{regression}
\keyword{neural}
\keyword{nonlinear}
