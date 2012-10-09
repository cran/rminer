# identical to e1071 package naiveBayes, except that I corrected some errors to get correct predictions
naivebayes <- function(x, ...)
  UseMethod("naivebayes")

naivebayes.default <- function(x, y, laplace = 0, ...) {
  call <- match.call()
  Yname <- deparse(substitute(y))

  ## estimation-function
  est <- function(var)
    if (is.numeric(var)) {
      cbind(tapply(var, y, mean, na.rm = TRUE),
            tapply(var, y, sd, na.rm = TRUE))
    } else {
      tab <- table(y, var)
      (tab + laplace) / (rowSums(tab) + laplace * nlevels(var))
    }
  
  ## create tables
  apriori <- table(y)
  tables <- lapply(x, est)
  
  ## fix dimname names
  for (i in 1:length(tables))
    names(dimnames(tables[[i]])) <- c(Yname, colnames(x)[i])
  names(dimnames(apriori)) <- Yname

  structure(list(apriori = apriori,
                 tables = tables,
                 levels = levels(y),
                 call   = call
                 ),
            
            class = "naivebayes"
            )
}

naivebayes.formula <- function(formula, data, laplace = 0, ...,
                               subset, na.action = na.pass) {
  call <- match.call()
  Yname <- as.character(formula[[2]])

  if (is.data.frame(data)) {
    ## handle formula
    m <- match.call(expand.dots = FALSE)
    m$... <- NULL
    m$laplace = NULL
    m$na.action <- na.action
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    if (any(attr(Terms, "order") > 1)) 
      stop("naivebayes cannot handle interaction terms")
    Y <- model.extract(m, "response")
    X <- m[,-attr(Terms, "response")]

    return(naivebayes(X, Y, laplace = laplace, ...))
  } else if (is.array(data)) {
    ## Find Class dimension
    Yind <- which(names(dimnames(data)) == Yname)

    ## Create Variable index
    deps <- strsplit(as.character(formula)[3], ".[+].")[[1]]
    if (length(deps) == 1 && deps == ".")
      deps <- names(dimnames(data))[-Yind]
    Vind <- which(names(dimnames(data)) %in% deps)
    
    ## create tables
    apriori <- margin.table(data, Yind)
    tables <- lapply(Vind,
                     function(i) (margin.table(data, c(Yind, i)) + laplace) / (as.numeric(apriori) + laplace * dim(data)[i]))

    structure(list(apriori = apriori,
                   tables = tables,
                   levels = names(apriori),
                   call   = call
                   ),
              
              class = "naivebayes"
              )
  } else stop("naivebayes formula interface handles data frames or arrays only")

}


print.naivebayes <- function(x, ...) {
  cat("\nNaive Bayes Classifier for Discrete Predictors\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nA-priori probabilities:\n")
  print(x$apriori / sum(x$apriori))
  
  cat("\nConditional probabilities:\n")
  for (i in x$tables) {print(i); cat("\n")}
    
}

predict.naivebayes <- function(object,
                               newdata,
                               type = c("class", "raw"),
                               threshold = 0.001,
                               ...) {
#object=NB@object
#newdata=D[(TRAIN+21),-1535]
#type="raw"
#threshold=0.001

  type <- match.arg(type)
  nattribs <- ncol(newdata)
#nattribs=10
  isnumeric <- sapply(newdata, is.numeric)
  newdata <- data.matrix(newdata)
  L <- sapply(1:nrow(newdata), function(i) {
    ndata <- newdata[i,]
    L <- log(object$apriori) + 
      apply(log(sapply(1:nattribs, function(v) {
        nd <- ndata[v]
        if(is.na(nd))
          rep(1, length(object$apriori))
        else {
          prob <- if (isnumeric[v]) {
            msd <- object$tables[[v]]
            msd[is.na(msd)]<-0 # new code
	    msd[,2][msd[,2]==0] <-  threshold
            dnorm(nd, msd[,1], msd[,2])
          } else
            object$tables[[v]][,nd]
          prob[prob == 0] <- threshold
          prob
        }
      })), 1, sum)
    if (type == "class")
      L
    else {
             # new code
             # warning: currently I do not know if this works when length(L)>2 ...
      if ( max(L)-min(L)> 600) {IM=which.max(L); L[IM]=300; L[-IM]=-300;} # avoid Inf or 0,0,0,...
      else if( sum(L< -745)==length(L) ) { Disp=min(L)+700; L=L-Disp;} # avoid L= 0, 0, 0, ...
      else if ( max(L)>709 ) {Disp=max(L)-700; L=L-Disp;} # avoid Inf
             # { IM=which.max(L); L[IM]=1; L[-IM]=0; L; } #1
             # end of new code
      L <- exp(L)# original
      L / sum(L) 
    }
  })
  if (type == "class")
    factor(object$levels[apply(L, 2, which.max)], levels = object$levels)
  else
    t(L)
}
