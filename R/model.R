#-------------------------------------------------------------------------------------------------
# "model.R" code by Paulo Cortez 2010@, Department of Information Systems, University of Minho
#
# This file deals will all Data Mining models
#
#  - Several parts of code were gratefully copied/inspired from the kernlab source code package :)
#  - The naiveBayes is a slightly different version of the e1071 package. It eliminates an error
#  that constantly appeared when a high number of inputs was used (spam e-mail classification)
#-------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------
# 17/11/2010 - deleted text mining and holdout sequence code. not sure if I will need this.
#-------------------------------------------------------------------------------------------------
# libraries that need to be installed:
library(nnet)   # get nnet: MLP and Multiple Logistic Regression
library(rpart)  # get rpart: Decision Trees
library(kknn)   # get kknn: k-nearest neighbours
library(kernlab)# get the ksvm

##library(neuralnet) # get neuralnet # still experimental stuff

# suggested libraries:
#try(library(MASS),silent=TRUE)
#try(library(mda),silent=TRUE)
#try(library(randomForest),silent=TRUE)

#----------------------------------------------------------------------------------------------
##library(klaR) # get NaiveBayes
##library(RWeka) # get Weka classifiers # uncomment if you need to use Weka classifiers
                 # uncomment all wmodels code and adapt to work...
#-------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------
# Things to do in the future: # - ensembles?
#-------------------------------------------------------------------------------------------------

# load other rminer files:
###source("estimate.R")
###source("preprocess.R")

#-------------------------------------------------------------------------------------------------
# save/load functions:
#

#-- save a fitted model into a binary (ascii=FALSE) of text (ascii=TRUE) file:
savemodel<-function(mm_model,file,ascii=FALSE)
{ #if(mmmodel@model=="wnaivebayes") { .jcache( (mmmodel@object)$classifier ) } # weka objects
  save(mm_model,file=file,ascii=ascii)
}
#-- load a file that was saved using savemodel function
loadmodel<-function(file)
{mm_model=FALSE;load(file=file); return(mm_model) }

#-- save a fitted mining into a binary (ascii=FALSE) of text (ascii=TRUE) file:
savemining<-function(mmm_mining,file,ascii=TRUE)
{ save(mmm_mining,file=file,ascii=ascii) }

#-- load a fitted mining file that was saved using savemining function
loadmining<-function(file)
{mmm_mining=FALSE;load(file=file);return(mmm_mining)}

#-------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------
# set the global object "model", contains the slots:
# @formula - the formula used to create the model
# @outindex - the output index of data
# @attributes - index with all attributes considered 
# @model - the data mining model: "naive", "lr", "mr", "dt", "naivebayes", "knn", "mlp", "svm" (optional: "lda", "qda", "randomforest", "bruto", "mars")
# @object - the data mining object (e.g. nnet, rpart, etc...) 
# @time - the estimation/training time
# @created - string with the date and time (strptime) when the model was build
# @mpar - vector with the parameters of the model (e.g. svm, mlp, knn)
#         other: metric
#         mlp: c(NR,MaxIt,validationmethod,validationpar,metric) 
#         mlp: c(NR,MaxIt,validationmethod,validationpar,metric,PAR) PAR = HN if decay is used for the search or DECAY if HN is used for the search
#         svm: c(C,Epsilon,validationmethod,validationpar,metric)  if C=NA or Epsilon=NA then heuristics are used
#         knn: c(validationmethod,validationpar,metric)
#              metric = "AUC", "ACC", "MAD", "SSE"
# @task - type of data mining task, one of the following: 
# "default", "class" (or "cla" or "classification"), "prob" (or "probability") - pure class, no probabilities!, "reg" (or "regression")
#     "default" means if is.factor(y) "class" else "reg"
# @scale - if the data needs to be scaled (to use specially with "mlp"). The options are:
#     "none", "inputs", "all", "default"
# @transform - if the output needs to be transformed ("reg"):
#     "log","logpositive","positive","scale"
# @levels - if the output is factor, store the levels. This is needed to recover the level order when predictive results are latter
#           analysed from files.
setClass("model",representation(formula="formula",model="character",task="character",mpar="ANY",attributes="numeric",scale="character",transform="character",created="character",time="numeric",object="ANY",outindex="numeric",levels="character")) 
#--------------------------------------------------------------------------------------
#--- generic and powerful fit method --------------------------------------------------
# most of these parameters are explained in the rminer.R file

#setGeneric("fit", function(x, ...) standardGeneric("fit"))
#setMethod("fit",signature(x="formula"),
fit=function(x, data=NULL, model="default",task="default", search="heuristic",mpar=NULL,feature="none",scale="default", transform="none",created=NULL,...)
{
 #print("<ini fit>")
 if(is.null(data)) 
   { data<-model.frame(x) # only formula is used 
     outindex=output_index(x,data)
     x<-as.formula(paste(names(data)[outindex]," ~ .",sep=""))
   }
 #else { if(feature[1]=="aries" && length(feature)>1) outindex=as.numeric(feature[2]) 
 else outindex=output_index(x,data) #}
 task=defaultask(task,model,data[1,outindex])
 if(model=="lr") model="logistic" else if(model=="default") model<-defaultmodel(task)
 metric=getmetric(mpar,model,task)
 mpar=defaultmpar(mpar,model,task,metric)

 if(task!="reg" && !is.factor(data[,outindex])) data[,outindex]=factor(data[,outindex])
 else if(task=="reg" && transform_needed(transform)) data[,outindex]=xtransform(data[,outindex],transform) #transform output?

 if(is.factor(data[,outindex])) levels<-levels(data[1,outindex]) else levels=""
 params=NULL

 PTM<- proc.time() # start clock

 feature=defaultfeature(feature) # process feature
 if(feature_needed(feature[1]))
  { 
#print(feature)
#print("--")
    fstop=as.numeric(feature[2]) # -1 or number of delections
    Runs=as.numeric(feature[3]) 
    vmethod=feature[4] 
    vpar=as.numeric(feature[5]) 
    if(length(feature)>5) { fsearch=suppressWarnings(as.numeric(feature[6])) # primary search parameter
                            if(is.na(fsearch)) fsearch=feature[6]
                          }
    else fsearch=search # this may take a long time...
       # perform both FS and parameter search at the same time...
    RES=bssearch(x,data,algorithm=feature[1],Runs=Runs,method=c(vmethod,vpar),model=model,task=task,search=fsearch,mpar=mpar,scale=scale,transform="none",fstop=fstop)
    attributes=RES$attributes
    outindex=output_index(x,data[,attributes]) # update the output index if it has changed due to the deletion of variables
    if(length(fsearch)>1) { LRS=length(RES$search) 
                            if(substr(model,1,3)=="mlp")
                            { if(length(mpar)==6) search=RES$search[2] # decay
                              else search=RES$search[1] # Hn
                            }
                            else{ search=RES$search[1];
                                  if(LRS>1) mpar[1:(LRS-1)]=RES$search[2:LRS]
                                }
                          }
    data=data[,attributes] # new
  }
 else attributes=1:NCOL(data); # set to all

 # --- if for all models ---
 if (model=="naive")
  { if(substr(task,1,3)=="reg") M<-mean(data[,outindex])
    else if(is.factor(data[,outindex])) 
    { M<-data[1,outindex]; M[1]<-levels(data[1,outindex])[mostcommon(data[,outindex])] }
  }
 else if(model=="dt") # decision tree: based in the rpart library
  { if(task=="reg") METHOD="anova" else METHOD="class"
    M=rpart(x,data=data,method=METHOD,...)
  }
 #else if(model=="randomforest") # random forest: based in the randomForest package
 # {
   #M<-randomForest(x,ntree=1000,data=data,importance=TRUE,...)
 #  M<-randomForest(x,data=data,importance=TRUE,mtry=search,...)
 # }
 else if(model=="mr") # multiple/linear regression: uses the nnet function
  { M=mlp.fit(x,data,HN=0,NR=1,maxit=100,task=task,scale="none") }
 else if(model=="mars") # mars
  { M=mars(data[,-outindex],data[,outindex]) }
 else if(model=="bruto") # bruto
  { M=bruto(data[,-outindex],data[,outindex]) }
 #else if(model=="lm") # multiple/linear regression: uses the R lm function
 # {
 #  if(substr(task,1,3)!="reg" || is.factor(data[,outindex])) 
 #	print("Error: lm should be used in regression with a numeric output!!!")
 #  suppressWarnings(M<-lm(x, data=data,...))
 # }
 else if(model=="lda" || model=="qda") # linear/quadratic discriminant analysis
  {
   if(substr(task,1,3)=="reg") cat("Error:",model,"should be used in classification\n")
   if(model=="lda") suppressWarnings(M<-lda(x,data=data,...)) else suppressWarnings(M<-qda(x,data=data,...))
  }
 else if(model=="naivebayes") # || model=="wnaivebayes") # naive bayes: based in the naiveBayes function
  {
   if(substr(task,1,3)=="reg") print("Error: naive bayes should be used in classification\n")
   #if(model=="wnaivebayes") { WNB <- make_Weka_classifier("weka/classifiers/bayes/NaiveBayes"); M<-WNB(x,data=data,...); }
   #else 
   M<-naivebayes(x,data=data,...)
   #if(length(attributes)==2) M<-NULL # model is equal to the unique input attribute
   #else M<-naiveBayes(x,data=data,...)
 }
 else if(model=="logistic") 
  {
   M<-multinom(x,data=data,trace=FALSE,...)
   #else suppressWarnings(M<-glm(x,data=data,family=binomial(logit),...)) # glm gives some errors when factors have levels with few elements.
  }
 else if(model=="knn" || substr(model,1,3)=="mlp" || model=="svm" || model=="randomforest") # methods with internal search for best parameter
  {# --- else if ----------------------------------------------------------------------------------------
   if(scale=="default" && substr(model,1,3)=="mlp") # removed "svm" from this if due to a "object scal not found" error
     { if(is.factor(data[,outindex])) scale="inputs" else scale="all" }
   # set default mpar and metric
#print("mpar:")
#print(mpar)
   if(substr(model,1,3)=="mlp" && length(mpar)==6) decay=TRUE else decay=FALSE
   search=defaultsearch(search,model,task,COL=NCOL(data),decay)
#print("search:")
#print(search)
   if(is.list(search) && length(search$convex)==1) convex=search$convex else convex=0 # perform convex search: much faster and simpler...

#cat("convex:",convex,"\n")
   if( (is.list(search)&&length(search$search)==1) || length(search)==1) K=search # do not perform an expensive best fit search
   else if(!is.list(search)) K=bestfit(x,data,model,task,scale,search,mpar,transform="none",convex=convex,...) # simple grid search
   else # is list 
   { 
          if(length(search$smethod)==1) smethod=search$smethod else smethod="normal"
          if(smethod=="normal") K=bestfit(x,data,model,task,scale,search=search$search,mpar,transform="none",convex=convex,...)
          else if(substr(smethod,1,2)=="UD" || smethod=="2L") 
           {
             LS=length(search$search)
             if(substr(smethod,1,2)=="UD"){if(smethod=="UD1") points=13 else points=9
                                           factors=LS/2;upar=2^(uniform_design(search$search,factors=factors,points=points))} # (gamma,C) or # (gamma,C,epsilon)
             else upar=search$search # vector, 2L
             K=bestfit(x,data,model,task,scale,search=upar,mpar,transform="none",error=TRUE,convex=convex,...)
             KI=log(K$K,2)
             # 2nd level:
             if(substr(smethod,1,2)=="UD")
                              {upar=rep(0,LS);if(smethod=="UD1") points=9 else points=5
                               for(i in 1:factors) {ini=i*2-1;end=ini+1;aux=diff(range(search$search[ini:end]))/4;
                                                    upar[ini]=KI[i]-aux; upar[end]=KI[i]+aux
                                                   }
                               upar=uniform_design(upar,points=points,center=KI)
                               upar=2^upar
                              }
             else { # vector, 2L
                    KI=which(search$search==K$K[1])
                    if(model=="svm") {upar=log(search$search,2);KK=log(K$K[1],2);} else {upar=search$search;KK=K$K[1];}
                    if(KI==1) Range=diff(range(upar[1:2])) else if(KI==LS) Range=diff(range(upar[(LS-1):LS])) else Range=diff(range(upar[(KI-1):(KI+1)]))/2
                    upar=seq(upar[KI]-Range/2,upar[KI]+Range/2,length.out=LS)
                    upar=upar[setdiff(1:length(upar),which(upar==KK[1]))]
                    if( (substr(model,1,3)=="mlp" && length(mpar)<6)||model=="randomforest"||model=="knn") upar=unique(round(upar)) # H, mtry, k

                    if(decay) upar=upar[which(upar<=1)]
                    if(substr(model,1,3)=="mlp") upar=upar[which(upar>=0)] # H or decay
                    else if(model=="randomforest"||model=="knn") { upar=upar[which(upar>=1)]} # mtry, k 
                    else if(model=="svm") upar=2^upar # gamma = sigma
                  }
             if(length(upar)>0) K2=bestfit(x,data,model,task,scale,search=upar,mpar,transform="none",error=TRUE,convex=convex,...)
             if(length(upar)>0 && isbest(K2$error,K$error,metric)) K=K2$K else K=K$K
           }
   }
   if(model=="knn") # k-nearest neighbour, uses kknn implementation
   { M<-knn.fit(x=x,data=data,k=K); params=c(K) }
   else if(substr(model,1,3)=="mlp")
   {
    if(model=="mlpe") type=2
    else if(model=="mlp2") type=3
    else type=1
    # read the internal parameters
    NR<-as.numeric(mpar[1]); maxit<-as.numeric(mpar[2]);
    if(decay) {hn=as.numeric(mpar[6]); DECAY=K} else {DECAY=0.0;hn=K}
    #print(paste("NR:",NR," maxit:",maxit))
    M<-mlp.fit(x,data,HN=hn,decay=DECAY,NR=NR,maxit=maxit,task,scale,type=type,metric=metric)
    params=c(M$mpar)
   }
   else if(model=="svm") # support vector machine, uses kernlab implementation
   { 
    # read the internal parameters
    sigma=K[1];
    if(length(K)>2)
      { C=K[2];epsilon=K[3] }
    else if(length(K)>1)
      { C=K[2];epsilon=NA}
    else
      { 
       LP<-length(mpar)
       if(LP>=1 && !is.na(mpar[1])) C=as.numeric(mpar[1]) else C=NULL;
       if(LP>=2 && !is.na(mpar[2])) epsilon=as.numeric(mpar[2]) else epsilon=NULL;
      }
     M<-svm.fit(x,data,sigma=sigma,C=C,epsilon=epsilon,task,scale)
     if(substr(task,1,3)=="reg") epsilon=M$mpar[3] else if(task=="class") epsilon=0 else epsilon=1
     if(is.null(C)) C=M$mpar[2]
     params=c(sigma,C,epsilon)
   }
   else if(model=="randomforest")
   {
    M<-randomForest(x,data=data,importance=TRUE,mtry=K,...);params=K
   }
  } # --- end else if -----------------------------------------------------------------------------------
 TME<-(proc.time()-PTM)[3] # computes time elapsed, in seconds
 if(is.null(created)) created=format(Sys.time(),"%Y-%m-%d %H:%M:%S")
#print(" < end fit>")
 if(!is.null(params)) {params=data.frame(t(params));names(params)=modelnamespar(model);}
 return(new("model",formula=x,model=model,task=task,mpar=params,attributes=attributes,scale=scale,transform=transform,created=created,time=TME,object=M,outindex=outindex,levels=levels))
}#)

# 
defaultsearch=function(search,model,task,COL,decay=FALSE)
{
 if(class(search)=="character")
 { if(search=="heuristic") 
   { if(model=="knn") search=3 # 3nn by default
     else if(model=="randomforest") { if(task!="reg") search=floor(sqrt(COL-1)) else search=max(floor((COL-1)/3),1) }
     else if(substr(model,1,3)=="mlp") 
      { search=round((COL-1)/2) # simple heuristic for hidden nodes
        if(task!="reg" && search>10) search=10 # classification mlp cannot have huge search, as R entropy minimization gives an error...
      }
     else if(model=="svm") search=2^-7 #"automatic" # use automatic heuristic for the kernel parameter
   }
   else if(search=="heuristic5") # simple heuristic
     { if(model=="knn") search=seq(1,9,2) else if(substr(model,1,3)=="mlp"){if(decay) search=seq(0,0.2,0.05) else search=seq(0,8,2)} else if(model=="svm") search=2^seq(-15,3,4)
       else if(model=="randomforest") search=1:5 
     }
   else if(search=="heuristic10")
     { if(model=="knn") search=seq(1,10,1) else if(substr(model,1,3)=="mlp"){if(decay) search=seq(0,0.2,length.out=10) else search=seq(0,9,1)} else if(model=="svm") search=2^seq(-15,3,2)
       else if(model=="randomforest") search=1:10 
     }
   else if(substr(search,1,2)=="UD" && model=="svm")
     {
       #if(task=="reg") search=c(-8,0,-1,6,-8,-1) # gama, C, epsilon in libsvm guide?
       smethod=search
       if(task=="reg") search=c(-8,0,-1,6,-8,-1) else search=c(-15,3,-5,15) # gama, C, epsilon
       #if(task=="reg") search=c(-15,3,-5,15,-11,2) else search=c(-15,3,-5,15) # gama, C, epsilon
       search=list(smethod=smethod,search=search)
     }
 }
 return (search)
}

defaultmpar=function(mpar,model,task,metric)
{ DLP=switch(model,knn=,randomforest=3,mlp=,mlpe=,svm=5,1)
  LP=length(mpar)
  if(LP>=DLP && sum(is.na(mpar))==0) return (mpar)
  else
  { dmpar=switch(model,knn=,randomforest=c("holdout",2/3,metric),
                    mlp=,mlpe=c(3,100,"holdout",2/3,metric), # NR, maxit,...
                    svm=c(NA,NA,"holdout",2/3,metric), # C, epsilon,...
                    metric)
    if(LP<DLP) mpar=c(mpar,rep(NA,(DLP-LP)))
    for(i in 1:DLP) { if(is.na(mpar[i])) mpar[i]=dmpar[i]
                      if(i>1 && !is.na(mpar[i-1]) && mpar[i-1]=="kfold" && !is.na(mpar[i]) && mpar[i]==2/3) mpar[i]=3
                    }
    return (mpar)
  }
}

defaultfeature=function(feature)
{ if(!is.na(feature[1]))
   { if(feature[1]=="s") return("simp")
     else if(substr(feature[1],1,4)=="sens" || substr(feature[1],1,4)=="simp" || feature[1]=="none") return(feature[1]) 
   }
  else if(is.na(feature[1]) && length(feature)==1) return("none")
  LF=length(feature);DLF=5
  if(LF<DLF || sum(is.na(feature))!=0) 
  {
   dfeature=c("sabs",-1,1,"holdout",2/3)
   for(i in 1:DLF) { if(is.na(feature[i])) feature[i]=dfeature[i]
                     if(i>1 && !is.na(feature[i-1]) && feature[i-1]=="kfold" && !is.na(feature[i]) && feature[i]==2/3) feature[i]=3
                   }
  }
  return(feature)
}

# ====== Feature Selection Methods (beta development =======================
bssearch=function(x,data,algorithm="sabs",Runs,method,model,task,search,mpar,scale,transform,smeasure="g",fstop=-1,debug=FALSE)
{ 
 #debug=TRUE; cat("debug:",debug,"\n")
 metric=getmetric(mpar,model,task)
 OUT=output_index(x,data) # output
 Z=1:NCOL(data);BZ=Z;
 LZ=length(Z); if(fstop==-1)fstop=(LZ-2);LZSTOP=min(LZ-2,fstop)
 t=0 # iteration
 NILIM=LZ #convex 0 # LZ # 2 # I need to program this better, transform this parameter into a input user parameter!!!
 JBest=worst(metric);BK=search; # worst error
 notimprove=0 # iterations with improvement
 if(algorithm=="sabs") imethod=switch(smeasure,g="sensg",v="sensv",r="sensr")
#cat("imethod:",imethod,"\n") 
#print(BK)
 stop=FALSE;t=0;
 while(!stop)
 {
   if(algorithm=="sabs") 
   { 
    M=mining(x,data[,Z],model=model,task=task,method=method,feature=c(imethod,1,"all"),search=search,mpar=mpar,scale=scale,transform=transform,Runs=Runs)
    J=mean(M$error); K=medianminingpar(M)
    LZ=length(Z);Imp=vector(length=LZ);for(i in 1:LZ) Imp[i]=mean(M$sen[,i]);
   }
   else if(algorithm=="sbs")
   { if(t==0) { M=mining(x,data[,Z],Runs=Runs,model=model,method=method,task=task,feature="none",search=search,mpar=mpar,scale=scale,transform=transform)
               J=mean(M$error); K=medianminingpar(M)
              }
   }
   if(isbest(J,JBest,metric)) {JBest=J;BZ=Z;notimprove=0;BK=K;BT=t;} else notimprove=notimprove+1;
   if(debug){ cat(">> t:",t,"cerr:",J,"Best:",JBest,"BT:",BT,"Z:",Z,"S:",K,"\n",sep=" ");}
   if((t==LZSTOP || notimprove==NILIM)) stop=TRUE
   else # select next attribute to be deleted 
    { 
     if(algorithm=="sabs"){IOUT=which(Z==OUT);Zimp=setdiff(1:length(Z),IOUT);IDel=which.min(Imp[Zimp]);Del=Zimp[IDel]; Z=setdiff(Z,Del) }
     else if(algorithm=="sbs"){
                               J=worst(metric);Zb=Z;
                               for(i in setdiff(Z,OUT))
                               {
                                Zi=setdiff(Z,i)
                                M=mining(x,data[,Zi],Runs=Runs,model=model,method=method,task=task,feature="none",search=search,mpar=mpar,scale=scale,transform=transform)
                                Ji=mean(M$error); Ki=medianminingpar(M); 
                                if(isbest(Ji,J,metric)) {J=Ji;K=Ki;Zb=Zi}
                                #cat("   >> i:",i,"Zi",Zi,"Ki:",Ki,"err:",Ji,"best:",J,"\n")
                               }
                               Z=Zb
                              }
    }
   t=t+1
 }
 if(debug) cat(" | t:",BT,"Best:",JBest,"BZ:",BZ,"BS:",BK,"\n",sep=" ")
 return(list(attributes=BZ,search=BK)) # attributes and search parameter
}
# ====== End of Feature Selection Methods ==================================

# fast uniform design: only works for some uniform design setups: 2, 3 factors and 5,9,13 points
uniform_design<-function(limits,factors=length(limits)/2,points,center=NULL)
{
 ud=matrix(nrow=points,ncol=factors) # gamma, C (and epsilon?)
 
 SEQ=vector("list",length=factors)
 for(i in 1:factors)
  { ini=2*(i-1)+1;SEQ[[i]]=seq(limits[ini],limits[ini+1],length=points) }

# 2, 5 
if(factors==2)
{ if(points==5) m=t(matrix(c(1,2,2,5,4,1,5,4,3,3),ncol=5)) 
  else if(points==9) m=t(matrix(c(5,5,1,4,7,8,2,7,3,2,9,6,8,3,6,1,4,9),ncol=9))
  else if(points==13) m=t(matrix(c(5,4,12,3,2,11,9,10,7,7,6,13,3,2,11,12,13,8,10,5,1,6,4,9,8,1),ncol=13))
}
else if(factors==3)
{ if(points==5) m=t(matrix(c(5,3,3,4,4,5,3,1,1,2,5,2,1,2,4),ncol=5))
  else if(points==9) m=t(matrix(c(3,9,6,9,4,7,7,1,4,2,2,8,1,6,3,8,8,2,4,3,1,5,5,5,6,7,9),ncol=9))
  else if(points==13) m=t(matrix(c(9,9,1,7,4,5,8,13,10,2,3,11,4,6,2,13,7,9,6,1,8,10,5,12,12,11,6,3,12,4,1,8,7,5,10,13,11,2,3),ncol=13))
}
for(i in 1:nrow(m))
for(j in 1:ncol(m))
 {
  val=SEQ[[j]][m[i,j]]
  ud[i,j]=val
 }
if(length(center)>0) # delete center point
{ ALL=1:NROW(ud)
  I=which(ud[ALL,1]==center[1])
  ud=ud[setdiff(ALL,I),]
} 
return (ud)
}

# auxiliary create grid for classification and regression as suggested by 
# Hastie et al 2001 and the LIBSVM authors
svmgrid<-function(task="prob")
{
 l=1
 if(substr(task,1,3)=="reg")
   {
     res<-matrix(ncol=3,nrow=576)
     for(i in 0:-8)
     	for(j in -1:6)
     		for(k in -8:-1) 
			{ res[l,1]=2^i; res[l,2]=2^j; res[l,3]=2^k; l=l+1;}
   }
 else # "class" || "prob"
   { res<-matrix(ncol=2,nrow=110)
     for(i in -15+(0:9)*2)
     	for(j in -5+(0:10)*2)
		{ res[l,1]=2^i; res[l,2]=2^j; l=l+1}
   }
 return (res)
}

# --- the best fit internal function ---------------------------------------------------
bestfit=function(x,data,model,task,scale,search,mpar,transform,convex=0,error=FALSE,imp="none",...)
{
 #cat("Imp:",imp,"\n")
#print(mpar)
 outindex<-output_index(x,data)
 y<-data[,outindex]

 if(model=="randomforest") search=search[which(search<=(NCOL(data)-1))] # clean elements higher than inputs
 VEC_SEARCH<-is.vector(search)

 if(VEC_SEARCH)
   { NP<-length(search); BEST<-search[1]; }
 else # matrix, used for svm or models that require more than 1 hyperparameter
   { NP<-nrow(search) # one parameter per line
     BEST<-search[1,] # the first line 
   }
 P<-vector(length=0)
 
 LP<-length(mpar)
 DECAY=0.0

 metric=getmetric(mpar,model,task)

 if( substr(model,1,3)=="mlp" || model=="svm") 
   {K1=as.numeric(mpar[1]); K2=as.numeric(mpar[2]); vmethod=mpar[3]; K3=as.numeric(mpar[4]);
    if(LP>5 && substr(model,1,3)=="mlp") {DECAY=mpar[6]}
   }
 else {vmethod=mpar[1]; K3=as.numeric(mpar[2]);}

 BESTVAL=worst(metric)

 stop=FALSE; i=1; # for(i in 1:NP)
 if(convex>0) NLIM=convex
 else NLIM=(NP+1)
 notimprove=0

 if(substr(imp,1,4)=="simp") { FSRESP=TRUE;imp="sens"}
 else {FSRESP=FALSE}

 if(substr(imp,1,4)=="sens")
   {
     if(vmethod=="kfold") MULT=K3
     else MULT=1
     SEN<-matrix(ncol=NCOL(data),nrow=MULT)
   }
 else SEN=NULL

 if(substr(vmethod,1,6)=="kfoldo") order=TRUE
 else order=FALSE

 while(!stop)
   {
     if(vmethod=="all")
       {
          if(VEC_SEARCH) M<-fit(x,data,model=model,task=task,scale=scale,search=c(search[i]),mpar=mpar,transform=transform,...)
          else
            { if(model=="svm") { if(NCOL(search)>=2) mpar[1]<-search[i,2]
                                 if(NCOL(search)>=3) mpar[2]<-search[i,3]
                               }
              M<-fit(x,data,model=model,task=task,scale=scale,search=c(search[i,1]),mpar=mpar,transform=transform,...)
              if(!is.null(SEN)) { 
                           SEN[1,]<-Importance(M,data,method=imp,responses=FSRESP)$imp # store sen
                         }
            }
          P<-predict(M,data)
          TS<-y
       }
     else if(substr(vmethod,1,7)=="holdout") #all holdout types 
       {  
          if(vmethod=="holdoutorder") H<-holdout(y,ratio=K3,mode="order") 
         else if(vmethod=="holdoutfor5") H<-holdout(y,ratio=K3-5,mode="order") 
          else H<-holdout(y,ratio=K3)

          if(VEC_SEARCH) M<-fit(x,data[H$tr,],model=model,task=task,scale=scale,search=c(search[i]),mpar=mpar,transform=transform,...) 
          else
            { if(model=="svm") { if(NCOL(search)>=2) mpar[1]<-search[i,2]
                                 if(NCOL(search)>=3) mpar[2]<-search[i,3]
                               }
              M<-fit(x,data[H$tr,],model=model,task=task,scale=scale,search=c(search[i,1]),mpar=mpar,transform=transform,...) 
            }
          if(!is.null(SEN)) { SEN[1,]<-Importance(M,data[H$tr,],method=imp,responses=FSRESP)$imp # store sen
                         }
          TS<-y[H$ts]
       #if(method=="holdoutfor") P=lforecast(M,data,start=(1+nrow(data)-K3),horizon=K3) 
          #else if(method=="holdoutfor5")
          # {  Eval=vector(length=5)
          #    for(kkk in 1:5)
          #    {
          #      beg=1+nrow(data)-K3-5+kkk;
          #      P=lforecast(M,data,start=beg,horizon=K3) 
          #      Eval[kkk]=Error(data[beg:(beg+horizon),ncol(data)],P,metric)
          #    }
          #    Eval=mean(Eval)
          # }
        #else 
          P<-predict(M,data[H$ts,])
       }
     else if(substr(vmethod,1,5)=="kfold")
       {
          if(model=="knn") kmpar=c("all",0,metric)
          else kmpar=c(K1,K2,"all",0,metric)
          if(DECAY>0) kmpar<-c(kmpar,DECAY)
         
          if(VEC_SEARCH) KF<-crossvaldata(x,data,fit,predict,ngroup=K3,order=order,model=model,task=task,feature=imp,
                                          scale=scale,search=c(search[i]),mpar=kmpar,transform=transform,...)
          else
            { if(model=="svm") { if(NCOL(search)>=2) kmpar[1]<-search[i,2]
                                 if(NCOL(search)>=3) kmpar[2]<-search[i,3]
                               }
              
              KF<-crossvaldata(x,data,fit,predict,ngroup=K3,order=order,model=model,task=task,feature=imp,
                                          scale=scale,search=c(search[i,1]),mpar=kmpar,transform=transform,...)
            }
          if(!is.null(SEN)) { 
                           #print(paste("B:",Begi,"E:",Endi,"Mult:",MULT,"ncol:",ncol(SEN)))
                           #print(KF$sen)
                           SEN<-KF$sen # store sen
                         }

          P<-KF$cv.fit
#cat(" <---- i:",i,"----\n")
#print(summary(P))
          TS<-y
       }
     # get the error
       #if(method!="holdoutfor5") 
       Eval=Error(TS,P,metric)
       if(isbest(Eval,BESTVAL,metric)) {BESTVAL<-Eval;notimprove=0;if(VEC_SEARCH) BEST<-search[i] else BEST<-search[i,]} 
       else notimprove=notimprove+1
#cat("i",i,"at:",(ncol(data)-1),"nr:",nrow(data),"s:",search[i],"val:",Eval,"b:",BESTVAL,"best:",BEST,"\n")
###cat(names(data),"\n")
#print(M@object)
       i<-i+1
       if(notimprove==NLIM) stop=TRUE
       else if(i==NP+1) stop=TRUE
    } # end while 
 #print(paste("Best:",BEST))
 # refit the best model with all training data: 
 if(error) return(list(K=BEST,error=BESTVAL,sen=SEN))
 else return (BEST)
}
# -------- end of best fit -----------

output_index<-function(x,data=NULL)
{
 T<-terms(x,data=data)
 Y<-names(attr(T,"factors")[,1])[1] 
 return(which(names(data)==Y))
}

transform_needed=function(transform) { return (switch(transform,log=,logpositive=,positive=,scale=TRUE,FALSE)) }
#param_needed=function(model) { return (switch(model,mlp=,mlpe=,svm=,knn=,randomforest=TRUE,FALSE)) }
feature_needed=function(feature) { return (switch(feature,sbs=,sabs=TRUE,FALSE)) }

defaultask=function(task="default",model="default",output=1)
{ if(task=="default") { task=switch(model,logistic=,lr=,lda=,qda=,naivebayes="prob",mr=,mars=,bruto="reg","default")
                        if(task=="default") { if (is.factor(output)) task="prob" else task="reg" }
                      }
  else if(substr(task,1,1)=="c") task="class" else if(substr(task,1,1)=="p") task="prob"
  return (task)
}

defaultmodel<-function(task)
{ 
  if(substr(task,1,3)=="reg") return ("mr")
  else return ("dt")
}
#--- end of generic fit ------------------------------------------------------

#--- generic predict method --------------------------------------------------
# @@@ search pattern for this function
# object - a model created with fit
# newdata - data.frame or matrix or vector or factor
#           if formula, please use: predict(...,model.frame(data))
setMethod("predict",signature(object= "model"),
function(object,newdata){

 if(NCOL(newdata)>length(object@attributes)) newdata=newdata[,object@attributes]
 
 # --- if code for the models ---
 if(substr(object@model,1,3)=="mlp"){ if(object@scale=="inputs" || object@scale=="all") newdata=scaleinputs2(newdata,object@object$cx,object@object$sx) }
 if(object@model=="naive")
   {
    if(object@task=="reg") P=rep(object@object,length=nrow(newdata))
    else { L=levels(object@object[1])
           P<-matrix(nrow=nrow(newdata),ncol=length(L)); P[,]=0
           P[,which(L==object@object)]=1
           if(task=="class") P=majorClass(P,L)
         }
   }
 else if(object@model=="dt")
   { type=switch(object@task,class="class",prob="prob","vector")
     P=predict(object@object,newdata,type=type)
   }
 else if(object@model=="randomforest")
   { type=switch(object@task,prob="prob","response")
     P=predict(object@object,newdata,type=type)
     #attr(P,"class")=NULL
   }
 else if(object@model=="mr")
    { P<-predict(object@object$mlp,newdata)[,1] }
 else if(object@model=="mars")
    { P<-predict(object@object,newdata[,-object@outindex])[,1] }
 else if(object@model=="bruto")
    { P<-predict(object@object,as.matrix(newdata[,-object@outindex]))[,1] }
 #else if(object@model=="lm") { suppressWarnings(P<-predict(object@object,newdata)) }
 else if(object@model=="lda" || object@model=="qda")
   { P<-predict(object@object,newdata)
     if(object@task=="class") P=P$class else P=P$posterior 
   }
 else if(object@model=="naivebayes") # || object@model=="wnaivebayes")
   { if(is.null(object@object)) P<-newdata[,-object@outindex] # only 1 input
     else { type=switch(object@task,class="class","raw")
            suppressWarnings(P<-predict(object@object,newdata[,-object@outindex],type=type))
          }
     #else  # wnaivebayes
     #{ suppressWarnings(P<-predict(object@object,newdata[,-object@outindex],type="prob")) }
   }
 else if(object@model=="logistic")
   { type=switch(object@task,class="class","probs")
     P=predict(object@object,newdata,type=type)
     if(length(object@levels)==2 && object@task=="prob") {MM=matrix(ncol=2,nrow=length(P));MM[,2]=P;MM[,1]=1-P;P=MM;}
     #else P<-1/(1+exp( - predict(object@object,newdata)) )
   }
 else if(object@model=="knn") # k-nearest neighbour, uses kknn implementation
   { if(object@task=="reg") P=(kknn(object@object$x,object@object$data,newdata,k=object@object$k))$fitted.values
     else { P=(kknn(object@object$x,object@object$data,newdata,k=object@object$k))
            P=switch(object@task,prob=P$prob,P$fitted.values)
          }
   }
 else if(object@model=="mlp") 
   { if(object@task=="reg") 
       { P<-predict(object@object$mlp,newdata)[,1]
         if(object@scale=="all") P<-invtransform(P,"scale",object@object$cy,object@object$sy)
       }
     else
      { type=switch(object@task,class="class","raw")
        P<-predict(object@object$mlp,newdata,type=type)
        if(length(object@levels)==2 && object@task=="prob") {MM=matrix(ncol=2,nrow=length(P));MM[,2]=P;MM[,1]=1-P;P=MM;}
      }
   }
 else if(object@model=="mlpe") # 
   {  MR=object@mpar[1,3]
      if(object@task!="reg"){L=length(object@levels); P=matrix(0,ncol=L,nrow=NROW(newdata))} 
      else P=rep(0,NROW(newdata))
      for(i in 1:MR)
      {
       if(object@task=="reg") 
        { P1=predict(object@object$mlp[[i]],newdata)[,1]
          if(object@scale=="all") P1<-invtransform(P1,"scale",object@object$cy,object@object$sy)
          P=P+P1
        }
       else
       { 
        P1=predict(object@object$mlp[[i]],newdata)
        if(L==2) {MM=matrix(ncol=2,nrow=length(P1));MM[,2]=P1;MM[,1]=1-P1;P1=MM;}
        P=P+P1
       }
      }
      P=P/MR
      if(object@task=="class") P=majorClass(P,L=object@levels)
   }
 #else if(object@model=="mlp2") # in development, does not work for class and prob...
 #  {
 #     if(object@task=="reg") 
 #       { 
 #         P<-compute(object@object$mlp,newdata[-object@outindex])$net.result
 #         if(object@scale=="all") P<-invtransform(P,"scale",object@object$cy,object@object$sy)
 #       }
 #     else # classification: to be done...
 #     { L=length(object@levels)
 #       P<-predict(object@object$mlp,newdata)
 #       if(L==2) {MM=matrix(ncol=2,nrow=length(P));MM[,2]=P;MM[,1]=1-P;P=MM;}
 #     }
 #  }
 else if(object@model=="svm") 
   {
      if(object@task=="prob") P=predict(object@object$svm,newdata,type="probabilities")
      else if(object@task=="class") P=predict(object@object$svm,newdata) 
      else P<-predict(object@object$svm,newdata,type="response")[,1]
   }
 if(substr(object@task,1,3)=="reg" && transform_needed(object@transform))P=invtransform(P,object@transform)

 #if(object@task=="reg") attr(P,"names")=NULL # new line: should i keeep this? lets see...
 #else attr(P,"dimnames")=NULL # classification

 return(P)
})

#-----------------------------------------------------------------------------
modelnpar<-function(model)
{return(switch(model,knn=,randomforest=1,mlp=,mlpe=4,svm=3,0))}
modelnamespar<-function(model)
{return(switch(model,knn="k",randomforest="mtry",mlp=,mlpe=c("H","decay","Nr","Me"),svm=c("gamma","C","epsilon"),""))}

knn.fit<- function(x,data,k)
{ return (list(x=x,data=data,k=k)) #why? kknn does not have predict function that works like predict for other methods
}

svm.fit <- function(x,data,sigma=2^-6,C=NULL,epsilon=NULL,task,scale=TRUE)
{ 
 # to avoid this error: object "scal" not found, I will adopt SCALED=TRUE
 SCALED=TRUE; cx=NULL; sx=NULL; 
#print(paste("svm.fit Sigma: ",sigma," C: ",C," Ep:",epsilon))
 #if(scale=="all") SCALED=TRUE # XXX warning: be careful about this option. scaled by rminer of ksvm? 
 # check if ksvm has some problems with numeric output??
 #else SCALED=FALSE

 outindex=output_index(x,data)
 #if(scale=="inputs") # scale the inputs 
 #  { S<-scaleinputs(data,outindex)
 #    data=S$data; cx=S$cx; sx=S$sx
 #  }
 #else{cx=NULL; sx=NULL}

 # this code is needed for the SVM heuristics only:
 if(substr(task,1,3)=="reg") # && (scale=="output" || scale=="all")) # scale the output
   { CY=mean(data[,outindex])
     SY=sd(data[,outindex])
     Y<-xtransform(data[,outindex],"scale",A=CY,B=SY)
   }
 else{CY=0; SY=0; Y=data[,outindex]}
 
 if(is.null(C)) # set the C value
  {
   if(is.factor(Y)) { #Y<-as.numeric(Y) # reasonable ???
                      if(FALSE){
                      L=levels(Y);NL=length(L)
                      if(NL==1) { MEANY=0;SDY=1;}
                      else if(NL==2) 
                                 { NY=as.numeric(Y==L[2])
                                   MEANY=mean(NY);SDY=sd(NY)
                                 }
                      else { MEANY=vector(length=NL);SDY=MEANY;
                             for(i in 1:NL) 
                                 { NY=as.numeric(Y==L[i])
                                   MEANY[i]=mean(NY);SDY[i]=sd(NY)
                                 }
                           }
                               }
                      MEANY=0;SDY=1; 
                    }
   else
   {
    MEANY<-mean(Y)
    SDY<-sd(Y)
   }
   C<-max(abs(MEANY+3*SDY),abs(MEANY-3*SDY))
   #cat("Scaled:",SCALED,"MEANY",MEANY,"SDY:",SDY,"C:",C,"\n")
  }
 
 if(is.null(epsilon) && substr(task,1,3)=="reg") # set the epsilon, in regression
 {
  data2=data; data2[,outindex]=Y
  N=nrow(data)
  
  K3<-(kknn(x,train=data2,test=data2,k=3))$fitted.values # 3NN, as argued by CheMa04
  SUM<-sum((Y-K3)^2)
  SD=sqrt(1.5/N * SUM)
  #print(paste("SD: ",SD))
  if(N>100) epsilon<-3*SD*sqrt(log(N)/N) # eq.14 of CheMa04: use for large N 
  else epsilon=SD/sqrt(N) # eq.13 of CheMa04 : ok if N small 
 }

#cat("s:",sigma,"C:",C,"e:",epsilon,"prob:",prob,"\n")
 if(task=="prob") { SVM=ksvm(x,data=data,scaled=SCALED,kernel="rbfdot",kpar=list(sigma=sigma),C=C,prob.model=TRUE); mpar=c(sigma,C,1) }
 else if(task=="class") { SVM=ksvm(x,data=data,scaled=SCALED,kernel="rbfdot",kpar=list(sigma=sigma),C=C,prob.model=FALSE); mpar=c(sigma,C,0) }
 else # regression
  {
   SVM=ksvm(x,data=data,scaled=SCALED,kernel="rbfdot",kpar=list(sigma=sigma),C=C,epsilon=epsilon)
   mpar=c(sigma,C,epsilon)
  }
if(is.null(epsilon)) epsilon<-NA
#return(list(svm=SVM,mpar=mpar,cx=cx,sx=sx,cy=CY,sy=SY,C=C,epsilon=epsilon))
return(list(svm=SVM,mpar=mpar,cx=cx,sx=sx,cy=0,sy=0,C=C,epsilon=epsilon))
}

#------------------------------------------------------------------------------------------
# Create and fit a MLP
# current version: performs NetRuns trainings and selects the MLP 
#                  with lowest minimum criterion (SSE+decay?) 
# HN - number of hidden nodes 
# decay - decay value (in 0 ... 1)
# NR - NetRuns - number of MLP trainings
# maxit   - maximum number of training epochs
# type = 1 mlp, 2 mlpe, 3 mlp neuralnet (rprop)
mlp.fit<- function(x,data,HN=2,decay=0.0,NR=3,maxit=100,task,scale,type=1,metric="SSE")
{
 if(type==2) { MinNet<-vector("list",NR) } # ensemble

#scale="inputs"
#print(summary(data)); NR=1;
#print(paste(" > mlp HN:",HN,"decay:",decay,"NR:",NR,"maxit:",maxit,"task:",task,"scale:",scale))
 outindex=output_index(x,data)
 if(scale=="inputs" || scale=="all") # scale the inputs 
   { S<-scaleinputs(data,outindex)
     data=S$data; cx=S$cx; sx=S$sx
   }
 else{cx=NULL; sx=NULL}
 if(substr(task,1,3)=="reg" && scale=="all") # scale the output
   { CY=mean(data[,outindex])
     SY=sd(data[,outindex])
     data[,outindex]<-xtransform(data[,outindex],"scale",A=CY,B=SY)
   }
 else{CY=0; SY=0}

 MinError<-Inf # maximum real number
 if(substr(task,1,3)=="reg") LINOUT<-TRUE
 else LINOUT=FALSE

# XXX 
#LINOUT=FALSE
 
 goon=TRUE;i=1;
 while(goon)#for(i in 1:NR)
    {
# XXX
      if(type==1)
      { if(HN>0) Net<-nnet(x,data=data,size=HN,decay=decay,trace=FALSE,skip=FALSE,maxit=maxit,linout=LINOUT)
        else     Net<-nnet(x,data=data,size=0,decay=decay,trace=FALSE,skip=TRUE,maxit=maxit,linout=LINOUT)# 0 hidden nodes
        err<-Net$value #---- this code used minimum penalized error (SSE+decay...) #err<-Sse(data[,outindex],Net$fitted.values)
      }
      else if(type==2)
      { if(HN>0) MinNet[[i]]<-nnet(x,data=data,size=HN,decay=decay,trace=FALSE,skip=FALSE,maxit=maxit,linout=LINOUT)
        else     MinNet[[i]]<-nnet(x,data=data,size=0,decay=decay,trace=FALSE,skip=TRUE,maxit=maxit,linout=LINOUT)# 0 hidden nodes
      }
      #else # type==3
      #{
      # #if(metric=="SSE") ERR="sse" 
      # #else ERR="sad" # sse2 #sad
      # NM=names(data); ALL=setdiff(1:ncol(data),outindex) 
      # fx=as.formula(paste(NM[outindex], paste(NM[ALL],collapse='+'), sep='~'))
      # Net=neuralnet(fx,data=data,hidden=HN,linear.output=LINOUT,rep=1) #,err.fct=ERR)
      # #cat(data[1,outindex],Net$response[1],"\n")
      # err<-Sse(data[,outindex],Net$response)
      #}
#print(paste("HN:",HN,"decay:",decay,"Err:",err,"Min:",MinError,sep=" "))
      if(type!=2 && MinError>err){ MinError<-err; MinNet<-Net;}
#cat("metric:",metric,"HN:",HN,"Err:",err,"Min:",MinError,"\n")
      i=i+1
      if(i>NR || MinError==0) goon=FALSE
    }
NNmodel<-list(mlp=MinNet,mpar=c(HN,decay,NR,maxit),cx=cx,sx=sx,cy=CY,sy=SY)
#print(" > end mlpfit <")
return (NNmodel)
}

#---------------------------------------------------------------------------------
# Performs the Data Mining stage.
# The dataset (x - inputs, y - output) is used to fit the model.
# It returns the performance of a model using all data, holdout or kfold.
# PARAMETERS:
# x - a formula (default) or input dataframe, matrix, vector or factor
# data - NULL or the dataframe (default) or the output vector or factor 
# Runs - number of runs (e.g. 1, 10, 20, 30, 100, 200)
# method - vector with 2 arguments (see point A) or list (point B)
#      A) 1st: estimation method: "all", "holdout", "holdoutorder", "holdoutinc", "kfold", "kfoldo"
#          - holdoutorder is similar to holdout, except that a sequencial ordered split is used
#          - holdoutfor is equal to holdoutorder but is used for time series forecasting, special holdout that does not use out-samples
#            instead of the random split. This is quite useful for temporal data.
#          - holdoutinc - incremental retraining method (vpar=batchsize)
#         2nd: the estimation method parameter (e.g. 2/3, 3/4 - holdout or 3, 5 or 10 - kfold)
# model - the type of model used. Currently: "lm", "svm","mlp","logistic","naivebayes", "knn", "dt","randomforest"
# task - "default", "class" (or "classification", "cla"), "reg" (or "regression")
# search - NULL or a vector of seaching hyperparameters (used for "mlp", "svm", "knn")
# mpar - vector with the internal searching parameters (for "mlp", "svm", "knn"):
# mpar[1] - knn: internal estimation method ("all", "kfold", "holdout", "holdoutorder") 
# mpar[2] - knn: internal split ratio or number of kfolds
# mpar[1] - svm or mlp: 1st internal parameter (svm - C, mlp - NR)
# mpar[2] - svm or mlp: 2nd internal parameter (svm - Epsilon, mlp - maxit) 
# mpar[3] - svm or mlp: internal estimation method ("all", "kfold", "holdout", "holdoutorder")
# mpar[4] - svm or mlp: internal split ratio or number of kfolds 
# scale   - "none", "inputs", "all" - if a 0 mean and 1 std scale is used
#           "default" - the best generic method is automatically applied
# transform - if "log" and task="reg" then a logistic transformation is used in the output  
#             if "logpositive" is equal to "log" except that negative output values are set to zero.
# feature - "default" - "none" or "sens" if model="mlp","svm",...
#           "sens" - sensitivity needs to be stored for each fit... ("sensg" - sensitivity with gradient)
#           "simpv" "simpg" - similar to above but stores more info (sresponses) (v -variance, g- gradient)
#           "none" - if no feature selection is used
#            when using a feature selection the following scheme should be adopted:
#            -  feature=c(FM,Runs,Method,MethodPar,Search), where
#                  FM - "sabs", "sfs", "sffs", "sbs", "sbfs"
#                       "sabsv" - sabs with variance, "sabsg" - sabs with gradient
#                  Runs/Fixed - number of feature selection runs OR "fN", where N is the number of fstop deletions
#                  Method - "holdout" or "kfold"
#                  MethodPar - splitratio or kfold (vpar)
#                  Search - primary search parameter (mlp, svm, knn)
#---------------------------------------------------------------------------------
mining=function(x,data=NULL,Runs=1,method=NULL,model="default",task="default",search="heuristic",mpar=NULL,feature="none",scale="default",transform="none",debug=FALSE,...)
{
 if(is.null(data)) 
   { data<-model.frame(x) # only formula is used 
     outindex=output_index(x,data)
     x<-as.formula(paste(names(data)[outindex]," ~ .",sep=""))
   }
 #else { if(feature[1]=="aries" && length(feature)>1) outindex=as.numeric(feature[2]) 
 else outindex=output_index(x,data) 

 firstm <- proc.time()
 previoustm<-firstm

 time<-vector(length=Runs)
 error<-vector(length=Runs)
 PRED<-vector("list",Runs); TEST<-vector("list",Runs)

 vpar=2/3; # default 
 if(is.null(method)) {vmethod="holdout"; vpar=2/3;}
 else if(length(method)>1) 
 { vpar=as.numeric(method[2]); vmethod=method[1];}
 else {vmethod=method[1]}

 if(vmethod=="all") vpar=1
 else if(vmethod=="kfold" && vpar<1) vpar<-10 # default
 else if(substr(vmethod,1,7)=="holdout" && vpar>=1 && !is.null(data)) #all holdout types 
   {  
      NR<-nrow(data)
      if(vmethod=="holdoutinc") Runs=(ceiling(NR/vpar)-1)
   }

 feature=defaultfeature(feature)
 task=defaultask(task,model,data[1,outindex])

 npar=modelnpar(model)
 if(npar>0) { if(vmethod=="kfold") par=matrix(nrow=Runs*vpar,ncol=npar) else par=matrix(nrow=Runs,ncol=npar) }
 else par=NULL

 metric=getmetric(mpar,model,task)
 #cat(" >> metric:",metric,"\n")

#cat(" >> feature:",feature,"\n")
 if(substr(feature[1],1,4)=="sens" || substr(feature[1],1,4)=="sabs" || substr(feature[1],1,4)=="simp") 
    { 
       if(vmethod=="kfold") MULT=vpar
       else MULT=1
       SEN<-matrix(ncol=NCOL(data),nrow=(Runs*MULT) )
       if(substr(feature[1],1,4)!="sens") { SRESP=TRUE;FSRESP=TRUE;}
       else {SRESP=NULL; FSRESP=FALSE;}
       imethod=switch(feature[1],sabsv=,simpv="sensv",sabs=,simp=,simpg="sensg",sabsr=,simpr="sensr",feature[1])
    }
 else {SEN<-NULL; SRESP=NULL;FSRESP=FALSE;imethod="none"}

 #cat("is.null sen?",is.null(SEN),"\n")
 if(feature[1]!="none" && substr(feature[1],1,4)!="sens") # store features?
  { if(vmethod=="kfold") attrib=vector("list",Runs*vpar)
    else attrib=vector("list",Runs)
  }
 else attrib=NULL 

 if(substr(vmethod,1,6)=="kfoldo") order=TRUE
 else order=FALSE

 for(i in 1:Runs)
 {
  #cat("i:",i,"of ",Runs,"\n")
  if(debug) { curtm <- proc.time(); print(paste("Run: ", i, " time:",(curtm-firstm)[3])) }

  if(vmethod=="all")
     { L<-fit(x=x,data=data,model=model,task=task,search=search,mpar=mpar,scale=scale,transform=transform,feature=feature,...)
       P<-predict(L,data)
       TS<-data[,outindex]
       if(!is.null(SEN)) IMPORTANCE<-Importance(L,data,method=imethod,responses=FSRESP) # store sen
     }
  else if(substr(vmethod,1,7)=="holdout") #all holdout types 
     { 
       if(vmethod=="holdoutorder"|| vmethod=="holdoutfor") H<-holdout(data[,outindex],ratio=vpar,mode="order")
       else if(vmethod=="holdoutinc") H=holdout(data[,outindex],ratio=vpar,mode="incremental",iter=i)
       else H<-holdout(data[,outindex],ratio=vpar)
       L<-fit(x=x,data=data[H$tr,],model=model,task=task,search=search,mpar=mpar,scale=scale,transform=transform,feature=feature,...)
       P<-predict(L,data[H$ts,])
     ##if(method!="holdoutfor") P<-predict(L,data[H$ts,])
     ##else P=lforecast(L,data,start=(1+nrow(data)-vpar),horizon=vpar) 
       TS<-data[H$ts,outindex]
       if(!is.null(SEN)) IMPORTANCE<-Importance(L,data[H$tr,],method=imethod,responses=FSRESP) # store sen
     }
  else if(substr(vmethod,1,5)=="kfold")
     { L=crossvaldata(x,data,fit,predict,ngroup=vpar,order=order,model=model,task=task,feature=feature,
                      scale=scale,search=search,mpar=mpar,transform=transform,...)
       P<-L$cv.fit
       TS<-data[,outindex]
     }

  if(!is.null(SEN))
     { 
       if(vmethod=="kfold")
         {  Begi<-(i-1)*MULT+1; Endi=Begi+vpar-1;
            SEN[Begi:Endi,]<-L$sen # store sen
            if(!is.null(SRESP)) 
                 { 
                  if(i==1) SRESP=L$sresponses # store sen
                  else{ for(j in 1:length(SRESP))
                           {  
                            if(!is.null(L$sresponses[[j]]) ) # && !is.null(IMPORTANCE$sresponses[[j]])$x) # $x)) 
                              { if(is.null(SRESP[[j]])) SRESP[[j]]=L$sresponses[[j]]
                                else{ SRESP[[j]]$x=c(SRESP[[j]]$x,L$sresponses[[j]]$x);
                                      if(task=="prob") SRESP[[j]]$y=rbind(SRESP[[j]]$y,L$sresponses[[j]]$y,deparse.level=0)
                                      else if(task=="class") SRESP[[j]]$y=addfactor(SRESP[[j]]$y,L$sresponses[[j]]$y)
                                      else SRESP[[j]]$y=c(SRESP[[j]]$y,L$sresponses[[j]]$y)
                                    }
                              }
                           }
                      }
                 }
         }
       else
         { 
           SEN[i,]<-IMPORTANCE$imp # store sen
           if(!is.null(SRESP)) 
            { if(i==1) SRESP=IMPORTANCE$sresponses # store sen
              else{ for(j in 1:length(SRESP)) 
                       { 
                         if(!is.null(IMPORTANCE$sresponses[[j]]) ) # && !is.null(IMPORTANCE$sresponses[[j]])$x) # $x)) 
                           { if(is.null(SRESP[[j]])) SRESP[[j]]=IMPORTANCE$sresponses[[j]]
                             else{ SRESP[[j]]$x=c(SRESP[[j]]$x,IMPORTANCE$sresponses[[j]]$x);
                                   if(task=="prob") SRESP[[j]]$y=rbind(SRESP[[j]]$y,IMPORTANCE$sresponses[[j]]$y,deparse.level=0)
                                   else if(task=="class") SRESP[[j]]$y=addfactor(SRESP[[j]]$y,IMPORTANCE$sresponses[[j]]$y)
                                   else SRESP[[j]]$y=c(SRESP[[j]]$y,IMPORTANCE$sresponses[[j]]$y)
                                 }
                           }
                       }
                  }
            }
         }
     }

  if(!is.null(attrib))
  { if(vmethod=="kfold") for(j in 1:vpar) attrib[[ (i-1)*vpar+j ]]<-L$attributes[[j]]
    else attrib[[i]]<-L@attributes
  }

  if(npar>0) { if(vmethod=="kfold"){ ini=(i-1)*vpar+1;end=ini+vpar-1;par[ini:end,]=L$mpar} else par[i,]=as.numeric(L@mpar[1,]) }

  # store P and TS
  PRED[[i]]<-P
  TEST[[i]]<-TS
#print(PRED[[i]])
#print(TEST[[i]])
#cat(" i:",i,"metric:",metric,"\n")
  error[i]=Error(TS,P,metric)
#cat(" e:",error[i],"\n")
  curtm<- proc.time()
  time[i]<-(curtm-previoustm)[3] # execution time for run i
  previoustm<-curtm
 }
 if(!is.null(SRESP)) { for(i in 1:length(SRESP)) if( !is.null(SRESP[[i]]) && is.factor(data[,i])) SRESP[[i]]$x=factor(SRESP[[i]]$x) }
 if(npar>0) { par=data.frame(par); names(par)=modelnamespar(model);}
 return(list(time=time,test=TEST,pred=PRED,error=error,mpar=par,model=model,task=task,method=c(vmethod,vpar),sen=SEN,sresponses=SRESP,runs=Runs,attributes=attrib,feature=feature))
}

medianminingpar=function(M,parorder=NULL)
{
 NP=switch(M$model,randomforest=,knn=1,mlp=,mlpe=4,svm=3,0) 
 if(NP>0)
 {
  MPAR=M$mpar
  if(is.null(parorder)) parorder=1:NP
  i=1;stop=FALSE
  while(!stop)
  {
   val=medianfirst(MPAR[,parorder[i]])$val
   ind=which(MPAR[,parorder[i]]==val);Lind=length(ind)
   if(Lind>0) MPAR=MPAR[ind,]
   if(Lind==1) stop=TRUE
   else { i=i+1
          if(i>NP) { stop=TRUE
                     if(is.vector(MPAR)) MPAR=MPAR[1] else if(Lind>1) MPAR=MPAR[1,]
                   }
        }
  }
  return (as.numeric(MPAR))
 }
 else return (NULL)
}


# 
getmetric=function(mpar,model,task)
{
 L=length(mpar)
 if(L==0) {metric=switch(task,prob="AUC",class="ACC","MAD")}
 else if(L>2 && (model=="knn" || model=="randomforest") ) metric=mpar[3]
 else if(L>4 && (substr(model,1,3)=="mlp" || model=="svm")) metric=mpar[5]
 else if(L>0 && (substr(model,1,3)!="mlp" && model!="svm" && model!="knn") ) metric=mpar[1]
 else metric=switch(task,prob="AUC",class="ACC","MAD")
 return (metric)
}

#-- example code of variant to allow direct use of models built with nnet, ksvm, randomForest, etc...
#-- needs more coding... still in beta phase
#-- type should be explicitly defined, as it depends on the predict function.
#-- for instance: nnet - type="raw" or "class", kvsm - type ="response", "votes, etc...

# method="randomforest" - Leo Breiman method; "sens","sensv", "sensg", "sensv", "sensi" - Embrechts method ; "data" - uses all data.frame
#        "sensi" - Embrechts method extended to i-D sensitivity (with interactions) 
# sampling - "regular" - grid spread from min to max or "quantile" - based on the distribution percentiles
# responses - if TRUE then the full response (y) values are stored
# interactions - default NULL else is a vector of one or two variables to explore 2 variable interactions (VECC and VEC3)...
# baseline - only used by "sens" method: "mean" (pure kewley method), "median" or row of base values

# To be used only if model is not of "model" type:
# outindex - the output column index 
# task - NULL or "class" or "prob" or "reg" 
# PRED - new prediction function, needs to return a vector for "reg" or matrix of probabilities for "prob"

# to do: think about I-D sensitivity, where I is the number of features? 
# measure - "variance", "range", "gradient"
Importance=function(M,data,RealL=6,method="sens",measure="gradient",sampling="regular",baseline="mean",responses=TRUE,
                    outindex=NULL,task="default",PRED=NULL,interactions=NULL)
{
#model=M;data=d;RealL=6;method="sens";measure="variance";sampling="regular";responses=TRUE;outindex=NULL;task=NULL;

 #model<<-model;DDD<<-data;RealL<<-RealL;method<<-method;measure<<-measure;sampling<<-sampling;responses<<-responses; 
 if(class(M)=="model"){outindex=M@outindex; task=M@task; Attributes=M@attributes}
 else # another R supervised learning model
 { 
  task=defaultask(task,model="default",data[1,outindex])
  if(is.null(PRED)) 
  { nclass=class(M);nclass=nclass[length(nclass)] 
    if(task=="prob")
     {PRED=switch(nclass,
                       lda=,qda=function(M,D){predict(M,D)$posterior},
                       randomForest=function(M,D){predict(M,D,type="prob")},
                       # and so on, to be completed...
                       NULL)
     }
    else if(task=="class")
     {PRED=switch(nclass,
                       lda=,qda=function(M,D){predict(M,D)$class},
                       randomForest=function(M,D){predict(M,D,type="response")},
                       # and so on, to be completed...
                       NULL)
     }
    else{ PRED=switch(nclass,
                       lm==function(M,D){predict(M,D)},
                       randomForest=function(M,D){predict(M,D)},
                       # and so on, to be completed...
                       NULL)}
  }
  Attributes=1:ncol(data)
 }

 # use the method proposed by Leo Breiman 2001:
 if(method=="randomforest") # user should be sure this is a randomforest!
 { if(class(M)=="model") Sv=randomForest::importance(M@object,type=1)
   else Sv=randomForest::importance(M,type=1) # true randomforest
   ASv<-abs(Sv)
   imp=100*abs(ASv)/sum(ASv) 
   RESP=NULL
 }
 else if(substring(method,1,4)=="sens" || method=="data") # set SPREAD
 {
  measure=switch(method,sensv="variance",sens=,sensg="gradient",sensr="range",measure)
  if(method=="sensi") INTERACT=TRUE else INTERACT=FALSE
  method=switch(method,sensi=,sensv=,sensg=,sensr="sens",method)

  Dsize=NCOL(data); SPREAD=vector("list",Dsize); MR=RealL;
  if(method=="data") Lsize=NROW(data)
  IRANGE=setdiff(Attributes,outindex)
  for(i in IRANGE) # set the SPREAD 
  { 

    if(is.ordered(data[1,i])) { AUX=levels(data[1,i]);LNAUX=length(AUX)
                                if(LNAUX<RealL) SPREAD[[i]]=AUX
                                else if(sampling=="regular") { NAUX=1:LNAUX; SPREAD[[i]]=AUX[unique(round(seq(NAUX[1],NAUX[LNAUX],length.out=RealL)))] }
                                else if(sampling=="quantile")SPREAD[[i]]=AUX[unique(quantile(as.numeric(data[,i]),seq(0,1,length.out=RealL)))]
                              }
    else if(is.factor(data[1,i])) SPREAD[[i]]=levels(data[1,i])
    else # numeric 
       { if(sampling=="regular"){ SPREAD[[i]]=seq(min(data[,i]),max(data[,i]),length.out=RealL);}
         else if(sampling=="quantile"){ SPREAD[[i]]=unique(quantile(data[,i],seq(0,1,length.out=RealL)));attr(SPREAD[[i]],"names")=NULL;}
       }
    MR=max(MR,length(SPREAD[[i]]))
  } #end of for cycle, SPREAD SET!

 if(!is.null(interactions))
 {
  LINT=length(interactions)
  i1=interactions[1]; L=length(SPREAD[[i1]]);
  if(LINT==2) {Dsize=1;JDOMAIN=interactions[2];}
  else if(LINT>2){Dsize=1; INTERACT=TRUE;}
  else JDOMAIN=setdiff(IRANGE,i1)
 }
 else LINT=0
 if(method=="sens")
 { 
   if(class(baseline)=="data.frame") v=baseline
   else #--- set the average/median input
   { v=data[1,]
    for(i in IRANGE) 
    { 
     if(is.factor(data[1,i]))  
       { if(is.ordered(data[1,i])) ModeClass=middleclass(data[,i],method=baseline) # "mean" or "median"
         else ModeClass=mostcommon(data[,i])
         v[1,i]=levels(data[1,i])[ModeClass]
       }
     else if(baseline=="mean") v[1,i]=mean(data[,i])
     else if(baseline=="median") v[1,i]=median(data[,i])
    } #end of for cycle 
   }
   # new data frame with MR rows, for speeding up the Importance computation:
   data=v[1,]
   if(LINT>0) 
   { L2=2;
     if(!INTERACT) {for(j in JDOMAIN) L2=max(L2,length(SPREAD[[j]]));MR=L*L2}
     else{ MR=1;for(j in 1:LINT) MR=MR*length(SPREAD[[interactions[j]]])
         }
   }
   data=v[rep(1,MR),]
#cat("MR:",MR,"\n")
 }
 if(is.null(interactions)) # method=="sens" || method=="sensgrad" || method=="data"
 {
  Sv=rep(0,Dsize)
  if(responses) RESP=vector("list",length=Dsize) else RESP=NULL
  for(i in IRANGE)
   { 
    L=length(SPREAD[[i]]);MEM=data[,i]
    if(responses) { xinp=SPREAD[[i]];xname=names(data)[i];}
    if(method=="data") 
      { Y=NULL 
        for(k in 1:L)
        {
         if(is.factor(data[1,i])) data[,i]=factor(rep(SPREAD[[i]][k],Lsize),levels=SPREAD[[i]])
         else data[,i]=rep(SPREAD[[i]][k],Lsize)
         if(class(M)=="model") y=predict(M,data) else y=PRED(M,data)
         data[,i]=MEM # restore previous values
         Y=c(Y,mean(y))
        }
      }
    else{ data[(1:L),i]=SPREAD[[i]]
          if(class(M)=="model") Y=predict(M,data[(1:L),]) else Y=PRED(M,data[(1:L),])
          data[,i]=MEM # restore previous values
        }
    Sv[i]=switch(measure,variance=variance_responses(Y),gradient=gradient_responses(Y),range=range_responses(Y))
    if(responses) RESP[[i]]=list(n=xname,l=L,x=xinp,y=Y)
   }
  Sum=sum(Sv)
  if(Sum==0) # no change in the model, equal importances to all attributes
  {imp=rep(1/(Dsize-1),length=Dsize);imp[outindex]=0;}
  else imp=Sv/Sum
  } #end of if null interactions
  else if(!INTERACT) # LINT<3) # if(!is.null(interactions))
  {
   Sv=rep(0,Dsize)
   if(responses) RESP=vector("list",length=Dsize) else RESP=NULL
   if(method=="data") MR=NROW(data) 
   for(j in JDOMAIN)
    {
     if(method=="sens")
     { MEM=data[,j] 
       L2=length(SPREAD[[j]]); MR=L*L2;
       for(k in 1:L2) {ini=1+(k-1)*L;end=ini+L-1;data[ini:end,i1]=SPREAD[[i1]];data[ini:end,j]=SPREAD[[j]][k];}
       if(responses) {xinp=data[1:MR,c(i1,j)];xname=c(names(data)[i1],names(data)[j]);}
       if(class(M)=="model") Y=predict(M,data[(1:MR),]) else Y=PRED(M,data[(1:MR),])
       data[,j]=MEM # restore previous values
     }
     else if(method=="data") 
     { L2=length(SPREAD[[j]]);Y=NULL
       MEM=data[,c(i1,j)];
       if(responses) {xname=c(names(data)[i1],names(data)[j]);xinp=data[1:(L*L2),c(i1,j)];o=1;}
       for(k in 1:L2) 
       for(l in 1:L) 
       { 
         if(responses){xinp[o,1]=SPREAD[[i1]][l];xinp[o,2]=SPREAD[[j]][k];o=o+1;}
         if(is.factor(data[1,i1])) data[,i1]=factor(SPREAD[[i1]][l],levels=levels(data[1,i1])) else data[,i1]=SPREAD[[i1]][l]
         if(is.factor(data[1,j])) data[,j]=factor(SPREAD[[j]][k],levels=levels(data[1,j])) else data[,j]=SPREAD[[j]][k]
         if(class(M)=="model") y=predict(M,data) else y=PRED(M,data)
         Y=c(Y,mean(y)) 
         #data[,i1]=MEM[,1];data[,j]=MEM[,2] # restore previous values
       }
       data[,i1]=MEM[,1];data[,j]=MEM[,2] # restore previous values
     }
     if(length(interactions)==2) k=1 else k=j
     Sv[k]=switch(measure,variance=variance_responses(Y),gradient=gradient_responses(Y),range=range_responses(Y))
     if(responses) RESP[[k]]=list(n=xname,l=L,x=xinp,y=Y) 
    }
#cat("SV:",Sv,"\n")
    Sum=sum(Sv)
    if(Sum==0) # no change in the model, equal importances to all attributes
    {imp=rep(1/(Dsize-1),length=Dsize);imp[outindex]=0;}
    else imp=Sv/Sum
  }
  else # LINT>=3, INTERACT
  {
   Sv=rep(0,1)
   if(responses) {RESP=vector("list",length=1);L=vector(length=LINT);} else RESP=NULL
   # interactions / SPREAD / data
   E=1;
#cat("LINT:",LINT,"MR:",MR,"\n")
   for(i in 1:LINT)
   { k=interactions[[i]];LS=length(SPREAD[[k]]);T=MR/(E*LS)
#cat("i:",i,"LS:",LS,"T:",T,"\n")
     if(is.factor(data[,k][1])) data[,k]=rep(factor(SPREAD[[k]],levels=SPREAD[[k]]),each=E,times=T) else data[,k]=rep(SPREAD[[k]],each=E,times=T)
     E=E*LS;
     if(responses) L[i]=LS
   }
   if(class(M)=="model") y=predict(M,data) else y=PRED(M,data)
   Sv=switch(measure,variance=variance_responses(y),gradient=gradient_responses(y),range=range_responses(y))
   if(responses) RESP[[1]]=list(n=names(data)[interactions],l=L,x=data[,interactions],y=y) 
   imp=1
  }
 }
 return (list(value=Sv,imp=imp,sresponses=RESP))
}

# -- 
avg_imp_1D=function(I,measure="variance")
{
#II<<-I
 x=I$sresponses[[1]]$x 
 NC=NCOL(x)
 
 R=list(responses=vector("list",length=NC),value=vector(length=NC),imp=rep(1/NC,NC))
 for(i in 1:ncol(x))
 {
  I1=avg_imp(I,i,measure=measure)
  R$sresponses[[i]]$x=I1$sresponses[[1]]$x
  y=I1$sresponses[[1]]$y
  R$sresponses[[i]]$y=y
#print(y)
  R$value[i]=switch(measure,variance=variance_responses(y),gradient=gradient_responses(y),range=range_responses(y))
 }
#print("RVALUE:")
#print(R$value)
 SUM=sum(R$value)
 if(SUM>0) R$imp=R$value/SUM 
 return(R)
}

avg_imp=function(I,AT,measure="variance")
{
 x=I$sresponses[[1]]$x
 y=I$sresponses[[1]]$y
 CY=NCOL(y)
 X1=unique(x[,AT[1]])

 if(length(AT)>1)
 { 
   X2=unique(x[,AT[2]])
   LX=length(X2)*length(X1)
   if(is.matrix(y)) my=matrix(ncol=NCOL(y),nrow=LX)
   else if(is.factor(y)) my=factor(rep(levels(y)[1],LX),levels=levels(y))
   else my=vector(length=LX)
   Im=vector(length=LX)
   k=1;
   for(i in X1)
     for(j in X2)
     {
      W=which(x[,AT[1]]==i & x[,AT[2]]==j)
      Im[k]=W[1]
      if(is.matrix(y)) my[k,]=mean_resp(y[W,])
      else my[k]=mean_resp(y[W])
      k=k+1
     }
 }
 else
 { LX=length(X1)
   if(is.matrix(y)) my=matrix(ncol=NCOL(y),nrow=LX)
   else if(is.factor(y)) my=factor(rep(levels(y)[1],LX),levels=levels(y))
   else my=vector(length=LX)
   Im=vector(length=LX); k=1;
   for(i in X1)
     {
      W=which(x[,AT[1]]==i)
      Im[k]=W[1]
      if(is.matrix(y)) my[k,]=mean_resp(y[W,])
      else my[k]=mean_resp(y[W])
      k=k+1
     }
 }
 I$sresponses[[1]]$x=x[Im,AT]
 I$sresponses[[1]]$y=my
 I$sresponses[[1]]$n=I$sresponses[[1]]$n[AT]
 I$sresponses[[1]]$l=I$sresponses[[1]]$l[AT]
 I$value=switch(measure,variance=variance_responses(my),gradient=gradient_responses(my),range=range_responses(my))
 return(I)
}

#----------------------------------------------------------------------------------------------------
#-- Auxiliary private functions, in principle you should not need to use these:
#----------------------------------------------------------------------------------------------------
mean_resp=function(y)
{
 if(is.matrix(y)) return (colMeans(y))
 else if(is.ordered(y)) return (levels(y)[middleclass(y,method="mean")])
 else if(is.factor(y)) return (levels(y)[mostcommon(y)])
 else return (mean(y))
}


range_responses=function(y)
{
  if(is.ordered(y)) return(range_responses(as.numeric(y)))
  else if(is.factor(y)) return(range_responses(one_of_c(y)))
  else
  { 
    LV=NCOL(y) 
    if(LV==1) return (abs(diff(range(y))))
    else { res=0; for(i in 1:LV) res=res+abs(diff(range(y[,i]))); return(res/LV) }
  }
}

#-- compute the gradient of response y: vector, matrix of numeric data (probabilities) or factor
#-- ordered is only used if y is factor, true if it is an ordered factor or else false
# XXX think here!
gradient_responses=function(y,ABS=TRUE)
{
 if(is.ordered(y)) return(gradient_responses(as.numeric(y)))
 else if(is.factor(y)) return(gradient_responses(one_of_c(y)))
 else{ if(ABS) FUNC=abs else FUNC=identity
       G=mean(FUNC(diff(y)))
     }
 return (G)
}
 
# compute the variance of response y (vector or matrix or data.frame)
variance_responses=function(y)
{ if(is.ordered(y)) return(variance_responses(as.numeric(y)))
  else if(is.factor(y)) return(variance_responses(one_of_c(y)))
  else
  { LV=NCOL(y)
    if(LV==1) return (var(y))
    else if(LV==2) return (var(y[,1]))
    else # LV>2 
    { res=0; for(i in 1:LV) res=res+var(y[,i]); return(res/LV) }
  }
}

# --- auxiliary functions, do not use directly:
# x is a $sresponses object
resp_to_list=function(x,TC=-1)
{
 LX=length(x$x);lines=LX/x$l;
 RES=vector("list",lines)
 for(i in 1:lines)
    { 
      ini=(i-1)*x$l+1;end=ini+x$l-1
      if(TC==-1) M=cbind(x$x[ini:end],x$y[ini:end])
      else M=cbind(x$x[ini:end],x$y[ini:end,TC])
      RES[[i]]=M   
    }
 return (RES)
}
# --- end of auxiliary functions --------------
