isbest=function(Cur,Best,metric)
{if(is.na(Cur)) return (FALSE) 
 else{ return (switch(metric,
                ALIFTATPERC=,NALIFT=,ALIFT=,ACC=,KAPPA=,COR=,TPR=,TNR=,R2=,R22=,NAREC=,NAUC=,TPRATFPR=,AUC=Cur>Best,
                Cur<Best))}
} 
worst=function(metric)
{ return (switch(metric, 
                 ALIFTATPERC=,NALIFT=,ALIFT=,ACC=,KAPPA=,COR=,TPR=,TNR=,R2=,R22=,NAREC=,NAUC=,TPRATFPR=,TOLERANCE=,AUC=-Inf,
                 Inf))} 

# --- Error function:
# if you want to choose a different internal model selection function, change here!
Error=function(y,x,metric,D=0.5,TC=-1,val=NULL)
{
 if(class(metric)=="function") return(metric(y,x))
 else 
 { if(TC==-1 && (metric=="ALIFTATPERC" || metric=="ALIFT" || metric=="NALIFT" || metric=="NAUC"||metric=="TPRATFPR"||metric=="TPR"||metric=="TNR") ) TC=length(levels(y[1])) # set TC to the last class index
   return (switch(metric,
                AUC=ROCcurve(y,x,TC=TC)$auc, # classification G
                ACC=Accuracy(y,x,D,TC=TC), # classification G
                CONF=Conf(y,x,D,TC=TC), # classification G
                ALIFT=LIFTcurve(y,x,TC=TC)$area, # classification
                TPR=metrics(y,x,D,AUC=FALSE,BRIER=FALSE)$tpr[TC], # classification
                TNR=metrics(y,x,D,AUC=FALSE,BRIER=FALSE)$tnr[TC], # classification
                KAPPA=Kappa(y,x,D,TC=TC), # classification G
                SAD=Sad(y,x), 
                MAD=,MAE=Mad(y,x),MdAE=Mdae(y,x),GMAE=,GMAD=Gmad(y,x),
                RMAD=,RAE=Rae(y,x),
                SSE=Sse(y,x),
                MSE=Mse(y,x),MdSE=Mdse(y,x),GMSE=Gmse(y,x),
                RSE=Rse(y,x),
                RMSE=Rmse(y,x),
                RRSE=Rrse(y,x),
                ME=Me(y,x),
                SMinkowski3=SMinkowski(y,x,3), # experimental
                MMinkowski3=MMinkowski(y,x,3), # experimental
                MdMinkowski3=MdMinkowski(y,x,3), # experimental
                COR=Correlation(y,x),
                R2=R2(y,x),
                R22=R22(y,x),
                BRIER=Tbrier(y,x,TC=TC), # classification G
                NAREC=Narec(y,x,val=val),
                TOLERANCE=Tolerance(y,x,val=val),
                NAUC=Nauc(y,x,TC=TC,val=val), # classification
                NALIFT=Nalift(y,x,TC=TC,val=val), # classification
                TPRATFPR=Tprfpr(y,x,val=val,TC=TC), # classification
                ALIFTATPERC=Aliftperc(y,x,val=val,TC=TC), # classification
                MdAPE=Mdape(y,x),RMSPE=Rmspe(y,x),RMdSPE=Rmdspe(y,x),
                MAPE=Mape(y,x), SMAPE=Smape(y,x), SMdAPE=Smdape(y,x), #SMAPE2=Smape2(y,x),
		MRAE=Mrae(y,x,ts=val),MdRAE=Mdrae(y,x,ts=val),GMRAE=Gmrae(y,x,ts=val),THEILSU2=TheilsU2(y,x,ts=val),MASE=Mase(y,x,ts=val),
                NA # not defined...
         ))
 }
}

#---------------------------------------------------------
# Brier Score: computes the brier score (SSE) for probabilitic output models ([0,1]) 
# y - target factor
# x - numeric predictions (in 0 to 1 probabilities) 
#---------------------------------------------------------
Brier=function(y,x,TC=-1)
{
  L=levels(y[1]);N=length(y); 
  if(TC>0){ NL=1;x=x[,TC];L=L[TC];MSE=1;} else{NL=length(L); PROP=table(y)[]/N; TBRIER=0;MSE=rep(FALSE,(NL+1));} 
  for(i in 1:NL) 
  {
      T=as.numeric(y==L[i])
      if(is.vector(x)) MSE[i]=Sse(T,x)/N 
      else { MSE[i]=Sse(T,x[,i])/N # Brier=MSE
             if(PROP[i]>0) TBRIER=TBRIER+PROP[i]*MSE[i]
           }
  }
  if(!is.vector(x)) MSE[(NL+1)]=TBRIER
  return (MSE)
}
Tbrier=function(y,x,TC=-1){ B=Brier(y,x,TC=TC);return(B[length(B)])}

#----------------------------------------
# RECurve by Paulo Cortez, 2006@
#
# following the article of Bi & Bennett 2003:
# J. Bi and K. Bennett, Regression Error Characteristic curve
# In Proceedings of 20th Int. Conf. on Machine Learning (ICML),  
# Washington DC, USA, 2003.
#
# vector.error - vector with the residuals or errors
#                vector.error = y (desired) - x (predicted)
#              or vector.error= y and x= predictions (2nd mode)
RECcurve<-function(vector.error,x=NULL)
{
#print(vector.error)
 if(!is.null(x)) vector.error=(vector.error-x)
 Correct=0; Eprev=0; 
 ES<-sort(abs(vector.error))
#print(ES)
 M<-length(ES)+1; M1=M-1;
 X<-matrix(nrow=M,ncol=2)
 M<-length(ES); k=1;i=1;notstop=TRUE;
 while(notstop)
  { a=0; while( (i+a)<M && ES[(i+a+1)]==ES[(i+a)] ) a=a+1;
    if(a>0) {i=i+a-1; Correct=Correct+a-1;}
#cat(" >> i:",i,"a:",a,"k:",k,"prev:",Eprev,"ESi:",ES[i],"\n")
    if(Eprev<ES[i])
      { X[k,1]<-Eprev; X[k,2]<-Correct/M1; Eprev<-ES[i]; k=k+1;}
    Correct=Correct+1
    i=i+1;
    if(i>M1) notstop=FALSE;
  }
  X[k,1]<-ES[M]
  X[k,2]<-Correct/M1
#print(X)
#cat("M:",M,"k:",k,"Cor:",Correct,"\n")
  #X=na.omit(X) #X[,2]<-100*X[,2] # put in percentage
  return (X[(1:k),])
}

# several metric functions (for fast access)
# sum of squared errors
# x - vector of predictions, y - vector of desired values  


SMinkowski=function(y,x,q=1){ return (sum( abs(y-x)^q )) } # minkowski loss function, Bishop 2006, pattern recognition...
MMinkowski=function(y,x,q=1){ return (mean( abs(y-x)^q )) } 
MdMinkowski=function(y,x,q=1){ return (median( abs(y-x)^q )) }

Me=function(y,x) { return (sum(y-x)) }
Sse=function(y,x) { return (sum((y-x)^2)) }
Mse=function(y,x) { return (mean((y-x)^2)) }
Mdse=function(y,x) { return (median((y-x)^2)) }
Gmse=function(y,x){ return (Gmean((y-x)^2)) }
Rse=function(y,x,ymean=mean(y)) { return (100* sum((y-x)^2)/sum((y-ymean)^2)) } # relative squared error
Rmse=function(y,x) { return (sqrt(mean((y-x)^2))) }
Rrse=function(y,x) { return ( 100*sqrt(sum((y-x)^2)/sum((y-mean(y))^2)) ) }
Sad=function(y,x) { return (sum(abs(y-x))) }
Mad=function(y,x) { return (mean(abs(y-x)))}
Gmad=function(y,x){ return (Gmean(abs(y-x)))}
Rae=function(y,x,ymean=mean(y)) { return (100*sum(abs(y-x))/sum(abs(y-ymean))) } # also known as CumRAE, RelMAE
#Smape=function(y,x) { return ( 100*mean(abs(y-x)/(x+y)))} # wikipedia def.

# forecasting specific, R. Hyndman IJF 2006 "Another Look at Measures of Forecast Accuracy" definition:
Mdae=function(y,x) { return (median(abs(y-x)))}
Mape=function(y,x) { return (100*mean(abs((y-x)/y)))} # R. Hyndman IJF 2006 def:
Mdape=function(y,x) { return (100*median(abs((y-x)/y)))}
Rmspe=function(y,x) { return (sqrt(100*mean(((y-x)/y)^2)))}
Rmdspe=function(y,x) { return (sqrt(100*median(((y-x)/y)^2)))}
#Smape=function(y,x) { return (100*mean(abs(y-x)/((abs(x)+abs(y))/2)))} # def. of http://www.neural-forecasting-competition.com/motivation.htm
Smape=function(y,x) { return (200*mean(abs(y-x)/(abs(x)+abs(y))))}
Smdape=function(y,x) { return (200*median(abs(y-x)/(abs(x)+abs(y))))}

# relative errors:
# b= benchmark naive method forecasts.

# random walk, with or without drift
# ts - time series in samples
# H - number of forecasts, horizon
randomwalk=function(ts,H,drift=TRUE)
{ if(drift) drift=mean(diff(ts)) else drift=0
  return (rep(ts[length(ts)],H)+1:H*drift)
}
Gmean=function(x){return(prod(x)^(1/(length(x))))} # auxiliar geometric mean

Mrae=function(y,x,ts,b=randomwalk(ts,length(y))){if(is.null(ts)&&is.null(b)) return(NA) else return(mean(abs((y-x)/(y-b))))}
Mdrae=function(y,x,ts,b=randomwalk(ts,length(y))) { if(is.null(ts) && is.null(b)) return(NA) else return (median(abs((y-x)/(y-b)))) }
Gmrae=function(y,x,ts,b=randomwalk(ts,length(y))) { if(is.null(ts) && is.null(b)) return(NA) else return (Gmean(abs((y-x)/(y-b)))) }
TheilsU2=function(y,x,ts,b=randomwalk(ts,length(y))) { if(is.null(ts) && is.null(b)) return(NA) else return (Rmse(y,x)/Rmse(y,b)) } # theils'U or U2
#Mase=function(y,x,ts) { N=length(ts); K=1/(N-1); SUM=sum(abs(ts[2:N]-ts[1:(N-1)])); return ( mean( abs( (y-x)/(K*SUM) )) ) }
# ts - time series in samples 
Mase=function(y,x,ts) { N=length(ts); MEAN=mean(abs(ts[2:N]-ts[1:(N-1)])); return ( mean(abs((y-x)/MEAN)) ) } # faster variant?

# wikipedia:
R2=function(y,x,ymean=mean(y)) { return (1-sum((y-x)^2)/sum((y-ymean)^2)) }
R22=function(y,x,ymean=mean(y)) { return (sum((x-mean(x))^2)/sum((y-ymean)^2)) }
Ss=function(y){return (sum((y-mean(y))^2))} # sum of squares


Correlation=function(y,x) { COR=suppressWarnings(cor(x,y)); if(is.na(COR)) COR=0; return(COR)}
Narec=function(y,x,val=1){ if(is.null(val)) val=1; R=RECcurve(y,x);if(R[nrow(R),1]>val)R=partialcurve(R,val) else val=R[nrow(R),1];return(curvearea(R,val))}
Tolerance=function(y,x,val=1){ if(is.null(val)) val=1; R=RECcurve(y,x);if(R[nrow(R),1]>val)R=partialcurve(R,val);if(is.vector(R)) return (R[2]) else return(R[nrow(R),2])}
Nauc=function(y,x,val=1,TC=-1){ if(is.null(val)) val=1; RR=ROCcurve(y,x,TC=TC); RR=partialcurve(RR$roc,val); return(curvearea(RR,val))}
Nalift=function(y,x,val=1,TC=-1){ if(is.null(val)) val=1; RR=LIFTcurve(y,x,TC=TC); RR=partialcurve(RR$alift,val); return(curvearea(RR,val))}
Tprfpr=function(y,x,val,TC=-1){ if(is.null(val)) val=0.01; RR=ROCcurve(y,x,TC=TC); RR=partialcurve(RR$roc,val);if(is.vector(RR)) return (RR[2]) else return(RR[nrow(RR),2])}
Aliftperc=function(y,x,val,TC=-1){ if(is.null(val)) val=0.1; RR=LIFTcurve(y,x,TC=TC); RR=partialcurve(RR$alift,val);if(is.vector(RR)) return (RR[2]) else return(RR[nrow(RR),2])}

# x - vector of predictions, y - vector of desired values  
# MEAN - it should be the mean value of the training set, XXX 
metrics=function(y,x=NULL,D=0.5,TC=-1,AUC=TRUE,BRIER=FALSE,Run=1,task="default")
{ 
  if(is.list(y) && is.null(x)) { x=y$pred[[Run]];y=y$test[[Run]];}
  if( is.factor(y)) # || (is.list(y) && !is.null(y$task) && (y$task=="class"||y$task=="prob"))) 
  {
    conf=Conf(y,x,D=D,TC=TC)
    Total=sum(conf);Diag=0;DiagR=0;C=NCOL(conf)
    for(i in 1:C) 
    { Diag=Diag+conf[i,i]
      DiagR=DiagR+(sum(conf[i,])*(sum(conf[,i])/Total))
    }
    Kap=100*(Diag-DiagR)/(Total-DiagR)
    Acc=c(Diag/Total)*100 # total accuracy, in percentage
    LC=C; if(C<3) LC=1
    acc_class=vector(length=LC);sen_class=acc_class;spe_class=acc_class;pre_class=acc_class;
    for(k in 1:C) # class 
 	   {
       	    TP<-conf[k,k]
	    FN<-0
 	    for(i in 1:C) # iterator?
		if(i!=k) FN<-FN+conf[k,i]
            FP<-0
 	    for(i in 1:C) # iterator?
		if(i!=k) FP<-FP+conf[i,k]
            TN<-Total-TP-FN-FP
            acc_class[k]<- 100*(TP+TN)/Total 
            if(TP!=0) sen_class[k]<- 100*TP/(FN+TP) 
            else sen_class[k]=0
            if(TN!=0) spe_class[k]<- 100*TN/(TN+FP)
            else spe_class[k]=0
            if(TP!=0) pre_class[k]<- 100*TP/(TP+FP)
            else pre_class[k]=0
           }
    tbrier=NULL;brier=NULL;auc=NULL;tauc=NULL;
    if(!is.factor(x)){
         if(BRIER) { brier=Brier(y,x,TC=TC);if(TC<1){tbrier=brier[C+1];brier=brier[1:C];} else {brier=c(brier,brier);tbrier=brier[1]}}
         if(AUC){ if(C>2) {roc=multiROC(y,x);tauc=roc$auc;auc=vector(length=C); for(i in 1:C) auc[i]=roc$roc[[i]]$auc;}
             else {tauc=twoclassROC(y,x[,TC],Positive=levels(y[1])[TC])$auc;auc=c(tauc,tauc);}
                     }
        }
    return (list(conf=conf,acc=Acc,kappa=Kap,acclass=acc_class,tpr=sen_class,tnr=spe_class,precision=pre_class,tauc=tauc,auc=auc,brier=brier,tbrier=tbrier))
  }
  else return (list(me=Me(y,x),mad=Mad(y,x),sse=Sse(y,x),mse=Mse(y,x),rmse=Rmse(y,x),rae=Rae(y,x),rrse=Rrse(y,x),rse=Rse(y,x),
                    mape=Mape(y,x),smape=Smape(y,x),cor=Correlation(y,x),r2=R2(y,x)
                   ))
}

# convert matrix or data.frame into factor with major class 
majorClass=function(x,L)
{
 if(is.vector(x)) return (factor(L[which.max(x)],levels=L))
 else 
 { NX=nrow(x)
   y=vector(length=NX)
   for(i in 1:NX) y[i]=L[which.max(x[i,])]
   return (factor(y,levels=L))
 }
}

# target - vector of factor with the desired values 
# predict - vector of factor with the predicted values
# D - decision thresold
# TC - target concept class, -1 not used
Conf=function(target,pred,D=0.5,TC=-1)
{
 L=levels(target[1])
 if(is.vector(pred)) 
 { if(length(L)>2) pred<-factor(pred,levels=L)
   else  { 
           if(TC==1) LB=c("TRUE","FALSE") else LB=c("FALSE","TRUE")
           pred=factor(pred>D,levels=LB); target=factor((target==L[TC]),levels=LB)
         }
 }
 else if(is.factor(pred) && TC>0) 
 {
   pred=factor(pred==L[TC],levels=c("FALSE","TRUE"))
   target=factor((target==L[TC]),levels=c("FALSE","TRUE"))
 }
 else if(!is.factor(pred)) #if(is.matrix(pred) || is.data.frame(pred)) 
  { if(TC>0) { pred=factor(pred[,TC]>D,levels=c("FALSE","TRUE")); target=factor((target==L[TC]),levels=c("FALSE","TRUE"));}
    else pred=majorClass(pred,L)
  }
 return(table(target,pred))
}

# classification: 
# y - vector of factor with the desired values 
# x - vector of factor with the predicted values
# D - decision thresold
# TC - target concept class, -1 not used
Accuracy<-function(y,x,D=0.5,TC=-1)
{
 conf=Conf(y,x,TC=TC,D=D)
 D=0; for(i in 1:NCOL(conf))D=D+conf[i,i]
 return (100*D/length(y))
}

Kappa=function(y,x,D=0.5,TC=-1)
{
 conf=Conf(y,x,D,TC);Total=sum(conf);Diag=0;DiagR=0
 for(i in 1:NCOL(conf)) 
    {
      Diag<-Diag+conf[i,i]
      DiagR<-DiagR+(sum(conf[i,])*(sum(conf[,i])/Total))
    }
 return (100*(Diag-DiagR)/(Total-DiagR))
}

# ------------------------------------------------------------------
# call of the ROC function: 
# - calls multiROCcurve: if x is matrix or data.frame!
# - calls ROCcurve: else.
# y - vector of factor or numeric (0,1) with the desired values 
# x - vector or matriz of numeric with the predicted values
# TC - target class
ROCcurve<-function(y,x,TC=-1) #,method="int")
{
 if(TC>0 && !is.vector(x) ) { x=x[,TC];} else TC=2
 NC=NCOL(x)
 if(NC>2) return (multiROC(y,x)) #,method=method))
 else{
       if(is.factor(y[1])) POS=levels(y[1])[TC] else POS=1
       if(NC==2) x=x[,2]
       #cat(" >> POS:",POS,"lev:",levels(y),"\n")
       return (twoclassROC(y,x,Positive=POS)) #,method=method))
     }
}

# ------------------------------------------------------------------
LIFTcurve<-function(y,x,TC)
{
 if(TC>0 && !is.vector(x) ) { x=x[,TC];} else TC=2
 NC=NCOL(x)
 if(is.factor(y[1])) POS=levels(y[1])[TC] else POS=1
 if(NC==2) x=x[,2]
 NR=NROW(x); if(NR>100) NR=100
 alift=twoclassLift(y,x,Positive=POS,STEPS=NR,type=3)
 return (list(alift=alift,area=curvearea(alift,1)))
}

# ------------------------------------------------------------------
# Provost and Domingos AUC formulation for Multiclass problems
# y - vector of factor with the desired values 
# x - matriz of numeric with the predicted values
#     Note: the sum of x[i,] should be 1 for all i!!!
multiROC<-function(y,x) #,method="int")
{
 C=NCOL(x) # number of classes
 ROC<-vector("list",C)
 # prevalence of each class:
 SUM<-length(y)
 Lev<-levels(y[1])
 p=table(y)[]/SUM
 aux=0.0
 for(i in 1:C)
   { #print(paste("i:",i))
     R<-twoclassROC(y,x[,i],Positive=Lev[i]) #,method=method)
     ROC[[i]]=R
     if(p[i]>0) aux=aux+R$auc*p[i]
     #cat("i:",i,"p:",p[i],"auc:",R$auc,"aux:",aux,"\n")
   }
  ROC<-list(roc=ROC,auc=aux)
  return (ROC) # use: ob$roc[[i]]$roc or ob$roc[[i]]$auc to access individual rocs for each class i
}

# ------------------------------------------------------------------
# practical efficient method for ROC and AUC value
# algorithm 2 of Fawcett 2003, algorithm 1 of Fawcett 2006
# notes: use only with 2 classes
#        this is 2nd implementation, where the <-c(,) was replaced
#        by a much faster [,1]<- and [,2]<- instructions 
#
# y - vector of factor/numeric with the desired values 
# x - vector of numeric with the predicted values
# Positive - a label or number that corresponds to a TRUE/positive class value
## method = "int" - interpolate between 2 points, "pes" - pessimistic Fawcett point, "opt" - optimistic Fawcett point
# ------------------------------------------------------------------
twoclassROC<-function(y, x, Positive=1) #, method="int")
{
 #print(method)
# if(is.factor(y)) {y=as.numeric(y)-1;Positive=1;}
#  print(summary(y))
#  print(summary(x))

#YY<<-y; XX<<-x

  Xsize<-length(y)
  Pos<-sum(y[]==Positive) # total actual positives
#PP<<-Positive
#cat("Pos:",Pos,"\n")
  Neg<-Xsize-Pos          # total actual negatives

  Ord<-order(x,decreasing=TRUE) # very fast sort of vector
 
  FP<- 0
  TP<- 0
  FPprev<- 0
  TPprev<- 0
  A<-0  
  fprev<- -Inf

#cat(" --- AUC:",A,"\n")
  R<-matrix(ncol=2,nrow=(Xsize+1))
  k<-1
  for(i in 1:Xsize) 
     {
       if (x[Ord[i]]!=fprev)  
            { 
              if(FP>0) R[k,1]<-FP/Neg else R[k,1]=0
              if(TP>0) R[k,2]<-TP/Pos else R[k,2]=0 # the ROC point
#cat("k:",k,"FP:",FP,"TP:",TP,"Neg:",Neg,"Pos:",Pos,"\n")
              if(!is.na(R[k,1])) 
              {k<-k+1
               A<-A+trap_area(FP,FPprev,TP,TPprev) #,method) # compute the AUC
               fprev <- x[Ord[i]]
               FPprev<- FP
               TPprev<- TP
              }
            }
       if (y[Ord[i]]==Positive) TP<-TP+1 # test[i]
       else FP<-FP+1 
     }
  if(FP>0) R[k,1]<-FP/Neg else R[k,1]=0
  if(TP>0) R[k,2]<-TP/Pos else R[k,2]=1 # the ROC point
#cat("k:",k,"FP:",FP,"TP:",TP,"Neg:",Neg,"Pos:",Pos,"\n")
  if(FP==0 && k<(Xsize+1)) {k=k+1; R[k,]=c(1,1)}
#cat(" --- AUC:",A,"pos:",Pos,"neg:",Neg,"TP",TP,"FP",FP,"\n")
#cat(" ---: FPprev:",FPprev,"TPprev:",TPprev,"trap:",trap_area(1,FPprev,Pos,Pos),"\n")
  if(Neg>0) A<-A+trap_area(Neg,FPprev,Pos,TPprev) #,method)
  else A<-A+trap_area(1,FPprev,Pos,Pos) #,method)
  if(Pos>0 && Neg>0) A<-A/((1.0*Pos)*Neg) 
  else if(Neg>0) A<-A/((1.0*Neg))
  else A<-A/((1.0*Pos))

#RR<<-R;RR=na.omit(RR);print(RR[])
 
  ROC<-list(roc=R[(1:k),],auc=A)
#cat(" --- AUC:",A,"\n")
  return (ROC)
}
# ------------------------------------------------------------------
# internal R function used by ROCcurve: do not use this
# ------------------------------------------------------------------
trap_area<-function(X1,X2,Y1,Y2) #,method="int")
{ 
  return ( (abs(X1-X2)) * ((Y1+Y2)/2) )
  #if(method=="int") return ( (abs(X1-X2)) * ((Y1+Y2)/2) )
  #else if(method=="opt") return ( (abs(X1-X2)) * ((Y1+Y1)/2) ) 
  #else if(method=="pes") return ( (abs(X1-X2)) * ((Y2+Y2)/2) ) 
}
# ------------------------------------------------------------------
xmiddle_point<-function(X1,X2,Y1,Y2,X3)
{ 
 m=(Y2-Y1)/(X2-X1); b=Y1-m*X1;
 return (m*X3+b)
}
#-------------------------------------------------------------------
# vertical averaging of ROC curves, algorithm 3 from Fawcett 2006
# samples - number of FP samples
# ROCS list with length(ROCS) ROC curves, each ROC is [,1] frp and [,2] tpr
#
# returns tpravg with samples+1 rows and 3+nrocs columns: fpr, tpr, mean, confint95, tpr_roc[[1]],...,tpr_roc[[nrocs]]
vaveraging<-function(samples,ROCS,min=0,max=1)
{
  s=1
  nrocs=length(ROCS)
  fprsamples=seq(min,max,length.out=samples)
#cat("min:",min,"max:",max,"\n")
  tpravg=matrix(ncol=(3+nrocs),nrow=length(fprsamples))
  for(k in fprsamples)
  {
    #tprsum=0
    tprsum=rep(0,nrocs)
    for(i in 1:nrocs)
    {
#cat("k:",k,"i:",i,"\n")
     #tprsum=tprsum+TPR_FOR_FPR(k,ROCS[[i]],nrow(ROCS[[i]]))
     tprsum[i]=tprsum[i]+TPR_FOR_FPR(k,ROCS[[i]],nrow(ROCS[[i]]))
#cat("tprsum[",i,"]=",tprsum[i],"\n")
    }
    #tpravg[s,]=c(k,tprsum/nrocs)
#cat("conf:\n")
#TPR<<-tprsum
    tpravg[s,]=c(k,mean(tprsum),conflevel(tprsum),tprsum)
#cat("conf done\n")
    s=s+1
  } 
  return(tpravg)
}
# internal R functions used by vaveraging: do not use this
TPR_FOR_FPR<-function(fprsample,ROC,npts)
{
 # error here, think later...
 #RRR<<-ROC
 #NPTS<<-npts
 #FPR<<-fprsample
 i=1
 NR=nrow(ROC)
 while(i<npts && ROC[(i+1),1]<=fprsample) i=i+1;
 if(i<NR && ROC[i,1]==fprsample) return (ROC[i,2])
 else if(i<NR) return (INTERPOLATE(ROC[i,],ROC[(i+1),],fprsample))
 else return (ROC[i,2])
# else return (INTERPOLATE(ROC[(NR-1),],ROC[NR,],fprsample))
}
INTERPOLATE<-function(ROCP1,ROCP2,X)
{
 #cat("rocp1:",ROCP1,"rocp2:",ROCP2,"x:",X,"\n")
 slope=(ROCP2[2]-ROCP1[2])/(ROCP2[1]-ROCP1[1])
 return (ROCP1[2]+slope*(X-ROCP1[1]))
}
# 95% confidence interval according to a t-student distribution
conflevel=function(x,level=0.95)
{
 RES=try( (t.test(x,conf.level=level)$conf[2]-t.test(x,conf.level=level)$conf[1])/2 , silent=TRUE)
 if(class(RES)=="numeric") return(RES) else return (0)
}

# partial curve (roc, rec, ...)
partialcurve=function(Curve,threshold=1) #,method="int")
{
 I=which(Curve[,1]<=threshold)
 IND=I[length(I)]
 if(Curve[IND,1]==threshold) M=Curve[(1:IND),]
 else 
 {
   M=Curve[(1:(IND+1)),]
   M[(IND+1),]=c(threshold,xmiddle_point(Curve[IND,1],Curve[(IND+1),1],Curve[IND,2],Curve[(IND+1),2],threshold))
   ##if(method=="int") M[(IND+1),]=c(threshold,xmiddle_point(Curve[IND,1],Curve[(IND+1),1],Curve[IND,2],Curve[(IND+1),2],threshold))
   ##else if(method=="pes") M[(IND+1),]=c(threshold,Curve[IND,2])
   ##else if(method=="opt") M[(IND+1),]=c(threshold,Curve[(IND+1),2])
 }
 return(M)
 #return(list(curve=M,auc=rocarea(M,threshold),tprfpr=M[(IND+1),2]))
}

# area of a curve using trapesoidal method
# examples: auc of a roc curve, area of rec, etc...
curvearea<-function(Curve,threshold=1.0) #,method="int")
{ 
  if(is.vector(Curve)) return (0)
  else 
  { if(Curve[nrow(Curve),1]>threshold) Curve=partialcurve(Curve,threshold) #,method=method)
    A=0
    for(i in 2:nrow(Curve)) 
       {
        A=A+trap_area(Curve[i,1],Curve[(i-1),1],Curve[i,2],Curve[(i-1),2]) #,method=method)
#cat("A",A,"T",threshold,"\n")
       }
#cat("A",A,"A/T",A/threshold,"\n")
    if(threshold>0) return (A/threshold)
    else return (0)
  }
}

#-------------------------------------------------------------------

# tolerance of a rec curve
tolerance<-function(REC,tol=0.5)
{
 stop=FALSE; i=1;N=nrow(REC)
 while(i<N && REC[i,1]< tol ) {i=i+1;}

 if(i==N || REC[i,1]==tol) return (REC[i,2])
 else if(i>1) return ( xmiddle_point(REC[(i-1),1],REC[i,1],REC[(i-1),2],REC[i,2],tol) )
}

# mean and confidence interval using t.test
meanint<-function(x,level=0.95)
{
 if(is.matrix(x)||is.data.frame(x))
 {
  C=ncol(x); M=rep(0,C); Conf=rep(0,C);
  for(i in 1:C)
  { M[i]=mean(x[,i])
    Conf[i]=conflevel(x[,i],level=level)
  }
 }
 else
 {
  M=mean(x)
  Conf=conflevel(x,level=level)
 }
 return(list(mean=M,int=Conf))
}
# --------

# m - is matrix or data.frame
mpairwise=function(m,p.adj="bonf",paired=TRUE)
{
 NC=NCOL(m);NR=NROW(m)
 x=vector(length=NC*NR); g=x;
 for(i in 1:NC) { ini=(i-1)*NR+1;end=ini+NR-1;
                  x[ini:end]=m[,i];g[ini:end]=rep(i,NR);}
 g=factor(g)
 P=pairwise.t.test(x,g,p.adj =p.adj,paired=paired)
 return(P)
}


# TC - target concept class, -1 not used
## method="int" # not used currently...
# M - mining or target
# X - NULL or predictions
# b - baseline (for forecasting errors)
mmetric=function(y,x=NULL,metric,D=0.5,TC=-1,val=NULL) #,method="int")
{
 if(!is.null(x)) RES=Error(y,x,metric,D=D,TC=TC,val=val) 
 else if(is.list(y)){ R=y$runs
        if(metric=="CONF") {  
 	                     for(i in 1:R)  
  	                      { #cat("i:",i,"\n")
   		                if(i==1) RES=Error(y$test[[i]],y$pred[[i]],metric,D=D,TC=TC,val=val)
                                else RES=RES+Error(y$test[[i]],y$pred[[i]],metric,D=D,TC=TC,val=val)
  	                      }
                           }
        else 
        { RES=vector(length=R) 
 	  for(i in 1:R)  
  	  { #cat("i:",i,"\n")
   		RES[i]=Error(y$test[[i]],y$pred[[i]],metric,D=D,TC=TC,val=val)
  	  }
        }
     }
 else RES=NA
 return(RES)
}

# experimental stuff:
# type 1 - normal lift
# type 2 - cumulative
# type 3 - cumulative percentage
twoclassLift<-function(y, x, Positive=1,STEPS=10,type=3)
{
#YY<<-y;XX<<-x;PP<<-Positive;
#y=YY;x=XX;Positive=PP;STEPS=10;type=3
  #STEPS=11;Positive="setosa";type=3
  DIV=STEPS^2
  Xsize<-length(y)
  APos<-sum(y[]==Positive) # total actual positives
  if(is.data.frame(x)) x=x[,1] 
  Ord<-order(x,decreasing=TRUE) # very fast sort of vector
  ALL=APos/Xsize
  if(type==3){R=matrix(ncol=2,nrow=STEPS+1)
              R[1,]=c(0,0)
             }
  else R=matrix(ncol=2,nrow=STEPS)
  Portion=Xsize/STEPS
  for(i in 1:STEPS) 
     {
       N=round(i*Portion)
       if(type==1 || type==3) IND=1:N # total actual positives
       else if(type==2) IND=(N-STEPS+1):N
       Pos=sum(y[Ord[IND]]==Positive)
       if(type==3) { R[(i+1),1]=i*STEPS/DIV; R[(i+1),2]=Pos/APos }
       else {R[i,1]=i*STEPS/DIV; R[i,2]=Pos/(ALL*N)}
     }
  return (R)
}
