#-------------------------------------------------------------------------------------------------
# "estimate.R" code by Paulo Cortez 2009@, Department of Information Systems, University of Minho
#
# This file deals will two estimation functions: crossvaldata - k-fold cross validation and holdout - houldout validation
#
#-------------------------------------------------------------------------------------------------

# adapted from the bootstrap library to use x - formula and data!!!
# order = if samples are timely ordered and this order should be kept:
crossvaldata<- function(x,data,theta.fit,theta.predict,ngroup=n,order=FALSE,model,task,feature="none",...)
{
    call <- match.call()
    outindex<-output_index(x,data)
    y<- data[,outindex]
    n <- length(y)
    ngroup <- trunc(ngroup)

    if( ngroup < 2){
        stop ("ngroup should be greater than or equal to 2")
    }
    if(ngroup > n){
        stop ("ngroup should be less than or equal to the number of observations")
    }
  
    if(ngroup==n) {groups <- 1:n; leave.out <- 1}
    if(ngroup<n){
     leave.out <- trunc(n/ngroup)
     groups <- vector("list",ngroup)

     if(is.factor(y)) # stratified crossvalidation
     { L<-levels(y)
       NL<-length(L)
       I<-vector("list",NL)
       g<-vector("list",NL)
       for(i in 1:NL)
          {
            I[[i]]<-which(y==L[i])
            g[[i]]<-crossfolds(y[ I[[i]] ],ngroup)
          }
       for(j in 1:ngroup)
          {
           f=NULL
           for(i in 1:NL) f=c(f, I[[i]][ g[[i]][[j]] ] )
           groups[[j]]=f
          }
      }
      else #normal crossvalidation
      {
        if(order) o=1:n
        else o <- sample(1:n)
        for(j in 1:(ngroup-1)){
            jj <- (1+(j-1)*leave.out)
            groups[[j]] <- (o[jj:(jj+leave.out-1)])
        }
        groups[[ngroup]] <- o[(1+(ngroup-1)*leave.out):n]
      }
    }
    u <- vector("list",ngroup)
    npar=modelnpar(model)
    par=matrix(nrow=ngroup,ncol=npar)
    if(substr(feature[1],1,4)=="sens" || substr(feature[1],1,4)=="sabs" || substr(feature[1],1,4)=="simp") 
    {
         SEN <- matrix(nrow=ngroup, ncol=ncol(data) )
         if(substr(feature[1],1,4)!="sens")
         { SRESP=TRUE;NCOL=ncol(data);FSRESP=TRUE;}
         else
         { SRESP=FALSE;FSRESP=FALSE;}
    }
    else {SEN<-NULL;SRESP=FALSE}

    if(feature[1]!="none" && substr(feature[1],1,4)!="sens") ATTRIB <- vector("list",ngroup)
    else ATTRIB<-NULL

    #cv.fit <- rep(NA,n) # this line was changed to:
    if(task=="prob") cv.fit <- matrix(ncol=length(levels(y)),nrow=n) 
    else if(task=="class") cv.fit=y
    else cv.fit <- rep(NA,n)

    imethod=switch(feature[1],sabsv=,simpv="sensv",sabs=,simp=,simpg="sensg",sabsr=,simpr="sensr",feature[1])

    for(j in 1:ngroup)
    {
        u <- theta.fit(x,data[-groups[[j]], ],task=task,model=model,feature=feature,...) # is this the correct line or add feature=feature?
        if(!is.null(SEN)) 
         {
            #cat("----- j:",j,"\n",sep=" ")
            #aux=(Importance(u,data[-groups[[j]], ],method="sensgrad"))$imp 
            #print(paste(" --- j:",j,"---"))
            #print(aux)
            IMPORTANCE=Importance(u,data[-groups[[j]], ],method=imethod,responses=FSRESP)
#cat(" >> L:",length(IMPORTANCE$sresponses),"NCOL:",NCOL,"\n")
            SEN[j,]=IMPORTANCE$imp 
            if(FSRESP)
            { 
             if(j==1) SRESP=IMPORTANCE$sresponses
             else{ for(l in 1:NCOL) 
                    {
                     if(!is.null(IMPORTANCE$sresponses[[l]]) ) # && !is.null(IMPORTANCE$sresponses[[j]])$x) # $x)) 
                      { if(is.null(SRESP[[l]])) SRESP[[l]]=IMPORTANCE$sresponses[[l]]
                        else{ SRESP[[l]]$x=c(SRESP[[l]]$x,IMPORTANCE$sresponses[[l]]$x);
                              if(task=="prob") SRESP[[l]]$y=rbind(SRESP[[l]]$y,IMPORTANCE$sresponses[[l]]$y,deparse.level=0)
                              else if(task=="class") SRESP[[l]]$y=addfactor(SRESP[[l]]$y,IMPORTANCE$sresponses[[l]]$y)
                              else SRESP[[l]]$y=c(SRESP[[l]]$y,IMPORTANCE$sresponses[[l]]$y)
                            }
                      }
                    }
                }
            }
         }
        if(!is.null(ATTRIB)) ATTRIB[[j]]=u@attributes
        #print(u)
        if(npar>0) par[j,]=as.numeric(u@mpar[1,])

#cat("cv fit class:",class(cv.fit),"\n")
        if(is.matrix(cv.fit)) cv.fit[ groups[[j]], ] <-  theta.predict(u,data[groups[[j]],]) # probabilities!
        else cv.fit[ groups[[j]] ] <-  theta.predict(u,data[groups[[j]],]) # regression or classification, 1 output
        #print(cv.fit[groups[[j]]])
    }
    if(leave.out==1) groups <- NULL
    return(list(cv.fit=cv.fit, 
                mpar=par,
                sen=SEN,sresponses=SRESP,
                attributes=ATTRIB,
                ngroup=ngroup, 
                leave.out=leave.out,
                groups=groups, 
                call=call)) 
}

#  auxiliar function, adapted from the bootstap library: only makes the groups, you should not need to use this:
crossfolds<-function(y,ngroup)
{
  n <- length(y)
  leave.out <- trunc(n/ngroup)
  groups <- vector("list",ngroup)
  o <- sample(1:n)
  for(j in 1:(ngroup-1)){
            jj <- (1+(j-1)*leave.out)
            groups[[j]] <- (o[jj:(jj+leave.out-1)])
  }
  groups[[ngroup]] <- o[(1+(ngroup-1)*leave.out):n]
  return(groups)
}


#---------------------------------------------------------
# holdout: create indexes for spliting the data into training and test datasets
#          the holdout is statified if the output is factor 
# a list is returned with:
#  $tr - indexes of all training examples
#  $ts - indexes of all test examples
#  if internalsplit is TRUE then another holdout if performed on tr:
#  $itr - indexes of tr for internal training 
#  $vtr - indexes of tr for internal testing (validation)
#
# 
# parameters:
# y - is the desired output vector or a vector with a numeric sequence (e.g. c(1,1,2,3,3,3,4,4) )
# ratio is the ratio of training set (in percentage)
# internalsplit if TRUE then another stratified holdout is used within the training set
# mode - the sampling mode used for the holdout
#      > "random" - is the default mode, uses standard random split
#      > "order" - uses the sequencial order of y, no random is used. 
#                  the first examples are used as tr while the lattest are used as ts
#      need to check in the future is this makes any sense at all ?
#      > "incremental" - incremental retraining, ratio=batch-size, iter=iterator
#---------------------------------------------------------
holdout<-function(y,ratio=2/3,internalsplit=FALSE,mode="random",iter=1)
{ 
  ALLITR=NULL; VAL=NULL;
  NSIZE=length(y)
 
 if(mode=="incremental")
 { batches=ceiling(NSIZE/ratio)-1
   aux=iter*ratio; ALLTR=1:aux;
   end=aux+ratio; if(end>NSIZE) end=NSIZE;
   TS=(length(ALLTR)+1):end
 }
 else 
 {
  if(ratio>=1) { ratio=(NSIZE-ratio)/NSIZE}

  if(mode=="order")
  { 
    TRS<-round(ratio*NSIZE)
    ALLTR<-1:TRS
    TS<-(TRS+1):NSIZE
    if(internalsplit)
      {
       TRIS<-round(ratio*TRS)
       ALLITR<-1:TRIS
       VAL<-(TRIS+1):TRS
      }
  }
  else # default random holdout
  {
   ALL<-1:NSIZE
   FACTOR<-is.factor(y)
   ALLITR<-vector(length=0)
    if(FACTOR) 
     { L<-levels(y)
       NL<-length(L)
       I<-vector("list",NL)
       TR<-vector("list",NL)
       ALLTR<-vector(length=0)

       for(i in 1:NL)
          {
            I[[i]]<- which(y==L[i])
            NI<- length( (I[[i]]) )
            TR[[i]]<-sample( (I[[i]]), (NI*ratio) )
            ALLTR<-c(ALLTR,TR[[i]])
 
            if(internalsplit==TRUE)
               {
                 NTR<-length(TR[[i]]) 
                 ITR<-sample(TR[[i]],(NTR*ratio))
                 ALLITR<-c(ALLITR,ITR)
               }
          }
      }
    else
     {
       ALLTR<-sample(ALL,(ratio*NSIZE)) 
       if(internalsplit)
         {
          NALLTR<-length(ALLTR)
          ALLITR<-sample(ALLTR,(ratio*NALLTR)) 
          
         }
     }
  TS<-setdiff(ALL,ALLTR)
  if(internalsplit==TRUE)VAL<-setdiff(ALLTR,ALLITR)       
  }
 }
  return(list(tr=ALLTR,itr=ALLITR,val=VAL,ts=TS))
}
#-----------------------------------------------------------
