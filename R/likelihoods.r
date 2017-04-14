########################## Likelihoods functions ##########################
# All of these functions assume that datamat has in its first columns Z1,Z2,..,ZK X1,X2,...,XJ,
# in that order and that data is sorted in ascending failure time

#' Loglikelihood L as a function of (xi,rho,phi) to get starting values
#'
#' Internal function used in estimation procedure
#' @export
#' @keywords internal

loglik=function(par,nxi,nX,nZ,pset,datamat,psetTT,txi,tt,datInd) {###par contains (xi,rho,phi)

  xi=par[1:nxi]
  if(nX==0){rho=as.matrix(0)}else{rho=as.matrix(c(rep(0,times=nZ),par[(nxi+1):(nxi+nX)]))}
  if(nZ==0){phi=as.matrix(0)}else{phi=as.matrix(c(par[(nX+nxi+1):(nZ+nX+nxi)],rep(0,nX)))}

  ll<-0
  for(pp in pset)
  {
    q<-which(psetTT==pp)

    datXi<-NULL
    for(ff in 1:nxi)
      datXi<-cbind(datXi,pp*exp(datamat[,1:(nZ+nX)]%*%rho)+(1-pp)*exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi))

    # get denominator: sum of contributions of individuals at risk for each uniqute type-q event time
    ttq<-unique(datamat[datamat[,"status"]==q,"AgeExit"])
    cum<-vector()
    for(ff in 1:nxi)
    {tind<-(tt%in%ttq)&(txi[ff]<=tt)&(tt<ifelse(is.na(txi[ff+1]),Inf,txi[ff+1]))
    datq<-datInd[,tind]
    cum<-c(cum,as.matrix(t(datXi[,ff])%*%datq))}

    # get power of denominator: number of ties for each uniqute type-q event time
    ttq<-datamat[datamat[,"status"]==q,"AgeExit"]
    tt2<-table(ttq)
    tt2<-data.frame(ID=1:nrow(tt2),ttq=as.numeric(dimnames(tt2)$ttq),tt2,row.names=NULL)[,c(1,2,4)]

    ll<-ll-sum(log((cum)^tt2$Freq))

    for(ff in 1:nxi)
    {
    dp<-(datamat[,"status"]==q&txi[ff]<=datamat[,"AgeExit"]&datamat[,"AgeExit"]<ifelse(is.na(txi[ff+1]),Inf,txi[ff+1]))
    ll<-ll+sum(log(datXi[dp,ff]))
    }
  }
  return(-ll) #ll is negative because optim performs minimization
}

#' Loglikelihood L as a function of (rho,phi) only, with xi known from previous iteration
#' Exactly the same as loglik except now xi is not an argument to the function so the first three
#' lines change
#'
#' Internal function used in estimation procedure
#' @export
#' @keywords internal

loglikopt=function(par,nxi,nX,nZ,pset,datamat,psetTT,txi,tt,datInd,opt,j){ ###par contains (rho,phi)

  xi=opt[j,1:nxi]
  if(nX==0){rho=as.matrix(0)}else{rho=as.matrix(c(rep(0,nZ),par[1:nX]))}
  if(nZ==0){phi=as.matrix(0)}else{phi=as.matrix(c(par[(nX+1):(nZ+nX)],rep(0,nX)))}

  ll<-0
  for(pp in pset)
  {
    q<-which(psetTT==pp)

    datXi<-NULL
    for(ff in 1:nxi)
      datXi<-cbind(datXi,pp*exp(datamat[,1:(nZ+nX)]%*%rho)+(1-pp)*exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi))

    # get denominator: sum of contributions of individuals at risk for each uniqute type-q event time
    ttq<-unique(datamat[datamat[,"status"]==q,"AgeExit"])
    cum<-vector()
    for(ff in 1:nxi)
    {tind<-(tt%in%ttq)&(txi[ff]<=tt)&(tt<ifelse(is.na(txi[ff+1]),Inf,txi[ff+1]))
     datq<-datInd[,tind]
     cum<-c(cum,as.matrix(t(datXi[,ff])%*%datq))}

    # get power of denominator: number of ties for each uniqute type-q event time
    ttq<-datamat[datamat[,"status"]==q,"AgeExit"]
    tt2<-table(ttq)
    tt2<-data.frame(ID=1:nrow(tt2),ttq=as.numeric(dimnames(tt2)$ttq),tt2,row.names=NULL)[,c(1,2,4)]

    ll<-ll-sum(log((cum)^tt2$Freq))

    for(ff in 1:nxi)
    {
      dp<-(datamat[,"status"]==q&txi[ff]<=datamat[,"AgeExit"]&datamat[,"AgeExit"]<ifelse(is.na(txi[ff+1]),Inf,txi[ff+1]))
      ll<-ll+sum(log(datXi[dp,ff]))
    }
  }
  return(-ll) #ll is negative because optim performs minimization
}

#' Loglikelihood L* as a function of xi only, with rho and phi known from previous iteration
#' Only difference with L is denominator
#'
#' Internal function used in estimation procedure
#' @export
#' @keywords internal
#'

logliksteropt=function(par,nxi,nX,nZ,pset,datamat,psetTT,txi,tt,datInd,opt,i){ ###par contains (xi)
  xi=par
  if(nX==0){rho=as.matrix(0)}else{rho=as.matrix(c(rep(0,nZ),opt[i,(nxi+1):(nxi+nX)]))}
  if(nZ==0){phi=as.matrix(0)}else{phi=as.matrix(c(opt[i,(nX+nxi+1):(nZ+nX+nxi)],rep(0,nX)))}

  ll<-0
  Q<-length(pset)
  datXiT<-NULL
  for(ff in 1:nxi)
    datXiT<-cbind(datXiT,sum(pset)*exp(datamat[,1:(nZ+nX)]%*%rho)+(Q-sum(pset))*exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi)) #aT

  for(pp in pset)
  {
    q<-which(psetTT==pp)


    datXi<-NULL
    for(ff in 1:nxi)
      datXi<-cbind(datXi,pp*exp(datamat[,1:(nZ+nX)]%*%rho)+(1-pp)*exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi))


    # get denominator: sum of contributions of individuals at risk for each uniqute type-q event time
    ttq<-unique(datamat[datamat[,"status"]==q,"AgeExit"])
    cum<-vector()
    for(ff in 1:nxi)
    {tind<-(tt%in%ttq)&(txi[ff]<=tt)&(tt<ifelse(is.na(txi[ff+1]),Inf,txi[ff+1]))
     datq<-datInd[,tind]
     cum<-c(cum,as.matrix(t(datXiT[,ff])%*%datq))}

    # get power of denominator: number of ties for each uniqute type-q event time
    ttq<-datamat[datamat[,"status"]==q,"AgeExit"]
    tt2<-table(ttq)
    tt2<-data.frame(ID=1:nrow(tt2),ttq=as.numeric(dimnames(tt2)$ttq),tt2,row.names=NULL)[,c(1,2,4)]

    ll<-ll-sum(log((cum)^tt2$Freq))

    for(ff in 1:nxi)
    {
      dp<-(datamat[,"status"]==q&txi[ff]<=datamat[,"AgeExit"]&datamat[,"AgeExit"]<ifelse(is.na(txi[ff+1]),Inf,txi[ff+1]))
      ll<-ll+sum(log(datXi[dp,ff]))
    }
  }
  return(-ll)

}


