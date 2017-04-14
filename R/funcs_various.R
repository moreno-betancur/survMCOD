#' Function to get entry indicators
#'
#' Internal function used to deal with delayed entry
#' @export
#' @keywords internal
getEntry<-function(x,tt){
  tt1<-tt-x["AgeEntry"]
  return(ifelse(length(tt1[tt1>0])==0,length(tt)+1,
                min(which(tt1==min(tt1[tt1>0])))))}

#' Function to get exit indicators
#'
#' Internal function used to deal with delayed entry
#' @export
#' @keywords internal
getExit<-function(x,tt){
  tt1<-x["AgeExit"]-tt
  return(ifelse(length(tt1[tt1>=0])==0,0,
                min(which(tt1==min(tt1[tt1>=0])))))}

#' Function to get at-risk indicators
#'
#' Internal function used to deal with delayed entry
#' @export
#' @keywords internal
atrisk<-function(x,tt){
  if(x["IndEntry"]==(length(tt)+1)|x["IndExit"]==0) st<-rep(0,length(tt)) else
    if(x["IndEntry"]!=1&x["IndExit"]!=length(tt)) st<-c(rep(0,x["IndEntry"]-1),rep(1,x["IndExit"]-x["IndEntry"]+1),rep(0,length(tt)-x["IndExit"])) else
      if(x["IndEntry"]==1&x["IndExit"]!=length(tt)) st<-c(rep(1,x["IndExit"]),rep(0,length(tt)-x["IndExit"])) else
        if(x["IndEntry"]!=1&x["IndExit"]==length(tt)) st<-c(rep(0,x["IndEntry"]-1),rep(1,x["IndExit"]-x["IndEntry"]+1)) else
          if(x["IndEntry"]==1&x["IndExit"]==length(tt)) st<-rep(1,length(tt))
          return(st)}

#' Non-parametric estimator of log-ratio of the baseline hazards
#'
#' Internal function used to generate an starting value for this parameter for optimization
#' @export
#' @keywords internal

stat<-function(datc,ind,nxi)
{ datm<-datc[ind,]
dat<-datm
ret<-rep((log((sum(dat$status==1))/
                (sum(dat$status==2)))),nxi)
return(ret)
}


#' Inverse of the total cumulative hazard in simulation model
#'
#' Internal function used used for simulating data
#' @export
#' @keywords internal

cumhazinv<-function(u,X1,Z1,xi,rho,phi,pgen,lambda,v){
  Qgen<-length(pgen)
  cumhinv<-(-log(u)/(lambda*(sum(pgen)*exp(rho*X1)+
                               (Qgen-sum(pgen))*exp(-xi+phi*Z1))))^(1/v)
  return(cumhinv)
}

