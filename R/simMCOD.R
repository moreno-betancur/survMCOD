#' Simulate survival data with multiple causes of death
#'
#' \code{simMCOD} simulates data from the multiple-cause model of Moreno-Betancur et al. (2017),
#' specifically under the same assumptions as in their simulation study.
#'
#' @export
#'
#' @param n Size of dataset to be generated
#' @param xi Log-ratio of the pure baseline hazards, assumed constant
#' @param rho Effect of the binary exposure X1(=Z1) on the pure hazard of the disease of interest
#' @param phi Effect of the binary exposure Z1(=X1) on the pure hazard of other diseases
#' @param pgen Vector of weights corresponding to each of the 6 states in the multistate model.
#' @param lambda Weibull scale parameter used to specify the pure baseline hazard of the disease of interest.
#' @param v Weibull shape parameter used to specify the pure baseline hazard of the disease of interest.
#' @param pUC Vector indicating the values in \code{pgen} corresponding to weight values assigned to the underlying cause of death.
#'
#' @details The function can be used to generate data from any of the scenarios considered in the the main simulation settings
#' of Moreno-Betancur et al. (2017). See that reference for details.
#' @return A data.frame as required by \code{\link[survMCOD]{survMCOD}} function.
#' Specifically, a dataset with one row per individual and the following variables:
#'
#'\describe{
#'   \item{X1}{Binary exposure (to include in model for disease of interest).}
#'   \item{Z1}{Binary exposure equal to X1 (to include in model for other diseases).}
#'   \item{TimeEntry}{Time of entry into the study}
#'   \item{TimeExit}{Time of exit from the study}
#'   \item{Status}{Indicator of vital status, with Status=1 if individual died and Status=0 if individual is censored.}
#'   \item{Pi}{Weight indicating the proportion of the death attributed to the disease of interest (missing, i.e. "NA", if Status=0).}
#'   \item{UC}{The Underlying Cause Indicator, which is 1 if the individual died and the disease of interest was selected as
#' underlying cause of death and 0 otherwise (i.e. 0 if individual censored or if individual died but the disease of interest
#' was not selected as underlying cause of death).}
#' }
#'@references
#'
#'Moreno-Betancur M, Sadaoui H, Piffaretti C, Rey G. Survival analysis with multiple causes of death:
#'Extending the competing risks model. Epidemiology 2017; 28(1): 12-19.
#'
#'
#'@examples
#'
#'   datEx<-simMCOD(n=400,xi=-1,rho=-2,phi=0,
#'           pgen=c(1,0,0.75,0.25,0.125,0.083),
#'           lambda=0.001,v=2,pUC=c(1,0.75))
#'
#'   head(datEx)
#'
#'@importFrom utils flush.console
#'@importFrom stats rbinom runif rmultinom terms model.matrix model.extract quantile as.formula optim pnorm rmultinom


simMCOD<-function(n,xi,rho,phi,pgen,lambda,v,pUC){

  # Generate balanced binary exposure
  X1<-rbinom(n,1,0.5)
  Z1<-X1

  # Generate failure times
  u<-runif(n,0,1)
  times<-cumhazinv(u,X1,Z1,xi,rho,phi,pgen,lambda,v)

  # Generate censoring times (results in around 30% censoring)
  cens<-runif(n,5,15)
  sum(1*(cens<=times)+0)/n

  # Observed right-censored time-to-event
  tt<-pmin(times, cens)

  # Determine event-type
  Qgen<-length(pgen)
  h1<-lambda*v*times^(v-1)*exp(rho*X1)
  h0<-lambda*v*times^(v-1)*exp(-xi+phi*Z1)
  allh<-sum(pgen)*h1+(Qgen-sum(pgen))*h0

  pr<-NULL
  for(ff in 1:Qgen)
    pr<-cbind(pr,(pgen[ff]*h1+(1-pgen[ff])*h0)/allh)

  event<-rep(NA,n)
  for( m in 1:n) {
    mtn<-rmultinom(1,size=1,prob =pr[m,])
    event[m]<-which(mtn==1)
  }
  case<-(times<=cens)*event+0
  Pi<-rep(NA,n)
  Pi[case!=0]<-pgen[case[case!=0]]
  Status<-1*(case!=0)+0
  #data
  dat<-data.frame(X1,Z1,rep(0,n),tt,Status,Pi)
  dat<-dat[order(dat[,"tt"]),]
  names(dat)[3:4]<-c("TimeEntry","TimeExit")
  dat$UC<-with(dat,ifelse(is.na(Pi),0,ifelse(Pi%in%pUC,1,0)))
  return(dat)}

