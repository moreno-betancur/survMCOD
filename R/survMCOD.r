#' Fit Cox regression models for multiple cause-of-death data
#'
#' \code{survMCOD} fits Cox regression models for the pure hazard of death due to a disease of interest based on multiple
#' cause-of-death data ("Multiple-cause analysis"), thus acknowledging that death may be caused by several disease processes acting
#' concurrently. The pure hazard is the rate of deaths caused exclusively by the disease of interest, and is thus a quantity that
#' is conceptually closer to the marginal "causal" hazard than the cause-specific hazard. The latter is the quantity modeled
#' when using competing risks Cox regression based on the so-called "underlying cause of death" and ignoring all other diseases
#' mentioned on the death certificate ("Single-cause analysis").
#'
#'
#' @export
#' @param formula A formula object with the response on the left of the ~ operator being an object as returned by the \code{SurvM} function
#' (see \strong{Examples} section below). The terms on the right are the regressors to include in the the model for the pure hazard of the
#' disease of interest.
#' @param formOther A formula object with empty response, i.e. nothing on the left of the ~ operator (if a response is given it is ignored).
#' The terms on the right are the regressors to include in the the model for the pure hazard of other diseases. If unspecified, this model will
#' include the same regressors as the model for the disease of interest as specified in \code{formula} above.
#' @param data A data.frame in which to interpret the variables named in \code{formula} and \code{formOther} above.
#' Specifically, the dataset must contain one row per individual and the following variables:
#' - The time-to-event (for type="right") or the entry and exit times (for type="counting")
#' - A status variable indicating whether the individual died (=1) or is censored (=0)
#' - A weight variable which is missing ("NA") if the individual is censored, and a proportion between 0 and 1 if the individual died.
#' The weight corresponds to the proportion of that death that is attributed to the disease of interest according to the chosen weight
#' attribution strategy (see \strong{Details} section below).
#' - All the regressors named in \code{formula} and \code{formOther}.
#' - The Underlying Cause Indicator, specified in \code{UC_indicator}, which is 1 if the individual died and the disease of interest was selected as
#' underlying cause of death and 0 otherwise (i.e. 0 if individual censored or if individual died but the disease of interest
#' was not selected as underlying cause of death).
#' @param UC_indicator Name of the Underlying Cause Indicator variable in the dataset.
#' @param Iter Number of iterations for iteration procedure. Default is 4 which is generally enough to achieve convergence. See \strong{Value} below
#' for more details.
#' @details The \code{survMCOD} function can be used to fit Cox regression models
#' for the pure hazard of death due to a disease of interest based on multiple
#' cause-of-death data. In addition to results from the multiple-cause analysis,
#' the function also provides the results from the single-cause analysis (i.e. from a
#' competing risks Cox regression based on the so-called "underlying cause of death").
#'
#' The key preliminary step to using this function to perform the multiple-cause analysis
#' is to assign a weight to each death that represents the proportion of the death
#' attributed to the disease of interest. The user is referred to Moreno-Betancur et al. (2017),
#' Piffaretti et al. (2016) and Rey et al. (2017) for descriptions and discussions of various
#' weight-attribution strategies.
#'
#' The assumptions of the multiple-cause model and details of the estimation procedure are provided in
#' Moreno-Betancur et al. (2017). A key feature of the model is that regression coefficients
#' need to be estimated simultaneaously for a Cox model for the disease of interest and a
#' Cox model for other causes, and deaths with a weight between 0 and 1 will contribute to both.
#' This is why the user needs to specify regressors for each of these models using the two
#' arguments \code{formula} and \code{formOther}.
#'
#' Another aspect of the multiple-cause model is that a fully parametric model for the log ratio of the
#' baseline pure hazards needs to be posited. The current default is to parametrise this
#' a piecewise constant function with cut-offs at the 25th, 50th and 75th percentile of the
#  failure time distribution of pure events (with weight=1). Future versions will provide
#' the user control over this.
#'
#' Convergence of the multiple-cause model fitting procedure should be checked
#' using \code{check.survMCOD}.
#'
#'@return This development version returns a list with three components. First, a list with the results of the multiple-cause analysis.
#'These are estimates of log hazard ratios for the models for the disease of interest and other diseases, and
#'estimates of the piecewise constant log ratio of the baseline pure hazards. Second, a list with the single-cause analysis
#'log hazard ratio estimates for the disease of interest and other diseases. All estimates are accompanied with their corresponding
#'standard errors, 95\%  confidence intervals and p-values. This will change in future versions when a proper class of objects and
#'summary and other such methods are developed. Third, a data.frame with Iter+2 estimates of each parameter, the first two corresponding
#'to a starting value and an improved starting value (see Moreno-Betancur et al. 2017 for details). This data.frame is used by function
#'\code{check.survMCOD}, which enables the user to check convergence (type ?check.survMCOD for details).
#'
#'@references
#'
#'Moreno-Betancur M, Sadaoui H, Piffaretti C, Rey G. Survival analysis with multiple causes of death:
#'Extending the competing risks model. Epidemiology 2017; 28(1): 12-19.
#'
#'Piffaretti C, Moreno-Betancur M, Lamarche-Vadel A, Rey G. Quantifying cause-related mortality by
#'weighting multiple causes of death. Bulletin of the World Health Organization 2016; 94:870-879B.
#'
#'Rey G, Piffaretti C, Rondet C, Lamarche-Vadel A, Moreno-Betancur M. Analyse de la mortalite par cause :
#'ponderation des causes multiples. Bulletin Epidemiologique Hebdomadaire, 2017; (1): 13-9.
#'
#'@examples
#'
#'   ## Example ##
#'
#'   # First we simulate data using the simMCOD function:
#'   datEx<-simMCOD(n=1000,xi=-1,rho=-2,phi=0,
#'          pgen=c(1,0,0.75,0.25,0.125,0.083),
#'          lambda=0.001,v=2,pUC=c(1,0.75))
#'
#'   # Run analysis
#'
#'   fitMCOD<-survMCOD(SurvM(time=TimeEntry,time2=TimeExit,status=Status,
#'                           weight=Pi)~X1,
#'                    formOther=~Z1,data=datEx,UC_indicator="UC")
#'
#'   # Multiple-cause analysis results
#'    fitMCOD[[1]]
#'
#'   # Single-cause analysis results
#'    fitMCOD[[2]]
#'
#'   # Check convergence of multiple-cause analysis
#'     check.survMCOD(fitMCOD)
#'
#' @import survival
#' @importFrom Matrix Matrix
#' @importFrom utils flush.console
#' @importFrom stats rbinom runif rmultinom terms model.matrix model.extract quantile as.formula optim pnorm rmultinom



survMCOD<-function(formula, formOther=formula[-2], data,
                   UC_indicator,Iter=4)
{
  ##################################### Preliminaries and checks #####################################

  call <- match.call()
  m <- match.call(expand.dots = F)
  Terms <- terms(formula, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m$formOther<-m$UC_indicator<-NULL
  m <- eval(m, parent.frame())
  S <- model.extract(m, "response")
  if (!inherits(S, "SurvM"))
    stop("Response must be a SurvM object")

  if(missing(UC_indicator))
    stop("Underlying cause indicator must be provided")
  # Recover names of all outcome variables in model
  tstart <- as.character(Terms[[2]][2])
  tstop <- as.character(Terms[[2]][3])
  sname <- as.character(Terms[[2]][4])
  wname <- as.character(Terms[[2]][5])

  # Create priliminary design matrix and count number of dummy covariates for each.
  matX<-as.data.frame(model.matrix(formula[-2],data)[,-1])
  namX<-dimnames(model.matrix(formula[-2],data))[[2]][-1]
  matZ<-as.data.frame(model.matrix(formOther,data)[,-1])
  namZ<-dimnames(model.matrix(formOther,data))[[2]][-1]


  nX<-ncol(matX)
  nZ<-ncol(matZ)
  if(nX==0|nZ==0)
    stop("Covariates must be included in models for both pure hazards")
  names(matX)<-paste("X",1:nX,sep="")
  names(matZ)<-paste("Z",1:nZ,sep="")


  ################################## PREPARE DATA AND ARGUMENTS FOR FUNCTIONS #####################################
  #### Preparation of "datInd" and "tt" to feed to functions.
  # This is to use age as time-scale
  dat<-cbind(matZ,matX,"AgeEntry"=data[, tstart],"AgeExit"=data[, tstop],"Uncens"=data[, sname])
  dat<-dat[order(dat$AgeExit),]
  datamat<-as.matrix(dat)
  datInd<-dat
  tt<-unique(datamat[datamat[,"Uncens"]!=0,"AgeExit"])
  IndEntry<-apply(datamat,1,getEntry,tt)
  IndExit<-apply(datamat,1,getExit,tt)
  datInd<-cbind(datInd,IndEntry,IndExit)
  datInd2<-apply(datInd,1,atrisk,tt)
  datInd<-cbind(datInd,t(datInd2))
  datInd<-datInd[,-(1:(ncol(datInd)-nrow(datInd2)))]
  datInd<-Matrix(as.matrix(datInd),sparse=T)
  n<-nrow(dat)

  #### Prepare dataset "datamat" to feed to functions
  # in particular, determine weights and new status variable
  dat2<-data
  dat2$pi<-data[, wname]
  dat2$AgeEntry<-dat2[, tstart]
  dat2$AgeExit<-dat2[, tstop]
  psetTT<-c(1,0,sort(setdiff(unique(dat2$pi),c(0,1)),decreasing = T))
  dat2$status<-apply(cbind(1:nrow(dat2),dat2$pi),1,function(x)return(ifelse(is.na(x[2]),0,which(round(psetTT,10)==round(x[2],10)))))
  pset<-psetTT
  statuses<-1:length(psetTT)
  dat<-dat2[,c(paste("Z",1:nZ,sep=""),paste("X",1:nX,sep=""),"AgeEntry","AgeExit","status")]
  dat<-dat[order(dat$AgeExit),]
  datamat<-as.matrix(dat)

  #### Determine age-groups with at least one "pure" event for piecewise constant estimation of baseline hazard
  #### This leads to txi and nxi to feed to functions.
  txi<-c(0,round(quantile(dat$AgeExit[dat$status==1])[2:4],2))
  nxi<-length(txi)
  namXi<-paste("xi.",txi,sep="")



  ################################## FIT MULTIPLE-CAUSE MODEL #####################################

  ################### Starting values

  #### rho and phi: Single-cause analysis: Underlying cause + two separate Cox models
  dat1<-dat2[,c(paste("Z",1:nZ,sep=""),paste("X",1:nX,sep=""),"AgeEntry","AgeExit",sname,UC_indicator)]
  dat1$NoUC<-ifelse(dat1[,sname]==0,0,ifelse(dat1[,UC_indicator]==1,0,1))
  dat1$UC<-dat1[,UC_indicator]

  fit<-summary(coxph(as.formula(paste("Surv(time=AgeEntry,time2=AgeExit,type=\"counting\",
                                      event=UC)~",paste("X",1:nX,sep="",collapse="+"))),data=dat1))
  if(nX==1)
  rhoSCOD<-c(fit$coef[,c(1,3)],log(fit$conf.int[,3:4]),fit$coef[,c(5)]) else
  rhoSCOD<-cbind(fit$coef[,c(1,3)],log(fit$conf.int[,3:4]),fit$coef[,c(5)])


  fit<-summary(coxph(as.formula(paste("Surv(time=AgeEntry,time2=AgeExit,type=\"counting\",
                                      event=NoUC)~",paste("Z",1:nZ,sep="",collapse="+"))), data=dat1))
  if(nZ==1)
  phiSCOD<-c(fit$coef[,c(1,3)],log(fit$conf.int[,3:4]),fit$coef[,c(5)]) else
  phiSCOD<-cbind(fit$coef[,c(1,3)],log(fit$conf.int[,3:4]),fit$coef[,c(5)])

  names(rhoSCOD)<-c("Coef","SE","CIupp","CIlow","pvalue")
  names(phiSCOD)<-c("Coef","SE","CIupp","CIlow","pvalue")

  resSCOD<-rbind(rhoSCOD,phiSCOD)

  #### xi:  Non-param estimator of xi
  xiNonpar<-as.data.frame(stat(dat,1:nrow(dat),nxi))

  ################### Optimization
  #loglik and the other likelihood functions use the following objects defined above: nxi, txi, nZ, nX, pset, pseTT, datamat, tt, datInd

  opt<-matrix(NA,ncol=nX+nZ+nxi,nrow=Iter+2) # opt contains (xi,rho,phi), in that order

  ## Get second starting values from L, using as first starting values the nonpar + classic analyses
  opt[1,]<-c(xiNonpar[,1],resSCOD[1:nX,1],resSCOD[(nX+1):(nX+nZ),1]) #first starting values
  opt[2,]<-optim(par=opt[1,],fn=loglik, gr = NULL,nxi,nX,nZ,pset,datamat,psetTT,txi,tt,datInd,method="BFGS")$par  # second (better) starting values

  ## Run iterative procedure
  for(j in 3:(Iter+2))
  {
    i<-(j-1)
    opt[j,1:nxi]<-optim(par=opt[i,1:nxi],fn=logliksteropt, gr = NULL,
                                      nxi,nX,nZ,pset,datamat,psetTT,txi,tt,datInd,opt,i,method = "BFGS")$par
    opt[j,(nxi+1):(nZ+nX+nxi)]<-optim(par=opt[i,(nxi+1):(nZ+nX+nxi)],fn=loglikopt, gr = NULL,
                                      nxi,nX,nZ,pset,datamat,psetTT,txi,tt,datInd,opt,j,method = "BFGS")$par
    print(paste("Iteration",j-2,"out of 4 completed",sep=" "))
    flush.console()
  }


  ################### Derive variance estimators
  H<-hess(opt[(Iter+2),],nxi,nX,nZ,pset,datamat,psetTT,txi,tt,datInd)
  if(inherits(try(solve(H),silent=T),"try-error"))
    SE<-NA else{
      D2=-H
      D2[1:nxi,(nxi+1):(nxi+nX+nZ)]<-D2[(nxi+1):(nxi+nX+nZ),1:nxi]
      V2<-solve(H)%*%D2%*%t(solve(H))
      SE<-sqrt(diag(V2))}


  tabres<-data.frame(Coef=opt[(Iter+2),],SE=SE,CIlow=opt[(Iter+2),]-1.96*SE,
                      CIupp=opt[(Iter+2),]+1.96*SE, "pvalue"=2*(1-pnorm(q=abs(opt[(Iter+2),]/SE))),
                      row.names=
                        c(namXi,namX,namZ))

  RESMCOD<-list()
  RESMCOD[["Disease of interest"]]<-  tabres[(nxi+1):(nxi+nX),]
  RESMCOD[["Other diseases"]]<-  tabres[(nxi+nX+1):(nxi+nX+nZ),]
  RESMCOD[["Piecewise constant log ratio of the baseline pure hazards"]]<-  tabres[1:nxi,]

  row.names(resSCOD)<-c(namX,namZ)
  resSCOD<-data.frame(resSCOD)
  RESSCOD<-list()
  RESSCOD[["Disease of interest"]]<-resSCOD[1:nX,]
  RESSCOD[["Other diseases"]]<-resSCOD[(nX+1):(nX+nZ),]

  RES<-list()
  RES[["Multiple-cause"]]<-RESMCOD
  RES[["Single-cause"]]<-RESSCOD
  opt<-data.frame(opt)
  names(opt)<-c(namXi,namX,namZ)
  RES[["opt"]]<-opt

  return(RES)

}
