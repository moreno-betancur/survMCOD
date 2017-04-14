#' Produce plots to check convergence of multiple-cause model
#'
#' @export
#' @param fit A model fit returned by \code{survMCOD}.
#' @details
#'The plots produced show Iter+2 estimates of each parameter, where Iter is the number of iterations specified by the user in the call
#'to \code{survMCOD}. The first two estimates correspond to a starting value and an improved starting value
#'(see Moreno-Betancur et al. 2017 for details). Healthy convergence is seen by curves showing variation in estimates across the first
#'three points, followed by a stabilisation of the curve around the final estimate.
#'
#'
#'@references
#'
#'Moreno-Betancur M, Sadaoui H, Piffaretti C, Rey G. Survival analysis with multiple causes of death:
#'Extending the competing risks model. Epidemiology 2017; 28(1): 12-19.
#'
#'
#'@examples
#'
#'   ## Example ## uncomment to run
#'
#'   # First we simulate data using the simMCOD function:
#'   # datEx<-simMCOD(n=1000,xi=-1,rho=-2,phi=0,
#'   #       pgen=c(1,0,0.75,0.25,0.125,0.083),
#'   #      lambda=0.001,v=2,pUC=c(1,0.75))
#'
#'   # Run analysis
#'
#'   # fitMCOD<-survMCOD(SurvM(time=TimeEntry,time2=TimeExit,status=Status,
#'   #                       weight=Pi)~X1,
#'   #                  formOther=~Z1,data=datEx,UC_indicator="UC")
#'
#'   # Check convergence of multiple-cause analysis
#'   # check.survMCOD(fitMCOD)
#'
#'@importFrom graphics par plot title

check.survMCOD<-function(fit)
{
  opt<-fit$opt
  nc<-ncol(opt)
  cc<-floor(sqrt(nc))
  if((cc*(cc+1))>=nc) siz<-c(cc,cc+1) else
    siz<-c(cc,cc+2)

  par(mfrow=siz,oma=c(0, 0, 3, 0))
  for(xx in 1:nc)
     plot(1:nrow(opt),opt[,xx],type="b",main=names(opt)[xx],xlab="Iteration",ylab="Estimate")

  title("Check convergence of estimates of each parameter",outer=T,line=1)
}
