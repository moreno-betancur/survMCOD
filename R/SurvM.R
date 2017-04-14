#' Create a multiple-cause survival object
#'
#' \code{SurvM} Creates a multiple-cause survival object, which is used as a response variable
#' in the model formula provided to \code{\link[survMCOD]{survMCOD}}.
#'
#' @export
#'
#'
#' @param time For right censored data, this is the follow up time. For counting process data,
#' the first argument is the entry time.
#' @param time2 For counting process data only, it is the exit time.
#' @param status The vital status indicator such that 0=alive and 1=dead.
#' @param weight The weight indicating the proportion of the death attributed to the disease of interest.
#' It should be NA if individual alive (status=0) and a number between 0 and 1 if dead (status=1).
#' @param type Character string specifying the type of censoring. Possible values are "right" and "counting".
#'
#' @details This function needs to be used to in the formula argument of the \code{\link[survMCOD]{survMCOD}} function.
#'
#' The user is referred to Moreno-Betancur et al. (2017), Piffaretti et al. (2016) and Rey et al. (2017)
#' for descriptions and discussions of various weight-attribution strategies to determine the \code{weight} argument.
#'
#' @return An object of class \code{SurvM}.
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
#'
#'@examples
#'
#'   datEx<-simMCOD(n=400,xi=-1,rho=-2,phi=0,pgen=c(1,0,0.75,0.25,0.125,0.083),
#'           lambda=0.001,v=2,pUC=c(1,0.75))
#'
#'   SurvM(time=datEx$TimeEntry, time2=datEx$TimeExit, status=datEx$Status,
#'   weight=datEx$Pi, type="counting")
#'
#'@importFrom utils flush.console
#'@importFrom stats rbinom runif rmultinom terms model.matrix model.extract quantile as.formula optim pnorm rmultinom

SurvM<-
function (time, time2, status, weight, type = c("right", "counting"))
{
  # Check time argument
  if (missing(time))
    stop("Must have a time argument")
  if (inherits(time, "difftime"))
    time <- unclass(time)
  if (!missing(time2) && class(time2) == "difftime")
    time2 <- as.numeric(time2)
  if (!is.numeric(time))
    stop("Time variable is not numeric")
  nn <- length(time)

  # Check status argument
  if (missing(status)) {
    stop("Must specify 'status' argument")
  }else{
  if(any(as.numeric(status)!=1&as.numeric(status)!=0))
    stop("Argument 'status' must be logical or numeric equal to 0 or 1")}
  if (length(status) != nn)
    stop("Time and status are different lengths")

  # Check weight argument
    if (!is.numeric(weight))
    stop("Argument 'weight' must be numeric")
  if (missing(weight)) {
    stop("Must specify 'weight' argument")
  } else
  if(any(weight>1|weight<0,na.rm=T))
      stop("Weights must be comprised between 0 and 1")
  if (length(weight) != nn)
    stop("Time and weight are different lengths")
  if (any(is.na(weight)&status==1))
    stop("A non-missing weight must be provided for all subjects with an event (status=1)")
  if (any(!is.na(weight)&status==0))
    stop("Non-missing weights for censored subjects (status=0) will be ignored")

  # Check type argument
  ng <- (!missing(time)) + (!missing(time2))
  mtype <- match.arg(type)

  if(!missing(type)&mtype!="right"&mtype!="counting")
    stop("Type must be 'right' or 'counting'")

  if (missing(type)) {
    if (ng == 1 )
      type <- "right"
    else if (ng == 2)
      type <- "counting"
    else stop("No time variable!")
  }else {
    type <- mtype
    if (ng < 2 && (type == "counting"))
      stop("Wrong number of args for this type of survival data")
    if (ng >1 && (type == "right"))
      stop("Wrong number of args for this type of survival data")
  }

 #collect key info in vector ss
 if (type == "right") {
    if (!is.numeric(time))
      stop("Time variable is not numeric")
    if (any(time<=0))
      stop("Time variable needs to be >0")
    ss <- cbind(start = rep(0,length(time)), stop = time, status = status, weight=weight)
  } else if (type == "counting") {
    if (length(time2) != nn)
      stop("Start and stop are different lengths")
    if (!is.numeric(time))
      stop("Start time is not numeric")
    if (!is.numeric(time2))
      stop("Stop time is not numeric")
    if (any(time2<=0))
      stop("Stop time needs to be >0")
    temp <- (time >= time2)
    if (any(temp & !is.na(temp))) {
      time[temp] <- NA
      warning("Stop time must be > start time, NA created")
    }

    ss <- cbind(start = time, stop = time2,
                status = status, weight=weight)
  }

  inputAttributes <- list()
  if (!is.null(attributes(time)))
    inputAttributes$time <- attributes(time)
  if (!missing(time2) && !is.null(attributes(time2)))
    inputAttributes$time2 <- attributes(time2)
  if (!missing(weight) && !is.null(attributes(weight)))
    inputAttributes$weight <- attributes(weight)
  cname <- dimnames(ss)[[2]]
  dimnames(ss) <- list(NULL, cname)
  attr(ss, "type") <- type
  if (length(inputAttributes) > 0)
    attr(ss, "inputAttributes") <- inputAttributes
  class(ss) <- "SurvM"
  ss
}
