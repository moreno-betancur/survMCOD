
#' Hessian matrix as a function of (xi,rho,phi)
#'
#' Internal function used in estimation procedure
#' @export
#' @keywords internal
#'

hess<-function(par,nxi,nX,nZ,pset,datamat,psetTT,txi,tt,datInd){
  #par contient (xi,rho,phi)

  xi=par[1:nxi]
  if(nX==0){rho=as.matrix(0)}else{rho=as.matrix(c(rep(0,nZ),par[(nxi+1):(nxi+nX)]))}
  if(nZ==0){phi=as.matrix(0)}else{phi=as.matrix(c(par[(nX+nxi+1):(nZ+nX+nxi)],rep(0,nX)))}

  Hessian<-matrix(0,nrow=nZ+nX+nxi,ncol=nZ+nX+nxi)
  Q<-length(pset)


  # This data frame has the "a_tot" term (denominator in l*).
  # Column j has the term corresponding to xi_j
  datXiT<-NULL

  # These data frames have the derivatives of the "a_tot" term (denominator in l*).
  # Column j has the term corresponding to xi_j
  daTrho1<-NULL
  daTrhorho1C<-NULL
  daTrhophiC<-NULL
  daTrhoxi<-NULL

  daTphi1<-NULL
  daTphiphi1C<-NULL
  daTphirhoC<-NULL
  daTphixi1<-NULL

  daTxi<-NULL
  daTxixi<-NULL
  daTxiphi1C<-NULL
  daTxirhoC<-NULL

  for(ff in 1:nxi)
  {

    datXiT<-cbind(datXiT,sum(pset)*exp(datamat[,1:(nZ+nX)]%*%rho)+(Q-sum(pset))*exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi))

    daTrho1<-cbind(daTrho1,sum(pset)*(as.vector(exp(datamat[,1:(nZ+nX)]%*%rho))))
    daTrhorho1C<-cbind(daTrhorho1C,sum(pset)*(as.vector(exp(datamat[,1:(nZ+nX)]%*%rho))))
    daTrhophiC<-cbind(daTrhophiC,rep(0,nrow(datamat)))
    daTrhoxi<-cbind(daTrhoxi,rep(0,nrow(datamat)))

    daTphi1<-cbind(daTphi1,(Q-sum(pset))*as.vector(exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi)))
    daTphiphi1C<-cbind(daTphiphi1C,(Q-sum(pset))*as.vector(exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi)))
    daTphirhoC<-cbind(daTphirhoC,rep(0,nrow(datamat)))
    daTphixi1<-cbind(daTphixi1,(Q-sum(pset))*as.vector(exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi))*(-1))

    daTxi<-cbind(daTxi,(Q-sum(pset))*as.vector(exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi))*(-1))
    daTxixi<-cbind(daTxixi,(Q-sum(pset))*as.vector(exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi))*(-1)*(-1))
    daTxiphi1C<-cbind(daTxiphi1C,(Q-sum(pset))*as.vector(exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi))*(-1))
    daTxirhoC<-cbind(daTxirhoC,rep(0,nrow(datamat)))
  }


  # Here we loop through the pi's to get the a_pi terms
  for(pp in pset)
  {
    q<-which(psetTT==pp)


    # This data frame has the "a_pi" term (denominator in l).
    # Column j has the term corresponding to xi_j
    datXi<-NULL

    # These data frames have the derivatives of the "a_pi" term (denominator in l).
    # Column j has the term corresponding to xi_j
    darho1<-NULL
    darhorho1C<-NULL
    darhophiC<-NULL
    darhoxi<-NULL

    daphi1<-NULL
    daphiphi1C<-NULL
    daphirhoC<-NULL
    daphixi1<-NULL

    daxi<-NULL
    daxixi<-NULL
    daxirhoC<-NULL
    daxiphi1C<-NULL

    for(ff in 1:nxi)
    {
      datXi<-cbind(datXi,pp*exp(datamat[,1:(nZ+nX)]%*%rho)+(1-pp)*exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi))

      darho1<-cbind(darho1,  pp*(as.vector(exp(datamat[,1:(nZ+nX)]%*%rho))))
      darhorho1C<-cbind(darhorho1C,pp*(as.vector(exp(datamat[,1:(nZ+nX)]%*%rho))))
      darhophiC<-cbind(darhophiC, rep(0,nrow(datamat)))
      darhoxi<-cbind(  darhoxi, rep(0,nrow(datamat)))

      daphi1<-cbind(daphi1,(1-pp)*as.vector(exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi)))
      daphiphi1C<-cbind(daphiphi1C,(1-pp)*as.vector(exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi)))
      daphirhoC<-cbind( daphirhoC, rep(0,nrow(datamat)))
      daphixi1<-cbind(daphixi1,(1-pp)*as.vector(exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi))*(-1))

      daxi<-cbind(daxi,(1-pp)*as.vector(exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi))*(-1))
      daxixi<-cbind(daxixi,(1-pp)*as.vector(exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi))*(-1)*(-1))
      daxirhoC<-cbind(daxirhoC,rep(0,nrow(datamat)))
      daxiphi1C<-cbind(daxiphi1C,(1-pp)*as.vector(exp(-xi[ff]+datamat[,1:(nZ+nX)]%*%phi))*(-1))
    }

    # get power of denominator: number of ties for each uniqute type-q event time
    ttq<-datamat[datamat[,"status"]==q,"AgeExit"]
    tt2<-table(ttq)
    tt2<-data.frame(ID=1:nrow(tt2),ttq=as.numeric(dimnames(tt2)$ttq),tt2,row.names=NULL)[,c(1,2,4)]



    for(rrow in 1:(nxi+nX+nZ)) {
      for(ccol in 1:(nZ+nX+nxi)) {


        #In this part: finish calculating the derivaties of both a_pi and a_tot terms,
        # which requires multiplication by the values of the covariates.

        if(rrow-nxi>0 & rrow-nxi-nX<=0) dxr<-datamat[,nZ+rrow-nxi] else dxr<-rep(NA,nrow(datamat))
        if(ccol-nxi>0 & ccol-nxi-nX<=0) dxc<-datamat[,nZ+ccol-nxi] else dxc<-rep(NA,nrow(datamat))
        if(rrow-nxi-nX>0) dzr<-datamat[,rrow-nxi-nX] else dzr<-rep(NA,nrow(datamat))
        if(ccol-nxi-nX>0) dzc<-datamat[,ccol-nxi-nX] else dzc<-rep(NA,nrow(datamat))

        darho<- darho1*dxr
        daTrho<-daTrho1*dxr
        darhoC<- darho1*dxc
        daTrhoC<-daTrho1*dxc

        darhorhoC<-darhorho1C*dxr*dxc
        daTrhorhoC<-daTrhorho1C*dxr*dxc

        daphi<-daphi1*dzr
        daTphi<-daTphi1*dzr
        daphiC<-daphi1*dzc
        daTphiC<-daTphi1*dzc

        daphiphiC<-daphiphi1C*dzr*dzc
        daTphiphiC<-daTphiphi1C*dzr*dzc
        daphixi<-daphixi1*dzr
        daTphixi<-daTphixi1*dzr

        daxiphiC<-daxiphi1C*dzc
        daTxiphiC<-daTxiphi1C*dzc


        if(rrow>nxi)
        {         nn<-1*(rrow<=(nxi+nX))+2*(rrow>(nxi+nX))
        mm<-1*(ccol<=nxi)+2*(nxi<ccol&ccol<=(nxi+nX))+3*(ccol>nxi+nX)
        theta1<-c("rho", "phi")[nn]
        theta2<-c("xi","rhoC","phiC")[mm]

        term1<-0
        term2<-0

        if(theta2=="xi") loopin<-ccol else loopin<-1:nxi  #only individuals with event in that frame contribute here

        for(ff in loopin)
        {
          dp<-(datamat[,"status"]==q&txi[ff]<=datamat[,"AgeExit"]&datamat[,"AgeExit"]<ifelse(is.na(txi[ff+1]),Inf,txi[ff+1]))

          term1<-term1+sum((datXi[dp,ff]*get(paste("da",theta1,theta2,sep=""))[dp,ff]-
                              get(paste("da",theta1,sep=""))[dp,ff]*get(paste("da",theta2,sep=""))[dp,ff])/(datXi[dp,ff]^2))

          tind<-(tt%in%ttq)&(txi[ff]<=tt)&(tt<ifelse(is.na(txi[ff+1]),Inf,txi[ff+1]))
          datq<-datInd[,tind]

          term2<-term2+sum(
            tt2$Freq[txi[ff]<=tt2$ttq&tt2$ttq<ifelse(is.na(txi[ff+1]),Inf,txi[ff+1])]*
              (((t(datXi[,ff])%*%datq)*(t(get(paste("da",theta1,theta2,sep=""))[,ff])%*%datq)-
                  (t(get(paste("da",theta1,sep=""))[,ff])%*%datq)*
                  (t(get(paste("da",theta2,sep=""))[,ff])%*%datq))/((t(datXi[,ff])%*%datq)^2)))

        }

        contr<-term1-term2

        Hessian[rrow ,ccol]<-Hessian[rrow,ccol]+ contr
        }else{

          mm<-1*(ccol<=nxi)+2*(nxi<ccol&ccol<=(nxi+nX))+3*(ccol>nxi+nX)
          theta1<-"xi"
          theta2<-c("xi","rhoC","phiC")[mm]

          term1<-0
          term2<-0

          if(theta2!="xi"|ccol==rrow)
          {
            ff<-(1:nxi)[rrow]  #only events occurring in that time frame contribute to that row

            dp<-(datamat[,"status"]==q&txi[ff]<=datamat[,"AgeExit"]&datamat[,"AgeExit"]<ifelse(is.na(txi[ff+1]),Inf,txi[ff+1]))

            term1<-term1+sum((datXi[dp,ff]*get(paste("da",theta1,theta2,sep=""))[dp,ff]-
                                get(paste("da",theta1,sep=""))[dp,ff]*get(paste("da",theta2,sep=""))[dp,ff])/(datXi[dp,ff]^2))

            tind<-(tt%in%ttq)&(txi[ff]<=tt)&(tt<ifelse(is.na(txi[ff+1]),Inf,txi[ff+1]))
            datq<-datInd[,tind]

            term2<-term2+ sum(
              tt2$Freq[txi[ff]<=tt2$ttq&tt2$ttq<ifelse(is.na(txi[ff+1]),Inf,txi[ff+1])]*
                (((t(datXiT[,ff])%*%datq)*(t(get(paste("daT",theta1,theta2,sep=""))[,ff])%*%datq)-
                    (t(get(paste("daT",theta1,sep=""))[,ff])%*%datq)*
                    (t(get(paste("daT",theta2,sep=""))[,ff])%*%datq))/((t(datXiT[,ff])%*%datq)^2)))

          }
          contr<-term1-term2

          Hessian[rrow ,ccol]<-Hessian[rrow,ccol]+ contr
          #print(Hessian[1,1])
          #flush.console()
        }
      }
    }}
  return(Hessian)
}

