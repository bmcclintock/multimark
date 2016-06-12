#' Simulate spatially-explicit capture-mark-recapture data from a (demographically) closed population with multiple non-invasive marks
#'
#' This function generates encounter histories from spatially-explicit capture-mark-recapture data consisting of multiple non-invasive marks. 
#'
#'
#' @param N True population size or abundance.
#' @param ntraps The number of traps. If \code{trapCoords=NULL}, the square root of \code{ntraps} must be a whole number in order to create a regular grid of trap coordinates on the unit square.
#' @param noccas Scaler indicating the number of sampling occasions per trap.
#' @param pbeta Complementary loglog-scale intercept term for detection probability (p). Must be a scaler or vector of length \code{noccas}.
#' @param tau Additive complementary loglog-scale behavioral effect term for recapture probability (c).
#' @param sigma2_scr Complementary loglog-scale variance for the detection function.
#' @param delta_1 Conditional probability of type 1 encounter, given detection.
#' @param delta_2 Conditional probability of type 2 encounter, given detection.
#' @param alpha Conditional probability of simultaneous type 1 and type 2 detection, given both types encountered. Only applies when \code{data.type="sometimes"}.
#' @param data.type Specifies the encounter history data type. All data types include non-detections (type 0 encounter), type 1 encounter (e.g., left-side), and type 2 encounters (e.g., right-side). When both type 1 and type 2 encounters occur for the same individual within a sampling occasion, these can either be "non-simultaneous" (type 3 encounter) or "simultaneous" (type 4 encounter). Three data types are currently permitted:
#' 
#'  \code{data.type="never"} indicates both type 1 and type 2 encounters are never observed for the same individual within a sampling occasion, and observed encounter histories therefore include only type 1 or type 2 encounters (e.g., only left- and right-sided photographs were collected). Observed encounter histories can consist of non-detections (0), type 1 encounters (1), and type 2 encounters (2). See \code{\link{bobcat}}. Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 3 encounters (3).
#'
#'  \code{data.type="sometimes"} indicates both type 1 and type 2 encounters are sometimes observed (e.g., both-sided photographs are sometimes obtained, but not necessarily for all individuals). Observed encounter histories can consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), type 3 encounters (3), and type 4 encounters (4). Type 3 encounters can only be observed when an individual has at least one type 4 encounter. Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), type 3 encounters (3), and type 4 encounters (4). 
#'
#'  \code{data.type="always"} indicates both type 1 and type 2 encounters are always observed, but some encounter histories may still include only type 1 or type 2 encounters. Observed encounter histories can consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 4 encounters (4). Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 4 encounters (4).
#'
#' @param detection Detection function for detection probability. Must be "\code{half-normal}" or "\code{exponential}".
#' @param spatialInputs A list of length 3 composed of objects named \code{trapCoords}, \code{studyArea}, and \code{centers}:
#' 
#'  \code{trapCoords} is a matrix of dimension \code{ntraps} x 2 indicating the Cartesian coordinates for the traps, where rows correspond to trap, the first column the x-coordinate, and the second column the y-coordinate. If \code{spatialInputs=NULL} (the default), the traps are placed in a regular grid on the unit square.
#'
#'  \code{studyArea} is a 3-column matrix defining the study area and available habitat. Each row corresponds to a grid cell. The first 2 columns indicate the Cartesian x- and y-coordinate for the centroid of each grid cell, and the third column indicates whether the cell is available habitat (=1) or not (=0). All cells must have the same resolution. If \code{spatialInputs=NULL} (the default), the study area is assumed to be composed of 4000 grid cells of available habitat over the unit square plus a buffer of 3*\code{sqrt(sigma2_scr)}.
#'   
#'  \code{centers} is a matrix containing the true (latent) coordinates of the activity centers for each individual in the population, where each row corresponds to an individual, the first column the x-coordinate, and the second column the y-coordinate. If \code{spatialInputs=NULL} (the default), the activity centers are randomly placed on the study area.
#'
#' @details Please be very careful when specifying your own \code{spatialInputs}; \code{\link{multimarkClosedSCR}} does not verify that these make sense during model fitting.  
#'
#' @return A list containing the following:
#' \item{Enc.Mat}{A matrix containing the observed encounter histories with rows corresponding to individuals and (\code{ntraps}*\code{noccas}) columns corresponding to traps and sampling occasions.  The first \code{noccas} columns correspond to trap 1, the second \code{noccas} columns corresopond to trap 2, etc.}
#' \item{trueEnc.Mat}{A matrix containing the true (latent) encounter histories with rows corresponding to individuals and (\code{ntraps}*\code{noccas}) columns corresponding to traps and sampling occasions.  The first \code{noccas} columns correspond to trap 1, the second \code{noccas} columns corresopond to trap 2, etc.}
#' \item{spatialInputs}{}
#' @author Brett T. McClintock 
#' @seealso \code{\link{processdata}}, \code{\link{multimarkClosedSCR}}
#' @references
#' Bonner, S. J., and Holmberg J. 2013. Mark-recapture with multiple, non-invasive marks. \emph{Biometrics} 69: 766-775.
#' 
#' King. R., B. T. McClintock, D. Kidney, and D. Borchers. 2016. Capture-recapture abundance estimation using a semi-complete data likelihood approach. The Annals of Applied Statistics 10: 264-285.
#' 
#' McClintock, B. T., Conn, P. B., Alonso, R. S., and Crooks, K. R. 2013. Integrated modeling of bilateral photo-identification data in mark-recapture analyses. \emph{Ecology} 94: 1464-1471.
#' @examples
#' #simulate data for data.type="sometimes" using defaults
#' data<-simdataClosedSCR(data.type="sometimes")
simdataClosedSCR <- function(N=100,ntraps=16,noccas=5,pbeta=log(0.4),tau=0,sigma2_scr=0.1,delta_1=0.4,delta_2=0.4,alpha=0.5,data.type="never",detection="half-normal",spatialInputs=NULL){
  
  if(length(pbeta)==1){
    pbeta=rep(pbeta,noccas)
  } else if(length(pbeta)!=noccas){
    stop(paste0("'pbeta' must be a scaler or vector of length ",noccas))
  }
  delta_B<-1-(delta_1+delta_2)
  if(delta_B<0) stop ("delta_1 and delta_2 must have sum less than 1")
  
  if(data.type=="never"){
    alpha<-0
  } else if(data.type=="always"){
    alpha<-1
  } else if(data.type!="sometimes"){
    stop("'data.type' must be 'never', 'sometimes', or 'always'")
  }
  
  if(is.null(spatialInputs)){
    if(sqrt(ntraps)%%1) stop("The square root of 'ntraps' must be a whole number")
    spatialInputs=list()
    spatialInputs$trapCoords<-as.matrix(expand.grid(seq(0,1,length=sqrt(ntraps)),seq(0,1,length=sqrt(ntraps)))) #trap coordinates on unit square
    buffer <- 3*sqrt(sigma2_scr)
    #spatialInputs$centers<-cbind(runif(N,0-buffer,1+buffer),runif(N,0-buffer,1+buffer))                    #activity center coordinates on unit square plus buffer
    studyArea<-as.matrix(expand.grid(seq(0-buffer,1+buffer,length=sqrt(3600)),seq(0-buffer,1+buffer,length=sqrt(3600)))) #study area grid
    spatialInputs$studyArea<-cbind(studyArea,rep(1,3600))
    spatialInputs$centers<-sample.int(nrow(studyArea),N,replace=TRUE)
  } else {
    if(!is.list(spatialInputs) | length(spatialInputs)!=3 | any(sort(names(spatialInputs))!=c("centers","studyArea","trapCoords"))) stop("'spatialInputs' must be a list of length 3 containing the object 'trapCoords', 'studyArea', and 'centers'")
    if(ntraps!=nrow(spatialInputs$trapCoords)){
      warning("'ntraps' not equal to the number of rows of spatialInputs$trapCoords; the latter will be used.")
      ntraps<-nrow(spatialInputs$trapCoords)
    }
    if(N!=length(spatialInputs$centers)){
      warning("'N' not equal to the length of spatialInputs$centers; the latter will be used.")
      N<-length(spatialInputs$centers)
    }
  }
  colnames(spatialInputs$trapCoords)<-c("x","y")
  rownames(spatialInputs$trapCoords)<-paste0("trap",1:ntraps)
  #colnames(spatialInputs$centers)<-c("x","y")
  #rownames(spatialInputs$centers)<-paste0("center",1:N)
  colnames(spatialInputs$studyArea)<-c("x","y","avail")
  rownames(spatialInputs$studyArea)<-paste0("cell",1:nrow(spatialInputs$studyArea))
  trapCoords<-spatialInputs$trapCoords
  centers<-spatialInputs$studyArea[spatialInputs$centers,] 
  
  tEnc.Mat<-matrix(0,nrow=N,ncol=noccas*ntraps)        #"true" latent histories
  
  dist2<-getdist(centers,trapCoords)
  
  if(detection=="half-normal") {
    dexp<-2
  } else if(detection=="exponential"){
    dexp<-1
  } else {
    stop("'detection' argument must be 'half-normal' or 'exponential'")
  }
  
  for(i in 1:N){
    for(k in 1:ntraps){
      ind<-0
      for(j in 1:noccas){
        p<-invcloglog(pbeta[j]-1./(2*sigma2_scr)*dist2[i,k]^dexp)
        c<-invcloglog(pbeta[j]+tau-1./(2*sigma2_scr)*dist2[i,k]^dexp)
        tEnc.Mat[i,(k-1)*noccas+j] <- rbinom(1,1,((1-ind)*p+ind*c) )       #"true" latent histories
        if(tEnc.Mat[i,(k-1)*noccas+j]==1){
          ind<-1
        }
      }
    }
  }
  Rand.Mat<-matrix(runif(N*noccas*ntraps,0,1),N,noccas*ntraps)
  tEnc.Mat[which(tEnc.Mat==1 & Rand.Mat<delta_2)] <- 2      # type 2 encounters
  tEnc.Mat[which(tEnc.Mat==1 & Rand.Mat>(1-delta_B))] <- 4  # type 1 and type 2 encounters
  tEnc.Mat[which(tEnc.Mat==4)] <- tEnc.Mat[which(tEnc.Mat==4)]-(runif(base::sum(tEnc.Mat==4))<(1-alpha))   # unobserved type 1 and type 2 encounters
  
  Enc.Mat <- get_Enc(tEnc.Mat,data.type)
  return(list(Enc.Mat=Enc.Mat,trueEnc.Mat=tEnc.Mat,spatialInputs=spatialInputs))
}

pstarintegrandSCR<-function(noccas,beta,sigma2,DM,spatialInputs,dexp){
  
  XB <- DM %*% beta
  dist2 <- spatialInputs$dist2
  detProb<-invcloglogtol(matrix(XB,nrow=dim(dist2)[1],ncol=dim(dist2)[2]*noccas,byrow=TRUE)-1/(2*sigma2)*dist2[,rep(1:dim(dist2)[2],each=noccas)]^dexp)
  pdot<-1.-apply(1.-detProb,1,function(x) max(prod(x),tol))
  esa<-sum(pdot)*spatialInputs$a
  esa/spatialInputs$A
}

loglikeClosedSCR<-function(parms,DM,noccas,ntraps,C,All.hists,spatialInputs){
  
  H <- parms$H
  pbeta <- parms$pbeta
  if(DM$mod.p.h){
    zp <- parms$zp
  } else {
    zp <- rep(0,length(H))
  }
  if(DM$mod.delta != ~NULL){
    if(DM$mod.delta==formula(~type)){
      delta_1 <- parms$delta_1
      delta_2 <- parms$delta_2
    } else {
      delta_1 <- delta_2 <- parms$delta
    }
    alpha <- parms$alpha
  } else {
    delta_1 <- 1.0
    delta_2 <- 0.0
    alpha <- 0.0
  }
  
  dexp<-ifelse(DM$mod.det=="half-normal",2,1)
  
  Hind <- H[which(H>1)]
  centers <- spatialInputs$studyArea[parms$centers[which(H>1)],]
  indhist <- All.hists[Hind,]
  n<-length(Hind)
  #firstcap<- (C[Hind]>=matrix(rep(1:noccas,each=n),nrow=n,ncol=noccas))

  msk <- matrix(1,nrow=noccas,ncol=ntraps) #currently assumes all traps are active on all occasions
  msk2 <- array(NA, c(n, noccas, ntraps))
  for(i in 1:n){
    msk2[i, 1:noccas, 1:ntraps] <- msk[1:noccas, 1:ntraps]
  }
  msk2 <- as.vector(msk2)
  
  Yaug <- array(0, dim=c(n, noccas, ntraps))
  for(j in 1:n){
    Yaug[j, 1:noccas, 1:ntraps] <- matrix(indhist[j,],nrow=noccas,ncol=ntraps)#byrow=TRUE))#Y[j, 1:nT, 1:ntraps]
  }  
  # create covariate of previous capture
  prevcap <- array(0,c(n,noccas,ntraps))
  #prevcap2 <- matrix(0,nrow=n,ncol=ntraps*noccas)
  #dist2mat <- matrix(0,nrow=n,ncol=ntraps*noccas)
  #dist2<-getdist(centers,spatialInputs$trapCoords)
  for(i in 1:n){
    for(j in 1:ntraps){
      tmp <- Yaug[i, 1:noccas, j]
      #tmp2 <- indhist[i,(j-1)*noccas+1:noccas]
      if(any(tmp > 0)){
        fst <- min( (1:noccas)[tmp > 0] )
        if(fst < noccas)
          prevcap[i, (fst+1):noccas, j] <- 1
      }
      #if(any(tmp2 > 0)){
      #  fst <- min( (1:noccas)[tmp2 > 0] )
      #  if(fst < noccas)
      #    prevcap2[i, (j-1)*noccas+(fst+1):noccas] <- 1
      #}
      #dist2mat[i,(j-1)*noccas+1:noccas]<-dist2[i,j]
    }
  }
  prevcap <- as.vector(prevcap)
  
  ## vectorize all the data objects
  arr.trues <- array(TRUE, c(n,noccas,ntraps))
  idx <- which(arr.trues, arr.ind = TRUE)
  y <- as.vector(Yaug)
  y <- y[msk2==1]
  prevcap <- prevcap[msk2==1]   #### + 1   ### add 1 if want 1/2 otherwise dummy
  indid <- idx[msk2==1,1] ## [AMG] Individual IDs
  repid <- idx[msk2==1,2] ## [AMG] Replicate/Sampling Occasion IDs
  trapid <- idx[msk2==1,3] ## [AMG] Trap IDs
  
  trapgridbig <- spatialInputs$trapCoords[trapid,]   # stretches out the trap coord matrix
  c1 <- (centers[indid,1] - trapgridbig[,1])^2
  c2 <- (centers[indid,2] - trapgridbig[,2])^2
  
  p <- invcloglogtol(rep(DM$p%*%pbeta,each=n)*(1-prevcap) + rep(DM$c%*%pbeta,each=n)*prevcap - 1./(2*parms$sigma2_scr)*sqrt(c1+c2)^dexp)
  #p2 <- invcloglogtol(matrix(rep(DM$p%*%pbeta,each=n)*(1-prevcap2)+rep(DM$c%*%pbeta,each=n)*prevcap2,nrow=n,ncol=noccas*ntraps)- 1./(2*parms$sigma2_scr)*(dist2mat)^dexp)
  
  loglike <- base::sum( log( (y==0) * (1. - p)
                             + (y==1) * p * delta_1  
                             + (y==2) * p * delta_2
                             + (y==3) * p * (1. - delta_1 - delta_2) * (1. - alpha)
                             + (y==4) * p * (1. - delta_1 - delta_2) * alpha ))
  
  #loglike2 <- base::sum( log( (indhist==0) * (1. - p2)
  #                           + (indhist==1) * p2 * delta_1  
  #                           + (indhist==2) * p2 * delta_2
  #                           + (indhist==3) * p2 * (1. - delta_1 - delta_2) * (1. - alpha)
  #                           + (indhist==4) * p2 * (1. - delta_1 - delta_2) * alpha ))
  
  #if(DM$mod.p.h){
    pstar <- pstarintegrandSCR(noccas,pbeta,parms$sigma2_scr,DM$p,spatialInputs,dexp)
  #} else {
  #  pstar <- 1-min(1.-tol,max(tol,prod(1-expit(DM$p%*%pbeta))))
  #}    
  loglike <- loglike + dbinom(n,parms$N,pstar,1) - n * log(pstar)  
  loglike
}

priorsClosedSCR<-function(parms,DM,priorparms,data_type,spatialInputs){
  
  priors <- (base::sum(dnorm(parms$pbeta,priorparms$mu0,sqrt(priorparms$sigma2_mu0),log=TRUE))
             + -log(parms$N))
  
  if(DM$mod.delta != ~NULL){
    if(DM$mod.delta==formula(~type)){
      priors <- priors + ddirichlet(c(parms$delta_1,parms$delta_2,1.-parms$delta_1-parms$delta_2),priorparms$a0delta)
    } else {
      priors <- priors + dbeta(2*parms$delta,priorparms$a0delta[1],priorparms$a0delta[2],log=TRUE)
    }
    if(data_type=="sometimes"){
      priors <- priors + dbeta(parms$alpha,priorparms$a0alpha,priorparms$b0alpha,log=TRUE)
    }
    priors <- priors + (base::sum(dbinom((parms$H>1),1,parms$psi,log=TRUE))
                         + dbeta(parms$psi,priorparms$a0psi,priorparms$b0psi,log=TRUE))
  }
  
  priors <- priors + log(2.0*dcauchy(sqrt(parms$sigma2_scr),0.0,priorparms$a,log=FALSE))
  
  priors <- priors + length(parms$centers)*log(1./spatialInputs$A)
  
  if(DM$mod.p.h){
    priors <- priors + (base::sum(dnorm(parms$zp,0.0,sqrt(parms$sigma2_zp),log=TRUE))
                        + log(2.0*dcauchy(sqrt(parms$sigma2_zp),0.0,priorparms$a,log=FALSE)))
  }        
  priors
}

posteriorClosedSCR<-function(parms,DM,mms,priorparms,spatialInputs){
  nchains<-length(parms)
  ntraps<-nrow(spatialInputs$trapCoords)
  noccas<-ncol(mms@Enc.Mat)/ntraps
  M<-nrow(mms@Enc.Mat)
  All.hists<-matrix(mms@vAll.hists,byrow=TRUE,ncol=noccas*ntraps)
  for(ichain in 1:nchains){
    temp<-parms[[ichain]]
    
    loglike <- loglikeClosedSCR(temp,DM,noccas,ntraps,mms@C,All.hists,spatialInputs)
    
    if(!is.finite(loglike)) {
      stop(paste0("initial model likelihood is ",loglike," for chain ",ichain,". Try different initial values."))
    }
    
    posterior <- loglike + priorsClosedSCR(temp,DM,priorparms,mms@data.type,spatialInputs)
    
    if(!is.finite(posterior)) {
      stop(paste("initial model posterior is",posterior,"for chain",ichain,". Try different initial values or prior parameters"))
    }
  }
}

checkClosedSCR<-function(parms,parmlist,mms,DM,iter,adapt,bin,thin,burnin,taccept,tuneadjust,maxnumbasis,a0delta,a0alpha,b0alpha,a,sigma2_mu0,a0psi,b0psi){
  
  if(mms@data.type!="sometimes" & any(parms=="alpha")) stop("Parameter 'alpha' only applies to models for the 'sometimes' data type")
  
  params<-parms
  if(any(parms=="all")){
    if(mms@data.type=="sometimes"){
      params<-parmlist
    } else {
      params<-parmlist[which(parmlist!="alpha")]
    }
  } else {
    if(!all(match(params,parmlist,nomatch=0))) stop(paste0("'",params[match(params,parmlist,nomatch=0)==0],"' is not a valid parameter\n  "))
  }  
  
  if(adapt<0) stop("'adapt' must be >=0") 
  if((bin<1 | bin>iter) & iter>0) stop("'bin' must be >0 and <",iter)
  if(thin>max(1,floor((iter-burnin+1)/2)) | thin<1) stop("'thin' must be >0 and <=",max(1,floor((iter-burnin+1)/2)))
  if(taccept<=0 | taccept>1) stop ("'taccept' must be >0 and <=1")
  if(tuneadjust<=0 | tuneadjust>1) stop ("'tuneadjust' must be >0 and <=1")
  if(mms@ncolbasis & (maxnumbasis<1 | maxnumbasis>mms@ncolbasis)) stop("'maxnumbasis' must be between 1 and ",mms@ncolbasis)
  if(!all(c(a0delta,a0alpha,b0alpha,a,sigma2_mu0,a0psi,b0psi)>0)) stop("'a0delta', 'a0alpha', 'b0alpha', 'a', 'sigma2_mu0', 'a0psi', and 'b0psi' must be >0")
  
  pdim<-ncol(DM$p)
  if(!pdim) stop("'mod.p' must include at least 1 parameter")
  
  params
}

mcmcClosedSCR<-function(ichain,mms,DM,params,inits,iter,adapt,bin,thin,burnin,taccept,tuneadjust,Prop.sd,Prop.center,spatialInputs,maxnumbasis,a0delta,a0alpha,b0alpha,a,mu0,sigma2_mu0,a0psi,b0psi,printlog){
  
  gq<-gauss.quad(500,kind="hermite")
  weights<-gq$weights
  nodes<-gq$nodes
  npoints<-length(weights)
  
  ntraps<-nrow(spatialInputs$trapCoords)
  noccas<-ncol(mms@Enc.Mat)/ntraps
  M<-nrow(mms@Enc.Mat)
  DMp<-DM$p
  DMc<-DM$c
  mod.p.h<-DM$mod.p.h
  pdim<-ncol(DMp)
  dexp<-ifelse(DM$mod.det=="half-normal",2,1)
  firstcap<-get_C(matrix(mms@vAll.hists,byrow=TRUE,ncol=noccas),type="SCR")
  
  #declare and initialize parameters
  pbeta<-rep(NA,max(1,floor(iter/thin))*(pdim))
  zp<-rep(NA,ifelse(any(params=="zp"),max(1,floor(iter/thin))*M,M))
  H<-rep(NA,ifelse(any(params=="H"),max(1,floor(iter/thin))*M,M))
  sigma2_zp<-rep(NA,max(1,floor(iter/thin)))
  centers<-rep(NA,ifelse(any(params=="centers"),max(1,floor(iter/thin))*M,M)) 
  sigma2_scr<-rep(NA,max(1,floor(iter/thin)))
  alpha<-rep(NA,max(1,floor(iter/thin)))
  delta_1<-rep(NA,max(1,floor(iter/thin)))
  delta_2<-rep(NA,max(1,floor(iter/thin)))
  N<-rep(NA,max(1,floor(iter/thin)))
  psi<-rep(NA,max(1,floor(iter/thin)))
  logPosterior<-rep(NA,max(1,floor(iter/thin)))
  
  pbeta[1:pdim] <- inits[[ichain]]$pbeta
  zp[1:M] <- inits[[ichain]]$zp
  H[1:M] <- inits[[ichain]]$H-1
  sigma2_zp[1] <- inits[[ichain]]$sigma2_zp
  centers[1:M] <- inits[[ichain]]$centers-1
  sigma2_scr[1] <- inits[[ichain]]$sigma2_scr
  alpha[1] <- inits[[ichain]]$alpha
  delta_1[1] <- inits[[ichain]]$delta_1
  delta_2[1] <- inits[[ichain]]$delta_2
  N[1] <- inits[[ichain]]$N
  psi[1] <- inits[[ichain]]$psi
  
  arate<-numeric(M+pdim+1)
  
  posterior <- .C(ClosedSCRC,as.integer(ichain),as.numeric(mu0), as.numeric(sigma2_mu0), as.numeric(pbeta), as.numeric(zp), as.numeric(sigma2_zp), as.numeric(sigma2_scr), as.numeric(delta_1),as.numeric(delta_2),as.numeric(alpha), as.integer(inits[[ichain]]$x), as.numeric(N), as.numeric(psi), as.integer(H), as.integer(centers),
                  as.integer(ntraps),as.integer(noccas), as.integer(M), as.numeric(a0delta), as.numeric(a0alpha), as.numeric(b0alpha), as.numeric(a), as.numeric(a0psi), as.numeric(b0psi),
                  as.numeric(Prop.sd),as.integer(Prop.center$NNvect),as.integer(Prop.center$numnn),as.numeric(arate),as.numeric(logPosterior),
                  as.integer(length(mms@vAll.hists)/(noccas*ntraps)),as.integer(mms@vAll.hists), as.integer(firstcap), as.integer(mms@indBasis-1), as.integer(mms@ncolbasis), as.integer(mms@knownx), as.numeric(as.vector(t(DMp))), as.numeric(as.vector(t(DMc))),as.integer(pdim),
                  as.integer(iter), as.integer(thin), as.integer(adapt), as.integer(bin), as.numeric(taccept),as.numeric(tuneadjust),as.integer(maxnumbasis),
                  as.integer(npoints),as.numeric(weights),as.numeric(nodes),as.integer(mod.p.h),as.integer(mms@data.type=="sometimes"),as.integer(any(params=="zp")),as.integer(any(params=="H")),as.integer(any(params=="centers")),as.integer(DM$mod.delta != ~NULL),as.integer(DM$mod.delta==formula(~type)),as.numeric(dexp),as.numeric(spatialInputs$dist2),as.integer(nrow(spatialInputs$studyArea)),as.numeric(spatialInputs$A),as.integer(printlog),NAOK = TRUE) 
  
  names(posterior) <- c("ichain","mu_0","sigma2_mu","pbeta", "zp", "sigma2_zp", "sigma2_scr", "delta_1","delta_2","alpha", "x", "N", "psi","H", "centers", "ntraps", "noccas", "M","a0delta", "a0alpha", "b0alpha","a","a0psi","b0psi","Prop.sd", "NNvect", "numnn", "arate","logPosterior","nHists","vAll.hists","firstcap", "indBasis", "ncolBasis","knownx","DMp","DMc","pdim","iter", "thin", "adapt", "bin", "taccept","tuneadjust","maxnumbasis","npoints","weights","nodes","mod.p.h","sometimes?","zp?","H?","centers?","updatedelta?","type?","dexp","dist2","ncells","Area","printlog?")
  
  g <- posterior$iter
  x <- posterior$x
  if(any(params=="zp")){
    temp<-cbind(matrix(posterior$pbeta[(floor(burnin/thin)*pdim+1):(max(1,floor(iter/thin))*pdim)],ncol=pdim,byrow=T),matrix(posterior$zp[(floor(burnin/thin)*M+1):(max(1,floor(iter/thin))*M)],ncol=M,byrow=T),posterior$sigma2_scr[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$delta_1[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$delta_2[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$alpha[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$N[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$psi[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))]) 
    zp <- NULL
  } else {
    zp <- posterior$zp
    temp<-cbind(matrix(posterior$pbeta[(floor(burnin/thin)*pdim+1):(max(1,floor(iter/thin))*pdim)],ncol=pdim,byrow=T),posterior$sigma2_scr[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$delta_1[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$delta_2[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$alpha[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$N[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$psi[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))])       
  }
  if(any(params=="H")){
    posterior<-cbind(temp,matrix(posterior$H[(floor(burnin/thin)*M+1):(max(1,floor(iter/thin))*M)]+1,ncol=M,byrow=T),posterior$logPosterior[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))]) 
    H <- NULL
  } else {
    H <- posterior$H+1
    posterior<-cbind(temp,posterior$logPosterior[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))])       
  }
  return(list(posterior=posterior,x=x,H=H,zp=zp,g=g))
}

#' Fit spatially-explicit population abundance models for capture-mark-recapture data consisting of multiple non-invasive marks
#'
#' This function fits spatially-explicit population abundance models for capture-mark-recapture data consisting of multiple non-invasive marks using Bayesian analysis methods. Markov chain Monte Carlo (MCMC) is used to draw samples from the joint posterior distribution. 
#'
#'
#' @param Enc.Mat A matrix containing the observed encounter histories with rows corresponding to individuals and (\code{ntraps}*\code{noccas}) columns corresponding to traps and sampling occasions.  The first \code{noccas} columns correspond to trap 1, the second \code{noccas} columns corresopond to trap 2, etc. Ignored unless \code{mms=NULL}.
#' @param data.type Specifies the encounter history data type. All data types include non-detections (type 0 encounter), type 1 encounter (e.g., left-side), and type 2 encounters (e.g., right-side). When both type 1 and type 2 encounters occur for the same individual within a sampling occasion, these can either be "non-simultaneous" (type 3 encounter) or "simultaneous" (type 4 encounter). Three data types are currently permitted:
#' 
#'  \code{data.type="never"} indicates both type 1 and type 2 encounters are never observed for the same individual within a sampling occasion, and observed encounter histories therefore include only type 1 or type 2 encounters (e.g., only left- and right-sided photographs were collected). Observed encounter histories can consist of non-detections (0), type 1 encounters (1), and type 2 encounters (2). See \code{\link{bobcat}}. Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 3 encounters (3).
#'
#'  \code{data.type="sometimes"} indicates both type 1 and type 2 encounters are sometimes observed (e.g., both-sided photographs are sometimes obtained, but not necessarily for all individuals). Observed encounter histories can consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), type 3 encounters (3), and type 4 encounters (4). Type 3 encounters can only be observed when an individual has at least one type 4 encounter. Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), type 3 encounters (3), and type 4 encounters (4). 
#'
#'  \code{data.type="always"} indicates both type 1 and type 2 encounters are always observed, but some encounter histories may still include only type 1 or type 2 encounters. Observed encounter histories can consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 4 encounters (4). Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 4 encounters (4).
#'
#' @param covs A data frame of temporal covariates for detection probabilities (ignored unless \code{mms=NULL}). The number of rows in the data frame must equal the number of sampling occasions. Covariate names cannot be "time", "age", or "h"; these names are reserved for temporal, behavioral, and individual effects when specifying \code{mod.p} and \code{mod.phi}.
#' @param mms An optional object of class \code{multimarksetup-class}; if \code{NULL} it is created. See \code{\link{processdata}}.
#' @param spatialInputs A list of length 3 composed of objects named \code{trapCoords}, \code{studyArea}, and \code{centers}:
#' 
#'  \code{trapCoords} is a matrix of dimension \code{ntraps} x 2 indicating the Cartesian coordinates for the traps, where rows correspond to trap, the first column the x-coordinate, and the second column the y-coordinate. 
#'
#'  \code{studyArea} is a 3-column matrix containing the coordinates for the centroids a contiguous grid of cells that define the study area and available habitat. Each row corresponds to a grid cell. The first 2 columns indicate the Cartesian x- and y-coordinate for the centroid of each grid cell, and the third column indicates whether the cell is available habitat (=1) or not (=0). All cells must have the same resolution.
#'   
#'  \code{centers} is an integer vector indicating the grid cell (i.e., the row of \code{studyArea}) that contains the true (latent) coordinates of the activity centers for each individual in the population.
#'
#' @param mod.p Model formula for detection probability as a function of distance from activity centers. For example, \code{mod.p=~1} specifies no effects (i.e., intercept only) other than distance, \code{mod.p~time} specifies temporal effects, \code{mod.p~c} specifies behavioral reponse (i.e., trap "happy" or "shy"), and \code{mod.p~time+c} specifies additive temporal and behavioral effects.
#' @param mod.delta Model formula for conditional probabilities of type 1 (delta_1) and type 2 (delta_2) encounters, given detection. Currently only \code{mod.delta=~1} (i.e., \eqn{\delta_1 = \delta_2}) and \code{mod.delta=~type} (i.e., \eqn{\delta_1 \ne \delta_2}) are implemented.
#' @param detection Detection function for detection probability. Must be "\code{half-normal}" or "\code{exponential}".
#' @param parms A character vector giving the names of the parameters and latent variables to monitor. Possible parameters are cloglog-scale detection probability parameters ("\code{pbeta}"), population abundance ("\code{N}"), conditional probability of type 1 or type 2 encounter, given detection ("\code{delta})", probability of simultaneous type 1 and type 2 detection, given both types encountered ("\code{alpha}"), cloglog-scale variance term for the detection function ("\code{sigma2_scr}"), and the probability that a randomly selected individual from the \code{M = nrow(Enc.Mat)} observed individuals belongs to the \eqn{n} unique individuals encountered at least once ("\code{psi}"). Individual activity centers ("\code{centers}"), encounter history indices ("\code{H}"), and the log posterior density ("\code{logPosterior}") may also be monitored. Setting \code{parms="all"} monitors all possible parameters and latent variables.
#' @param nchains The number of parallel MCMC chains for the model.
#' @param iter The number of MCMC iterations.
#' @param adapt The number of iterations for proposal distribution adaptation. If \code{adapt = 0} then no adaptation occurs.
#' @param bin Bin length for calculating acceptance rates during adaptive phase (\code{0 < bin <= iter}).
#' @param thin Thinning interval for monitored parameters.
#' @param burnin Number of burn-in iterations (\code{0 <= burnin < iter}).
#' @param taccept Target acceptance rate during adaptive phase (\code{0 < taccept <= 1}). Acceptance rate is monitored every \code{bin} iterations. Default is \code{taccept = 0.44}.
#' @param tuneadjust Adjustment term during adaptive phase (\code{0 < tuneadjust <= 1}). If acceptance rate is less than \code{taccept}, then proposal term (\code{proppbeta}, \code{propzp}, or \code{propsigma}) is multiplied by \code{tuneadjust}. If acceptance rate is greater than or equal to \code{taccept}, then proposal term is divided by \code{tuneadjust}. Default is \code{tuneadjust = 0.95}.
#' @param proppbeta Scaler or vector (of length k) specifying the initial standard deviation of the Normal(pbeta[j], proppbeta[j]) proposal distribution. If \code{proppbeta} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{proppbeta = 0.1}.
#' @param propzp Scaler or vector (of length M) specifying the initial standard deviation of the Normal(zp[i], propzp[i]) proposal distribution. If \code{propzp} is a scaler, then this value is used for all i = 1, ..., M individuals. Default is \code{propzp = 1}.
#' @param propsigma Scaler specifying the initial Gamma(shape = 1/\code{propsigma}, scale = sigma_scr * \code{propsigma}) proposal distribution for sigma_scr = sqrt(sigma2_scr). Default is \code{propsigma=1}.
#' @param propcenter Scaler specifying the neighborhood distance (on re-scaled coordinates between 0 and \code{scalemax}) when proposing updates to activity centers. When \code{propcenter=NULL} (the default), then propcenter = a*10, where a is the re-scaled cell size for the study area grid.  When \code{propcenter=NULL} and \code{scalemax=10} (the default), each cell has approximately 300 neighbors (at most). 
#' @param maxnumbasis Maximum number of basis vectors to use when proposing latent history frequency updates. Default is \code{maxnumbasis = 1}, but higher values can potentially improve mixing.
#' @param a0delta Scaler or vector (of length d) specifying the prior for the conditional (on detection) probability of type 1 (delta_1), type 2 (delta_2), and both type 1 and type 2 encounters (1-delta_1-delta_2). If \code{a0delta} is a scaler, then this value is used for all a0delta[j] for j = 1, ..., d. For \code{mod.delta=~type}, d=3 with [delta_1, delta_2, 1-delta_1-delta_2] ~ Dirichlet(a0delta) prior. For \code{mod.delta=~1}, d=2 with [tau] ~ Beta(a0delta[1],a0delta[2]) prior, where (delta_1,delta_2,1-delta_1-delta_2) = (tau/2,tau/2,1-tau). See McClintock et al. (2013) for more details.
#' @param a0alpha Specifies "shape1" parameter for [alpha] ~ Beta(a0alpha, b0alpha) prior. Only applicable when \code{data.type = "sometimes"}. Default is \code{a0alpha = 1}. Note that when \code{a0alpha = 1} and \code{b0alpha = 1}, then [alpha] ~ Unif(0,1).
#' @param b0alpha Specifies "shape2" parameter for [alpha] ~ Beta(a0alpha, b0alpha) prior. Only applicable when \code{data.type = "sometimes"}. Default is \code{b0alpha = 1}. Note that when \code{a0alpha = 1} and \code{b0alpha = 1}, then [alpha] ~ Unif(0,1).
#' @param a Scale parameter for [sigma_z] ~ half-Cauchy(a) prior for the detection function term sigma_scr = sqrt(sigma2_scr). Default is ``uninformative'' \code{a = 25}.
#' @param mu0 Scaler or vector (of length k) specifying mean of pbeta[j] ~ Normal(mu0[j], sigma2_mu0[j]) prior. If \code{mu0} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{mu0 = 0}.
#' @param sigma2_mu0 Scaler or vector (of length k) specifying variance of pbeta[j] ~ Normal(mu0[j], sigma2_mu0[j]) prior. If \code{sigma2_mu0} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{sigma2_mu0 = 1.75}.
#' @param a0psi Specifies "shape1" parameter for [psi] ~ Beta(a0psi,b0psi) prior. Default is \code{a0psi = 1}.
#' @param b0psi Specifies "shape2" parameter for [psi] ~ Beta(a0psi,b0psi) prior. Default is \code{b0psi = 1}.
#' @param initial.values Optional list of \code{nchain} list(s) specifying initial values for parameters and latent variables. Default is \code{initial.values = NULL}, which causes initial values to be generated automatically. In addition to the parameters ("\code{pbeta}", "\code{N}", "\code{delta_1}", "\code{delta_2}", "\code{alpha}", "\code{sigma2_scr}", "\code{zp}", and "\code{psi}"), initial values can be specified for the initial latent history frequencies ("\code{x}") and initial individual encounter history indices ("\code{H}").
#' @param known Optional integer vector indicating whether the encounter history of an individual is known with certainty (i.e., the observed encounter history is the true encounter history). Encounter histories with at least one type 4 encounter are automatically assumed to be known, and \code{known} does not need to be specified unless there exist encounter histories that do not contain a type 4 encounter that happen to be known with certainty (e.g., from independent telemetry studies). If specified, \code{known = c(v_1,v_2,...,v_M)} must be a vector of length \code{M = nrow(Enc.Mat)} where \code{v_i = 1} if the encounter history for individual \code{i} is known (\code{v_i = 0} otherwise). Note that known all-zero encounter histories (e.g., `000') are ignored.
#' @param scalemax Upper bound for internal re-scaling of grid cell centroid coordinates. Default is \code{scalemax=10}, which re-scales the centroids to be between 0 and 10.  Re-scaling is done internally to avoid numerical overflows during model fitting.
#' @param printlog Logical indicating whether to print the progress of chains and any errors to a log file in the working directory. Ignored when \code{nchains=1}. Updates are printed to log file as 1\% increments of \code{iter} of each chain are completed. With >1 chains, setting \code{printlog=TRUE} is probably most useful for Windows users because progress and errors are automatically printed to the R console for "Unix-like" machines (i.e., Mac and Linux) when \code{printlog=FALSE}. Default is \code{printlog=FALSE}.
#' @param ... Additional "\code{parameters}" arguments for specifying \code{mod.p}. See \code{\link[RMark]{make.design.data}}.
#'
#' @details The first time \code{multimarkClosed} (or \code{\link{multimarkCJS}}) is called, it will likely produce a firewall warning alerting users that R has requested the ability to accept incoming network connections. Incoming network connections are required to use parallel processing as implemented in \code{multimarkClosed}. Note that setting \code{parms="all"} is required for any \code{multimarkClosed} model output to be used in \code{\link{multimodelClosed}}.
#' @return A list containing the following:
#' \item{mcmc}{Markov chain Monte Carlo object of class \code{\link[coda]{mcmc.list}}.}
#' \item{mod.p}{Model formula for detection probability (as specified by \code{mod.p} above).}
#' \item{mod.delta}{Model formula for conditional probability of type 1 or type 2 encounter, given detection (as specified by \code{mod.delta} above).}
#' \item{DM}{A list of design matrices for detection probability generated for model \code{mod.p}, where DM$p is the design matrix for initial capture probability (p) and DM$c is the design matrix for recapture probability (c).}
#' \item{initial.values}{A list containing the parameter and latent variable values at iteration \code{iter} for each chain. Values are provided for "\code{pbeta}", "\code{N}", "\code{delta_1}", "\code{delta_2}", "\code{alpha}", "\code{sigma2_scr}", "\code{zp}", "\code{psi}", "\code{x}", and "\code{H}".}
#' \item{mms}{An object of class \code{multimarksetup}}
#' @author Brett T. McClintock
#' @seealso \code{\link{bobcat}}, \code{\link{processdata}}, \code{\link{multimodelClosed}}
#' @references
#' Bonner, S. J., and Holmberg J. 2013. Mark-recapture with multiple, non-invasive marks. \emph{Biometrics} 69: 766-775.
#' 
#' McClintock, B. T., Conn, P. B., Alonso, R. S., and Crooks, K. R. 2013. Integrated modeling of bilateral photo-identification data in mark-recapture analyses. \emph{Ecology} 94: 1464-1471.
#' 
#' McClintock, B. T., Bailey, L. L., Dreher, B. P., and Link, W. A. 2014. Probit models for capture-recapture data subject to imperfect detection, individual heterogeneity and misidentification. \emph{The Annals of Applied Statistics} 8: 2461-2484.
#' @examples
#' \dontshow{
#' test<-multimarkClosed(Enc.Mat=bobcat,data.type="never",iter=10,burnin=0,bin=5)}
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # Example uses unrealistically low values for nchain, iter, and burnin
#' 
#' #Run single chain using the default model for bobcat data
#' bobcat.dot<-multimarkClosed(bobcat)
#' 
#' #Posterior summary for monitored parameters
#' summary(bobcat.dot$mcmc)
#' plot(bobcat.dot$mcmc)}
multimarkClosedSCR<-function(Enc.Mat,data.type="never",covs=data.frame(),mms=NULL,spatialInputs=NULL,mod.p=~1,mod.delta=~type,detection="half-normal",parms=c("pbeta","delta","N","sigma2_scr"),nchains=1,iter=12000,adapt=1000,bin=50,thin=1,burnin=2000,taccept=0.44,tuneadjust=0.95,proppbeta=0.1,propzp=1,propsigma=1,propcenter=NULL,maxnumbasis=1,a0delta=1,a0alpha=1,b0alpha=1,a=25,mu0=0,sigma2_mu0=1.75,a0psi=1,b0psi=1,initial.values=NULL,known=integer(),scalemax=10,printlog=FALSE,...){
  
  if(is.null(mms)) mms <- processdata(Enc.Mat,data.type,covs,known)
  if(class(mms)!="multimarksetup") stop("'mms' must be an object of class 'multimarksetup'")
  validObject(mms)
  
  if(class(mod.p)!="formula") stop("'mod.p' must be an object of class 'formula'")
  if(class(mod.delta)!="formula") stop("'mod.delta' must be an object of class 'formula'")
  DM<-get_DMClosed(mod.p,mod.delta,mms@Enc.Mat,covs=mms@covs,ntraps=nrow(spatialInputs$trapCoords),detection=detection,...)
  
  if(iter>0){
    if(iter<=burnin) stop(paste("'burnin' must be less than ",iter))
  } else {
    burnin<-0
  }
  
  if(mod.delta != ~NULL) {
    parmlist<-c("pbeta","delta","N","sigma2_scr","alpha","psi","H","logPosterior","centers")
  } else {
    parmlist<-c("pbeta","N","sigma2_scr","logPosterior","centers")    
  }
  params <- checkClosedSCR(parms,parmlist,mms,DM,iter,adapt,bin,thin,burnin,taccept,tuneadjust,maxnumbasis,a0delta,a0alpha,b0alpha,a,sigma2_mu0,a0psi,b0psi)
  
  data.type<-mms@data.type
  Enc.Mat<-mms@Enc.Mat
  M<-nrow(Enc.Mat)
  noccas<-ncol(Enc.Mat)
  covs<-mms@covs
  pdim<-ncol(DM$p)
  
  mu0 <- checkvecs(mu0,pdim,"mu0")
  sigma2_mu0 <- checkvecs(sigma2_mu0,pdim,"sigma2_mu0")
  a0delta <- checkvecs(a0delta,ifelse(mod.delta==formula(~type),3,2),"a0delta")
  
  S <- spatialInputs$studyArea # total study area
  A <- subset(S,"avail">0,c("x","y"))  # available habitat study area
  minCoord <- apply(A, 2, min) 
  Arange <- max(apply(A, 2, max) - minCoord)/scalemax 
  
  availSpatialInputs=list()
  availSpatialInputs$studyArea <- scale(A, center=minCoord, scale=rep(Arange,2))
  availSpatialInputs$trapCoords <- scale(spatialInputs$trapCoords, center=minCoord, scale=rep(Arange,2))
  availSpatialInputs$a <- sp::points2grid(sp::SpatialPoints(availSpatialInputs$studyArea[,1:2]))@cellsize[1]
  availSpatialInputs$A <- availSpatialInputs$a * sum(spatialInputs$studyArea[,"avail"])
  availSpatialInputs$dist2 <- getdist(availSpatialInputs$studyArea,availSpatialInputs$trapCoords)
  
  inits<-get_initsSCR(mms,nchains,initial.values,M,data.type,a0alpha,b0alpha,a0delta,a0psi,b0psi,DM,availSpatialInputs)
  
  priorparms <-list(a0delta=a0delta,a0alpha=a0alpha,b0alpha=b0alpha,a=a,mu0=mu0,sigma2_mu0=sigma2_mu0,a0psi=a0psi,b0psi=b0psi)
  
  message("\nFitting spatial abundance model with cloglog link\n")
  if(mod.delta != ~NULL) message("data type = \"",data.type,"\"\n")
  message("p model = ",as.character(mod.p))
  if(mod.delta != ~NULL) message("delta model = ",as.character(mod.delta))
  message("\nInitializing model \n")
  posteriorClosedSCR(inits,DM,mms,priorparms,availSpatialInputs)
  
  propzp <- checkvecs(propzp,M,"propzp")
  proppbeta <- checkvecs(proppbeta,pdim,"proppbeta")
  if(length(propsigma)!=1) stop("'propsigma' must be a scaler")
  
  NN<-list()
  RAD <- ifelse(is.null(propcenter),availSpatialInputs$a*10,propcenter) # Change propcenter to get more or fewer neighbors
  for(i in 1:nrow(availSpatialInputs$studyArea)){
    od <- sqrt( (availSpatialInputs$studyArea[i,1]-availSpatialInputs$studyArea[,1])^2  +  (availSpatialInputs$studyArea[i,2]-availSpatialInputs$studyArea[,2])^2  )
    od <- (1:length(od))[od < RAD]
    NN[[i]]<-od
  }
  Prop.center=list(NNvect=unlist(NN),numnn=lapply(NN,length))
  
  Prop.sd <- c(propzp,proppbeta,propsigma)
  
  message("Updating...",ifelse(printlog | nchains==1,"","set 'printlog=TRUE' to follow progress of chains in a working directory log file"),"\n",sep="")
  if(printlog & nchains==1) printlog<-FALSE
  
  if(nchains>1){
    if(nchains>detectCores()) warning("Number of parallel chains (nchains) is greater than number of cores \n")
    modlog <- ifelse(mod.delta != ~NULL,"multimarkClosed","markClosed")
    cl <- makeCluster( nchains ,outfile=ifelse(printlog,paste0(modlog,"_log_",format(Sys.time(), "%Y-%b-%d_%H%M.%S"),".txt"),""))
    clusterExport(cl,list("mcmcClosed"),envir=environment())  
    chains <- parLapply(cl,1:nchains, function(ichain) mcmcClosedSCR(ichain,mms,DM,params,inits,iter,adapt,bin,thin,burnin,taccept,tuneadjust,Prop.sd,Prop.center,availSpatialInputs,maxnumbasis,a0delta,a0alpha,b0alpha,a,mu0,sigma2_mu0,a0psi,b0psi,printlog))
    stopCluster(cl)
    gc()
  } else {
    chains <- vector('list',nchains)
    chains[[nchains]] <- mcmcClosedSCR(nchains,mms,DM,params,inits,iter,adapt,bin,thin,burnin,taccept,tuneadjust,Prop.sd,Prop.center,availSpatialInputs,maxnumbasis,a0delta,a0alpha,b0alpha,a,mu0,sigma2_mu0,a0psi,b0psi,printlog)
    gc()
  }
  
  chains <- processClosedchains(chains,params,DM,M,noccas,nchains,iter,burnin,thin)
  return(list(mcmc=chains$chains,mod.p=mod.p,mod.delta=mod.delta,DM=list(p=DM$p,c=DM$c),initial.values=chains$initial.values,priorparms=priorparms,mms=mms))
}