#' Bobcat data
#' 
#' Example bobcat data for \code{multimark} package.
#' @name bobcat
#' @docType data
#' @format The data are summarized in a 46x8 matrix containing observed encounter histories for 46 bobcats across 8 sampling occasions. Bobcats are bilaterially asymmetrical, and sampling was conducted using camera stations consisting of a single camera. 
#' 
#' Because the left-side cannot be reconciled with the right-side, the two types of ``marks'' in this case are the pelage patterns on the left- and right-side of each individual. Encounter type 0 corresponds to non-detection, encounter type 1 corresponds to left-sided detection, encounter type 2 corresponds to right-sided detection. 
#' 
#' Both-sided encounters were never observed in this dataset, hence the most appropriate \code{multimark} data type is \code{data.type="never".}
#' @seealso \code{\link{multimarkClosed}}, \code{\link{processdata}}
#' @source 
#' McClintock, B. T., Conn, P. B., Alonso, R. S., and Crooks, K. R. 2013. Integrated modeling of bilateral photo-identification data in mark-recapture analyses. \emph{Ecology} 94: 1464-1471.
#' @examples
#' data(bobcat)
#' @keywords data datasets
NULL

#' Class \code{"multimarksetup"}
#'
#' A class of 'mulitmark' model inputs
#'
#'
#' @section Objects from the Class:
#' Objects can be created by calls of the form \code{processdata(Enc.Mat, ...)} or \code{new("multimarksetup", ...)}.
#' @slot Enc.Mat Object of class \code{"matrix"}. The observed encounter histories (with rows corresponding to individuals and columns corresponding to sampling occasions).
#' @slot data.type Object of class \code{"character"}. The encounter history data type ("never", "sometimes", or "always") .
#' @slot vAll.hists Object of class \code{"integer"}. An ordered vector containing all possible encounter histories in sequence .
#' @slot Aprime Object of class \code{"sparseMatrix"}. Transpose of the A matrix mapping latent encounter histories to observed histories .
#' @slot indBasis Object of class \code{"numeric"}.An ordered vector of the indices of the three encounter histories updated by each basis vector .
#' @slot ncolbasis Object of class \code{"integer"}. The number of needed basis vectors .
#' @slot knownx Object of class \code{"integer"}. Frequencies of known encounter histories .
#' @slot C Object of class \code{"integer"}. Sampling occasion of first capture for each encounter history .
#' @slot L Object of class \code{"integer"}. Sampling occasion of last capture for each encounter history .
#' @slot naivex Object of class \code{"integer"}. ``Naive'' latent history frequencies assuming a one-to-one mapping with \code{Enc.Mat} .
#' @slot covs Object of class \code{"data.frame"}. Temporal covariates for detection probability (the number of rows in the data frame must equal the number of sampling occasions).
#' @section Methods:
#' No methods defined with class "multimarksetup".
#' @author Brett T. McClintock
#' @seealso \code{\link{processdata}}
#' @examples
#' showClass("multimarksetup")
#' @keywords classes
setClass("multimarksetup", representation=list(Enc.Mat="matrix",data.type="character",vAll.hists="integer",Aprime="sparseMatrix",indBasis="integer",ncolbasis="integer",knownx="integer",C="integer",L="integer",naivex="integer",covs="data.frame"),
         prototype=list(Enc.Mat=matrix(0,0,0),data.type=character(),vAll.hists=integer(),Aprime=Matrix(0,0,0),indBasis=integer(),ncolbasis=integer(),knownx=integer(),C=integer(),L=integer(),naivex=integer(),covs=data.frame()),
         package="multimark")

tol <- 1.e-6

#' Generate encounter history frequencies
#'
#' This function generates an ordered vector of encounter history frequencies 
#' 
#' @param Enc.Mat A matrix of observed encounter histories with rows corresponding to individuals and columns corresponding to sampling occasions.
#' @param data.type Specifies the encounter history data type. All data types include non-detections (type 0 encounter), type 1 encounter (e.g., left-side), and type 2 encounters (e.g., right-side). When both type 1 and type 2 encounters occur for the same individual within a sampling occasion, these can either be "non-simultaneous" (type 3 encounter) or "simultaneous" (type 4 encounter). Three data types are currently permitted:
#' 
#'  \code{data.type="never"} indicates both type 1 and type 2 encounters are never observed for the same individual within a sampling occasion, and observed encounter histories therefore include only type 1 or type 2 encounters (e.g., only left- and right-sided photographs were collected). Observed encounter histories can consist of non-detections (0), type 1 encounters (1), and type 2 encounters (2). See \code{\link{bobcat}}. Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 3 encounters (3).
#'
#'  \code{data.type="sometimes"} indicates both type 1 and type 2 encounters are sometimes observed (e.g., both-sided photographs are sometimes obtained, but not necessarily for all individuals). Observed encounter histories can consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), type 3 encounters (3), and type 4 encounters (4). Type 3 encounters can only be observed when an individual has at least one type 4 encounter. Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), type 3 encounters (3), and type 4 encounters (4). 
#'
#'  \code{data.type="always"} indicates both type 1 and type 2 encounters are always observed, but some encounter histories may still include only type 1 or type 2 encounters. Observed encounter histories can consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 4 encounters (4). Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 4 encounters (4).
#'
#' @details This function calculates frequencies for the latent encounter history frequencies by assuming a 1-to-1 mapping of the observed histories (\code{Enc.Mat}) to latent histories. 
#' 
#' For the \code{"never"} and \code{"always"} data types, the index of latent history \eqn{i = 1 + \sum_{t=1}^T \omega_t*4^{T-t}}{i = 1 + \sum_{t=1}^T \omega_t*4^{T-t}}, with \eqn{\omega_t \in {0, 1, 2, 3}} to represent a non-encounter, a type 1 encounter, a type 2 encounter, and both type 1 and type 2 encounter (i.e., type 3 encounter for \code{"never"} and type 4 encounter for \code{"always"}), respectively, for sampling occasion \eqn{t=1,\ldots,T}. 
#' 
#' For the "sometimes" data type, the index of latent history \eqn{i = 1 + \sum_{t=1}^T \omega_t*5^{T-t}}{i = 1 + \sum_{t=1}^T \omega_t*5^{T-t}}, with \eqn{\omega_t \in {0, 1, 2, 3,  4}} to represent a non-encounter, a type 1 encounter, a type 2 encounter, type 3 encounter, and type 4 encounter, respectively, for sampling occasion \eqn{t=1,\ldots,T}.
#' @return A vector of encounter history frequencies of length \eqn{4^T} for the "never" and "always" data types or of length \eqn{5^T} for the "sometimes" data type.
#' @author Brett T. McClintock
#' @references 
#' Bonner, S. J., and Holmberg J. 2013. Mark-recapture with multiple, non-invasive marks. \emph{Biometrics} 69: 766-775.
#' 
#' McClintock, B. T., Conn, P. B., Alonso, R. S., and Crooks, K. R. 2013. Integrated modeling of bilateral photo-identification data in mark-recapture analyses. \emph{Ecology} 94: 1464-1471.
#' @examples
#' #simulate closed capture data with 2 sampling occasions for data.type="sometimes"
#' data_type <- "sometimes"
#' data <- simdataClosed(noccas=2,data.type=data_type)
#'  
#' #generate initial latent history freqencies from observed encounter histories
#' getfreq(data$Enc.Mat,data.type=data_type)
#'  
#' #calculate freqencies for true encounter histories
#' getfreq(data$trueEnc.Mat,data.type=data_type)
getfreq <- function(Enc.Mat,data.type){
  
  temp.Enc.Mat <- Enc.Mat
  
  if(data.type=="never"){
    ind<-4
  } else if(data.type=="always"){
    ind<-4
    temp.Enc.Mat[which(Enc.Mat==4)]=3
  } else if(data.type=="sometimes"){
    ind<-5
  } 
  M<-nrow(temp.Enc.Mat)
  temp.noccas<-ncol(temp.Enc.Mat)
  Hist.num<-integer(M)
  Value<-ind^c((temp.noccas-1):0)
  Hist.num <- 1 + temp.Enc.Mat %*% Value
  
  temp.x<-integer(ind^temp.noccas)
  temp.x[sort(unique(Hist.num))]=table(Hist.num)
  temp.x[1] <- M-base::sum(temp.x[-1])
  names(temp.x)=paste0("x[",1:ind^temp.noccas,"]")
  as.vector(temp.x,mode="integer")
}

get_A <- function(noccas,data.type){
  
  if(data.type=="never"){
    ind<-4
    r<-4^noccas-2*(2^noccas-1)          #r is the number of free variables (i.e., the number of basis vectors)
  } else if(data.type=="always"){
    ind<-4
    r<-3^noccas-2^(noccas+1)+2
  } else if(data.type=="sometimes"){
    ind<-5
    r<-4^noccas-2*(2^noccas-1) 
  }
  
  cat("Dimension of basis for null space of A is ",ind^noccas," x ",r," = ",ind^noccas*r,"\n")
  if(noccas>6){
    cat("  This might be a while...that's a darned big null space! \n")
  }
  
  #first, create correponding matrix w/all possible latent histories using recursive algorithm
  dimAllhists <- (ind^noccas)*noccas
  nonzeroAllhists<-ind^(noccas-1)*noccas*(ind-1)
  indexAllhists<-integer(nonzeroAllhists)
  valueAllhists<-integer(nonzeroAllhists)
  
  indexcount<-0
  for(i in 1:noccas){
    if(data.type=="never"){
      series<-rep(c(rep(0,ind^(noccas-i)),rep(1,ind^(noccas-i)),rep(2,ind^(noccas-i)),rep(3,ind^(noccas-i))),ind^(i-1))     
    } else if(data.type=="sometimes"){
      series<-rep(c(rep(0,ind^(noccas-i)),rep(1,ind^(noccas-i)),rep(2,ind^(noccas-i)),rep(3,ind^(noccas-i)),rep(4,ind^(noccas-i))),ind^(i-1))      
    } else if(data.type=="always"){
      series<-rep(c(rep(0,ind^(noccas-i)),rep(1,ind^(noccas-i)),rep(2,ind^(noccas-i)),rep(4,ind^(noccas-i))),ind^(i-1))
    }
    seriesind<-which(series>0)
    indexAllhists[indexcount+1:length(seriesind)] <- seriesind+(i-1)*dimAllhists/noccas;
    valueAllhists[indexcount+1:length(seriesind)] <- series[seriesind];
    indexcount<-indexcount+length(seriesind)
  }
  
  All.hists<-Matrix(0,nrow=(ind^noccas),ncol=noccas)
  All.hists[indexAllhists] <- valueAllhists
  
  Value<-ind^c((noccas-1):0)
  
  # Construct A matrix   
  ivect<-which(((rowSums(All.hists==1)>0 & rowSums(All.hists==2)>0) | rowSums(All.hists==3)>0) & (rowSums(All.hists==4)==0))  
  temp.hist<-All.hists[ivect,]
  temp.1<-Matrix(0,nrow=length(ivect),ncol=noccas)
  temp.2<-Matrix(0,nrow=length(ivect),ncol=noccas)
  temp.1[which(temp.hist==1 | temp.hist>2)] <- 1
  temp.2[which(temp.hist>1)] <- 2
  places1<-as.vector(temp.1 %*% Value) 
  places2<-as.vector(temp.2 %*% Value)
  
  A<-sparseMatrix(i=c(ivect,ivect),j=c(places1+1,places2+1),dims=c(ind^noccas,ind^noccas),x=1)
  diag(A)[-ivect] <- 1
  A<-A[,-1]
  A<-A[,-which(colSums(A)==0)]
  
  gc()
  out<-list(Aprime=t(A),vAll.hists=as.vector(t(All.hists),mode="integer"),ivect=ivect)
}

get_basis_vectors <- function(noccas,tA,ivect,div,data.type){
  #This function caculates basis vectors based on data type (data.type) and latent frequencies (x).  
  #Function returns a matrix of the relevant basis vectors for the null space of A'.
  
  # Arguments: 
  # noccas = number of sampling occasions (T in paper)
  # A and ivect are objects returned by "get_A" above.
  # div = integer for dividing matrix manipulation workload for very large null spaces. If experiencing memory issues, try larger values for div.  If memory not an issue and computation is too slow, try lower values for div.
  # data.type = data type that determines mapping of recorded histories to latent histories (see Table 1 in paper). 
  #   Data type "never" indicates simultaneous type 1 and type 2 detections are never observed, "sometimes" indicates simultaneous type 1 and type 2 detections are sometimes observed, and "always" indicates simultaneous type 1 and type 2 detections are always observed
  
  #determine basis vectors
  if(data.type=="never"){
    ind<-4
    r<-4^noccas-2*(2^noccas-1)                   #r is the number of free variables (i.e., the number of basis vectors)
  } else if(data.type=="always"){
    ind<-4
    r<-3^noccas-2^(noccas+1)+2
  } else if(data.type=="sometimes"){
    ind<-5
    r<-4^noccas-2*(2^noccas-1) 
  }
  
  cat("Constructing basis vectors \n")
  
  free<-c(1,ivect)   # indices for the free variables (i.e., the latent histories that spawn >1 recorded history)
  
  bound<-seq(1:(ind^noccas))[-free]  # indices for the bound variables (i.e., the latent histories that spawn only 1 recorded history)
  
  k<-length(bound)
  
  Basis<-Diagonal(ind^noccas)[,-bound]
  
  if(div>1 & noccas>3){
    if(div>ceiling(r/2)){stop(paste0("'divBasis' must be less than ",ceiling(r/2)+1," when 'data.type'=",data.type," and 'noccas'=",noccas))}
    div_r  <-  ceiling(seq(2,r+1,length=div))
    pb <- txtProgressBar(min=1,max=length(div_r),char="+",width=100,style=3)
    setTxtProgressBar(pb, 1)
    for(i in 1:(length(div_r)-1)){
      Basis[bound,(div_r[i]):(div_r[i+1]-1)] <- -tA[,free[(div_r[i]):(div_r[i+1]-1)]] 
      setTxtProgressBar(pb, i+1)
      gc()
    }
    close(pb)
  } else {
    Basis[bound,-1] <- -tA[,free[-1]]     
  }
  gc()
  return(Basis)
}

reduced_basis<-function(noccas,Basis,x,div,data.type){
  #Given the initial frequency vector 'x', eliminate Basis vectors which always produce negative frequencies.  
  #For example, with T=3, if latent history \omega_2 = '001' has frequency x_2 = 0, then any basis vector with a '-1' in the second row will produce negative frequencies and can be eliminated. For the never and sometimes data types, one would eliminate basis vectors 2, 5, 6, 23, 24, 29, and 30. For the always data type, one would eliminate basis vectors 3, 9, and 13.
  
  # Arguments: 
  # noccas = number of sampling occasions (T in paper)
  # Basis is basis vectors object returned by "get_basis_vectors" above.
  # x = initial frequencies for latent histories return by "getfreq" above.
  # div = integer for dividing matrix manipulation workload for very large null spaces. If experiencing memory issues, try larger values for div.  If memory not an issue and computation is too slow, try lower values for div.
  # data.type = data type that determines mapping of recorded histories to latent histories (see Table 1 in paper). 
  #   Data type "never" indicates simultaneous type 1 and type 2 detections are never observed, "sometimes" indicates simultaneous type 1 and type 2 detections are sometimes observed, and "always" indicates simultaneous type 1 and type 2 detections are always observed
  
  if(data.type=="never"){
    ind<-4
    r<-4^noccas-2*(2^noccas-1)                   #r is the number of free variables (i.e., the number of basis vectors)
  } else if(data.type=="always"){
    ind<-4
    r<-3^noccas-2^(noccas+1)+2
  } else if(data.type=="sometimes"){
    ind<-5
    r<-4^noccas-2*(2^noccas-1) 
  }
  
  cat("Reducing dimension of basis based on observed histories \n")
  x<-c(1,which(x[-1]>0)+1)
  freqvect<-numeric(ind^noccas)
  freqvect[x] <- 1
  if(div>1 & noccas>3){
    if(div>ceiling(r/2)){stop(paste0("'divredBasis' must be less than ",ceiling(r/2)+1," when 'data.type'=",data.type," and 'noccas'=",noccas))}
    temp<-numeric(r)
    temp[1]<-base::sum((freqvect+Basis[,1])<0)
    div_r  <-  ceiling(seq(2,r+1,length=div))
    pb <- txtProgressBar(min=1,max=length(div_r),char="+",width=100,style=3)
    setTxtProgressBar(pb, 1)
    for(i in 1:(length(div_r)-1)){
      sequ <- (div_r[i]):(div_r[i+1]-1)
      temp[sequ] <- colSums((freqvect+Basis[,sequ])<0)
      setTxtProgressBar(pb, i+1)
      gc()
    }
    close(pb)
  } else {
    temp<-colSums((freqvect+Basis)<0)
  }
  gc()
  return(Basis[,which(temp<1)])
}

get_Enc <- function(tEnc.Mat,data.type){
  noccas <- ncol(tEnc.Mat)
  Enc.Mat <- tEnc.Mat
  z.Enc.Mat<-which(rowSums(Enc.Mat)==0)
  if(length(z.Enc.Mat)) {
    Enc.Mat<-Enc.Mat[-z.Enc.Mat,]        #remove all zero simulated histories
  }
  
  if(data.type=="never"){
    for(i in 1:nrow(Enc.Mat)){      #histories with 1 and 2 or 3 result in 2 'ghost' histories; change 2's to 0 and 3's to 1 to get 1 ghosts, then change 3's to 2 and append with 2 ghosts to end of array
      if((base::sum(Enc.Mat[i,]==1)>0 & base::sum(Enc.Mat[i,]==2)>0) | base::sum(Enc.Mat[i,]==3)>0) {
        cur2<-which(Enc.Mat[i,]==2)
        cur3<-which(Enc.Mat[i,]==3)
        Enc.Mat[i,cur2] <- 0
        Enc.Mat[i,cur3] <- 1
        new.hist<-rep(0,noccas)
        new.hist[cur2] <- 2
        new.hist[cur3] <- 2
        Enc.Mat<-rbind(Enc.Mat,new.hist,deparse.level=0)
      }
    }
  } else if(data.type=="always"){
    for(i in 1:nrow(Enc.Mat)){      #histories with 1 and 2 but no 4 result in 2 'ghost' histories; change 2's to 0 to get 1 ghosts and append 2 ghosts to end of array
      if(base::sum(Enc.Mat[i,]==1)>0 & base::sum(Enc.Mat[i,]==2)>0 & base::sum(Enc.Mat[i,]==4)==0) {
        cur2<-which(Enc.Mat[i,]==2)
        Enc.Mat[i,cur2] <- 0
        new.hist<-rep(0,noccas)
        new.hist[cur2] <- 2
        Enc.Mat<-rbind(Enc.Mat,new.hist,deparse.level=0)
      }
    }
  } else if(data.type=="sometimes"){
    for(i in 1:nrow(Enc.Mat)){      #histories with 1 and 2 or 3 but no 4 result in 2 'ghost' histories; change 2's to 0 and 3's to 1 to get 1 ghosts, then change 3's to 2 and append with 2 ghosts to end of array
      if(base::sum(Enc.Mat[i,]==1)>0 & base::sum(Enc.Mat[i,]==2)>0 & base::sum(Enc.Mat[i,]==4)==0 | (base::sum(Enc.Mat[i,]==3)>0 & base::sum(Enc.Mat[i,]==4)==0)) {
        cur2<-which(Enc.Mat[i,]==2)
        cur3<-which(Enc.Mat[i,]==3)
        Enc.Mat[i,cur2] <- 0
        Enc.Mat[i,cur3] <- 1
        new.hist<-rep(0,noccas)
        new.hist[cur2] <- 2
        new.hist[cur3] <- 2
        Enc.Mat<-rbind(Enc.Mat,new.hist,deparse.level=0)
      }
    }
  }
  Enc.Mat
}

get_H <- function(mms,x){
  
  if(all(x==mms@naivex)){
    
    temp.Enc.Mat <- mms@Enc.Mat
    
    if(mms@data.type=="never"){
      ind<-4
    } else if(mms@data.type=="always"){
      ind<-4
      tmp.simult <- which(mms@Enc.Mat==4)
      if(length(tmp.simult)){
        temp.Enc.Mat[tmp.simult] <- 3
      }
    } else if(mms@data.type=="sometimes"){
      ind<-5
    }
    
    noccas <- ncol(temp.Enc.Mat)
    Value<-ind^c((noccas-1):0)
    H <- 1 + temp.Enc.Mat %*% Value
    
  } else {
    H<-integer(base::sum(x))
    cumx<-c(0,cumsum(x))
    for(i in which(x>0)){
      H[cumx[i]+1:x[i]] <- i
    }
  }
  H
}

get_C <-function(noccas,data.type){
  if(data.type=="never"){
    ind<-4
  } else if(data.type=="always"){
    ind<-4
  } else if(data.type=="sometimes"){
    ind<-5
  }
  Cind<-ind^(seq(1,noccas))-ind^(seq(0,noccas-1))
  c(as.integer(noccas+1),rep(seq(noccas,1),times=Cind))
}

get_L <-function(noccas,data.type){
  if(data.type=="never"){
    ind<-4
  } else if(data.type=="always"){
    ind<-4
  } else if(data.type=="sometimes"){
    ind<-5
  }
  
  L<-rep(noccas,ind^noccas)
  for(i in 1:(noccas-1)){
    L[seq(1,ind^noccas,ind^i)]=noccas-i
  }
  L[1]=0
  as.integer(L)
}
  
expit<-function(x){
  1/(1+exp(-x))
}

expittol<-function(x){
  l <- dim(x)
  expittol <- pmin(pmax(tol,expit(x)),1.-tol)
  if(!is.null(l)){
    expittol <- matrix(expittol,nrow=l[1],ncol=l[2])
  }
  return(expittol)
}

rdirichlet<-function (n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
}

ddirichlet <- function(x, alpha) {
  dirichlet1 <- function(x, alpha) {
    logD <- base::sum(lgamma(alpha)) - lgamma(base::sum(alpha))
    s <- base::sum((alpha - 1) * log(x))
    exp(base::sum(s) - logD)
  }
  if (!is.matrix(x)) 
    if (is.data.frame(x)) 
      x <- as.matrix(x)
  else x <- t(x)
  if (!is.matrix(alpha)) 
    alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x), 
                    byrow = TRUE)
  if (any(dim(x) != dim(alpha))) 
    stop("Mismatch between dimensions of x and alpha in ddirichlet().\n")
  pd <- vector(length = nrow(x))
  for (i in 1:nrow(x)) pd[i] <- dirichlet1(x[i, ], alpha[i, 
                                                         ])
  pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
  pd[apply(x, 1, function(z) all.equal(base::sum(z), 1) != TRUE)] <- 0
  return(log(pd))
}

integrandlogit<-function(x,beta,DM){
  dim<-length(DM)
  expit(DM%*%beta[1:dim]+x)*dnorm(x,0,sqrt(beta[dim+1]))
}

integrandprobit<-function(x,beta,DM){
  dim<-length(DM)
  pnorm(DM%*%beta[1:dim]+x)*dnorm(x,0,sqrt(beta[dim+1]))
}

isValidx<-function(mms,x,M){
  ifelse(ncol(mms@Aprime)==length(x),
         ifelse((all((mms@Aprime%*%x)==(mms@Aprime%*%mms@naivex)) & all(x>=mms@knownx) & base::sum(x)==M),TRUE,FALSE),
         FALSE)
}

get_inits<-function(mms,nchains,initial.values,M,data.type,a0alpha,b0alpha,a0delta,DM,gq=NULL){
  
  inits<-vector("list",nchains)
  
  if(!is.null(initial.values)){
    for(ichain in 1:length(initial.values)){
      if(!is.list(initial.values[[ichain]])) stop(paste0("'initial.values' for chain ",ichain," must be a list"))
    }
    if(length(initial.values)<nchains){
      for(ichain in (length(initial.values)+1):nchains){
        initial.values[[ichain]]<-list()
      }
    }
  } else {
    initial.values<-inits
  }
  
  mod.p.h<-DM$mod.p.h
  pdim<-ncol(DM$p)
  
  for(ichain in 1:nchains){
    
    if(length(initial.values[[ichain]]$H)){
      tab<-table(initial.values[[ichain]]$H)
      if(length(initial.values[[ichain]]$x)){
        if(isValidx(mms,initial.values[[ichain]]$x,M)){
          inits[[ichain]]$x<-initial.values[[ichain]]$x
        } else {
          stop(paste("impermissible initial latent frequencies (x) for chain",ichain))
        }
      } else {
        initial.values[[ichain]]$x<-rep(0,length(mms@naivex))
        initial.values[[ichain]]$x[as.integer(names(tab))]<-as.vector(tab)
        if(isValidx(mms,initial.values[[ichain]]$x,M)){
          inits[[ichain]]$x<-initial.values[[ichain]]$x
        } else {
          stop(paste("impermissible initial indiviual histories (H) for chain",ichain))
        }
      }
      if(any(tab!=inits[[ichain]]$x[as.integer(names(tab))])){
        stop(paste("initial individual histories (H) not compatible with initial latent frequencies (x) for chain",ichain)) 
      } else {
        inits[[ichain]]$H<-initial.values[[ichain]]$H         
      }
    } else {
      if(length(initial.values[[ichain]]$x)){
        if(isValidx(mms,initial.values[[ichain]]$x,M)){
          inits[[ichain]]$x<-initial.values[[ichain]]$x
        } else {
          stop(paste("impermissible initial latent frequencies (x) for chain",ichain))
        }
      } else {
        inits[[ichain]]$x<-mms@naivex
      }
      inits[[ichain]]$H<-get_H(mms,inits[[ichain]]$x)
    }
    
    if(length(initial.values[[ichain]]$pbeta)){
      if(length(initial.values[[ichain]]$pbeta)==pdim){
        inits[[ichain]]$pbeta<-initial.values[[ichain]]$pbeta
      } else {
        stop(paste("initial values for pbeta must be vector of length",pdim))
      }
    } else {
      inits[[ichain]]$pbeta<-rnorm(pdim)
    }
    if(mod.p.h){
      if(length(initial.values[[ichain]]$sigma2_zp)){
        if(length(initial.values[[ichain]]$sigma2_zp)==1 & initial.values[[ichain]]$sigma2_zp>0){
          inits[[ichain]]$sigma2_zp<-max(initial.values[[ichain]]$sigma2_zp,tol)
        } else {
          stop("initial value for sigma2_zp must be a positive scalar")
        }
      } else {
        inits[[ichain]]$sigma2_zp<-runif(1,tol,5)
      }
      if(length(initial.values[[ichain]]$zp)){
        if(length(initial.values[[ichain]]$zp)==M){
          inits[[ichain]]$zp<-initial.values[[ichain]]$zp
        } else {
          stop(paste("initial values for zp must be vector of length",M))
        }
      } else {
        inits[[ichain]]$zp<-rnorm(M,0.0,sqrt(inits[[ichain]]$sigma2_zp))
      }
    } else {
      inits[[ichain]]$sigma2_zp<-0.0
      inits[[ichain]]$zp<-rep(0.0,M)
    }
    if(!is.null(DM$phi)){
      
      mod.phi.h<-DM$mod.phi.h
      phidim<-ncol(DM$phi)  
      
      if(length(initial.values[[ichain]]$phibeta)){
        if(length(initial.values[[ichain]]$phibeta)==phidim){
          inits[[ichain]]$phibeta<-initial.values[[ichain]]$phibeta
        } else {
          stop(paste("initial values for phibeta must be vector of length",phidim))
        }
      } else {
        inits[[ichain]]$phibeta<-rnorm(phidim)
      }
      if(mod.phi.h){
        if(length(initial.values[[ichain]]$sigma2_zphi)){
          if(length(initial.values[[ichain]]$sigma2_zphi)==1 & initial.values[[ichain]]$sigma2_zphi>0){
            inits[[ichain]]$sigma2_zphi<-max(initial.values[[ichain]]$sigma2_zphi,tol)
          } else {
            stop("initial value for sigma2_zphi must be a positive scalar")
          }
        } else {
          inits[[ichain]]$sigma2_zphi<-runif(1,tol,5)
        }
        if(length(initial.values[[ichain]]$zphi)){
          if(length(initial.values[[ichain]]$zphi)==M){
            inits[[ichain]]$zphi<-initial.values[[ichain]]$zphi
          } else {
            stop(paste("initial values for zphi must be vector of length",M))
          }
        } else {
          inits[[ichain]]$zphi<-rnorm(M,0.0,sqrt(inits[[ichain]]$sigma2_zphi))
        }
      } else {
        inits[[ichain]]$sigma2_zphi<-0.0
        inits[[ichain]]$zphi<-rep(0.0,M)
      }
      if(length(initial.values[[ichain]]$z)){
        if(all(dunif(initial.values[[ichain]]$z)==1) & is.matrix(initial.values[[ichain]]$z) & all(dim(initial.values[[ichain]]$z)==dim(mms@Enc.Mat))){
          inits[[ichain]]$z<-initial.values[[ichain]]$z
        } else {
          stop(paste0("initial values for 'z' must be a ",M,"x",ncol(mms@Enc.Mat)," binary (0,1) matrix"))
        }
      } else {
        inits[[ichain]]$z<-get_z(mms,DM,inits[[ichain]]$H,inits[[ichain]]$pbeta,inits[[ichain]]$zp,inits[[ichain]]$phibeta,inits[[ichain]]$zphi)
      }
    } else {
      if(length(initial.values[[ichain]]$N)){
        if(length(initial.values[[ichain]]$N)==1 & initial.values[[ichain]]$N>=base::sum(inits[[ichain]]$H>1) & (!abs(initial.values[[ichain]]$N-round(initial.values[[ichain]]$N))>0)){
          inits[[ichain]]$N<-initial.values[[ichain]]$N
        } else {
          stop(paste0("initial value for N for chain ",ichain," must be positive integer >=",base::sum(inits[[ichain]]$H>1)))
        }
      } else {
        if(mod.p.h){
          pstar <- pstarintegrand(inits[[ichain]]$pbeta,sqrt(inits[[ichain]]$sigma2_zp),DM$p,gq)
        } else {
          pstar<-1-min(1.-tol,max(tol,prod(1-expit(DM$p%*%inits[[ichain]]$pbeta))))
        }
        inits[[ichain]]$N<-base::sum(inits[[ichain]]$H>1)+rnbinom(1,base::sum(inits[[ichain]]$H>1),pstar)
      }      
    }
    if(mms@data.type=="never"){    
      inits[[ichain]]$alpha <- 0.0
    } else if(mms@data.type=="sometimes"){
      if(length(initial.values[[ichain]]$alpha)){
        if(length(initial.values[[ichain]]$alpha)==1 & initial.values[[ichain]]$alpha>0 & initial.values[[ichain]]$alpha<1){
          inits[[ichain]]$alpha <- initial.values[[ichain]]$alpha
        } else {
          stop(paste("chain",ichain,"initial value for alpha must be a scalar between 0 and 1"))
        }
      } else {
        inits[[ichain]]$alpha <- rbeta(1,a0alpha,b0alpha)
      }    
    } else {
      inits[[ichain]]$alpha <- 1.0
    }
    if(length(initial.values[[ichain]]$delta_1) & length(initial.values[[ichain]]$delta_2)){
      if(length(initial.values[[ichain]]$delta_1)==1 & length(initial.values[[ichain]]$delta_2)==1 & (dunif(initial.values[[ichain]]$delta_1+initial.values[[ichain]]$delta_2))){
        inits[[ichain]]$delta_1 <- initial.values[[ichain]]$delta_1
        inits[[ichain]]$delta_2 <- initial.values[[ichain]]$delta_2
      } else {
        stop("initial values for delta_1 and delta_2 must be positive scalars with sum less than 1")
      }
    } else if(length(initial.values[[ichain]]$delta_1) | length(initial.values[[ichain]]$delta_2)){
      stop("initial values for delta_1 and delta_2 must be positive scalars with sum less than 1")    
    } else {
      delta<-rdirichlet(1,a0delta)
      inits[[ichain]]$delta_1<-delta[1]
      inits[[ichain]]$delta_2<-delta[2]    
    }   
    if(length(initial.values[[ichain]]$psi)){
      if(length(initial.values[[ichain]]$psi)==1 & initial.values[[ichain]]$psi>0 & initial.values[[ichain]]$psi<1){
        inits[[ichain]]$psi<-initial.values[[ichain]]$psi
      } else {
        stop("initial value for psi must be a scalar between 0 and 1")
      }
    } else {
      inits[[ichain]]$psi <- rbeta(1,tol+base::sum(inits[[ichain]]$H>1),1+M-base::sum(inits[[ichain]]$H>1))
    }
  }
  return(inits)
}

get_known<-function(known,Enc.Mat,naivex,data.type){
  M <- nrow(Enc.Mat)
  if(length(known) & base::sum(known)>0){
    if(length(known)!=M | base::sum(known)>M){
      stop(paste0("'known' must be an integer vector of length ",M," with sum between 0 and ",M))
    } else {
      knownx <- getfreq(Enc.Mat[which(known>0),],data.type)
      if(base::sum(apply(Enc.Mat==3 | Enc.Mat==4,1,base::sum)>0)>base::sum(knownx)) stop("'known' vector misspecified. Encounter histories containing simultaneous encounters are known")     
    }
  } else {
    knownx <- integer(length(naivex))
  }
  knownx
}

#' Generate model inputs for fitting 'multimark' models
#'
#' This function generates an object of class \code{multimarksetup} that is required to fit `multimark' models. 
#'
#'
#' @param Enc.Mat A matrix of observed encounter histories with rows corresponding to individuals and columns corresponding to sampling occasions (ignored unless \code{mms=NULL}).
#' @param data.type Specifies the encounter history data type. All data types include non-detections (type 0 encounter), type 1 encounter (e.g., left-side), and type 2 encounters (e.g., right-side). When both type 1 and type 2 encounters occur for the same individual within a sampling occasion, these can either be "non-simultaneous" (type 3 encounter) or "simultaneous" (type 4 encounter). Three data types are currently permitted:
#' 
#'  \code{data.type="never"} indicates both type 1 and type 2 encounters are never observed for the same individual within a sampling occasion, and observed encounter histories therefore include only type 1 or type 2 encounters (e.g., only left- and right-sided photographs were collected). Observed encounter histories can consist of non-detections (0), type 1 encounters (1), and type 2 encounters (2). See \code{\link{bobcat}}. Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 3 encounters (3).
#'
#'  \code{data.type="sometimes"} indicates both type 1 and type 2 encounters are sometimes observed (e.g., both-sided photographs are sometimes obtained, but not necessarily for all individuals). Observed encounter histories can consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), type 3 encounters (3), and type 4 encounters (4). Type 3 encounters can only be observed when an individual has at least one type 4 encounter. Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), type 3 encounters (3), and type 4 encounters (4). 
#'
#'  \code{data.type="always"} indicates both type 1 and type 2 encounters are always observed, but some encounter histories may still include only type 1 or type 2 encounters. Observed encounter histories can consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 4 encounters (4). Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 4 encounters (4).
#'
#' @param covs A data frame of temporal covariates for detection probabilities (ignored unless \code{mms=NULL}). The number of rows in the data frame must equal the number of sampling occasions. Covariate names cannot be "time", "age", or "h"; these names are reserved for temporal, behavioral, and individual effects when specifying \code{mod.p} and \code{mod.phi}.
#' @param known Optional integer vector indicating whether the encounter history of an individual is known with certainty (i.e., the observed encounter history is the true encounter history). Encounter histories with at least one type 4 encounter are automatically assumed to be known, and \code{known} does not need to be specified unless there exist encounter histories that do not contain a type 4 encounter that happen to be known with certainty. If specified, \code{known = c(v_1,v_2,...,v_M)} must be a vector of length \code{M = nrow(Enc.Mat)} where \code{v_i = 1} if the encounter history for individual \code{i} is known (\code{v_i = 0} otherwise). Note that known all-zero encounter histories (e.g., `000') are ignored.
#' @param divBasis Integer scaler for dividing up the matrix calculations when determining the possible set of basis vectors for updating the latent encounter histories (based on \code{data.type} and number of sampling occasions). Default is \code{divBasis=100}. For very large problems, increasing \code{divBasis} can help reduce memory requirements and speed up computations. For smaller problems, reducing \code{divBasis} can speed up computation.
#' @param divredBasis Integer scaler for dividing up the matrix calculations when reducing the possible set of basis vectors to those that are needed (based on \code{Enc.Mat}) for updating the latent encounter histories. Default is \code{divredBasis=100}. For very large problems, increasing \code{divredBasis} can help reduce memory requirements and speed up computations. For smaller problems, reducing \code{divredBasis} can speed up computation.
#' 
#' @return An object of class \code{multimarksetup}.
#' @author Brett T. McClintock
#' @seealso \code{\link{multimarksetup-class}}, \code{\link{multimarkClosed}}, \code{\link{bobcat}}
#' @references
#' Bonner, S. J., and Holmberg J. 2013. Mark-recapture with multiple, non-invasive marks. \emph{Biometrics} 69: 766-775.
#' 
#' McClintock, B. T., Conn, P. B., Alonso, R. S., and Crooks, K. R. 2013. Integrated modeling of bilateral photo-identification data in mark-recapture analyses. \emph{Ecology} 94: 1464-1471.
#' @examples
#' \donttest{
#' #Generate object of class "multimarksetup"
#' setup <- processdata(bobcat)
#' 
#' #Run two parallel chains using the default model for bobcat data
#' bobcat.dot<-multimarkClosed(mms=setup,nchains=2)
#' 
#' #Run two parallel chains for bobcat data with temporal effects (i.e., mod.p=~time)
#' bobcat.time <- multimarkClosed(mms=setup,mod.p=~time,nchains=2)}
processdata<-function(Enc.Mat,data.type="never",covs=data.frame(),known=integer(),divBasis=100,divredBasis=100){
  
  if(!is.matrix(Enc.Mat)) stop("'Enc.Mat' must be a matrix")
  
  if(data.type=="never"){
    if(!all(match(unique(c(Enc.Mat)),c(0,1,2,3),nomatch=0))) stop("Encounter histories for 'never' data type can only include 0, 1, 2, and 3 entries")
    if(any(Enc.Mat==3)) warning("Data type is 'never' but includes type 3 encounters")
    ind<-4
  } else if(data.type=="always"){
    if(!all(match(unique(c(Enc.Mat)),c(0,1,2,4),nomatch=0))) stop("Encounter histories for 'always' data type can only include 0, 1, 2, and 4 entries")
    if(!any(Enc.Mat==4)) warning("Encounter histories contain no simulataneous encounters -- should you be using the 'never' data type?")
    ind<-4
  } else if(data.type=="sometimes"){
    if(!all(match(unique(c(Enc.Mat)),c(0,1,2,3,4),nomatch=0))) stop("Encounter histories for 'sometimes' data type can only include 0, 1, 2, 3, and 4 entries")
    temp.check <- which(rowSums(Enc.Mat==3)>0 & rowSums(Enc.Mat==4)==0)
    if(length(temp.check)) warning(paste("Encounter history",temp.check,"includes a type 3 encounter but no type 4 encounter. \n"))
    if(!any(Enc.Mat==4)) warning("Encounter histories contain no simulataneous encounters -- should you be using the 'never' data type?")
    ind<-5
  } else {
    stop("Data type ('data.type') must be 'never', 'sometimes', or 'always'")
  }
  
  M<-nrow(Enc.Mat)
  noccas<-ncol(Enc.Mat)
  
  if(!is.data.frame(covs)) stop("covariates ('covs') must be a data frame")
  
  if(length(covs)){
    if(nrow(covs)!=noccas){
      stop(paste("covariates (covs) must contain an entry for each occasion"))
    }
    nonames <- c("group","time","Time","c","h")
    if(any(match(colnames(covs),nonames,nomatch=0)>0)) stop(paste0("'",nonames[match(colnames(covs),nonames,nomatch=0)],"' cannot be used for covariate ('covs') names. "))
  }
  
  naivex<-getfreq(Enc.Mat,data.type)
  A<- get_A(noccas,data.type)
  gc()
  Basis<-get_basis_vectors(noccas,A$Aprime,A$ivect,div=divBasis,data.type=data.type)
  gc()
  C<-get_C(noccas,data.type=data.type)
  L<-get_L(noccas,data.type=data.type)
  redBasis<-reduced_basis(noccas,Basis,naivex,div=divredBasis,data.type=data.type)[,-1]
  gc()
  ncolbasis<-ncol(redBasis)
  indBasis<-as.vector(which(redBasis!=0)-(ind^noccas)*rep(seq(0,ncolbasis-1),each=3),mode="integer")
  knownx<-get_known(known,Enc.Mat,naivex,data.type)
  mms<-new(Class="multimarksetup",Enc.Mat=Enc.Mat,data.type=data.type,vAll.hists=A$vAll.hists,Aprime=A$Aprime,indBasis=indBasis,ncolbasis=ncolbasis,knownx=knownx,C=C,L=L,naivex=naivex,covs=covs)  
  return(mms)
}

checkvecs<-function(vec,dim,parm){
  if(length(vec)==1){
    vec <- rep(vec,dim) 
  } else if(length(vec)!=dim){
    stop(paste(parm,"must be scaler or vector of length",dim))
  }
  vec
}

checkmats<-function(mat,dim,parm){
  if(!is.matrix(mat)){
    if(length(mat)==1){
      mat <- diag(mat,dim) 
    }
  } else if(!all(dim(mat)==dim)){
    stop(paste(parm,"must be scaler or matrix of dimension",dim,"x",dim))
  }
  mat
}

inverseXB<-function(ichain,iout,betanames,mod.h,DM,noccas,varind,vars,parm,sigparm,link){
  
  if(!varind){
    beta<-as.matrix(iout$mcmc[[ichain]][,betanames])   
  } else {
    beta<-matrix(iout$mcmc[[ichain]][betanames],nrow=1) 
  }
  
  if(!mod.h){
    if(link=="logit"){
      invbeta<-matrix(expit(apply(beta,1,function(x) DM%*%x)),byrow=T,ncol=noccas)
    } else if(link=="probit"){
      invbeta<-matrix(pnorm(apply(beta,1,function(x) DM%*%x)),byrow=T,ncol=noccas)  
    } else {
      stop("link function must be 'logit' or 'probit'")
    }
  } else {
    if(!any(match(sigparm,vars,nomatch=0))) stop(paste(sigparm,"parameter not found"))
    if(!varind){
      sigma2<-as.matrix(iout$mcmc[[ichain]][,sigparm])
    } else {
      sigma2<-matrix(iout$mcmc[[ichain]][sigparm],nrow=1)  
    }
    cat("Performing numerical integration over ",parm," individual effects for chain ",ichain,". This might take a while. \n",sep="")
    
    invbeta<-matrix(0,nrow=length(sigma2),ncol=noccas)
    pb <- txtProgressBar(min=0,max=noccas,char="+",width=100,style=3)
    setTxtProgressBar(pb, 1)
    if(link=="logit"){
      for(j in 1:noccas){  
        invbeta[,j] <- apply(cbind(beta,sigma2),1,function(x) integrate(integrandlogit,-Inf,Inf,beta=x,DM=DM[j,])$value)
        setTxtProgressBar(pb, j)
      }
    } else if(link=="probit"){
      for(j in 1:noccas){  
        invbeta[,j] <- apply(cbind(beta,sigma2),1,function(x) integrate(integrandprobit,-Inf,Inf,beta=x,DM=DM[j,])$value)
        setTxtProgressBar(pb, j)
      }      
    } else {
      stop("link function must be 'logit' or 'probit'")
    }
    close(pb)
  }
  invbeta
}

rinvgamma <- function(n, shape, scale){
  return(1/rgamma(n = n, shape = shape, rate = scale))
}

dinvgamma<-function(x, shape, scale){
  if (shape <= 0 | scale <= 0) stop("Shape or scale parameter negative in dinvgamma().\n")
  alpha <- shape
  beta <- scale
  log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - (beta/x)
  return(log.density)
}

get_missingparms <- function(parms,multiparms){
  multiparms[-match(parms,multiparms)]
}

derivedlogitfun<-function(parms,var){
  if(any(parms==var)){
    getlogit<-function(h,DM,beta,sig2){ 
      if(h){
        temp <- apply(DM,1,function(x) integrate(integrandlogit,-Inf,Inf,beta=c(beta,sig2),DM=x)$value)
      } else {
        temp <- expit(DM %*% beta)
      }
      return(temp)
    }    
  } else {
    getlogit <- function(h,DM,beta,sig2){
      NULL
    }
  }
  return(getlogit)
}

extractmissingparms<-function(missingparms,parm){
  lapply(missingparms,function(x) unlist(x,use.names=FALSE)[which(substr(unlist(x),1,nchar(parm))==parm)])
}

derivedprobitfun<-function(parms,var){
  if(any(parms==var)){
    getprobit<-function(h,DM,beta,sig2){ 
      if(h){
        temp <- apply(DM,1,function(x) integrate(integrandprobit,-Inf,Inf,beta=c(beta,sig2),DM=x)$value)
      } else {
        temp <- pnorm(DM %*% beta)
      }
      return(temp)
    }    
  } else {
    getprobit <- function(h,DM,beta,sig2){
      NULL
    }
  }
  return(getprobit)
}