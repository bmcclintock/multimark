#' Simulate spatially-explicit capture-mark-recapture data from a (demographically) closed population with multiple non-invasive marks
#'
#' This function generates encounter histories from spatially-explicit capture-mark-recapture data consisting of multiple non-invasive marks. 
#'
#'
#' @param N True population size or abundance.
#' @param ntraps The number of traps. If \code{trapCoords=NULL}, the square root of \code{ntraps} must be a whole number in order to create a regular grid of trap coordinates on a square.
#' @param noccas Scaler indicating the number of sampling occasions per trap.
#' @param pbeta Complementary loglog-scale intercept term for detection probability (p). Must be a scaler or vector of length \code{noccas}.
#' @param tau Additive complementary loglog-scale behavioral effect term for recapture probability (c).
#' @param sigma2_scr Complementary loglog-scale term for effect of distance in the ``half-normal'' detection function. Ignored unless \code{detection=``half-normal''}.
#' @param lambda Complementary loglog-scale term for effect of distance in the ``exponential'' detection function. Ignored unless \code{detection=``exponential''}.
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
#' @param detection Model for detection probability as a function of distance from activity centers. Must be "\code{half-normal}" (of the form \eqn{\exp{(-d^2 / (2*\sigma^2))}}, where \eqn{d} is distance) or "\code{exponential}" (of the form \eqn{\exp{(-d / \lambda)}}).
#' @param spatialInputs A list of length 3 composed of objects named \code{trapCoords}, \code{studyArea}, and \code{centers}:
#' 
#'  \code{trapCoords} is a matrix of dimension \code{ntraps} x (2 + \code{noccas}) indicating the Cartesian coordinates and operating occasions for the traps, where rows correspond to trap, the first column the x-coordinate (``x''), and the second column the y-coordinate (``y''). The last \code{noccas} columns indicate whether or not the trap was operating on each of the occasions, where `1' indciates the trap was operating and `0' indicates the trap was not operating.  
#'
#'  \code{studyArea} is a 3-column matrix defining the study area and available habitat. Each row corresponds to a grid cell. The first 2 columns (``x'' and ``y'') indicate the Cartesian x- and y-coordinate for the centroid of each grid cell, and the third column (``avail'') indicates whether the cell is available habitat (=1) or not (=0). All grid cells must have the same resolution. Note that rows should be ordered in raster cell order (raster cell numbers start at 1 in the upper left corner, and increase from left to right, and then from top to bottom).
#'   
#'  \code{centers} is a \code{N}-vector indicating the grid cell (i.e., the row of \code{studyArea}) that contains the true (latent) activity centers for each individual in the population. 
#'
#' If \code{spatialInputs=NULL} (the default), then all traps are assumed to be operating on all occasions, the study area is assumed to be composed of \code{ncells} grid cells, grid cells within \code{buffer} of the trap array are assumed to be available habitat, and the activity centers are randomly assigned to grid cells of available habitat.
#' @param buffer A scaler indicating the buffer around the bounding box of \code{trapCoords} for defining the study area and available habitat when \code{spatialInputs=NULL}.  Default is \code{buffer=3*sqrt(sigma2_scr)}. Ignored unless \code{spatialInputs=NULL}.
#' @param ncells The number of grid cells in the study area when \code{studyArea=NULL}. The square root of \code{ncells} must be a whole number. Default is \code{ncells=1024}. Ignored unless \code{spatialInputs=NULL}.
#' @param scalemax Upper bound for grid cell centroid x- and y-coordinates. Default is \code{scalemax=10}, which scales the x- and y-coordinates to be between 0 and 10.  Ignored unless \code{spatialInputs=NULL}.
#' @param plot Logical indicating whether to plot the simulated trap coordinates, study area, and activity centers using \code{\link{plotSpatialData}}.  Default is \code{plot=TRUE}
#' @details Please be very careful when specifying your own \code{spatialInputs}; \code{\link{multimarkClosedSCR}} and \code{\link{markClosedSCR}} do little to verify that these make sense during model fitting.  
#'
#' @return A list containing the following:
#' \item{Enc.Mat}{Matrix containing the observed encounter histories with rows corresponding to individuals and (\code{ntraps}*\code{noccas}) columns corresponding to traps and sampling occasions.  The first \code{noccas} columns correspond to trap 1, the second \code{noccas} columns corresopond to trap 2, etc.}
#' \item{trueEnc.Mat}{Matrix containing the true (latent) encounter histories with rows corresponding to individuals and (\code{ntraps}*\code{noccas}) columns corresponding to traps and sampling occasions.  The first \code{noccas} columns correspond to trap 1, the second \code{noccas} columns corresopond to trap 2, etc.}
#' \item{spatialInputs}{List of length 2 with objects named \code{trapCoords} and \code{studyArea}:
#' 
#'  \code{trapCoords} is a matrix of dimension \code{ntraps} x (2 + \code{noccas}) indicating the Cartesian coordinates and operating occasions for the traps, where rows correspond to trap, the first column the x-coordinate, and the second column the y-coordinate. The last \code{noccas} columns indicate whether or not the trap was operating on each of the occasions, where `1' indciates the trap was operating and `0' indicates the trap was not operating.
#'
#'  \code{studyArea} is a 3-column matrix containing the coordinates for the centroids a contiguous grid of cells that define the study area and available habitat. Each row corresponds to a grid cell. The first 2 columns indicate the Cartesian x- and y-coordinate for the centroid of each grid cell, and the third column indicates whether the cell is available habitat (=1) or not (=0). All cells must have the same resolution.}
#' \item{centers}{\code{N}-vector indicating the grid cell (i.e., the row of \code{spatialInputs$studyArea}) that contains the true (latent) activity centers for each individual in the population.}
#' @author Brett T. McClintock 
#' @seealso \code{\link{processdataSCR}}, \code{\link{multimarkClosedSCR}}, \code{\link{markClosedSCR}}
#' @references
#' Bonner, S. J., and Holmberg J. 2013. Mark-recapture with multiple, non-invasive marks. \emph{Biometrics} 69: 766-775.
#' 
#' McClintock, B. T., Conn, P. B., Alonso, R. S., and Crooks, K. R. 2013. Integrated modeling of bilateral photo-identification data in mark-recapture analyses. \emph{Ecology} 94: 1464-1471.
#' 
#' Royle, J.A., Karanth, K.U., Gopalaswamy, A.M. and Kumar, N.S. 2009. Bayesian inference in camera trapping studies for a class of spatial capture-recapture models.  \emph{Ecology} 90: 3233-3244.
#' @examples
#' #simulate data for data.type="sometimes" using defaults
#' data<-simdataClosedSCR(data.type="sometimes")
#' 
#' @export
simdataClosedSCR <- function(N=30,ntraps=9,noccas=5,pbeta=0.25,tau=0,sigma2_scr=0.75,lambda=0.75,delta_1=0.4,delta_2=0.4,alpha=0.5,data.type="never",detection="half-normal",spatialInputs=NULL,buffer=3*sqrt(sigma2_scr),ncells=1024,scalemax=10,plot=TRUE){
  
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
    spatialInputs=list()
    if(sqrt(ntraps)%%1) stop("The square root of 'ntraps' must be a whole number")
    if(sqrt(ncells)%%1) stop("The square root of 'ncells' must be a whole number")
    if(buffer>=scalemax/2) stop("'buffer' must be <",scalemax/2)
    spatialInputs$trapCoords<-as.matrix(expand.grid(seq(0+buffer,scalemax-buffer,length=sqrt(ntraps)),seq(0+buffer,scalemax-buffer,length=sqrt(ntraps)))) #trap coordinates on square
    spatialInputs$trapCoords<-cbind(spatialInputs$trapCoords,matrix(1,nrow=ntraps,ncol=noccas)) # assumes all traps are operational on all occasions
    studyArea<-as.matrix(expand.grid(seq(0,scalemax,length=sqrt(ncells)),seq(0,scalemax,length=sqrt(ncells)))) #study area grid
    studyArea <- studyArea[order(-studyArea[,2],studyArea[,1]),] # put in raster cell order (left to right, top left to bottom right)
    NN<-numeric()
    for(i in 1:nrow(spatialInputs$trapCoords)){
      od <- sqrt( (studyArea[,1]-spatialInputs$trapCoords[i,1])^2  +  (studyArea[,2]-spatialInputs$trapCoords[i,2])^2  )
      od <- (1:length(od))[od <= buffer]
      NN<-c(NN,od)
    }
    spatialInputs$studyArea<-cbind(studyArea,rep(0,ncells))
    spatialInputs$studyArea[unique(NN),3]<-1
    if(!sum(spatialInputs$studyArea[,3])) stop("No available habitat; 'buffer' probably too small")
    spatialInputs$centers<-sample.int(nrow(studyArea),N,replace=TRUE,prob=spatialInputs$studyArea[,3])
  } else {
    if(!is.list(spatialInputs) | length(spatialInputs)!=3 | any(sort(names(spatialInputs))!=c("centers","studyArea","trapCoords"))) stop("'spatialInputs' must be a list of length 3 containing the object 'trapCoords', 'studyArea', and 'centers'")
    if(ntraps!=nrow(spatialInputs$trapCoords)){
      warning("'ntraps' not equal to the number of rows of spatialInputs$trapCoords; the latter will be used.")
      ntraps<-nrow(spatialInputs$trapCoords)
    }
    if(noccas!=ncol(spatialInputs$trapCoords[,-c(1,2)])){
      stop("'spatialInputs$trapCoords' must have ",2+noccas," columns")
    }
    if(N!=length(spatialInputs$centers)){
      warning("'N' not equal to the length of spatialInputs$centers; the latter will be used.")
      N<-length(spatialInputs$centers)
    }
  }
  checkSpatialInputs(spatialInputs$trapCoords,spatialInputs$studyArea,centers=spatialInputs$centers)
  
  colnames(spatialInputs$trapCoords)<-c("x","y",paste0("occ",1:noccas))
  rownames(spatialInputs$trapCoords)<-paste0("trap",1:ntraps)
  colnames(spatialInputs$studyArea)<-c("x","y","avail")
  rownames(spatialInputs$studyArea)<-paste0("cell",1:nrow(spatialInputs$studyArea))
  
  trapCoords<-spatialInputs$trapCoords
  centers<-spatialInputs$studyArea[spatialInputs$centers,1:2] 
  
  tEnc.Mat<-matrix(0,nrow=N,ncol=noccas*ntraps)        #"true" latent histories
  
  dist2<-getdist(centers,trapCoords)
  
  if(detection=="half-normal") {
    dexp<-2
    sigma2<-sigma2_scr
  } else if(detection=="exponential"){
    dexp<-1
    sigma2<-lambda
  } else {
    stop("'detection' argument must be 'half-normal' or 'exponential'")
  }
  
  for(i in 1:N){
    for(k in 1:ntraps){
      ind<-0
      for(j in 1:noccas){
        p<-invcloglog(pbeta[j]-1./(dexp*sigma2)*dist2[i,k]^dexp)*trapCoords[k,2+j]
        c<-invcloglog(pbeta[j]+tau-1./(dexp*sigma2)*dist2[i,k]^dexp)*trapCoords[k,2+j]
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
  if(plot) plotSpatialData(tEnc.Mat,trapCoords=spatialInputs$trapCoords,studyArea=spatialInputs$studyArea,centers=spatialInputs$centers,trapLines=TRUE)
  return(list(Enc.Mat=Enc.Mat,trueEnc.Mat=tEnc.Mat,spatialInputs=spatialInputs))
}

#' Plot spatial capture-mark-recapture data
#'
#' This function plots the study area grid, available habitat, and trap coordinates for spatial capture-recapture studies.  Activity centers and capture locations can also be plotted. 
#'
#'
#' @param mms An optional object of class \code{multimarkSCRsetup-class} from which the (re-scaled) study area and trap coordinates are plotted.
#' @param trapCoords A matrix of dimension \code{ntraps} x (2 + \code{noccas}) indicating the Cartesian coordinates and operating occasions for the traps, where rows correspond to trap, the first column the x-coordinate, and the second column the y-coordinate. The last \code{noccas} columns indicate whether or not the trap was operating on each of the occasions, where `1' indciates the trap was operating and `0' indicates the trap was not operating. Ignored unless \code{mms=NULL}.  
#' @param studyArea A 3-column matrix defining the study area and available habitat. Each row corresponds to a grid cell. The first 2 columns indicate the Cartesian x- and y-coordinate for the centroid of each grid cell, and the third column indicates whether the cell is available habitat (=1) or not (=0). All cells must have the same resolution. Ignored unless \code{mms=NULL}.
#' Note that rows should be ordered in raster cell order (raster cell numbers start at 1 in the upper left corner, and increase from left to right, and then from top to bottom).
#' @param centers An optional vector indicating the grid cell (i.e., the row of \code{studyArea}) that contains the true (latent) activity centers for each individual. If \code{mms} is provided, then \code{centers} must be of length \code{nrow(Enc.Mat)} (i.e., a center must be provided for each observed individual). 
#' @param trapLines Logical indicating whether to draw lines from activity centers to respective traps at which each individual was captured. Default is \code{trapLines=FALSE}. Ignored when \code{mms=NULL} or \code{centers=NULL}. 
#' @author Brett T. McClintock 
#' @examples
#' #Plot the tiger example data
#' plotSpatialData(trapCoords=tiger$trapCoords,studyArea=tiger$studyArea)
#' 
#' @export
#' @importFrom graphics hist legend lines par points
#' @importFrom sp points2grid SpatialPoints bbox gridded SpatialGrid over coordinates
#' @importFrom raster extent raster rasterize
#' @importFrom grDevices rainbow gray
#' @importFrom utils capture.output
plotSpatialData<-function(mms=NULL,trapCoords,studyArea,centers=NULL,trapLines=FALSE){
  
  cur.par<-par(no.readonly=TRUE)
  
  if(!is.null(mms)){
    if(inherits(mms,"multimarkSCRsetup")){
      studyArea<-mms@spatialInputs$studyArea
      trapCoords<-mms@spatialInputs$trapCoords
      Enc.Mat<-mms@Enc.Mat
    } else {
      Enc.Mat<-mms
    }
  }
  checkSpatialInputs(trapCoords,studyArea,centers)
  
  gridtol <- checkGridTol(studyArea,warn=FALSE)
  
  if(!isTRUE(all.equal(1:nrow(studyArea),order(-studyArea[,"y"], studyArea[,"x"])))) {
    warning("studyArea rows do not appear to be in raster cell order; these will be reordered accordingly")
    studyArea <- studyArea[order(-studyArea[,"y"], studyArea[,"x"]),] # put in raster cell order (left to right, top left to bottom right)
  }
  
  sAgrid<-sp::SpatialGrid(sp::points2grid(sp::SpatialPoints(studyArea[,c("x","y")]),tolerance = gridtol))
  e<-extent(sAgrid)
  cells.dim<-attributes(sAgrid)$grid@cells.dim
  r<-raster(e,nrow=cells.dim["y"],ncol=cells.dim["x"])
  x<-rasterize(coordinates(sAgrid),r,studyArea[,3],fun=max)
  breakpoints <- c(0,0.5,1)
  colors <- c("white","green","green")
  sp::plot(x,breaks=breakpoints,col=colors,legend=FALSE)
  sp::plot(sAgrid,add=TRUE,col=gray(.75))
  
  if(is.null(centers)){
    tmp<-legend(x="top",legend=c("traps                          ","available habitat"),pch=c(17,15),col=c(1,"green"),box.col="white",bty="o",bg="white",cex=.8)
  } else {
    M<-ifelse(!is.null(mms),nrow(Enc.Mat),length(centers))
    rainColors<-grDevices::rainbow(M)
    ntraps<-nrow(trapCoords)
    noccas<-ncol(trapCoords[,-c(1,2),drop=FALSE])
    if(length(centers)!=M) stop("'centers' must be of length ",M)
    for(i in 1:M){
      points(studyArea[centers[i],"x"],studyArea[centers[i],"y"],pch=20,col=rainColors[i])
      if(!is.null(mms) & trapLines){
        for(k in 1:ntraps){
          for(j in 1:noccas){
            if(Enc.Mat[i,(k-1)*noccas+j]>0){
              wig1<-rnorm(1,0,0.01)
              wig2<-rnorm(1,0,0.01)
              points(studyArea[centers[i],"x"],studyArea[centers[i],"y"],pch=20,col=rainColors[i])
              #points(trapCoords[k,"x"]+wig1,trapCoords[k,"y"]+wig2,2,pch=17)
              lines(rbind(trapCoords[k,c("x","y")]+c(wig1,wig2),studyArea[centers[i],c("x","y")]),col=rainColors[i],lty=2,lwd=.5)
            }
          }
        }
      }
    }
    tmp<-legend(x="top",legend=c("traps                          ","available habitat", "activity centers"),pch=c(17,15,1),col=c(1,"green"),box.col="white",bty="o",bg="white",cex=.8)
  }
  points(trapCoords[,c("x","y")],pch=17)
  
  par(cur.par)
}

pstarintegrandSCR<-function(noccas,beta,sigma2,DM,spatialInputs,dexp){
  msk2 <- which(c(t(spatialInputs$msk))==1)
  XB <- DM[msk2,,drop=FALSE] %*% beta
  dist2 <- spatialInputs$dist2
  detProb<-invcloglogtol(matrix(XB,nrow=dim(dist2)[1],ncol=length(msk2),byrow=TRUE)-1/(dexp*sigma2)*dist2[,rep(1:dim(dist2)[2],times=apply(spatialInputs$msk,1,sum))]^dexp)
  #detProb<-exp(matrix(0,nrow=dim(dist2)[1],ncol=length(msk2),byrow=TRUE)-1/(dexp*sigma2)*dist2[,rep(1:dim(dist2)[2],times=apply(spatialInputs$msk,1,sum))]^dexp)
  pdot<-1.-apply(1.-detProb,1,function(x) max(prod(x),tol))
  esa<-sum(pdot)*spatialInputs$a
  esa/spatialInputs$A
}

loglikeClosedSCR<-function(parms,DM,noccas,ntraps,C,All.hists,spatialInputs){
  
  H <- parms$H
  pbeta <- parms$pbeta
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
  if(DM$mod.det=="half-normal"){
    dexp <- 2
    sigma2 <- parms$sigma2_scr
  } else {
    dexp <-1
    sigma2 <- parms$lambda
  }
  
  Hind <- H[which(H>1)]
  centers <- spatialInputs$studyArea[parms$centers[which(H>1)],]
  indhist <- All.hists[Hind,]
  n<-length(Hind)
  #firstcap<- (C[Hind]>=matrix(rep(1:noccas,each=n),nrow=n,ncol=noccas))
  
  msk <- t(spatialInputs$msk) #matrix(1,nrow=noccas,ncol=ntraps) 
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
  prevcap <- prevcap[msk2==1]   
  indid <- idx[msk2==1,1] 
  repid <- idx[msk2==1,2] 
  trapid <- idx[msk2==1,3] 
  
  trapgridbig <- spatialInputs$trapCoords[trapid,c(1,2)]   
  c1 <- (centers[indid,1] - trapgridbig[,1])^2
  c2 <- (centers[indid,2] - trapgridbig[,2])^2
  
  p <- invcloglogtol(rep(DM$p%*%pbeta,each=n)[msk2==1]*(1-prevcap) + rep(DM$c%*%pbeta,each=n)[msk2==1]*prevcap - 1./(dexp*sigma2)*sqrt(c1+c2)^dexp)
  #p <- exp(- 1./(dexp*sigma2)*sqrt(c1+c2)^dexp)
  #p2 <- invcloglogtol(matrix(rep(DM$p%*%pbeta,each=n)*(1-prevcap2)+rep(DM$c%*%pbeta,each=n)*prevcap2,nrow=n,ncol=noccas*ntraps)- 1./(dexp*sigma2)*(dist2mat)^dexp)
  
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
  
  pstar <- pstarintegrandSCR(noccas,pbeta,sigma2,DM$p,spatialInputs,dexp)
  
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
  
  #priors <- priors + log(2.0*dcauchy(sqrt(parms$sigma2_scr),0.0,priorparms$a,log=FALSE))
  priors <- priors + dunif(sqrt(ifelse(DM$mod.det=="half-normal",parms$sigma2_scr,parms$lambda)),priorparms$sigma_bounds[1],priorparms$sigma_bounds[2],log=TRUE)
  
  priors <- priors + length(parms$centers)*log(1./spatialInputs$A)
  
  priors
}

posteriorClosedSCR<-function(parms,DM,mms,priorparms,spatialInputs){
  nchains<-length(parms)
  ntraps<-nrow(spatialInputs$trapCoords)
  noccas<-ncol(mms@Enc.Mat)/ntraps
  M<-nrow(mms@Enc.Mat)
  All.hists<-matrix(mms@vAll.hists,byrow=TRUE,ncol=noccas*ntraps)
  
  priorparms$sigma_bounds <- priorparms$sigma_bounds / mms@spatialInputs$Srange
  
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

checkGridTol <- function(studyArea,warn=TRUE){
  gridtol <- sqrt(.Machine$double.eps)
  checkTol <- tryCatch(sp::points2grid(sp::SpatialPoints(studyArea[,1:2])),error=function(e) e)
  if(inherits(checkTol,"error")){
    tmp <- utils::capture.output(tryCatch(sp::points2grid(sp::SpatialPoints(studyArea[,1:2])),error=function(e){} ))
    tmp <- tmp[-which(tmp=="NULL")]
    gridtol <- as.numeric(gsub("suggested tolerance minimum:","",tmp))
    if(warn) warning(gsub("Error in ","",checkTol),"   tolerance set to ",gridtol," -- if this number isn't very small then check studyArea coordinates!")
  }
  gridtol
}

checkSpatialInputs<-function(trapCoords,studyArea,centers=NULL){
  
  if(ncol(studyArea)!=3){
    stop("'studyArea' must have 3 columns")
  }
  if(!all(studyArea[,3] %in% c(0,1))){
    stop("Indicators for available habitat in 'studyArea' must be 0 or 1")
  }
  if(!all(trapCoords[,-c(1,2)] %in% c(0,1))){
    stop("Indicators for each occasion in 'trapCoords' must be 0 or 1")
  }
  
  gridtol <- checkGridTol(studyArea,warn=TRUE)
  
  if(!sp::gridded(SpatialGrid(sp::points2grid(sp::SpatialPoints(studyArea[,1:2]),tolerance = gridtol)))){
    stop("'studyArea' must be a regular grid ")
  } 
  
  cellsize<-sp::points2grid(sp::SpatialPoints(studyArea[,1:2]),tolerance = gridtol)@cellsize
  if(!(diff(range(cellsize)) < .Machine$double.eps ^ 0.5)) stop("studyArea grid cells must be square")
  
  if(any(is.na(sp::over(sp::SpatialPoints(trapCoords[,1:2]),sp::SpatialGrid(sp::points2grid(sp::SpatialPoints(studyArea[,1:2]),tolerance = gridtol)))))){
    stop("'trapCoords' must be within 'studyArea'")
  }
  if(!is.null(centers)){
    if(any(centers<1) | any(centers>nrow(studyArea))){
      stop("activity center cell indices must have values between 1 and ",nrow(studyArea))
    }  
    if(any(studyArea[centers,3]==0)){
      stop("activity centers must be in cells of available habitat")
    }
  }
  
}

mapCenters<-function(centers,centermap1,centermap2,origtoavail=TRUE){
  if(origtoavail){
    ncenters <- centers - centermap1[centers]
  } else {
    ncenters <- match(centers,centermap2)
  }
  as.integer(ncenters)
}

checkClosedSCR<-function(parms,parmlist,mms,DM,iter,adapt,bin,thin,burnin,taccept,tuneadjust,maxnumbasis,a0delta,a0alpha,b0alpha,sigma_bounds,sigma2_mu0,a0psi,b0psi){
  
  if(mms@data.type!="sometimes" & any(parms=="alpha")) stop("Parameter 'alpha' only applies to models for the 'sometimes' data type")
  if(DM$mod.det=="half-normal" & any(parms=="lambda")) stop("Parameter 'lambda' only applies to `exponential' detection function")
  if(DM$mod.det=="exponential" & any(parms=="sigma2_scr")) stop("Parameter 'sigma2_scr' only applies to `half-normal' detection function")
  
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
  if(!all(c(a0delta,a0alpha,b0alpha,sigma_bounds,sigma2_mu0,a0psi,b0psi)>0)) stop("'a0delta', 'a0alpha', 'b0alpha', 'sigma_bounds', 'sigma2_mu0', 'a0psi', and 'b0psi' must be >0")
  
  pdim<-ncol(DM$p)
  if(!pdim) stop("'mod.p' must include at least 1 parameter")
  
  params
}

mcmcClosedSCR<-function(ichain,mms,DM,params,inits,iter,adapt,bin,thin,burnin,taccept,tuneadjust,Prop.sd,Prop.center,spatialInputs,maxnumbasis,a0delta,a0alpha,b0alpha,sigma_bounds,mu0,sigma2_mu0,a0psi,b0psi,printlog){
  
  ntraps<-nrow(spatialInputs$trapCoords)
  noccas<-ncol(mms@Enc.Mat)/ntraps
  M<-nrow(mms@Enc.Mat)
  DMp<-DM$p
  DMc<-DM$c
  pdim<-ncol(DMp)
  dexp<-ifelse(DM$mod.det=="half-normal",2,1)
  firstcap<-mms@C
  
  #declare and initialize parameters
  pbeta<-rep(NA,max(1,floor(iter/thin))*(pdim))
  H<-rep(NA,ifelse(any(params=="H"),max(1,floor(iter/thin))*M,M))
  centers<-rep(NA,ifelse(any(params=="centers"),max(1,floor(iter/thin))*M,M)) 
  sigma2_scr<-rep(NA,max(1,floor(iter/thin)))
  alpha<-rep(NA,max(1,floor(iter/thin)))
  delta_1<-rep(NA,max(1,floor(iter/thin)))
  delta_2<-rep(NA,max(1,floor(iter/thin)))
  N<-rep(NA,max(1,floor(iter/thin)))
  psi<-rep(NA,max(1,floor(iter/thin)))
  logPosterior<-rep(NA,max(1,floor(iter/thin)))
  
  pbeta[1:pdim] <- inits[[ichain]]$pbeta
  H[1:M] <- inits[[ichain]]$H-1
  centers[1:M] <- inits[[ichain]]$centers-1
  sigma2_scr[1] <- ifelse(DM$mod.det=="half-normal",inits[[ichain]]$sigma2_scr,inits[[ichain]]$lambda)
  alpha[1] <- inits[[ichain]]$alpha
  delta_1[1] <- inits[[ichain]]$delta_1
  delta_2[1] <- inits[[ichain]]$delta_2
  N[1] <- inits[[ichain]]$N
  psi[1] <- inits[[ichain]]$psi
  
  arate<-numeric(M+pdim+1)
  
  cummind<-c(0,cumsum(apply(spatialInputs$msk,1,sum)))
  mind<-c(unlist(apply(spatialInputs$msk,1,function(x) which(x>0))))-1
  
  posterior <- .C(ClosedSCRC,as.integer(ichain),as.numeric(mu0), as.numeric(sigma2_mu0), as.numeric(pbeta), as.numeric(sigma2_scr), as.numeric(delta_1),as.numeric(delta_2),as.numeric(alpha), as.integer(inits[[ichain]]$x), as.numeric(N), as.numeric(psi), as.integer(H), as.integer(centers),
                  as.integer(ntraps),as.integer(noccas), as.integer(M), as.numeric(a0delta), as.numeric(a0alpha), as.numeric(b0alpha), as.numeric(sigma_bounds/mms@spatialInputs$Srange), as.numeric(a0psi), as.numeric(b0psi),
                  as.numeric(Prop.sd),as.integer(Prop.center$NNvect-1),as.integer(Prop.center$numnn),as.integer(Prop.center$cumnumnn),as.numeric(arate),as.numeric(logPosterior),
                  as.integer(length(mms@vAll.hists)/(noccas*ntraps)),as.integer(mms@vAll.hists), as.integer(firstcap), as.integer(mms@indBasis-1), as.integer(mms@ncolbasis), as.integer(mms@knownx), as.numeric(as.vector(t(DMp))), as.numeric(as.vector(t(DMc))),as.integer(pdim),
                  as.integer(iter), as.integer(thin), as.integer(adapt), as.integer(bin), as.numeric(taccept),as.numeric(tuneadjust),as.integer(maxnumbasis),
                  as.integer(mms@data.type=="sometimes"),as.integer(any(params=="H")),as.integer(any(params=="centers")),as.integer(DM$mod.delta != ~NULL),as.integer(DM$mod.delta==formula(~type)),as.numeric(dexp),as.numeric(spatialInputs$dist2^dexp),as.integer(nrow(spatialInputs$studyArea)),as.numeric(spatialInputs$A),as.integer(c(t(spatialInputs$msk))),as.integer(cummind),as.integer(mind),as.integer(printlog),NAOK = TRUE) 
  
  names(posterior) <- c("ichain","mu_0","sigma2_mu","pbeta","sigma2_scr", "delta_1","delta_2","alpha", "x", "N", "psi","H", "centers", "ntraps", "noccas", "M","a0delta", "a0alpha", "b0alpha","sigma_bounds","a0psi","b0psi","Prop.sd", "NNvect", "numnn","cumnumnn", "arate","logPosterior","nHists","vAll.hists","firstcap", "indBasis", "ncolBasis","knownx","DMp","DMc","pdim","iter", "thin", "adapt", "bin", "taccept","tuneadjust","maxnumbasis","sometimes?","H?","centers?","updatedelta?","type?","dexp","dist2","ncells","Area","msk","cummind","mind","printlog?")
  
  g <- posterior$iter
  x <- posterior$x
  posterior$centers<-mapCenters(posterior$centers+1,spatialInputs$centermap1,spatialInputs$centermap2,0)
  posterior$sigma2_scr <- posterior$sigma2_scr * mms@spatialInputs$Srange^2
  
  if(any(params=="centers")){
    temp<-cbind(matrix(posterior$pbeta[(floor(burnin/thin)*pdim+1):(max(1,floor(iter/thin))*pdim)],ncol=pdim,byrow=T),matrix(posterior$centers[(floor(burnin/thin)*M+1):(max(1,floor(iter/thin))*M)],ncol=M,byrow=T),posterior$sigma2_scr[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$delta_1[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$delta_2[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$alpha[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$N[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$psi[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))]) 
    centers <- NULL
  } else {
    centers <- posterior$centers
    temp<-cbind(matrix(posterior$pbeta[(floor(burnin/thin)*pdim+1):(max(1,floor(iter/thin))*pdim)],ncol=pdim,byrow=T),posterior$sigma2_scr[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$delta_1[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$delta_2[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$alpha[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$N[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))],posterior$psi[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))])       
  }
  if(any(params=="H")){
    posterior<-cbind(temp,matrix(posterior$H[(floor(burnin/thin)*M+1):(max(1,floor(iter/thin))*M)]+1,ncol=M,byrow=T),posterior$logPosterior[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))]) 
    H <- NULL
  } else {
    H <- posterior$H+1
    posterior<-cbind(temp,posterior$logPosterior[(floor(burnin/thin)+1):(max(1,floor(iter/thin)))])       
  }
  return(list(posterior=posterior,x=x,H=H,centers=centers,g=g))
}

processClosedSCRchains<-function(chains,params,DM,M,nchains,iter,burnin,thin){
  
  parms<-params
  if(any(parms=="pbeta")){
    parms<-c(paste0("pbeta[",colnames(DM$p),"]"),params[which(params!="pbeta")])
  }
  if(any(parms=="delta")){
    if(DM$mod.delta==formula(~type)){
      deltaname<-c("delta_1","delta_2")   
      parms<-c(parms[which(parms!="delta")],"delta_1","delta_2") 
    } else {
      deltaname<-c("delta")   
      parms<-c(parms[which(parms!="delta")],"delta_1") 
    }
  } else {
    deltaname<-NULL
  }
  if(any(parms=="centers")){
    centersname<-paste0("center[",1:M,"]")
    parms<-c(centersname,parms[which(parms!="centers")])
  } else {
    centersname<-NULL
  }
  if(any(parms=="H")){
    Hname<-paste0("H[",1:M,"]")
    parms<-c(Hname,parms[which(parms!="H")])
  } else {
    Hname<-NULL
  }
  sig2name<-ifelse(DM$mod.det=="half-normal","sigma2_scr","lambda")
  
  initial.values <- list()
  
  for(ichain in 1:nchains){
    checkend <- chains[[ichain]]$g
    if(checkend<iter | !is.finite(chains[[ichain]]$posterior[nrow(chains[[ichain]]$posterior),ncol(chains[[ichain]]$posterior)])) {
      warning(paste0("chain ",ichain," terminated at iteration ",checkend,"; check log for more information"))
      if(!checkend & burnin<1){
        initstemp <- chains[[ichain]]$posterior[1,]
      } else if(floor(checkend/thin)>floor(burnin/thin)){
        initstemp <- chains[[ichain]]$posterior[floor(checkend/thin)-floor(burnin/thin),]  
      } else {
        initstemp <- chains[[ichain]]$posterior[nrow(chains[[ichain]]$posterior),]
      }
    } else {
      initstemp <- chains[[ichain]]$posterior[nrow(chains[[ichain]]$posterior),]
    }
    names(initstemp) <- c(paste0("pbeta[",colnames(DM$p),"]"),centersname,sig2name,"delta_1","delta_2","alpha","N","psi",Hname,"logPosterior")
    if(any(params=="centers")){
      initial.values[[ichain]] <- list(pbeta=initstemp[paste0("pbeta[",colnames(DM$p),"]")],centers=initstemp[centersname],sigma2_scr=initstemp[sig2name],delta_1=initstemp["delta_1"],delta_2=initstemp["delta_2"],alpha=initstemp["alpha"],N=initstemp["N"],psi=initstemp["psi"],x=chains[[ichain]]$x,H=chains[[ichain]]$H)
    } else {
      initial.values[[ichain]] <- list(pbeta=initstemp[paste0("pbeta[",colnames(DM$p),"]")],centers=chains[[ichain]]$centers,sigma2_scr=initstemp[sig2name],delta_1=initstemp["delta_1"],delta_2=initstemp["delta_2"],alpha=initstemp["alpha"],N=initstemp["N"],psi=initstemp["psi"],x=chains[[ichain]]$x,H=chains[[ichain]]$H)
      names(initial.values[[ichain]]$centers) <- paste0("center[",1:M,"]")
    }
    if(any(params=="H")){
      initial.values[[ichain]]$H <- initstemp[Hname]
    } else {
      initial.values[[ichain]]$H <- chains[[ichain]]$H
      names(initial.values[[ichain]]$H) <- paste0("H[",1:M,"]")
    }
    names(initial.values[[ichain]]$x) <- paste0("x[",1:length(initial.values[[ichain]]$x),"]")
    names(initial.values[[ichain]])[match("sigma2_scr",names(initial.values[[ichain]]))]<-sig2name
    chains[[ichain]] <- chains[[ichain]]$posterior
    colnames(chains[[ichain]]) <- names(initstemp)   
    chains[[ichain]] <- chains[[ichain]][,parms]
    if(!is.null(deltaname)){
      if(!is.null(nrow(chains[[ichain]]))) {
        colnames(chains[[ichain]])[which(substr(colnames(chains[[ichain]]),1,nchar("delta"))=="delta")] <- deltaname
      } else {
        names(chains[[ichain]])[which(substr(names(chains[[ichain]]),1,nchar("delta"))=="delta")] <- deltaname     
      }
    }
    chains[[ichain]] <- mcmc(chains[[ichain]],start=1,thin=1)
    if(burnin<thin){
      temp=seq(thin,max(1,iter),thin)
    } else {
      temp=seq(thin*(floor(burnin/thin)+1),iter,thin)
    }
    attributes(chains[[ichain]])$mcpar <- c(head(temp,n=1),tail(temp,n=1),thin)  
  }
  chains<-as.mcmc.list(chains)
  return(list(chains=chains,initial.values=initial.values))  
}

getSpatialInputs<-function(mms){
  spatialInputs=list()
  spatialInputs$studyArea <- mms@spatialInputs$studyArea[which(mms@spatialInputs$studyArea[,"avail"]==1),c("x","y")]  # available habitat study area
  spatialInputs$trapCoords <- mms@spatialInputs$trapCoords[,c(1,2)]
  gridtol <- checkGridTol(mms@spatialInputs$studyArea,warn=FALSE)
  spatialInputs$a <- sp::points2grid(sp::SpatialPoints(mms@spatialInputs$studyArea[,1:2]),tolerance = gridtol)@cellsize[1]
  spatialInputs$A <- spatialInputs$a * nrow(spatialInputs$studyArea)
  spatialInputs$dist2 <- getdist(spatialInputs$studyArea,spatialInputs$trapCoords)
  spatialInputs$msk <- mms@spatialInputs$trapCoords[,-c(1,2),drop=FALSE]
  spatialInputs$centermap1 <- cumsum(mms@spatialInputs$studyArea[,"avail"]==0)
  spatialInputs$centermap2 <- cumsum(mms@spatialInputs$studyArea[,"avail"]==1)
  spatialInputs
}

getPropCenter<-function(spatialInputs,propcenter){
  NN<-list()
  RAD <- ifelse(is.null(propcenter),spatialInputs$a*10,propcenter) # Change propcenter to get more or fewer neighbors
  for(i in 1:nrow(spatialInputs$studyArea)){
    od <- sqrt( (spatialInputs$studyArea[i,1]-spatialInputs$studyArea[,1])^2  +  (spatialInputs$studyArea[i,2]-spatialInputs$studyArea[,2])^2  )
    od <- (1:length(od))[od < RAD]
    NN[[i]]<-od
  }
  Prop.center<-list(NNvect=unlist(NN),numnn=unlist(lapply(NN,length)))
  Prop.center$cumnumnn<-c(0,cumsum(Prop.center$numnn))[1:length(Prop.center$numnn)]
  Prop.center
}

#' Fit spatial population abundance models for ``traditional'' capture-mark-recapture data consisting of a single mark type
#'
#' This function fits spatial population abundance models for ``traditional'' capture-mark-recapture data consisting of a single mark type using Bayesian analysis methods. Markov chain Monte Carlo (MCMC) is used to draw samples from the joint posterior distribution. 
#'
#'
#' @param Enc.Mat A matrix containing the observed encounter histories with rows corresponding to individuals and (\code{ntraps}*\code{noccas}) columns corresponding to traps and sampling occasions.  The first \code{noccas} columns correspond to trap 1, the second \code{noccas} columns corresopond to trap 2, etc.
#' @param trapCoords A matrix of dimension \code{ntraps} x (2 + \code{noccas}) indicating the Cartesian coordinates and operating occasions for the traps, where rows correspond to trap, the first column the x-coordinate (``x''), and the second column the y-coordinate (``y''). The last \code{noccas} columns indicate whether or not the trap was operating on each of the occasions, where `1' indciates the trap was operating and `0' indicates the trap was not operating. Ignored unless \code{mms=NULL}.
#' @param studyArea is a 3-column matrix containing the coordinates for the centroids a contiguous grid of cells that define the study area and available habitat. Each row corresponds to a grid cell. The first 2 columns (``x'' and ``y'') indicate the Cartesian x- and y-coordinate for the centroid of each grid cell, and the third column (``avail'') indicates whether the cell is available habitat (=1) or not (=0). All cells must have the same resolution. If \code{studyArea=NULL} (the default) and  \code{mms=NULL}, then a square study area grid composed of \code{ncells} cells of available habitat is drawn around the bounding box of \code{trapCoords} based on \code{buffer}. Ignored unless \code{mms=NULL}.
#' Note that rows should be ordered by raster cell order (raster cell numbers start at 1 in the upper left corner, and increase from left to right, and then from top to bottom).
#' @param buffer A scaler in same units as \code{trapCoords} indicating the buffer around the bounding box of \code{trapCoords} for defining the study area when \code{studyArea=NULL}.  Ignored unless \code{studyArea=NULL}.
#' @param ncells The number of grid cells in the study area when \code{studyArea=NULL}. The square root of \code{ncells} must be a whole number. Default is \code{ncells=1024}. Ignored unless \code{studyArea=NULL} and \code{mms=NULL}.
#' @param covs A data frame of time- and/or trap-dependent covariates for detection probabilities (ignored unless \code{mms=NULL}). The number of rows in the data frame must equal the number of traps times the number of sampling occasions (\code{ntraps*noccas}), where the first \code{noccas} rows correspond to trap 1, the \code{noccas} rows correspond to trap 2, etc. Covariate names cannot be "time", "age", or "h"; these names are reserved for temporal, behavioral, and individual effects when specifying \code{mod.p} and \code{mod.phi}.
#' @param mod.p Model formula for detection probability. For example, \code{mod.p=~1} specifies no effects (i.e., intercept only), \code{mod.p~time} specifies temporal effects, \code{mod.p~c} specifies behavioral reponse (i.e., trap "happy" or "shy"), \code{mod.p~trap} specifies trap effects, and \code{mod.p~time+c} specifies additive temporal and behavioral effects.
#' @param detection Model for detection probability as a function of distance from activity centers . Must be "\code{half-normal}" (of the form \eqn{\exp{(-d^2 / (2*\sigma^2))}}, where \eqn{d} is distance) or "\code{exponential}" (of the form \eqn{\exp{(-d / \lambda)}}).
#' @param parms A character vector giving the names of the parameters and latent variables to monitor. Possible parameters are cloglog-scale detection probability parameters ("\code{pbeta}"), population abundance ("\code{N}"), and cloglog-scale distance term for the detection function ("\code{sigma2_scr}" when \code{detection=``half-normal''} or "\code{lambda}" when \code{detection=``exponential''}). Individual activity centers ("\code{centers}") and the log posterior density ("\code{logPosterior}") may also be monitored. Setting \code{parms="all"} monitors all possible parameters and latent variables.
#' @param nchains The number of parallel MCMC chains for the model.
#' @param iter The number of MCMC iterations.
#' @param adapt The number of iterations for proposal distribution adaptation. If \code{adapt = 0} then no adaptation occurs.
#' @param bin Bin length for calculating acceptance rates during adaptive phase (\code{0 < bin <= iter}).
#' @param thin Thinning interval for monitored parameters.
#' @param burnin Number of burn-in iterations (\code{0 <= burnin < iter}).
#' @param taccept Target acceptance rate during adaptive phase (\code{0 < taccept <= 1}). Acceptance rate is monitored every \code{bin} iterations. Default is \code{taccept = 0.44}.
#' @param tuneadjust Adjustment term during adaptive phase (\code{0 < tuneadjust <= 1}). If acceptance rate is less than \code{taccept}, then proposal term (\code{proppbeta} or \code{propsigma}) is multiplied by \code{tuneadjust}. If acceptance rate is greater than or equal to \code{taccept}, then proposal term is divided by \code{tuneadjust}. Default is \code{tuneadjust = 0.95}.
#' @param proppbeta Scaler or vector (of length k) specifying the initial standard deviation of the Normal(pbeta[j], proppbeta[j]) proposal distribution. If \code{proppbeta} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{proppbeta = 0.1}.
#' @param propsigma Scaler specifying the initial Gamma(shape = 1/\code{propsigma}, scale = sigma_scr * \code{propsigma}) proposal distribution for sigma_scr = sqrt(sigma2_scr). Default is \code{propsigma=1}.
#' @param propcenter Scaler specifying the neighborhood distance when proposing updates to activity centers. When \code{propcenter=NULL} (the default), then propcenter = a*10, where a is the cell size for the study area grid, and each cell has (at most) approximately 300 neighbors. 
#' @param sigma_bounds Positive vector of length 2 for the lower and upper bounds for the [sigma_scr] ~ Uniform(sigma_bounds[1], sigma_bounds[2]) (or [sqrt(lambda)] when \code{detection=``exponential''}) prior for the detection function term sigma_scr = sqrt(sigma2_scr) (or sqrt(lambda)). When \code{sigma_bounds = NULL} (the default), then \code{sigma_bounds = c(1.e-6,max(diff(range(studyArea[,"x"])),diff(range(studyArea[,"y"]))))}.
#' @param mu0 Scaler or vector (of length k) specifying mean of pbeta[j] ~ Normal(mu0[j], sigma2_mu0[j]) prior. If \code{mu0} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{mu0 = 0}.
#' @param sigma2_mu0 Scaler or vector (of length k) specifying variance of pbeta[j] ~ Normal(mu0[j], sigma2_mu0[j]) prior. If \code{sigma2_mu0} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{sigma2_mu0 = 1.75}.
#' @param initial.values Optional list of \code{nchain} list(s) specifying initial values for "\code{pbeta}", "\code{N}", "\code{sigma2_scr}", and "\code{centers}". Default is \code{initial.values = NULL}, which causes initial values to be generated automatically.
#' @param scalemax Upper bound for internal re-scaling of grid cell centroid coordinates. Default is \code{scalemax=10}, which re-scales the centroids to be between 0 and 10.  Re-scaling is done internally to avoid numerical overflows during model fitting.
#' @param printlog Logical indicating whether to print the progress of chains and any errors to a log file in the working directory. Ignored when \code{nchains=1}. Updates are printed to log file as 1\% increments of \code{iter} of each chain are completed. With >1 chains, setting \code{printlog=TRUE} is probably most useful for Windows users because progress and errors are automatically printed to the R console for "Unix-like" machines (i.e., Mac and Linux) when \code{printlog=FALSE}. Default is \code{printlog=FALSE}.
#' @param ... Additional "\code{parameters}" arguments for specifying \code{mod.p}. See \code{\link[RMark]{make.design.data}}.
#'
#' @details The first time \code{markClosedSCR} is called, it will likely produce a firewall warning alerting users that R has requested the ability to accept incoming network connections. Incoming network connections are required to use parallel processing as implemented in \code{markClosed}. Note that setting \code{parms="all"} is required for any \code{markClosed} model output to be used in \code{\link{multimodelClosed}}.
#' @return A list containing the following:
#' \item{mcmc}{Markov chain Monte Carlo object of class \code{\link[coda]{mcmc.list}}.}
#' \item{mod.p}{Model formula for detection probability (as specified by \code{mod.p} above).}
#' \item{mod.delta}{Formula always \code{NULL}; only for internal use in \code{\link{multimodelClosedSCR}}.}
#' \item{mod.det}{Model formula for detection function (as specified by \code{detection} above).}
#' \item{DM}{A list of design matrices for detection probability generated for model \code{mod.p}, where DM$p is the design matrix for initial capture probability (p) and DM$c is the design matrix for recapture probability (c).}
#' \item{initial.values}{A list containing the parameter and latent variable values at iteration \code{iter} for each chain. Values are provided for "\code{pbeta}", "\code{N}", "\code{sigma2_scr}", and "\code{centers}".}
#' \item{mms}{An object of class \code{multimarkSCRsetup}}
#' @author Brett T. McClintock
#' @seealso \code{\link{multimodelClosedSCR}}
#' @references
#' Gopalaswamy, A.M., Royle, J.A., Hines, J.E., Singh, P., Jathanna, D., Kumar, N. and Karanth, K.U. 2012. Program SPACECAP: software for estimating animal density using spatially explicit capture-recapture models. \emph{Methods in Ecology and Evolution} 3:1067-1072.
#'
#' King, R., McClintock, B. T., Kidney, D., and Borchers, D. L. 2016. Capture-recapture abundance estimation using a semi-complete data likelihood approach. \emph{The Annals of Applied Statistics} 10: 264-285 
#' 
#' Royle, J.A., Karanth, K.U., Gopalaswamy, A.M. and Kumar, N.S. 2009. Bayesian inference in camera trapping studies for a class of spatial capture-recapture models.  \emph{Ecology} 90: 3233-3244.
#'
#' @examples
#' \dontshow{
#' sim.data<-simdataClosedSCR(delta_1=1,delta_2=0)
#' Enc.Mat<-sim.data$Enc.Mat
#' trapCoords<-sim.data$spatialInputs$trapCoords
#' studyArea<-sim.data$spatialInputs$studyArea
#' test<-markClosedSCR(Enc.Mat,trapCoords,studyArea,iter=10,burnin=0,bin=5)}
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # Example uses unrealistically low values for nchain, iter, and burnin
#' 
#' #Run single chain using the default model for ``traditional'' tiger data of Royle et al (2009)
#' Enc.Mat<-tiger$Enc.Mat
#' trapCoords<-tiger$trapCoords
#' studyArea<-tiger$studyArea
#' tiger.dot<-markClosedSCR(Enc.Mat,trapCoords,studyArea,iter=100,adapt=50,burnin=50)
#' 
#' #Posterior summary for monitored parameters
#' summary(tiger.dot$mcmc)
#' plot(tiger.dot$mcmc)}
#' 
#' @export
markClosedSCR<-function(Enc.Mat,trapCoords,studyArea=NULL,buffer=NULL,ncells=1024,covs=data.frame(),mod.p=~1,detection="half-normal",parms=c("pbeta","N"),nchains=1,iter=12000,adapt=1000,bin=50,thin=1,burnin=2000,taccept=0.44,tuneadjust=0.95,proppbeta=0.1,propsigma=1,propcenter=NULL,sigma_bounds=NULL,mu0=0,sigma2_mu0=1.75,initial.values=NULL,scalemax=10,printlog=FALSE,...){
  if(any(Enc.Mat>1 | Enc.Mat<0)) stop("With a single mark type, encounter histories can only contain 0's (non-detections) and 1's (detections)")
  mms <- processdataSCR(Enc.Mat,trapCoords,studyArea,buffer,ncells,covs=covs,known=rep(1,nrow(Enc.Mat)),scalemax=scalemax)
  out <- multimarkClosedSCR(mms=mms,mod.p=mod.p,mod.delta=~NULL,detection=detection,parms=parms,nchains=nchains,iter=iter,adapt=adapt,bin=bin,thin=thin,burnin=burnin,taccept=taccept,tuneadjust=tuneadjust,proppbeta=proppbeta,propsigma=propsigma,propcenter=propcenter,sigma_bounds=sigma_bounds,mu0=mu0,sigma2_mu0=sigma2_mu0,initial.values=initial.values,scalemax=scalemax,printlog=printlog,...)
  out$initial.values <- lapply(out$initial.values,function(x) list(pbeta=x$pbeta,sigma2_scr=x$sigma2_scr,N=x$N,centers=x$centers))
  return(out)
}

#' Fit spatially-explicit population abundance models for capture-mark-recapture data consisting of multiple non-invasive marks
#'
#' This function fits spatially-explicit population abundance models for capture-mark-recapture data consisting of multiple non-invasive marks using Bayesian analysis methods. Markov chain Monte Carlo (MCMC) is used to draw samples from the joint posterior distribution. 
#'
#'
#' @param Enc.Mat A matrix containing the observed encounter histories with rows corresponding to individuals and (\code{ntraps}*\code{noccas}) columns corresponding to traps and sampling occasions.  The first \code{noccas} columns correspond to trap 1, the second \code{noccas} columns corresopond to trap 2, etc. Ignored unless \code{mms=NULL}.
#' @param trapCoords A matrix of dimension \code{ntraps} x (2 + \code{noccas}) indicating the Cartesian coordinates and operating occasions for the traps, where rows correspond to trap, the first column the x-coordinate (``x''), and the second column the y-coordinate (``y''). The last \code{noccas} columns indicate whether or not the trap was operating on each of the occasions, where `1' indciates the trap was operating and `0' indicates the trap was not operating. Ignored unless \code{mms=NULL}.
#' @param studyArea is a 3-column matrix containing the coordinates for the centroids of a contiguous grid of cells that define the study area and available habitat. Each row corresponds to a grid cell. The first 2 columns (``x'' and ``y'') indicate the Cartesian x- and y-coordinate for the centroid of each grid cell, and the third column (``avail'') indicates whether the cell is available habitat (=1) or not (=0). All cells must be square and have the same resolution. If \code{studyArea=NULL} (the default) and  \code{mms=NULL}, then a square study area grid composed of \code{ncells} cells of available habitat is drawn around the bounding box of \code{trapCoords} based on \code{buffer}. Ignored unless \code{mms=NULL}.
#' Note that rows should be ordered in raster cell order (raster cell numbers start at 1 in the upper left corner, and increase from left to right, and then from top to bottom).
#' @param buffer A scaler in same units as \code{trapCoords} indicating the buffer around the bounding box of \code{trapCoords} for defining the study area when \code{studyArea=NULL}.  Ignored unless \code{studyArea=NULL} and \code{mms=NULL}.
#' @param ncells The number of grid cells in the study area when \code{studyArea=NULL}. The square root of \code{ncells} must be a whole number. Default is \code{ncells=1024}. Ignored unless \code{studyArea=NULL} and \code{mms=NULL}.
#' @param data.type Specifies the encounter history data type. All data types include non-detections (type 0 encounter), type 1 encounter (e.g., left-side), and type 2 encounters (e.g., right-side). When both type 1 and type 2 encounters occur for the same individual within a sampling occasion, these can either be "non-simultaneous" (type 3 encounter) or "simultaneous" (type 4 encounter). Three data types are currently permitted:
#' 
#'  \code{data.type="never"} indicates both type 1 and type 2 encounters are never observed for the same individual within a sampling occasion, and observed encounter histories therefore include only type 1 or type 2 encounters (e.g., only left- and right-sided photographs were collected). Observed encounter histories can consist of non-detections (0), type 1 encounters (1), and type 2 encounters (2). See \code{\link{bobcat}}. Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 3 encounters (3).
#'
#'  \code{data.type="sometimes"} indicates both type 1 and type 2 encounters are sometimes observed (e.g., both-sided photographs are sometimes obtained, but not necessarily for all individuals). Observed encounter histories can consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), type 3 encounters (3), and type 4 encounters (4). Type 3 encounters can only be observed when an individual has at least one type 4 encounter. Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), type 3 encounters (3), and type 4 encounters (4). 
#'
#'  \code{data.type="always"} indicates both type 1 and type 2 encounters are always observed, but some encounter histories may still include only type 1 or type 2 encounters. Observed encounter histories can consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 4 encounters (4). Latent encounter histories consist of non-detections (0), type 1 encounters (1), type 2 encounters (2), and type 4 encounters (4).
#'
#' @param covs A data frame of time- and/or trap-dependent covariates for detection probabilities (ignored unless \code{mms=NULL}). The number of rows in the data frame must equal the number of traps times the number of sampling occasions (\code{ntraps*noccas}), where the first \code{noccas} rows correspond to trap 1, the second \code{noccas} rows correspond to trap 2, etc. Covariate names cannot be "time", "age", or "h"; these names are reserved for temporal, behavioral, and individual effects when specifying \code{mod.p} and \code{mod.phi}.
#' @param mms An optional object of class \code{multimarkSCRsetup-class}; if \code{NULL} it is created. See \code{\link{processdataSCR}}.
#' @param mod.p Model formula for detection probability as a function of distance from activity centers. For example, \code{mod.p=~1} specifies no effects (i.e., intercept only) other than distance, \code{mod.p~time} specifies temporal effects, \code{mod.p~c} specifies behavioral reponse (i.e., trap "happy" or "shy"), \code{mod.p~trap} specifies trap effects, and \code{mod.p~time+c} specifies additive temporal and behavioral effects.
#' @param mod.delta Model formula for conditional probabilities of type 1 (delta_1) and type 2 (delta_2) encounters, given detection. Currently only \code{mod.delta=~1} (i.e., \eqn{\delta_1 = \delta_2}) and \code{mod.delta=~type} (i.e., \eqn{\delta_1 \ne \delta_2}) are implemented.
#' @param detection Model for detection probability as a function of distance from activity centers . Must be "\code{half-normal}" (of the form \eqn{\exp{(-d^2 / (2*\sigma^2))}}, where \eqn{d} is distance) or "\code{exponential}" (of the form \eqn{\exp{(-d / \lambda)}}).
#' @param parms A character vector giving the names of the parameters and latent variables to monitor. Possible parameters are cloglog-scale detection probability parameters ("\code{pbeta}"), population abundance ("\code{N}"), conditional probability of type 1 or type 2 encounter, given detection ("\code{delta})", probability of simultaneous type 1 and type 2 detection, given both types encountered ("\code{alpha}"), cloglog-scale distance term for the detection function ("\code{sigma2_scr}" when \code{detection=``half-normal''} or "\code{lambda}" when \code{detection=``exponential''}), and the probability that a randomly selected individual from the \code{M = nrow(Enc.Mat)} observed individuals belongs to the \eqn{n} unique individuals encountered at least once ("\code{psi}"). Individual activity centers ("\code{centers}"), encounter history indices ("\code{H}"), and the log posterior density ("\code{logPosterior}") may also be monitored. Setting \code{parms="all"} monitors all possible parameters and latent variables.
#' @param nchains The number of parallel MCMC chains for the model.
#' @param iter The number of MCMC iterations.
#' @param adapt The number of iterations for proposal distribution adaptation. If \code{adapt = 0} then no adaptation occurs.
#' @param bin Bin length for calculating acceptance rates during adaptive phase (\code{0 < bin <= iter}).
#' @param thin Thinning interval for monitored parameters.
#' @param burnin Number of burn-in iterations (\code{0 <= burnin < iter}).
#' @param taccept Target acceptance rate during adaptive phase (\code{0 < taccept <= 1}). Acceptance rate is monitored every \code{bin} iterations. Default is \code{taccept = 0.44}.
#' @param tuneadjust Adjustment term during adaptive phase (\code{0 < tuneadjust <= 1}). If acceptance rate is less than \code{taccept}, then proposal term (\code{proppbeta} or \code{propsigma}) is multiplied by \code{tuneadjust}. If acceptance rate is greater than or equal to \code{taccept}, then proposal term is divided by \code{tuneadjust}. Default is \code{tuneadjust = 0.95}.
#' @param proppbeta Scaler or vector (of length k) specifying the initial standard deviation of the Normal(pbeta[j], proppbeta[j]) proposal distribution. If \code{proppbeta} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{proppbeta = 0.1}.
#' @param propsigma Scaler specifying the initial Gamma(shape = 1/\code{propsigma}, scale = sigma_scr * \code{propsigma}) proposal distribution for sigma_scr = sqrt(sigma2_scr) (or sqrt(lambda) = lambda if \code{detection=``exponential''}). Default is \code{propsigma=1}.
#' @param propcenter Scaler specifying the neighborhood distance when proposing updates to activity centers. When \code{propcenter=NULL} (the default), then propcenter = a*10, where a is the cell size for the study area grid, and each cell has (at most) approximately 300 neighbors. 
#' @param maxnumbasis Maximum number of basis vectors to use when proposing latent history frequency updates. Default is \code{maxnumbasis = 1}, but higher values can potentially improve mixing.
#' @param a0delta Scaler or vector (of length d) specifying the prior for the conditional (on detection) probability of type 1 (delta_1), type 2 (delta_2), and both type 1 and type 2 encounters (1-delta_1-delta_2). If \code{a0delta} is a scaler, then this value is used for all a0delta[j] for j = 1, ..., d. For \code{mod.delta=~type}, d=3 with [delta_1, delta_2, 1-delta_1-delta_2] ~ Dirichlet(a0delta) prior. For \code{mod.delta=~1}, d=2 with [tau] ~ Beta(a0delta[1],a0delta[2]) prior, where (delta_1,delta_2,1-delta_1-delta_2) = (tau/2,tau/2,1-tau). See McClintock et al. (2013) for more details.
#' @param a0alpha Specifies "shape1" parameter for [alpha] ~ Beta(a0alpha, b0alpha) prior. Only applicable when \code{data.type = "sometimes"}. Default is \code{a0alpha = 1}. Note that when \code{a0alpha = 1} and \code{b0alpha = 1}, then [alpha] ~ Unif(0,1).
#' @param b0alpha Specifies "shape2" parameter for [alpha] ~ Beta(a0alpha, b0alpha) prior. Only applicable when \code{data.type = "sometimes"}. Default is \code{b0alpha = 1}. Note that when \code{a0alpha = 1} and \code{b0alpha = 1}, then [alpha] ~ Unif(0,1).
#' @param sigma_bounds Positive vector of length 2 for the lower and upper bounds for the [sigma_scr] ~ Uniform(sigma_bounds[1], sigma_bounds[2]) (or [sqrt(lambda)] when \code{detection=``exponential''}) prior for the detection function term sigma_scr = sqrt(sigma2_scr) (or sqrt(lambda)). When \code{sigma_bounds = NULL} (the default), then \code{sigma_bounds = c(1.e-6,max(diff(range(studyArea[,"x"])),diff(range(studyArea[,"y"]))))}.
#' @param mu0 Scaler or vector (of length k) specifying mean of pbeta[j] ~ Normal(mu0[j], sigma2_mu0[j]) prior. If \code{mu0} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{mu0 = 0}.
#' @param sigma2_mu0 Scaler or vector (of length k) specifying variance of pbeta[j] ~ Normal(mu0[j], sigma2_mu0[j]) prior. If \code{sigma2_mu0} is a scaler, then this value is used for all j = 1, ..., k. Default is \code{sigma2_mu0 = 1.75}.
#' @param a0psi Specifies "shape1" parameter for [psi] ~ Beta(a0psi,b0psi) prior. Default is \code{a0psi = 1}.
#' @param b0psi Specifies "shape2" parameter for [psi] ~ Beta(a0psi,b0psi) prior. Default is \code{b0psi = 1}.
#' @param initial.values Optional list of \code{nchain} list(s) specifying initial values for parameters and latent variables. Default is \code{initial.values = NULL}, which causes initial values to be generated automatically. In addition to the parameters ("\code{pbeta}", "\code{N}", "\code{delta_1}", "\code{delta_2}", "\code{alpha}", "\code{sigma2_scr}", "\code{centers}", and "\code{psi}"), initial values can be specified for the initial latent history frequencies ("\code{x}") and initial individual encounter history indices ("\code{H}").
#' @param known Optional integer vector indicating whether the encounter history of an individual is known with certainty (i.e., the observed encounter history is the true encounter history). Encounter histories with at least one type 4 encounter are automatically assumed to be known, and \code{known} does not need to be specified unless there exist encounter histories that do not contain a type 4 encounter that happen to be known with certainty (e.g., from independent telemetry studies). If specified, \code{known = c(v_1,v_2,...,v_M)} must be a vector of length \code{M = nrow(Enc.Mat)} where \code{v_i = 1} if the encounter history for individual \code{i} is known (\code{v_i = 0} otherwise). Note that known all-zero encounter histories (e.g., `000') are ignored.
#' @param scalemax Upper bound for internal re-scaling of grid cell centroid coordinates. Default is \code{scalemax=10}, which re-scales the centroids to be between 0 and 10.  Re-scaling is done internally to avoid numerical overflows during model fitting. Ignored unless \code{mms=NULL}.
#' @param printlog Logical indicating whether to print the progress of chains and any errors to a log file in the working directory. Ignored when \code{nchains=1}. Updates are printed to log file as 1\% increments of \code{iter} of each chain are completed. With >1 chains, setting \code{printlog=TRUE} is probably most useful for Windows users because progress and errors are automatically printed to the R console for "Unix-like" machines (i.e., Mac and Linux) when \code{printlog=FALSE}. Default is \code{printlog=FALSE}.
#' @param ... Additional "\code{parameters}" arguments for specifying \code{mod.p}. See \code{\link[RMark]{make.design.data}}.
#'
#' @details The first time \code{multimarkSCRClosed} is called, it will likely produce a firewall warning alerting users that R has requested the ability to accept incoming network connections. Incoming network connections are required to use parallel processing as implemented in \code{multimarkClosed}. Note that setting \code{parms="all"} is required for any \code{multimarkClosed} model output to be used in \code{\link{multimodelClosed}}.
#' @return A list containing the following:
#' \item{mcmc}{Markov chain Monte Carlo object of class \code{\link[coda]{mcmc.list}}.}
#' \item{mod.p}{Model formula for detection probability (as specified by \code{mod.p} above).}
#' \item{mod.delta}{Model formula for conditional probability of type 1 or type 2 encounter, given detection (as specified by \code{mod.delta} above).}
#' \item{mod.det}{Model formula for detection function (as specified by \code{detection} above).}
#' \item{DM}{A list of design matrices for detection probability generated for model \code{mod.p}, where DM$p is the design matrix for initial capture probability (p) and DM$c is the design matrix for recapture probability (c).}
#' \item{initial.values}{A list containing the parameter and latent variable values at iteration \code{iter} for each chain. Values are provided for "\code{pbeta}", "\code{N}", "\code{delta_1}", "\code{delta_2}", "\code{alpha}", "\code{sigma2_scr}", "\code{centers}", "\code{psi}", "\code{x}", and "\code{H}".}
#' \item{mms}{An object of class \code{multimarkSCRsetup}}
#' @author Brett T. McClintock
#' @seealso \code{\link{processdataSCR}}.
#' @references
#' Bonner, S. J., and Holmberg J. 2013. Mark-recapture with multiple, non-invasive marks. \emph{Biometrics} 69: 766-775.
#' 
#' Gopalaswamy, A.M., Royle, J.A., Hines, J.E., Singh, P., Jathanna, D., Kumar, N. and Karanth, K.U. 2012. Program SPACECAP: software for estimating animal density using spatially explicit capture-recapture models. \emph{Methods in Ecology and Evolution} 3:1067-1072.
#'
#' King, R., McClintock, B. T., Kidney, D., and Borchers, D. L. 2016. Capture-recapture abundance estimation using a semi-complete data likelihood approach. \emph{The Annals of Applied Statistics} 10: 264-285 
#' 
#' McClintock, B. T., Conn, P. B., Alonso, R. S., and Crooks, K. R. 2013. Integrated modeling of bilateral photo-identification data in mark-recapture analyses. \emph{Ecology} 94: 1464-1471.
#' 
#' McClintock, B. T., Bailey, L. L., Dreher, B. P., and Link, W. A. 2014. Probit models for capture-recapture data subject to imperfect detection, individual heterogeneity and misidentification. \emph{The Annals of Applied Statistics} 8: 2461-2484.
#' 
#' Royle, J.A., Karanth, K.U., Gopalaswamy, A.M. and Kumar, N.S. 2009. Bayesian inference in camera trapping studies for a class of spatial capture-recapture models.  \emph{Ecology} 90: 3233-3244.
#'
#' @examples
#' \dontshow{
#' sim.data<-simdataClosedSCR(N=30,noccas=5,ntraps=4)
#' Enc.Mat <- sim.data$Enc.Mat
#' trapCoords <- sim.data$spatialInputs$trapCoords
#' studyArea <- sim.data$spatialInputs$studyArea
#' test<-multimarkClosedSCR(Enc.Mat,trapCoords,studyArea,iter=10,burnin=0,bin=5)}
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # Example uses unrealistically low values for nchain, iter, and burnin
#' 
#' #Generate object of class "multimarkSCRsetup" from simulated data
#' sim.data<-simdataClosedSCR()
#' Enc.Mat <- sim.data$Enc.Mat
#' trapCoords <- sim.data$spatialInputs$trapCoords
#' studyArea <- sim.data$spatialInputs$studyArea
#' 
#' #Run single chain using the default model for simulated data
#' example.dot<-multimarkClosedSCR(Enc.Mat,trapCoords,studyArea)
#' 
#' #Posterior summary for monitored parameters
#' summary(example.dot$mcmc)
#' plot(example.dot$mcmc)}
#' 
#' @export
#' @importFrom methods validObject
multimarkClosedSCR<-function(Enc.Mat,trapCoords,studyArea=NULL,buffer=NULL,ncells=1024,data.type="never",covs=data.frame(),mms=NULL,mod.p=~1,mod.delta=~type,detection="half-normal",parms=c("pbeta","delta","N"),nchains=1,iter=12000,adapt=1000,bin=50,thin=1,burnin=2000,taccept=0.44,tuneadjust=0.95,proppbeta=0.1,propsigma=1,propcenter=NULL,maxnumbasis=1,a0delta=1,a0alpha=1,b0alpha=1,sigma_bounds=NULL,mu0=0,sigma2_mu0=1.75,a0psi=1,b0psi=1,initial.values=NULL,known=integer(),scalemax=10,printlog=FALSE,...){
  
  if(is.null(mms)) mms <- processdataSCR(Enc.Mat,trapCoords,studyArea,buffer,ncells,data.type,covs,known,scalemax)
  if(!inherits(mms,"multimarkSCRsetup")) stop("'mms' must be an object of class 'multimarkSCRsetup'")
  validObject(mms)
  
  match.arg(detection,c("half-normal","exponential"))
  if(is.null(detection)) stop("detection function cannot be NULL")
  
  if(!inherits(mod.p,"formula")) stop("'mod.p' must be an object of class 'formula'")
  if(!inherits(mod.delta,"formula")) stop("'mod.delta' must be an object of class 'formula'")
  DM<-get_DMClosed(mod.p,mod.delta,mms@Enc.Mat,covs=mms@covs,ntraps=nrow(mms@spatialInputs$trapCoords),detection=detection,...)
  
  if(iter>0){
    if(iter<=burnin) stop(paste("'burnin' must be less than ",iter))
  } else {
    burnin<-0
  }
  
  detParm<-ifelse(DM$mod.det=="half-normal","sigma2_scr","lambda")
  if(mod.delta != ~NULL) {
    parmlist<-c("pbeta","delta","N",detParm,"alpha","psi","H","centers","logPosterior")
  } else {
    parmlist<-c("pbeta","N",detParm,"centers","logPosterior")    
  }
  
  if(is.null(sigma_bounds)){
    sigma_bounds<-c(1.e-6,max(diff(range(mms@spatialInputs$studyArea[,"x"])),diff(range(mms@spatialInputs$studyArea[,"y"])))*mms@spatialInputs$Srange)
  }
  params <- checkClosedSCR(parms,parmlist,mms,DM,iter,adapt,bin,thin,burnin,taccept,tuneadjust,maxnumbasis,a0delta,a0alpha,b0alpha,sigma_bounds,sigma2_mu0,a0psi,b0psi)
  
  data.type<-mms@data.type
  Enc.Mat<-mms@Enc.Mat
  M<-nrow(Enc.Mat)
  covs<-mms@covs
  pdim<-ncol(DM$p)
  
  mu0 <- checkvecs(mu0,pdim,"mu0")
  sigma2_mu0 <- checkvecs(sigma2_mu0,pdim,"sigma2_mu0")
  a0delta <- checkvecs(a0delta,ifelse(mod.delta==formula(~type),3,2),"a0delta")
  
  spatialInputs <- getSpatialInputs(mms)
  
  inits<-get_initsSCR(mms,nchains,initial.values,M,data.type,a0alpha,b0alpha,a0delta,sigma_bounds/mms@spatialInputs$Srange,a0psi,b0psi,DM,spatialInputs)
  
  priorparms <-list(a0delta=a0delta,a0alpha=a0alpha,b0alpha=b0alpha,sigma_bounds=sigma_bounds,mu0=mu0,sigma2_mu0=sigma2_mu0,a0psi=a0psi,b0psi=b0psi)
  
  message("\nFitting spatial abundance model with cloglog link\n")
  if(mod.delta != ~NULL) message("data type = \"",data.type,"\"\n")
  message("detection function = \"",DM$mod.det,"\"\n")
  message("p model = ",as.character(mod.p))
  if(mod.delta != ~NULL) message("delta model = ",as.character(mod.delta))
  message("\nInitializing model \n")
  posteriorClosedSCR(inits,DM,mms,priorparms,spatialInputs)
  
  proppbeta <- checkvecs(proppbeta,pdim,"proppbeta")
  if(length(propsigma)!=1) stop("'propsigma' must be a scaler")
  
  if(!is.null(propcenter)) propcenter <- propcenter/mms@spatialInputs$Srange
  Prop.center<-getPropCenter(spatialInputs,propcenter)
  Prop.sd <- c(proppbeta,propsigma)
  
  message("Updating...",ifelse(printlog | nchains==1,"","set 'printlog=TRUE' to follow progress of chains in a working directory log file"),"\n",sep="")
  if(printlog & nchains==1) printlog<-FALSE
  
  if(nchains>1){
    if(nchains>detectCores()) warning("Number of parallel chains (nchains) is greater than number of cores \n")
    modlog <- ifelse(mod.delta != ~NULL,"multimarkClosedSCR","markClosedSCR")
    cl <- makeCluster( nchains ,outfile=ifelse(printlog,paste0(modlog,"_log_",format(Sys.time(), "%Y-%b-%d_%H%M.%S"),".txt"),""))
    clusterExport(cl,list("mcmcClosedSCR"),envir=environment())  
    clusterSetRNGStream(cl)
    chains <- parLapply(cl,1:nchains, function(ichain) 
      mcmcClosedSCR(ichain,mms,DM,params,inits,iter,adapt,bin,thin,burnin,taccept,tuneadjust,Prop.sd,Prop.center,spatialInputs,maxnumbasis,a0delta,a0alpha,b0alpha,sigma_bounds,mu0,sigma2_mu0,a0psi,b0psi,printlog))
    stopCluster(cl)
    gc()
  } else {
    chains <- vector('list',nchains)
    chains[[nchains]] <- mcmcClosedSCR(nchains,mms,DM,params,inits,iter,adapt,bin,thin,burnin,taccept,tuneadjust,Prop.sd,Prop.center,spatialInputs,maxnumbasis,a0delta,a0alpha,b0alpha,sigma_bounds,mu0,sigma2_mu0,a0psi,b0psi,printlog)
    gc()
  }
  
  chains <- processClosedSCRchains(chains,params,DM,M,nchains,iter,burnin,thin)
  return(list(mcmc=chains$chains,mod.p=mod.p,mod.delta=mod.delta,mod.det=detection,DM=list(p=DM$p,c=DM$c),initial.values=chains$initial.values,priorparms=priorparms,mms=mms))
}

#' Calculate posterior capture and recapture probabilities
#'
#' This function calculates posterior spatial capture (\eqn{p}) and recapture (\eqn{c}) probabilities (at zero distance from an activity center) for each sampling occasion from \code{\link{multimarkClosedSCR}} output. 
#'
#'
#' @param out List of output returned by \code{\link{multimarkClosedSCR}}.
#' @param link Link function for detection probability. Must be "\code{cloglog}". Note that \code{\link{multimarkClosedSCR}} is currently implemented for the cloglog link only.
#' @return An object of class \code{\link[coda]{mcmc.list}} containing the following:
#' \item{p}{Posterior samples for capture probability (\eqn{p}) for each sampling occasion (first index) and trap (second index).}
#' \item{c}{Posterior samples for recapture probability (\eqn{c}) for each sampling occasion (first index) and trap (second index).}
#' @author Brett T. McClintock
#' @seealso \code{\link{multimarkClosedSCR}}
#' @examples
#' \dontshow{
#' sim.data<-simdataClosedSCR()
#' Enc.Mat<-sim.data$Enc.Mat
#' trapCoords<-sim.data$spatialInputs$trapCoords
#' studyArea<-sim.data$spatialInputs$studyArea
#' test<-getprobsClosedSCR(multimarkClosedSCR(Enc.Mat,trapCoords,studyArea,iter=10,burnin=0,bin=5))}
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # Example uses unrealistically low values for nchain, iter, and burnin
#' 
#' #Run behavior model for simulated data with constant detection probability (i.e., mod.p=~c)
#' sim.data<-simdataClosedSCR()
#' Enc.Mat<-sim.data$Enc.Mat
#' trapCoords<-sim.data$spatialInputs$trapCoords
#' studyArea<-sim.data$spatialInputs$studyArea
#' example.c <- multimarkClosedSCR(Enc.Mat,trapCoords,studyArea,mod.p=~c,
#'                                 iter=1000,adapt=500,burnin=500)
#'   
#' #Calculate capture and recapture probabilities
#' pc <- getprobsClosedSCR(example.c)
#' summary(pc)}
#' 
#' @export
getprobsClosedSCR<-function(out,link="cloglog"){
  
  DMp<-out$DM$p
  DMc<-out$DM$c
  
  ntraps<-nrow(out$mms@spatialInputs$trapCoords)
  noccas<-ncol(out$mms@Enc.Mat)/ntraps
  
  if(noccas*ntraps<2) stop("must have >1 sampling occasion or >1 trap")
  
  pbetanames<-paste0("pbeta[",colnames(DMp),"]")
  nchains<-length(out$mcmc)
  
  pc<-vector("list",nchains)
  
  varind <- is.null(varnames(out$mcmc))
  if(!varind){
    vars <- varnames(out$mcmc)
  } else {
    vars <- names(out$mcmc[[1]])    
  }
  if(!any(match(pbetanames,vars,nomatch=0))) stop("'pbeta' parameters not found")
  
  for(ichain in 1:nchains){
    
    if(!varind){
      pbeta<-out$mcmc[[ichain]][,pbetanames,drop=FALSE]   
    } else {
      pbeta<-matrix(out$mcmc[[ichain]][pbetanames],nrow=1) 
    }
    
    if(link=="cloglog"){
      p <- matrix(invcloglog(apply(pbeta,1,function(x) DMp%*%x)),byrow=T,ncol=noccas*ntraps)
      rc <- matrix(invcloglog(apply(pbeta,1,function(x) DMc%*%x)),byrow=T,ncol=noccas*ntraps)
    } else {
      stop("link function must be 'cloglog'")
    }
    
    if(dim(rc)[1]==1){
      rc <- matrix(rc[,-seq(1,ntraps*noccas,noccas)],nrow=1)
    } else if(ntraps*noccas<3){
      rc <- matrix(rc[,-1],ncol=1)      
    } else {
      rc <- rc[,-seq(1,ntraps*noccas,noccas)]
    }
    colnames(p)  <- paste0("p[",1:noccas,",",rep(1:ntraps,each=noccas),"]")
    colnames(rc)  <- paste0("c[",2:noccas,",",rep(1:ntraps,each=noccas-1),"]")
    pc[[ichain]]<- mcmc(cbind(p,rc),start=start(out$mcmc),end=end(out$mcmc),thin=attributes(out$mcmc[[ichain]])$mcpar[3])
  }
  return(as.mcmc.list(pc))
}

#' Calculate population density estimates
#'
#' This function calculates posterior population density estimates from \code{\link{multimarkClosedSCR}} output as D = N/A, where D is density, N is abundance, and A is the area of available habitat within the study area. 
#'
#'
#' @param out List of output returned by \code{\link{multimarkClosedSCR}}.
#' @return An object of class \code{\link[coda]{mcmc.list}} containing the following:
#' \item{D}{Posterior samples for density.}
#' @author Brett T. McClintock
#' @seealso \code{\link{multimarkClosedSCR}}
#' @examples
#' \dontshow{
#' sim.data<-simdataClosedSCR()
#' Enc.Mat<-sim.data$Enc.Mat
#' trapCoords<-sim.data$spatialInputs$trapCoords
#' studyArea<-sim.data$spatialInputs$studyArea
#' test<-getdensityClosedSCR(multimarkClosedSCR(Enc.Mat,trapCoords,studyArea,iter=10,burnin=0,bin=5))}
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # Example uses unrealistically low values for nchain, iter, and burnin
#' 
#' #Run behavior model for simulated data with constant detection probability (i.e., mod.p=~c)
#' sim.data<-simdataClosedSCR()
#' Enc.Mat<-sim.data$Enc.Mat
#' trapCoords<-sim.data$spatialInputs$trapCoords
#' studyArea<-sim.data$spatialInputs$studyArea
#' example.dot <- multimarkClosedSCR(Enc.Mat,trapCoords,studyArea,mod.p=~1)
#'   
#' #Calculate capture and recapture probabilities
#' D <- getdensityClosedSCR(example.dot)
#' summary(D)}
#' 
#' @export
getdensityClosedSCR<-function(out){
  
  nchains<-length(out$mcmc)
  
  spatialInputs<-getSpatialInputs(out$mms)
  
  A <- out$mms@spatialInputs$origCellRes^2 * nrow(spatialInputs$studyArea)
  
  varind <- is.null(varnames(out$mcmc))
  if(!varind){
    vars <- varnames(out$mcmc)
  } else {
    vars <- names(out$mcmc[[1]])    
  }
  if(!any(match("N",vars,nomatch=0))) stop("'N' parameter not found")
  
  D <- vector("list",nchains)
  
  for(ichain in 1:nchains){
    iD <- out$mcmc[[ichain]][,"N"]/A
    D[[ichain]]<- mcmc(iD,start=start(out$mcmc),end=end(out$mcmc),thin=attributes(out$mcmc[[ichain]])$mcpar[3])
  }
  return(as.mcmc.list(D))
}

monitorparmsClosedSCR <- function(parms,parmlist,noccas,ntraps){ 
  
  if(!all(match(parms,parmlist,nomatch=0))) stop(paste0("monitored parameters ('monparms') can only include: ",paste(parmlist[-length(parmlist)],collapse=", "),", or ",parmlist[length(parmlist)]))
  
  commonparms <- parms
  
  if(any(parms=="D")){
    namesD <- "D"
  } else {
    namesD <- NULL
  }
  getcloglogp <- derivedcloglogfun(parms,"p")
  getcloglogc <- derivedcloglogfun(parms,"c")  
  
  if(any(parms=="p")){
    namesp <- paste0("p[",1:noccas,",",rep(1:ntraps,each=noccas),"]")
    commonparms <- commonparms[-which(parms=="p")]
    parms <- parms[-which(parms=="p")]
    parms <- c(parms,namesp)
  } else {
    namesp <- NULL
  }
  if(any(parms=="c")){
    namesc <- paste0("c[",2:noccas,",",rep(1:ntraps,each=noccas-1),"]") 
    commonparms <- commonparms[-which(parms=="c")]
    parms <- parms[-which(parms=="c")]
    parms <- c(parms,namesc)
  } else {
    namesc <- NULL
  }
  list(commonparms=commonparms,parms=parms,namesD=namesD,namesp=namesp,namesc=namesc,getcloglogp=getcloglogp,getcloglogc=getcloglogc)
}

getcurClosedSCRparmslist<-function(cur.parms,DM,M,noccas,data_type,alpha,centermap1,centermap2,Srange2){
  
  parmslist=vector('list',1)
  parmslist[[1]]$H<-cur.parms[paste0("H[",1:M,"]")]
  parmslist[[1]]$N <- cur.parms["N"]
  parmslist[[1]]$pbeta <- cur.parms[paste0("pbeta[",colnames(DM$p),"]")]
  parmslist[[1]]$centers <- mapCenters(cur.parms[paste0("center[",1:M,"]")],centermap1,centermap2)
  names(parmslist[[1]]$centers) <- paste0("center[",1:M,"]")
  parmslist[[1]]$sigma2_scr <- cur.parms["sigma2_scr"] / Srange2
  parmslist[[1]]$lambda <- cur.parms["lambda"] / Srange2
  
  parmslist[[1]]$psi <- cur.parms["psi"]
  parmslist[[1]]$delta_1 <- cur.parms["delta_1"]
  parmslist[[1]]$delta_2 <- cur.parms["delta_2"]
  parmslist[[1]]$delta <- cur.parms["delta"]
  
  if(data_type=="sometimes"){
    parmslist[[1]]$alpha <- cur.parms["alpha"]
  } else {
    parmslist[[1]]$alpha <- alpha   
  }
  parmslist
}

#' @importFrom utils flush.console
#' @importFrom Brobdingnag brob as.brob sum
rjmcmcClosedSCR <- function(ichain,mms,M,noccas,ntraps,spatialInputs,data_type,alpha,C,All.hists,modlist,DMlist,deltalist,detlist,priorlist,mod.p.h,iter,miter,mburnin,mthin,modprior,M1,monitorparms,missing,pbetapropsd,sigpropmean,sigpropsd,pmodnames,deltamodnames,printlog){
  
  multimodel <- matrix(0,nrow=(max(1,floor(miter/mthin)))-(floor(mburnin/mthin)),ncol=length(monitorparms$parms)+1,dimnames=list(NULL,c(monitorparms$parms,"M")))
  
  nmod <- length(modlist)
  mod.prob.brob <- as.brob(numeric(nmod))
  
  commonparms <- monitorparms$commonparms

  if(any(unlist(lapply(deltalist,function(x) {x== ~NULL })))){
    H<-get_H(mms,mms@naivex)
    names(H)<-paste0("H[",1:M,"]")
  } else {
    H<-NULL
  }
  
  A <- mms@spatialInputs$origCellRes^2 * nrow(spatialInputs$studyArea)
  Srange2 <- mms@spatialInputs$Srange^2
  for(imod in 1:nmod){
    priorlist[[imod]]$sigma_bounds<-priorlist[[imod]]$sigma_bounds/mms@spatialInputs$Srange
  }
  
  sigppropshape <- sigpropmean^2/(sigpropsd^2) + 2
  sigppropscale <- sigpropmean*(sigppropshape-1)/mms@spatialInputs$Srange
  
  M.cur<- M1
  
  modmissingparms <- drawmissingClosed(M.cur,missing,pbetapropsd,sigppropshape,sigppropscale)
  cur.parms <- c(modlist[[M.cur]][sample(iter,1),],modmissingparms,H)
  
  DM <- DMlist[[M.cur]]
  DM$mod.delta <- deltalist[[M.cur]]
  DM$mod.det <- detlist[[M.cur]]
  DM$mod.p.h <- mod.p.h[[M.cur]]
  
  cur.parms.list <- getcurClosedSCRparmslist(cur.parms,DM,M,noccas,data_type,alpha,spatialInputs$centermap1,spatialInputs$centermap2,Srange2)  
  
  for(iiter in 1:miter){
    
    mod.prob.brob[M.cur] <- getbrobprobClosed(M.cur,modprior,cur.parms["logPosterior"],cur.parms,missing,pbetapropsd,sigppropshape,sigppropscale)
    
    for(imod in (1:nmod)[-M.cur]){ 
      
      DM <- DMlist[[imod]]
      DM$mod.delta <- deltalist[[imod]]
      DM$mod.det <- detlist[[imod]]
      DM$mod.p.h <- mod.p.h[[imod]]
      
      cur.parms.list[[1]]$pbeta <- cur.parms[paste0("pbeta[",colnames(DM$p),"]")]
      
      loglike <- loglikeClosedSCR(cur.parms.list[[1]],DM,noccas,ntraps,C,All.hists,spatialInputs)
      
      posterior <- loglike + priorsClosedSCR(cur.parms.list[[1]],DM,priorlist[[imod]],data_type,spatialInputs)
      
      mod.prob.brob[imod] <- getbrobprobClosed(imod,modprior,posterior,cur.parms,missing,pbetapropsd,sigppropshape,sigppropscale)
    }
    
    if(any(is.na(as.numeric(mod.prob.brob)))){
      warning(paste0("'NA' posterior for model '","p(",pmodnames[is.na(as.numeric(mod.prob.brob))],")delta(",deltamodnames[is.na(as.numeric(mod.prob.brob))],")' at iteration ",iiter,"; model move rejected."))
      flush.console()
    } else {       
      mod.prob <- as.numeric(mod.prob.brob/Brobdingnag::sum(mod.prob.brob))
      M.cur <- (1:nmod)[rmultinom(1, 1, mod.prob)==1]
    }
    
    modmissingparms <- drawmissingClosed(M.cur,missing,pbetapropsd,sigppropshape,sigppropscale)
    cur.parms <- c(modlist[[M.cur]][sample(iter,1),],modmissingparms,H)
    
    DM <- DMlist[[M.cur]]
    DM$mod.delta <- deltalist[[M.cur]]
    DM$mod.det <- detlist[[M.cur]]
    DM$mod.p.h <- mod.p.h[[M.cur]]
    
    cur.parms.list <- getcurClosedSCRparmslist(cur.parms,DM,M,noccas,data_type,alpha,spatialInputs$centermap1,spatialInputs$centermap2,Srange2)
    
    if(iiter>mburnin & !iiter%%mthin){
      multimodel[iiter/mthin-floor(mburnin/mthin),"M"] <- M.cur
      multimodel[iiter/mthin-floor(mburnin/mthin),commonparms] <- cur.parms[commonparms]
      multimodel[iiter/mthin-floor(mburnin/mthin),monitorparms$namesD] <- cur.parms["N"] / A
      multimodel[iiter/mthin-floor(mburnin/mthin),monitorparms$namesp] <- monitorparms$getcloglogp(DM$p,cur.parms.list[[1]]$pbeta)
      multimodel[iiter/mthin-floor(mburnin/mthin),monitorparms$namesc] <- monitorparms$getcloglogc(DM$c,cur.parms.list[[1]]$pbeta)[-seq(1,ntraps*noccas,noccas)]
    }
    
    if(!(iiter%%(miter/ min(miter,100)))) {
      if(printlog){
        cat("Chain ",ichain," is ",100*(iiter/miter),"% complete \n",sep="")        
      } else{
        cat("\rChain ",ichain," is ",100*(iiter/miter),"% complete",sep="")
      }
    }
  }
  return(multimodel)
}

#' Multimodel inference for 'multimark' spatial population abundance models
#'
#' This function performs Bayesian multimodel inference for a set of 'multimark' spatial population abundance models using the reversible jump Markov chain Monte Carlo (RJMCMC) algorithm proposed by Barker & Link (2013).
#'
#'
#' @param modlist A list of individual model output lists returned by \code{\link{multimarkClosedSCR}} or \code{\link{markClosedSCR}}. The models must have the same number of chains and MCMC iterations.
#' @param modprior Vector of length \code{length(modlist)} containing prior model probabilities. Default is \code{modprior = rep(1/length(modlist), length(modlist))}.
#' @param monparms Parameters to monitor. Only parameters common to all models can be monitored (e.g., "\code{pbeta[(Intercept)]}", "\code{N}", "\code{sigma2_scr}"), but derived density ("\code{D}") as well as capture ("\code{p}") and recapture ("\code{c}") probabilities (at distance zero from activity centers) can also be monitored. Default is \code{monparms = "N"}.
#' @param miter The number of RJMCMC iterations per chain. If \code{NULL}, then the number of MCMC iterations for each individual model chain is used.
#' @param mburnin Number of burn-in iterations (\code{0 <= mburnin < miter}).
#' @param mthin Thinning interval for monitored parameters.
#' @param M1 Integer vector indicating the initial model for each chain, where \code{M1_j=i} initializes the RJMCMC algorithm for chain j in the model corresponding to \code{modlist[[i]]} for i=1,...,  \code{length(modlist)}. If \code{NULL}, the algorithm for all chains is initialized in the most general model. Default is \code{M1=NULL}.
#' @param pbetapropsd Scaler specifying the standard deviation of the Normal(0, pbetapropsd) proposal distribution for "\code{pbeta}"  parameters. Default is \code{pbetapropsd=1}. See Barker & Link (2013) for more details.
#' @param sigpropmean Scaler specifying the mean of the inverse Gamma proposal distribution for \code{sigma2_scr} (or \code{lambda} if \code{detection=``exponential''}). Only applies if models do not have the same detection function (i.e., ``half-normal'' or ``exponential''). Default is \code{sigpropmean=0.8}. See Barker & Link (2013) for more details.
#' @param sigpropsd Scaler specifying the standard deviation of the inverse Gamma proposal distribution for \code{sigma2_scr} (or \code{lambda} if \code{detection=``exponential''}). Only applies if models do not have the same detection function (i.e., ``half-normal'' or ``exponential''). Default is \code{sigpropsd=0.4}. See Barker & Link (2013) for more details.
#' @param printlog Logical indicating whether to print the progress of chains and any errors to a log file in the working directory. Ignored when \code{nchains=1}. Updates are printed to log file as 1\% increments of \code{iter} of each chain are completed. With >1 chains, setting \code{printlog=TRUE} is probably most useful for Windows users because progress and errors are automatically printed to the R console for "Unix-like" machines (i.e., Mac and Linux) when \code{printlog=FALSE}. Default is \code{printlog=FALSE}.
#' @details Note that setting \code{parms="all"} is required when fitting individual \code{\link{multimarkClosedSCR}} or \code{\link{markClosedSCR}} models to be included in \code{modlist}.
#' @return A list containing the following:
#' \item{rjmcmc}{Reversible jump Markov chain Monte Carlo object of class \code{\link[coda]{mcmc.list}}. Includes RJMCMC output for monitored parameters and the current model at each iteration ("\code{M}").}
#' \item{pos.prob}{A list of calculated posterior model probabilities for each chain, including the overall posterior model probabilities across all chains.}
#' @author Brett T. McClintock
#' @seealso \code{\link{multimarkClosedSCR}}, \code{\link{processdataSCR}}
#' @references
#' Barker, R. J. and Link. W. A. 2013. Bayesian multimodel inference by RJMCMC: a Gibbs sampling approach. The American Statistician 67: 150-156.
#' @examples
#' \dontshow{
#' sim.data<-simdataClosedSCR()
#' Enc.Mat<-sim.data$Enc.Mat
#' trapCoords<-sim.data$spatialInputs$trapCoords
#' studyArea<-sim.data$spatialInputs$studyArea
#' setup<-processdataSCR(Enc.Mat,trapCoords,studyArea)
#' test.dot<-multimarkClosedSCR(mms=setup,parms="all",iter=10,burnin=0,bin=5)
#' test<-multimodelClosedSCR(modlist=list(mod1=test.dot,mod2=test.dot))
#' }
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # Example uses unrealistically low values for nchain, iter, and burnin
#' 
#' #Generate object of class "multimarkSCRsetup"
#' sim.data<-simdataClosedSCR()
#' Enc.Mat<-sim.data$Enc.Mat
#' trapCoords<-sim.data$spatialInputs$trapCoords
#' studyArea<-sim.data$spatialInputs$studyArea
#' setup<-processdataSCR(Enc.Mat,trapCoords,studyArea)
#'  
#' #Run single chain using the default model for simulated data. Note parms="all".
#' example.dot <- multimarkClosedSCR(mms=setup,parms="all",iter=1000,adapt=500,burnin=500)
#' 
#' #Run single chain for simulated data with behavior effects. Note parms="all".
#' example.c <- multimarkClosedSCR(mms=setup,mod.p=~c,parms="all",iter=1000,adapt=500,burnin=500)
#' 
#' #Perform RJMCMC using defaults
#' modlist <- list(mod1=example.dot,mod2=example.c)
#' example.M <- multimodelClosedSCR(modlist=modlist,monparms=c("N","D","sigma2_scr"))
#' 
#' #Posterior model probabilities
#' example.M$pos.prob
#'  
#' #multimodel posterior summary for abundance and density
#' summary(example.M$rjmcmc[,c("N","D")])}
#' 
#' @export
multimodelClosedSCR<-function(modlist,modprior=rep(1/length(modlist),length(modlist)),monparms="N",miter=NULL,mburnin=0,mthin=1,M1=NULL,pbetapropsd=1,sigpropmean=0.8,sigpropsd=0.4,printlog=FALSE){
  
  nmod <- length(modlist)
  iter <- unlist(unique(lapply(modlist,function(x) unique(lapply(x$mcmc,nrow)))))
  nchains <- unlist(unique(lapply(modlist,function(x) length(x$mcmc))))
  mmslist <- unlist(unique(lapply(modlist, function(x) {x$mms@covs<-data.frame();x$mms})))
  
  params <- lapply(modlist,function(x) varnames(x$mcmc))
  
  if(is.null(M1)) M1 <- rep(which.max(lapply(params,length))[1],nchains)
  
  if(is.null(miter)) miter <- iter
  
  mms<-checkmmClosedinput(mmslist,modlist,nmod,nchains,iter,miter,mburnin,mthin,modprior,M1,type="SCR")
  
  spatialInputs<-getSpatialInputs(mms)
  
  noccas<-ncol(mms@spatialInputs$trapCoords[,-c(1,2),drop=FALSE])
  ntraps<-nrow(mms@spatialInputs$trapCoords)
  M<-nrow(mms@Enc.Mat)
  All.hists<-matrix(mms@vAll.hists,byrow=TRUE,ncol=noccas*ntraps)
  C<-mms@C
  
  checkparmsClosed(mms,modlist,params,parmlist=c("pbeta[(Intercept)]","N","logPosterior"),M,type="SCR")
  
  pmodnames <- lapply(modlist,function(x) x$mod.p) 
  deltamodnames <- lapply(modlist,function(x) x$mod.delta)
  detlist <- lapply(modlist,function(x) x$mod.det)
  detmodnames<-unlist(detlist)
  
  message("\nPerforming spatial population abundance Bayesian multimodel inference by RJMCMC \n")
  if(all(unlist(lapply(deltamodnames,function(x) {x!= ~NULL })))) {
    if(length(unique(detmodnames))>1){
      message(paste0("mod",1:nmod,": ","p(",pmodnames,")delta(",deltamodnames,") ",detmodnames,"\n"))
    } else {
      message(paste0("mod",1:nmod,": ","p(",pmodnames,")delta(",deltamodnames,") \n"))      
    }
  } else if(all(unlist(lapply(deltamodnames,function(x) {x== ~NULL})))){
    message(paste0("mod",1:nmod,": ","p(",pmodnames,")  (",detmodnames,")\n"))
  }
  
  missing <- missingparmnamesClosed(params,M,noccas,NULL) 
  
  monitorparms <- monitorparmsClosedSCR(monparms,c(missing$commonparms,"p","c","D"),noccas,ntraps)
  
  DMlist <- lapply(modlist,function(x) x$DM)
  deltalist <- lapply(modlist,function(x) x$mod.delta)
  priorlist <- lapply(modlist,function(x) x$priorparms) 
  mod.p.h <- unlist(lapply(modlist,function(x) any("h"==attributes(terms(x$mod.p))$term.labels)))
  
  data_type <- mms@data.type
  if(data_type=="never"){
    alpha <- 0
  } else if(data_type=="always"){
    alpha <- 1
  } else {
    alpha <- numeric(0)
  }
  
  
  message("Updating...",ifelse(printlog | nchains==1,"","set 'printlog=TRUE' to follow progress of chains in a working directory log file"),"\n",sep="")
  if(printlog & nchains==1) printlog<-FALSE
  
  if(nchains>1){
    if(nchains>detectCores()) warning("Number of parallel chains (nchains) is greater than number of cores \n")
    cl <- makeCluster( nchains ,outfile=ifelse(printlog,paste0("multimodelClosedSCR_log_",format(Sys.time(), "%Y-%b-%d_%H%M.%S"),".txt"),""))
    clusterExport(cl,list("rjmcmcClosedSCR"),envir=environment())
    clusterSetRNGStream(cl)
    multimodel <- parLapply(cl,1:nchains, function(ichain) 
        rjmcmcClosedSCR(ichain,mms,M,noccas,ntraps,spatialInputs,data_type,alpha,C,All.hists,lapply(modlist,function(x) x$mcmc[[ichain]]),DMlist,deltalist,detlist,priorlist,mod.p.h,iter,miter,mburnin,mthin,modprior,M1[ichain],monitorparms,missing,pbetapropsd,sigpropmean,sigpropsd,pmodnames,deltamodnames,printlog))
    stopCluster(cl)
    gc()
  } else {
    multimodel <- vector('list',nchains)
    multimodel[[nchains]] <- rjmcmcClosedSCR(nchains,mms,M,noccas,ntraps,spatialInputs,data_type,alpha,C,All.hists,lapply(modlist,function(x) x$mcmc[[nchains]]),DMlist,deltalist,detlist,priorlist,mod.p.h,iter,miter,mburnin,mthin,modprior,M1,monitorparms,missing,pbetapropsd,sigpropmean,sigpropsd,pmodnames,deltamodnames,printlog)
    gc()
  }
  
  if(mburnin<mthin){
    temp=seq(mthin,max(1,miter),mthin)
  } else {
    temp=seq(mthin*(floor(mburnin/mthin)+1),miter,mthin)
  }
  
  pos.prob <- vector('list',nchains)
  for(ichain in 1:nchains){
    pos.prob[[ichain]] <-hist(multimodel[[ichain]][,"M"],plot=F,breaks=0:nmod)$density
    if(all(unlist(lapply(deltamodnames,function(x) {x!= ~NULL })))){
      if(length(unique(detmodnames))>1){
        names(pos.prob[[ichain]]) <- paste0("mod",1:nmod,": ","p(",pmodnames,")delta(",deltamodnames,") ",detmodnames) 
      } else {
        names(pos.prob[[ichain]]) <- paste0("mod",1:nmod,": ","p(",pmodnames,")delta(",deltamodnames,")")       
      }
    } else {
      names(pos.prob[[ichain]]) <- paste0("mod",1:nmod,": ","p(",pmodnames,")  (",detmodnames,")")
    }
    multimodel[[ichain]] <- mcmc(multimodel[[ichain]])
    attributes(multimodel[[ichain]])$mcpar <- c(head(temp,n=1),tail(temp,n=1),mthin)
  }  
  
  multimodel <- as.mcmc.list(multimodel)
  names(pos.prob) <- paste0("chain",1:nchains)
  pos.prob[["overall"]]<- hist(unlist(multimodel[, "M"]),plot = F, breaks = 0:nmod)$density
  if(all(unlist(lapply(deltamodnames,function(x) {x!= ~NULL })))){
    if(length(unique(detmodnames))>1){
      names(pos.prob$overall) <- paste0("mod",1:nmod,": ","p(",pmodnames,")delta(",deltamodnames,") ",detmodnames)
    } else {
      names(pos.prob$overall) <- paste0("mod",1:nmod,": ","p(",pmodnames,")delta(",deltamodnames,")") 
    }
  } else {
    names(pos.prob$overall) <- paste0("mod",1:nmod,": ","p(",pmodnames,")  (",detmodnames,")")
  }
  list(rjmcmc=multimodel,pos.prob=pos.prob) 
}