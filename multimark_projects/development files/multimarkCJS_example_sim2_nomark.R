library(multimark)
library(RMark)

nchains <- 10
niter <- 1000000
nburnin <- 500000
niterCJS <- 50000
nburninCJS <- 10000

useInitial <- TRUE

noccas=10

N=250
pbeta = -0.25
phibeta= 1.2
sigma2_zphi <- 0
delta_combs <- matrix(c(0.6,0.2,0.2,
                        1/3,1/3,1/3,
                        .1,.1,.8,
                        .5,.5,0,
                        .2,.4,.4,
                        0,  0,  1,
                        1/3,0,2/3,
                        2/3,0,1/3,
                        1,  0,  0),ncol=3,byrow=TRUE)
truealpha=seq(0,1,length=5)

nsims=40

#out=outCJS<-NULL

#for(delta_1 in truedelta_1){
#out[[delta_1]] <- outCJS[[delta_1]] <- list()
#  for(delta_2 in truedelta_2){
#out[[delta_1]][[delta_2]] <- outCJS[[delta_1]][[delta_2]] <- list()
for(alpha in truealpha){
  
  for(delt in 1:nrow(delta_combs)){
    delta_1 <- delta_combs[delt,1]
    delta_2 <- delta_combs[delt,2]
    
    #out[[delta_1]][[delta_2]][[alpha]] <- outCJS[[delta_1]][[delta_2]][[alpha]] <- list()
    if(delta_combs[delt,3] | (!delta_combs[delt,3] & alpha==0)){
      #if(alpha==0 | alpha==1) {
      #  if(alpha==0) data_type <- "never"
      #  else data_type <- "always"
      #  parms <- c("pbeta","phibeta","delta","psi","loglike")
      if(alpha==0){
        data_type <- "never"
        parms <- c("pbeta","phibeta","delta","psi","loglike")
      } else if(alpha==1){
        data_type <- "always"
        parms <- c("pbeta","phibeta","delta","psi","loglike")
      } else {
        data_type <- "sometimes"
        parms <- c("pbeta","phibeta","delta","alpha","psi","loglike")
      }
      if(sigma2_zphi){ 
        mod.phi = ~h
        parms = c(parms,"sigma2_zphi")
      } else {
        mod.phi = ~ 1
      }
      data <- list()
      set.seed(1,kind="Mersenne-Twister",normal.kind = "Inversion")
      for(isim in 1:nsims){#isim=1#
        data[[isim]]=simdataCJS(N=N,noccas=noccas,pbeta=pbeta,sigma2_zphi=sigma2_zphi,phibeta=phibeta,delta_1=delta_1,delta_2=delta_2,alpha=alpha,data.type=data_type)##fun_sim_data(N,noccas,truemu,truetau,truesigma2_z,truedelta_1,truedelta_2,truealpha,data.type=data_type,link="logit")#simdata[[data_type]][[isim]]#
      }
      
      out <- outCJS <- NULL
      
      phibetahat=phibetahatCJS<-matrix(NA,nrow=nsims,ncol=7)
      pbetahat=pbetahatCJS<-matrix(NA,nrow=nsims,ncol=7)
      if(sigma2_zphi) sig2zphihat=sig2zphihatm0<-matrix(NA,nrow=nsims,ncol=7)
      phihat=phihatCJS<-matrix(NA,nrow=nsims,ncol=7)
      phat=phatCJS<-matrix(NA,nrow=nsims,ncol=7)
      
      psihat<-matrix(NA,nrow=nsims,ncol=6)
      if(data_type=="sometimes") alphahat<-matrix(NA,nrow=nsims,ncol=6)
      delta1hat<-matrix(NA,nrow=nsims,ncol=6)
      delta2hat<-matrix(NA,nrow=nsims,ncol=6)
      
      psihatCJS<-matrix(NA,nrow=nsims,ncol=6)
      if(data_type=="sometimes") alphahatCJS<-matrix(NA,nrow=nsims,ncol=6)
      delta1hatCJS<-matrix(NA,nrow=nsims,ncol=6)
      delta2hatCJS<-matrix(NA,nrow=nsims,ncol=6)
      
      for(isim in 1:nsims){#isim=1#
        
        #data=simdataCJS(N=N,noccas=noccas,pbeta=pbeta,phibeta=phibeta,data.type=data_type)##fun_sim_data(N,noccas,truemu,truetau,truesigma2_z,truedelta_1,truedelta_2,truealpha,data.type=data_type,link="logit")#simdata[[data_type]][[isim]]#
        Enc.Mat=data[[isim]]$Enc.Mat
        tEnc.Mat=data[[isim]]$trueEnc.Mat
        setup=processdata(Enc.Mat,data.type=data_type)
        
        
        tEnc.Mat[which(tEnc.Mat>0)] <- 1
        tEnc.Mat=tEnc.Mat[apply(tEnc.Mat>0,1,sum)>0,]
        ch<-as.character(as.matrix( apply(tEnc.Mat, 1, paste, collapse = ""), ncol=1 ))
        CH<-RMark::process.data(data.frame(ch,stringsAsFactors = FALSE),model="CJS")
        #tEnc.Mat1=tEnc.Mat[apply(tEnc.Mat>0,1,sum)>0,] # true encounter histories with at least 1 detection
        #N1 <- nrow(tEnc.Mat1) # true number of individuals with at least 1 detection
        #inits1 <- list(Freq=getfreq(rbind(tEnc.Mat1,matrix(0,ncol=noccas,nrow=max(0,nrow(Enc.Mat)-N1))),data.type=data_type))
        #markout<-RMark::mark(CH,output=FALSE,profile.int=TRUE,delete=TRUE,silent=TRUE)
        #phihatMARK[isim,]<-as.numeric(markout$results$real[1,1:4])
        #phatMARK[isim,]<-as.numeric(markout$results$real[2,1:4])  
        xvec <- multimark:::getfreq(data[[isim]]$trueEnc.Mat,setup@vAll.hists,data_type)
        xvec[1] <- nrow(data[[isim]]$Enc.Mat)-nrow(data[[isim]]$trueEnc.Mat)
        
        initial.values <- vector("list",nchains)
        if((1-(delta_1+delta_2))==0){
          if(delta_1==delta_2){
            tmpdelta_1 <- 0.5 - 1.e-12
            tmpdelta_2 <- 0.5 - 1.e-12
          } else {
            tmpdelta_1 <- min(1-1.e-6-1.e-6,max(1.e-6,delta_1))
            tmpdelta_2 <- min(1-1.e-6-1.e-6,max(1.e-6,delta_2))
          }
        } else {
          tmpdelta_1 <- min(1-1.e-6,max(1.e-6,delta_1))
          tmpdelta_2 <- min(1-1.e-6,max(1.e-6,delta_2))
        }
        initial.values[[1]]=list(pbeta=pbeta,phibeta=phibeta,delta_1=tmpdelta_1,delta_2=tmpdelta_2,alpha=min(1-1.e-6,max(1.e-6,alpha)),x=xvec,H=multimark:::get_H(setup,xvec))
        if(sigma2_zphi) initial.values[[1]]$sigma2_zphi <- sigma2_zphi
        for(i in 2:nchains){
          initial.values[[i]] <-  list(pbeta=pbeta,phibeta=phibeta,delta_1=tmpdelta_1,delta_2=tmpdelta_2,alpha=min(1-1.e-6,max(1.e-6,alpha)),x=xvec,H=multimark:::get_H(setup,xvec))
        }
        out[[isim]] <- grb <- tryCatch(stop(),error=function(e) e)
        while(inherits(out[[isim]],"error") | inherits(grb,"error")){
          out[[isim]]=tryCatch(suppressMessages(multimarkCJS(mms=setup,mod.p=~1,mod.phi=mod.phi,useInitial=useInitial,nchains=nchains,iter=niter,burnin=nburnin,parms=parms,initial.values = initial.values)),error=function(e) e)
          #print(summary(out[[isim]]))
          grb=tryCatch(coda::gelman.diag(out[[isim]]$mcmc,transform=TRUE,autoburnin=FALSE),error=function(e) e)
        }
        #print(grb)
        posphibeta=unlist(out[[isim]]$mcmc[,"phibeta[(Intercept)]"])
        posphi=pnorm(posphibeta)
        pospbeta<-unlist(out[[isim]]$mcmc[,"pbeta[(Intercept)]"])
        posp=pnorm(pospbeta)
        if(sigma2_zphi) possig2zp=unlist(out[[isim]]$mcmc[,"sigma2_zphi"])
        pospsi= unlist(out[[isim]]$mcmc[,"psi"])
        if(data_type=="sometimes") posalpha= unlist(out[[isim]]$mcmc[,"alpha"])
        posdelta1= unlist(out[[isim]]$mcmc[,"delta_1"])
        posdelta2= unlist(out[[isim]]$mcmc[,"delta_2"])
        phibetahat[isim,]=c(mean(posphibeta),median(posphibeta),sd(posphibeta),quantile(posphibeta,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(posphibeta))  
        phihat[isim,]=c(mean(posphi),median(posphi),sd(posphi),quantile(posphi,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(posphi))  
        #rm(posphi)
        pbetahat[isim,]=c(mean(pospbeta),median(pospbeta),sd(pospbeta),quantile(pospbeta,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(pospbeta))  
        phat[isim,]=c(mean(posp),median(posp),sd(posp),quantile(posp,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(posp))  
        if(sigma2_zphi) sig2zphihat[isim,]=c(mean(possig2zp),median(possig2zp),sd(possig2zp),quantile(possig2zp,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(possig2zp))
        psihat[isim,]=c(mean(pospsi),median(pospsi),sd(pospsi),quantile(pospsi,c(.025,.975),names=F),coda::effectiveSize(pospsi))  
        if(data_type=="sometimes") alphahat[isim,]=c(mean(posalpha),median(posalpha),sd(posalpha),quantile(posalpha,c(.025,.975),names=F),coda::effectiveSize(posalpha))  
        delta1hat[isim,]=c(mean(posdelta1),median(posdelta1),sd(posdelta1),quantile(posdelta1,c(.025,.975),names=F),coda::effectiveSize(posdelta1))  
        delta2hat[isim,]=c(mean(posdelta2),median(posdelta2),sd(posdelta2),quantile(posdelta2,c(.025,.975),names=F),coda::effectiveSize(posdelta2))  
        
        outCJS[[isim]] <- grb <- tryCatch(stop(),error=function(e) e)
        while(inherits(outCJS[[isim]],"error") | inherits(grb,"error")){
          outCJS[[isim]] <- tryCatch(suppressMessages(multimarkCJS(Enc.Mat = rbind(data[[isim]]$trueEnc.Mat,matrix(0,nrow(data[[isim]]$Enc.Mat)-nrow(data[[isim]]$trueEnc.Mat),ncol=noccas)),data.type=data_type,mod.p=~1,mod.phi=mod.phi,useInitial=useInitial,nchains=2,iter=niterCJS,burnin=nburninCJS,parms=parms,known=c(rep(1,nrow(data[[isim]]$trueEnc.Mat)),rep(0,nrow(data[[isim]]$Enc.Mat)-nrow(data[[isim]]$trueEnc.Mat))))),error=function(e) e)
          #print(summary(out[[isim]]))
          grb=tryCatch(coda::gelman.diag(outCJS[[isim]]$mcmc,transform=TRUE,autoburnin=FALSE),error=function(e) e)
        }
        posphibeta=unlist(outCJS[[isim]]$mcmc[,"phibeta[(Intercept)]"])
        posphi=pnorm(posphibeta)
        pospbeta<-unlist(outCJS[[isim]]$mcmc[,"pbeta[(Intercept)]"])
        posp=pnorm(pospbeta)
        if(sigma2_zphi) possig2zp=unlist(outCJS[[isim]]$mcmc[,"sigma2_zphi"])
        pospsi= unlist(outCJS[[isim]]$mcmc[,"psi"])
        if(data_type=="sometimes") posalpha= unlist(outCJS[[isim]]$mcmc[,"alpha"])
        posdelta1= unlist(outCJS[[isim]]$mcmc[,"delta_1"])
        posdelta2= unlist(outCJS[[isim]]$mcmc[,"delta_2"])
        #print(grb)
        phibetahatCJS[isim,]=c(mean(posphibeta),median(posphibeta),sd(posphibeta),quantile(posphibeta,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(posphibeta))
        phihatCJS[isim,]=c(mean(posphi),median(posphi),sd(posphi),quantile(posphi,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(posphi))  
        pbetahatCJS[isim,]=c(mean(pospbeta),median(pospbeta),sd(pospbeta),quantile(pospbeta,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(pospbeta))  
        phatCJS[isim,]=c(mean(posp),median(posp),sd(posp),quantile(posp,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(posp))  
        if(sigma2_zphi) sig2zphihatm0[isim,]=c(mean(possig2zp),median(possig2zp),sd(possig2zp),quantile(possig2zp,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(possig2zp))
        psihatCJS[isim,]=c(mean(pospsi),median(pospsi),sd(pospsi),quantile(pospsi,c(.025,.975),names=F),coda::effectiveSize(pospsi))  
        if(data_type=="sometimes") alphahatCJS[isim,]=c(mean(posalpha),median(posalpha),sd(posalpha),quantile(posalpha,c(.025,.975),names=F),coda::effectiveSize(posalpha))  
        delta1hatCJS[isim,]=c(mean(posdelta1),median(posdelta1),sd(posdelta1),quantile(posdelta1,c(.025,.975),names=F),coda::effectiveSize(posdelta1))  
        delta2hatCJS[isim,]=c(mean(posdelta2),median(posdelta2),sd(posdelta2),quantile(posdelta2,c(.025,.975),names=F),coda::effectiveSize(posdelta2))  
        
        #rm(posphi)
        
        cat("delta_1 ",delta_1," delta_2 ",delta_2," alpha ",alpha," Simulation ",isim,"\n")
        #cat("Simulation ",isim," phi ",phihat[isim,],phihatCJS[isim,],"\n")
        #cat("Simulation ",isim," phi ",phihat[isim,],phihatCJS[isim,],"\n")
        #cat("Simulation ",isim," phibeta ",phibetahat[isim,],phibetahatCJS[isim,],"\n")
        #cat("Simulation ",isim," p ",phat[isim,],phatCJS[isim,],"\n")
        #cat("Simulation ",isim," pbeta ",pbetahat[isim,],pbetahatCJS[isim,],"\n\n")
        cat("Simulation psi ",apply(psihat[1:isim,,drop=FALSE],2,mean),apply(psihatCJS[1:isim,,drop=FALSE],2,mean),"\n")
        if(data_type=="sometimes") cat("Simulation alpha ",apply(alphahat[1:isim,,drop=FALSE],2,mean),apply(alphahatCJS[1:isim,,drop=FALSE],2,mean),"\n")
        cat("Simulation delta_1 ",apply(delta1hat[1:isim,,drop=FALSE],2,mean),apply(delta1hatCJS[1:isim,,drop=FALSE],2,mean),"\n")
        cat("Simulation delta_2",apply(delta2hat[1:isim,,drop=FALSE],2,mean),apply(delta2hatCJS[1:isim,,drop=FALSE],2,mean),"\n")
        cat("Simulation phi ",apply(phihat[1:isim,,drop=FALSE],2,mean),apply(phihatCJS[1:isim,,drop=FALSE],2,mean),"\n")
        cat("Simulation p ",apply(phat[1:isim,,drop=FALSE],2,mean),apply(phatCJS[1:isim,,drop=FALSE],2,mean),"\n")
        cat("Simulation phibeta ",apply(phibetahat[1:isim,,drop=FALSE],2,mean),apply(phibetahatCJS[1:isim,,drop=FALSE],2,mean),"\n")
        if(sigma2_zphi) cat("Simulation sigma2_zphi ",apply(sig2zphihat[1:isim,,drop=FALSE],2,mean),apply(sig2zphihatm0[1:isim,,drop=FALSE],2,mean),"\n")
        cat("Simulation pbeta ",apply(pbetahat[1:isim,,drop=FALSE],2,mean),apply(pbetahatCJS[1:isim,,drop=FALSE],2,mean),"\n\n")
        
      }
      
      if(data_type=="sometimes") {
        if(sigma2_zphi){
          save(alphahat,alphahatCJS,delta1hat,delta1hatCJS,delta2hat,delta2hatCJS,pbetahat,pbetahatCJS,phat,phatCJS,phibetahat,phibetahatCJS,phihat,phihatCJS,psihat,psihatCJS,sig2zphihat,sig2zphihatm0,
               noccas,N,pbeta,phibeta,
               alpha,delta_1,delta_2,sigma2_zphi,
               data_type,niter,nburnin,niterCJS,nburninCJS,nchains,nsims,
               file=paste0(getwd(),"/multimark_projects/development files/CJSims/sims_noccas",noccas,"_N",N,"_pbeta",pbeta,"_sigma2_zphi",sigma2_zphi,"_phibeta",phibeta,"_delta1",round(delta_1,2),"_delta2",round(delta_2,2),"_alpha",round(alpha,2),"_prodComb",ifelse(useInitial,"","_nouseInitial"),".RData"))
        } else {
          save(alphahat,alphahatCJS,delta1hat,delta1hatCJS,delta2hat,delta2hatCJS,pbetahat,pbetahatCJS,phat,phatCJS,phibetahat,phibetahatCJS,phihat,phihatCJS,psihat,psihatCJS,
               noccas,N,pbeta,phibeta,
               alpha,delta_1,delta_2,sigma2_zphi,
               data_type,niter,nburnin,niterCJS,nburninCJS,nchains,nsims,
               file=paste0(getwd(),"/multimark_projects/development files/CJSims/sims_noccas",noccas,"_N",N,"_pbeta",pbeta,"_sigma2_zphi",sigma2_zphi,"_phibeta",phibeta,"_delta1",round(delta_1,2),"_delta2",round(delta_2,2),"_alpha",round(alpha,2),"_prodComb",ifelse(useInitial,"","_nouseInitial"),".RData"))
        }
      } else {
        if(sigma2_zphi){
          save(delta1hat,delta1hatCJS,delta2hat,delta2hatCJS,pbetahat,pbetahatCJS,phat,phatCJS,phibetahat,phibetahatCJS,phihat,phihatCJS,psihat,psihatCJS,sig2zphihat,sig2zphihatm0,
               noccas,N,pbeta,phibeta,
               alpha,delta_1,delta_2,sigma2_zphi,
               data_type,niter,nburnin,niterCJS,nburninCJS,nchains,nsims,
               file=paste0(getwd(),"/multimark_projects/development files/CJSims/sims_noccas",noccas,"_N",N,"_pbeta",pbeta,"_sigma2_zphi",sigma2_zphi,"_phibeta",phibeta,"_delta1",round(delta_1,2),"_delta2",round(delta_2,2),"_alpha",round(alpha,2),"_prodComb",ifelse(useInitial,"","_nouseInitial"),".RData"))
        } else {
          save(delta1hat,delta1hatCJS,delta2hat,delta2hatCJS,pbetahat,pbetahatCJS,phat,phatCJS,phibetahat,phibetahatCJS,phihat,phihatCJS,psihat,psihatCJS,
               noccas,N,pbeta,phibeta,
               alpha,delta_1,delta_2,sigma2_zphi,
               data_type,niter,nburnin,niterCJS,nburninCJS,nchains,nsims,
               file=paste0(getwd(),"/multimark_projects/development files/CJSims/sims_noccas",noccas,"_N",N,"_pbeta",pbeta,"_sigma2_zphi",sigma2_zphi,"_phibeta",phibeta,"_delta1",round(delta_1,2),"_delta2",round(delta_2,2),"_alpha",round(alpha,2),"_prodComb",ifelse(useInitial,"","_nouseInitial"),".RData"))
        }        
      }
      
    }
  }
  #  }
}

for(delt in 1:nrow(delta_combs)){
  delta_1 <- delta_combs[delt,1]
  delta_2 <- delta_combs[delt,2]

  for(alpha in truealpha){
    
    if(delta_combs[delt,3] | (!delta_combs[delt,3] & alpha==0)){
      
      if(alpha==0){
        data_type <- "never"
        parms <- c("pbeta","phibeta","delta","psi")
      } else {
        data_type <- "sometimes"
        parms <- c("pbeta","phibeta","delta","alpha","psi")
      }
      
      if(data_type=="sometimes") {
        load(paste0(getwd(),"/multimark_projects/development files/CJSims/sims_noccas",noccas,"_N",N,"_pbeta",pbeta,"_sigma2_zphi",sigma2_zphi,"_phibeta",phibeta,"_delta1",round(delta_1,2),"_delta2",round(delta_2,2),"_alpha",round(alpha,2),"_prodComb",ifelse(useInitial,"","_nouseInitial"),".RData"))
      } else {
        load(paste0(getwd(),"/multimark_projects/development files/CJSims/sims_noccas",noccas,"_N",N,"_pbeta",pbeta,"_sigma2_zphi",sigma2_zphi,"_phibeta",phibeta,"_delta1",round(delta_1,2),"_delta2",round(delta_2,2),"_alpha",round(alpha,2),"_prodComb",ifelse(useInitial,"","_nouseInitial"),".RData"))
        
      }
      
      cat("delta_1 ",delta_1," delta_2 ",delta_2," alpha ",alpha,"\n")
      #cat("Simulation ",isim," phi ",phihat[isim,],phihatCJS[isim,],"\n")
      #cat("Simulation ",isim," phi ",phihat[isim,],phihatCJS[isim,],"\n")
      #cat("Simulation ",isim," phibeta ",phibetahat[isim,],phibetahatCJS[isim,],"\n")
      #cat("Simulation ",isim," p ",phat[isim,],phatCJS[isim,],"\n")
      #cat("Simulation ",isim," pbeta ",pbetahat[isim,],pbetahatCJS[isim,],"\n\n")
      cat("Simulation psi ",apply(psihat,2,mean),apply(psihatCJS,2,mean),"\n")
      if(data_type=="sometimes") cat("Simulation alpha ",apply(alphahat,2,mean),apply(alphahatCJS,2,mean),"\n")
      cat("Simulation delta_1 ",apply(delta1hat,2,mean),apply(delta1hatCJS,2,mean),"\n")
      cat("Simulation delta_2",apply(delta2hat,2,mean),apply(delta2hatCJS,2,mean),"\n")
      cat("Simulation phi ",apply(phihat,2,mean),apply(phihatCJS,2,mean),"\n")
      cat("Simulation p ",apply(phat,2,mean),apply(phatCJS,2,mean),"\n")
      if(sigma2_zphi) cat("Simulation sigma2_zphi ",apply(sig2zphihat,2,mean),apply(sig2zphihatm0,2,mean),"\n")
      cat("Simulation phibeta ",apply(phibetahat,2,mean),apply(phibetahatCJS,2,mean),"\n")
      cat("Simulation pbeta ",apply(pbetahat,2,mean),apply(pbetahatCJS,2,mean),"\n\n")
      if(sigma2_zphi) cat("Overall PRB sig2zp ",mean((sig2zphihat[,1]-sigma2_zphi)/sigma2_zphi*100),mean((sig2zphihatm0[,1]-sigma2_zphi)/sigma2_zphi*100),"\n")
      cat("Overall bias phibeta ",mean((phibetahat[,1]-phibeta)),mean((phibetahatCJS[,1]-phibeta)),"\n")
      cat("Overall PRB phi ",mean((phihat[,1]-pnorm(phibeta))/pnorm(phibeta)*100),mean((phihatCJS[,1]-pnorm(phibeta))/pnorm(phibeta)*100),"\n\n")
      
    }
  }
}
