library(multimark)
library(RMark)

nchains <- 10
niter <- 1000000
nadapt <- 100000
nbin <- 1000
nburnin <- 500000
niterm0 <- 50000
nadaptm0 <- 10000
nbinm0 <- 100
nburninm0 <- 10000

noccas=7

N=100
pbeta = -0.25
sigma2_zp = 0.5
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

#out=outm0<-NULL

#for(delta_1 in truedelta_1){
#out[[delta_1]] <- outm0[[delta_1]] <- list()
#  for(delta_2 in truedelta_2){
#out[[delta_1]][[delta_2]] <- outm0[[delta_1]][[delta_2]] <- list()
for(delt in 1:nrow(delta_combs)){
  delta_1 <- delta_combs[delt,1]
  delta_2 <- delta_combs[delt,2]
  for(alpha in truealpha){
    
    #out[[delta_1]][[delta_2]][[alpha]] <- outm0[[delta_1]][[delta_2]][[alpha]] <- list()
    if(delta_combs[delt,3] | (!delta_combs[delt,3] & alpha==0)){
      #if(alpha==0 | alpha==1) {
      #  if(alpha==0) data_type <- "never"
      #  else data_type <- "always"
      #  parms <- c("pbeta","N","delta","psi","logPosterior")
      if(alpha==0){
        data_type <- "never"
        parms <- c("pbeta","N","delta","psi","logPosterior")
      } else if(alpha==1){
        data_type <- "always"
        parms <- c("pbeta","N","delta","psi","logPosterior")
      } else {
        data_type <- "sometimes"
        parms <- c("pbeta","N","delta","alpha","psi","logPosterior")
      }
      if(sigma2_zp){ 
        mod.p = ~h
        parms = c(parms,"sigma2_zp")
      } else {
        mod.p = ~ 1
      }
      data <- list()
      set.seed(1,kind="Mersenne-Twister",normal.kind = "Inversion")
      for(isim in 1:nsims){#isim=1#
        data[[isim]]=simdataClosed(N=N,noccas=noccas,pbeta=pbeta,sigma2_zp=sigma2_zp,delta_1=delta_1,delta_2=delta_2,alpha=alpha,data.type=data_type)##fun_sim_data(N,noccas,truemu,truetau,truesigma2_z,truedelta_1,truedelta_2,truealpha,data.type=data_type,link="logit")#simdata[[data_type]][[isim]]#
      }
      
      out <- outm0 <- NULL
      
      Nhat=Nhatm0<-matrix(NA,nrow=nsims,ncol=7)
      pbetahat=pbetahatm0<-matrix(NA,nrow=nsims,ncol=7)
      if(sigma2_zp) sig2zphat=sig2zphatm0<-matrix(NA,nrow=nsims,ncol=7)
      phat=phatm0<-matrix(NA,nrow=nsims,ncol=7)
      
      psihat<-matrix(NA,nrow=nsims,ncol=6)
      if(data_type=="sometimes") alphahat<-matrix(NA,nrow=nsims,ncol=6)
      delta1hat<-matrix(NA,nrow=nsims,ncol=6)
      delta2hat<-matrix(NA,nrow=nsims,ncol=6)
      
      psihatm0<-matrix(NA,nrow=nsims,ncol=6)
      if(data_type=="sometimes") alphahatm0<-matrix(NA,nrow=nsims,ncol=6)
      delta1hatm0<-matrix(NA,nrow=nsims,ncol=6)
      delta2hatm0<-matrix(NA,nrow=nsims,ncol=6)
      
      for(isim in 1:nsims){#isim=1#
        
        Enc.Mat=data[[isim]]$Enc.Mat
        tEnc.Mat=data[[isim]]$trueEnc.Mat
        setup=processdata(Enc.Mat,data.type=data_type)
        
        
        tEnc.Mat[which(tEnc.Mat>0)] <- 1
        tEnc.Mat=tEnc.Mat[apply(tEnc.Mat>0,1,sum)>0,]
        ch<-as.character(as.matrix( apply(tEnc.Mat, 1, paste, collapse = ""), ncol=1 ))
        xvec <- multimark:::getfreq(data[[isim]]$trueEnc.Mat,setup@vAll.hists,data_type)
        xvec[1] <- nrow(data[[isim]]$Enc.Mat)-nrow(data[[isim]]$trueEnc.Mat)+xvec[1]
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
        initial.values[[1]]=list(pbeta=pbeta,N=N,delta_1=tmpdelta_1,delta_2=tmpdelta_2,alpha=min(1-1.e-6,max(1.e-6,alpha)),x=xvec,H=multimark:::get_H(setup,xvec))
        for(i in 2:nchains){
          initial.values[[i]] <-  list(pbeta=pbeta,N=N,delta_1=tmpdelta_1,delta_2=tmpdelta_2,alpha=min(1-1.e-6,max(1.e-6,alpha)),x=xvec,H=multimark:::get_H(setup,xvec))
        }
        
        out[[isim]]=suppressMessages(multimarkClosed(mms=setup,mod.p=mod.p,nchains=nchains,iter=niter,adapt=nadapt,bin=nbin,burnin=nburnin,parms=parms))
        posN=unlist(out[[isim]]$mcmc[,"N"])
        pospbeta<-unlist(out[[isim]]$mcmc[,"pbeta[(Intercept)]"])
        posp=multimark:::expit(pospbeta)
        if(sigma2_zp) possig2zp=unlist(out[[isim]]$mcmc[,"sigma2_zp"])
        pospsi= unlist(out[[isim]]$mcmc[,"psi"])
        if(data_type=="sometimes") posalpha= unlist(out[[isim]]$mcmc[,"alpha"])
        posdelta1= unlist(out[[isim]]$mcmc[,"delta_1"])
        posdelta2= unlist(out[[isim]]$mcmc[,"delta_2"])
        #print(summary(out[[isim]]))
        grb=coda::gelman.diag(out[[isim]]$mcmc,transform=TRUE,autoburnin=FALSE)
        #print(grb)
        Nhat[isim,]=c(mean(posN),median(posN),sd(posN),quantile(posN,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(posN))  
        pbetahat[isim,]=c(mean(pospbeta),median(pospbeta),sd(pospbeta),quantile(pospbeta,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(pospbeta))  
        phat[isim,]=c(mean(posp),median(posp),sd(posp),quantile(posp,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(posp))  
        if(sigma2_zp) sig2zphat[isim,]=c(mean(possig2zp),median(possig2zp),sd(possig2zp),quantile(possig2zp,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(possig2zp))
        psihat[isim,]=c(mean(pospsi),median(pospsi),sd(pospsi),quantile(pospsi,c(.025,.975),names=F),coda::effectiveSize(pospsi))  
        if(data_type=="sometimes") alphahat[isim,]=c(mean(posalpha),median(posalpha),sd(posalpha),quantile(posalpha,c(.025,.975),names=F),coda::effectiveSize(posalpha))  
        delta1hat[isim,]=c(mean(posdelta1),median(posdelta1),sd(posdelta1),quantile(posdelta1,c(.025,.975),names=F),coda::effectiveSize(posdelta1))  
        delta2hat[isim,]=c(mean(posdelta2),median(posdelta2),sd(posdelta2),quantile(posdelta2,c(.025,.975),names=F),coda::effectiveSize(posdelta2))  
        
        outm0[[isim]] <- suppressMessages(multimarkClosed(Enc.Mat = rbind(data[[isim]]$trueEnc.Mat[which(apply(data[[isim]]$trueEnc.Mat,1,sum)>0),],matrix(0,max(0,nrow(data[[isim]]$Enc.Mat)-nrow(data[[isim]]$trueEnc.Mat[which(apply(data[[isim]]$trueEnc.Mat,1,sum)>0),])),ncol=noccas)),data.type=data_type,mod.p=mod.p,nchains=2,iter=niterm0,adapt=nadaptm0,bin=nbinm0,burnin=nburninm0,parms=parms,known=c(rep(1,nrow(data[[isim]]$trueEnc.Mat[which(apply(data[[isim]]$trueEnc.Mat,1,sum)>0),])),rep(0,max(0,nrow(data[[isim]]$Enc.Mat)-nrow(data[[isim]]$trueEnc.Mat[which(apply(data[[isim]]$trueEnc.Mat,1,sum)>0),]))))))
        posN=unlist(outm0[[isim]]$mcmc[,"N"])
        pospbeta<-unlist(outm0[[isim]]$mcmc[,"pbeta[(Intercept)]"])
        posp=multimark:::expit(pospbeta)
        if(sigma2_zp) possig2zp=unlist(outm0[[isim]]$mcmc[,"sigma2_zp"])
        pospsi= unlist(outm0[[isim]]$mcmc[,"psi"])
        if(data_type=="sometimes") posalpha= unlist(outm0[[isim]]$mcmc[,"alpha"])
        posdelta1= unlist(outm0[[isim]]$mcmc[,"delta_1"])
        posdelta2= unlist(outm0[[isim]]$mcmc[,"delta_2"])
        
        #print(summary(out[[isim]]))
        grb=coda::gelman.diag(outm0[[isim]]$mcmc,transform=TRUE,autoburnin=FALSE)
        #print(grb)
        Nhatm0[isim,]=c(mean(posN),median(posN),sd(posN),quantile(posN,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(posN))
        pbetahatm0[isim,]=c(mean(pospbeta),median(pospbeta),sd(pospbeta),quantile(pospbeta,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(pospbeta))  
        phatm0[isim,]=c(mean(posp),median(posp),sd(posp),quantile(posp,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(posp))  
        if(sigma2_zp) sig2zphatm0[isim,]=c(mean(possig2zp),median(possig2zp),sd(possig2zp),quantile(possig2zp,c(.025,.975),names=F),grb$mpsrf,coda::effectiveSize(possig2zp))
        psihatm0[isim,]=c(mean(pospsi),median(pospsi),sd(pospsi),quantile(pospsi,c(.025,.975),names=F),coda::effectiveSize(pospsi))  
        if(data_type=="sometimes") alphahatm0[isim,]=c(mean(posalpha),median(posalpha),sd(posalpha),quantile(posalpha,c(.025,.975),names=F),coda::effectiveSize(posalpha))  
        delta1hatm0[isim,]=c(mean(posdelta1),median(posdelta1),sd(posdelta1),quantile(posdelta1,c(.025,.975),names=F),coda::effectiveSize(posdelta1))  
        delta2hatm0[isim,]=c(mean(posdelta2),median(posdelta2),sd(posdelta2),quantile(posdelta2,c(.025,.975),names=F),coda::effectiveSize(posdelta2))  
        
        #rm(posphi)
        
        cat("delta_1 ",delta_1," delta_2 ",delta_2," alpha ",alpha," Simulation ",isim,"\n")
        cat("Simulation psi ",apply(psihat[1:isim,,drop=FALSE],2,mean),apply(psihatm0[1:isim,,drop=FALSE],2,mean),"\n")
        if(data_type=="sometimes") cat("Simulation alpha ",apply(alphahat[1:isim,,drop=FALSE],2,mean),apply(alphahatm0[1:isim,,drop=FALSE],2,mean),"\n")
        cat("Simulation delta_1 ",apply(delta1hat[1:isim,,drop=FALSE],2,mean),apply(delta1hatm0[1:isim,,drop=FALSE],2,mean),"\n")
        cat("Simulation delta_2",apply(delta2hat[1:isim,,drop=FALSE],2,mean),apply(delta2hatm0[1:isim,,drop=FALSE],2,mean),"\n")
        cat("Simulation p ",apply(phat[1:isim,,drop=FALSE],2,mean),apply(phatm0[1:isim,,drop=FALSE],2,mean),"\n")
        cat("Simulation pbeta ",apply(pbetahat[1:isim,,drop=FALSE],2,mean),apply(pbetahatm0[1:isim,,drop=FALSE],2,mean),"\n")
        if(sigma2_zp) cat("Simulation sigma2_zp ",apply(sig2zphat[1:isim,,drop=FALSE],2,mean),apply(sig2zphatm0[1:isim,,drop=FALSE],2,mean),"\n")
        cat("Simulation N ",apply(Nhat[1:isim,,drop=FALSE],2,mean),apply(Nhatm0[1:isim,,drop=FALSE],2,mean),"\n\n")
        
        
      }
      
      if(data_type=="sometimes") {
        if(sigma2_zp){
          save(alphahat,alphahatm0,delta1hat,delta1hatm0,delta2hat,delta2hatm0,pbetahat,pbetahatm0,phat,phatm0,Nhat,Nhatm0,psihat,psihatm0,sig2zphat,sig2zphatm0,
               noccas,N,pbeta,
               alpha,delta_1,delta_2,sigma2_zp,
               data_type,niter,nburnin,niterm0,nburninm0,nchains,nsims,
               file=paste0(getwd(),"/multimark_projects/development files/closedSims/sims_noccas",noccas,"_N",N,"_pbeta",pbeta,"_sig2zp",sigma2_zp,"_delta1",round(delta_1,2),"_delta2",round(delta_2,2),"_alpha",round(alpha,2),"_prod.RData"))
        } else {
          save(alphahat,alphahatm0,delta1hat,delta1hatm0,delta2hat,delta2hatm0,pbetahat,pbetahatm0,phat,phatm0,Nhat,Nhatm0,psihat,psihatm0,
               noccas,N,pbeta,
               alpha,delta_1,delta_2,sigma2_zp,
               data_type,niter,nburnin,niterm0,nburninm0,nchains,nsims,
               file=paste0(getwd(),"/multimark_projects/development files/closedSims/sims_noccas",noccas,"_N",N,"_pbeta",pbeta,"_sig2zp",sigma2_zp,"_delta1",round(delta_1,2),"_delta2",round(delta_2,2),"_alpha",round(alpha,2),"_prod.RData"))
        }
      } else {
        if(sigma2_zp){
          save(delta1hat,delta1hatm0,delta2hat,delta2hatm0,pbetahat,pbetahatm0,phat,phatm0,Nhat,Nhatm0,psihat,psihatm0,sig2zphat,sig2zphatm0,
               noccas,N,pbeta,
               alpha,delta_1,delta_2,sigma2_zp,
               data_type,niter,nburnin,niterm0,nburninm0,nchains,nsims,
               file=paste0(getwd(),"/multimark_projects/development files/closedSims/sims_noccas",noccas,"_N",N,"_pbeta",pbeta,"_sig2zp",sigma2_zp,"_delta1",round(delta_1,2),"_delta2",round(delta_2,2),"_alpha",round(alpha,2),"_prod.RData"))
        } else {
          save(delta1hat,delta1hatm0,delta2hat,delta2hatm0,pbetahat,pbetahatm0,phat,phatm0,Nhat,Nhatm0,psihat,psihatm0,
               noccas,N,pbeta,
               alpha,delta_1,delta_2,sigma2_zp,
               data_type,niter,nburnin,niterm0,nburninm0,nchains,nsims,
               file=paste0(getwd(),"/multimark_projects/development files/closedSims/sims_noccas",noccas,"_N",N,"_pbeta",pbeta,"_sig2zp",sigma2_zp,"_delta1",round(delta_1,2),"_delta2",round(delta_2,2),"_alpha",round(alpha,2),"_prod.RData"))
          
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
        parms <- c("pbeta","N","delta","psi")
      } else {
        data_type <- "sometimes"
        parms <- c("pbeta","N","delta","alpha","psi")
      }
      
      if(data_type=="sometimes") {
        load(paste0(getwd(),"/multimark_projects/development files/closedSims/sims_noccas",noccas,"_N",N,"_pbeta",pbeta,"_sig2zp",sigma2_zp,"_delta1",round(delta_1,2),"_delta2",round(delta_2,2),"_alpha",round(alpha,2),"_prod.RData"))
      } else {
        load(paste0(getwd(),"/multimark_projects/development files/closedSims/sims_noccas",noccas,"_N",N,"_pbeta",pbeta,"_sig2zp",sigma2_zp,"_delta1",round(delta_1,2),"_delta2",round(delta_2,2),"_alpha",round(alpha,2),"_prod.RData"))
        
      }
      
      cat("delta_1 ",delta_1," delta_2 ",delta_2," alpha ",alpha,"\n")
      cat("Simulation psi ",apply(psihat,2,mean),apply(psihatm0,2,mean),"\n")
      if(data_type=="sometimes") cat("Simulation alpha ",apply(alphahat,2,mean),apply(alphahatm0,2,mean),"\n")
      cat("Simulation delta_1 ",apply(delta1hat,2,mean),apply(delta1hatm0,2,mean),"\n")
      cat("Simulation delta_2",apply(delta2hat,2,mean),apply(delta2hatm0,2,mean),"\n")
      cat("Simulation p ",apply(phat,2,mean),apply(phatm0,2,mean),"\n")
      cat("Simulation pbeta ",apply(pbetahat,2,mean),apply(pbetahatm0,2,mean),"\n")
      if(sigma2_zp) cat("Simulation sigma2_zp ",apply(sig2zphat,2,mean),apply(sig2zphatm0,2,mean),"\n")
      cat("Simulation N ",apply(Nhat,2,mean),apply(Nhatm0,2,mean),"\n")
      cat("Overall PRB p ",mean((phat[,1]-multimark:::expit(pbeta))/multimark:::expit(pbeta)*100),mean((phatm0[,1]-multimark:::expit(pbeta))/multimark:::expit(pbeta)*100),"\n")
      cat("Overall PRB N ",mean((Nhat[,1]-N)/N*100),mean((Nhatm0[,1]-N)/N*100),"\n")
      if(sigma2_zp) cat("Overall PRB sig2zp ",mean((sig2zphat[,1]-sigma2_zp)/sigma2_zp*100),mean((sig2zphatm0[,1]-sigma2_zp)/sigma2_zp*100),"\n\n")
    }
  }
}
