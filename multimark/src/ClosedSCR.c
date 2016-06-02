#include "ranlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#define tol 1.e-6

#ifndef min
  #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

// Define function ClosedC to draw samples from the posterior distribution

void ClosedSCRC(int *ichain, double *mu0, double *sigma2_mu0, double *beta, double *z, double *sigma2_z, double *sigma2_scr, double *delta_1, double *delta_2, double *alpha, int *x, double *N, double *psi, int *H, int *center,
              int *ntraps, int *noccas, int *M, double *a0delta, double *a0alpha, double *b0alpha, double *A, double *a0psi, double *b0psi,
              double *Propsd, int *NNvect, int *numnn, double *accept, double *posterior,
              int *nHists, int *Allhists, int *C, int *indBasis, int *ncolBasis, int *knownx, double *DMp, double *DMc, int *pdim,
              int *iter, int *thin, int *adapt, int *bin, double *taccept, double *tuneadjust, int *numbasis,
              int *npoints, double *weights, double *nodes, int *mod_h, int *data_type, int *zind, int *Hind, int *centerind, int *updatedelta, int *delta_type, double *dexp, double *dist2, int *ncell, int *printlog)
{
  
  GetRNGstate(); 

  /* Declare functions */
  double invcloglog();
  double LIKESCR();
  double POSTERIORSCR();
  double FREQSUM();
  int GETCK();
  int sample();
  void PROPFREQ();
  void GETDELTA();
  double GETPSTARSCR();

  int T = *noccas;
  int K = *ntraps;
  int supN = *M; 
  int Mk = supN*T*K;      
  int datatype = *data_type;
  int deltatype = *delta_type;
    
  int niter, th, ada, n;
  int J = *nHists;
  int indhist;
  
  niter = *iter;              /* Number of iterations in the Markov chain */
  th = *thin;                 /* Number of iterations for thinning */
  ada = *adapt;               /* Number of iterations for pilot tuning */
  n = *bin;                   /* Number of iterations between each pilot tune */
  
  int dimp = *pdim;               /* length of beta vector */
  int ncells = *ncell;
  
  double betas[dimp], zs[supN], zstar, sigma2_zs, sigma_zstar, sigma2_scrs, sigma2_scrstar, sigma_scrstar, alphas, delta_1s, delta_2s, psis;
  double betastar[dimp];
  int xs[J], xnew[J], knownxs[J], Hs[supN], Hnew[supN], centers[supN];
  int qs[supN], qnew[supN]; 
  double Ns;
  double p[Mk], c[Mk];
  double sha, sca;
  int base, nbasis;
  double sqrttol = sqrt(tol);

  int t, i, j, k, g=0;
  double op, np;
  int dimrem, dimadd;
  
  for(k=0; k<dimp; k++){
    betas[k]=beta[k];
    betastar[k]=betas[k];
  }
  sigma2_zs= ( *mod_h ? sigma2_z[0] : (double) 0.0);
  double sigma_zs = sqrt(sigma2_zs);
  sigma2_scrs = sigma2_scr[0];
  sigma2_scrstar = sigma2_scrs;
  double sigma_scrs = sqrt(sigma2_scrs);
  sigma_scrstar = sigma_scrs;
  delta_1s= (*updatedelta ? delta_1[0] : 1.);
  delta_2s= (*updatedelta ? delta_2[0] : 0.);
  alphas= (*updatedelta ? alpha[0] : 0.);
  psis= (*updatedelta ? psi[0] : 1.);

  double ns=0.;
  for(i=0; i< supN; i++)  {
    zs[i] = ( *mod_h ? z[i] : (double) 0.0);
    Hs[i] = H[i];
    qs[i] = ((Hs[i]) ? 1 : 0);
    ns += (double) qs[i];
    qnew[i] = qs[i];
    Hnew[i] = Hs[i];
    centers[i] = center[i];
  }
  Ns = N[0];
  double Nstar = Ns;
  double nstar = ns;
  
  for(j=0; j<J; j++) {
    xs[j]=x[j];
    xnew[j]=xs[j];
    knownxs[j]= knownx[j];
  }
  int xind;
  
  int npts = *npoints;
  double weight[npts],node[npts];
  for(i=0; i<npts; i++){
    weight[i]=weights[i];
    node[i]=nodes[i];
  }
  
  double cloglogp[T*K], cloglogpstar[T*K], cloglogc[T*K], cloglogcstar[T*K];
  double propp[Mk], propc[Mk];
  for(k=0; k<K; k++){
    for(t=0; t<T; t++){
      cloglogp[k*T+t]=0.;
      cloglogc[k*T+t]=0.;
      for(j=0; j<dimp; j++){
        cloglogp[k*T+t]+=DMp[(k*T+t)*dimp+j]*betas[j];  
        cloglogc[k*T+t]+=DMc[(k*T+t)*dimp+j]*betas[j];  
      }
      cloglogpstar[k*T+t]=cloglogp[k*T+t];
      cloglogcstar[k*T+t]=cloglogc[k*T+t];
      for(i=0; i<supN; i++){
        p[i*T*K+k*T+t]=invcloglog((cloglogp[k*T+t]-1.0/(2.0*sigma2_scrs)*pow(dist2[k*ncells+centers[i]],*dexp)+zs[i]));
        c[i*T*K+k*T+t]=invcloglog((cloglogc[k*T+t]-1.0/(2.0*sigma2_scrs)*pow(dist2[k*ncells+centers[i]],*dexp)+zs[i]));
        //Rprintf("k %d t %d i %d indhist %d d2 %f p %f c %f \n",k,t,i,Allhists[Hs[i] * T*K + k*T+t],pow(dist2[k*ncells+centers[i]],*dexp),p[i*T*K+k*T+t],c[i*T*K+k*T+t]);
        propp[i*T*K+k*T+t]=p[i*T*K+k*T+t];
        propc[i*T*K+k*T+t]=c[i*T*K+k*T+t];
      }
    }
  }
  double pstar=GETPSTARSCR(dist2, cloglogp, sigma2_scrs, T, K, ncells, *dexp);
  Rprintf("pstar %f \n",pstar);
  double proppstar=pstar;
  
  double a0_delta[3], deltavect[3];
  for(i=0; i<3; i++){
    a0_delta[i]=a0delta[i];
  }
  deltavect[0]=delta_1s;
  deltavect[1]=delta_2s;
  deltavect[2]=1.-delta_1s-delta_2s;

  /* Vector to store numerator of acceptance rate for each parameter  */   
  double arate[(dimp+supN+1)];       
  for (j=0; j<(dimp+supN+1); j++) {
    arate[j]=0.0;
  }
  
  double indbase[*ncolBasis];
  for(j=0; j< *ncolBasis; j++){
    indbase[j]=0.0;
  }
  int c_k, temp;
  int ind, obasesum, nbasesum;
    
  /* Calculate the log-likelihood */  
  double ll=LIKESCR(p,c,qs,delta_1s,delta_2s,alphas,Allhists,Hs,T,K,supN,C,Ns,pstar);
  posterior[0]=POSTERIORSCR(ll,betas,qs,zs,deltavect,alphas,sigma_zs,Ns,psis,mu0,sigma2_mu0,a0_delta,*a0alpha,*b0alpha,*A,*a0psi,*b0psi,supN,dimp,*mod_h,datatype,*updatedelta,deltatype);
  if(!R_FINITE(ll)) {
    Rprintf("Fatal error in chain %d: initial likelihood is '%f'. \n",*ichain,ll);
    *iter = g;
    return;
  } else if(!R_FINITE(posterior[0])) {
    Rprintf("Fatal error in chain %d: initial posterior is '%f'. \n",*ichain,posterior[0]);
    *iter = g;
    return;
  }
  double nl;
  
  Rprintf("Chain %d initialized \n",*ichain);     
  
  /* Begin Markov chain */  
  for (g=1; g < (niter+1); g++)  {
  
    /* Update betas  */
    for(j=0; j<dimp; j++)  {
      betastar[j]=betas[j]+rnorm(0.0,Propsd[supN+j]); 
      for(t=0; t<T; t++){
        cloglogpstar[t]=0.;
        cloglogcstar[t]=0.;
        for(k=0; k<dimp; k++){
          cloglogpstar[t]+=DMp[t*dimp+k]*betastar[k];
          cloglogcstar[t]+=DMc[t*dimp+k]*betastar[k];
        }
        for(i=0; i<supN; i++){
          propp[i*T+t]=invcloglog(cloglogpstar[t]+zs[i]);
          propc[i*T+t]=invcloglog(cloglogcstar[t]+zs[i]);
        }
      }
      proppstar=GETPSTARSCR(dist2, cloglogpstar, sigma2_scrs, T, K, ncells, *dexp);
      np=LIKESCR(propp,propc,qs,delta_1s,delta_2s,alphas,Allhists,Hs,T,supN,C,Ns,proppstar);
      if(runif(0.0,1.0)<exp(np+dnorm(betastar[j],mu0[j],sqrt(sigma2_mu0[j]),1)-ll-dnorm(betas[j],mu0[j],sqrt(sigma2_mu0[j]),1))){
        betas[j]=betastar[j];
        for(t=0; t<T; t++){
          for(i=0; i<supN; i++){
            p[i*T+t]=propp[i*T+t];
            c[i*T+t]=propc[i*T+t];
          }
          cloglogp[t]=cloglogpstar[t];
          cloglogc[t]=cloglogcstar[t];
        }
        pstar=proppstar;
        ll=np;
        accept[supN+j]+=1;
      } else {
        betastar[j]=betas[j];
      }
    }
    
    if(*mod_h){
    /* update z_i  (cloglog scale)  */
      for(i=0; i<supN; i++)  {
        op=0.0;
        np=0.0;
        zstar=zs[i]+rnorm(0.0,Propsd[i]);
        temp = C[Hs[i]];
        for(t=0; t<T; t++) {
          propp[i*T+t]=invcloglog(cloglogp[t]+zstar);
          propc[i*T+t]=invcloglog(cloglogc[t]+zstar);
          if(qs[i]){
            indhist = Allhists[Hs[i] * T + t];
            if(t<temp){
              op += (indhist==0)*log(1.0-p[i*T+t])     + (indhist>=1)*log(p[i*T+t]);
              np += (indhist==0)*log(1.0-propp[i*T+t]) + (indhist>=1)*log(propp[i*T+t]);
            } else {
              op += (indhist==0)*log(1.0-c[i*T+t])     + (indhist>=1)*log(c[i*T+t]);   
              np += (indhist==0)*log(1.0-propc[i*T+t]) + (indhist>=1)*log(propc[i*T+t]);             
            }
          }
        }
        op+=dnorm(zs[i],0.0,sigma_zs,1);
        np+=dnorm(zstar,0.0,sigma_zs,1); 
        if(runif(0.0,1.0)<exp(np-op)){
          zs[i]=zstar;
          for(t=0; t<T; t++){
            p[i*T+t]=propp[i*T+t];
            c[i*T+t]=propc[i*T+t];
          }
          accept[i]+=1;
        }
      }  
      sha=1.0/Propsd[supN+dimp];
      sca=sigma_zs/sha;
      sigma_zstar = rgamma(sha,sca);
      if(sigma_zstar>sqrttol){
        op=dgamma(sigma_zstar,sha,sca,1) - pgamma(sqrttol,sha,sca,0,1);
        sca=sigma_zstar/sha;
        np=dgamma(sigma_zs,sha,sca,1) - pgamma(sqrttol,sha,sca,0,1);
        op+=log(2.0*dcauchy(sigma_zs,0.0,*A,0));
        np+=log(2.0*dcauchy(sigma_zstar,0.0,*A,0));
        for(i=0; i<supN; i++){
          op+=dnorm(zs[i],0.0,sigma_zs,1);
          np+=dnorm(zs[i],0.0,sigma_zstar,1);
        }
        proppstar=GETPSTARSCR(dist2, cloglogp, sigma2_scrstar, T, K, ncells, *dexp);
        op+=dbinom((double) ns,(double) Ns,pstar,1) - ns * log(pstar);
        np+=dbinom((double) ns,(double) Ns,proppstar,1) - ns * log(proppstar);
        if(runif(0.0,1.0)<exp(np-op)){
          sigma_zs=sigma_zstar;
          sigma2_zs=sigma_zstar*sigma_zstar;
          pstar=proppstar;
          accept[supN+dimp]+=1;
        }
      }
    }
    
    if(*updatedelta){
      /* update alpha */
      if(datatype){
        sha = *a0alpha+FREQSUM(xs,Allhists,T,J,4);
        sca = *b0alpha+FREQSUM(xs,Allhists,T,J,3);
        alphas = rbeta(sha,sca);
      }
      /* update delta_1 and delta_2 */
      if(deltatype){
        GETDELTA(deltavect, xs, Allhists, T, J, 3, a0_delta); 
        delta_1s=deltavect[0];
        delta_2s=deltavect[1];   
      } else {
        sha = a0_delta[0] + FREQSUM(xs,Allhists,T,J,1) + FREQSUM(xs,Allhists,T,J,2);
        sca = a0_delta[1] + FREQSUM(xs,Allhists,T,J,3) + FREQSUM(xs,Allhists,T,J,4);
        delta_1s = rbeta(sha,sca) / 2.0;
        delta_2s = delta_1s;
      }
      /* Update psi */
      sha = (double) *a0psi + ns;
      sca= (double) *b0psi + supN - ns;
      psis = rbeta(sha,sca);
    }
    
    ll=LIKESCR(p,c,qs,delta_1s,delta_2s,alphas,Allhists,Hs,T,supN,C,Ns,pstar);
  
    /* update x and H (latent history frequencies) */
    op=0.0;
    np=0.0;
    nbasis = (*ncolBasis ? GETCK(*numbasis,0) : 0);
    obasesum=0;
    ind=0;
    for(k=0; k< *ncolBasis; k++){
      if(min(xs[indBasis[k*3+2]]-knownxs[indBasis[k*3+2]],xs[0])+min(xs[indBasis[k*3]]-knownxs[indBasis[k*3]],xs[indBasis[k*3+1]]-knownxs[indBasis[k*3+1]])){
        indbase[k]=1.0;
        obasesum+=1;
        ind=1;
      } else {
        indbase[k]=0.0;
      }
    }
    for(j=0; j< nbasis; j++){

      if(obasesum){
        base=sample(*ncolBasis, indbase);
        xind = min(xnew[indBasis[base*3+2]] - knownxs[indBasis[base*3+2]],xnew[0]);
        c_k=GETCK(xind+min(xnew[indBasis[base*3]]-knownxs[indBasis[base*3]],xnew[indBasis[base*3+1]]-knownxs[indBasis[base*3+1]]),xind) - xind;
        if(c_k>0){
          dimadd=c_k;
          dimrem=2*c_k;
        } else {
          dimadd=2*(-c_k);
          dimrem= -c_k;
        }
        double nprop[dimadd+dimrem];
        double oprop[dimadd+dimrem];
        
        PROPFREQ(base,c_k,Hnew,indBasis,J,xnew,supN,T,p,c,C,delta_1s,delta_2s,alphas,Allhists,nprop,oprop);
        
        for(i=0; i<(dimrem+dimadd); i++){
          if(!R_FINITE(nprop[i])) {Rprintf("PROPFREQ nprop %f fatal error in chain %d: please report to <brett.mcclintock@noaa.gov> \n",nprop[i],*ichain); *iter = g; return;}
          if(!R_FINITE(oprop[i])) {Rprintf("PROPFREQ oprop %f fatal error in chain %d: please report to <brett.mcclintock@noaa.gov> \n",oprop[i],*ichain); *iter = g; return;}
          op += nprop[i];
          np += oprop[i]; 
        }
        nbasesum=0;
        for(i=0; i< *ncolBasis; i++){
          if(min(xnew[indBasis[i*3+2]]-knownxs[indBasis[i*3+2]],xnew[0])+min(xnew[indBasis[i*3]]-knownxs[indBasis[i*3]],xnew[indBasis[i*3+1]]-knownxs[indBasis[i*3+1]])){
            indbase[i]=1.0;
            nbasesum+=1;
          } else {
            indbase[i]=0.0;
          }
        }
        op += -log((double) obasesum);
        np += -log((double) nbasesum);
        obasesum=nbasesum;
      }
    }

    if(ind){
      op += ll;
      
      nstar=0.;
      for(i=0; i<supN; i++){
        qnew[i] = ((Hnew[i]) ? (int) 1 : (int) 0);
        nstar += (double) qnew[i];
        op += dbinom((double) qs[i],1.0,psis,1);
        np += dbinom((double) qnew[i],1.0,psis,1);
      }
      Nstar = nstar+rnbinom((double) nstar,pstar);
  
      nl = LIKESCR(p,c,qnew,delta_1s,delta_2s,alphas,Allhists,Hnew,T,supN,C,Nstar,pstar);
      np += nl;     
      
      op += dnbinom((double) Nstar - nstar,(double) nstar,pstar,1.0) - log((double) Ns);
      np += dnbinom((double) Ns - ns,(double) ns,pstar,1.0) - log((double) Nstar);
      
      if(runif(0.0,1.0)<exp(np-op)){
        for(i=0; i<J; i++){
          xs[i]=xnew[i];
        }
        Ns=Nstar;
        ns=nstar;
        for(i=0; i<supN; i++){
          Hs[i]=Hnew[i];
          qs[i]=qnew[i];
        }
        ll=nl;
      }
    }
    for(i=0; i<J; i++){
      xnew[i]=xs[i];
    }
    for(i=0; i<supN; i++){
      Hnew[i]=Hs[i];
      qnew[i]=qs[i];
    }
    Nstar=Ns;
    nstar=ns;
    
    /* Update N */
    Ns=ns+rnbinom((double) ns,pstar);
    
    ll=LIKESCR(p,c,qs,delta_1s,delta_2s,alphas,Allhists,Hs,T,supN,C,Ns,pstar);   
    
     /* Save draws according to thin specification */ 
    if((g % th)==0)  {
            
      N[(g/th - 1)]=Ns;

      for(j=0; j<(dimp); j++){
        beta[(g/th - 1)*(dimp)+j]=betas[j];
      }
      sigma2_z[(g/th - 1)]=sigma2_zs;
      delta_1[(g/th - 1)]=delta_1s;
      delta_2[(g/th - 1)]=delta_2s;
      alpha[(g/th - 1)]=alphas;
      psi[(g/th - 1)]=psis;
      
      if(*zind){
        for(i=0; i<supN; i++){
          z[(g/th - 1)*supN+i]=zs[i];
        }
      } else {
        for(i=0; i<supN; i++){
          z[i]=zs[i];
        }        
      }

      if(*Hind){
        for(i=0; i<supN; i++){
          H[(g/th - 1)*supN+i]=Hs[i];
        }
      } else {
        for(i=0; i<supN; i++){
          H[i]=Hs[i];
        }
      }
  
      for(j=0; j<J; j++){
        x[j]=xs[j];
      }
      
      posterior[(g/th - 1)]=POSTERIORSCR(ll,betas,qs,zs,deltavect,alphas,sigma_zs,Ns,psis,mu0,sigma2_mu0,a0_delta,*a0alpha,*b0alpha,*A,*a0psi,*b0psi,supN,dimp,*mod_h,datatype,*updatedelta,deltatype); 
      if(!R_FINITE(posterior[(g/th - 1)])) {Rprintf("Fatal error in chain %d: please report to <brett.mcclintock@noaa.gov> \n",*ichain); *iter = g; return;}
      
    }

    if(!((g+1) % n)){    
      if((g+1)<ada){
        if(*mod_h){
          for(i=0; i<supN+dimp+1; i++){
            arate[i]= accept[i]/n;
            Propsd[i]=Propsd[i]*(arate[i] < *taccept) * *tuneadjust+Propsd[i]*(arate[i] >= *taccept) / *tuneadjust;
            accept[i]=0;
          }
        } else {
          for(i=supN; i<supN+dimp; i++){
            arate[i]= accept[i]/n;
            Propsd[i]=Propsd[i]*(arate[i] < *taccept)* *tuneadjust+Propsd[i]*(arate[i] >= *taccept) / *tuneadjust;
            accept[i]=0;
          }          
        }
      } 
    }
    if(!(g%(niter/ min(niter,100)))){
      if(*printlog){
        Rprintf("Chain %d is %.0f%% complete\n",*ichain,(double) 100*g/niter);
      } else {
        Rprintf("\rChain %d is %.0f%% complete",*ichain,(double) 100*g/niter);
      }
    }
  }
  /* End Markov chain */
  if(! *printlog) Rprintf("\n");
  PutRNGstate(); 
}     
/* End function MCMCloop */

double invcloglog(double x)
{
  double invcloglog = fmin(1.-tol,fmax(tol,1.0-exp(-exp(x))));
  return(invcloglog);
}

double LIKESCR(double *p, double *c, int *qs, double delta_1, double delta_2, double alpha, int *Allhists, int *Hs, int T, int K, int supN, int *C, double Ns, double pstar)
{
  int i,k,t;
  double logdens=0.;
  int indhist;
  int temp;
  double n=0.;
  for(i=0; i<supN; i++)  {
    if(qs[i]){
      n+=1.;
      for(k=0; k<K; k++){
        temp = C[Hs[i]*K+k];
        for(t=0; t<temp; t++) {
          indhist = Allhists[Hs[i]*T*K+k*T+t];
          logdens += log( (indhist==0) * (1.0-p[i*T*K+k*T+t])
                        + (indhist==1) * p[i*T*K+k*T+t] * delta_1  
                        + (indhist==2) * p[i*T*K+k*T+t] * delta_2
                        + (indhist==3) * p[i*T*K+k*T+t] * (1.-delta_1-delta_2) * (1.-alpha)
                        + (indhist==4) * p[i*T*K+k*T+t] * (1.-delta_1-delta_2) * alpha );
          //Rprintf("i %d k %d t %d H %d indhist %d p %f delta_1 %f delta_2 %f alpha %f logdens %f \n",i,k,t,Hs[i],indhist,p[i*T*K+k*T+t],delta_1,delta_2,alpha,logdens);
        }
        for(t=temp; t<T; t++) {
          indhist = Allhists[Hs[i]*T*K+k*T+t];      
          logdens += log( (indhist==0) * (1.0-c[i*T*K+k*T+t])
                        + (indhist==1) * c[i*T*K+k*T+t] * delta_1  
                        + (indhist==2) * c[i*T*K+k*T+t] * delta_2
                        + (indhist==3) * c[i*T*K+k*T+t] * (1.-delta_1-delta_2) * (1.-alpha)
                        + (indhist==4) * c[i*T*K+k*T+t] * (1.-delta_1-delta_2) * alpha );
          //Rprintf("i %d k %d t %d H %d indhist %d c %f delta_1 %f delta_2 %f alpha %f logdens %f \n",i,k,t,Hs[i],indhist,c[i*T*K+k*T+t],delta_1,delta_2,alpha,logdens);
        } 
      }
      //Rprintf("i %d H[i] %d logdens %f \n",i,Hs[i],logdens);
    }
  }
  logdens += dbinom(n,Ns,pstar,1) - n * log(pstar);
  Rprintf("pstar %f logdens %f \n",pstar,logdens);
  return(logdens);   
}

double POSTERIORSCR(double ll, double *beta, int *qs, double *z, double *deltavect, double alpha, double sigma_z, double Ns, double psi, double *mu0, double *sigma2_mu0, double *a0_delta, double a0_alpha, double b0_alpha, double A, double a0psi, double b0psi, int supN, int pdim, int modh, int datatype, int updatedelta, int deltatype)
{
  double pos=ll;
  int i,j;
  for(j=0; j<pdim; j++){
    pos += dnorm(beta[j],mu0[j],sqrt(sigma2_mu0[j]),1);
  }
  if(modh){
    for(i=0; i<supN; i++){
      pos += dnorm(z[i],0.0,sigma_z,1);
    }
    pos += log(2.0*dcauchy(sigma_z,0.0,A,0));
  }
  if(updatedelta){
    if(deltatype){
      pos += DDIRICHLET(deltavect,a0_delta,3);
    } else {
      pos += dbeta((deltavect[0]+deltavect[1]),a0_delta[0],a0_delta[1],1);
    }
    if(datatype){
      pos += dbeta(alpha,a0_alpha,b0_alpha,1);
    }
    for(i=0; i<supN; i++){
      pos += dbinom((double) qs[i],1.0,psi,1);
    }
    pos += dbeta(psi,a0psi,b0psi,1);
  }
  pos += -log(Ns);
  return(pos);
}

double GETPSTARSCR(double *dist2, double *cloglogp, double sigma2_scr, int T, int K, int ncells, double dexp)
{
  int i,t, k;
  double oneminuspstar[ncells], esa=0.;
  for(i=0; i<ncells; i++){
    oneminuspstar[i]=1.;
    for(k=0; k<K; k++){
      for(t=0; t<T; t++){
        oneminuspstar[i] *= 1. - invcloglog(cloglogp[k*T+t]-1.0/(2.0*sigma2_scr)*pow(dist2[k*ncells+i],dexp));
      }
    }
    esa += 1. - fmax(oneminuspstar[i],tol);
    Rprintf("cell %d oneminuspstar %f esa %f \n",i,oneminuspstar[i],esa);
  }
  double pstar = esa / ncells;
  Rprintf("pstar %f esa %f dexp %f ncells %d sigma2_scr %f \n",pstar,esa,dexp,ncells,sigma2_scr);
  return(pstar);
}