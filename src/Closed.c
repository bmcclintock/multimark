#include "ranlib.h"
#include "func.h"
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

/* FUNCTION DEFINITIONS */

double EXPIT(double x)
{
  double expit = fmin(1.-tol,fmax(tol,1.0/(1.0+exp(-x))));
  return(expit);
}

double LIKE(double *p, double *c, int *qs, double delta_1, double delta_2, double alpha, int *Allhists, int *Hs, int T, int supN, int *C, double Ns, double pstar)
{
  int i,t;
  double logdens=0.;
  int indhist;
  int temp;
  double n=0.;
  for(i=0; i<supN; i++)  {
    if(qs[i]){
      n+=1.;
      temp = C[Hs[i]];
      for(t=0; t<temp; t++) {
        indhist = Allhists[Hs[i] * T + t];
        logdens += log( (indhist==0) * (1.0-p[i*T+t])
                          + (indhist==1) * p[i*T+t] * delta_1  
                          + (indhist==2) * p[i*T+t] * delta_2
                          + (indhist==3) * p[i*T+t] * (1.-delta_1-delta_2) * (1.-alpha)
                          + (indhist==4) * p[i*T+t] * (1.-delta_1-delta_2) * alpha );
      }
      for(t=temp; t<T; t++) {
        indhist = Allhists[Hs[i] * T + t];      
        logdens += log( (indhist==0) * (1.0-c[i*T+t])
                          + (indhist==1) * c[i*T+t] * delta_1  
                          + (indhist==2) * c[i*T+t] * delta_2
                          + (indhist==3) * c[i*T+t] * (1.-delta_1-delta_2) * (1.-alpha)
                          + (indhist==4) * c[i*T+t] * (1.-delta_1-delta_2) * alpha );
      } 
    }
  }
  logdens += dbinom(n,Ns,pstar,1) - n * log(pstar);
  return(logdens);   
}

double POSTERIOR(double ll, double *beta, int *qs, double *z, double *deltavect, double alpha, double sigma_z, double Ns, double psi, double *mu0, double *sigma2_mu0, double *a0_delta, double a0_alpha, double b0_alpha, double A, double a0psi, double b0psi, int supN, int pdim, int modh, int datatype, int updatedelta, int deltatype)
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

double GETprodh(int *Allhists, double *p, double *c, int *C, double delta_1, double delta_2, double alpha, int j,int T, int i)
{
  int t;
  double logdens=0.0;
  int indhist;
  
  int temp = C[j];
  for(t=0; t<temp; t++) {
    indhist = Allhists[j * T + t];
    logdens += log( (indhist==0) * (1.0-p[i*T+t]) 
                      + (indhist==1) * p[i*T+t] * delta_1 
                      + (indhist==2) * p[i*T+t] * delta_2
                      + (indhist==3) * p[i*T+t] * (1.-delta_1-delta_2) * (1.-alpha)
                      + (indhist==4) * p[i*T+t] * (1.-delta_1-delta_2) * alpha );  
  }
  for(t=temp; t<T; t++) {
    indhist = Allhists[j * T + t];      
    logdens += log( (indhist==0) * (1.0-c[i*T+t]) 
                      + (indhist==1) * c[i*T+t] * delta_1 
                      + (indhist==2) * c[i*T+t] * delta_2
                      + (indhist==3) * c[i*T+t] * (1.-delta_1-delta_2) * (1.-alpha)
                      + (indhist==4) * c[i*T+t] * (1.-delta_1-delta_2) * alpha );  
  }
  double dens = exp(logdens);
  if(dens<tol) dens = tol;
  return(dens);
}

void PROPFREQ(int icol,int c_k,int *Hnew, int *indBasis, int J, int *xnew, int supN, int T, double *p, double *c, int *C, double delta_1, double delta_2, double alpha, int *Allhists, double *nprop, double *oprop)
{  
  int remove_xi[3];
  int add_xi[3];
  int j,i,k;
  int absc_k=abs(c_k);
  int remove[absc_k], add[absc_k];
  double prodz[supN], prodh[supN];
  double prodzsum, prodhsum;  
  double temp = runif(0.0,1.0);
  if(c_k > 0){
    remove_xi[0]=1;
    remove_xi[1]=1;
    remove_xi[2]=0;
    add_xi[0]=0;
    add_xi[1]=0;
    add_xi[2]=1;
    if(xnew[0]<(c_k)) temp = 0.0;
  } else if(c_k < 0){
    add_xi[0]=1;
    add_xi[1]=1;
    add_xi[2]=0;
    remove_xi[0]=0;
    remove_xi[1]=0;
    remove_xi[2]=1;
    if(xnew[0]<(-2*c_k)) temp = 0.0;
  } 
  int count=0;
  if(temp<0.5) {
    goto S10;
  } else {
    goto S20;
  }
  S10:
    for(j=0; j<3; j++){
      prodzsum=0.0;
      prodhsum=0.0;
      if(remove_xi[j]){
        for(i=0; i<supN; i++) {
          prodz[i] = -1.0;
          prodh[i] = -1.0;
          if(Hnew[i]==indBasis[icol*3+j]){
            prodz[i] = 1. - GETprodh(Allhists,p,c,C,delta_1,delta_2,alpha,indBasis[icol*3+j],T,i);
            prodzsum+=prodz[i];
          } else if(!Hnew[i]){
            prodh[i] = GETprodh(Allhists,p,c,C,delta_1,delta_2,alpha,indBasis[icol*3+j],T,i);
            prodhsum+=prodh[i];
          }          
        }
        ProbSampleNoReplace(supN, prodz, absc_k, remove); 
        for(k=0; k<absc_k; k++){
          Hnew[remove[k]]=0;
          prodh[remove[k]] = 1. - prodz[remove[k]];
          prodhsum+=prodh[remove[k]];    
        }
        for(k=0; k<absc_k; k++){
          nprop[count] = log(prodz[remove[k]])-log(prodzsum);
          oprop[count] = log(prodh[remove[k]])-log(prodhsum);
          prodzsum -=  prodz[remove[k]];
          prodhsum -=  prodh[remove[k]];
          count+=1;
        }
      }
    }
    if(temp<0.5) goto S20;
    else goto S30;
    S20:
      for(j=0; j<3; j++){
        prodzsum=0.0;
        prodhsum=0.0;
        if(add_xi[j]){
          for(i=0; i<supN; i++) {
            prodz[i] = -1.0;
            prodh[i] = -1.0;
            if(!Hnew[i]){
              prodh[i] = GETprodh(Allhists,p,c,C,delta_1,delta_2,alpha,indBasis[icol*3+j],T,i);
              prodhsum+=prodh[i];
            } else if(Hnew[i]==indBasis[icol*3+j]){
              prodz[i] = 1. - GETprodh(Allhists,p,c,C,delta_1,delta_2,alpha,indBasis[icol*3+j],T,i);
              prodzsum+=prodz[i];
            }
          }
          ProbSampleNoReplace(supN, prodh, absc_k, add);
          for(k=0; k<absc_k; k++){
            Hnew[add[k]]=indBasis[icol*3+j];
            prodz[add[k]] = 1. - prodh[add[k]];
            prodzsum+=prodz[add[k]];
          }
          for(k=0; k<absc_k; k++){
            nprop[count] = log(prodh[add[k]])-log(prodhsum);  
            oprop[count] = log(prodz[add[k]])-log(prodzsum); 
            prodzsum -=  prodz[add[k]];
            prodhsum -=  prodh[add[k]];
            count+=1;
          }
        }
      }
      if(temp>=0.5) goto S10;
      else goto S30;
      S30:
        xnew[indBasis[icol*3]]-=c_k;
      xnew[indBasis[icol*3+1]]-=c_k;
      xnew[indBasis[icol*3+2]]+=c_k;  
      xnew[0]+=c_k;
}

double GETPSTAR(int npts, double *weight, double *node, double *logitp, double sigma_zs, int dimp, int T)
{
  int i,t;
  double temp[npts];
  double oneminuspstar=0.;
  for(i=0; i<npts; i++){
    temp[i]=1.;
    for(t=0; t<T; t++){
      temp[i]*=(1. - 1. / (1. + exp(-(sqrt(2.0)*sigma_zs*node[i]+logitp[t]))));
    }
    oneminuspstar += 1.0/sqrt(3.14159265359)*weight[i]*temp[i]; 
  }
  double pstar = 1.-fmin(1.-tol,fmax(tol,oneminuspstar));
  return(pstar);
}

// Define function ClosedC to draw samples from the posterior distribution

void ClosedC(int *ichain, double *mu0, double *sigma2_mu0, double *beta, double *z, double *sigma2_z, double *delta_1, double *delta_2, double *alpha, int *x, double *N, double *psi, int *H, 
              int *noccas, int *M, double *a0delta, double *a0alpha, double *b0alpha, double *A, double *a0psi, double *b0psi,
              double *Propsd, double *accept, double *posterior,
              int *nHists, int *Allhists, int *C, int *indBasis, int *ncolBasis, int *knownx, double *DMp, double *DMc, int *pdim,
              int *iter, int *thin, int *adapt, int *bin, double *taccept, double *tuneadjust, int *numbasis,
              int *npoints, double *weights, double *nodes, int *mod_h, int *data_type, int *zind, int *Hind, int *updatedelta, int *delta_type, int *printlog)
{
  
  GetRNGstate(); 

  int T = *noccas;
  int supN = *M; 
  int Mk = supN*(T);      
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
  
  double betas[dimp], zs[supN], zstar, sigma2_zs, sigma_zstar, alphas, delta_1s, delta_2s, psis;
  double betastar[dimp];
  int xs[J], xnew[J], knownxs[J], Hs[supN], Hnew[supN];
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
  
  double logitp[T], logitpstar[T], logitc[T], logitcstar[T];
  double propp[Mk], propc[Mk];
  for(t=0; t<T; t++){
    logitp[t]=0.;
    logitc[t]=0.;
    for(k=0; k<dimp; k++){
      logitp[t]+=DMp[t*dimp+k]*betas[k];  
      logitc[t]+=DMc[t*dimp+k]*betas[k];  
    }
    logitpstar[t]=logitp[t];
    logitcstar[t]=logitc[t];
    for(i=0; i<supN; i++){
      p[i*T+t]=EXPIT((logitp[t]+zs[i]));
      c[i*T+t]=EXPIT((logitc[t]+zs[i]));
      propp[i*T+t]=p[i*T+t];
      propc[i*T+t]=c[i*T+t];
    }
  }
  double pstar=GETPSTAR(npts, weight, node, logitp, sigma_zs, dimp, T);
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
  
  double indbase[*ncolBasis ? *ncolBasis : 1];
  for(j=0; j< *ncolBasis; j++){
    indbase[j]=0.0;
  }
  int c_k, temp;
  int ind, obasesum, nbasesum;
    
  /* Calculate the log-likelihood */  
  double ll=LIKE(p,c,qs,delta_1s,delta_2s,alphas,Allhists,Hs,T,supN,C,Ns,pstar);
  posterior[0]=POSTERIOR(ll,betas,qs,zs,deltavect,alphas,sigma_zs,Ns,psis,mu0,sigma2_mu0,a0_delta,*a0alpha,*b0alpha,*A,*a0psi,*b0psi,supN,dimp,*mod_h,datatype,*updatedelta,deltatype);
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
        logitpstar[t]=0.;
        logitcstar[t]=0.;
        for(k=0; k<dimp; k++){
          logitpstar[t]+=DMp[t*dimp+k]*betastar[k];
          logitcstar[t]+=DMc[t*dimp+k]*betastar[k];
        }
        for(i=0; i<supN; i++){
          propp[i*T+t]=EXPIT(logitpstar[t]+zs[i]);
          propc[i*T+t]=EXPIT(logitcstar[t]+zs[i]);
        }
      }
      proppstar=GETPSTAR(npts, weight, node, logitpstar, sigma_zs, dimp, T);
      np=LIKE(propp,propc,qs,delta_1s,delta_2s,alphas,Allhists,Hs,T,supN,C,Ns,proppstar);
      if(runif(0.0,1.0)<exp(np+dnorm(betastar[j],mu0[j],sqrt(sigma2_mu0[j]),1)-ll-dnorm(betas[j],mu0[j],sqrt(sigma2_mu0[j]),1))){
        betas[j]=betastar[j];
        for(t=0; t<T; t++){
          for(i=0; i<supN; i++){
            p[i*T+t]=propp[i*T+t];
            c[i*T+t]=propc[i*T+t];
          }
          logitp[t]=logitpstar[t];
          logitc[t]=logitcstar[t];
        }
        pstar=proppstar;
        ll=np;
        accept[supN+j]+=1;
      } else {
        betastar[j]=betas[j];
      }
    }
    
    if(*mod_h){
    /* update z_i  (logit scale)  */
      for(i=0; i<supN; i++)  {
        op=0.0;
        np=0.0;
        zstar=zs[i]+rnorm(0.0,Propsd[i]);
        temp = C[Hs[i]];
        for(t=0; t<T; t++) {
          propp[i*T+t]=EXPIT(logitp[t]+zstar);
          propc[i*T+t]=EXPIT(logitc[t]+zstar);
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
        proppstar=GETPSTAR(npts, weight, node, logitp, sigma_zstar, dimp, T);
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
    
    ll=LIKE(p,c,qs,delta_1s,delta_2s,alphas,Allhists,Hs,T,supN,C,Ns,pstar);
  
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
  
      nl = LIKE(p,c,qnew,delta_1s,delta_2s,alphas,Allhists,Hnew,T,supN,C,Nstar,pstar);
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
    
    ll=LIKE(p,c,qs,delta_1s,delta_2s,alphas,Allhists,Hs,T,supN,C,Ns,pstar);   
    
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
      
      posterior[(g/th - 1)]=POSTERIOR(ll,betas,qs,zs,deltavect,alphas,sigma_zs,Ns,psis,mu0,sigma2_mu0,a0_delta,*a0alpha,*b0alpha,*A,*a0psi,*b0psi,supN,dimp,*mod_h,datatype,*updatedelta,deltatype); 
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
