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

double invcloglog(double x)
{
  double invcloglog = fmin(1.-tol,fmax(tol,1.0-exp(-exp(x))));
  return(invcloglog);
}

double getdist(double x1, double y1, double x2, double y2)
{
  double dist = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
  return(dist);
}

double GETcellprob(int indhist, double prob, double delta_1, double delta_2, double da3, double da4)
{
  double cellprob=0.;
  switch(indhist)
  {
  case 0 :
    cellprob = log(1.0-prob);
    break;
  case 1 :
    cellprob = log(prob * delta_1);
    break;
  case 2 :
    cellprob = log(prob * delta_2);
    break;
  case 3 :
    cellprob = log(prob * da3);
    break;
  case 4 :
    cellprob = log(prob * da4);
    break;
  }
  return(cellprob);
}

double GETprodhSCR(int *Allhists, double *p, double *c, int *C, double delta_1, double delta_2, double alpha, int j,int T, int K, int i)
{
  int t, k;
  double logdens=0.0;
  //int indhist;
  int firstcap;
  double da3 = (1.-delta_1-delta_2) * (1.-alpha);
  double da4 = (1.-delta_1-delta_2) * alpha; 
  
  for(k=0; k<K; k++){
    firstcap = C[j*K+k];
    for(t=0; t<firstcap; t++) {
      logdens += GETcellprob(Allhists[j*T*K+k*T+t], p[i*T*K+k*T+t], delta_1, delta_2, da3, da4);
      //indhist = Allhists[j*T*K+k*T+t];
      //logdens += log( (indhist==0) * (1.0-p[i*T*K+k*T+t]) 
      //                  + (indhist==1) * p[i*T*K+k*T+t] * delta_1 
      //                  + (indhist==2) * p[i*T*K+k*T+t] * delta_2
      //                  + (indhist==3) * p[i*T*K+k*T+t] * da3
      //                  + (indhist==4) * p[i*T*K+k*T+t] * da4);  
    }
    for(t=firstcap; t<T; t++) {
      logdens += GETcellprob(Allhists[j*T*K+k*T+t], c[i*T*K+k*T+t], delta_1, delta_2, da3, da4);
      //indhist = Allhists[j*T*K+k*T+t];   
      //logdens += log( (indhist==0) * (1.0-c[i*T*K+k*T+t]) 
      //                  + (indhist==1) * c[i*T*K+k*T+t] * delta_1 
      //                  + (indhist==2) * c[i*T*K+k*T+t] * delta_2
      //                  + (indhist==3) * c[i*T*K+k*T+t] * da3
      //                  + (indhist==4) * c[i*T*K+k*T+t] * da4);  
    }
  }
  double dens = exp(logdens);
  if(dens<tol) dens = tol;
  return(dens);
}

double LIKESCR(double *p, double *c, int *qs, double delta_1, double delta_2, double alpha, int *Allhists, int *Hs, int T, int K, int supN, int *C, double Ns, double pstar)
{
  int i, k, t;
  double logdens=0.;
  //int indhist;
  int temp;
  double n=0.;
  double da3 = (1.-delta_1-delta_2) * (1.-alpha);
  double da4 = (1.-delta_1-delta_2) * alpha; 
  
  for(i=0; i<supN; i++)  {
    if(qs[i]){
      n+=1.;
      for(k=0; k<K; k++){
        temp = C[Hs[i]*K+k];
        for(t=0; t<temp; t++) {
          logdens += GETcellprob(Allhists[Hs[i]*T*K+k*T+t], p[i*T*K+k*T+t], delta_1, delta_2, da3, da4);
        }
        for(t=temp; t<T; t++) {
          logdens += GETcellprob(Allhists[Hs[i]*T*K+k*T+t], c[i*T*K+k*T+t], delta_1, delta_2, da3, da4);
        }
      }
    }
  }
  logdens += dbinom(n,Ns,pstar,1) - n * log(pstar);
  //Rprintf("pstar %f logdens %f \n",pstar,logdens);
  return(logdens);   
}

double POSTERIORSCR(double ll, double *beta, int *qs, double *deltavect, double alpha, double sigma2_scr, double Ns, double psi, double *mu0, double *sigma2_mu0, double *a0_delta, double a0_alpha, double b0_alpha, double *sigma_bounds, double a0psi, double b0psi, int supN, int pdim, int datatype, int updatedelta, int deltatype, double Area)
{
  double pos=ll;
  int i,j;
  for(j=0; j<pdim; j++){
    pos += dnorm(beta[j],mu0[j],sqrt(sigma2_mu0[j]),1);
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
  pos += dunif(sqrt(sigma2_scr),sigma_bounds[0],sigma_bounds[1],1);//log(2.0*dcauchy(sqrt(sigma2_scr),0.0,A,0));
  pos += supN * log(1./Area);
  pos += -log(Ns);
  return(pos);
}

void PROPFREQSCR(int icol,int c_k,int *Hnew, int *indBasis, int J, int *xnew, int supN, int T, int K, double *p, double *c, int *C, double delta_1, double delta_2, double alpha, int *Allhists, double *nprop, double *oprop)
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
            prodz[i] = 1. - GETprodhSCR(Allhists,p,c,C,delta_1,delta_2,alpha,indBasis[icol*3+j],T,K,i);
            prodzsum+=prodz[i];
          } else if(!Hnew[i]){
            prodh[i] = GETprodhSCR(Allhists,p,c,C,delta_1,delta_2,alpha,indBasis[icol*3+j],T,K,i);
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
              prodh[i] = GETprodhSCR(Allhists,p,c,C,delta_1,delta_2,alpha,indBasis[icol*3+j],T,K,i);
              prodhsum+=prodh[i];
            } else if(Hnew[i]==indBasis[icol*3+j]){
              prodz[i] = 1. - GETprodhSCR(Allhists,p,c,C,delta_1,delta_2,alpha,indBasis[icol*3+j],T,K,i);
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

void GETPC(double *p, double *c, double *cloglogp, double *cloglogc, double *beta, double sigma2, double *DMp, double *DMc, double *dist2centers, int dimp, int supN, int T, int K, int *msk, int *cummind, int *mind, double dexp)
{
  int k, t, j, i;
  double tmp = 1.0/(dexp * sigma2);
  //double cloglogc[T*K];
  for(k=0; k<K; k++){
    for(t=cummind[k]; t<cummind[k+1]; t++){
      //for(t=0; t<T; t++){
      cloglogp[k*T+mind[t]]=0.;
      cloglogc[k*T+mind[t]]=0.;
      for(j=0; j<dimp; j++){
        cloglogp[k*T+mind[t]]+=DMp[(k*T+mind[t])*dimp+j]*beta[j];  
        cloglogc[k*T+mind[t]]+=DMc[(k*T+mind[t])*dimp+j]*beta[j];  
      }
      for(i=0; i<supN; i++){
        p[i*T*K+k*T+mind[t]]=invcloglog((cloglogp[k*T+mind[t]]-tmp*dist2centers[k*supN+i]));
        c[i*T*K+k*T+mind[t]]=invcloglog((cloglogc[k*T+mind[t]]-tmp*dist2centers[k*supN+i]));
        //p[i*T*K+k*T+t]=exp(-tmp*dist2centers[k*supN+i]) * msk[k*T+t];
        //c[i*T*K+k*T+t]=exp(-tmp*dist2centers[k*supN+i]) * msk[k*T+t];
      }
      //}
    }
  }    
}

double GETPSTARSCR(double *dist2, double *cloglogp, double sigma2_scr, int T, int K, int ncells, int *msk, int *cummind, int *mind, double dexp)
{
  int i,t, k;
  double tmp = 1.0/(dexp * sigma2_scr);
  double oneminuspstar[ncells], esa=0., clog;
  double ncell = (double) ncells;
  for(i=0; i<ncells; i++){
    oneminuspstar[i]=1.;
  }
  for(k=0; k<K; k++){
    for(t=cummind[k]; t<cummind[k+1]; t++){
      clog = cloglogp[k*T+mind[t]];
      for(i=0; i<ncells; i++){
        oneminuspstar[i] *= (1. - invcloglog(clog - tmp * dist2[k*ncells+i]));
      }
    }
  }
  for(i=0; i<ncells; i++){
    esa += 1. - fmax(oneminuspstar[i],tol);
  }
  double pstar = esa / ncell;
  //Rprintf("pstar %f esa %f dexp %f ncells %d sigma2_scr %f \n",pstar,esa,dexp,ncells,sigma2_scr);
  return(pstar);
}

// Define function ClosedC to draw samples from the posterior distribution

void ClosedSCRC(int *ichain, double *mu0, double *sigma2_mu0, double *beta, double *sigma2_scr, double *delta_1, double *delta_2, double *alpha, int *x, double *N, double *psi, int *H, int *center,
                int *ntraps, int *noccas, int *M, double *a0delta, double *a0alpha, double *b0alpha, double *sigma_bounds, double *a0psi, double *b0psi,
                double *Propsd, int *NNvect, int *numnn, int *cumnumnn, double *accept, double *posterior,
                int *nHists, int *Allhists, int *C, int *indBasis, int *ncolBasis, int *knownx, double *DMp, double *DMc, int *pdim,
                int *iter, int *thin, int *adapt, int *bin, double *taccept, double *tuneadjust, int *numbasis,
                int *data_type, int *Hind, int *centerind, int *updatedelta, int *delta_type, double *dexp, double *dist2, int *ncell, double *Area, int *msk, int *cummind, int *mind, int *printlog)
{
  
  GetRNGstate(); 
  
  int T = *noccas;
  int K = *ntraps;
  int supN = *M; 
  int Mk = supN*T*K;      
  int datatype = *data_type;
  int deltatype = *delta_type;
  
  int niter, th, ada, n;
  int J = *nHists;
  int firstcap;
  
  niter = *iter;              /* Number of iterations in the Markov chain */
  th = *thin;                 /* Number of iterations for thinning */
  ada = *adapt;               /* Number of iterations for pilot tuning */
  n = *bin;                   /* Number of iterations between each pilot tune */
  
  int dimp = *pdim;               /* length of beta vector */
  int ncells = *ncell;
  
  double betas[dimp], sigma2_scrs, sigma2_scrstar, sigma_scrstar, alphas, delta_1s, delta_2s, psis;
  double betastar[dimp];
  int xs[J], xnew[J], knownxs[J], Hs[supN], Hnew[supN], centers[supN], centerstar;
  int qs[supN], qnew[supN]; 
  double Ns;
  double p[Mk], c[Mk];
  double sha, sca;
  int base, nbasis;
  double sqrttol = sqrt(tol);
  
  int t, i, j, k, l, g=0;
  double op, np;
  int dimrem, dimadd;
  
  for(k=0; k<dimp; k++){
    betas[k]=beta[k];
    betastar[k]=betas[k];
  }
  
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
  
  //int NNvect[(cumnumnn[ncells-1]+numnn[ncells-1])];
  
  //double dist2[K*ncells], dist2powcells[K*ncells];
  //double dist2powcells[K*ncells];
  //for(i=0; i<ncells; i++){
  //  for(k=0; k<K; k++){
  //    //dist2[k*ncells+i] = getdist(studyArea[i],studyArea[ncells+i],trapCoords[k],trapCoords[K+k]);
  //    dist2powcells[k*ncells+i] = pow(dist2[k*ncells+i],*dexp);
  //  }
  //}
  
  double cloglogp[T*K], cloglogc[T*K];
  double propp[T*K], propc[T*K];
  //double cloglogp[T*K], cloglogpstar[T*K], cloglogc[T*K], cloglogcstar[T*K];
  //double propp[Mk], propc[Mk];
  
  double dist2centers[K*supN], diststar[K];
  for(k=0; k<K; k++){
    for(i=0; i<supN; i++){
      dist2centers[k*supN+i] = dist2[k*ncells+centers[i]];
      for(t=0; t<T; t++){
        p[i*T*K+k*T+t]=0.;
        c[i*T*K+k*T+t]=0.;
      }
    }
  }
  for(k=0; k<K; k++){
    for(t=0; t<T; t++){
      cloglogp[k*T+t]=0.;
      cloglogc[k*T+t]=0.;
      propp[k*T+t]=0.;
      propc[k*T+t]=0.;
    }
  }
  
  GETPC(p,c,cloglogp,cloglogc,betas,sigma2_scrs,DMp,DMc,dist2centers,dimp,supN,T,K,msk, cummind, mind,*dexp);
  
  double pstar=GETPSTARSCR(dist2, cloglogp, sigma2_scrs, T, K, ncells, msk, cummind, mind, *dexp);
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
  int c_k, tmp;
  int ind, obasesum, nbasesum;
  double temp;
  
  /* Calculate the log-likelihood */  
  double ll=LIKESCR(p,c,qs,delta_1s,delta_2s,alphas,Allhists,Hs,T,K,supN,C,Ns,pstar);
  posterior[0]=POSTERIORSCR(ll,betas,qs,deltavect,alphas,sigma2_scrs,Ns,psis,mu0,sigma2_mu0,a0_delta,*a0alpha,*b0alpha,sigma_bounds,*a0psi,*b0psi,supN,dimp,datatype,*updatedelta,deltatype,*Area);
  //Rprintf("ll %f posterior %f \n",ll,posterior[0]);
  if(!R_FINITE(ll)) {
    Rprintf("Fatal error in chain %d: initial likelihood is '%f'. \n",*ichain,ll);
    *iter = g;
    return;
  } else if(!R_FINITE(posterior[0])) {
    Rprintf("Fatal error in chain %d: initial posterior is '%f'. \n",*ichain,posterior[0]);
    *iter = g;
    return;
  }
  double nl, ol;
  
  Rprintf("Chain %d initialized \n",*ichain);     
  
  /* Begin Markov chain */  
  for (g=1; g < (niter+1); g++)  {
    
    /* Update betas  */
    for(l=0; l<dimp; l++)  {
      betastar[l]=betas[l]+rnorm(0.0,Propsd[l]); 
      GETPC(p,c,cloglogp,cloglogc,betastar,sigma2_scrs,DMp,DMc,dist2centers,dimp,supN,T,K,msk, cummind, mind,*dexp);
      //for(k=0; k<K; k++){
      //  for(t=0; t<T; t++){
      //    cloglogpstar[k*T+t]=0.;
      //    cloglogcstar[k*T+t]=0.;
      //    for(j=0; j<dimp; j++){
      //      cloglogpstar[k*T+t]+=DMp[(k*T+t)*dimp+j]*betastar[j];  
      //      cloglogcstar[k*T+t]+=DMc[(k*T+t)*dimp+j]*betastar[j];  
      //    }
      //    for(i=0; i<supN; i++){
      //      propp[i*T*K+k*T+t]=invcloglog((cloglogpstar[k*T+t]-1.0/(*dexp * sigma2_scrs)*dist2centers[k*supN+i]));
      //      propc[i*T*K+k*T+t]=invcloglog((cloglogcstar[k*T+t]-1.0/(*dexp * sigma2_scrs)*dist2centers[k*supN+i]));
      //    }
      //  }
      //}      
      proppstar=GETPSTARSCR(dist2, cloglogp, sigma2_scrs, T, K, ncells, msk, cummind, mind, *dexp);
      np=LIKESCR(p,c,qs,delta_1s,delta_2s,alphas,Allhists,Hs,T,K,supN,C,Ns,proppstar);
      //Rprintf("g %d l %d betastar %f beta %f np %f ll %f nprior %f oprior %f R %f \n",g,l,betastar[l],betas[l],np,ll,dnorm(betastar[l],mu0[l],sqrt(sigma2_mu0[l]),1),dnorm(betas[l],mu0[l],sqrt(sigma2_mu0[l]),1),exp(np+dnorm(betastar[l],mu0[l],sqrt(sigma2_mu0[l]),1)-ll-dnorm(betas[l],mu0[l],sqrt(sigma2_mu0[l]),1)));
      if(runif(0.0,1.0)<exp(np+dnorm(betastar[l],mu0[l],sqrt(sigma2_mu0[l]),1)-ll-dnorm(betas[l],mu0[l],sqrt(sigma2_mu0[l]),1))){
        betas[l]=betastar[l];
        //for(k=0; k<K; k++){
        //  for(t=0; t<T; t++){
        //    //for(i=0; i<supN; i++){
        //    //  p[i*T*K+k*T+t]=propp[i*T*K+k*T+t];
        //    //  c[i*T*K+k*T+t]=propc[i*T*K+k*T+t];
        //    //}
        //    cloglogp[k*T+t]=cloglogpstar[k*T+t];
        //    cloglogc[k*T+t]=cloglogcstar[k*T+t];
        //  }
        //}
        pstar=proppstar;
        ll=np;
        accept[l]+=1;
      } else {
        betastar[l]=betas[l];
      }
    }
    
    //GETPC(p,c,cloglogp,cloglogc,betas,sigma2_scrs,DMp,DMc,dist2centers,dimp,supN,T,K,msk, cummind, mind,*dexp);
    
    // Update sigma2_scr
    sha=1.0/Propsd[dimp];
    sca=sigma_scrs/sha;
    sigma_scrstar = rgamma(sha,sca);
    sigma2_scrstar = sigma_scrstar * sigma_scrstar;
    if(sigma_scrstar>sqrttol){
      op=dgamma(sigma_scrstar,sha,sca,1) - pgamma(sqrttol,sha,sca,0,1);
      sca=sigma_scrstar/sha;
      np=dgamma(sigma_scrs,sha,sca,1) - pgamma(sqrttol,sha,sca,0,1);
      op+=dunif(sigma_scrs,sigma_bounds[0],sigma_bounds[1],1);//log(2.0*dcauchy(sigma_scrs,0.0,*A,0));
      np+=dunif(sigma_scrstar,sigma_bounds[0],sigma_bounds[1],1);//log(2.0*dcauchy(sigma_scrstar,0.0,*A,0));
      
      GETPC(p,c,cloglogp,cloglogc,betas,sigma2_scrstar,DMp,DMc,dist2centers,dimp,supN,T,K,msk, cummind, mind,*dexp);      
      proppstar=GETPSTARSCR(dist2, cloglogp, sigma2_scrstar, T, K, ncells, msk, cummind, mind, *dexp);
      
      //op+=dbinom((double) ns,(double) Ns,pstar,1) - ns * log(pstar);
      //np+=dbinom((double) ns,(double) Ns,proppstar,1) - ns * log(proppstar);
      
      //for(k=0; k<K; k++){
      //  for(t=0; t<T; t++){
      //    for(i=0; i<supN; i++){
      //      propp[i*T*K+k*T+t]=invcloglog((cloglogp[k*T+t]-1.0/(*dexp * sigma2_scrstar)*dist2centers[k*supN+i]));
      //      propc[i*T*K+k*T+t]=invcloglog((cloglogc[k*T+t]-1.0/(*dexp * sigma2_scrstar)*dist2centers[k*supN+i]));
      //    }
      //  }
      //}      
      
      op+=ll;
      nl=LIKESCR(p,c,qs,delta_1s,delta_2s,alphas,Allhists,Hs,T,K,supN,C,Ns,proppstar);      
      np+=nl;
      //Rprintf("g %d sigma2_scrstar %f sigma2_scr %f nl %f ll %f np %f op %f R %f \n",g,sigma2_scrstar,sigma2_scrs,nl,ll,np,op,exp(np-op));
      
      if(runif(0.0,1.0)<exp(np-op)){
        sigma_scrs=sigma_scrstar;
        sigma2_scrs=sigma_scrs*sigma_scrs;
        pstar=proppstar;
        //for(k=0; k<K; k++){
        //  for(t=0; t<T; t++){
        //    for(i=0; i<supN; i++){
        //      p[i*T*K+k*T+t]=propp[i*T*K+k*T+t];
        //      c[i*T*K+k*T+t]=propc[i*T*K+k*T+t];
        //    }
        //  }
        //}      
        ll=nl;
        accept[dimp]+=1;
      }
    }
    
    GETPC(p,c,cloglogp,cloglogc,betas,sigma2_scrs,DMp,DMc,dist2centers,dimp,supN,T,K,msk, cummind, mind,*dexp);   
    
    //Update centers
    //int tmpl=1;
    for(i=0; i<supN; i++){
      nl=0.;
      ol=0.;
      tmp = (int) floor(runif(0.0,numnn[centers[i]]));
      centerstar = NNvect[cumnumnn[centers[i]]+tmp];
      //Rprintf("g %d i %d tmp %d centerstar %d center %d numnn[center] %d cumnumnn[center] %d \n",g,i,tmp,centerstar,centers[i],numnn[centers[i]],cumnumnn[centers[i]]);
      
      for(k=0; k<K; k++){
        diststar[k] = dist2[k*ncells+centerstar];
        temp = 1.0/(*dexp * sigma2_scrs) * diststar[k];
        if(qs[i]){
          firstcap = C[Hs[i]*K+k];
          for(t=0; t<firstcap; t++) {
            propp[k*T+t] = invcloglog(cloglogp[k*T+t]-temp) * msk[k*T+t];
            propc[k*T+t] = invcloglog(cloglogc[k*T+t]-temp) * msk[k*T+t];
            //propp[k*T+t] = exp(-temp) * msk[k*T+t];
            //propc[k*T+t] = exp(-temp) * msk[k*T+t];
            nl += GETcellprob(Allhists[Hs[i]*T*K+k*T+t], propp[k*T+t],       delta_1s, delta_2s, ((1.-delta_1s-delta_2s) * (1.-alphas)), ((1.-delta_1s-delta_2s) * alphas));
            ol += GETcellprob(Allhists[Hs[i]*T*K+k*T+t], p[i*T*K+k*T+t],     delta_1s, delta_2s, ((1.-delta_1s-delta_2s) * (1.-alphas)), ((1.-delta_1s-delta_2s) * alphas));
            //if(fabs(nl-ol)>1.e-6 && tmpl) {
            //  Rprintf("g %d i %d k %d t %d firstcap %d p %f propp %f nl %f ol %f dist2star %f dist2 %f temp %f sigma2_scr %f cloglogp %f cloglogc %f \n",g,i,k,t,firstcap,p[i*T*K+k*T+t],propp[i*T*K+k*T+t],nl,ol,diststar[k],dist2centers[k*supN+i],temp,sigma2_scrs,cloglogp[k*T+t],cloglogc[k*T+t]);
            //  tmpl=0;
            //}
          }
          for(t=firstcap; t<T; t++) {
            propp[k*T+t] = invcloglog(cloglogp[k*T+t]-temp) * msk[k*T+t];
            propc[k*T+t] = invcloglog(cloglogc[k*T+t]-temp) * msk[k*T+t];
            //propp[k*T+t] = exp(-temp) * msk[k*T+t];
            //propc[k*T+t] = exp(-temp) * msk[k*T+t];
            nl += GETcellprob(Allhists[Hs[i]*T*K+k*T+t], propc[k*T+t],       delta_1s, delta_2s, ((1.-delta_1s-delta_2s) * (1.-alphas)), ((1.-delta_1s-delta_2s) * alphas));
            ol += GETcellprob(Allhists[Hs[i]*T*K+k*T+t], c[i*T*K+k*T+t],     delta_1s, delta_2s, ((1.-delta_1s-delta_2s) * (1.-alphas)), ((1.-delta_1s-delta_2s) * alphas));
            //if(fabs(nl-ol)>1.e-6 && tmpl) {
            //  Rprintf("g %d i %d k %d t %d firstcap %d c %f propc %f nl %f ol %f dist2star %f dist2 %f tmp %f sigma2_scr %f cloglogp %f cloglogc %f \n",g,i,k,t,firstcap,c[i*T*K+k*T+t],propc[i*T*K+k*T+t],nl,ol,diststar[k],dist2centers[k*supN+i],temp,sigma2_scrs,cloglogp[k*T+t],cloglogc[k*T+t]);
            //  tmpl=0;
            //}
          }
        } else {
          for(t=0; t<T; t++){
            propp[k*T+t]=invcloglog((cloglogp[k*T+t]-temp)) * msk[k*T+t];
            propc[k*T+t]=invcloglog((cloglogc[k*T+t]-temp)) * msk[k*T+t];   
            //p[i*T*K+k*T+t]=exp(-temp) * msk[k*T+t];
            //c[i*T*K+k*T+t]=exp(-temp) * msk[k*T+t]; 
          }
        }
      }
      //Rprintf("g %d i %d nl %f ol %f nprop %f oprop %f R %f \n",g,i,nl,ol,-log(numnn[centers[i]]),-log(numnn[centerstar]),exp(nl-log(numnn[centers[i]])-ol+log(numnn[centerstar])));
      if(runif(0.0,1.0)<exp(nl-log(numnn[centerstar])-ol+log(numnn[centers[i]]))){
        centers[i] = centerstar;
        for(k=0; k<K; k++){
          dist2centers[k*supN+i] = diststar[k];
          for(t=0; t<T; t++){
            p[i*T*K+k*T+t]=propp[k*T+t];
            c[i*T*K+k*T+t]=propc[k*T+t];              
          }
        }
        accept[dimp+1+i]+=1;
      }
    }
    
    if(*updatedelta){
      /* update alpha */
      if(datatype){
        sha = *a0alpha+FREQSUM(xs,Allhists,T*K,J,4);
        sca = *b0alpha+FREQSUM(xs,Allhists,T*K,J,3);
        alphas = rbeta(sha,sca);
      }
      /* update delta_1 and delta_2 */
      if(deltatype){
        GETDELTA(deltavect, xs, Allhists, T*K, J, 3, a0_delta); 
        delta_1s=deltavect[0];
        delta_2s=deltavect[1];   
      } else {
        sha = a0_delta[0] + FREQSUM(xs,Allhists,T*K,J,1) + FREQSUM(xs,Allhists,T*K,J,2);
        sca = a0_delta[1] + FREQSUM(xs,Allhists,T*K,J,3) + FREQSUM(xs,Allhists,T*K,J,4);
        delta_1s = rbeta(sha,sca) / 2.0;
        delta_2s = delta_1s;
      }
      /* Update psi */
      sha = (double) *a0psi + ns;
      sca= (double) *b0psi + supN - ns;
      psis = rbeta(sha,sca);
    }
    
    /* Update N */
    Ns=ns+rnbinom((double) ns,pstar);
    Nstar=Ns;
    
    ll=LIKESCR(p,c,qs,delta_1s,delta_2s,alphas,Allhists,Hs,T,K,supN,C,Ns,pstar);
    //Rprintf("g %d delta_1s %f delta_2s %f ll %f \n",g,delta_1s,delta_2s,ll);
    
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
        
        PROPFREQSCR(base,c_k,Hnew,indBasis,J,xnew,supN,T,K,p,c,C,delta_1s,delta_2s,alphas,Allhists,nprop,oprop);
        
        for(i=0; i<(dimrem+dimadd); i++){
          if(!R_FINITE(nprop[i])) {Rprintf("PROPFREQSCR nprop %f fatal error in chain %d: please report to <brett.mcclintock@noaa.gov> \n",nprop[i],*ichain); *iter = g; return;}
          if(!R_FINITE(oprop[i])) {Rprintf("PROPFREQSCR oprop %f fatal error in chain %d: please report to <brett.mcclintock@noaa.gov> \n",oprop[i],*ichain); *iter = g; return;}
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
      
      nl = LIKESCR(p,c,qnew,delta_1s,delta_2s,alphas,Allhists,Hnew,T,K,supN,C,Nstar,pstar);  
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
    
    //ll=LIKESCR(p,c,qs,delta_1s,delta_2s,alphas,Allhists,Hs,T,K,supN,C,Ns,pstar);   
    
    /* Save draws according to thin specification */ 
    if((g % th)==0)  {
      
      N[(g/th - 1)]=Ns;
      
      for(j=0; j<(dimp); j++){
        beta[(g/th - 1)*(dimp)+j]=betas[j];
      }
      sigma2_scr[(g/th - 1)]=sigma2_scrs;
      delta_1[(g/th - 1)]=delta_1s;
      delta_2[(g/th - 1)]=delta_2s;
      alpha[(g/th - 1)]=alphas;
      psi[(g/th - 1)]=psis;
      
      if(*centerind){
        for(i=0; i<supN; i++){
          center[(g/th - 1)*supN+i]=centers[i];
        }
      } else {
        for(i=0; i<supN; i++){
          center[i]=centers[i];
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
      
      posterior[(g/th - 1)]=POSTERIORSCR(ll,betas,qs,deltavect,alphas,sigma2_scrs,Ns,psis,mu0,sigma2_mu0,a0_delta,*a0alpha,*b0alpha,sigma_bounds,*a0psi,*b0psi,supN,dimp,datatype,*updatedelta,deltatype,*Area);
      if(!R_FINITE(posterior[(g/th - 1)])) {Rprintf("Fatal error in chain %d: please report to <brett.mcclintock@noaa.gov> \n",*ichain); *iter = g; return;}
      
    }
    
    if(!((g+1) % n)){    
      if((g+1)<ada){
        for(i=0; i<dimp+1; i++){
          arate[i]= accept[i]/n;
          Propsd[i]=Propsd[i]*(arate[i] < *taccept)* *tuneadjust+Propsd[i]*(arate[i] >= *taccept) / *tuneadjust;
          accept[i]=0;
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
    //if(!g%256) R_CheckUserInterrupt();
  }
  /* End Markov chain */
  if(! *printlog) Rprintf("\n");
  PutRNGstate(); 
}     
/* End function MCMCloop */
