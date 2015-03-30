#include "ranlib.h"
#include "matrix.h"
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

// Define function ProbitCJSC to draw samples from the posterior distribution

void ProbitCJSC(int *ichain, double *pbeta0, double *pprec0, double *pbeta, double *phibeta0, double *phiprec0, double *phibeta, double *zp, double *sigma2_zp, double *zphi, double *sigma2_zphi, double *delta_1, double *delta_2, double *alpha, int *x, double *psi, int *H, int *z,
              int *noccas, int *M, double *a0delta, double *a0alpha, double *b0alpha, double *l0p, double *d0p, double *l0phi, double *d0phi,
              double *loglike,
              int *nHists, int *Allhists, int *C, int *L, int *indBasis, int *ncolBasis, int *knownx, double *DMp, double *DMphi, int *pdim, int *phidim,
              int *iter, int *thin, int *numbasis,
              int *modp_h, int *modphi_h, int *data_type, int *zpind, int *zphiind, int *zind, int *Hind, int *printlog)
{
  
  GetRNGstate(); 

  /* Declare functions */
  double LIKEProbitCJS();
  double FREQSUMCJS();
  double INVPROBIT();
  int GETCK();
  int sample();
  void BETAUPDATE();
  void GETZ();
  void GETTILDE();
  void PROPFREQProbitCJS();
  void GETDELTACJS();

  int T = *noccas-1;
  double Ttmp = (double) T;
  int supN = *M; 
  int Mk = supN*(T);      
  int datatype = *data_type;
    
  int niter, th;
  int J = *nHists;
  
  niter = *iter;              /* Number of iterations in the Markov chain */
  th = *thin;                 /* Number of iterations for thinning */
  
  int dimp = pdim[1];               /* length of p beta vector */
  int dimphi = phidim[1];           /* length of phi beta vector */
  
  double pbetas[dimp], phibetas[dimphi], zps[supN], zphis[supN], sigma2_zps, sigma2_zphis, alphas, delta_1s, delta_2s, psis;
  int xs[J], xnew[J], knownxs[J], Hs[supN], Hnew[supN];
  int qs[supN], qnew[supN]; 
  
  double sha, sca;
  int base, nbasis;

  int t, i, j, k, g=0;
  double op, np;
  int dimrem, dimadd;
  
  for(k=0; k<dimp; k++){
    pbetas[k]=pbeta[k];
  }
  for(k=0; k<dimphi; k++){
    phibetas[k]=phibeta[k];
  }
  
  sigma2_zps= ( *modp_h ? sigma2_zp[0] : (double) 0.0);
  sigma2_zphis= ( *modphi_h ? sigma2_zphi[0] : (double) 0.0);
  double preczps= ( *modp_h ? 1./sigma2_zps : (double) 0.0);
  double preczphis= ( *modphi_h ? 1./sigma2_zphis : (double) 0.0);
  
  double zps2=0.;
  double zphis2=0.;
  delta_1s=delta_1[0];
  delta_2s=delta_2[0];
  psis=psi[0];

  double ns=0.;
  for(i=0; i< supN; i++)  {
    zps[i] = ( *modp_h ? zp[i] : (double) 0.0);
    zps2 += (zps[i]*zps[i]);
    zphis[i] = ( *modphi_h ? zphi[i] : (double) 0.0);
    zphis2 += (zphis[i]*zphis[i]);
    Hs[i] = H[i];
    qs[i] = (((C[Hs[i]]-1)<T) ? 1 : 0);
    ns += (double) qs[i];
    qnew[i] = qs[i];
    Hnew[i] = Hs[i];
  }
  double nstar = ns;
  
  for(j=0; j<J; j++) {
    xs[j]=x[j];
    xnew[j]=xs[j];
    knownxs[j]= knownx[j];
  }
  
  double up[Mk], uphi[Mk];
  int zs[supN*(T+1)], znew[supN*(T+1)];
  double propz[supN], newpropz[supN];
  
  for(i=0; i<supN; i++){
    for(t=0; t<(T+1); t++){
      zs[i*(T+1)+t] = z[i*(T+1)+t];
    }
  }
  
  double probitp[T*T], probitphi[T*T];
  
  int cohort;
  int cohortseq[T];
  cohortseq[0]=0;
  for(i=1; i<T; i++){
    cohortseq[i]=cohortseq[i-1]+T-i+1; 
  }

  for(cohort=0; cohort<T; cohort++){
    for(t=0; t<T; t++){
      probitp[cohort*T+t]=0.;
      probitphi[cohort*T+t]=0.;
    }
  }  
  for(cohort=0; cohort<T; cohort++){
    for(t=cohort; t<T; t++){
      for(k=0; k<dimp; k++){
        probitp[cohort*T+t]+=DMp[(cohortseq[cohort]+t-cohort)*dimp+k]*pbetas[k];
      }
      for(k=0; k<dimphi; k++){
        probitphi[cohort*T+t]+=DMphi[(cohortseq[cohort]+t-cohort)*dimphi+k]*phibetas[k];
      }
    }
  }
  
  GETTILDE(up,uphi,probitp,probitphi,zps,zphis,zs,T,supN,C,Hs,Allhists);
  
  double vbp,mbp;
  double vbphi,mbphi;
  
  double pprec0s[dimp*dimp];
  for(i=0; i<dimp; i++){
    for(j=0; j<dimp; j++){
      pprec0s[i*dimp+j]=pprec0[i*(dimp)+j];
    }
  }

  double phiprec0s[dimphi*dimphi];
  for(i=0; i<dimphi; i++){
    for(j=0; j<dimphi; j++){
      phiprec0s[i*dimphi+j]=phiprec0[i*(dimphi)+j];
    }
  }

  double a0_delta[3], deltavect[3];
  for(i=0; i<3; i++){
    a0_delta[i]=a0delta[i];
  }
  deltavect[0]=delta_1s;
  deltavect[1]=delta_2s;
  deltavect[2]=1.-delta_1s-delta_2s;
  
  alphas=alpha[0];
  
  double indbase[*ncolBasis];
  for(j=0; j< *ncolBasis; j++){
    indbase[j]=0.0;
  }
  int c_k;
  int ind, obasesum, nbasesum;
  
  /* Calculate the log-likelihood */  
  double ll=LIKEProbitCJS(zs,probitp,probitphi,zps,zphis,delta_1s,delta_2s,alphas,Allhists,Hs,T,supN,C);
  loglike[0]=ll;
  if(!R_FINITE(ll)) {
    Rprintf("Fatal error in chain %d: initial likelihood is '%f'. \n",*ichain,ll);
    *iter = g;
    return;
  }
  double nl;
  
  if(*printlog){
    Rprintf("Chain %d initialized \n",*ichain);     
    if(niter) Rprintf("Chain %d is %.0f%% complete \n",*ichain,(double) 100*g/niter);
  } else {
    Rprintf("Chain %d initialized \r",*ichain);     
    if(niter) Rprintf("Chain %d is %.0f%% complete \r",*ichain,(double) 100*g/niter);
  }
  
  /* Begin Markov chain */  
  for (g=1; g < (niter+1); g++)  {
  
    /* Update pbetas */
    BETAUPDATE(pbetas,up,zps,DMp,zs,dimp,T,supN,C,Hs,pbeta0,pprec0s,cohortseq,1);
   
    /* Update phibetas */
    BETAUPDATE(phibetas,uphi,zphis,DMphi,zs,dimphi,T,supN,C,Hs,phibeta0,phiprec0s,cohortseq,0);
    
    for(cohort=0; cohort<T; cohort++){
      for(t=cohort; t<T; t++){
        probitp[cohort*T+t]=0.;
        probitphi[cohort*T+t]=0.;
        for(k=0; k<dimp; k++){
          probitp[cohort*T+t]+=DMp[(cohortseq[cohort]+t-cohort)*dimp+k]*pbetas[k];
        }
        for(k=0; k<dimphi; k++){
          probitphi[cohort*T+t]+=DMphi[(cohortseq[cohort]+t-cohort)*dimphi+k]*phibetas[k];
        }
      }
    }
    
    if(*modp_h){
      /* Update zp */

      zps2 = 0.0;
      for(i=0; i<supN; i++){
        mbp = 0.0;
        Ttmp = 0.;
        for(t=(C[Hs[i]]-1); t<T; t++){  
          mbp += (up[i*T+t] - probitp[(C[Hs[i]]-1)*T+t]) * zs[i*(T+1)+t+1]; 
          Ttmp += (double) zs[i*(T+1)+t+1];
        }
        vbp = 1. / (preczps + Ttmp);
        zps[i]=rnorm((vbp*mbp),sqrt(vbp));
        if(qs[i]) zps2 += zps[i]*zps[i];
      }
      
      /* Update sigma2_zp */
      sha = *l0p + ns/2.0;
      sca = *d0p + zps2/2.0;
      preczps = rgamma((double) sha, (double) (1.0/sca));
      sigma2_zps = 1. / preczps;
    }
    
    if(*modphi_h){
      /* Update zphi */
      zphis2 = 0.0;
      for(i=0; i<supN; i++){
        mbphi = 0.0;
        Ttmp = 0.;
        for(t=(C[Hs[i]]-1); t<T; t++){  
          mbphi += (uphi[i*T+t] - probitphi[(C[Hs[i]]-1)*T+t]) * zs[i*(T+1)+t]; 
          Ttmp += (double) zs[i*(T+1)+t];
        }
        vbphi = 1. / (preczphis + Ttmp);
        zphis[i]=rnorm((vbphi*mbphi),sqrt(vbphi));
        if(qs[i]) zphis2 += zphis[i]*zphis[i];
      }
      
      /* Update sigma2_zphi */
      sha = *l0phi + ns/2.0;
      sca = *d0phi + zphis2/2.0;
      preczphis = rgamma((double) sha, (double) (1.0/sca));
      sigma2_zphis = 1. / preczphis;
    }
  
    /* update alpha */
    if(datatype){
      sha = *a0alpha+FREQSUMCJS(xs,Allhists,T,J,4,C);
      sca = *b0alpha+FREQSUMCJS(xs,Allhists,T,J,3,C);
      alphas = rbeta(sha,sca);
    }
    
    /* update delta_1 and delta_2 */
    GETDELTACJS(deltavect, xs, Allhists, T, J, 3, a0_delta,C); 
    delta_1s=deltavect[0];
    delta_2s=deltavect[1];    
      
    /* update z */
    for(i=0; i<supN; i++){
      GETZ(i,zs,T,probitp,probitphi,zps,zphis,C,L,Hs[i],propz);
      for(t=0; t<(T+1); t++){
        znew[i*(T+1)+t] = zs[i*(T+1)+t];
      }
      newpropz[i]=propz[i];
    }
    
    ll=LIKEProbitCJS(zs,probitp,probitphi,zps,zphis,delta_1s,delta_2s,alphas,Allhists,Hs,T,supN,C);
  
    /* update x and H (latent history frequencies) */
    op=0.0;
    np=0.0;
    nbasis = GETCK(*numbasis,0);
    obasesum=0;
    ind=1;
    for(k=0; k< *ncolBasis; k++){
      if((xs[indBasis[k*3+2]]+min(xs[indBasis[k*3]],xs[indBasis[k*3+1]])) && !(xs[indBasis[k*3+2]]==knownxs[indBasis[k*3+2]] && min(xs[indBasis[k*3]],xs[indBasis[k*3+1]])==min(knownxs[indBasis[k*3]],knownxs[indBasis[k*3+1]]))){
        indbase[k]=1.0;
        obasesum+=1;
      } else {
        indbase[k]=0.0;
      }
    }
    
    for(j=0; j< nbasis; j++){

      if(obasesum){
        base=sample(*ncolBasis, indbase);
        c_k=GETCK(xnew[indBasis[base*3+2]]+min(xnew[indBasis[base*3]],xnew[indBasis[base*3+1]]),xnew[indBasis[base*3+2]]) - xnew[indBasis[base*3+2]];
        if(c_k>0){
          dimadd=c_k;
          dimrem=2*c_k;
        } else {
          dimadd=2*(-c_k);
          dimrem= -c_k;
        }
        double nprop[dimadd+dimrem];
        double oprop[dimadd+dimrem];
        
        PROPFREQProbitCJS(base,c_k,Hnew,indBasis,J,xnew,supN,T,znew,probitp,probitphi,zps,zphis,C,L,delta_1s,delta_2s,alphas,Allhists,nprop,oprop,newpropz);
        
        for(i=0; i<(dimrem+dimadd); i++){
          if(!R_FINITE(nprop[i])) {Rprintf("PROPFREQ nprop %f fatal error in chain %d: please report to <brett.mcclintock@noaa.gov> \n",nprop[i],*ichain); *iter = g; return;}
          if(!R_FINITE(oprop[i])) {Rprintf("PROPFREQ oprop %f fatal error in chain %d: please report to <brett.mcclintock@noaa.gov> \n",oprop[i],*ichain); *iter = g; return;}
          op += nprop[i];
          np += oprop[i]; 
        }
        nbasesum=0;
        for(i=0; i< *ncolBasis; i++){
          if(xnew[indBasis[i*3+2]]+min(xnew[indBasis[i*3]],xnew[indBasis[i*3+1]])){
            indbase[i]=1.0;
            nbasesum+=1;
          } else {
            indbase[i]=0.0;
          }
        }
        for(k=0; k<3; k++){
          ind = ( (xnew[indBasis[base*3+k]]<knownxs[indBasis[base*3+k]]) ? (int) 0 : ind);
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
        qnew[i] = (((C[Hnew[i]]-1)<T) ? (int) 1 : (int) 0);//((Hnew[i]) ? (int) 1 : (int) 0);//
        nstar += (double) qnew[i];
        op += dbinom((double) qs[i],1.0,psis,1);
        np += dbinom((double) qnew[i],1.0,psis,1);
        GETZ(i,znew,T,probitp,probitphi,zps,zphis,C,L,Hnew[i],newpropz);
        op += newpropz[i];
        np += propz[i];
      }
      
      nl = LIKEProbitCJS(znew,probitp,probitphi,zps,zphis,delta_1s,delta_2s,alphas,Allhists,Hnew,T,supN,C);
      np += nl;
      
      if(runif(0.0,1.0)<exp(np-op)){
        for(i=0; i<J; i++){
          xs[i]=xnew[i];
        }
        ns=nstar;
        for(i=0; i<supN; i++){
          Hs[i]=Hnew[i];
          qs[i]=qnew[i];
          for(t=0; t<(T+1); t++){
            zs[i*(T+1)+t] = znew[i*(T+1)+t];
          }
          propz[i]=newpropz[i];
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
    nstar=ns;
    
    /* Update psi */
    sha = (double) 1.e-6 + ns;
    sca= (double) 1.0 + supN - ns;
    psis = rbeta(sha,sca);
    
    /* update up and uphi */    
    GETTILDE(up,uphi,probitp,probitphi,zps,zphis,zs,T,supN,C,Hs,Allhists);
        
     /* Save draws according to thin specification */ 
    if((g % th)==0)  {

      for(j=0; j<(dimp); j++){
        pbeta[(g/th - 1)*(dimp)+j]=pbetas[j];
      }
      for(j=0; j<(dimphi); j++){
        phibeta[(g/th - 1)*(dimphi)+j]=phibetas[j];
      }
      sigma2_zp[(g/th - 1)]=sigma2_zps;
      sigma2_zphi[(g/th - 1)]=sigma2_zphis;
      delta_1[(g/th - 1)]=delta_1s;
      delta_2[(g/th - 1)]=delta_2s;
      alpha[(g/th - 1)]=alphas;
      psi[(g/th - 1)]=psis;
      
      if(*zpind){
        for(i=0; i<supN; i++){
          zp[(g/th - 1)*supN+i]=zps[i];
        }
      } else {
        for(i=0; i<supN; i++){
          zp[i]=zps[i];
        }        
      }
      if(*zphiind){
        for(i=0; i<supN; i++){
          zphi[(g/th - 1)*supN+i]=zphis[i];
        }
      } else {
        for(i=0; i<supN; i++){
          zphi[i]=zphis[i];
        }        
      }
      
      if(*zind){
        for(i=0; i<supN; i++){
          for(t=0; t<(T+1); t++){
            z[(g/th - 1)*supN*(T+1)+i*(T+1)+t]=zs[i*(T+1)+t];
          }
        }
      } else {
        for(i=0; i<supN; i++){
          for(t=0; t<(T+1); t++){
            z[i*(T+1)+t]=zs[i*(T+1)+t];
          }
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
      
      loglike[(g/th - 1)]=ll; 
      if(!R_FINITE(loglike[(g/th - 1)])) {Rprintf("Fatal error in chain %d: please report to <brett.mcclintock@noaa.gov> \n",*ichain); *iter = g; return;}
      
    }

    if(!(g%(niter/ min(niter,100)))){
      if(*printlog){
        Rprintf("Chain %d is %.0f%% complete \n",*ichain,(double) 100*g/niter);      
      } else {
        Rprintf("Chain %d is %.0f%% complete \r",*ichain,(double) 100*g/niter);
      }
    }
  }
  /* End Markov chain */
  if(! *printlog) Rprintf("\n");
  PutRNGstate(); 
}     
/* End function MCMCloop */



/* FUNCTION DEFINITIONS */

double INVPROBIT(double q, double mean, double sd, int lower_tail, int log_p)
{
  double invprobit = fmin(1.-tol,fmax(tol,pnorm(q,mean,sd,lower_tail,log_p)));
  return(invprobit);
}

void GETTILDE(double *up, double *uphi, double *probitp, double *probitphi, double *zps, double *zphis, int *zs, int T, int supN, int *C, int *Hs, int *Allhists)
{  
  int i,t;
  double muzp,muzphi;
  int indhist;
  double trunmuzp,trunmuzphi;
  
  for(i=0; i<supN; i++){
    for(t=0; t<min(C[Hs[i]]-1,T); t++){
      muzp=probitp[t]+zps[i];
      muzphi=probitphi[t]+zphis[i];
      up[i*T+t] = rnorm(muzp,1.0);
      uphi[i*T+t] = rnorm(muzphi,1.0);   
    }
    for(t=(C[Hs[i]]-1); t<T; t++){
      indhist=Allhists[Hs[i] * (T+1) + t +1];
      muzp=probitp[(C[Hs[i]]-1)*T+t]+zps[i];
      muzphi=probitphi[(C[Hs[i]]-1)*T+t]+zphis[i];
      trunmuzp=INVPROBIT(0.0,muzp,1.0,1,0); 
      trunmuzphi=INVPROBIT(0.0,muzphi,1.0,1,0);  
      up[i*T+t] = ( (indhist>0) ? qnorm(runif(trunmuzp,1.0),muzp,1.0,1,0) :  qnorm(runif(0.0,trunmuzp),muzp,1.0,1,0) );
      uphi[i*T+t] = ( (zs[i*(T+1)+t+1]>0) ? qnorm(runif(trunmuzphi,1.0),muzphi,1.0,1,0) :  qnorm(runif(0.0,trunmuzphi),muzphi,1.0,1,0) );
    }
  }
}  

void BETAUPDATE(double *beta, double *u, double *z, double *DM, int *zs, int dim, int T, int supN, int *C, int *Hs, double *beta0, double *prec0, int *cohortseq, int pind)
{
  void mvnorm();
  void matrix_invert();
  Matrix *matrix_alloc();
  
  int i,j,k,t;
  double vbeta[dim*dim];
  Matrix *invvbeta = matrix_alloc(dim,dim); 
  double tempmbeta[dim];
  double mbeta[dim];

  for(i=0; i<dim; i++){
    for(j=0; j<dim; j++){
      invvbeta->matrix_entry[i][j]=0.;
      for(k=0; k<supN; k++){
        for(t=(C[Hs[k]]-1); t<T; t++){
          invvbeta->matrix_entry[i][j] += DM[(cohortseq[C[Hs[k]]-1]+t-(C[Hs[k]]-1))*dim+i] * DM[(cohortseq[C[Hs[k]]-1]+t-(C[Hs[k]]-1))*dim+j] * zs[k*(T+1)+t+pind] * (1. + prec0[i*dim+j]);
        }
      }
    }
  } 
  matrix_invert(invvbeta);
  for(i=0; i<dim; i++){
    for(j=0; j<dim; j++){
      vbeta[i*dim+j]= invvbeta->matrix_entry[i][j];
    }
  }
  for(i=0; i<dim; i++){
    tempmbeta[i]=0.;
    for(k=0; k<supN; k++){
      for(t=(C[Hs[k]]-1); t<T; t++){
        tempmbeta[i] += DM[(cohortseq[C[Hs[k]]-1]+t-(C[Hs[k]]-1))*dim+i] * (u[k*T+t]-z[k])  * zs[k*(T+1)+t+pind];
        for(j=0; j<dim; j++){
          tempmbeta[i] +=  DM[(cohortseq[C[Hs[k]]-1]+t-(C[Hs[k]]-1))*dim+i] * DM[(cohortseq[C[Hs[k]]-1]+t-(C[Hs[k]]-1))*dim+j] * zs[k*(T+1)+t+pind] * prec0[i*dim+j] * beta0[j];
        }
      }
    }
  }
  for(i=0; i<dim; i++){
    mbeta[i]=0.;
    for(j=0; j<dim; j++){
      mbeta[i]+=vbeta[i*dim+j]*tempmbeta[j];
    }
  }
  mvnorm(dim, mbeta, vbeta, beta); 
}

void GETZ(int i, int *z, int T, double *probitp, double *probitphi, double *zp, double *zphi, int *C, int *L, int Hind, double *propz)
{
  int t;
  int firstcap, lastcap;
  double p,phi,num,probz;

  propz[i]=0.;
  firstcap = C[Hind]-1;
  lastcap = L[Hind]-1;
  for(t=0; t<firstcap; t++){
    z[i*(T+1)+t] = 0;
  }
  if(Hind){
    z[i*(T+1)+firstcap] = 1;
    if(firstcap<T){
      for(t=firstcap; t<lastcap+1; t++){
        z[i*(T+1)+t] = 1;
      }
      for(t=lastcap+1; t<T+1; t++){
        p = INVPROBIT((probitp[firstcap*T+t-1] + zp[i]),0.0,1.0,1,0);
        phi = INVPROBIT((probitphi[firstcap*T+t-1] + zphi[i]),0.0,1.0,1,0); 
        num = phi * z[i*(T+1)+t-1] * (1.-p);
        if(t<T) num *= (1.-INVPROBIT((probitphi[firstcap*T+t] + zphi[i]),0.0,1.0,1,0));
        probz = num / (num + (1.-phi*z[i*(T+1)+t-1]));
        if(t<T){
          if(z[i*(T+1)+t+1]) probz = 1.;
        }
        z[i*(T+1)+t] = (int) rbinom(1.0,probz);
        propz[i] += dbinom((double) z[i*(T+1)+t],1.0,probz,1);
      }
    }
  }
}

double FREQSUMCJS(int *x, int *Allhists, int T, int J, int ind, int *C)
{
  int j, t;
  int freqsum=0;
  for(j=0; j<J; j++){
    for(t=(C[j]-1); t<T; t++){
      freqsum += (Allhists[j*(T+1)+t+1]==ind)*x[j];
    }
  }
  return(freqsum); 
}

/* Define function GETDELTA for updating delta by drawing from the full conditional posterior distribution */
void GETDELTACJS(double *deltavect, int *xs, int *Allhists, int T, int J, int dim, double *a0delta, int *C) 
{
  int k, kk;
  double nu[dim];
  
  for (kk=0; kk < dim; kk++)  {
    nu[kk]=FREQSUMCJS(xs,Allhists,T,J,kk+1,C);
  }
  nu[dim-1]+=FREQSUMCJS(xs,Allhists,T,J,dim+1,C);
      
  double xx[dim], sumx=0.0;
  for (k = 0; k < dim; k++) {
    xx[k]=rgamma((nu[k]+a0delta[k]),1.0);
    sumx += xx[k];
  }

  for (k = 0; k < dim; k++) {
    deltavect[k]=xx[k]/sumx;
  }
}

double LIKEProbitCJS(int *z, double *probitp, double *probitphi, double *zp, double *zphi, double delta_1, double delta_2, double alpha, int *Allhists, int *Hs, int T, int supN, int *C)
{
  int i,t;
  double logdens=0.;
  int indhist;
  int firstcap;
  double p, phi;
  for(i=0; i<supN; i++)  {
    firstcap = C[Hs[i]]-1;
    for(t=firstcap; t<T; t++){
      if(z[i*(T+1)+t]){
        indhist = Allhists[Hs[i] * (T+1) + t + 1];
        p = INVPROBIT((probitp[firstcap*T+t] + zp[i]),0.0,1.0,1,0);
        phi = INVPROBIT((probitphi[firstcap*T+t] + zphi[i]),0.0,1.0,1,0);
        logdens += log( (indhist==0) * ((1.-p) * phi * z[i*(T+1)+t+1] + (1.-phi)*(1.-z[i*(T+1)+t+1]))
                      + (indhist==1) * p * delta_1 * phi
                      + (indhist==2) * p * delta_2 * phi
                      + (indhist==3) * p * (1.-delta_1-delta_2) * (1.-alpha) * phi
                      + (indhist==4) * p * (1.-delta_1-delta_2) * alpha * phi );
      }
    }     
  }
  return(logdens); 
}

double GETprodhProbitCJS(int *Allhists, int *z, double *probitp, double *probitphi, double *zp, double *zphi, int *C, double delta_1, double delta_2, double alpha, int j,int T, int i, double propz)
{
  int t;
  double logdens=0.0, dens;
  int indhist;
  int firstcap = C[j]-1;
  double p, phi;
  
  for(t=firstcap; t<T; t++){
    if(z[i*(T+1)+t]){
      indhist = Allhists[j * (T+1) + t + 1];
      p = INVPROBIT((probitp[firstcap*T+t] + zp[i]),0.0,1.0,1,0);
      phi = INVPROBIT((probitphi[firstcap*T+t] + zphi[i]),0.0,1.0,1,0);
      logdens += log( (indhist==0) * ((1.-p) * phi * z[i*(T+1)+t+1] + (1.-phi)*(1.-z[i*(T+1)+t+1]))
                    + (indhist==1) * p * delta_1 * phi
                    + (indhist==2) * p * delta_2 * phi
                    + (indhist==3) * p * (1.-delta_1-delta_2) * (1.-alpha) *phi
                    + (indhist==4) * p * (1.-delta_1-delta_2) * alpha *phi );
    }
  }  
  dens = exp(propz + logdens);
  if(dens<tol) dens = tol;
  return(dens);
}

void PROPFREQProbitCJS(int icol,int c_k,int *Hnew, int *indBasis, int J, int *xnew, int supN, int T, int *zs, double *probitp, double *probitphi, double *zp, double *zphi, int *C, int *L, double delta_1, double delta_2, double alpha, int *Allhists, double *nprop, double *oprop, double *propz)
{  
  void ProbSampleNoReplace();
  
  int remove_xi[3];
  int add_xi[3];
  int j,i,k,t;
  int absc_k=abs(c_k);
  int remove[absc_k], add[absc_k];
  double prodz[supN], prodh[supN], temppropz[supN];
  double prodzsum, prodhsum;  
  double temp = runif(0.0,1.0);
  int ztemp[(T+1)*supN];
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
          temppropz[i]=propz[i];
          if(Hnew[i]==indBasis[icol*3+j]){
            prodz[i] =  1. - GETprodhProbitCJS(Allhists,zs,probitp,probitphi,zp,zphi,C,delta_1,delta_2,alpha,indBasis[icol*3+j],T,i,propz[i])*((C[indBasis[icol*3+j]]-1)<T);
            prodzsum+=prodz[i];
          } else if(!Hnew[i]){
            for(t=0; t<(T+1); t++){
              ztemp[i*(T+1)+t] = 0;
            }
            GETZ(i,ztemp,T,probitp,probitphi,zp,zphi,C,L,indBasis[icol*3+j],temppropz);
            prodh[i] = GETprodhProbitCJS(Allhists,ztemp,probitp,probitphi,zp,zphi,C,delta_1,delta_2,alpha,indBasis[icol*3+j],T,i,temppropz[i]);
            prodhsum+=prodh[i];
          }          
        }
        ProbSampleNoReplace(supN, prodz, absc_k, remove); 
        for(k=0; k<absc_k; k++){
          Hnew[remove[k]]=0;
          prodh[remove[k]] = 1. - prodz[remove[k]]*((C[indBasis[icol*3+j]]-1)<T);
          for(t=0; t<T+1; t++){
            zs[remove[k]*(T+1)+t] = 0;
          }
          prodhsum+=prodh[remove[k]];    
        }
        for(k=0; k<absc_k; k++){
          nprop[count] = log(prodz[remove[k]])-log(prodzsum);
          oprop[count] = log(prodh[remove[k]])-log(prodhsum);// + propz[remove[k]];
          propz[remove[k]] = 0.;
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
          temppropz[i]=propz[i];
          if(!Hnew[i]){
            for(t=0; t<(T+1); t++){
              ztemp[i*(T+1)+t] = 0;
            }
            GETZ(i,ztemp,T,probitp,probitphi,zp,zphi,C,L,indBasis[icol*3+j],temppropz);
            prodh[i] = GETprodhProbitCJS(Allhists,ztemp,probitp,probitphi,zp,zphi,C,delta_1,delta_2,alpha,indBasis[icol*3+j],T,i,temppropz[i]);
            prodhsum+=prodh[i];
          } else if(Hnew[i]==indBasis[icol*3+j]){
            prodz[i] = 1. - GETprodhProbitCJS(Allhists,zs,probitp,probitphi,zp,zphi,C,delta_1,delta_2,alpha,indBasis[icol*3+j],T,i,propz[i])*((C[indBasis[icol*3+j]]-1)<T);
            prodzsum+=prodz[i];
          }
        }
        ProbSampleNoReplace(supN, prodh, absc_k, add);
        for(k=0; k<absc_k; k++){
          Hnew[add[k]]=indBasis[icol*3+j];
          prodz[add[k]] = 1. - prodh[add[k]]*((C[indBasis[icol*3+j]]-1)<T);
          for(t=0; t<T+1; t++){
            zs[add[k]*(T+1)+t] = ztemp[add[k]*(T+1)+t];
          }
          propz[add[k]] = temppropz[add[k]];
          prodzsum+=prodz[add[k]];
        }
        for(k=0; k<absc_k; k++){
          nprop[count] = log(prodh[add[k]])-log(prodhsum);// + propz[add[k]];  
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

void mvnorm(int dim, double *mean, double *covar, double *answer)
{

void setgmn();
void genmn();

long i,j,p=dim;
float mmean[p],ccovar[p * p],param[(p *(p+3)/2+1)],temp[p],work[p];

    for(i=0; i< p; i++) mmean[i] = (float) mean[i];
    for(i=0; i< p*p; i++) ccovar[i] = (float) covar[i];

    setgmn(mmean,ccovar,p,param);
    genmn(param,work,temp);
    for(j=0; j<p; j++) *(answer+j) = *(work+j);
}


void spofa(float *a,long lda,long n,long *info)
/*
     SPOFA FACTORS A REAL SYMMETRIC POSITIVE DEFINITE MATRIX.
     SPOFA IS USUALLY CALLED BY SPOCO, BUT IT CAN BE CALLED
     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
     (TIME FOR SPOCO) = (1 + 18/N)*(TIME FOR SPOFA) .
     ON ENTRY
        A       REAL(LDA, N)
                THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE
                DIAGONAL AND UPPER TRIANGLE ARE USED.
        LDA     INTEGER
                THE LEADING DIMENSION OF THE ARRAY  A .
        N       INTEGER
                THE ORDER OF THE MATRIX  A .
     ON RETURN
        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R
                WHERE  TRANS(R)  IS THE TRANSPOSE.
                THE STRICT LOWER TRIANGLE IS UNALTERED.
                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
        INFO    INTEGER
                = 0  FOR NORMAL RETURN.
                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
     LINPACK.  THIS VERSION DATED 08/14/78 .
     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
     SUBROUTINES AND FUNCTIONS
     BLAS SDOT
     FORTRAN SQRT
     INTERNAL VARIABLES
*/
{
extern float sdot(long n,float *sx,long incx,float *sy,long incy);
static long j,jm1,k;
static float t,s;
/*
     BEGIN BLOCK WITH ...EXITS TO 40
*/
    for(j=1; j<=n; j++) {
        *info = j;
        s = 0.0;
        jm1 = j-1;
        if(jm1 < 1) goto S20;
        for(k=0; k<jm1; k++) {
            t = *(a+k+(j-1)*lda)-sdot(k,(a+k*lda),1L,(a+(j-1)*lda),1L);
            t /=  *(a+k+k*lda);
            *(a+k+(j-1)*lda) = t;
            s += (t*t);
        }
S20:
        s = *(a+j-1+(j-1)*lda)-s;
/*
     ......EXIT
*/
        if(s <= 0.0) goto S40;
        *(a+j-1+(j-1)*lda) = sqrt(s);
    }
    *info = 0;
S40:
    return;
}

float sdot(long n,float *sx,long incx,float *sy,long incy)
{
static long i,ix,iy,m,mp1;
static float sdot,stemp;
    stemp = sdot = 0.0;
    if(n <= 0) return sdot;
    if(incx == 1 && incy == 1) goto S20;
    ix = iy = 1;
    if(incx < 0) ix = (-n+1)*incx+1;
    if(incy < 0) iy = (-n+1)*incy+1;
    for(i=1; i<=n; i++) {
        stemp += (*(sx+ix-1)**(sy+iy-1));
        ix += incx;
        iy += incy;
    }
    sdot = stemp;
    return sdot;
S20:
    m = n % 5L;
    if(m == 0) goto S40;
    for(i=0; i<m; i++) stemp += (*(sx+i)**(sy+i));
    if(n < 5) goto S60;
S40:
    mp1 = m+1;
    for(i=mp1; i<=n; i+=5) stemp += (*(sx+i-1)**(sy+i-1)+*(sx+i)**(sy+i)+*(sx+i
      +1)**(sy+i+1)+*(sx+i+2)**(sy+i+2)+*(sx+i+3)**(sy+i+3));
S60:
    sdot = stemp;
    return sdot;
}

void setgmn(float *meanv,float *covm,long p,float *parm)
/*
**********************************************************************
     void setgmn(float *meanv,float *covm,long p,float *parm)
            SET Generate Multivariate Normal random deviate
                              Function
      Places P, MEANV, and the Cholesky factoriztion of COVM
      in GENMN.
                              Arguments
     meanv --> Mean vector of multivariate normal distribution.
     covm   <--> (Input) Covariance   matrix    of  the  multivariate
                 normal distribution
                 (Output) Destroyed on output
     p     --> Dimension of the normal, or length of MEANV.
     parm <-- Array of parameters needed to generate multivariate norma
                deviates (P, MEANV and Cholesky decomposition of
                COVM).
                1 : 1                - P
                2 : P + 1            - MEANV
                P+2 : P*(P+3)/2 + 1  - Cholesky decomposition of COVM
               Needed dimension is (p*(p+3)/2 + 1)
**********************************************************************
*/
{
extern void spofa(float *a,long lda,long n,long *info);
//static long T1;
static long i,icount,info,j,D2,D3,D4,D5;
//    T1 = p*(p+3)/2+1;
/*
     TEST THE INPUT
*/
    if(!(p <= 0)) goto S10;
    Rprintf("P nonpositive in SETGMN: please report to <brett.mcclintock@noaa.gov> \n");
    return;
S10:
    *parm = p;
/*
     PUT P AND MEANV INTO PARM
*/
    for(i=2,D2=1,D3=(p+1-i+D2)/D2; D3>0; D3--,i+=D2) *(parm+i-1) = *(meanv+i-2);
/*
      Cholesky decomposition to find A s.t. trans(A)*(A) = COVM
*/
    spofa(covm,p,p,&info);
    if(!(info != 0)) goto S30;
    Rprintf(" COVM not positive definite in SETGMN: please report to <brett.mcclintock@noaa.gov> \n");
    return;
S30:
    icount = p+1;
/*
     PUT UPPER HALF OF A, WHICH IS NOW THE CHOLESKY FACTOR, INTO PARM
          COVM(1,1) = PARM(P+2)
          COVM(1,2) = PARM(P+3)
                    :
          COVM(1,P) = PARM(2P+1)
          COVM(2,2) = PARM(2P+2)  ...
*/
    for(i=1,D4=1,D5=(p-i+D4)/D4; D5>0; D5--,i+=D4) {
        for(j=i-1; j<p; j++) {
            icount += 1;
            *(parm+icount-1) = *(covm+i-1+j*p);
        }
    }
}

float ranf(void)
/*
**********************************************************************
     float ranf(void)
                RANDom number generator as a Function
     Returns a random floating point number from a uniform distribution
     over 0 - 1 (endpoints of this interval are not returned) using the
     current generator
     This is a transcription from Pascal to Fortran of routine
     Uniform_01 from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
**********************************************************************
*/
{
static float ranf;
/*
     4.656613057E-10 is 1/M1  M1 is set in a data statement in IGNLGI
      and is currently 2147483563. If M1 changes, change this also.
*/
    ranf = ignlgi()*4.656613057E-10;
    return ranf;
}

float snorm(void)
/*
**********************************************************************


     (STANDARD-)  N O R M A L  DISTRIBUTION


**********************************************************************
**********************************************************************

     FOR DETAILS SEE:

               AHRENS, J.H. AND DIETER, U.
               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM
               SAMPLING FROM THE NORMAL DISTRIBUTION.
               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.

     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'
     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)

     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
     SUNIF.  The argument IR thus goes away.

**********************************************************************
     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
*/
{
static float a[32] = {
    0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
    0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
    0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
    1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
    1.862732,2.153875
};
static float d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
    0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
    0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
    0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
};
static float t[31] = {
    7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
    1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
    2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
    4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
    9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
};
static float h[31] = {
    3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
    4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
    4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
    5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
    8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
};
static long i;
static float snorm,u,s,ustar,aa,w,y,tt;
    u = ranf();
    s = 0.0;
    if(u > 0.5) s = 1.0;
    u += (u-s);
    u = 32.0*u;
    i = (long) (u);
    if(i == 32) i = 31;
    if(i == 0) goto S100;
/*
                                START CENTER
*/
    ustar = u-(float)i;
    aa = *(a+i-1);
S40:
    if(ustar <= *(t+i-1)) goto S60;
    w = (ustar-*(t+i-1))**(h+i-1);
S50:
/*
                                EXIT   (BOTH CASES)
*/
    y = aa+w;
    snorm = y;
    if(s == 1.0) snorm = -y;
    return snorm;
S60:
/*
                                CENTER CONTINUED
*/
    u = ranf();
    w = u*(*(a+i)-aa);
    tt = (0.5*w+aa)*w;
    goto S80;
S70:
    tt = u;
    ustar = ranf();
S80:
    if(ustar > tt) goto S50;
    u = ranf();
    if(ustar >= u) goto S70;
    ustar = ranf();
    goto S40;
S100:
/*
                                START TAIL
*/
    i = 6;
    aa = *(a+31);
    goto S120;
S110:
    aa += *(d+i-1);
    i += 1;
S120:
    u += u;
    if(u < 1.0) goto S110;
    u -= 1.0;
S140:
    w = u**(d+i-1);
    tt = (0.5*w+aa)*w;
    goto S160;
S150:
    tt = u;
S160:
    ustar = ranf();
    if(ustar > tt) goto S50;
    u = ranf();
    if(ustar >= u) goto S150;
    u = ranf();
    goto S140;
}

void gsrgs(long getset,long *qvalue)
/*
**********************************************************************
     void gsrgs(long getset,long *qvalue)
               Get/Set Random Generators Set
     Gets or sets whether random generators set (initialized).
     Initially (data statement) state is not set
     If getset is 1 state is set to qvalue
     If getset is 0 state returned in qvalue
**********************************************************************
*/
{
static long qinit = 0;

    if(getset == 0) *qvalue = qinit;
    else qinit = *qvalue;
}

void gssst(long getset,long *qset)
/*
**********************************************************************
     void gssst(long getset,long *qset)
          Get or Set whether Seed is Set
     Initialize to Seed not Set
     If getset is 1 sets state to Seed Set
     If getset is 0 returns T in qset if Seed Set
     Else returns F in qset
**********************************************************************
*/
{
static long qstate = 0;
    if(getset != 0) qstate = 1;
    else  *qset = qstate;
}

void gscgn(long getset,long *g)
/*
**********************************************************************
     void gscgn(long getset,long *g)
                         Get/Set GeNerator
     Gets or returns in G the number of the current generator
                              Arguments
     getset --> 0 Get
                1 Set
     g <-- Number of the current random number generator (1..32)
**********************************************************************
*/
{
#define numg 32L
static long curntg = 1;
    if(getset == 0) *g = curntg;
    else  {
        if(*g < 0 || *g > numg) {
            Rprintf(" Generator number out of range in GSCGN: please report to <brett.mcclintock@noaa.gov> \n");
            return;
        }
        curntg = *g;
    }
#undef numg
}

long Xm1,Xm2,Xa1,Xa2,Xcg1[32],Xcg2[32],Xa1w,Xa2w,Xig1[32],Xig2[32],Xlg1[32],
    Xlg2[32],Xa1vw,Xa2vw;
long Xqanti[32];

void inrgcm(void)
/*
**********************************************************************
     void inrgcm(void)
          INitialize Random number Generator CoMmon
                              Function
     Initializes common area  for random number  generator.  This saves
     the  nuisance  of  a  BLOCK DATA  routine  and the  difficulty  of
     assuring that the routine is loaded with the other routines.
**********************************************************************
*/
{
#define numg 32L
extern void gsrgs(long getset,long *qvalue);
extern long Xm1,Xm2,Xa1,Xa2,Xa1w,Xa2w,Xa1vw,Xa2vw;
extern long Xqanti[];
static long T1;
static long i;
/*
     V=20;                            W=30;
     A1W = MOD(A1**(2**W),M1)         A2W = MOD(A2**(2**W),M2)
     A1VW = MOD(A1**(2**(V+W)),M1)    A2VW = MOD(A2**(2**(V+W)),M2)
   If V or W is changed A1W, A2W, A1VW, and A2VW need to be recomputed.
    An efficient way to precompute a**(2*j) MOD m is to start with
    a and square it j times modulo m using the function MLTMOD.
*/
    Xm1 = 2147483563L;
    Xm2 = 2147483399L;
    Xa1 = 40014L;
    Xa2 = 40692L;
    Xa1w = 1033780774L;
    Xa2w = 1494757890L;
    Xa1vw = 2082007225L;
    Xa2vw = 784306273L;
    for(i=0; i<numg; i++) *(Xqanti+i) = 0;
    T1 = 1;
/*
     Tell the world that common has been initialized
*/
    gsrgs(1L,&T1);
#undef numg
}

void initgn(long isdtyp)
/*
**********************************************************************
     void initgn(long isdtyp)
          INIT-ialize current G-e-N-erator
     Reinitializes the state of the current generator
     This is a transcription from Pascal to Fortran of routine
     Init_Generator from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     isdtyp -> The state to which the generator is to be set
          isdtyp = -1  => sets the seeds to their initial value
          isdtyp =  0  => sets the seeds to the first value of
                          the current block
          isdtyp =  1  => sets the seeds to the first value of
                          the next block
**********************************************************************
*/
{
#define numg 32L
extern void gsrgs(long getset,long *qvalue);
extern void gscgn(long getset,long *g);
extern long Xm1,Xm2,Xa1w,Xa2w,Xig1[],Xig2[],Xlg1[],Xlg2[],Xcg1[],Xcg2[];
static long g;
static long qrgnin;
/*
     Abort unless random number generator initialized
*/
    gsrgs(0L,&qrgnin);
    if(qrgnin) goto S10;
    Rprintf(" INITGN called before random number generator  initialized -- abort!: please report to <brett.mcclintock@noaa.gov> \n");
    return;
S10:
    gscgn(0L,&g);
    if(-1 != isdtyp) goto S20;
    *(Xlg1+g-1) = *(Xig1+g-1);
    *(Xlg2+g-1) = *(Xig2+g-1);
    goto S50;
S20:
    if(0 != isdtyp) goto S30;
    goto S50;
S30:
/*
     do nothing
*/
    if(1 != isdtyp) goto S40;
    *(Xlg1+g-1) = mltmod(Xa1w,*(Xlg1+g-1),Xm1);
    *(Xlg2+g-1) = mltmod(Xa2w,*(Xlg2+g-1),Xm2);
    goto S50;
S40:
    Rprintf("isdtyp not in range in INITGN: please report to <brett.mcclintock@noaa.gov> \n");
    return;
S50:
    *(Xcg1+g-1) = *(Xlg1+g-1);
    *(Xcg2+g-1) = *(Xlg2+g-1);
#undef numg
}

long mltmod(long a,long s,long m)
/*
**********************************************************************
     long mltmod(long a,long s,long m)
                    Returns (A*S) MOD M
     This is a transcription from Pascal to Fortran of routine
     MULtMod_Decompos from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     a, s, m  -->
**********************************************************************
*/
{
#define h 32768L
static long mltmod,a0,a1,k,p,q,qh,rh;
/*
     H = 2**((b-2)/2) where b = 32 because we are using a 32 bit
      machine. On a different machine recompute H
*/
    if(!(a <= 0 || a >= m || s <= 0 || s >= m)) goto S10;
    Rprintf(" a, m, s out of order in mltmod - ABORT!: please report to <brett.mcclintock@noaa.gov> \n");
    Rprintf(" mltmod requires: 0 < a < m; 0 < s < m: please report to <brett.mcclintock@noaa.gov> \n");
    return(NA_REAL);
S10:
    if(!(a < h)) goto S20;
    a0 = a;
    p = 0;
    goto S120;
S20:
    a1 = a/h;
    a0 = a-h*a1;
    qh = m/h;
    rh = m-h*qh;
    if(!(a1 >= h)) goto S50;
    a1 -= h;
    k = s/qh;
    p = h*(s-k*qh)-k*rh;
S30:
    if(!(p < 0)) goto S40;
    p += m;
    goto S30;
S40:
    goto S60;
S50:
    p = 0;
S60:
/*
     P = (A2*S*H)MOD M
*/
    if(!(a1 != 0)) goto S90;
    q = m/a1;
    k = s/q;
    p -= (k*(m-a1*q));
    if(p > 0) p -= m;
    p += (a1*(s-k*q));
S70:
    if(!(p < 0)) goto S80;
    p += m;
    goto S70;
S90:
S80:
    k = p/qh;
/*
     P = ((A2*H + A1)*S)MOD M
*/
    p = h*(p-k*qh)-k*rh;
S100:
    if(!(p < 0)) goto S110;
    p += m;
    goto S100;
S120:
S110:
    if(!(a0 != 0)) goto S150;
/*
     P = ((A2*H + A1)*H*S)MOD M
*/
    q = m/a0;
    k = s/q;
    p -= (k*(m-a0*q));
    if(p > 0) p -= m;
    p += (a0*(s-k*q));
S130:
    if(!(p < 0)) goto S140;
    p += m;
    goto S130;
S150:
S140:
    mltmod = p;
    return mltmod;
#undef h
}

void setall(long iseed1,long iseed2)
/*
**********************************************************************
     void setall(long iseed1,long iseed2)
               SET ALL random number generators
     Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
     initial seeds of the other generators are set accordingly, and
     all generators states are set to these seeds.
     This is a transcription from Pascal to Fortran of routine
     Set_Initial_Seed from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     iseed1 -> First of two integer seeds
     iseed2 -> Second of two integer seeds
**********************************************************************
*/
{
#define numg 32L
extern void gsrgs(long getset,long *qvalue);
extern void gssst(long getset,long *qset);
extern void gscgn(long getset,long *g);
extern long Xm1,Xm2,Xa1vw,Xa2vw,Xig1[],Xig2[];
static long T1;
static long g,ocgn;
static long qrgnin;
    T1 = 1;
/*
     TELL IGNLGI, THE ACTUAL NUMBER GENERATOR, THAT THIS ROUTINE
      HAS BEEN CALLED.
*/
    gssst(1,&T1);
    gscgn(0L,&ocgn);
/*
     Initialize Common Block if Necessary
*/
    gsrgs(0L,&qrgnin);
    if(!qrgnin) inrgcm();
    *Xig1 = iseed1;
    *Xig2 = iseed2;
    initgn(-1L);
    for(g=2; g<=numg; g++) {
        *(Xig1+g-1) = mltmod(Xa1vw,*(Xig1+g-2),Xm1);
        *(Xig2+g-1) = mltmod(Xa2vw,*(Xig2+g-2),Xm2);
        gscgn(1L,&g);
        initgn(-1L);
    }
    gscgn(1L,&ocgn);
#undef numg
}

long ignlgi(void)
/*
**********************************************************************
     long ignlgi(void)
               GeNerate LarGe Integer
     Returns a random integer following a uniform distribution over
     (1, 2147483562) using the current generator.
     This is a transcription from Pascal to Fortran of routine
     Random from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
**********************************************************************
*/
{
#define numg 32L
extern void gsrgs(long getset,long *qvalue);
extern void gssst(long getset,long *qset);
extern void gscgn(long getset,long *g);
extern void inrgcm(void);
extern long Xm1,Xm2,Xa1,Xa2,Xcg1[],Xcg2[];
extern long Xqanti[];
static long ignlgi,curntg,k,s1,s2,z;
static long qqssd,qrgnin;
/*
     IF THE RANDOM NUMBER PACKAGE HAS NOT BEEN INITIALIZED YET, DO SO.
     IT CAN BE INITIALIZED IN ONE OF TWO WAYS : 1) THE FIRST CALL TO
     THIS ROUTINE  2) A CALL TO SETALL.
*/
    gsrgs(0L,&qrgnin);
    if(!qrgnin) inrgcm();
    gssst(0,&qqssd);
    if(!qqssd) setall(1234567890L,123456789L);
/*
     Get Current Generator
*/
    gscgn(0L,&curntg);
    s1 = *(Xcg1+curntg-1);
    s2 = *(Xcg2+curntg-1);
    k = s1/53668L;
    s1 = Xa1*(s1-k*53668L)-k*12211;
    if(s1 < 0) s1 += Xm1;
    k = s2/52774L;
    s2 = Xa2*(s2-k*52774L)-k*3791;
    if(s2 < 0) s2 += Xm2;
    *(Xcg1+curntg-1) = s1;
    *(Xcg2+curntg-1) = s2;
    z = s1-s2;
    if(z < 1) z += (Xm1-1);
    if(*(Xqanti+curntg-1)) z = Xm1-z;
    ignlgi = z;
    return ignlgi;
#undef numg
}

void genmn(float *parm,float *x,float *work)
/*
**********************************************************************
     void genmn(float *parm,float *x,float *work)
              GENerate Multivariate Normal random deviate
                              Arguments
     parm --> Parameters needed to generate multivariate normal
               deviates (MEANV and Cholesky decomposition of
               COVM). Set by a previous call to SETGMN.
               1 : 1                - size of deviate, P
               2 : P + 1            - mean vector
               P+2 : P*(P+3)/2 + 1  - upper half of cholesky
                                       decomposition of cov matrix
     x    <-- Vector deviate generated.
     work <--> Scratch array
                              Method
     1) Generate P independent standard normal deviates - Ei ~ N(0,1)
     2) Using Cholesky decomposition find A s.t. trans(A)*A = COVM
     3) trans(A)E + MEANV ~ N(MEANV,COVM)
**********************************************************************
*/
{
static long i,icount,j,p,D1,D2,D3,D4;
static float ae;

    p = (long) (*parm);
/*
     Generate P independent normal deviates - WORK ~ N(0,1)
*/
    for(i=1; i<=p; i++) *(work+i-1) = snorm();
    for(i=1,D3=1,D4=(p-i+D3)/D3; D4>0; D4--,i+=D3) {
/*
     PARM (P+2 : P*(P+3)/2 + 1) contains A, the Cholesky
      decomposition of the desired covariance matrix.
          trans(A)(1,1) = PARM(P+2)
          trans(A)(2,1) = PARM(P+3)
          trans(A)(2,2) = PARM(P+2+P)
          trans(A)(3,1) = PARM(P+4)
          trans(A)(3,2) = PARM(P+3+P)
          trans(A)(3,3) = PARM(P+2-1+2P)  ...
     trans(A)*WORK + MEANV ~ N(MEANV,COVM)
*/
        icount = 0;
        ae = 0.0;
        for(j=1,D1=1,D2=(i-j+D1)/D1; D2>0; D2--,j+=D1) {
            icount += (j-1);
            ae += (*(parm+i+(j-1)*p-icount+p)**(work+j-1));
        }
        *(x+i-1) = ae+*(parm+i);
    }
}

/*
* Copyright (C) 2013- Acho Arnold
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

void matrix_print (Matrix *matrix)
{
  int i,j;
  Rprintf("\n");
  for(i = 0; i < matrix->row_size;i++)
  {
    Rprintf("\t\t");
    for(j = 0; j < matrix->col_size;j++)
    {
      Rprintf("%9.2f", matrix->matrix_entry[i][j]);
    }
    Rprintf("\n");
  }
  Rprintf("\n");
}

void matrix_print_part(Matrix *matrix, int start_index)
{
  int j,i;
  for (i = 0; i < matrix->row_size; ++i)\
  {
    for (j = start_index; j < matrix->col_size; ++j)
    {
      Rprintf("\t\t%9.2f", matrix->matrix_entry[i][j]);
    }
    Rprintf("\n");
  }
}

/*Function to create an identity matrix */
Matrix * matrix_callalloc(int matrix_size)
{
  Matrix *result = matrix_alloc(matrix_size, matrix_size);
  int i,j;
  for (i = 0; i < matrix_size; i += 1)
  {
    for (j = 0; j < matrix_size; j += 1)
    {
      if (j == i)
      {
        result->matrix_entry[i][j] = 1;
      }
      else
      {
        result->matrix_entry[i][j] = 0;
      }
    }
  }
  return result;
}

Matrix *matrix_alloc(int row_size, int col_size)
{
  int j;
  Matrix *new_matrix = malloc(sizeof(Matrix));
  //Allocating memory for the new matrix structure
  new_matrix->row_size = row_size;
  new_matrix->col_size = col_size;
  new_matrix->matrix_entry = malloc( new_matrix->row_size *sizeof(float *));
  for(j = 0 ; j < new_matrix->row_size ; j++)
  {
    new_matrix->matrix_entry[j] = malloc( new_matrix->col_size*sizeof(float) );
  }
  return new_matrix;
}

/*Copies Matrix1 into matrix2 */
void matrix_copy(Matrix *matrix1, Matrix *matrix2)
{
  int i, j;
  for (i = 0; i < matrix1->row_size; i += 1)
  {
    for (j = 0; j < matrix1->col_size; j += 1)
    {
      matrix2->matrix_entry[i][j] = matrix1->matrix_entry[i][j];
    }
  }
}

Matrix *matrix_multiply(const Matrix *matrix1, const Matrix *matrix2)
{
  int i, j,k, sum;
  if (matrix1->col_size != matrix2->row_size)
  {
    Rprintf("ERROR: The number columns of matrix1 != number of rows in matrix2!: please report to <brett.mcclintock@noaa.gov> \n");
  }
  Matrix *result = matrix_alloc( matrix1->row_size,matrix2->col_size);
  for (i = 0; i < matrix1->row_size; i += 1)
    {
    for (k = 0; k < matrix2->col_size; k += 1)
    {
      sum = 0;
      for (j = 0; j < matrix1->col_size; j += 1)
      {
        sum += matrix1->matrix_entry[i][j] * matrix2->matrix_entry[j][k];
      }
      result->matrix_entry[i][k] = sum;
    }
  }
  return result;
}

Matrix * matrix_pow(Matrix *matrix, int index)
{
  if(index == 1)
  {
    Matrix *result = matrix_alloc (matrix->row_size, matrix->col_size);
    matrix_copy(matrix, result);
    return result;
  }
  else
  {
    int i, j,k,l,sum,count;
    Matrix *temp = matrix_alloc (matrix->row_size, matrix->col_size); //Allocating space for a temporal matrix
    Matrix *result = matrix_alloc (matrix->row_size, matrix->col_size); //Allocating space for the result matrix
    matrix_copy(matrix, temp);
    count = index/2 -1;
    if (count < 1)
    {
      matrix_copy(matrix, result);
    }
    else
    {
      for (l = 0; l < count; l += 1)
      {
        for (i = 0; i < matrix->row_size; i += 1)
        {
          for (k = 0; k < matrix->col_size; k += 1)
          {
            sum = 0;
            for (j = 0; j < matrix->col_size; j += 1)
            {
              sum += (temp->matrix_entry[i][j] * matrix->matrix_entry[j][k]);
            }
            result->matrix_entry[i][k] = sum;
          }
        }
        /* Copying the result matrix into the temp matrix for further
        * multiplication */
        matrix_copy(result, temp);
      }
    }
    /* Freeing the temp matrix */
    matrix_free(temp);
    if (index%2 == 0)
    {
      Matrix *result_final = matrix_multiply(result, result);
      /* Freeing the result Matrix */
      matrix_free(result);
      return result_final;
    }
    else
    {
      Matrix *temp = matrix_multiply(matrix, result);
      Matrix *result_final = matrix_multiply(temp, result);
      /* Freeing the temp matrix */
      matrix_free(temp);
      /* Freeing the result Matrix */
      matrix_free(result);
      return result_final;
    }//End of else statement
  }
}

void matrix_free( Matrix *matrix)
{
  int j;
  for(j = 0 ; j < matrix->row_size ; j++)
  {
    free(matrix->matrix_entry[j]);
  }
  free(matrix->matrix_entry);
  free(matrix);
}

/*Function which divides all row entries by the value of a the diagonal */
void row_divide(Matrix *matrix, int pivot)
{
  int j;
  float divisor = matrix->matrix_entry[pivot][pivot],
  result;
  for(j = pivot; j < matrix->col_size; j++)
  {
    result = (matrix->matrix_entry[pivot][j] / divisor);
    matrix->matrix_entry[pivot][j] = result;
  }
}

/*Function to carry out row operations*/
void row_operation(Matrix *multiplier_matrix,Matrix *matrix, int pivot, int row_index)
{
  int j;
  float multiplier = (matrix->matrix_entry[row_index][pivot] / matrix->matrix_entry[pivot][pivot]);
  //Loop which checks if matrix is provided to store the multiplier
  if(multiplier_matrix != NULL)
  {
    multiplier_matrix ->matrix_entry[row_index][pivot] = multiplier;
  }
  for(j=0; j < matrix->col_size; j++)
  {
    matrix->matrix_entry[row_index][j] -= multiplier * matrix->matrix_entry[pivot][j];
  }
}

void matrix_row_reduce( Matrix *matrix, int zero_control )
{
  int pivot, row_index;
  //float multiplier;
  for( pivot = 0; pivot < matrix->row_size ; pivot++)
  {
    error_zeros(matrix, zero_control); //Function checks if there are too many zeros in a single row
    if(  (matrix->matrix_entry[pivot][pivot] != 1) && (matrix->matrix_entry[pivot][pivot] != 0)  )
    {
      row_divide(matrix, pivot);
    }
    for (row_index = pivot+1; row_index < matrix->row_size; row_index++)
    {
      if (matrix->matrix_entry[pivot][pivot] != 0)
      {
        row_operation(NULL,matrix, pivot, row_index);
      }
    }
    for(row_index = pivot-1; row_index >=0; row_index --)
    {
      if (matrix->matrix_entry[pivot][pivot] != 0)
      {
      row_operation(NULL,matrix, pivot, row_index);
      }
    }
  }
}

void LU_decompose(Matrix *upper_triangular, Matrix *lower_triangular)
{
  int pivot, row_index;
  //float multiplier;
  for( pivot = 0; pivot < upper_triangular->row_size ; pivot++)
  {
    error_zeros(upper_triangular, upper_triangular->col_size); //Function checks if there are too many zeros in a single row
    for (row_index = pivot+1; row_index < upper_triangular->row_size; row_index++)
    {
      if ( upper_triangular->matrix_entry[pivot][pivot] != 0)
      {
        row_operation(lower_triangular,upper_triangular, pivot, row_index);
      }
    }
  }
}

void matrix_subtract(Matrix *result, Matrix *matrix1, Matrix *matrix2)
{
  int i, j;
  if ( !(matrix_equal_size(matrix1, matrix2 )) || \
  !(matrix_equal_size(matrix2, result)))
  {
    Rprintf("ERROR: The matrices you are trying to subtract have different sizes: please report to <brett.mcclintock@noaa.gov> \n");
    return;
  }
  for(i = 0; i < matrix1->row_size; i += 1)
  {
    for (j = 0; j < matrix1->col_size; j += 1)
    {
      result->matrix_entry[i][j] = matrix1->matrix_entry[i][j] - matrix2->matrix_entry[i][j];
    }
  }
}

void matrix_add(Matrix *result, Matrix *matrix1, Matrix *matrix2)
{
  int i, j;
  if ( !(matrix_equal_size(matrix1, matrix2 )) || \
  !(matrix_equal_size(matrix2, result)))
  {
    Rprintf("ERROR: The matrices you are trying to add have different sizes: please report to <brett.mcclintock@noaa.gov> \n");
    return;
  }
  for (i = 0; i < matrix1->row_size; i += 1)
  {
    for (j = 0; j < matrix1->col_size; j += 1)
    {
    result->matrix_entry[i][j] = matrix1->matrix_entry[i][j] + matrix2->matrix_entry[i][j];
    }
  }
}

void matrix_invert(Matrix *inverse_matrix)
{
  int j,k;
  /*Temporal matrix used in this function */
  Matrix *temp_matrix = matrix_alloc(inverse_matrix->row_size, inverse_matrix->col_size *2);
  matrix_copy(inverse_matrix, temp_matrix);
  /* Adding an identity matrix at the end of the temporal matrix */
  for(j = 0; j< temp_matrix->row_size; j++)
  {
    for(k = inverse_matrix->col_size; k < temp_matrix->col_size; k++)
    {
      if( (j+inverse_matrix->col_size) == k)
      {
        temp_matrix->matrix_entry[j][k] = 1;
      }
      else
      {
        temp_matrix->matrix_entry[j][k] = 0;
      }
    }
  }
  matrix_row_reduce(temp_matrix, temp_matrix->row_size);
  /*  for(j = 0; j< inverse_matrix->row_size; j++)
  {
    for(k = 0; k < inverse_matrix->col_size; k++)
    {
      Rprintf("matrix %d %d %f \n",j+1,k+1,inverse_matrix->matrix_entry[j][k]);
    }
    for(k = 0; k < temp_matrix->col_size; k++)
    {
      Rprintf("temp_matrix %d %d %f \n",j+1,k+1,temp_matrix->matrix_entry[j][k]);
    }
  }*/
  /* Copying the inverse matrix from the temp_matrix to the invse_matrix */
  for(j = 0; j< temp_matrix->row_size; j++)
  {
    for(k = inverse_matrix->col_size; k < temp_matrix->col_size; k++)
    {
      inverse_matrix->matrix_entry[j][k-inverse_matrix->col_size] = temp_matrix->matrix_entry[j][k];
      //Rprintf("matrix %d %d %f \n",j+1,k-inverse_matrix->col_size+1,temp_matrix->matrix_entry[j][k]);
    }
  }/*
      for(j = 0; j< inverse_matrix->row_size; j++)
  {
    for(k = 0; k < inverse_matrix->col_size; k++)
    {
      Rprintf("matrix %d %d %f \n",j+1,k+1,inverse_matrix->matrix_entry[j][k]);
    }
  }*/
  matrix_free(temp_matrix);
}

int matrix_equal_size( Matrix *matrix1, Matrix *matrix2)
{
  return (matrix1->row_size == matrix2->row_size && \
  matrix1->col_size == matrix2->col_size);
}

/*
This function checks if there is a line containing too many zero's and it exits
if such a line is found
*/
void error_zeros( Matrix *matrix, int control_index)
{
  int i,j,count;
  for(i=0; i<matrix->row_size; i++)
  {
    count=0;
    for(j = 0; j < matrix->col_size; j++)
    {
      if( matrix->matrix_entry[i][j] == 0)
      {
        count++;
      }
      else
      {
        return;
      }
      if(count == control_index)
      {
        Rprintf("Process fail because row %d contains %d zeros: please report to <brett.mcclintock@noaa.gov> \n",i+1,control_index);
        matrix_print(matrix);
        return;
      }
    }
  }
}