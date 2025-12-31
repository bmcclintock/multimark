#include "ranlib.h"
#include "matrix.h"
#include "func.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#define tol 1.e-200

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

/* FUNCTION DEFINITIONS */

double INVPROBIT(double x, double mean, double sd, int lower_tail, int log_p)
{
  double invprobit = fmin(1.-tol,fmax(tol,pnorm(x,mean,sd,lower_tail,log_p)));
  return(invprobit);
}

void GETTILDE(double *up, double *uphi, double *probitp, double *probitphi, double *zps, double *zphis, int *qs, int T, int supN, int *C, int *Hs, int *Allhists)
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
      up[i*T+t] = ( (indhist>0) ? qnorm(fmin(1.-tol,fmax(runif(trunmuzp,1.0),tol)),muzp,1.0,1,0) :  qnorm(fmin(1.-tol,fmax(runif(0.0,trunmuzp),tol)),muzp,1.0,1,0) );
      uphi[i*T+t] = ( (qs[i*(T+1)+t+1]>0) ? qnorm(fmin(1.-tol,fmax(runif(trunmuzphi,1.0),tol)),muzphi,1.0,1,0) :  qnorm(fmin(1.-tol,fmax(runif(0.0,trunmuzphi),tol)),muzphi,1.0,1,0) );
    }
  }
}  

void mvnorm(int dim, double *mean, double *covar, double *answer)
{
  
  long i,j,p=dim;
  float mmean[p],ccovar[p * p],param[(p *(p+3)/2+1)],temp[p],work[p];
  
  for(i=0; i< p; i++) mmean[i] = (float) mean[i];
  for(i=0; i< p*p; i++) ccovar[i] = (float) covar[i];
  
  setgmn(mmean,ccovar,p,param);
  genmn(param,work,temp);
  for(j=0; j<p; j++) *(answer+j) = *(work+j);
}

void BETAUPDATE(double *beta, double *u, double *z, double *DM, int *qs, int dim, int T, int supN, int *C, int *Hs, double *beta0, double *prec0, int *cohortseq, int pind)
{
  int i,j,k,t;
  double vbeta[dim*dim];
  Matrix *invvbeta = matrix_alloc(dim,dim); 
  double tempmbeta[dim];
  double mbeta[dim];
  
  for(i=0; i<dim; i++){
    for(j=0; j<dim; j++){
      invvbeta->matrix_entry[i][j] = prec0[i*dim+j];
      for(k=0; k<supN; k++){
        for(t=(C[Hs[k]]-1); t<T; t++){
          invvbeta->matrix_entry[i][j] += DM[(cohortseq[C[Hs[k]]-1]+t-(C[Hs[k]]-1))*dim+i] * DM[(cohortseq[C[Hs[k]]-1]+t-(C[Hs[k]]-1))*dim+j] * qs[k*(T+1)+t+pind];
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
    for(j=0; j<dim; j++){
      tempmbeta[i] += prec0[i*dim+j] * beta0[j];
    }
    for(k=0; k<supN; k++){
      for(t=(C[Hs[k]]-1); t<T; t++){
        tempmbeta[i] += DM[(cohortseq[C[Hs[k]]-1]+t-(C[Hs[k]]-1))*dim+i] * (u[k*T+t]-z[k])  * qs[k*(T+1)+t+pind];
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

void NEWZ(int i, int *q, int T, double *probitphi, double *zphi, int *C, int *L, int Hind, double *propz, int update)
{
  int t;
  int firstcap, lastcap;
  double phi,probz;
  
  propz[i]=0.;
  firstcap = C[Hind]-1;
  lastcap = L[Hind]-1;
  if(update){
    for(t=0; t<firstcap; t++){
      q[i*(T+1)+t] = 0;
    }
  }
  if(Hind){
    if(update) q[i*(T+1)+firstcap] = 1;
    if(firstcap<T){
      if(update) {
        for(t=firstcap; t<lastcap+1; t++){
          q[i*(T+1)+t] = 1;
        }
      }
      for(t=lastcap+1; t<T+1; t++){
        phi = INVPROBIT((probitphi[firstcap*T+t-1] + zphi[i]),0.0,1.0,1,0); 
        probz = q[i*(T+1)+t-1] * phi;
        if(update) q[i*(T+1)+t] = (int) rbinom(1.0,probz);
        propz[i] += dbinom((double) q[i*(T+1)+t],1.0,probz,1);
      }
    }
  }
}

void GETZ(int i, int *q, int T, double *probitp, double *probitphi, double *zp, double *zphi, int *C, int *L, int Hind, double *propz)
{
  int t;
  int firstcap, lastcap;
  double p,phi,num,probz;
  
  propz[i]=0.;
  firstcap = C[Hind]-1;
  lastcap = L[Hind]-1;
  for(t=0; t<firstcap; t++){
    q[i*(T+1)+t] = 0;
  }
  if(Hind){
    q[i*(T+1)+firstcap] = 1;
    if(firstcap<T){
      for(t=firstcap; t<lastcap+1; t++){
        q[i*(T+1)+t] = 1;
      }
      for(t=lastcap+1; t<T+1; t++){
        p = INVPROBIT((probitp[firstcap*T+t-1] + zp[i]),0.0,1.0,1,0);
        phi = INVPROBIT((probitphi[firstcap*T+t-1] + zphi[i]),0.0,1.0,1,0); 
        num = phi * q[i*(T+1)+t-1] * (1.-p);
        if(t<T) num *= (1.-INVPROBIT((probitphi[firstcap*T+t] + zphi[i]),0.0,1.0,1,0));
        probz = num / (num + (1.-phi*q[i*(T+1)+t-1]));
        if((q[i*(T+1)+t-1]) & (t<T)){
          if(q[i*(T+1)+t+1]) probz = 1.;
        }
        q[i*(T+1)+t] = (int) rbinom(1.0,probz);
        propz[i] += dbinom((double) q[i*(T+1)+t],1.0,probz,1);
      }
    }
  }
}

double FREQSUMCJS(int *x, int *Allhists, int T, int J, int ind, int *C, int useInitial)
{
  int j, t;
  int freqsum=0;
  for(j=0; j<J; j++){
    for(t=(C[j]-1); t<(T+useInitial); t++){
      freqsum += (Allhists[j*(T+1)+t+1-useInitial]==ind)*x[j];
    }
  }
  return(freqsum); 
}

int FIRSTCAP(int *x, int T, int J, int ind, int *C)
{
  int j;
  int freqsum=0;
  for(j=0; j<J; j++){
    freqsum += (((C[j]-1)==ind)*x[j]);
  }
  return(freqsum); 
}

/* Define function GETDELTA for updating delta by drawing from the full conditional posterior distribution */
void GETDELTACJS(double *deltavect, int *xs, int *Allhists, int T, int J, int dim, double *a0delta, int *C, int useInitial) 
{
  int k, kk;
  double nu[dim];
  
  for (kk=0; kk < dim; kk++)  {
    nu[kk]=FREQSUMCJS(xs,Allhists,T,J,kk+1,C,useInitial);
  }
  nu[dim-1]+=FREQSUMCJS(xs,Allhists,T,J,dim+1,C,useInitial);
  
  double xx[dim], sumx=0.0;
  for (k = 0; k < dim; k++) {
    xx[k]=rgamma((nu[k]+a0delta[k]),1.0);
    sumx += xx[k];
  }
  
  for (k = 0; k < dim; k++) {
    deltavect[k]=xx[k]/sumx;
  }
}

double GETprodhProbitCJS(int *Allhists, double *probitp, double *probitphi, double *zp, double *zphi, int *C, int *L, double delta_1, double delta_2, double alpha, int j,int T, int i, int useInitial)//, double psi)
{
  int t;
  int indhist;
  int firstcap = C[j]-1;
  int lastcap = L[j]-1;
  double p, phi;
  double alphavec[2];
  double alp0;
  double alp1;
  double alphavecsum;
  double alphasum=0.0;
  double Gamma[4];
  
  alphavec[0] = 1.;
  alphavec[1] = 0.;
  
  if(firstcap<(T+1)){
    if(useInitial){
      indhist = Allhists[j * (T+1) + firstcap];
      alphavec[0] = ( (indhist==1) * delta_1
                        + (indhist==2) * delta_2 
                        + (indhist==3) * (1.-delta_1-delta_2) * (1.-alpha)
                        + (indhist==4) * (1.-delta_1-delta_2) * alpha );
                        alphavecsum = alphavec[0]+alphavec[1];
                        alphasum += log(alphavecsum);
                        alphavec[0] = alphavec[0]/(alphavecsum);
                        alphavec[1] = alphavec[1]/(alphavecsum);
    }
    for(t=(firstcap+1); t<(lastcap+1); t++){
      indhist = Allhists[j * (T+1) + t];
      p = INVPROBIT((probitp[firstcap*T+t-1] + zp[i]),0.0,1.0,1,0);
      phi = INVPROBIT((probitphi[firstcap*T+t-1] + zphi[i]),0.0,1.0,1,0);
      alphavec[0] *= ( (indhist==0) * (1.-p) * phi
                         + (indhist==1) * p * delta_1 * phi
                         + (indhist==2) * p * delta_2 * phi
                         + (indhist==3) * p * (1.-delta_1-delta_2) * (1.-alpha) * phi
                         + (indhist==4) * p * (1.-delta_1-delta_2) * alpha * phi );
                         alphavecsum = alphavec[0]+alphavec[1];
                         alphasum += log(alphavecsum);
                         alphavec[0] = alphavec[0]/(alphavecsum);
                         alphavec[1] = alphavec[1]/(alphavecsum);
    }   
    for(t=(lastcap+1); t<(T+1); t++){
      p = INVPROBIT((probitp[firstcap*T+t-1] + zp[i]),0.0,1.0,1,0);
      phi = INVPROBIT((probitphi[firstcap*T+t-1] + zphi[i]),0.0,1.0,1,0);
      Gamma[0] = (1.-p) * phi;
      Gamma[1] = 0.0;
      Gamma[2] = (1.-phi);
      Gamma[3] = 1.0;
      alp0 = alphavec[0];
      alp1 = alphavec[1];
      alphavec[0] = alp0 * Gamma[0] + alp1 * Gamma[1];
      alphavec[1] = alp0 * Gamma[2] + alp1 * Gamma[3];
      alphavecsum = alphavec[0]+alphavec[1];
      alphasum += log(alphavecsum);
      alphavec[0] = alphavec[0]/alphavecsum;
      alphavec[1] = alphavec[1]/alphavecsum;
    }
  }
  double dens = exp(alphasum);
  dens = fmin(1.0, fmax(dens, tol));
  return dens;
}

//double LIKEProbitCJS(int *q, double *probitp, double *probitphi, double *zp, double *zphi, double delta_1, double delta_2, double alpha, double psi, int *Allhists, int *Hs, int T, int supN, int *C, int *L, int *x, int J, int *firstcapvec, int useInitial)
//{
//  int i,t,j;
//  double logdens=0.;
//  int n=0;
//  for(i=0; i<supN; i++)  {
//    n += (Hs[i] ? (int) 1 : (int) 0);
//    logdens += log(GETprodhProbitCJS(Allhists, probitp, probitphi, zp, zphi,
//                      C, L, delta_1, delta_2, alpha,
//                      Hs[i], T, i, useInitial));
//  }
//  //logdens += lgamma((double) supN + 1.) - lgamma((double) n + 1.) - lgamma((double) supN - n + 1.);
//  logdens += n * log(psi) + (supN - n) * log(1.0 - psi);
//  //if (!useInitial) {
//  //  for (t = 0; t < T; t++) {
//  //    logdens += lgamma((double) firstcapvec[t] + 1.0);
//  //  }
//  //}
//  //for(t=0; t<T+useInitial; t++){
//  //  logdens += lgamma((double) firstcapvec[t]+1.);
//  //}
//  //for(j=1; j<J; j++){
//  //  if(((C[j]-1)<T) | useInitial) logdens -= lgamma((double) x[j]+1.);
//  //}
//  //Rprintf("logdens %f pbeta %f phibeta %f alpha %f delta_1 %f delta_2 %f x0 %d \n",logdens,probitp[firstcap*T+t],probitphi[firstcap*T+t],alpha,delta_1,delta_2,x[0]);
//  return(logdens); 
//}

double LIKEProbitCJS(int *q, double *probitp, double *probitphi, double *zp, double *zphi, 
                     double delta_1, double delta_2, double alpha, 
                     int *Allhists, int *Hs, int T, int supN, int *C, int *x, int J, 
                     int *firstcapvec, int useInitial)
{
  int i, t;
  double logdens = 0.0;
  
  // 1. Calculate Individual Likelihoods
  for(i = 0; i < supN; i++) {
    int firstcap = C[Hs[i]] - 1;
    if(Hs[i] > 0 && firstcap < (T + 1)) {
      if(useInitial) {
        int indhist = Allhists[Hs[i] * (T + 1) + firstcap];
        logdens += log((indhist == 1) * delta_1 + (indhist == 2) * delta_2 + 
          (indhist == 3) * (1. - delta_1 - delta_2) * (1. - alpha) + 
          (indhist == 4) * (1. - delta_1 - delta_2) * alpha);
      }
      for(t = firstcap; t < T; t++) {
        if(q[i * (T + 1) + t]) {
          int indhist = Allhists[Hs[i] * (T + 1) + t + 1];
          double p = INVPROBIT((probitp[firstcap * T + t] + zp[i]), 0.0, 1.0, 1, 0);
          double phi = INVPROBIT((probitphi[firstcap * T + t] + zphi[i]), 0.0, 1.0, 1, 0);
          logdens += log((indhist == 0) * ((1. - p) * phi * q[i * (T + 1) + t + 1] + (1. - phi) * (1. - q[i * (T + 1) + t + 1])) + 
            (indhist == 1) * p * delta_1 * phi + (indhist == 2) * p * delta_2 * phi + 
            (indhist == 3) * p * (1. - delta_1 - delta_2) * (1. - alpha) * phi + 
            (indhist == 4) * p * (1. - delta_1 - delta_2) * alpha * phi);
        }
      }
    }
  }
  
  // 2. Branched Combinatorial Correction
  if (useInitial) {
    // Model: Random Type -> Stratify by Time AND Type
    int count[T + 1][5]; memset(count, 0, sizeof(count));
    for (i = 0; i < supN; i++) {
      if (Hs[i] > 0) {
        int tc = C[Hs[i]] - 1;
        int type = Allhists[Hs[i] * (T + 1) + tc];
        if(tc < T+1 && type >=1 && type <=4) count[tc][type]++;
      }
    }
    for (t = 0; t < T + 1; t++) {
      for (int m = 1; m <= 4; m++) {
        if (count[t][m] > 1) logdens -= lgamma((double)count[t][m] + 1.0);
      }
    }
  } else {
    // Model: Fixed Release -> Stratify by Time Only
    int count[T + 1]; memset(count, 0, sizeof(count));
    for (i = 0; i < supN; i++) {
      if (Hs[i] > 0) {
        int tc = C[Hs[i]] - 1;
        if(tc < T+1) count[tc]++;
      }
    }
    for (t = 0; t < T + 1; t++) {
      if (count[t] > 1) logdens -= lgamma((double)count[t] + 1.0);
    }
  }
  return(logdens); 
}

void PROPFREQProbitCJS(int icol, int c_k,
                       int *Hnew, int *indBasis, int J, int *xnew,
                       int supN, int T,
                       double *probitp, double *probitphi, double *zp, double *zphi,
                       int *C, int *L,
                       double delta_1, double delta_2, double alpha,
                       int *Allhists,
                       double *nprop, double *oprop,
                       int useInitial)
{
  const double mytol = tol;
  int remove_xi[3] = {0,0,0}, add_xi[3] = {0,0,0}, bInd[3];
  int absc_k = abs(c_k);
  if (absc_k == 0) return;
  
  if (c_k > 0) {
    remove_xi[0]=1; remove_xi[1]=1; remove_xi[2]=0;
    add_xi[0]=0;   add_xi[1]=0;    add_xi[2]=1;
    bInd[0]=0; bInd[1]=1; bInd[2]=2;
  } else {
    add_xi[0]=1; add_xi[1]=1; add_xi[2]=0;
    remove_xi[0]=0; remove_xi[1]=0; remove_xi[2]=1;
    bInd[0]=2; bInd[1]=1; bInd[2]=0;
  }
  
  int count = 0;
  double *prodz = (double *)calloc(supN, sizeof(double));
  double *prodh = (double *)calloc(supN, sizeof(double));
  int *sel_remove = (int *)calloc(absc_k, sizeof(int));
  int *sel_add = (int *)calloc(absc_k, sizeof(int));
  
  for (int jm = 0; jm < 3; jm++) {
    int jhist = indBasis[icol*3 + bInd[jm]];
    //int firstcap = C[jhist] - 1;
    //int type = Allhists[jhist * (T + 1) + firstcap];
    
    // 1. Compute n_old for the exchangeable pool (type-time cohort)
    //int n_old = 0;
    //for (int i = 0; i < supN; i++) {
    //  if (Hnew[i] > 0 && (C[Hnew[i]]-1) == firstcap && Allhists[Hnew[i]*(T+1)+firstcap] == type) n_old++;
    //}
    
    // 2. Compute individual weights
    for (int i = 0; i < supN; i++) {
      double h = GETprodhProbitCJS(Allhists, probitp, probitphi, zp, zphi, C, L, delta_1, delta_2, alpha, jhist, T, i, useInitial);
      prodz[i] = fmin(1.0, fmax(h, mytol));
      prodh[i] = fmin(1.0, fmax(1.0 - h, mytol));
    }
    
    double *wh_sample = (double *)calloc(supN, sizeof(double));
    double *wz_sample = (double *)calloc(supN, sizeof(double));
    double sum_rem_fwd = 0.0, sum_add_fwd = 0.0;
    
    for (int i = 0; i < supN; i++) {
      if (Hnew[i] == jhist) { wh_sample[i] = prodh[i]; sum_rem_fwd += prodh[i]; }
      if (Hnew[i] == 0)     { wz_sample[i] = prodz[i]; sum_add_fwd += prodz[i]; }
    }
    
    // 3. Selection and Standard Proposal Logs
    if (remove_xi[bInd[jm]]) {
      ProbSampleNoReplace(supN, wh_sample, absc_k, sel_remove);
      double cur_sum = sum_rem_fwd;
      for (int k = 0; k < absc_k; k++) {
        nprop[count + k] = log(prodh[sel_remove[k]]) - log(cur_sum);
        cur_sum -= prodh[sel_remove[k]];
      }
    } else {
      ProbSampleNoReplace(supN, wz_sample, absc_k, sel_add);
      double cur_sum = sum_add_fwd;
      for (int k = 0; k < absc_k; k++) {
        nprop[count + k] = log(prodz[sel_add[k]]) - log(cur_sum);
        cur_sum -= prodz[sel_add[k]];
      }
    }
    
    // --- COMBINATORIAL SYMMETRY CORRECTION ---
    // If removing individuals from the active pool, add log(choose(n_old, absc_k))
    //if (!useInitial && remove_xi[bInd[jm]]) {
    //  double log_choose = lgamma((double)n_old + 1.0) - lgamma((double)absc_k + 1.0) - lgamma((double)(n_old - absc_k) + 1.0);
    //  for (int k = 0; k < absc_k; k++) nprop[count + k] += (log_choose / (double)absc_k);
    //}
    
    // 4. Update Hnew
    for (int k = 0; k < absc_k; k++) {
      if (add_xi[bInd[jm]]) Hnew[sel_add[k]] = jhist;
      else Hnew[sel_remove[k]] = 0;
    }
    
    // 5. Compute Reverse Proposal
    //int n_new = 0;
    //for (int i = 0; i < supN; i++) {
    //  if (Hnew[i] > 0 && (C[Hnew[i]]-1) == firstcap && Allhists[Hnew[i]*(T+1)+firstcap] == type) n_new++;
    //}
    
    double sum_rem_rev = 0.0, sum_add_rev = 0.0;
    for (int i = 0; i < supN; i++) {
      if (Hnew[i] == jhist) sum_rem_rev += prodh[i];
      if (Hnew[i] == 0) sum_add_rev += prodz[i];
    }
    
    if (remove_xi[bInd[jm]]) { // Reverse is add
      double cur_sum = sum_add_rev;
      for (int k = 0; k < absc_k; k++) {
        oprop[count + k] = log(prodz[sel_remove[k]]) - log(cur_sum);
        cur_sum -= prodz[sel_remove[k]];
      }
    } else { // Reverse is remove
      double cur_sum = sum_rem_rev;
      for (int k = 0; k < absc_k; k++) {
        oprop[count + k] = log(prodh[sel_add[k]]) - log(cur_sum);
        cur_sum -= prodh[sel_add[k]];
      }
      // Add log(choose(n_new, absc_k)) for reverse remove move
      //if(!useInitial){
      //  double log_choose_rev = lgamma((double)n_new + 1.0) - lgamma((double)absc_k + 1.0) - lgamma((double)(n_new - absc_k) + 1.0);
      //  for (int k = 0; k < absc_k; k++) oprop[count + k] += (log_choose_rev / (double)absc_k);
      //}
    }
    
    count += absc_k;
    free(wh_sample); free(wz_sample);
  }
  
  free(prodz); free(prodh); free(sel_remove); free(sel_add);
  xnew[indBasis[icol*3+0]] -= c_k;
  xnew[indBasis[icol*3+1]] -= c_k;
  xnew[indBasis[icol*3+2]] += c_k;
  xnew[0] += c_k;
}

// Define function ProbitCJSC to draw samples from the posterior distribution

void ProbitCJSC(int *ichain, double *pbeta0, double *pprec0, double *pbeta, double *phibeta0, double *phiprec0, double *phibeta, double *zp, double *sigma2_zp, double *zphi, double *sigma2_zphi, double *delta_1, double *delta_2, double *alpha, int *x, double *psi, int *H, int *q,
                int *noccas, int *M, double *a0delta, double *a0alpha, double *b0alpha, double *l0p, double *d0p, double *l0phi, double *d0phi, double *a0psi, double *b0psi,
                double *loglike,
                int *nHists, int *Allhists, int *C, int *L, int *indBasis, int *ncolBasis, int *knownx, double *DMp, double *DMphi, int *pdim, int *phidim,
                int *iter, int *thin, int *numbasis,
                int *modp_h, int *modphi_h, int *data_type, int *zpind, int *zphiind, int *zind, int *Hind, int *updatedelta, int *delta_type, int *printlog, int *useInitial)
{
  GetRNGstate();
  
  int T = *noccas-1;
  double Ttmp = (double) T;
  int supN = *M;
  int Mk = supN*(T);
  int datatype = *data_type;
  int deltatype = *delta_type;
  
  int niter, th;
  int J = *nHists;
  
  niter = *iter;
  th = *thin;
  
  int dimp = pdim[1];
  int dimphi = phidim[1];
  
  double pbetas[dimp], phibetas[dimphi], zps[supN], zphis[supN], sigma2_zps, sigma2_zphis, alphas, delta_1s, delta_2s, psis;
  int xs[J], xnew[J], knownxs[J], Hs[supN], Hnew[supN];
  int ws[supN], wnew[supN];
  int firstcapvec[T+1], firstcapvecstar[T+1];
  
  double sha, sca;
  int base, nbasis;
  
  int t, i, j, k, g=0;
  double op, np;
  int dimrem, dimadd;
  
  for(k=0; k<dimp; k++){ pbetas[k]=pbeta[k]; }
  for(k=0; k<dimphi; k++){ phibetas[k]=phibeta[k]; }
  
  sigma2_zps= ( *modp_h ? sigma2_zp[0] : (double) 0.0);
  sigma2_zphis= ( *modphi_h ? sigma2_zphi[0] : (double) 0.0);
  double preczps= ( *modp_h ? 1./sigma2_zps : (double) 0.0);
  double preczphis= ( *modphi_h ? 1./sigma2_zphis : (double) 0.0);
  
  double zps2=0.;
  double zphis2=0.;
  delta_1s= (*updatedelta ? delta_1[0] : 1.);
  delta_2s= (*updatedelta ? delta_2[0] : 0.);
  alphas= (*updatedelta ? alpha[0] : 0.);
  psis= (*updatedelta ? psi[0] : 1.);
  
  double ns=0.;
  for(i=0; i< supN; i++)  {
    zps[i] = ( *modp_h ? zp[i] : (double) 0.0);
    zps2 += (zps[i]*zps[i]);
    zphis[i] = ( *modphi_h ? zphi[i] : (double) 0.0);
    zphis2 += (zphis[i]*zphis[i]);
    Hs[i] = H[i];
    //if(*useInitial) ws[i] = (Hs[i] ? (int) 1 : (int) 0);
    //else ws[i] = (((C[Hs[i]]-1)<T) ? (int) 1 : (int) 0);
    ws[i] = (Hs[i] ? (int) 1 : (int) 0);
    ns += (double) ws[i];
    wnew[i] = ws[i];
    Hnew[i] = Hs[i];
  }
  double nstar = ns;
  
  for(j=0; j<J; j++) {
    xs[j]=x[j];
    xnew[j]=xs[j];
    knownxs[j]= knownx[j];
  }
  int xind;
  
  double up[Mk], uphi[Mk];
  int qs[supN*(T+1)], qnew[supN*(T+1)];
  double propz[supN], newpropz[supN];
  
  for(i=0; i<supN; i++){
    for(t=0; t<(T+1); t++){
      qs[i*(T+1)+t] = q[i*(T+1)+t];
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
  
  GETTILDE(up,uphi,probitp,probitphi,zps,zphis,qs,T,supN,C,Hs,Allhists);
  
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
  for(i=0; i<3; i++){ a0_delta[i]=a0delta[i]; }
  deltavect[0]=delta_1s;
  deltavect[1]=delta_2s;
  deltavect[2]=1.-delta_1s-delta_2s;
  
  double indbase[*ncolBasis ? *ncolBasis : 1];
  for(j=0; j< *ncolBasis; j++){ indbase[j]=0.0; }
  int c_k;
  int ind, obasesum, nbasesum;
  
  for(t=0; t<(T+1); t++){
    firstcapvec[t] = FIRSTCAP(xs, T, J, t, C);
  }
  
  /* Calculate the log-likelihood */
  double ll=LIKEProbitCJS(qs,probitp,probitphi,zps,zphis,delta_1s,delta_2s,alphas,Allhists,Hs,T,supN,C,xs,J,firstcapvec,*useInitial);
  loglike[0]=ll;
  
  if(!R_FINITE(ll)) {
    Rprintf("FATAL ERROR in chain %d: initial likelihood is '%f'. \n",*ichain,ll);
    *iter = g;
    return;
  }
  double nl;
  
  /* Begin Markov chain */
  for (g=1; g < (niter+1); g++)  {
    
    /* Update pbetas */
    BETAUPDATE(pbetas,up,zps,DMp,qs,dimp,T,supN,C,Hs,pbeta0,pprec0s,cohortseq,1);
    
    /* Update phibetas */
    BETAUPDATE(phibetas,uphi,zphis,DMphi,qs,dimphi,T,supN,C,Hs,phibeta0,phiprec0s,cohortseq,0);
    
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
        double mbp = 0.0;
        Ttmp = 0.;
        for(t=(C[Hs[i]]-1); t<T; t++){
          mbp += (up[i*T+t] - probitp[(C[Hs[i]]-1)*T+t]) * qs[i*(T+1)+t+1];
          Ttmp += (double) qs[i*(T+1)+t+1];
        }
        double vbp = 1. / (preczps + Ttmp);
        zps[i]=rnorm((vbp*mbp),sqrt(vbp));
        zps2 += zps[i]*zps[i];
      }
      /* Update sigma2_zp */
      sha = *l0p + supN/2.0;
      sca = *d0p + zps2/2.0;
      preczps = rgamma((double) sha, (double) (1.0/sca));
      sigma2_zps = 1. / preczps;
    }
    
    if(*modphi_h){
      /* Update zphi */
      zphis2 = 0.0;
      for(i=0; i<supN; i++){
        double mbphi = 0.0;
        Ttmp = 0.;
        for(t=(C[Hs[i]]-1); t<T; t++){
          mbphi += (uphi[i*T+t] - probitphi[(C[Hs[i]]-1)*T+t]) * qs[i*(T+1)+t];
          Ttmp += (double) qs[i*(T+1)+t];
        }
        double vbphi = 1. / (preczphis + Ttmp);
        zphis[i]=rnorm((vbphi*mbphi),sqrt(vbphi));
        zphis2 += zphis[i]*zphis[i];
      }
      /* Update sigma2_zphi */
      sha = *l0phi + supN/2.0;
      sca = *d0phi + zphis2/2.0;
      preczphis = rgamma((double) sha, (double) (1.0/sca));
      sigma2_zphis = 1. / preczphis;
    }
    
    if(*updatedelta){
      /* update alpha */
      if(datatype){
        sha = *a0alpha+FREQSUMCJS(xs,Allhists,T,J,4,C,*useInitial);
        sca = *b0alpha+FREQSUMCJS(xs,Allhists,T,J,3,C,*useInitial);
        alphas = rbeta(sha,sca);
      }
      /* update delta_1 and delta_2 */
      if(deltatype){
        GETDELTACJS(deltavect, xs, Allhists, T, J, 3, a0_delta,C,*useInitial);
        delta_1s=deltavect[0];
        delta_2s=deltavect[1];
      } else {
        sha = a0_delta[0] + FREQSUMCJS(xs,Allhists,T,J,1,C,*useInitial) + FREQSUMCJS(xs,Allhists,T,J,2,C,*useInitial);
        sca = a0_delta[1] + FREQSUMCJS(xs,Allhists,T,J,3,C,*useInitial) + FREQSUMCJS(xs,Allhists,T,J,4,C,*useInitial);
        delta_1s = rbeta(sha,sca) / 2.0;
        delta_2s = delta_1s;
      }
      /* Update psi */
      sha = (double) *a0psi + ns;
      sca= (double) *b0psi + supN - ns;
      psis = rbeta(sha,sca);
    }
    
    ll=LIKEProbitCJS(qs,probitp,probitphi,zps,zphis,delta_1s,delta_2s,alphas,Allhists,Hs,T,supN,C,xs,J,firstcapvec,*useInitial);
    
    /* update q */
    for(i=0; i<supN; i++){
      GETZ(i,qs,T,probitp,probitphi,zps,zphis,C,L,Hs[i],propz);
      for(t=0; t<(T+1); t++){ qnew[i*(T+1)+t] = qs[i*(T+1)+t]; }
      newpropz[i]=propz[i];
    }
    
    ll=LIKEProbitCJS(qs,probitp,probitphi,zps,zphis,delta_1s,delta_2s,alphas,Allhists,Hs,T,supN,C,xs,J,firstcapvec,*useInitial);
    
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
        
        /* Forward proposal: basis selection */
        xind = min(xnew[indBasis[base*3+2]] - knownxs[indBasis[base*3+2]],xnew[0]);
        c_k=GETCK(xind+min(xnew[indBasis[base*3]]-knownxs[indBasis[base*3]],xnew[indBasis[base*3+1]]-knownxs[indBasis[base*3+1]]),xind) - xind;
        if(c_k>0){ dimadd=c_k; dimrem=2*c_k; }
        else { dimadd=2*(-c_k); dimrem= -c_k; }
        
        double nprop[dimadd+dimrem];
        double oprop[dimadd+dimrem];
        
        /* Forward proposal: c_k selection */
        op += -log((double) xind + min(xnew[indBasis[base*3]]-knownxs[indBasis[base*3]],xnew[indBasis[base*3+1]]-knownxs[indBasis[base*3+1]]));
        
        /* Call PROPFREQ */
        PROPFREQProbitCJS(base,c_k,Hnew,indBasis,J,xnew,supN,T,probitp,probitphi,zps,zphis,C,L,delta_1s,delta_2s,alphas,Allhists,nprop,oprop,*useInitial);
        
        /* Reverse proposal: c_k selection */
        xind = min(xnew[indBasis[base*3+2]] - knownxs[indBasis[base*3+2]],xnew[0]);
        np += -log((double) xind + min(xnew[indBasis[base*3]]-knownxs[indBasis[base*3]],xnew[indBasis[base*3+1]]-knownxs[indBasis[base*3+1]]));
        
        for(t=0; t<(T+1); t++){
          firstcapvecstar[t] = FIRSTCAP(xnew, T, J, t, C);
        }
        
        /* Accumulate PROPFREQ contributions */
        for(i=0; i<(dimrem+dimadd); i++){
          if(!R_FINITE(nprop[i])) {
            Rprintf("ERROR: PROPFREQ nprop[%d]=%f is not finite at iter %d, chain %d\n",i,nprop[i],g,*ichain); 
            *iter = g; 
            return;
          }
          if(!R_FINITE(oprop[i])) {
            Rprintf("ERROR: PROPFREQ oprop[%d]=%f is not finite at iter %d, chain %d\n",i,oprop[i],g,*ichain); 
            *iter = g; 
            return;
          }
          op += nprop[i];
          np += oprop[i];
        }
        
        /* Recompute basis pool */
        nbasesum=0;
        for(i=0; i< *ncolBasis; i++){
          if(min(xnew[indBasis[i*3+2]]-knownxs[indBasis[i*3+2]],xnew[0])+min(xnew[indBasis[i*3]]-knownxs[indBasis[i*3]],xnew[indBasis[i*3+1]]-knownxs[indBasis[i*3+1]])){
            indbase[i]=1.0;
            nbasesum+=1;
          } else {
            indbase[i]=0.0;
          }
        }
        
        /* Basis selection probabilities */
        op += -log((double) nbasesum);
        np += -log((double) obasesum);
        
        obasesum=nbasesum;
      }
    }
    
    if(ind){
      
      nstar=0.;
      for(i=0; i<supN; i++){
        //if(*useInitial) wnew[i] = (Hnew[i] ? (int) 1 : (int) 0);
        //else wnew[i] = (((C[Hnew[i]]-1)<T) ? (int) 1 : (int) 0);
        wnew[i] = (Hnew[i] ? (int) 1 : (int) 0);
        nstar += (double) wnew[i];
        
        if(Hnew[i]!=Hs[i]) {
          NEWZ(i,qnew,T,probitphi,zphis,C,L,Hnew[i],newpropz,1);
          NEWZ(i,qs,T,probitphi,zphis,C,L,Hs[i],propz,0);
          //GETZ(i,qnew,T,probitp,probitphi,zps,zphis,C,L,Hnew[i],newpropz);
          op += newpropz[i];
          np += propz[i];
        }
      }
      
      /* Attach likelihoods */
      op += ll;
      nl = LIKEProbitCJS(qnew,probitp,probitphi,zps,zphis,delta_1s,delta_2s,alphas,Allhists,Hnew,T,supN,C,xnew,J,firstcapvecstar,*useInitial);
      np += nl;
      
      if(runif(0.0,1.0)<exp(np-op)){
        //Rprintf("bam %d nl %f ll %f np %f op %f R %f newx0 %d oldx0 %d psi %f alpha %f delta_1 %f delta_2 %f phi %f p %f sigma2_zphi %f \n",g,nl,ll,np,op,exp(np-op),xnew[0],xs[0],psis,alphas,delta_1s,delta_2s,probitphi[0],probitp[0],sigma2_zphis);
        for(i=0; i<J; i++){ xs[i]=xnew[i]; }
        ns=nstar;
        for(i=0; i<supN; i++){
          if(Hnew[i]!=Hs[i]){
            for(t=0; t<(T+1); t++){
              qs[i*(T+1)+t] = qnew[i*(T+1)+t];
            }
            propz[i]=newpropz[i];
          }
          Hs[i]=Hnew[i];
          ws[i]=wnew[i];
        }
        for(t=0; t<(T+1); t++){
          firstcapvec[t] = firstcapvecstar[t];
        }
        ll=nl;
        
        /* update q */
        //for(i=0; i<supN; i++){
        //  GETZ(i,qs,T,probitp,probitphi,zps,zphis,C,L,Hs[i],propz);
        //  for(t=0; t<(T+1); t++){ qnew[i*(T+1)+t] = qs[i*(T+1)+t]; }
        //  newpropz[i]=propz[i];
        //}
      } 
    }
    
    for(i=0; i<J; i++){ xnew[i]=xs[i]; }
    for(i=0; i<supN; i++){ Hnew[i]=Hs[i]; wnew[i]=ws[i]; }
    nstar=ns;
    
    /* update up and uphi */
    GETTILDE(up,uphi,probitp,probitphi,zps,zphis,qs,T,supN,C,Hs,Allhists);
    
    /* Save draws according to thin specification */
    if((g % th)==0)  {
      
      for(j=0; j<(dimp); j++){ pbeta[(g/th - 1)*(dimp)+j]=pbetas[j]; }
      for(j=0; j<(dimphi); j++){ phibeta[(g/th - 1)*(dimphi)+j]=phibetas[j]; }
      sigma2_zp[(g/th - 1)]=sigma2_zps;
      sigma2_zphi[(g/th - 1)]=sigma2_zphis;
      delta_1[(g/th - 1)]=delta_1s;
      delta_2[(g/th - 1)]=delta_2s;
      alpha[(g/th - 1)]=alphas;
      psi[(g/th - 1)]=psis;
      
      if(*zpind){
        for(i=0; i<supN; i++){ zp[(g/th - 1)*supN+i]=zps[i]; }
      } else {
        for(i=0; i<supN; i++){ zp[i]=zps[i]; }
      }
      if(*zphiind){
        for(i=0; i<supN; i++){ zphi[(g/th - 1)*supN+i]=zphis[i]; }
      } else {
        for(i=0; i<supN; i++){ zphi[i]=zphis[i]; }
      }
      
      if(*zind){
        for(i=0; i<supN; i++){
          for(t=0; t<(T+1); t++){
            q[(g/th - 1)*supN*(T+1)+i*(T+1)+t]=qs[i*(T+1)+t];
          }
        }
      } else {
        for(i=0; i<supN; i++){
          for(t=0; t<(T+1); t++){
            q[i*(T+1)+t]=qs[i*(T+1)+t];
          }
        }
      }
      
      if(*Hind){
        for(i=0; i<supN; i++){ H[(g/th - 1)*supN+i]=Hs[i]; }
      } else {
        for(i=0; i<supN; i++){ H[i]=Hs[i]; }
      }
      
      for(j=0; j<J; j++){ x[j]=xs[j]; }
      
      loglike[(g/th - 1)]=ll;
      if(!R_FINITE(loglike[(g/th - 1)])) {
        Rprintf("FATAL ERROR in chain %d at iter %d: loglike not finite\n",*ichain,g); 
        *iter = g; 
        return;
      }
    }
    
    
    if(!(g%(niter/ min(niter,100)))){
      if(*printlog){
        Rprintf("Chain %d is %.0f%% complete\n",*ichain,(double) 100*g/niter);
      } else {
        Rprintf("\rChain %d is %.0f%% complete",*ichain,(double) 100*g/niter);
      }
    }
      
  } /* End Markov chain */
    
    PutRNGstate();
}

/* End function MCMCloop */
