#include "ranlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

double FREQSUM(int *x, int *Allhists, int T, int J, int ind)
{
  int j, t;
  int freqsum=0;
  for(j=0; j<J; j++){
    for(t=0; t<T; t++){
      freqsum += (Allhists[j*T+t]==ind)*x[j];
    }
  }
  return(freqsum); 
}

double DDIRICHLET(double *x, double *alpha, int dim)
{
  int i;
  double sumalpha=0.0;
  double sumx=0.0;
  double s=0.0;
  double logD=0.0;
  double logdens;
  int ind=1;
  for(i=0; i<dim; i++){
    sumx += x[i];
    sumalpha += alpha[i];
    s += (alpha[i]-1.)*log(x[i]);
    logD += lgamma(alpha[i]);
    if(x[i]<0 || x[i]>1) ind=0;
  }
  logD +=  - lgamma(sumalpha);
  if(fabs(sumx-1.)>1.e-7) ind=0;
  logdens=( ind ? (s-logD) : log(0.0));
  return(logdens);
}

void ProbSampleNoReplace(int n, double *prob, int nans, int *ans)
{
  double rT, mass, totalmass, pz[n], pzsum=0.0;
  int i, j, k, n1;
  int perm[n];
  
  /* Record element identities */
  for (i = 0; i < n; i++) {perm[i] = i + 1; pz[i]=prob[i];}
  
  /* Sort probabilities into descending order */
  /* Order element identities in parallel */
  revsort(pz, perm, n);
  
  /* set impermissible probs to zero */
  for (i = 0; i < n; i++) {if(pz[i]<0) pz[i]=0.; pzsum+=pz[i];}
  
  /* Compute the sample */
  totalmass = pzsum;
  for (i = 0, n1 = n-1; i < nans; i++, n1--) {
    rT = totalmass * unif_rand();
    mass = 0;
    for (j = 0; j < n1; j++) {
      mass += pz[j];
      if (rT <= mass) break;
    }
    ans[i] = perm[j] - 1;
    totalmass -= pz[j];
    for(k = j; k < n1; k++) {
      pz[k] = pz[k + 1];
      perm[k] = perm[k + 1];
    }
  }
}

/* Define function GETDELTA for updating delta by drawing from the full conditional posterior distribution */
void GETDELTA(double *deltavect, int *xs, int *Allhists, int T, int J, int dim, double *a0delta) 
{
  int k, kk;
  double nu[dim];
  
  for (kk=0; kk < dim; kk++)  {
    nu[kk]=FREQSUM(xs,Allhists,T,J,kk+1);
  }
  nu[dim-1]+=FREQSUM(xs,Allhists,T,J,dim+1);
  
  double xx[dim], sumx=0.0;
  for (k = 0; k < dim; k++) {
    xx[k]=rgamma((nu[k]+a0delta[k]),1.0);
    sumx += xx[k];
  }
  
  for (k = 0; k < dim; k++) {
    deltavect[k]=xx[k]/sumx;
  }
}

int GETCK(int populationSize, int nogo)
{
  int n = 1;
  int N = populationSize;
  int sample;
  
  int t = 0; 
  int m = 0; 
  double u;
  
  while (m < n)
  {
    u = runif(0.0,1.0); 
    
    if ( ((N - t)*u > n - m))
    {
      t++;
    }
    else
    {
      sample = (t>= nogo) ? t+1 : t;
      t++; m++;
    }
  }
  return(sample);
}

int sample(int n, double *prob)
{
  double rT, mass, totalmass, pz[n], pzsum=0.0;
  int i, j, k, n1;
  int perm[n];
  int nans = 1;
  int ans;
  
  /* Record element identities */
  for (i = 0; i < n; i++) {perm[i] = i + 1; pz[i]=prob[i]; pzsum+=pz[i];}
  
  /* Sort probabilities into descending order */
  /* Order element identities in parallel */
  revsort(pz, perm, n);
  
  /* Compute the sample */
  totalmass = pzsum;
  for (i = 0, n1 = n-1; i < nans; i++, n1--) {
    rT = totalmass * unif_rand();
    mass = 0;
    for (j = 0; j < n1; j++) {
      mass += pz[j];
      if (rT <= mass) break;
    }
    ans = perm[j] - 1;
    totalmass -= pz[j];
    for(k = j; k < n1; k++) {
      pz[k] = pz[k + 1];
      perm[k] = perm[k + 1];
    }
  }
  return(ans);
}
