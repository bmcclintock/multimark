#include "ranlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

double FREQSUM(int *x, int *Allhists, int T, int J, int ind);
double DDIRICHLET(double *x, double *alpha, int dim);
void ProbSampleNoReplace(int n, double *prob, int nans, int *ans);
void GETDELTA(double *deltavect, int *xs, int *Allhists, int T, int J, int dim, double *a0delta); 
int GETCK(int populationSize, int nogo);
int sample(int n, double *prob);
