/* ranlib.h is free software, released under the GNU LGPL 
* Original FORTRAN77 version by Barry Brown, James Lovato; original C version by John Burkardt, Cleve Moler
* Modified by Brett McClintock (March 2015) for use in the 'multimark' R package.
*/

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <R.h>
# include <Rmath.h>
# include <Rinternals.h>

# include "ranlib.h"

static float sdot(long n,float *sx,long incx,float *sy,long incy)
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
        *(a+j-1+(j-1)*lda) = sqrtf(s);
    }
    *info = 0;
S40:
    return;
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
    *parm = (float) p;
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
    ranf = (float) ignlgi() * 4.656613057E-10f;
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
    0.0,3.917609E-2f,7.841241E-2f,0.11777f,0.1573107f,0.1970991f,0.2372021f,0.2776904f,
    0.3186394f,0.36013f,0.4022501f,0.4450965f,0.4887764f,0.5334097f,0.5791322f,
    0.626099f,0.6744898f,0.7245144f,0.7764218f,0.8305109f,0.8871466f,0.9467818f,
    1.00999f,1.077516f,1.150349f,1.229859f,1.318011f,1.417797f,1.534121f,1.67594f,
    1.862732f,2.153875f
};
static float d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843f,0.2425085f,0.2255674f,0.2116342f,0.1999243f,
    0.1899108f,0.1812252f,0.1736014f,0.1668419f,0.1607967f,0.1553497f,0.1504094f,
    0.1459026f,0.14177f,0.1379632f,0.1344418f,0.1311722f,0.128126f,0.1252791f,
    0.1226109f,0.1201036f,0.1177417f,0.1155119f,0.1134023f,0.1114027f,0.1095039f
};
static float t[31] = {
    7.673828E-4f,2.30687E-3f,3.860618E-3f,5.438454E-3f,7.0507E-3f,8.708396E-3f,
    1.042357E-2f,1.220953E-2f,1.408125E-2f,1.605579E-2f,1.81529E-2f,2.039573E-2f,
    2.281177E-2f,2.543407E-2f,2.830296E-2f,3.146822E-2f,3.499233E-2f,3.895483E-2f,
    4.345878E-2f,4.864035E-2f,5.468334E-2f,6.184222E-2f,7.047983E-2f,8.113195E-2f,
    9.462444E-2f,0.1123001f,0.136498f,0.1716886f,0.2276241f,0.330498f,0.5847031f
};
static float h[31] = {
    3.920617E-2f,3.932705E-2f,3.951E-2f,3.975703E-2f,4.007093E-2f,4.045533E-2f,
    4.091481E-2f,4.145507E-2f,4.208311E-2f,4.280748E-2f,4.363863E-2f,4.458932E-2f,
    4.567523E-2f,4.691571E-2f,4.833487E-2f,4.996298E-2f,5.183859E-2f,5.401138E-2f,
    5.654656E-2f,5.95313E-2f,6.308489E-2f,6.737503E-2f,7.264544E-2f,7.926471E-2f,
    8.781922E-2f,9.930398E-2f,0.11556f,0.1404344f,0.1836142f,0.2790016f,0.7010474f
};
static long i;
static float snorm,u,s,ustar,aa,w,y,tt;
    u = ranf();
    s = 0.0;
    if(u > 0.5) s = 1.0;
    u += (u-s);
    u = 32.0f*u;
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
    tt = (0.5f*w+aa)*w;
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
    u -= 1.0f;
S140:
    w = u**(d+i-1);
    tt = (0.5f*w+aa)*w;
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
    return(NA_INTEGER);
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
