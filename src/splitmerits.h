#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#define sqr(x) ((x)*(x))

//Gini index
void gi(int N,double *score,int *dec,double *minGI,double *bestThre){
  if(N<2){
  minGI[0]=INFINITY;
  bestThre[0]=score[0];//disfunctional split
  return;
 }
 int totT=0;
 double *x=(double*)R_alloc(N,sizeof(double));
 for(int e=0;e<N;e++) x[e]=score[e];
 int *y=(int*)R_alloc(N,sizeof(int));
 for(int e=0;e<N;e++){
  y[e]=dec[e];
  if(y[e]) totT++;
 }
 int totF=N-totT;

 rsort_with_index(x,y,N);
 //x is now sorted in increasing order
 // and y permuted the same way x was.

 int befT=0;
 int befF=0;

 minGI[0]=INFINITY;
 bestThre[0]=x[0];//disfunctional split

 for(int e=0;e<N-1;e++){
  befT+=y[e];
  befF+=!(y[e]);
  int aftT=totT-befT;
  int aftF=totF-befF;
  int befN=e+1;
  int aftN=N-e-1;
  if(x[e+1]!=x[e]){
   //GI calculation.
   //Warning: this may never happen
   double gi=((double)befN)/((double)N)*(
    1-
    sqr(((double)befT)/((double)befN))
    -
    sqr(((double)befF)/((double)befN)))
    +((double)aftN)/((double)N)*(
    1-
    sqr(((double)aftT)/((double)aftN))
    -
    sqr(((double)aftF)/((double)aftN)));
    if(gi<minGI[0]){
     minGI[0]=gi;
     bestThre[0]=(x[e]+x[e+1])/2.;
    }
  }
 }
 return;
}

SEXP C_colGI(SEXP scores,SEXP y){
 int N=length(y);
 int nCases=length(scores)/N;

 SEXP ans;
 PROTECT(ans=allocVector(REALSXP,2*nCases));
 SEXP ansDim;
 PROTECT(ansDim=allocVector(INTSXP,2));
 INTEGER(ansDim)[0]=2;
 INTEGER(ansDim)[1]=nCases;
 setAttrib(ans,R_DimSymbol,ansDim);
 double *theAns=REAL(ans);
 double *theScores=REAL(scores);
 int *theY=INTEGER(y);

 for(int e=0;e<nCases;e++)
  gi(N,theScores+(N*e),theY,theAns+(2*e),theAns+(2*e+1));


 UNPROTECT(2);
 return(ans);
}
