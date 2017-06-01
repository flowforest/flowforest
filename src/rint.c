#include <R.h>
#include <Rinternals.h>
#include "tools.h"
#include "flowset.h"
#include "splitmerits.h"

//R interface for flowForest C part

#define EMERGE_R_FROM_R \
 GetRNGstate(); \
 uint32_t a=(uint32_t)(((double)(~((uint32_t)0)))*unif_rand()); \
 uint32_t b=(uint32_t)(((double)(~((uint32_t)0)))*unif_rand()); \
 PutRNGstate(); \
 rng_t rngdata; \
 rng_t *rng=&rngdata; \
 SETSEED(a,b)

static void tubesetFin(SEXP ptr){
 if(!R_ExternalPtrAddr(ptr)) return;
 tubeset *t=R_ExternalPtrAddr(ptr);
 destroyTubeset(t);
 free(t);
 R_ClearExternalPtr(ptr);
}

SEXP C_loadTubeset(SEXP fileNames){
 tubeset *t;
 int numTubes=length(fileNames);
 if(numTubes<1) error("Ough, no tubes!");
 t=malloc(sizeof(tubeset));
 t->K=numTubes;
 t->tubes=malloc((t->K)*sizeof(flowset));
 double avgEvents=0.;
 for(int e=0;e<numTubes;e++){
  const char *fN=CHAR(STRING_ELT(fileNames,e));
  Rprintf("Loading tube %d/%d...\n", e+1, numTubes);
  if(!readFlowMany(fN,&(t->tubes[e]))) error("Cannot load tube.");

  avgEvents+=t->tubes[e].avgEvents;
 }
 t->N=t->tubes[0].N;
 avgEvents/=(double)numTubes;
 SEXP ans,ptr;
 PROTECT(ans=allocVector(INTSXP,4+t->K));
 PROTECT(ptr=R_MakeExternalPtr(t,install("tubesetPtr"),R_NilValue));
 R_RegisterCFinalizerEx(ptr,tubesetFin,TRUE);

 int *Ans=INTEGER(ans);
 Ans[0]=t->K;
 if(t->K>0){
  Ans[1]=t->tubes[0].N;
  Ans[2]=t->tubes[0].M;
  Ans[3]=round(avgEvents);
  for(int e=0;e<t->K;e++)
   Ans[4+e]=t->tubes[e].sumS;
 }else{
  Ans[1]=Ans[2]=Ans[3]=NA_INTEGER;
 }

 setAttrib(ans,install("tubeset_ptr"),ptr);
 UNPROTECT(2);
 return(ans);
}

tubeset* r2tubeset(SEXP rTubeset){
 return(
  R_ExternalPtrAddr(
   getAttrib(rTubeset,install("tubeset_ptr"))));
}

gate* r2gate(SEXP rGate){
 double *c=REAL(rGate);
 int dim=(length(rGate)-1)/3;
 gate *ans=malloc(sizeof(gate));
 ans->tube=(uint)floor(c[0]);
 ans->dim=dim;
 ans->dims=malloc(sizeof(gated)*dim);
 for(int e=0;e<dim;e++){
  ans->dims[e].param=(uint)floor(c[1+e*3+0]);
  ans->dims[e].dir=(uint)floor(c[1+e*3+1]);
  ans->dims[e].thresh=c[1+e*3+2];
 }
 return(ans);
}

SEXP gate2r(gate *g){
 SEXP ans;
 PROTECT(ans=allocVector(REALSXP,1+3*g->dim));
 double *Ans=REAL(ans);
 Ans[0]=(double)g->tube;
 for(int e=0;e<g->dim;e++){
  Ans[1+e*3+0]=(double)(g->dims[e].param);
  Ans[1+e*3+1]=(double)(g->dims[e].dir);
  Ans[1+e*3+2]=(double)(g->dims[e].thresh);
 }
 UNPROTECT(1);
 return(ans);
}

SEXP C_extractEvents(SEXP rTubeset,SEXP rTube,SEXP rIdx){
 tubeset *t=r2tubeset(rTubeset);
 flowset *f=&(t->tubes[INTEGER(rTube)[0]-1]);

 int *idx=INTEGER(rIdx);
 int idxN=length(rIdx);

 SEXP ans;
 PROTECT(ans=allocVector(REALSXP,idxN*(f->M)));
 double *Ans=REAL(ans);

 for(int e=0;e<idxN;e++){
  float *source=flowsetEvent(f,idx[e]-1);
  for(int ee=0;ee<f->M;ee++)
   Ans[(f->M)*e+ee]=source[ee];
 }

 UNPROTECT(1);
 return(ans);
}

SEXP C_inGate(SEXP rTubeset,SEXP rGate,SEXP rIdx){
 tubeset *t=r2tubeset(rTubeset);
 gate *g=r2gate(rGate);

 int *idx=INTEGER(rIdx);
 int idxN=length(rIdx);

 SEXP ans;
 PROTECT(ans=allocVector(LGLSXP,idxN));
 int *Ans=LOGICAL(ans);

 inGate(g,t,idxN,idx,Ans);

 UNPROTECT(1);
 return(ans);
}

SEXP C_randomGate(SEXP rTubeset,SEXP dModeN,SEXP dParam){
 tubeset *t=r2tubeset(rTubeset);
 gate *g;
 int modeN=INTEGER(dModeN)[0];
 int param=INTEGER(dParam)[0];
 EMERGE_R_FROM_R;

 if(modeN==0){
  g=makeRandomGatePseudosphere(t,_R);
 }else{
  error("Bad randomGate creation mode.");
 }

 SEXP ans;
 PROTECT(ans=gate2r(g));
 destroyGate(g);
 free(g);
 UNPROTECT(1);
 return(ans);
}

SEXP C_dotGateMasked(SEXP rTubeset,SEXP rMask,SEXP rGate){
 tubeset *t=r2tubeset(rTubeset);
 gate *g=r2gate(rGate);

 uint N=length(rMask);

 SEXP ans;
 PROTECT(ans=allocVector(REALSXP,
  (t->tubes[g->tube].M+1)*(N)
 ));

 SEXP dimV;
 PROTECT(dimV=allocVector(INTSXP,2));
 INTEGER(dimV)[1]=N;
 INTEGER(dimV)[0]=t->tubes[g->tube].M+1;

 int *mask=INTEGER(rMask);
 double *Ans=REAL(ans);
 dotGateMasked(g,t,mask,N,Ans);
 destroyGate(g);
 free(g);

 setAttrib(ans,R_DimSymbol,dimV);

 UNPROTECT(2);
 return(ans);
}

//Registation

static const R_CallMethodDef callMethods[]={
 //herein
 {"C_loadTubeset",(DL_FUNC)&C_loadTubeset,1},
 {"C_extractEvents",(DL_FUNC)&C_extractEvents,3},
 {"C_inGate",(DL_FUNC)&C_inGate,3},
 {"C_randomGate",(DL_FUNC)&C_randomGate,3},
 {"C_dotGateMasked",(DL_FUNC)&C_dotGateMasked,3},
 //splitmerits
 {"C_colGI",(DL_FUNC)&C_colGI,2},
 {NULL,NULL,0}
};

void R_init_flowForest(DllInfo *dll){
 R_registerRoutines(dll,NULL,callMethods,NULL,NULL);
 R_useDynamicSymbols(dll,FALSE);
 R_forceSymbols(dll,TRUE);
}
