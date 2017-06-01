#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>

#ifdef _WIN32
 #define FOR_WINDOWS 7
#endif
#ifdef _WIN64
 #define FOR_WINDOWS 7
#endif

#ifdef FOR_WINDOWS
 //Windows MapViewOfFile
 #include <Windows.h>
 #include <WinBase.h>
 #include <io.h>
#else
 //POSIX mmap
 #include <sys/mman.h>
 #include <err.h>
#endif

struct s_flowset{
 uint N; //flowframes in flowset
 uint M; //parameters
 uint *S; //rows in each flowframe
 uint sumS; //total rows in flowset, sum of S
 float *range; //ranges for each parameter in each flowframe
 float *quantile; //.05th and .95th quantiles in each flowframe
 float **q; //events
 double avgEvents; //average number of events in flowframe
 void *buf; //mmap'ed file
 size_t bufS;
};
typedef struct s_flowset flowset;

#define FLOWSET_EVENT(flowset,flowframe,index) ((flowset)->q[(flowframe)][(index)*((flowset)->M)])

float *flowsetEvent(flowset *fs,uint idx){
 uint iF=0;
 uint idxf;
 for(idxf=idx;idxf>fs->S[iF];idxf-=fs->S[iF++]);
 return(&FLOWSET_EVENT(fs,iF,idxf));
}

struct s_tubeset{
 uint K; //tubes
 uint N; //flowframes in flowset
 flowset *tubes;
};
typedef struct s_tubeset tubeset;

struct s_gated{
 uint param;
 uint dir;
 float thresh;
};
typedef struct s_gated gated;

struct s_gate{
 uint tube;
 uint dim;
 gated *dims;
};
typedef struct s_gate gate;

void dotGateMasked(gate *G,tubeset *ts,int *mask,int N,double *ans){
 gated *g=G->dims;
 uint gN=G->dim;
 flowset *fs=&(ts->tubes[G->tube]);
 uint P=1+(fs->M); //count+means
 for(uint eInMask=0;eInMask<N;eInMask++){
  uint eSet=mask[eInMask]-1;
  for(uint e=0;e<P;e++) ans[eInMask*P+e]=0.;
  for(uint eCell=0;eCell<(fs->S[eSet]);eCell++){
   float *cell=&FLOWSET_EVENT(fs,eSet,eCell);
   uint inGate=1;
   for(uint eGate=0;eGate<gN;eGate++)
    inGate=inGate&&
     ((cell[g[eGate].param]>g[eGate].thresh)==g[eGate].dir);
   if(inGate){
    ans[eInMask*P]++;
    for(uint e=0;e<(fs->M);e++)
     ans[eInMask*P+e+1]+=cell[e];
   }
  }
  for(uint e=0;e<(fs->M);e++)
   //                "Dirichlet prior"----------v
   ans[eInMask*P+e+1]=ans[eInMask*P+e+1]/(ans[eInMask*P]+1);
 }
}

//Fast check if a subset of flowset events is or is not in a given gate
void inGate(gate *G,tubeset *ts,int N,int  *restrict vIdx,int *restrict ans){
 gated *g=G->dims;
 uint gN=G->dim;
 flowset *fs=&(ts->tubes[G->tube]);
 for(int e=0;e<N;e++){
  float *cell=flowsetEvent(fs,vIdx[e]-1);
  uint inGate=1;
  for(uint eGate=0;eGate<gN;eGate++)
   inGate=inGate&&
    ((cell[g[eGate].param]>g[eGate].thresh)==g[eGate].dir);
  ans[e]=inGate;
 }
}

gate *makeRandomGatePseudosphere(tubeset *ts,R_){
 gate *G=malloc(sizeof(gate));
 G->tube=RINDEX(ts->K);
 flowset *fs=&(ts->tubes[G->tube]);
 G->dim=(fs->M)*2;
 gated *g=malloc((G->dim)*sizeof(gated));
 G->dims=g;
 for(uint e=0;e<fs->M;e++){
    g[2*e].param=e; g[2*e].dir=1; //up
  g[2*e+1].param=e; g[2*e+1].dir=0; //down
 }

 uint setI=RINDEX(fs->N);
 float *qs=fs->q[setI];
 uint cellO=RINDEX(fs->S[setI]);

 for(uint e=0;e<fs->M;e++){
  float range=fs->quantile[2*e+1]-fs->quantile[2*e];
  float S=qs[cellO*(fs->M)+e];
  float size;
    size=MAX(RUNIF_CLOSED*2., .1);
    g[2*e].thresh=S-range*size;
  g[2*e+1].thresh=S+range*size;
 }
 return(G);
}

void destroyGate(gate *G){
 if(G->dims) free(G->dims);
}

int readFlowMany(const char *fName,flowset *fs){
 fs->S=NULL;
 fs->q=NULL;
 fs->range=NULL;
 fs->quantile=NULL;

 struct stat st;
 stat(fName,&st);
 fs->bufS=st.st_size;

 int fd=-1;
 fd=open(fName,O_RDWR,0);
 if(fd==-1){
   error("Error while opening %s, exiting.\n",fName);
   return(-1);
 }

 #ifdef FOR_WINDOWS
  HANDLE h=(HANDLE)_get_osfhandle(fd);
  HANDLE fm=CreateFileMapping(h,NULL,PAGE_READONLY,0,0,NULL);
  if(!fm){
   close(fd);
   return(-1);
  }
  char *buf=MapViewOfFile(fm,FILE_MAP_READ,0,0,fs->bufS);
  CloseHandle(fm);
 #else
  char *buf=(char*)mmap(NULL,
   fs->bufS,
   PROT_READ,
   MAP_PRIVATE,
   fd,0);
 #endif

 close(fd);
 if(!buf)
  return(-1);
 fs->buf=buf;

 fs->N=((uint*)buf)[0]; buf+=sizeof(uint);
 fs->M=((uint*)buf)[0]; buf+=sizeof(uint);
 fs->S=malloc(fs->N*sizeof(uint));
 fs->sumS=0;
 fs->q=malloc(fs->N*sizeof(float*));
 for(int e=0;e<fs->N;e++) fs->q[e]=NULL;
 fs->range=malloc(sizeof(float)*2*(fs->M));
 fs->quantile=malloc(sizeof(float)*2*(fs->M));
 fs->avgEvents=0.;
 for(int e=0;e<fs->N;e++){
  fs->S[e]=((uint*)buf)[0]; buf+=sizeof(uint);
  fs->sumS+=fs->S[e];
  fs->avgEvents+=fs->S[e];
  fs->q[e]=((float*)buf); buf+=(fs->S[e])*(fs->M)*sizeof(float);
 }
 
 //Load ranges and quantiles from inside tube files
 for(uint e=0;e<fs->M;e++){
  fs->range[2*e]=((float*)buf)[e];
  fs->range[2*e+1]=((float*)buf)[e+(fs->M)];
  fs->quantile[2*e]=((float*)buf)[e+2*(fs->M)];
  fs->quantile[2*e+1]=((float*)buf)[e+3*(fs->M)];
 }
 buf+=sizeof(float)*4*(fs->M);
 Rprintf("\n");
 for(int e=0;e<fs->M;e++){
    Rprintf("Parameter %d. Range: %0.3f--%0.3f; ",e+1,fs->range[2*e],fs->range[2*e+1]);
    Rprintf(".05th quantile: %0.3f; ",fs->quantile[2*e]);
    Rprintf(".95th quantile: %0.3f;\n",fs->quantile[2*e+1]);
 }
 if(((uint*)buf)[0]!=1337){
  return(0);
 }
 fs->avgEvents/=((double)fs->N);
 return(1);
}

void destroyFlowset(flowset *fs){
 if(fs->buf){
  //This flowset is directly made by readFlowMany,
  // so it has a memory map of a bin file which we must release.
  #ifdef FOR_WINDOWS
   UnmapViewOfFile(fs->buf);
  #else
   munmap(fs->buf,fs->bufS);
  #endif
 }else{
  //This flowset just points to an other flowset
 }
 if(fs->q) free(fs->q);
 if(fs->S) free(fs->S);
 if(fs->range) free(fs->range);
 if(fs->quantile) free(fs->quantile);
}

void destroyTubeset(tubeset *ts){
 for(uint e=0;e<ts->K;e++)
  destroyFlowset(&(ts->tubes[e]));
 if(ts->tubes) free(ts->tubes);
}
