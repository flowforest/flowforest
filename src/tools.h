/*   Shared C code

  Copyright 2014--2016 Miron B. Kursa
*/

#include <stdint.h>

typedef uint32_t uint;
typedef uint32_t mask;
typedef int32_t sint;
typedef double score_t;

//Internal PRNG -- based on G. Marsaglia's not-better-idea generator
struct rng{
 uint32_t z;
 uint32_t w;
};

typedef struct rng rng_t;
#define znew (((rng->z)=36969*((rng->z)&65535)+((rng->z)>>16))<<16)
#define wnew (((rng->w)=18000*((rng->w)&65535)+((rng->w)>>16))&65535)
#define RINTEGER_MAX (~((uint32_t)0))
#define RINTEGER (znew+wnew)
//Gives a number from [0,1]
#define RUNIF_CLOSED (((double)RINTEGER)*(1./4294967295.))
//Gives a number from [0,1)
#define RUNIF_OPEN (((double)RINTEGER)*(1./4294967296.))
//Gives a number from 0 to upTo-1 (each value equally probable provided upTo is safely smaller than 1<<32)
#define RINDEX(upTo) ((uint32_t)(RUNIF_OPEN*((double)upTo)))
//Setting seed; give it two uint32s you like
#define SETSEED(a,b) rng->z=a; rng->w=b
#define R_ rng_t *rng
#define _R rng

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MIN_3(a,b,c) MIN(MIN(a,b),c)
#define MAX_3(a,b,c) MAX(MAX(a,b),c)
#define MIN_4(a,b,c,d) MIN(MIN_3(a,b,c),d)
#define MAX_4(a,b,c,d) MAX(MAX_3(a,b,c),d)
#define MIN_5(a,b,c,d,e) MIN(MIN_4(a,b,c,d),e)
#define MAX_5(a,b,c,d,e) MAX(MAX_4(a,b,c,d),e)
#define MIN_6(a,b,c,d,e,f) MIN(MIN_5(a,b,c,d,e),f)
#define MAX_6(a,b,c,d,e,f) MAX(MAX_5(a,b,c,d,e),f)
#define MIN_7(a,b,c,d,e,f,g) MIN(MIN_6(a,b,c,d,e,f),g)
#define MAX_7(a,b,c,d,e,f,g) MAX(MAX_6(a,b,c,d,e,f),g)


int floatCompare (const void * a, const void * b)
{
  float fa = *(const float*) a;
  float fb = *(const float*) b;
  return (fa > fb) - (fa < fb);
}

float MIN_ARRAY(const float *arr, size_t length) {
    // returns the minimum value of array
    size_t i;
    float minimum = arr[0];
    for (i = 1; i < length; ++i) {
        if (minimum > arr[i]) {
            minimum = arr[i];
        }
    }
    return minimum;
}

float MAX_ARRAY(const float *arr, size_t length) {
    // returns the maximum value of array
    size_t i;
    float maximum = arr[0];
    for (i = 1; i < length; ++i) {
        if (maximum < arr[i]) {
            maximum = arr[i];
        }
    }
    return maximum;
}
