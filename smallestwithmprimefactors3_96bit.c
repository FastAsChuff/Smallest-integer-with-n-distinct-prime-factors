//===================================================================================================================================
// smallestwithmprimefactors3.c: Searches for smallest integer in [a,b] with at least n distinct prime factors.
//===================================================================================================================================
// Author: Simon Goater Dec 2025
// Motivated by question https://math.stackexchange.com/questions/2430594/what-is-the-smallest-20-digit-number-with-at-least-13-distinct-prime-factors/5112563#5112563
// It's n distinct prime factors in the code, m in the Stack Exchange answer and program name, just to confuse.
// This code is given as is and is not rigorously tested. No liability for consequences of its use is recognised.
// 
// COPYRIGHT NOTICE: Copying, modifying, and distributing with conspicuous attribution for any purpose is permitted.
//===================================================================================================================================
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <gmp.h>
#include "/home/simon/Mairsonsprimesieve.c" // https://github.com/FastAsChuff/Primes-List/blob/main/Mairsonsprimesieve.c

// gcc smallestwithmprimefactors3_96bit.c -o smallestwithmprimefactors3_96bit.bin -O3 -Wall -mssse3 -lm -lpthread -lgmp

#define U128 unsigned __int128

uint64_t isqrtu128(U128 n) {
  if (n < 2) return n;
  U128 ai = sqrt(n);
  while (!((ai <= n/ai) && ((ai+1) > n/(ai+1)))) {    
    ai = (ai + n/ai)/2;
  }
  return ai;
}

_Bool isprimeu128(U128 n) {
  mpz_t mpzn;
  mpz_init(mpzn);
  mpz_set_ui(mpzn, (uint64_t)(n >> 64));
  mpz_mul_2exp(mpzn, mpzn, 64);
  mpz_add_ui(mpzn, mpzn, (uint64_t)(n & 0xffffffffffffffffULL));
  int res = mpz_probab_prime_p(mpzn,15);
  mpz_clear(mpzn);
  return res >= 1;
}

uint64_t icbrtu128(U128 n) {
  if (n < 2) return n;
  U128 ai = cbrt(n);
  while (!((ai <= n/(ai*ai)) && ((ai+1) > n/((ai+1)*(ai+1))))) {    
    ai = (2*ai + n/(ai*ai))/3;
  }
  return ai;
}

U128 atou128(char *in) {
  U128 res = 0;
  while (*in) {
    res *= 10;
    res += *in - '0';
    in++;
  }
  return res;
}

#ifndef MIN
  #define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif

void getcounts(U128 *remainingfactors, uint8_t *primefactorcounts, uint32_t primefactorcountssize, U128 a, uint32_t *primes, uint32_t numprimes, uint32_t n) {
  memset(primefactorcounts, 0, primefactorcountssize);
  // Only counts >= n guaranteed to be correct.
  for (uint32_t x=0; x<primefactorcountssize; x++) remainingfactors[x] = x+a;
  for (uint32_t primeix=0; primeix<numprimes; primeix++) {
    uint32_t p = primes[primeix];
    uint64_t x = a % p;
    if (x) x = p-x;
    while (true) {
      if (x >= primefactorcountssize) break;
      primefactorcounts[x]++;
      remainingfactors[x] /= p;
      while (0 == (remainingfactors[x] % p)) remainingfactors[x] /= p;
      x += p;
      if (x < p) break;
    }
  }
  for (uint32_t j=0; j<primefactorcountssize; j++) {
    if (1 != remainingfactors[j]) {
      if (2+primefactorcounts[j] >= n) {
        uint64_t q = isqrtu128(remainingfactors[j]);
        if ((U128)q*q == remainingfactors[j]) {
          primefactorcounts[j]++;
        } else if (isprimeu128(remainingfactors[j])) {
          primefactorcounts[j]++;
        } else {          
          primefactorcounts[j] += 2;
        }
      }
    }
  }
}
pthread_mutex_t intervalstartmutex = PTHREAD_MUTEX_INITIALIZER;
U128 result = 0;
U128 intervalstart = 0;

typedef struct {
  U128 b;
  uint32_t *primes;
  uint32_t numprimes;
  uint32_t n;
  uint32_t primefactorcountssize;
} searchintervalargs_s;

U128 nextintervalstart(U128 intervalstart, U128 b, uint32_t *primefactorcountssize) {
  if (intervalstart + *primefactorcountssize >= *primefactorcountssize) {
    if (intervalstart + *primefactorcountssize <= b) {
      return intervalstart + *primefactorcountssize;
    } else {
      *primefactorcountssize = 1+(b-intervalstart);
      return 0;
    }
  } else {
    *primefactorcountssize = 1+(b-intervalstart);
    return 0;
  }
}

void *searchinterval(void *args) {
  searchintervalargs_s *searchintervalargs = (searchintervalargs_s *)args;
  U128 b = searchintervalargs->b;
  uint32_t *primes = searchintervalargs->primes;
  uint32_t numprimes = searchintervalargs->numprimes;
  uint32_t n = searchintervalargs->n;
  uint32_t primefactorcountssize = searchintervalargs->primefactorcountssize;  
  U128 i;
  _Bool finished = false;
  pthread_mutex_lock(&intervalstartmutex);
  if (result || (0 == intervalstart)) {
    finished = true;
  } else {
    i = intervalstart;
    intervalstart = nextintervalstart(intervalstart, b, &primefactorcountssize);
  }
  pthread_mutex_unlock(&intervalstartmutex);
  uint8_t *primefactorcounts = calloc(primefactorcountssize,sizeof(uint8_t));
  assert(primefactorcounts);
  U128 *remainingfactors = malloc(primefactorcountssize*sizeof(U128));
  assert(remainingfactors);
  uint32_t j = 0;
  while (!finished) {
    //printf("[%lu, %lu]\n", i, i+primefactorcountssize-1);
    getcounts(remainingfactors, primefactorcounts, primefactorcountssize, i, primes, numprimes, n);
    for (j=0; j<primefactorcountssize; j++) {
      //printf("%u ", primefactorcounts[j]);
      if (primefactorcounts[j] >= n) {
        pthread_mutex_lock(&intervalstartmutex);
        if ((result == 0) || (i+j < result)) result = i+j;
        pthread_mutex_unlock(&intervalstartmutex);
        finished = true;
        break;
      }
    }
    pthread_mutex_lock(&intervalstartmutex);
    if (result || (0 == intervalstart)) {
      finished = true;
    } else {
      i = intervalstart;
      intervalstart = nextintervalstart(intervalstart, b, &primefactorcountssize);
    }
    pthread_mutex_unlock(&intervalstartmutex);
  }
  free(remainingfactors);
  free(primefactorcounts);
  return NULL;
}

int main(int argc, char* argv[]) {
  _Bool validinput = true;
  U128 a,b;
  uint32_t n;
  uint16_t numthreads;
  uint8_t maxprimorial = 20;
  if (argc < 5) {
    validinput = false;
  } else {
    numthreads = atoi(argv[1]);
    a = atou128(argv[2]);
    b = atou128(argv[3]);
    n = atoi(argv[4]);
    if (numthreads < 1) validinput = false;
    if (n < 1) validinput = false;
    if (n > maxprimorial) validinput = false;
    if (a < 1) validinput = false;
    if (a > b) validinput = false;
    if (a >> 96) validinput = false;
    if (b >> 96) validinput = false;
  }
  if (!validinput) {
    printf("This program searches for the smallest integer x with a <= x <= b with x having at least n distinct prime factors.\nAuthor: Simon Goater Dec 2025\nUsage: %s numthreads a b n \n1 <= a <= b < 2^96\n1 <= n <= %u\n", argv[0], maxprimorial);
    exit(0);
  }
  uint32_t primefactorcountssize = 1000000;
  primefactorcountssize = MIN(primefactorcountssize, 1+(b-a));
  uint32_t numprimes;
  uint32_t icbrtb = icbrtu128(b);
  uint32_t smallprimes[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71};
  assert(maxprimorial <= (sizeof(smallprimes)/sizeof(uint32_t)));
  U128 nthprimorial = 1;
  for (uint32_t i = 0; i<n; i++) nthprimorial *= smallprimes[i];
  if (b < nthprimorial) {
    printf("There are no integers in [%lu%018lu, %lu%018lu] with %u or more distinct prime factors.\n", (uint64_t)(a / 1000000000000000000ULL), (uint64_t)(a % 1000000000000000000ULL),(uint64_t)(b / 1000000000000000000ULL), (uint64_t)(b % 1000000000000000000ULL),n);
    exit(0);
  }
  if (a <= nthprimorial) {
    printf("%lu%018lu is the smallest integer in [%lu%018lu, %lu%018lu] with at least %u distinct prime factors.\n",(uint64_t)(nthprimorial / 1000000000000000000ULL), (uint64_t)(nthprimorial % 1000000000000000000ULL), (uint64_t)(a / 1000000000000000000ULL), (uint64_t)(a % 1000000000000000000ULL),(uint64_t)(b / 1000000000000000000ULL), (uint64_t)(b % 1000000000000000000ULL), n);
    exit(0);
  }
  uint32_t *primes = Mairsonsprimesieve(icbrtb, &numprimes);
  assert(primes);
  intervalstart = a;
  pthread_t threads[numthreads];
  searchintervalargs_s searchintervalargs[numthreads];
  for (uint32_t i=0; i<numthreads; i++) {
    searchintervalargs[i].b = b;
    searchintervalargs[i].n = n;
    searchintervalargs[i].primes = primes;
    searchintervalargs[i].numprimes = numprimes;
    searchintervalargs[i].primefactorcountssize = primefactorcountssize;
    pthread_create(&threads[i], NULL, searchinterval, &searchintervalargs[i]);
  }
  for (uint32_t i=0; i<numthreads; i++) pthread_join(threads[i], NULL);
  if (result) {
    printf("%lu%018lu is the smallest integer in [%lu%018lu, %lu%018lu] with at least %u distinct prime factors.\n",(uint64_t)(result / 1000000000000000000ULL), (uint64_t)(result % 1000000000000000000ULL), (uint64_t)(a / 1000000000000000000ULL), (uint64_t)(a % 1000000000000000000ULL), (uint64_t)(b / 1000000000000000000ULL), (uint64_t)(b % 1000000000000000000ULL), n);
  } else {
    printf("There are no integers in [%lu%018lu, %lu%018lu] with %u or more distinct prime factors.\n", (uint64_t)(a / 1000000000000000000ULL), (uint64_t)(a % 1000000000000000000ULL), (uint64_t)(b / 1000000000000000000ULL), (uint64_t)(b % 1000000000000000000ULL),n);
  }
  free(primes);
}
/*
1000000000000000014921840 is the smallest integer in [1000000000000000000000000, 1100000000000000000000000] with at least 13 distinct prime factors.
1000000000000000393986060 is the smallest integer in [1000000000000000000000000, 1100000000000000000000000] with at least 14 distinct prime factors.
1000000000000000000103482380 is the smallest integer in [1000000000000000000000000000, 1100000000000000000000000000] with at least 14 distinct prime factors.
1000000000000000001659399050 is the smallest integer in [1000000000000000000000000000, 1100000000000000000000000000] with at least 15 distinct prime factors.
*/
