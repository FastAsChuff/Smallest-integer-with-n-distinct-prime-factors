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

// gcc smallestwithmprimefactors3.c -o smallestwithmprimefactors3.bin -O3 -Wall -mssse3 -lm -lpthread -lgmp
 
uint32_t isqrtu64(uint64_t n) {
  if (n < 2) return n;
  uint64_t ai = sqrt(n);
  while (!((ai <= n/ai) && ((ai+1) > n/(ai+1)))) {    
    ai = (ai + n/ai)/2;
  }
  return ai;
}

_Bool isprimeu64(uint64_t n) {
  mpz_t mpzn;
  mpz_init(mpzn);
  mpz_set_ui(mpzn, n);
  int res = mpz_probab_prime_p(mpzn,15);
  mpz_clear(mpzn);
  return res >= 1;
}

uint32_t icbrtu64(uint64_t n) {
  if (n < 2) return n;
  uint64_t ai = cbrt(n);
  while (!((ai <= n/(ai*ai)) && ((ai+1) > n/((ai+1)*(ai+1))))) {    
    ai = (2*ai + n/(ai*ai))/3;
  }
  return ai;
}

uint64_t atou64(char *in) {
  uint64_t res = 0;
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

void getcount(uint8_t *primefactorcounts, uint32_t primefactorcountssize, uint64_t a, uint32_t *primes, uint32_t numprimes) {
  memset(primefactorcounts, 0, primefactorcountssize);
  for (uint32_t x=0; x<primefactorcountssize; x++) {
    uint64_t n=x+a;
    for (uint32_t primeix=0; primeix<numprimes; primeix++) {
      if ((n % primes[primeix]) == 0) {
        primefactorcounts[x]++;
        n /= primes[primeix];
        while (((n % primes[primeix]) == 0)) n /= primes[primeix];
      } 
      if (n == 1) break;
    }
    if (n > 1) primefactorcounts[x]++;
  }
}

void getcounts(uint64_t *remainingfactors, uint8_t *primefactorcounts, uint32_t primefactorcountssize, uint64_t a, uint32_t *primes, uint32_t numprimes, uint32_t n) {
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
        uint32_t q = isqrtu64(remainingfactors[j]);
        if ((uint64_t)q*q == remainingfactors[j]) {
          primefactorcounts[j]++;
        } else if (isprimeu64(remainingfactors[j])) {
          primefactorcounts[j]++;
        } else {          
          primefactorcounts[j] += 2;
        }
      }
    }
  }
}
pthread_mutex_t intervalstartmutex = PTHREAD_MUTEX_INITIALIZER;
uint64_t result = 0;
uint64_t intervalstart = 0;

typedef struct {
  uint64_t a,b;
  uint32_t *primes;
  uint32_t numprimes;
  uint32_t n;
  uint32_t primefactorcountssize;
} searchintervalargs_s;

uint64_t nextintervalstart(uint64_t intervalstart, uint64_t b, uint32_t *primefactorcountssize) {
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
  uint64_t a = searchintervalargs->a;
  uint64_t b = searchintervalargs->b;
  uint32_t *primes = searchintervalargs->primes;
  uint32_t numprimes = searchintervalargs->numprimes;
  uint32_t n = searchintervalargs->n;
  uint32_t primefactorcountssize = searchintervalargs->primefactorcountssize;  
  uint64_t i;
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
  uint64_t *remainingfactors = malloc(primefactorcountssize*sizeof(uint64_t));
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
  uint64_t a,b;
  uint32_t n;
  uint16_t numthreads;
  if (argc < 5) {
    validinput = false;
  } else {
    numthreads = atoi(argv[1]);
    a = atou64(argv[2]);
    b = atou64(argv[3]);
    n = atoi(argv[4]);
    if (numthreads < 1) validinput = false;
    if (n < 1) validinput = false;
    if (n > 30) validinput = false;
    if (a < 1) validinput = false;
    if (a > b) validinput = false;
  }
  if (!validinput) {
    printf("This program searches for the smallest integer x with a <= x <= b with x having at least n distinct prime factors.\nAuthor: Simon Goater Dec 2025\nUsage: %s numthreads a b n \n1 <= a <= b < 2^64\n1 <= n <= 30\n", argv[0]);
    exit(0);
  }
  uint32_t primefactorcountssize = 100000;
  primefactorcountssize = MIN(primefactorcountssize, 1+(b-a));
  uint32_t numprimes;
  uint32_t icbrtb = icbrtu64(b);
  uint32_t *primes = Mairsonsprimesieve(icbrtb, &numprimes);
  assert(primes);
  intervalstart = a;
  pthread_t threads[numthreads];
  searchintervalargs_s searchintervalargs[numthreads];
  for (uint32_t i=0; i<numthreads; i++) {
    searchintervalargs[i].a = a;
    searchintervalargs[i].b = b;
    searchintervalargs[i].n = n;
    searchintervalargs[i].primes = primes;
    searchintervalargs[i].numprimes = numprimes;
    searchintervalargs[i].primefactorcountssize = primefactorcountssize;
    pthread_create(&threads[i], NULL, searchinterval, &searchintervalargs[i]);
  }
  for (uint32_t i=0; i<numthreads; i++) pthread_join(threads[i], NULL);
  if (result) {
    printf("%lu is the smallest integer in [%lu, %lu] with at least %u distinct prime factors.\n",result, a,b,n);
  } else {
    printf("There are no integers in [%lu, %lu] with %u prime factors.\n", a,b,n);
  }
  free(primes);
}
/* i7-6700 3.4GHz 4c8t
simon@Bob:~$ time ./smallestwithmprimefactors3.bin 8 10000000000000000000 10000000118047509570 14
10000000118047509570 is the smallest integer in [10000000000000000000, 10000000118047509570] with at least 14 distinct prime factors.

real	22m39.842s
user	181m8.797s
sys	0m1.266s
simon@Bob:~$ time ./smallestwithmprimefactors3.bin 8 9999999939613775911 10000000000000000000  14
There are no integers in [9999999939613775911, 10000000000000000000] with 14 prime factors.

real	11m36.548s
user	92m48.156s
sys	0m0.287s
*/

