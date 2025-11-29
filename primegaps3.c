//===================================================================================================================================
// primegaps3.c: Searches for prime number gaps p_{i+1} - p_i >= gapsize with p_i in range [rangestart, rangeend].
//===================================================================================================================================
// Author: Simon Goater Nov 2025
// Motivated by deleted question https://math.stackexchange.com/questions/5108199/promegaps-kilogaps#comment11002118_5108199
// Question was looking for a list of all prime gaps >= 1000 in the range [5*10^16, 10^17]
// thus extending the results documented in Dr. Thomas R. Nicely's https://oeis.org/A000101/a000101.pdf.
// This program uses a segmented sieve of Eratosthenes. 
// Approx. 2.5 i7-6700 @3.4GHz (all core) CPU years required to cover that interval.
// I recommended a GPU accelerated search instead. Tomás Oliveira e Silva's fast_seive.c performance is similar on my CPU, 
// but is a single threaaded program deploying some obscure algorithm, suggesting a multi-threaded adaptation of it might be 
// significantly faster than both. Primegaps3.c's algorithm is simpler so arguably more easily auditable if anyone wanted to use 
// it for anything serious. As is, this code is given not rigorously tested. 
// 
// COPYRIGHT NOTICE: Copying, modifying, and distributing with conspicuous attribution for any purpose is permitted.
//===================================================================================================================================

// Relevant https://sweet.ua.pt/tos/software/prime_sieve.html

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include "/home/simon/bitarray.c"
#include "/home/simon/Mairsonsprimesieve.c" // https://github.com/FastAsChuff/Primes-List/blob/main/Mairsonsprimesieve.c

// gcc primegaps3.c -o primegaps3.bin -lm -O3 -mssse3 -Wall -std=c11 -pthread


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

//#define MAX64BITPRIMEGAPSIZE 1550U
#define MAX64BITPRIMEGAPSIZE 2048U // Use this if you don't want to assume 1550 is correct.

uint64_t bitarraymod6(uint64_t x) {
  // 6y + 1 -> 2y
  // 6y + 5 -> 2y + 1
  uint64_t y = x/6;
  return 2*y + ((x-6*y) == 5);
}

typedef struct {
  uint8_t *bitarray;
  uint64_t divisionstart;
  uint64_t minbitarrayix;
  uint64_t maxbitarrayix;
  uint32_t *smallprimes;
  uint32_t smallprimessize;
} updatebitarraythreadargs_s;

_Bool updatebitarraythreadcore(uint32_t pix, uint32_t *smallprimes, uint8_t *bitarray, uint64_t divisionstart, uint64_t minbitarrayix, uint64_t maxbitarrayix) {  
    uint32_t p = smallprimes[pix];
    uint64_t p2 = (uint64_t)p*p;
    if (p2 > maxbitarrayix + divisionstart) return false;
    uint64_t ptimes2 = (uint64_t)2*p;
    _Bool pis1mod6 = ((p % 6) == 1);
    uint64_t i = (minbitarrayix + divisionstart >= p2 - p + 1) ? (minbitarrayix + divisionstart - p2 + p - 1)/p : 0;
    //assert(p2 + i*p >= divisionstart + minbitarrayix);
    //if (i) assert(p2 + i*p < p + divisionstart + minbitarrayix);
    // minbitarrayix <= p2 + i*p - divisionstart <= maxbitarrayix
    // i <= (maxbitarrayix + divisionstart - p2)/p
    uint64_t maxi = (maxbitarrayix + divisionstart - p2)/p;
    //assert(p2 + maxi*p >= divisionstart + minbitarrayix);
    //assert(p2 + maxi*p <= divisionstart + maxbitarrayix);    
    // Want to update bitarray only if p2 + i*p - divisionstart = +-1 mod 6.
    // This occurs only when i*p = 0 or 4 mod 6.
    // If pis1mod6, i = 0 or 4 mod 6
    // otherwise, i = 0 or 2 mod 6    
    while (true) {
      if (((i % 6) == 0) || ((i % 6) == (pis1mod6 ? 4 : 2))) break;
      i++;
    }
    uint64_t bitarrayix = p2 + i*p - divisionstart; // i >= (minbitarrayix + divisionstart - p2 + p - 1)/p
    int32_t deltai = 4;
    if ((i % 6) == (pis1mod6 ? 4 : 0)) deltai = 2;
    for (; i<=maxi;) {
      setbitarray(bitarraymod6(bitarrayix), bitarray);
      bitarrayix += (ptimes2 << (deltai == 4));
      i += deltai;
      deltai = 6 - deltai;
    }
    return true;
}

void *updatebitarraythread(void *args) {
  updatebitarraythreadargs_s *updatebitarraythreadargs = (updatebitarraythreadargs_s*)args;
  uint8_t *bitarray = updatebitarraythreadargs->bitarray;
  uint64_t divisionstart = updatebitarraythreadargs->divisionstart;
  uint64_t minbitarrayix;
  uint64_t maxbitarrayix;
  uint32_t *smallprimes = updatebitarraythreadargs->smallprimes;
  uint32_t smallprimessize = updatebitarraythreadargs->smallprimessize;
  // Function Description:
  // =====================
  // Let z(i,p) = p^2 + ip - divisionstart.
  // For all smallprimes p >= 5, set as composite bitarraymod6(z(i,p)) in bitarray
  // such that minbitarrayix <= z(i,p) <= maxbitarrayix only when z(i,p) = +-1 mod 6 and i >= 0.
  // Note: minbitarrayix and maxbitarrayix represent uncompressed indices.
  // ==========================================================================================
  uint32_t primesperblock;
  uint32_t bitarrayblocksize; // update (bitarrayblocksize/24)k block of bitarray at a time.
  uint32_t pix = 2;
  primesperblock = 65536; // Optimise primesperblock and bitarrayblocksize for your hardware & input range.
  bitarrayblocksize = 8*1572864; 
  for (; pix+primesperblock<smallprimessize; pix+=primesperblock) {
    if (smallprimes[pix] > bitarrayblocksize) break;
    minbitarrayix = updatebitarraythreadargs->minbitarrayix;
    maxbitarrayix = MIN(updatebitarraythreadargs->maxbitarrayix, minbitarrayix + bitarrayblocksize - 1);
    while (true) {
      for (uint32_t i=0; i<primesperblock; i++) {
        updatebitarraythreadcore(i+pix, smallprimes, bitarray, divisionstart, minbitarrayix, maxbitarrayix);
      }
      if (maxbitarrayix == updatebitarraythreadargs->maxbitarrayix) break;
      minbitarrayix = 1 + maxbitarrayix;
      maxbitarrayix = MIN(updatebitarraythreadargs->maxbitarrayix, minbitarrayix + bitarrayblocksize - 1);
    }
  }
  minbitarrayix = updatebitarraythreadargs->minbitarrayix;
  maxbitarrayix = updatebitarraythreadargs->maxbitarrayix;
  for (; pix<smallprimessize; pix++) {
    if (!updatebitarraythreadcore(pix, smallprimes, bitarray, divisionstart, minbitarrayix, maxbitarrayix)) break;
  }
  return NULL;
}
/*
void *updatebitarraythreadold(void *args) {
  updatebitarraythreadargs_s *updatebitarraythreadargs = (updatebitarraythreadargs_s*)args;
  uint8_t *bitarray = updatebitarraythreadargs->bitarray;
  uint64_t divisionstart = updatebitarraythreadargs->divisionstart;
  uint64_t minbitarrayix = updatebitarraythreadargs->minbitarrayix;
  uint64_t maxbitarrayix = updatebitarraythreadargs->maxbitarrayix;
  uint32_t *smallprimes = updatebitarraythreadargs->smallprimes;
  uint32_t smallprimessize = updatebitarraythreadargs->smallprimessize;
  // Function Description:
  // =====================
  // Let z(i,p) = p^2 + ip - divisionstart.
  // For all smallprimes p >= 5, set as composite bitarraymod6(z(i,p)) in bitarray
  // such that minbitarrayix <= z(i,p) <= maxbitarrayix only when z(i,p) = +-1 mod 6 and i >= 0.
  // Note: minbitarrayix and maxbitarrayix represent uncompressed indices.
  // ==========================================================================================
  for (uint32_t pix=2; pix<smallprimessize; pix++) {
    uint32_t p = smallprimes[pix];
    uint64_t p2 = (uint64_t)p*p;
    if (p2 > maxbitarrayix + divisionstart) break;
    uint64_t ptimes2 = (uint64_t)2*p;
    _Bool pis1mod6 = ((p % 6) == 1);
    uint64_t i = (minbitarrayix + divisionstart >= p2 - p + 1) ? (minbitarrayix + divisionstart - p2 + p - 1)/p : 0;
    //assert(p2 + i*p >= divisionstart + minbitarrayix);
    //if (i) assert(p2 + i*p < p + divisionstart + minbitarrayix);
    // minbitarrayix <= p2 + i*p - divisionstart <= maxbitarrayix
    // i <= (maxbitarrayix + divisionstart - p2)/p
    uint64_t maxi = (maxbitarrayix + divisionstart - p2)/p;
    //assert(p2 + maxi*p >= divisionstart + minbitarrayix);
    //assert(p2 + maxi*p <= divisionstart + maxbitarrayix);    
    // Want to update bitarray only if p2 + i*p - divisionstart = +-1 mod 6.
    // This occurs only when i*p = 0 or 4 mod 6.
    // If pis1mod6, i = 0 or 4 mod 6
    // otherwise, i = 0 or 2 mod 6    
    while (true) {
      if (((i % 6) == 0) || ((i % 6) == (pis1mod6 ? 4 : 2))) break;
      i++;
    }
    uint64_t bitarrayix = p2 + i*p - divisionstart; // i >= (minbitarrayix + divisionstart - p2 + p - 1)/p
    int32_t deltai = 4;
    if ((i % 6) == (pis1mod6 ? 4 : 0)) deltai = 2;
    for (; i<=maxi;) {
      setbitarray(bitarraymod6(bitarrayix), bitarray);
      bitarrayix += ptimes2*(1 + (deltai == 4));
      i += deltai;
      deltai = 6 - deltai;
    }
  }
  return NULL;
}
*/
typedef struct {
  uint8_t *bitarray;
  uint64_t bitarraynumbytes;
  uint32_t *smallprimes;
  uint32_t smallprimessize;
  uint64_t divisionstart;
  uint64_t divisionend;
  uint32_t numthreads;
} updatebitarrayargs_s;

void *updatebitarray(void *args) {  
  updatebitarrayargs_s *updatebitarrayargs = (updatebitarrayargs_s*)args;
  uint8_t *bitarray = updatebitarrayargs->bitarray;
  uint64_t bitarraynumbytes = updatebitarrayargs->bitarraynumbytes;
  uint32_t *smallprimes = updatebitarrayargs->smallprimes;
  uint32_t smallprimessize = updatebitarrayargs->smallprimessize;
  uint64_t divisionstart = updatebitarrayargs->divisionstart;
  uint64_t divisionend = updatebitarrayargs->divisionend;
  uint32_t numthreads = updatebitarrayargs->numthreads;
  assert(numthreads);
  memset(bitarray, 0, bitarraynumbytes);
  updatebitarraythreadargs_s updatebitarraythreadargs[numthreads];
  uint64_t numsperthread = (divisionend - divisionstart) / numthreads;
  numsperthread -= (numsperthread % 24);
  if ((numthreads == 1) || (numsperthread == 0)) {
    updatebitarraythreadargs[0].bitarray = bitarray;
    updatebitarraythreadargs[0].minbitarrayix = 0;
    updatebitarraythreadargs[0].maxbitarrayix = divisionend - divisionstart + MAX64BITPRIMEGAPSIZE;
    updatebitarraythreadargs[0].divisionstart = divisionstart;
    updatebitarraythreadargs[0].smallprimes = smallprimes;
    updatebitarraythreadargs[0].smallprimessize = smallprimessize;
    updatebitarraythread(updatebitarraythreadargs);
    return NULL;
  }
  pthread_t threads[numthreads];
  for (uint32_t i=0; i<numthreads; i++) {
    updatebitarraythreadargs[i].bitarray = bitarray;
    updatebitarraythreadargs[i].minbitarrayix = i*numsperthread;
    updatebitarraythreadargs[i].maxbitarrayix = (i == (numthreads-1) ? divisionend - divisionstart + MAX64BITPRIMEGAPSIZE : ((i+1)*numsperthread) - 1);
    updatebitarraythreadargs[i].divisionstart = divisionstart;
    updatebitarraythreadargs[i].smallprimes = smallprimes;
    updatebitarraythreadargs[i].smallprimessize = smallprimessize;
    pthread_create(&threads[i], NULL, updatebitarraythread, &updatebitarraythreadargs[i]);
    //updatebitarraythread(&updatebitarraythreadargs[i]);
  }
  for (uint32_t i=0; i<numthreads; i++) pthread_join(threads[i], NULL);
  return NULL;
}

#define popcountu8(x) __builtin_popcount((x) & 0xff)

void printresults(uint8_t *bitarray, uint64_t bitarraynumbytes, uint64_t divisionstart, uint64_t divisionend, uint64_t *gapcount, uint32_t gapsize, uint64_t *primecount, _Bool printprimes) {
  uint64_t localprimecount = 0;
  uint64_t i;
  _Bool neargapsexist = false;
  uint64_t conseq0xff = 0;
  uint64_t minconseq0xff = (gapsize/24) - 1;  
  if (gapsize >= 48) { // Do quick scan to see if there are definately no gaps to count/print.
    for (i=0; i < bitarraynumbytes; i++) {
      uint8_t x = ~bitarray[i];
      if (i < (1+divisionend-divisionstart)/24) localprimecount += popcountu8(x);
      if (x) {
        neargapsexist |= (conseq0xff >= minconseq0xff);
        conseq0xff = 0;
      } else { // All composites
        conseq0xff++;
      }
    }
  }
  neargapsexist |= (conseq0xff >= minconseq0xff);
  if (neargapsexist) {
    uint64_t prevprime = -1;
    uint32_t deltai = 4;
    for (i=1; i <= (divisionend-divisionstart)+MAX64BITPRIMEGAPSIZE;) {
      if (!getbitarray(bitarraymod6(i),bitarray)) {
        if (i <= (divisionend-divisionstart)) (*primecount)++;
        if (prevprime != (uint64_t)-1) {
          if (i+divisionstart - prevprime >= gapsize) {
            if (prevprime <= divisionend) {
              if (printprimes) {
                printf("%lu %lu Diff = %lu\n", prevprime, i+divisionstart, i+divisionstart-prevprime);
                fflush(stdout);
              }
              (*gapcount)++;
            }
          }
        }
        prevprime = i+divisionstart;
      }
      i += deltai;
      deltai = 6 - deltai;
    }
  } else {
    (*primecount) += localprimecount;
  }
  return;
}
/*
void printresultsold(uint8_t *bitarray, uint64_t bitarraynumbytes, uint64_t divisionstart, uint64_t divisionend, uint64_t *gapcount, uint32_t gapsize, uint64_t *primecount, _Bool printprimes) {
  //return;
  uint64_t prevprime = -1;
  uint32_t deltai = 4;
  for (uint64_t i=1; i <= (divisionend-divisionstart)+MAX64BITPRIMEGAPSIZE;) {
    if (!getbitarray(bitarraymod6(i),bitarray)) {
      if (i <= (divisionend-divisionstart)) (*primecount)++;
      if (prevprime != (uint64_t)-1) {
        if (i+divisionstart - prevprime >= gapsize) {
          if (prevprime <= divisionend) {
            if (printprimes) {
              printf("%lu %lu Diff = %lu\n", prevprime, i+divisionstart, i+divisionstart-prevprime);
              fflush(stdout);
            }
            (*gapcount)++;
          }
        }
      }
      prevprime = i+divisionstart;
    }
    i += deltai;
    deltai = 6 - deltai;
  }
}
*/
uint32_t isqrtu64(uint64_t n) {
  if (n < 2) return n;
  uint64_t ai = sqrt(n);
  while (!((ai <= n/ai) && ((ai+1) > n/(ai+1)))) {    
    ai = (ai + n/ai)/2;
  }
  return ai;
}

int main(int argc, char* argv[]) {
  uint64_t rangeendmax = ((uint64_t)-1) - MAX64BITPRIMEGAPSIZE;
  rangeendmax -= (1+(rangeendmax % 24));
  if (argc < 4) {
    printf("This program searches an interval [rangestart, rangeend] for prime gaps >= gapsize, with the lower prime in the interval.\nUsage:- %s rangestart rangeend gapsize [--noprintgaps]\nRAM Requirement = ~3GB\nrangestart >= 24\nrangeend <= %lu\nrangestart, rangeend = 0, 23 mod 24 or search interval will be widened.\nAuthor: Simon Goater Nov 2025\n", argv[0], rangeendmax);
    exit(0);
  }
  uint64_t rangestartinput = atou64(argv[1]);
  uint64_t rangeendinput = atou64(argv[2]);
  uint64_t rangestart = rangestartinput - (rangestartinput % 24);
  uint64_t rangeend = rangeendinput + (23 - (rangeendinput % 24));
  uint32_t gapsize = atou64(argv[3]);
  _Bool printprimes = (argc <= 4 ? true : 0 != strcmp("--noprintgaps", argv[4]));
  assert(rangestart >= 24);
  assert((rangestart % 24) == 0);
  assert((rangeend % 24) == 23);
  assert(gapsize >= 2);
  assert(rangestart < rangeend);
  assert(rangeend < -(uint64_t)MAX64BITPRIMEGAPSIZE);
  uint64_t memmaxbits = 12000000000ULL;
  uint32_t smallprimessize;
  uint32_t numthreads = 8;
  uint32_t upperbound = isqrtu64(rangeend+MAX64BITPRIMEGAPSIZE);
  printf("Building an array of primes <= %u...\n", upperbound);
  uint32_t *smallprimes = Mairsonsprimesieve(upperbound, &smallprimessize);
  assert(smallprimes);
  // Divide the input interval into smaller divisions that can fit in RAM.
  uint64_t starttime = time(0);
  uint64_t numsperdivision = MIN(1 + (rangeend - rangestart), 3*memmaxbits - ((3*memmaxbits) % 24)); // Multiple of 24.
  uint64_t bitarraynumbytes;
  uint8_t *bitarrayupdating = makebitarray((numsperdivision/3)+MAX64BITPRIMEGAPSIZE, &bitarraynumbytes);
  assert(bitarrayupdating);
  uint8_t *bitarray = makebitarray((numsperdivision/3)+MAX64BITPRIMEGAPSIZE, &bitarraynumbytes);
  assert(bitarray);
  uint64_t divisionstart = rangestart;  
  uint64_t divisionstartprev = 0;  
  uint64_t divisionendprev = 0;  
  uint64_t gapcount = 0;
  uint64_t primecount = 0;
  while (divisionstart <= rangeend) {
    // BitArray is [divisionstart,divisionend] + [1+divisionend, divisionend+MAX64BITPRIMEGAPSIZE]
    uint64_t divisionendtemp = divisionstart + numsperdivision - 1;  
    uint64_t divisionend = (divisionendtemp < divisionstart ? rangeend : MIN(rangeend, divisionstart + numsperdivision - 1));  
    printf("Checking interval [%lu, %lu]...\n", divisionstart, divisionend);
    fflush(stdout);
    pthread_t updatebitarrayt;
    updatebitarrayargs_s updatebitarrayargs;
    updatebitarrayargs.bitarray = bitarrayupdating;
    updatebitarrayargs.bitarraynumbytes = bitarraynumbytes;
    updatebitarrayargs.smallprimes = smallprimes;
    updatebitarrayargs.smallprimessize = smallprimessize;
    updatebitarrayargs.divisionstart = divisionstart;
    updatebitarrayargs.divisionend = divisionend;
    updatebitarrayargs.numthreads = numthreads;
    pthread_create(&updatebitarrayt, NULL, updatebitarray, &updatebitarrayargs);
    if (divisionstart > rangestart) {
      printresults(bitarray, bitarraynumbytes, divisionstartprev, divisionendprev, &gapcount, gapsize, &primecount, printprimes);
      memset(bitarray, 0, bitarraynumbytes);
    }
    pthread_join(updatebitarrayt, NULL);
    uint8_t *temp = bitarrayupdating;
    bitarrayupdating = bitarray;
    bitarray = temp;
    if (1+divisionend > rangeend) break;
    divisionstartprev = divisionstart;
    divisionendprev = divisionend;
    divisionstart = 1+divisionend;
  }
  uint64_t endtime = time(0);
  printresults(bitarray, bitarraynumbytes, divisionstart, rangeend, &gapcount, gapsize, &primecount, printprimes);
  printf("Found %lu gaps >= %u, and %lu primes in the interval [%lu, %lu].\n", gapcount, gapsize, primecount, rangestart, rangeend);
  if (endtime - starttime > 5) printf("(%f billion ints/s in segmented sieve)\n", (rangeend - rangestart)/(1000000000.0*(endtime - starttime)));
  free(bitarray);
  free(bitarrayupdating);
  free(smallprimes);
}
// i7-6700 3.4GHz kilogap search performance around 10^17 ~650 million int/s.
// Relevant https://cs.stanford.edu/~knuth/programs/prime-sieve-sparse.w https://primegap-list-project.github.io/lists/prime-gaps-first-occurrences/
// Smallest diff >= 950 >= 5*10^16 is (50067121909334911 50067121909335901 Diff = 990 (found using TOS))...  cont from 50011591999999992
// 990 50067121909335901
/* E.g. ./primegaps3.bin 24000000000000000 24000019200000023 600 
Building an array of primes <= 154919395...
Checking interval [24000000000000000, 24000019200000023]...
24000002051017043 24000002051017747 Diff = 704
24000002699920847 24000002699921467 Diff = 620
24000003851474897 24000003851475521 Diff = 624
24000006105221303 24000006105221917 Diff = 614
24000006511104523 24000006511105127 Diff = 604
24000007004810983 24000007004811599 Diff = 616
24000009479809397 24000009479809997 Diff = 600
24000010262411849 24000010262412493 Diff = 644
24000013645709471 24000013645710079 Diff = 608
24000013913352049 24000013913352661 Diff = 612
24000015732242881 24000015732243493 Diff = 612
24000018376777357 24000018376777991 Diff = 634
Found 12 gaps >= 600, and 509040954 primes in the interval [24000000000000000, 24000019200000023].
*/
