/*
    Copyright (C) 2021-2023 Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <t.flouris@ucl.ac.uk>,
    Department of Genetics, Evolution and Environment,
    University College London, Gower Street, London WC1E 6BT, England
*/

#include "bpp-tools.h"

static unsigned int * z_rndu = NULL;

void rnd_init()
{
   unsigned int seed;
   long i;

   /* z_rndu = (unsigned int)opt_seed; */
   if (sizeof(int) != 4)
      fatal("oh-oh, we are in trouble.  int not 32-bit?  rndu() assumes 32-bit int.");

   if (opt_seed > 0) {
      seed = (unsigned int)opt_seed;
   }
   else {
      unsigned int raw_seed = 0;
      FILE *frand = fopen("/dev/urandom", "r");
      if (frand) {
         if (fread(&raw_seed, sizeof(raw_seed), 1, frand) != 1)
            fatal("failure to read white noise...");
         fclose(frand);
      }
      else {
         raw_seed = (unsigned int)time(NULL);
      }
      seed = raw_seed * 2u + 1u;
      if (!seed)
        seed = 1u;
   }

   assert(opt_threads >= 1);

   if (z_rndu) free(z_rndu);
   z_rndu = (unsigned int *)xmalloc((size_t)opt_threads * sizeof(unsigned int));
   for (i = 0; i < opt_threads; ++i)
     z_rndu[i] = seed;
}

void rnd_fini()
{
  free(z_rndu);
}

double rndu(long index)
{
/* 32-bit integer assumed.
   From Ripley (1987) p. 46 or table 2.4 line 2. 
   This may return 0 or 1, which can be a problem.
*/

   /* the below random number generator is the one used until v4.0.6.
      Change if 0 to if 1 to use it */
   #if 0
   z_rndu[index] = z_rndu[index]*69069 + 1;
   if(z_rndu[index] == 0 || z_rndu[index] == 4294967295)  z_rndu[index] = 13;
   return z_rndu[index]/4294967295.0;
   #else
   z_rndu[index] = z_rndu[index] * 69069 + 1;
   if (z_rndu[index] == 0)  z_rndu[index] = 12345671;
   return ldexp((double)(z_rndu[index]), -32);
   #endif
}


int MultiNomialAliasSetTable(int ncat, double* prob, double* F, int* L)
{
  /* This sets up the tables F and L for the alias algorithm for generating samples from the
     multinomial distribution MN(ncat, p) (Walker 1974; Kronmal & Peterson 1979).

     F[i] has cutoff probabilities, L[i] has aliases.
     I[i] is an indicator: -1 for F[i]<1; +1 for F[i]>=1; 0 if the cell is now empty.

     Should perhaps check whether prob[] sums to 1.
  */
  signed char* I = (signed char*)xmalloc((size_t)ncat * sizeof(signed char));
  int i, j, k, nsmall;

  for (i = 0; i < ncat; i++)  L[i] = -9;
  for (i = 0; i < ncat; i++)  F[i] = ncat * prob[i];
  for (i = 0, nsmall = 0; i < ncat; i++) {
    if (F[i] >= 1)  I[i] = 1;
    else { I[i] = -1; nsmall++; }
  }
  for (i = 0; nsmall > 0; i++) {
    for (j = 0; j < ncat; j++)  if (I[j] == -1) break;
    for (k = 0; k < ncat; k++)  if (I[k] == 1)  break;
    if (k == ncat)  break;

    L[j] = k;
    F[k] -= 1 - F[j];
    if (F[k] < 1) { I[k] = -1; nsmall++; }
    I[j] = 0;  nsmall--;
  }

  free(I);
  return(0);
}

int MultiNomialAlias(long index, int n, int ncat, double* F, int* L, int* nobs)
{
  /* This generates multinomial samples using the F and L tables set up before,
     using the alias algorithm (Walker 1974; Kronmal & Peterson 1979).

     F[i] has cutoff probabilities, L[i] has aliases.
     I[i] is an indicator: -1 for F[i]<1; +1 for F[i]>=1; 0 if the cell is now empty.
  */
  int i, k;
  double r;

  for (i = 0; i < ncat; i++)  nobs[i] = 0;
  for (i = 0; i < n; i++) {
    r = rndu(index) * ncat;
    k = (int)r;
    r -= k;
    if (r <= F[k]) nobs[k]++;
    else           nobs[L[k]]++;
  }
  return (0);
}
