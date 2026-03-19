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
   int seed = (int)opt_seed;
   long i;

   /* z_rndu = (unsigned int)opt_seed; */
   if (sizeof(int) != 4)
      fatal("oh-oh, we are in trouble.  int not 32-bit?  rndu() assumes 32-bit int.");

   if (seed <= 0) {
      FILE *frand = fopen("/dev/urandom", "r");
      if (frand) {
         if (fread(&seed, sizeof(int), 1, frand) != 1)
            fatal("failure to read white noise...");
         fclose(frand);
         seed = abs(seed * 2 - 1);
      }
      else {
         seed = abs(1234 * (int)time(NULL) + 1);
      }
   }

   assert(opt_threads >= 1);

   if (z_rndu) free(z_rndu);
   z_rndu = (unsigned int *)xmalloc((size_t)opt_threads * sizeof(unsigned int));
   for (i = 0; i < opt_threads; ++i)
     z_rndu[i] = (unsigned int)seed;
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
