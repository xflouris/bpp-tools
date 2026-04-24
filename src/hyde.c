/*
    Copyright (C) 2022-2023 Tomas Flouri

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

#define NT_CHARS 4
#define SAVINGS

static double * abba_tbl = NULL;
static double * baba_tbl = NULL;
static long * patt_count = NULL;

static const long thread_index_zero = 0;

/* biallelic only */
static const unsigned int map_solid_and_biallelic[256] =
 {
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  1,  0,  1,  0,  0,  0,  1,  0,  0,  0,  1,  0,  1,  0,  0,
   0,  0,  1,  1,  1,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,
   0,  1,  0,  1,  0,  0,  0,  1,  0,  0,  0,  1,  0,  1,  0,  0,
   0,  0,  1,  1,  1,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
 };

static const double allele_weights[256][4] =
 {
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  1,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   1,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   1,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0, 0.5, 0.5},
   {  0,   0,   0,   0},
   {0.5, 0.5,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {0.5,   0, 0.5,   0},
   {  0, 0.5, 0.5,   0},
   {  0,   0,   0,   1},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {0.5,   0,   0, 0.5},
   {  0,   0,   0,   0},
   {  0, 0.5,   0, 0.5},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  1,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   1,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   1,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0, 0.5, 0.5},
   {  0,   0,   0,   0},
   {0.5, 0.5,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {0.5,   0, 0.5,   0},
   {  0, 0.5, 0.5,   0},
   {  0,   0,   0,   1},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {0.5,   0,   0, 0.5},
   {  0,   0,   0,   0},
   {  0, 0.5,   0, 0.5},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0},
   {  0,   0,   0,   0}
 };

/* TF: 7/12/2023 */
#if 0
static long perms[24][4] =
 {
   {0,1,2,3},
   {0,1,3,2},
   {0,2,1,3},
   {0,2,3,1},
   {0,3,1,2},
   {0,3,2,1},
   {1,0,2,3},
   {1,0,3,2},
   {1,2,0,3},
   {1,2,3,0},
   {1,3,0,2},
   {1,3,2,0},
   {2,0,1,3},
   {2,0,3,1},
   {2,1,0,3},
   {2,1,3,0},
   {2,3,0,1},
   {2,3,1,0},
   {3,0,1,2},
   {3,0,2,1},
   {3,1,0,2},
   {3,1,2,0},
   {3,2,0,1},
   {3,2,1,0}
 };
#else
static long perms[6][4] =
 {
   {0,1,2,3},
   {0,2,1,3},
   {1,0,2,3},
   {1,2,0,3},
   {2,0,1,3},
   {2,1,0,3}
 };
#endif

#if 1
static long getnumdigits(long n)
{
  if (n == 0) return 1;
  return (long)(floor(log10(labs(n))+1));
}

static const char * extract_seqtag(const char * label, long msa_index)
{
  char * tag = strchr(label, '^');

  if (!tag)
    fatal("Cannot find species tag on sequence %s of locus %ld",
          label, msa_index);
  
  /* skip the '^' mark */
  ++label;
  if (!(*label))
    fatal("Sequence %s of locus %ld contains no label", label, msa_index);

  return tag;
}
#endif

static int cb_cmp_double(const void * a, const void * b)
{
  double * x = (double *)a;
  double * y = (double *)b;

  if (*x > *y) return 1;
  if (*x < *y) return -1;

  return 0;
}

#if 0
static long msa_map_compat(msa_t * msa, list_t * maps)
{
  char * label;
  long i;

  for (i = 0; i < msa->count; ++i)
  {
    /* get sequence */

    list_item_t * li = list->head;
    while (li)
    {
      mapping_t * map = (mapping_t *)(li->data);
      if (map->species)
      li = li->next;
    }

    
  }
    
}
#endif

/* given the indices of four species (species), we condense the alignment (msa)
   by collecting all sequences from those species using the index structure
   and the maparray

   Params:
     msa      : the concatenated alignment 
     maparray : the list of mappings (individual->species)
     index    : points to species element for given sequence in msa (0 to 3)
     species  : the 4 species to be used for d, points to maparray indices
     vec      : stores prob of each base for each site (size: 4x4*sites doubles)

*/
/* TODO: Get species tags */

static char invcharmap_acgt[16] =
{
  '-', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
  'T', 'W', 'Y', 'H', 'K', 'D', 'B', '-'
};

static long ends_with_hat_str(const char * x, const char * suffix)
{
  size_t xlen, slen;
  const char * p;
  const char * hat;

  xlen = strlen(x);
  slen = strlen(suffix);

  if (slen >= xlen)
    return 0;

  p = x + xlen - slen;
  hat = p-1;
  if (!strcmp(p,suffix) && *hat == '^')
    return 1;

  return 0;
}

static char ** map_sequences_to_species(msa_t * msa,
                                        list_t * maplist,
                                        long ** indexptr,
                                        long ** seq_per_species,
                                        long * sp_count)
{
  long i,j;
  long * index;
  char * lbl;
  list_item_t * li;
  char ** labels;

  if (!maplist || !msa) return NULL;

  labels = (char **)xmalloc((size_t)msa->count * sizeof(char *));

  for (i = 0; i < msa->count; ++i)
  {
    lbl= msa->label[i];

    li = maplist->head;
    while (li)
    {
      mapping_t * map = (mapping_t *)(li->data);

      if (ends_with_hat_str(lbl, map->individual))
      {
        labels[i] = xstrdup(map->species);
        break;
      }
      li = li->next;
    }

    if (!li)
    {
      free(labels);
      return NULL;
    }
  }

  printf("Came here\n");
  index = (long *)xmalloc((size_t)msa->count * sizeof(long));
  for (i = 0; i < msa->count; ++i)
    index[i] = -1;

  long unique = 0;
  /* count unique species labels */
  for (i = 0; i < msa->count; ++i)
  {
    if (index[i] != -1) continue;
    unique++;

    for (j = i+1; j < msa->count; ++j)
    {
      if (index[j] != -1) continue;

      if (!strcmp(labels[i], labels[j]))
        index[j] = i;
    }
  }

  char ** sp_labels = (char **)xmalloc((size_t)unique * sizeof(char *));
  long * sp_seqs = (long *)xcalloc((size_t)unique, sizeof(long));
  for (j=0, i=0; i < msa->count; ++i)
  {
    if (index[i] == -1)
    {
      sp_labels[j] = xstrdup(labels[i]);
      sp_seqs[j]++;
      index[i] = j++;
    }
    else
    {
      assert(index[index[i]] < j &&
             !strcmp(labels[i],labels[index[i]]) &&
             !strcmp(labels[i],sp_labels[index[index[i]]]));

      index[i] = index[index[i]];
      sp_seqs[index[i]]++;
    }
  }
  assert(j == unique);

  *indexptr = index;
  *sp_count = unique;
  *seq_per_species = sp_seqs;
  if (opt_debug)
    xdebug("Unique species labels: %ld\n", unique);

  for (i = 0; i < msa->count; ++i)
    free(labels[i]);
  free(labels);

  return sp_labels;
}

#if 1
static msa_t * condense_msa(msa_t * msa,
                            long * sp_seqassign,
                            char ** sp_labels,
                            long species_count,
//                            long * sp_seqcount,
                            double ** vecptr)
{
  /* TODO: Assert that we have 4 species */
  long i,j,k,m;
  long sp_index;
  msa_t * newmsa;
  double * sitevec;
  double * vec;
  double * pvec;
  unsigned int * sp_allele_code;

  long * sp_seqcount;

  sp_seqcount = (long *)xcalloc((size_t)species_count*msa->length,sizeof(long));

  sp_allele_code = (unsigned int *)xmalloc((size_t)species_count *
                                           sizeof(unsigned int));

  unsigned int keep;

  sitevec = (double *)xmalloc((size_t)species_count*4*sizeof(double));

  /* calculate amount of space needed for condensed alignment */
  for (i = 0; i < msa->length; ++i)
  {
    /* reset character */
    for (j = 0; j < species_count; ++j)
      sp_allele_code[j] = 0;

    for (j = 0; j < msa->count; ++j)
    {
      printf("HERE!\n");
      printf("j=%ld\n", j);
      sp_index = sp_seqassign[j]; 

      /* upto biallelic characters */
      if (map_solid_and_biallelic[(int)msa->sequence[j][i]])
      {
        /* update allele character for species sp_index */
        printf("%ld %ld\n", j, i);
        sp_allele_code[sp_index] |= pll_map_nt[(int)msa->sequence[j][i]];
      }
    }
  }

  /* allocate msa */
  newmsa = (msa_t *)xcalloc(1,sizeof(msa_t));
  newmsa->count = species_count;
  newmsa->length = msa->length;

  newmsa->sequence = (char **)xmalloc((size_t)species_count*sizeof(char *));
  newmsa->label = (char **)xmalloc((size_t)species_count*sizeof(char *));
  for (i = 0; i < species_count; ++i)
  {
    newmsa->sequence[i] = (char *)xmalloc((size_t)(msa->length+1)*sizeof(char));
    newmsa->label[i] = xstrdup(sp_labels[i]);
    newmsa->sequence[i][newmsa->length] = 0;
  }

  /* allocate probability vector (size: sites * states * species) */
  vec = (double *)xmalloc((size_t)(species_count*4*msa->length)*sizeof(double));
  pvec = vec;

  /* go through sites */
  for (m=0, i = 0; i < msa->length; ++i)
  {
    /* reset character */
    for (j = 0; j < species_count; ++j)
      sp_allele_code[j] = 0;

    for (j = 0; j < species_count*4; ++j)
      sitevec[j] = 0;

    for (j = 0; j < msa->count; ++j)
    {
      sp_index = sp_seqassign[j];

      /* upto biallelic characters */
      if (map_solid_and_biallelic[(int)msa->sequence[j][i]])
      {
        char ntchar = msa->sequence[j][i];
        for (k = 0; k < 4; ++k)
          sitevec[sp_index*4+k] += allele_weights[(int)ntchar][k];
        #if 0
        printf("[DBG] %s   %c -> %f %f %f %f     offset: %ld\n",
               sp_labels[sp_index],
               ntchar,
               allele_weights[(int)ntchar][0],
               allele_weights[(int)ntchar][1],
               allele_weights[(int)ntchar][2],
               allele_weights[(int)ntchar][3],
               species_count*sp_index);
        #endif

        /* update allele character for species sp_index */
        sp_allele_code[sp_index] |= pll_map_nt[(int)msa->sequence[j][i]];

        sp_seqcount[i*species_count+sp_index]++;
      }
    }

    #if 0
    for (keep=1, j = 0; j < species_count; ++j)
      keep = sp_allele_code[j] && keep;
    #endif

    #if 0
    printf("Writing data...\n");
    #endif
    double * sitepvec = pvec;
    for (j = 0; j < species_count; ++j)
    {
      newmsa->sequence[j][m] = invcharmap_acgt[sp_allele_code[j]];
      for (k = 0; k < 4; ++k)
      {
        pvec[k] = sitevec[4*j+k];
        #if 0
        printf(" %f", pvec[k]);
        #endif
      }
      pvec += 4;
    }
    #if 0
    printf("\n");
    printf("Repeat...\n");
    double * xxx = vec + species_count*4*i;
    for (k = 0; k < species_count*4; ++k)
      printf(" %f", xxx[k]);
    printf("\n");
    #endif

    ++m;
  }


  free(sitevec);

  *vecptr = vec;

  #if 0
  printf("[DBG] printing vec\n");
  for (i = 0; i < species_count; ++i)
  {
    printf("%s\n", sp_labels[i]);
    for (j = 0; j < newmsa->length; ++j)
    {
      double * v = vec + 4*(species_count*j+i);
      printf(" %ld:([%f,%f,%f,%f)", j+1, v[0], v[1], v[2], v[3]);
    }
    printf("\n");
  }
  for (i = 0; i < species_count; ++i)
    printf(" %s -> %ld\n", sp_labels[i], sp_seqcount[i]);
  #endif

  /* normalize */
  /* TODO: this is wrong, we ened sp_seqcount to be over sites*species, e.g. store # sequences at each site */
  for (i = 0; i < msa->length; ++i)
  {
    for (j = 0; j < species_count; ++j)
    {
      if (!sp_seqcount[i*species_count+j]) continue;

      for (k = 0; k < 4; ++k)
        vec[i*species_count*4 + 4*j+k] /= sp_seqcount[i*species_count+j];
    }
  }

  free(sp_seqcount);
  free(sp_allele_code);

  /* propagate pattern-compression metadata: condense_msa collapses rows
     (individuals -> species), not sites, so pattern_weights carry over */
  if (msa->pattern_weights)
  {
    newmsa->pattern_weights =
      (unsigned int *)xmalloc((size_t)msa->length * sizeof(unsigned int));
    memcpy(newmsa->pattern_weights, msa->pattern_weights,
           (size_t)msa->length * sizeof(unsigned int));
    newmsa->compress_model = msa->compress_model;
  }

  return newmsa;
}
#endif

static void show_tags(msa_t * msa)
{
  long i;

  for (i = 0; i < msa->count; ++i)
  {
    const char * x = extract_seqtag(msa->label[i], 1);
    printf("Tag: %s\n", x);
  }
}

static msa_t ** resample_msa_list(msa_t ** msa_list, long count)
{
  long i;
  msa_t ** new_msa_list;

  new_msa_list = (msa_t **)xmalloc((size_t)count * sizeof(msa_t *));

  for (i = 0; i < count; ++i)
  {
    new_msa_list[i] = msa_list[(long)(rndu(0)*count)];
  }

  return new_msa_list;
}

static msa_t * resample_sites(msa_t * msa)
{
  /* resample sites in msa with replacement */
  long i,j;
  msa_t * repl;

  long * indices = (long *)xmalloc((size_t)msa->length * sizeof(long));

  /* resample indices */
  for (i = 0; i < msa->length; ++i)
  {
    indices[i] = (long)(rndu(0)*msa->length);
  }

  /* create alignment */
  repl = (msa_t *)xcalloc(sizeof(msa_t),1);
  repl->count = msa->count;
  repl->length = msa->length;
  repl->label = (char **)xmalloc((size_t)msa->count * sizeof(char *));
  repl->sequence = (char **)xmalloc((size_t)msa->count * sizeof(char *));
  for (i = 0; i < msa->count; ++i)
  {
    repl->sequence[i] = (char *)xmalloc((size_t)msa->length * sizeof(char));
    repl->label[i] = xstrdup(msa->label[i]);
  }

  for (i = 0; i < msa->length; ++i)
  {
    for (j = 0; j < msa->count; ++j)
      repl->sequence[j][i] = msa->sequence[j][indices[i]];
  }
  free(indices);

  return repl;
}

static void bootstrap(msa_t * msa)
{
  long i;

  for (i = 0; i < 100; ++i)
  {
    
  }
  
}

static void abba_baba_score(unsigned int * s,
                            double * abbaptr,
                            double * babaptr,
                            long * patsptr)
{
  int ii[4] = {0,0,0,0};
  unsigned int code;
  long pats = 0;
  long abba = 0;
  long baba = 0;

  for (ii[0] = 0; ii[0] < 4; ++ii[0])
  {
    for (ii[1] = 0; ii[1] < 4; ++ii[1])
    {
      for (ii[2] = 0; ii[2] < 4; ++ii[2])
      {
        for (ii[3] = 0; ii[3] < 4; ++ii[3])
        {
          code = (((s[0] >> ii[0]) & 1) << 3) |
                 (((s[1] >> ii[1]) & 1) << 2) |
                 (((s[2] >> ii[2]) & 1) << 1) |
                 ((s[3] >> ii[3]) & 1);
          

          if (code == 0xf)
          {
            pats++;

            if (ii[0] == ii[3] && ii[1] == ii[2] && ii[0] != ii[1]) abba++;
            if (ii[0] == ii[2] && ii[1] == ii[3] && ii[0] != ii[1]) baba++;
          }
        }
      }
    }
  }
  assert(pats);
  *abbaptr = abba/(double)pats;
  *babaptr = baba/(double)pats;
  *patsptr = pats;
}

static void print_score_table()
{
  char nt[15] = "ACGTRYSWKMBDHVN";
  int ii[4] = {0,0,0,0};
  char site[5];
  double abba;
  double baba;
  long index;


  printf("Pattern   Code   f(abba)   f(baba)  pats\n");
  printf("----------------------------------------\n");
  
  site[4] = 0;
  for (ii[0] = 0; ii[0] < 0x1; ++ii[0])
  {
    for (ii[1] = 1; ii[1] < 0x0f; ++ii[1])
    {
      for (ii[2] = 0; ii[2] < 0x1; ++ii[2])
      {
        for (ii[3] = 0; ii[3] < 0x2; ++ii[3])
        {

          site[0] = nt[ii[0]];
          site[1] = nt[ii[1]];
          site[2] = nt[ii[2]];
          site[3] = nt[ii[3]];
          unsigned int s[4] = {pll_map_nt[(int)site[0]],
                               pll_map_nt[(int)site[1]],
                               pll_map_nt[(int)site[2]],
                               pll_map_nt[(int)site[3]]};
          //abba_baba_score(s,&abba,&baba,&pats);
          
          index = s[0] | (s[1] << 4) | (s[2] << 8) | (s[3] << 12);
          abba = abba_tbl[index];
          baba = baba_tbl[index];
          printf("  %c%c%c%c ", site[0], site[1], site[2], site[3]);

          printf("  %5ld", index);
          printf("  %f", abba);
          printf("  %f", baba);
          printf("  %4ld\n", patt_count[index]);
        }
      }
    }
  }
}

static void debug_print_vec(double * vec, msa_t * msa)
{
  long i,j,k;
  double * p = vec;
  long sites;
  long species;

  sites = msa->length;
  species = msa->count;

  #if 0
  for (i = 0; i < sites; ++i)
  {
    printf(" %ld:", i+1);
    for (j = 0; j < species; ++j)
    {
      for (k = 0; k < NT_CHARS; ++k)
        printf(" %f", *p++);
      printf(" ");
    }
    printf("\n");
  }
  printf("\n");
  #endif

  xdebug("Per site nucleotide character probability vectors for each species:");
  for (i = 0; i < species; ++i)
  {
    xdebug("%s  ", msa->label[i]);
    for (j = 0; j < msa->length; ++j)
    {
      double * v = vec + 4*(species*j+i);
      xdebug(" %ld: %f,%f,%f,%f", j+1, v[0], v[1], v[2], v[3]);
    }
  }
}

static msa_t * subsample_msa(msa_t * msa,
                             char ** labels,
                             double * vec,
                             double * abba_vec,
                             double * baba_vec,
                             double * bbaa_vec,
                             double **outvec,
                             long labels_count)
{
  long i,j,k;
  double * v;
  msa_t * newmsa;
  long span;

  v = (double *)xmalloc((size_t)(labels_count*NT_CHARS*msa->length)*sizeof(double));
  span = NT_CHARS*labels_count;

  /* allocate new alignment data structure */
  newmsa = (msa_t *)xcalloc(1,sizeof(msa_t));
  newmsa->count = labels_count;
  newmsa->length = msa->length;
  newmsa->sequence = (char **)xmalloc((size_t)labels_count*sizeof(char *));
  newmsa->label = (char **)xmalloc((size_t)labels_count*sizeof(char *));

  for (i = 0; i < labels_count; ++i)
  {
    for (j = 0; j < msa->count; ++j)
    {
      if (!strcmp(labels[i], msa->label[j]))
      {
        newmsa->sequence[i] = (char *)xmalloc((size_t)(msa->length+1)*sizeof(char));
        newmsa->label[i] = xstrdup(labels[i]);
        memcpy(newmsa->sequence[i], msa->sequence[j], msa->length*sizeof(char));
        newmsa->sequence[i][newmsa->length] = 0;

        for (k = 0; k < msa->length; ++k)
          memcpy(v+k*NT_CHARS*labels_count+NT_CHARS*i,
                 vec+k*NT_CHARS*msa->count+NT_CHARS*j,
                 NT_CHARS*sizeof(double));

        break;
      }
    }
    assert(j < msa->count);
  }
  *outvec = v;

  /* new code to precalculate per-site per-species abba-baba scores */
  #if defined(SAVINGS)
  double * sitevec = v;
  assert(newmsa->count == 4);
  memset(abba_vec,0,msa->length*sizeof(double));
  memset(baba_vec,0,msa->length*sizeof(double));
  memset(bbaa_vec,0,msa->length*sizeof(double));
  for (i = 0; i < msa->length; ++i)
  {
    double * sp1 = sitevec+0x00;
    double * sp2 = sitevec+0x04;
    double * sp3 = sitevec+0x08;
    double * sp4 = sitevec+0x0c;

    for (j = 0; j < 4; ++j)
    {
      for (k = 0; k < 4; ++k)
      {
        if (j == k)
          continue;
        
        abba_vec[i] += sp1[j]*sp2[k]*sp3[k]*sp4[j];
        baba_vec[i] += sp1[k]*sp2[j]*sp3[k]*sp4[j];
        bbaa_vec[i] += sp1[k]*sp2[k]*sp3[j]*sp4[j];
      }
    }
    sitevec += 16;
  }
  if (opt_debug)
  {
    /* abba */
    printf("abba_vec\n--------\n");
    printf("  [%.6f", abba_vec[0]);
    for (i = 1; i < msa->length; ++i)
    {
      printf(",%.6f", abba_vec[i]);
    }
    printf("]\n");

    /* baba */
    printf("baba_vec\n--------\n");
    printf("  [%.6f", baba_vec[0]);
    for (i = 1; i < msa->length; ++i)
    {
      printf(",%.6f", baba_vec[i]);
    }
    printf("]\n");

    /* bbaa */
    printf("bbaa_vec\n--------\n");
    printf("  [%.6f", bbaa_vec[0]);
    for (i = 1; i < msa->length; ++i)
    {
      printf(",%.6f", bbaa_vec[i]);
    }
    printf("]\n");
  }


  #endif

  /* propagate pattern-compression metadata (site count is unchanged) */
  if (msa->pattern_weights)
  {
    newmsa->pattern_weights =
      (unsigned int *)xmalloc((size_t)msa->length * sizeof(unsigned int));
    memcpy(newmsa->pattern_weights, msa->pattern_weights,
           (size_t)msa->length * sizeof(unsigned int));
    newmsa->compress_model = msa->compress_model;
  }

  return newmsa;
}

static void calculate_d(msa_t * msa)
{
  long i;
  long * pats = NULL;
  double abba = 0;
  double baba = 0;
  
  printf("--------\n");

  /* change order according to CSV options */
  char * s1 = msa->sequence[0];
  char * s2 = msa->sequence[1];
  char * s3 = msa->sequence[2];
  char * s4 = msa->sequence[3];

  pats = (long *)xcalloc(msa->length,sizeof(long));
  for (i = 0; i < msa->length; ++i)
  {
    long i1 = pll_map_nt[(int)s1[i]];
    long i2 = pll_map_nt[(int)s2[i]];
    long i3 = pll_map_nt[(int)s3[i]];
    long i4 = pll_map_nt[(int)s4[i]];

    pats[i] = i1 | (i2 << 4) | (i3 << 8) | (i4 << 12);
  }

  #if 1
  printf("Number of sites: %d\n", msa->length);
  for (i = 0; i < msa->length; ++i)
  {
    long s = pats[i];
    #if 1
    printf(" pattern: %c%c%c%c code: %5ld abba: %f baba: %f  pats: %ld\n",
           msa->sequence[0][i],
           msa->sequence[1][i],
           msa->sequence[2][i],
           msa->sequence[3][i],
           s, abba_tbl[s], baba_tbl[s], patt_count[s]);
    #endif
    abba += abba_tbl[s];
    baba += baba_tbl[s];
  }
  #endif

  printf("abba: %f     (%.2f%%)\n", abba, (abba /  msa->length)*100);
  printf("baba: %f     (%.2f%%)\n", baba, (baba /  msa->length)*100);

  free(pats);
  /* 
  decode_site(s);

  printf("Checking from precomputed table:\n");
  index = s[0] | (s[1] << 4) | (s[2] << 8) | (s[3] << 12);
  printf("abba: %f\n", abba_tbl[index]);
  printf("baba: %f\n", baba_tbl[index]);
  */
}

static double calculate_d_vec(msa_t * msa, double * vec, double * fabbaptr, double * fbabaptr, double * fbbaaptr)
{
  long i;
  long abba_count = 0;
  long baba_count = 0;
  long bbaa_count = 0;
  double fabba = 0;
  double fbaba = 0;
  double fbbaa = 0;
  double fabbabababbaa = 0;
  double fabbababa = 0;
  double patscore = 0;
  double dscore = 0;
  double gamma = 0;

  /* change order according to CSV options */
  char * s1 = msa->sequence[0];
  char * s2 = msa->sequence[1];
  char * s3 = msa->sequence[2];
  char * s4 = msa->sequence[3];

  for (i = 0; i < msa->length; ++i)
  {
    #if 0
    printf("Site %ld : %c%c%c%c\n", i+1, s1[i],s2[i],s3[i],s4[i]);
    #endif
    long i1 = pll_map_nt[(int)s1[i]];
    long i2 = pll_map_nt[(int)s2[i]];
    long i3 = pll_map_nt[(int)s3[i]];
    long i4 = pll_map_nt[(int)s4[i]];

    unsigned int s[4] = {i1,i2,i3,i4};
    int ii[4] = {0,0,0,0};

    for (ii[0] = 0; ii[0] < 4; ++ii[0])
    {
      for (ii[1] = 0; ii[1] < 4; ++ii[1])
      {
        for (ii[2] = 0; ii[2] < 4; ++ii[2])
        {
          for (ii[3] = 0; ii[3] < 4; ++ii[3])
          {
            unsigned int code,j;
            code = (((s[0] >> ii[0]) & 1) << 3) |
                   (((s[1] >> ii[1]) & 1) << 2) |
                   (((s[2] >> ii[2]) & 1) << 1) |
                   ((s[3] >> ii[3]) & 1);

            if (code == 0xf)
            {
              long found = 0;

              if (ii[0] == ii[3] && ii[1] == ii[2] && ii[0] != ii[1])
              {
                #if 0
                printf(ANSI_COLOR_BLUE);
                #endif
                abba_count++;
                found = 1;
              }
              if (ii[0] == ii[2] && ii[1] == ii[3] && ii[0] != ii[1])
              {
                #if 0
                printf(ANSI_COLOR_RED);
                #endif
                baba_count++;
                found = 1;
              }
              if (ii[0] == ii[1] && ii[2] == ii[3] && ii[0] != ii[2])
              {
                #if 0
                printf(ANSI_COLOR_RED);
                #endif
                bbaa_count++;
                found = 1;
              }

              #if 0
              for (j = 0; j < 4; ++j)
              {
                switch(ii[j])
                {
                  case 0:
                    printf("A");
                    break;
                  case 1:
                    printf("C");
                    break;
                  case 2:
                    printf("G");
                    break;
                  case 3:
                    printf("T");
                    break;
                  default:
                    fatal("Internal error");
                }
              }
              printf(" ");
              long j;
              for (j = 0; j < 4; ++j)  /* site position */
              {
                printf("%d",ii[j]);
              }
              #endif
              double * v = vec + i*msa->count*4;
              double patscore = 0;
              #if 0
              printf("  ");
              for (j = 0; j < 4; ++j)
              {
                printf("  %f", v[j*4+ii[j]]);
              }
              #endif

              if (found)
              {
                double patscore = 1;
                for (j = 0; j < 4; ++j)
                  patscore *= v[j*4+ii[j]];
                #if 0
                printf("  %f", patscore);
                #endif

                /* weight by pattern multiplicity if input is compressed */
                if (msa->pattern_weights)
                  patscore *= (double)msa->pattern_weights[i];

                fabbabababbaa += patscore;

                if (ii[0] == ii[3] && ii[1] == ii[2] && ii[0] != ii[1])
                {
                  /* ABBA */
                  fabba += patscore;
                  fabbababa += patscore;
                }
                else if (ii[0] == ii[2] && ii[1] == ii[3] && ii[0] != ii[1])
                {
                  /* BABA */
                  fbaba += patscore;
                  fabbababa += patscore;
                }
                else
                {
                  /* BBAA */
                  fbbaa += patscore;
                }
              }

              #if 0
              if (ii[0] == ii[3] && ii[1] == ii[2] && ii[0] != ii[1])
              {
                printf(ANSI_COLOR_RESET);
              }
              if (ii[0] == ii[2] && ii[1] == ii[3] && ii[0] != ii[1])
              {
                printf(ANSI_COLOR_RESET);
              }
              printf("\n");
              #endif
            }
          }
        }
      }
    }

  }

  if (abba_count+baba_count != 0)
    dscore = (fabba-fbaba)/fabbababa;

  if (fbbaa-2*fbaba+fabba != 0)
    gamma = (fbbaa-fbaba)/(fbbaa-2*fbaba+fabba);
  else
   gamma = -1;  /* undefined */

  if (fabbaptr)
    *fabbaptr = fabba;
  if (fbabaptr)
    *fbabaptr = fbaba;

  if (opt_debug)
  {
    //printf("\n");
    xdebug(ANSI_COLOR_BLUE
           "c(abba): %ld (%.2f%%)   f(abba): %f (%.2f%%)" ANSI_COLOR_RESET,
           abba_count, ((double)abba_count / msa->length)*100,
           fabba, (fabba/(fabba+fbaba))*100);
    xdebug(ANSI_COLOR_RED
           "c(baba): %ld (%.2f%%)   f(baba): %f (%.2f%%)" ANSI_COLOR_RESET,
           baba_count, ((double)baba_count / msa->length)*100,
           fbaba, (fbaba/(fabba+fbaba))*100);
    xdebug("f(abba)-f(baba): %f", fabba-fbaba);
    xdebug("f(abba)+f(baba): %f", fabbababa);
    if (abba_count+baba_count == 0)
      xdebug("d-score: 0 -- No abba/baba patterns found!");
    else
      xdebug("d-score: %f", (fabba-fbaba)/fabbababa);

    if (fbbaa-2*fbaba+fabba != 0)
      xdebug("gamma: %f", gamma);
  }
  return dscore;
}

static void precompute_table()
{
  char nt[15] = "ACGTRYSWKMBDHVN";
  int ii[4] = {0,0,0,0};
  char site[5];
  double abba;
  double baba;
  long index;
  long pats;

  abba_tbl = (double *)xcalloc(65536,sizeof(double));
  baba_tbl = (double *)xcalloc(65536,sizeof(double));
  patt_count = (long *)xcalloc(65536,sizeof(long));

  site[4] = 0;
  for (ii[0] = 0; ii[0] < 0xf; ++ii[0])
  {
    for (ii[1] = 0; ii[1] < 0xf; ++ii[1])
    {
      for (ii[2] = 0; ii[2] < 0xf; ++ii[2])
      {
        for (ii[3] = 0; ii[3] < 0xf; ++ii[3])
        {
          site[0] = nt[ii[0]];
          site[1] = nt[ii[1]];
          site[2] = nt[ii[2]];
          site[3] = nt[ii[3]];
          unsigned int s[4] = {pll_map_nt[(int)site[0]],
                               pll_map_nt[(int)site[1]],
                               pll_map_nt[(int)site[2]],
                               pll_map_nt[(int)site[3]]};
          abba_baba_score(s,&abba,&baba,&pats);

          index = s[0] | (s[1] << 4) | (s[2] << 8) | (s[3] << 12);
          abba_tbl[index]   = abba;
          baba_tbl[index]   = baba;
          patt_count[index] = pats;
        }
      }
    }
  }

  #if 0
  long i;
  double abba_sum = 0;
  double baba_sum = 0;
  for (i = 0; i < 65536; ++i)
  {
    abba_sum += abba_tbl[i];
    baba_sum += baba_tbl[i];
  }
  printf("abba sum: %f\n", abba_sum);
  printf("baba sum: %f\n", baba_sum);
  #endif
}

static int debug_decode_site(unsigned int * s)
{
  int ii[4] = {0,0,0,0};
  unsigned int code,j;
  long total = 0;
  long abba = 0;
  long baba = 0;

  for (ii[0] = 0; ii[0] < 4; ++ii[0])
  {
    for (ii[1] = 0; ii[1] < 4; ++ii[1])
    {
      for (ii[2] = 0; ii[2] < 4; ++ii[2])
      {
        for (ii[3] = 0; ii[3] < 4; ++ii[3])
        {
          code = (((s[0] >> ii[0]) & 1) << 3) |
                 (((s[1] >> ii[1]) & 1) << 2) |
                 (((s[2] >> ii[2]) & 1) << 1) |
                 ((s[3] >> ii[3]) & 1);

          if (code == 0xf)
          {
            total++;
            for (j = 0; j < 4; ++j)
            {
              switch(ii[j])
              {
                case 0:
                  printf("A");
                  break;
                case 1:
                  printf("C");
                  break;
                case 2:
                  printf("G");
                  break;
                case 3:
                  printf("T");
                  break;
                default:
                  fatal("Internal error");
              }
            }
            printf("\n");

            if (ii[0] == ii[3] && ii[1] == ii[2] && ii[0] != ii[1]) abba++;
            if (ii[0] == ii[2] && ii[1] == ii[3] && ii[0] != ii[1]) baba++;
          }
        }
      }
    }
  }
  printf("abba: %ld\n", abba);
  printf("baba: %ld\n", baba);
  printf("Total: %ld\n", total);
  printf("abba score: %f\n", abba/(double)total);
  printf("baba score: %f\n", baba/(double)total);
  return 0;
}

static msa_t * msa_collapse_to_quartet(msa_t * msa)
{
  long i,j;
  msa_t * newmsa = NULL;

  /* allocate new msa */
  newmsa = (msa_t *)xmalloc(sizeof(msa_t));
  
  newmsa->count = 4;
  //newmsa->length = XXX;

  /* find sites to be removed */
  for (i = 0; i < msa->length; i++)
  {
    for (j = 0; j < msa->count; ++j)
    {
    }
  }

  return newmsa;
}

#if 1
static char ** split4(const char * s)
{
  long i,k;
  long commas_count = 0;
  char ** taxa = NULL;

  /* compute number of commas in list of tips */
  for (i = 0; i < (long)strlen(s); ++i)
    if (s[i] == ',')
      ++commas_count;
  
  if (commas_count+1 != 4)
    fatal("ABBA-BABA test requires exactly four taxa");

  taxa = (char **)xmalloc((size_t)(commas_count+1) * sizeof(char *));

  k = 0;
  while (*s)
  {
    /* get next taxon */
    size_t taxon_len = strcspn(s,",");
    if (!taxon_len)
      fatal("Erroneous format in --dstat (taxon missing)");

    taxa[k++] = xstrndup(s, taxon_len);

    s += taxon_len;
    assert(*s == ',' || *s == '\0');
    if (*s == ',')
      ++s;
  }

  return taxa;
}
#endif

static void slice_condensed(msa_t * msa_src,
                            msa_t * msa_dst,
                            double * vec_src,
                            double * vec_dst,
                            long src_start,
                            long src_end,
                            long dst_start)
{
  long i,j;
  long length;
  long * indices;

  assert(msa_src->count == msa_dst->count);

  /* copy alignment */
  for (i = 0; i < msa_src->count; ++i)
  {
    memcpy(msa_dst->sequence[i]+dst_start,
           msa_src->sequence[i]+src_start,
           (src_end-src_start)*sizeof(char));
  }

  /* copy probability vectors */
  for (i = 0; i < msa_src->count; ++i)
  {
    memcpy(vec_dst+msa_src->count*NT_CHARS*dst_start,
           vec_src+msa_src->count*NT_CHARS*src_start,
           msa_src->count*NT_CHARS*(src_end-src_start)*sizeof(double));
  }

  length = src_end - src_start;

  /* resmaple sites */
  indices = (long *)xmalloc((size_t)length * sizeof(long));

  /* TF: 28.11.2023 */
  #if 1
  for (i = 0; i < length; ++i)
    indices[i] = (long)(rndu(0)*length);
  #else
  /* TODO: Delete, debug only */
  for (i = 0; i < length; ++i)
    indices[i] = i;
  #endif

  #if 1
  if (opt_debug)
  {
    for (i = 0; i < length; ++i) printf(" %ld", indices[i]);
    printf("\n");
  }
  #endif

  for (i = 0; i < length; ++i)
  {
    for (j = 0; j < msa_src->count; ++j)
      msa_dst->sequence[j][dst_start+i] = msa_src->sequence[j][src_start+indices[i]];

    memcpy(vec_dst+msa_dst->count*NT_CHARS*(dst_start+i),
           vec_src+msa_src->count*NT_CHARS*(src_start+indices[i]),
           msa_src->count*NT_CHARS*sizeof(double));
  }

  free(indices);
}

static double resample_condensed_save(msa_t ** msa_list,
                                      msa_t * ss,
                                      double * abba_vec,
                                      double * baba_vec,
                                      long count)
{
  long i,j;
  long indices;
  long * spos = NULL;
  long * rindex = NULL;
  double abba_count = 0;
  double baba_count = 0;

  rindex = (long *)xmalloc((size_t)count*sizeof(long));
  /* calculate starting position of each locus in condensed subsampled msa */
  spos = (long *)xmalloc((size_t)count*sizeof(long));
  spos[0] = 0;

  spos[0] = 0;
  for (i = 1; i < count; ++i)
    spos[i] = spos[i-1] + msa_list[i-1]->length;

  for (i = 0; i < count; ++i)
    rindex[i] = (long)(rndu(0)*count);

  for (i = 0; i < count; ++i)
  {
    long length = msa_list[rindex[i]]->length;

    for (j = 0; j < length; ++j)
    {
      indices = (long)(rndu(0)*length);

      abba_count += abba_vec[spos[rindex[i]] + indices];
      baba_count += baba_vec[spos[rindex[i]] + indices];
    }
  }

  free(rindex);
  free(spos);
  if (abba_count + baba_count == 0) return 0;
  return ((abba_count - baba_count) / (abba_count + baba_count));
}

static double * jackknife(msa_t ** msa_list,
                          msa_t * ss,
                          double * abba_vec,
                          double * baba_vec,
                          long count)
{
  long i,j,k;
  long * spos = NULL;
  double abba_count = 0;
  double baba_count = 0;
  double * D;

  /* calculate starting position of each locus in condensed subsampled msa */
  spos = (long *)xmalloc((size_t)count*sizeof(long));
  D = (double *)xmalloc((size_t)count*sizeof(double));
  spos[0] = 0;

  spos[0] = 0;
  for (i = 1; i < count; ++i)
    spos[i] = spos[i-1] + msa_list[i-1]->length;

  /* create count jacknife samples */
  for (i = 0; i < count; ++i)
  {
    /* one jackknife sample */
    abba_count = 0;
    baba_count = 0;
    for (j = 0; j < count; ++j)
    {
      if (i == j) continue;

      #if 0
      /* TF: 28.11.2023 */
      long dbg_index = (long)(rndu(0)*msa_list[j]->length);
      #endif

      for (k = 0; k < msa_list[j]->length; ++k)
      {
        #if 0
        /* TF: 28.11.2023 */
        if (k == dbg_index)
          continue;
        #endif

        abba_count += abba_vec[spos[j] + k];
        baba_count += baba_vec[spos[j] + k];
      }
    }
    if (abba_count + baba_count == 0)
      D[i] = 0;
    else
      D[i] = (abba_count - baba_count) / (abba_count + baba_count);
  }

  free(spos);

  return D;
}

static msa_t * resample_condensed(msa_t ** msa_list,
                                  msa_t * ss,
                                  double * vec,
                                  double **outvec,
                                  long count)
{
  long i;
  long * spos = NULL;
  long * rindex = NULL;
  msa_t * newmsa;
  double * newvec;

  newmsa = (msa_t *)xcalloc(1,sizeof(msa_t));
  rindex = (long *)xmalloc((size_t)count*sizeof(long));

  newmsa->count = ss->count;
  newmsa->length = 0;


  /* calculate starting position of each locus in condensed subsampled msa */
  spos = (long *)xmalloc((size_t)count*sizeof(long));
  spos[0] = 0;

  #if 0
  printf("spos[0] = %ld\n", spos[0]);
  #endif
  for (i = 1; i < count; ++i)
  {
    spos[i] = spos[i-1] + msa_list[i-1]->length;
    #if 0
    printf("spos[%ld] = %ld\n", i, spos[i]);
    #endif
  }

  /* loci indices resampled with replacement */
  #if 1
  if (opt_debug)
    xdebug_noendl("Resampled loci(sites): ");
  #endif
  for (i = 0; i < count; ++i)
  {
    rindex[i] = (long)(rndu(0)*count);
  #if 1
    if (opt_debug)
      printf(" %ld(%d)", rindex[i], msa_list[rindex[i]]->length);
  #endif
    newmsa->length += msa_list[rindex[i]]->length;
  }
  #if 1
  if (opt_debug)
    printf("\n");
  #endif
  newmsa->sequence = (char **)xmalloc((size_t)ss->count * sizeof(char *));
  newmsa->label    = (char **)xmalloc((size_t)ss->count * sizeof(char *));
  newvec = (double *)xmalloc(ss->count*NT_CHARS*newmsa->length*sizeof(double));

  #if 0
  xdebug("new alignment length: %d\n", newmsa->length);
  #endif
  for (i = 0; i < newmsa->count; ++i)
  {
    newmsa->label[i] = xstrdup(ss->label[i]);
    newmsa->sequence[i] = (char *)xcalloc((size_t)(newmsa->length+1),sizeof(char));
  }

  long newmsa_pos = 0;
  for (i = 0; i < count; ++i)
  {
    if (opt_debug)
      xdebug_noendl("Resampled site indices from locus %ld: ", rindex[i]);
    slice_condensed(ss,
                    newmsa,
                    vec,
                    newvec,
                    spos[rindex[i]],
                    spos[rindex[i]] + msa_list[rindex[i]]->length,
                    newmsa_pos);
    newmsa_pos += msa_list[rindex[i]]->length;
  }

  free(rindex);

  *outvec = newvec;

  free(spos);

  return newmsa;
}

#if 0
msa_t * phylip_concat(msa_t ** msa_list, long msa_count)
{
  long i,j,k,m;
  msa_t * msa = NULL;
  long total_length = 0;

  /* calculate total alignment length */
  for (i = 0; i < msa_count; ++i)
  {
    total_length += msa_list[i]->length;
  }

  msa = (msa_t *)xcalloc(1, sizeof(msa_t));
  msa->count  = 4;
  msa->length = total_length;
  msa->sequence = (char **)xmalloc(4*sizeof(char *));
  msa->label = (char **)xmalloc(4*sizeof(char *));
  for (i = 0; i < 4; ++i)
    msa->sequence[i] = (char *)xmalloc((size_t)(total_length+1) * sizeof(char));

  /* naive check of species labels */
  for (k=0, i = 0; i < msa_count; ++i)
  {
    if (msa_list[i]->count > 4)
      fatal("More than 4 sequences in alignment %ld", i);

    for (j = 0; j < msa_list[i]->count; ++j)
    {
      /* check that label is in the concenated alignment structure */
      for (m = 0; m < k; ++m)
        if (!strcmp(msa_list[i]->label[j],msa->label[m]))
          break;
      if (m == k)
      {
        if (k == 4)
          fatal("More than 4 sequences in full alignment");
        msa->label[k++] = xstrdup(msa_list[i]->label[j]);
      }
    }
  }

  if (k != 4)
    fatal("Error: only %ld sequences in alignments. Need 4 sequences.", k);

  /* create concatenated alignment */
  long offset = 0;
  for (i = 0; i < msa_count; ++i)
  {
    for (m = 0; m < k; ++m)
    {
      for (j = 0; j < msa_list[i]->count; ++j)
        if (!strcmp(msa_list[i]->label[j],msa->label[m]))
          break;

      if (j == msa_list[i]->count)
      {
        /* Sequence not found in current alignment, fill with missing data */
        memset(msa->sequence[m]+offset,'?',msa_list[i]->length);
      }
      else
      {
        /* Sequence found, copy data */
        memcpy(msa->sequence[m]+offset,
               msa_list[i]->sequence[j],
               msa_list[i]->length);
      }

    }
    offset += msa_list[i]->length;
  }
  assert(offset == total_length);

  for (i = 0; i < msa->count; ++i)
    msa->sequence[i][total_length] = 0;
  return msa;
}
#endif

void cmd_hyde()
{
  long i,j; //,k,m;
  long msa_count;
  phylip_t * fd;
  msa_t ** msa_list;

  if (opt_ci_alpha <= 0 || opt_ci_alpha >= 1)
    fatal("Confidence interval alpha must be in (0,1)");
  #if 0
  if (!opt_quiet)
    printf("Pre-computing table for site scores...\n");
  #endif
  /* this is probably redundant */
  precompute_table();



  #if 0
  char site[5] = "NRRN\0";
  unsigned int s[4] = {pll_map_nt[(int)site[0]],pll_map_nt[(int)site[1]],pll_map_nt[(int)site[2]],pll_map_nt[(int)site[3]]};
  printf("Decoding site: %s\n", site);
  decode_site(s);

  printf("Checking from precomputed table:\n");
  index = s[0] | (s[1] << 4) | (s[2] << 8) | (s[3] << 12);
  printf("abba: %f\n", abba_tbl[index]);
  printf("baba: %f\n", baba_tbl[index]);
  #endif

  /* open phylip file */
  fd = phylip_open(opt_msafile, pll_map_fasta);
  if (!fd)
    fatal("Cannot open file %s", opt_msafile);

  /* read alignment */
  msa_list = phylip_parse_multisequential(fd, &msa_count);
  assert(msa_list);

  #if 0
  for (i = 0; i < msa_count; ++i)
  {
    printf("MSA %ld\n", i);
    for (j = 0; j < msa_list[i]->count; ++j)
    {
      printf("  %s  %s\n", msa_list[i]->label[j], msa_list[i]->sequence[j]);
    }
  }
  #endif

  phylip_close(fd);

  /* reject JC69-compressed input (see cmd_dstat for the rationale) */
  for (i = 0; i < msa_count; ++i)
    if (msa_list[i]->compress_model == COMPRESS_JC69)
      fatal("--hyde does not support JC69 pattern-compressed alignments; "
            "use uncompressed input or GTR-compressed input");

  #if 0
  /* TODO: For now we only allow one alignment */
  assert(msa_count == 1);
  #endif

  if (!opt_mapfile)
    fatal("A map file needs to be specified with the --map option....");

  list_t * maplist = parse_mapfile(opt_mapfile);
  #if 0
  if (!opt_quiet)
  {
    printf("List of tag -> species mappings:\n");
    maplist_print(maplist);
  }
  #endif

  /* check that all tags are present in the map file */


  char ** taxa_list = split4(opt_hyde);

  #if 0
  printf("Tree: (((%s,%s),%s),%s);\n", taxa[0], taxa[1], taxa[2], taxa[3]);
  printf("Testing introgression between %s and %s, and between %s and %s\n",
         taxa[0], taxa[2], taxa[1], taxa[2]);
  #endif


  #if 0
  /* check for duplicate taxa */
  char ** labels = (char **)xcalloc(4,sizeof(char *));
  m = 0;
  for (i = 0; i < msa_count; ++i)
  {
    if (msa_list[i]->count > 4)
      fatal("More than 4 sequences found in alignment %ld.", i);
    for (j = 0; j < msa_list[i]->count; ++j)
    {
      for (k = 0; k < m; ++k)
        if (!strcmp(msa_list[i]->label[j],labels[k]))
          break;
      if (k == m)
      {
        if (m == 4)
          fatal("More than 4 sequences found in the dataset");

        labels[m++] = xstrdup(msa_list[i]->label[j]);
      }
    }
  }
  assert(m == 4);
  for (i = 0; i < 4; ++i)
    free(labels[i]);
  free(labels);
  #endif

  /* concatenate (possibly) multiple alignments and fill in missing data */
  #if 0
  msa_t * concat = phylip_concat(msa_list, msa_count);
  #else
  msa_t * concat = concatenate(msa_list, msa_count);
  #endif

  #if 0
  show_tags(concat);
  #endif

  if (opt_debug)
  {
    xdebug(ANSI_COLOR_RED "1. PHYLIP alignment loaded from %s" ANSI_COLOR_RESET, opt_msafile);
    phylip_print(stdout, concat);
  }

  long * index = NULL;
  long sp_count = 0;
  long * sp_seqcount = NULL;
  char ** sp_labels;
  sp_labels = map_sequences_to_species(concat, maplist, &index, &sp_seqcount, &sp_count);

  double * vec;
  msa_t * condmsa = condense_msa(concat, index, sp_labels, sp_count, &vec);

  #if 1
  if (opt_debug)
  {
    xdebug(ANSI_COLOR_RED "2. Condensed alignment:" ANSI_COLOR_RESET);
    phylip_print(stdout, condmsa);

    xdebug("Per site nucleotide character probability vectors for each species:");
    for (i = 0; i < sp_count; ++i)
    {
      xdebug("%s  ", sp_labels[i]);
      for (j = 0; j < condmsa->length; ++j)
      {
        double * v = vec + 4*(sp_count*j+i);
        xdebug(" %ld:([%f,%f,%f,%f)", j+1, v[0], v[1], v[2], v[3]);
      }
      printf("\n");
    }
  }
  #endif
  

  char * cilabel = NULL;
  if (opt_ci_alpha == 0.05)
    xasprintf(&cilabel, "95%% CI");
  else
    xasprintf(&cilabel, "%.1f%% CI", (1-opt_ci_alpha)*100);

  const char * colnames[] = {"D", "f(abba)", "f(baba)", cilabel, "tree"};
  long colsize[5] = {0,0,0,0,0};

  long t;
  #if 0
  long treestrsize = 9;  // 6 opar/cpar and 3 commas
  for (t = 0; t < 4; ++t)
    treestrsize += strlen(taxa_list[t]);
  printf("    D       f(abba)     f(baba)   p-value    tree\n");
  printf("=============================================");
  for (t = 0; t < treestrsize; ++t)
    printf("=");
  printf("\n");
  #endif

  char ** taxa = (char **)xmalloc((size_t)4*sizeof(char *));

  double * dilist = (double *)xmalloc((size_t)opt_bscount * sizeof(double));
  

  double * abba_vec = (double *)xmalloc((size_t)condmsa->length * sizeof(double));
  double * baba_vec = (double *)xmalloc((size_t)condmsa->length * sizeof(double));
  double * bbaa_vec = (double *)xmalloc((size_t)condmsa->length * sizeof(double));

  /* TF: 7/12/2023 
  changed to iterate 6 permutaions */
  #if 0
  for (t = 0; t < 24; ++t)
  #else
  for (t = 0; t < 6; ++t)
  #endif
  {
    taxa[0] = taxa_list[perms[t][0]];
    taxa[1] = taxa_list[perms[t][1]];
    taxa[2] = taxa_list[perms[t][2]];
    taxa[3] = taxa_list[perms[t][3]];

    /* subsampling */
    double * outvec = NULL;
    msa_t * ss = subsample_msa(condmsa,
                               taxa,
                               vec,
                               abba_vec,
                               baba_vec,
                               bbaa_vec,
                               &outvec,
                               4);

    if (opt_debug)
    {
      xdebug(ANSI_COLOR_RED "3. Filtered alignment for 4 species: %s" ANSI_COLOR_RESET, opt_hyde);
      phylip_print(stdout, ss);
    }

    #if 1
    if (opt_debug)
    {
      xdebug("Per site nucleotide character probability vectors for each species:");
      for (i = 0; i < ss->count; ++i)
      {
        xdebug("%s  ", ss->label[i]);
        for (j = 0; j < ss->length; ++j)
        {
          double * v = outvec + 4*(ss->count*j+i);
          xdebug(" %ld:([%f,%f,%f,%f)", j+1, v[0], v[1], v[2], v[3]);
        }
        //printf("\n");
      }
      xdebug("");
    }
    #endif

    #if 0
    printf("Calculating...\n");
    #endif
    double fabba = 0;
    double fbaba = 0;
    double fbbaa = 0;
    double dscore = calculate_d_vec(ss,outvec,&fabba, &fbaba, &fbbaa);

    #if defined(SAVINGS)
    double dbg_abba = 0;
    double dbg_baba = 0;
    double dbg_bbaa = 0;
    if (condmsa->pattern_weights)
    {
      for (i = 0; i < ss->length; ++i)
      {
        double w = (double)condmsa->pattern_weights[i];
        dbg_abba += abba_vec[i] * w;
        dbg_baba += baba_vec[i] * w;
        dbg_bbaa += bbaa_vec[i] * w;
      }
    }
    else
    {
      for (i = 0; i < ss->length; ++i)
      {
        dbg_abba += abba_vec[i];
        dbg_baba += baba_vec[i];
        dbg_bbaa += bbaa_vec[i];
      }
    }
    #if 1
    printf(ANSI_COLOR_RED "ABBA: %f    BABA: %f    BBAA: %f" ANSI_COLOR_RESET "\n", dbg_abba, dbg_baba, dbg_bbaa);
    #endif
    #endif

    #if 1
    if (t == 0)
    {
      colsize[0] = getnumdigits((long)dscore) + 1+1+1+6; // extra,sign,dot,floats
      colsize[1] = getnumdigits((long)fabba) + 4+1+6;  // extra,dot,floats
      colsize[2] = getnumdigits((long)fbaba) + 4+1+6;  // extra,dot,floats
      colsize[3] = 2*(getnumdigits((long)dscore)+2)+ 2+1+2+2+12; // extra,pars,comma,sign,dots,float
      colsize[4] = 10; //3 commas, 3 opars, 3 cpars, semicolon
      for (i = 0; i < 4; ++i)
        colsize[4] += strlen(taxa[i]);

      long totalsize = 0;
      for (i = 0; i < 5; ++i)
      {
        int lpadd = (int)((colsize[i] - strlen(colnames[i]))/2);
        int rpadd = (int)(colsize[i] - strlen(colnames[i]) - lpadd);
        printf("%*s%s%*s",
               lpadd, "",
               colnames[i],
               rpadd, "");
        if (i != 4)
        {
          printf(" ");
          totalsize++;
        }
        totalsize += colsize[i];
      }
      printf("\n");
      for (i = 0; i < totalsize; ++i)
      {
        printf("=");
      }
      printf("\n");

    }

    #endif

    if (opt_debug)
    {
      xdebug("");
      xdebug(ANSI_COLOR_RED "4. Resampling and bootstraping..." ANSI_COLOR_RESET);
      xdebug("");
    }
    exit(0);

    rnd_init();
    for (i = 0; i < opt_bscount; ++i)
    {
      #if defined(SAVINGS)
      dilist[i] = resample_condensed_save(msa_list, ss, abba_vec, baba_vec, msa_count);
      #else
      double * rsvec = NULL;
      msa_t * rsmsa = resample_condensed(msa_list, ss, outvec, &rsvec, msa_count);
      if (opt_debug)
      {
        phylip_print(stdout, rsmsa);
        debug_print_vec(rsvec, rsmsa);
      }
      double d = calculate_d_vec(rsmsa,rsvec, NULL, NULL);
      dilist[i] = d;
      msa_destroy(rsmsa);
      free(rsvec);
      #endif
    }

    qsort(dilist, opt_bscount, sizeof(double), cb_cmp_double);

    char * ci = NULL;
    double ci_lo = dilist[(long)(opt_bscount*(opt_ci_alpha/2))];
    double ci_hi = dilist[(long)(opt_bscount*(1-opt_ci_alpha/2))];

    /* TF 7/12/2023
      removed jackknife */
    #if 1
    double * djack = jackknife(msa_list, ss, abba_vec, baba_vec, msa_count);
    qsort(djack, msa_count, sizeof(double), cb_cmp_double);
    char * jack_ci = NULL;
    //double jack_ci_lo = djack[(long)(msa_count*(opt_ci_alpha/2))];
    double jack_ci_lo = djack[0];
    //double jack_ci_hi = djack[(long)(msa_count*(1-opt_ci_alpha/2))];
    double jack_ci_hi = djack[msa_count-1];
    xasprintf(&jack_ci, "(%.6f,%.6f)", jack_ci_lo, jack_ci_hi);

    xasprintf(&ci, "(%.6f,%.6f)", ci_lo, ci_hi);
    if (opt_ansi && (ci_lo > 0 || ci_hi < 0))
      printf(ANSI_COLOR_RED);

    /* TF 7/12/2023 removed jackknife from prinout */
    #if 0
    printf("%*.6f %*.6f %*.6f %*s (((%s,%s),%s),%s); %s\n", 
           (int)colsize[0], dscore,
           (int)colsize[1], fabba,
           (int)colsize[2], fbaba,
           (int)colsize[3], ci,
           taxa[0],taxa[1],taxa[2],taxa[3], jack_ci);
    #else
    printf("%*.6f %*.6f %*.6f %*s (((%s,%s),%s),%s);\n", 
           (int)colsize[0], dscore,
           (int)colsize[1], fabba,
           (int)colsize[2], fbaba,
           (int)colsize[3], ci,
           taxa[0],taxa[1],taxa[2],taxa[3]);
    #endif
    if (opt_ansi && (ci_lo > 0 || ci_hi < 0))
      printf(ANSI_COLOR_RESET);
    free(ci);

    free(jack_ci);
    #endif

    free(outvec);
    msa_destroy(ss);

    #if 0
    /* TF: 28.11.2023 */
    /* bootstrap variance */
    double di_bs_mean = 0;
    double di_bs_var = 0;
    for (i = 0; i < opt_bscount; ++i)
      di_bs_mean += dilist[i];
    di_bs_mean /= opt_bscount;

    for (i = 0; i < opt_bscount; ++i)
      di_bs_var += (dilist[i] - di_bs_mean)*(dilist[i] - di_bs_mean);
    di_bs_var /= (opt_bscount-1);

    /* jackknife variance */
    double bi_jk_mean = 0;
    double bi_jk_var = 0;
    for (i = 0; i < msa_count; ++i)
      bi_jk_mean += djack[i];
    bi_jk_mean /= msa_count;

    for (i = 0; i < msa_count; ++i)
      bi_jk_var += (djack[i] - bi_jk_mean)*(djack[i] - bi_jk_mean);
    di_bs_var /= (msa_count-1);

    printf("BS JK stdev: %f %f\n", sqrt(di_bs_var), sqrt(bi_jk_var));
    #endif
         

    /* TF: 7/12/2023
      removed jackknife */
    #if 1
    free(djack);
    //break;
    #endif
  }
  free(dilist);


  for (i = 0; i < sp_count; ++i)
    free(sp_labels[i]);
  free(sp_labels);
  free(index);
  free(sp_seqcount);

  //print_score_table();

  #if 0
  calculate_d(concat);

  printf("\nResampling and bootstraping...\n");
  rnd_init();

  #if 1
  msa_t ** resamp_msa_list = resample_msa_list(msa_list, msa_count);

  for (i = 0; i < msa_count; ++i)
  {
    resamp_msa_list[i] = resample_sites(resamp_msa_list[i]);
    
    msa_t * resamp_concat = concatenate(resamp_msa_list, msa_count);
    phylip_print(stdout, resamp_concat);
    
  }
  #endif

  //for (i = 0; i < 30; ++i)
  for (i = 0; i < 5; ++i)
  {
    msa_t * repl = resample_sites(concat);
    calculate_d(repl);
    msa_destroy(repl);
  }
  #endif

  rnd_fini();

  msa_destroy(concat);
  msa_destroy(condmsa);

  for (i = 0; i < 4; ++i)
    free(taxa_list[i]);
  free(taxa_list);
  free(taxa);
  free(cilabel);

  #if 0
  #if 1
  unsigned int dbg_site[4] = {pll_map_nt[(int)'V'],
                              pll_map_nt[(int)'W'],
                              pll_map_nt[(int)'W'],
                              pll_map_nt[(int)'N']};
  printf("DEBUG SITE:\n");
  debug_decode_site(dbg_site);
  #endif
  #endif

  for (i = 0; i < msa_count; ++i)
    msa_destroy(msa_list[i]);
  free(msa_list);

  if (abba_vec)
    free(abba_vec);
  if (baba_vec)
    free(baba_vec);
  if (abba_tbl)
    free(abba_tbl);
  if (baba_tbl)
    free(baba_tbl);
  if (patt_count)
    free(patt_count);

  list_clear(maplist, map_dealloc);
  free(maplist);
  free(vec);

  //assert(0);
  #if 0
  free(resamp_msa_list);
  #endif

}
