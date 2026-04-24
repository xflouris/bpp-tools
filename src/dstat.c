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
   {  0,   0,   0,   0}
 };

static long perms[6][4] =
 {
   {0,1,2,3},
   {0,2,1,3},
   {1,0,2,3},
   {1,2,0,3},
   {2,0,1,3},
   {2,1,0,3}
 };

static long getnumdigits(long n)
{
  if (n == 0) return 1;
  return (long)(floor(log10(labs(n))+1));
}

static int cb_cmp_double(const void * a, const void * b)
{
  double * x = (double *)a;
  double * y = (double *)b;

  if (*x > *y) return 1;
  if (*x < *y) return -1;

  return 0;
}

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
      for (j = 0; j < i; ++j)
        free(labels[j]);
      free(labels);
      return NULL;
    }
  }

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

static msa_t * condense_msa(msa_t * msa,
                            long * sp_seqassign,
                            char ** sp_labels,
                            long species_count,
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

  sitevec = (double *)xmalloc((size_t)species_count*4*sizeof(double));

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

        /* update allele character for species sp_index */
        sp_allele_code[sp_index] |= pll_map_nt[(int)msa->sequence[j][i]];

        sp_seqcount[i*species_count+sp_index]++;
      }
    }

    for (j = 0; j < species_count; ++j)
    {
      newmsa->sequence[j][m] = invcharmap_acgt[sp_allele_code[j]];
      for (k = 0; k < 4; ++k)
      {
        pvec[k] = sitevec[4*j+k];
      }
      pvec += 4;
    }

    ++m;
  }


  free(sitevec);

  *vecptr = vec;

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

static msa_t * subsample_msa(msa_t * msa,
                             char ** labels,
                             double * vec,
                             double * abba_vec,
                             double * baba_vec,
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

  /* pre-calculate per-site per-species abba-baba scores */
  #if defined(SAVINGS)
  double * sitevec = v;
  assert(newmsa->count == 4);
  memset(abba_vec,0,msa->length*sizeof(double));
  memset(baba_vec,0,msa->length*sizeof(double));
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
      }
    }
    sitevec += 16;
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

/* Build a prefix-sum array from pattern weights.  cum[j] is the total
   number of original sites represented by patterns [0..j], i.e.
   cum[j] = w[0] + w[1] + ... + w[j].  cum[n-1] is the original alignment
   length.  Caller owns the returned array. */
static unsigned long * build_cum_weights(const unsigned int * w, long n)
{
  long i;
  unsigned long s = 0;
  unsigned long * cum = (unsigned long *)xmalloc((size_t)n *
                                                 sizeof(unsigned long));
  for (i = 0; i < n; ++i)
  {
    s += w[i];
    cum[i] = s;
  }
  return cum;
}

/* Given a cumulative-weight array cum[0..n-1] and a uniform random integer
   r in [0, cum[n-1]-1], return the pattern index whose weight "owns" that
   original-site index.  That is, the smallest j such that r < cum[j]. */
static long cum_weight_lookup(const unsigned long * cum, long n,
                              unsigned long r)
{
  long lo = 0, hi = n - 1;
  while (lo < hi)
  {
    long mid = lo + (hi - lo) / 2;
    if (r < cum[mid]) hi = mid;
    else               lo = mid + 1;
  }
  return lo;
}

static double resample_condensed_save(msa_t ** msa_list,
                                      msa_t * ss,
                                      double * abba_vec,
                                      double * baba_vec,
                                      long count,
                                      unsigned long ** cum_weights)
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
  for (i = 1; i < count; ++i)
    spos[i] = spos[i-1] + msa_list[i-1]->length;

  for (i = 0; i < count; ++i)
    rindex[i] = (long)(rndu(0)*count);

  for (i = 0; i < count; ++i)
  {
    msa_t * m = msa_list[rindex[i]];
    long pat_count = m->length;
    long draws;
    unsigned long total;

    if (m->pattern_weights)
    {
      /* pattern-compressed: draw 'total' original sites with replacement,
         each draw produces a pattern index via weighted lookup. This is
         statistically equivalent to sampling sites from the uncompressed
         alignment. */
      total = cum_weights[rindex[i]][pat_count - 1];
      draws = (long)total;
    }
    else
    {
      total = (unsigned long)pat_count;
      draws = pat_count;
    }

    for (j = 0; j < draws; ++j)
    {
      if (m->pattern_weights)
      {
        unsigned long r = (unsigned long)(rndu(0) * (double)total);
        if (r >= total) r = total - 1;  /* clamp: rndu can return 1.0 */
        indices = cum_weight_lookup(cum_weights[rindex[i]], pat_count, r);
      }
      else
      {
        indices = (long)(rndu(0) * (double)pat_count);
        if (indices >= pat_count) indices = pat_count - 1;
      }

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

      if (msa_list[j]->pattern_weights)
      {
        for (k = 0; k < msa_list[j]->length; ++k)
        {
          double w = (double)msa_list[j]->pattern_weights[k];
          abba_count += abba_vec[spos[j] + k] * w;
          baba_count += baba_vec[spos[j] + k] * w;
        }
      }
      else
      {
        for (k = 0; k < msa_list[j]->length; ++k)
        {
          abba_count += abba_vec[spos[j] + k];
          baba_count += baba_vec[spos[j] + k];
        }
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

void cmd_dstat()
{
  long i,j;
  long msa_count;
  phylip_t * fd;
  msa_t ** msa_list;

  if (opt_ci_alpha <= 0 || opt_ci_alpha >= 1)
    fatal("Confidence interval alpha must be in (0,1)");

  /* open phylip file */
  fd = phylip_open(opt_msafile, pll_map_fasta);
  if (!fd)
    fatal("Cannot open file %s", opt_msafile);

  /* read alignment */
  msa_list = phylip_parse_multisequential(fd, &msa_count);
  assert(msa_list);

  phylip_close(fd);

  /* JC69 pattern-compression re-encodes alleles site-locally, which collides
     with the specific-allele assumptions of downstream pipelines; the
     frequency-based ABBA/BABA formula here is actually invariant under such
     relabeling, but the argument is subtle (and ambiguity sites bypass the
     JC69 path entirely, mixing two encodings in one buffer). Rather than
     silently produce results whose correctness depends on an invariant that
     future code changes might break, reject JC69 input explicitly. GTR
     compression is fine. */
  for (i = 0; i < msa_count; ++i)
    if (msa_list[i]->compress_model == COMPRESS_JC69)
      fatal("--dstat does not support JC69 pattern-compressed alignments; "
            "use uncompressed input or GTR-compressed input");

  if (!opt_mapfile)
    fatal("A map file needs to be specified with the --map option....");

  list_t * maplist = parse_mapfile(opt_mapfile);

  /* check that all tags are present in the map file */


  char ** taxa_list = split4(opt_dstat);

  /* concatenate (possibly) multiple alignments and fill in missing data */
  msa_t * concat = concatenate(msa_list, msa_count);

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
  if (!sp_labels)
    fatal("Failed to map all sequences to species. "
          "Check that the map file contains entries for all individuals.");

  double * vec;
  msa_t * condmsa = condense_msa(concat, index, sp_labels, sp_count, &vec);

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


  char * cilabel = NULL;
  if (opt_ci_alpha == 0.05)
    xasprintf(&cilabel, "95%% CI");
  else
    xasprintf(&cilabel, "%.1f%% CI", (1-opt_ci_alpha)*100);

  const char * colnames[] = {"D", "f(abba)", "f(baba)", cilabel, "tree"};
  long colsize[5] = {0,0,0,0,0};

  long t;

  char ** taxa = (char **)xmalloc((size_t)4*sizeof(char *));

  double * dilist = (double *)xmalloc((size_t)opt_bscount * sizeof(double));


  double * abba_vec = (double *)xmalloc((size_t)condmsa->length * sizeof(double));
  double * baba_vec = (double *)xmalloc((size_t)condmsa->length * sizeof(double));

  /* For pattern-compressed input, precompute cumulative weight arrays so
     the bootstrap can draw one site at a time by binary-searching the
     cumulative weights. One array per locus (only msa_list[0] matters
     today since concatenate() rejects multi-locus compressed input). */
  unsigned long ** cum_weights = NULL;
  if (condmsa->pattern_weights)
  {
    cum_weights = (unsigned long **)xcalloc((size_t)msa_count,
                                            sizeof(unsigned long *));
    for (i = 0; i < msa_count; ++i)
      if (msa_list[i]->pattern_weights)
        cum_weights[i] = build_cum_weights(msa_list[i]->pattern_weights,
                                           msa_list[i]->length);
  }

  /* TF: 7/12/2023
  changed to iterate 6 permutaions */
  for (t = 0; t < 6; ++t)
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
                               &outvec,
                               4);

    if (opt_debug)
    {
      xdebug(ANSI_COLOR_RED "3. Filtered alignment for 4 species: %s" ANSI_COLOR_RESET, opt_dstat);
      phylip_print(stdout, ss);
    }

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
      }
      xdebug("");
    }

    /* weighted sum when input is pattern-compressed (weights come from
       condmsa; ss has the same sites/patterns) */
    double fabba = 0;
    double fbaba = 0;
    if (condmsa->pattern_weights)
    {
      for (i = 0; i < ss->length; ++i)
      {
        double w = (double)condmsa->pattern_weights[i];
        fabba += abba_vec[i] * w;
        fbaba += baba_vec[i] * w;
      }
    }
    else
    {
      for (i = 0; i < ss->length; ++i)
      {
        fabba += abba_vec[i];
        fbaba += baba_vec[i];
      }
    }
    double dscore = (fabba + fbaba != 0) ? (fabba - fbaba) / (fabba + fbaba) : 0;

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

    if (opt_debug)
    {
      xdebug("");
      xdebug(ANSI_COLOR_RED "4. Resampling and bootstraping..." ANSI_COLOR_RESET);
      xdebug("");
    }

    char * ci = NULL;
    char * jack_ci = NULL;
    double ci_lo = 0, ci_hi = 0;
    double * djack = NULL;

    rnd_init();
    for (i = 0; i < opt_bscount; ++i)
    {
      dilist[i] = resample_condensed_save(msa_list, ss, abba_vec,
                                          baba_vec, msa_count,
                                          cum_weights);
    }

    qsort(dilist, opt_bscount, sizeof(double), cb_cmp_double);

    ci_lo = dilist[(long)(opt_bscount*(opt_ci_alpha/2))];
    ci_hi = dilist[(long)(opt_bscount*(1-opt_ci_alpha/2))];

    djack = jackknife(msa_list, ss, abba_vec, baba_vec, msa_count);
    qsort(djack, msa_count, sizeof(double), cb_cmp_double);
    double jack_ci_lo = djack[0];
    double jack_ci_hi = djack[msa_count-1];
    xasprintf(&jack_ci, "(%.6f,%.6f)", jack_ci_lo, jack_ci_hi);

    xasprintf(&ci, "(%.6f,%.6f)", ci_lo, ci_hi);
    if (opt_ansi && (ci_lo > 0 || ci_hi < 0))
      printf(ANSI_COLOR_RED);

    printf("%*.6f %*.6f %*.6f %*s (((%s,%s),%s),%s);\n",
           (int)colsize[0], dscore,
           (int)colsize[1], fabba,
           (int)colsize[2], fbaba,
           (int)colsize[3], ci,
           taxa[0],taxa[1],taxa[2],taxa[3]);
    if (opt_ansi && (ci_lo > 0 || ci_hi < 0))
      printf(ANSI_COLOR_RESET);
    free(ci);

    free(jack_ci);

    free(outvec);
    msa_destroy(ss);

    if (djack) free(djack);
  }
  free(dilist);


  for (i = 0; i < sp_count; ++i)
    free(sp_labels[i]);
  free(sp_labels);
  free(index);
  free(sp_seqcount);

  rnd_fini();

  msa_destroy(concat);
  msa_destroy(condmsa);

  for (i = 0; i < 4; ++i)
    free(taxa_list[i]);
  free(taxa_list);
  free(taxa);
  free(cilabel);

  for (i = 0; i < msa_count; ++i)
    msa_destroy(msa_list[i]);
  free(msa_list);

  free(abba_vec);
  free(baba_vec);

  if (cum_weights)
  {
    for (i = 0; i < msa_count; ++i)
      if (cum_weights[i]) free(cum_weights[i]);
    free(cum_weights);
  }

  list_clear(maplist, map_dealloc);
  free(maplist);
  free(vec);
}
