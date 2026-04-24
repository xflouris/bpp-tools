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
   {  0,   0,   0,   0},
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

  sitevec = (double *)xmalloc((size_t)species_count*4*sizeof(double));

  /* calculate amount of space needed for condensed alignment */
  for (i = 0; i < msa->length; ++i)
  {
    /* reset character */
    for (j = 0; j < species_count; ++j)
      sp_allele_code[j] = 0;

    for (j = 0; j < msa->count; ++j)
    {
      sp_index = sp_seqassign[j]; 

      /* upto biallelic characters */
      if (map_solid_and_biallelic[(int)msa->sequence[j][i]])
      {
        /* update allele character for species sp_index */
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

static void debug_print_vec(double * vec, msa_t * msa)
{
  long i,j;
  long species;

  species = msa->count;

  xdebug("Per site nucleotide character probability vectors for each species:");
  for (i = 0; i < species; ++i)
  {
    xdebug("%s  ", msa->label[i]);
    for (j = 0; j < msa->length; ++j)
    {
      double * v = vec + 4*(species*j+i);
      xdebug(" %ld:[%f,%f,%f,%f]", j+1, v[0], v[1], v[2], v[3]);
    }
    printf("\n");
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

  span = NT_CHARS*labels_count;
  v = (double *)xmalloc((size_t)(span*msa->length)*sizeof(double));

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
    printf("abba_vec:");
    printf("  [%.6f", abba_vec[0]);
    for (i = 1; i < msa->length; ++i)
    {
      printf(",%.6f", abba_vec[i]);
    }
    printf("]\n");

    /* baba */
    printf("baba_vec:");
    printf("  [%.6f", baba_vec[0]);
    for (i = 1; i < msa->length; ++i)
    {
      printf(",%.6f", baba_vec[i]);
    }
    printf("]\n");

    /* bbaa */
    printf("bbaa_vec:");
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

static long within_clade(rnode_t * node, rnode_t * clade_root)
{
  rnode_t * p = node;

  while (p && p != clade_root)
    p = p->parent;

  return !!p;
}

static int cb_cmp_double_asc(const void * a, const void * b)
{
  double * x = (double *)a;
  double * y = (double *)b;

  if ( *x > *y) return 1;
  if ( *x < *y) return -1;
  return 0;
}

static double median_A(msa_t * condmsa,
                       rnode_t ** a_tips,
                       long a_tip_count,
                       rnode_t ** b_tips,
                       long b_tip_count,
                       rnode_t * C,
                       rnode_t * O,
                       double * vec,
                       double *abba_vec,
                       double *baba_vec,
                       double *bbaa_vec)
{
  long i,j,k;
  double mgamma;
  double abba = 0;
  double baba = 0;
  double bbaa = 0;
  char * taxa[4];
  double * vgamma;

  /* vector for storing per-A minimum gammas */
  vgamma = (double *)xcalloc((size_t)a_tip_count,sizeof(double));

  /* go through each A population */
  for (i = 0; i < a_tip_count; ++i)
  {
    rnode_t * A = a_tips[i];
    double min_gamma = 0;
    int first = 1;

    /* find minimum gamma over all B tips for this A */
    for (j = 0; j < b_tip_count; ++j)
    {
      rnode_t * B = b_tips[j];

      /* build quartet */
      taxa[0] = A->label;
      taxa[1] = B->label;
      taxa[2] = C->label;
      taxa[3] = O->label;

      /* extract alignment and pattern counts for quartet */
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
        xdebug("Quartet: (((%s,%s),%s),%s)",
               taxa[0],taxa[1],taxa[2],taxa[3]);
      }

      /* compute pattern counts (weighted when input is pattern-compressed;
         weights come from condmsa since ss has the same sites/patterns) */
      abba = baba = bbaa = 0;
      if (condmsa->pattern_weights)
      {
        for (k = 0; k < ss->length; ++k)
        {
          double w = (double)condmsa->pattern_weights[k];
          abba += abba_vec[k] * w;
          baba += baba_vec[k] * w;
          bbaa += bbaa_vec[k] * w;
        }
      }
      else
      {
        for (k = 0; k < ss->length; ++k)
        {
          abba += abba_vec[k];
          baba += baba_vec[k];
          bbaa += bbaa_vec[k];
        }
      }

      /* calculate gamma with guard */
      double gamma = (bbaa <= baba || abba <= baba) ?
                       0 : (bbaa - baba)/(bbaa - 2*baba + abba);

      /* track minimum over B */
      if (first || gamma < min_gamma)
      {
        min_gamma = gamma;
        first = 0;
      }

      /* dealloc */
      msa_destroy(ss);
      free(outvec);
    }

    vgamma[i] = min_gamma;
  }

  /* sort gammas in ascending order */
  if (a_tip_count>1)
    qsort(vgamma, a_tip_count,sizeof(double),cb_cmp_double_asc);

  assert(a_tip_count);

  /* get gamma median value */
  mgamma = (a_tip_count % 2) ?
    vgamma[a_tip_count/2] : (vgamma[a_tip_count/2-1]+vgamma[a_tip_count/2])/2;

  if (opt_debug)
  {
    printf("    vgamma: [ %f", vgamma[0]);
    for (i = 1; i < a_tip_count; ++i)
      printf(",%f",vgamma[i]);
    printf("] gamma: %f\n", mgamma);
  }

  free(vgamma);
  return mgamma;
}

static void fill_tips(rnode_t * root, rnode_t ** outvec, long * index)
{
  if (!root->left)
  {
    /* tip node */
    assert(!root->right);

    outvec[*index] = root;
    (*index)++;
    return;
  }

  assert(root->left && root->right);
  fill_tips(root->left,  outvec, index);
  fill_tips(root->right, outvec, index);
}

static char * center(const char * s, int space)
{
  char * r;
  int len = (int)strlen(s);
  int left,right;

  left  = (space-len)/2;
  right = space-len-left;

  xasprintf(&r, "%*s%s%*s", left, "", s, right, "");

  return r;
}

static void print_matrix(double ** m, long r, long c, rnode_t * outgroup,
                         rtree_t * rtree)
{
  long i,j;
  int * digits;
  int prec = 3;
  char na[] = "N/A";

  /* calculate max number of digits per column (for correct alignment) */
  digits = (int *)xcalloc((size_t)c,sizeof(int));
  for (i = 0; i < c; ++i)
    digits[i] = MAX((int)floor(log10(i+1)+1),strlen(rtree->nodes[i]->label));

  assert(r > 0 && c > 0);

  for (i = 0; i < c; ++i)
  {
    for (j = 0; j < r; ++j)
    {
      if (m[j][i] == -1)
        digits[i] = MAX(digits[i],strlen(na));
      else
      {
        int d = (int)floor(log10(fabs(m[j][i])+1)+1);
        d += prec+1;
        if (m[j][i] < 0) d++;
        digits[i] = MAX(digits[i],d);
      }
    }
  }

  char * x = NULL;
  printf("    ");
  if (rtree->nodes[0] != outgroup)
  {
    x = center(rtree->nodes[0]->label,digits[0]);
    printf("%s",x);
    free(x);
  }
  for (i = 1; i < c; ++i)
  {
    /* skip outgroup column */
    if (rtree->nodes[i] == outgroup) continue;

    x = center(rtree->nodes[i]->label,digits[i]);
    printf("  %s",x);
    free(x);
  }
  printf("\n");

  for (i = 0; i < r; ++i)
  {
    /* skip outgroup row */
    if (rtree->nodes[i]->label && rtree->nodes[i] == outgroup) continue;

    /* skip these 4 branches */
    if (!rtree->nodes[i]->parent) continue;
    if (rtree->nodes[i]->parent)
    {
      if (!rtree->nodes[i]->parent->parent) continue;
      else if (!rtree->nodes[i]->parent->parent->parent)
        continue;
    }

    printf("%2ld  ", i);
    if (rtree->nodes[0] != outgroup)
    {
      if (m[i][0] == -1)
        printf("%*s",digits[0],na);
      else
      {
        xasprintf(&x,"%.*f",prec,m[i][0]);
        printf("%*s", digits[0],x);
        free(x);
      }
    }
    for (j = 1; j < c; ++j)
    {
      if (rtree->nodes[j] == outgroup) continue;

      if (m[i][j] == -1)
        printf("  %*s",digits[j],na);
      else
      {
        xasprintf(&x,"%.*f",prec,m[i][j]);
        printf("  %*s", digits[j],x);
        free(x);
      }
    }
    printf("\n");
  }

  free(digits);
}

static void set_height(rnode_t * root)
{
  if (!root->left)
  {
    /* tip */
    assert(!root->right);
    root->height = 1;
    return;
  }

  set_height(root->left);
  set_height(root->right);

  root->height = MAX(root->left->height, root->right->height)+1;
}

static void set_tau_uniformly(rtree_t * rtree)
{
  long i;

  long maxheight = rtree->root->height;
  rtree->root->tau = 1;
  for (i = 0; i < rtree->tip_count+rtree->inner_count; ++i)
  {
    rnode_t * x = rtree->nodes[i];
    if (!x->parent) continue;

    x->tau = x->height*(rtree->root->tau/(double)maxheight);
  }
}

static long is_emptyline(const char * line)
{
  size_t ws = strspn(line, " \t\r\n");
  if (!line[ws] || line[ws] == '*' || line[ws] == '#') return 1;
  return 0;
}

void cmd_fbranch()
{
  long i,j;
  long msa_count;
  char * newick;
  FILE * fp_tree;
  phylip_t * fd;
  msa_t ** msa_list;

  if (!opt_msafile)
    fatal("Specify alignment using --msafile");
  if (!opt_outgroup)
    fatal("Specify outgroup using --outgroup");

  /* open phylip file */
  fd = phylip_open(opt_msafile, pll_map_fasta);
  if (!fd)
    fatal("Cannot open file %s", opt_msafile);

  /* read alignment */
  msa_list = phylip_parse_multisequential(fd, &msa_count);
  assert(msa_list);
  phylip_close(fd);

  /* reject JC69-compressed input (see cmd_dstat for the rationale) */
  for (i = 0; i < msa_count; ++i)
    if (msa_list[i]->compress_model == COMPRESS_JC69)
      fatal("--fbranch does not support JC69 pattern-compressed alignments; "
            "use uncompressed input or GTR-compressed input");

  if (!opt_mapfile)
    fatal("A map file needs to be specified with the --map option....");

  list_t * maplist = parse_mapfile(opt_mapfile);

  /* concatenate (possibly) multiple alignments and fill in missing data */
  msa_t * concat = concatenate(msa_list, msa_count);

  if (opt_debug)
  {
    xdebug(ANSI_COLOR_RED "1. PHYLIP alignment loaded from %s" ANSI_COLOR_RESET,
           opt_msafile);
    phylip_print(stdout, concat);
  }

  /* deallocate individual MSAs */
  for (i = 0; i < msa_count; ++i)
    msa_destroy(msa_list[i]);
  free(msa_list);

  /* map sequences to species. Returns "index" which contains one integer per sequence
     sequence indicating its species. sp_seqcount contains the number of sequences for
     each species */
  long * index = NULL;
  long sp_count = 0;
  long * sp_seqcount = NULL;
  char ** sp_labels;
  sp_labels = map_sequences_to_species(concat,
                                       maplist,
                                       &index,
                                       &sp_seqcount,
                                       &sp_count);

  /* now create a "condensed" MSA where all sequences inside a population
     are represented by just one sequence. Ambiguous characters are used to
     represent polymorphism. We also construct a vector (vec) which is filled
     with per-site nucleodite weights (size = msa->length * 4). */
  double * vec;
  msa_t * condmsa = condense_msa(concat, index, sp_labels, sp_count, &vec);
  free(index);
  free(sp_seqcount);

  #if 1
  if (opt_debug)
  {
    xdebug(ANSI_COLOR_RED "2. Condensed alignment:" ANSI_COLOR_RESET);
    phylip_print(stdout, condmsa);
    debug_print_vec(vec,condmsa);
  }
  #endif

  double * abba_vec = (double *)xmalloc((size_t)condmsa->length*sizeof(double));
  double * baba_vec = (double *)xmalloc((size_t)condmsa->length*sizeof(double));
  double * bbaa_vec = (double *)xmalloc((size_t)condmsa->length*sizeof(double));


  fp_tree = xopen(opt_treefile, "r");
  long lineno = 0;

  /* storage space for tip nodes within the a and b clades */
  rnode_t ** a_tips = NULL;
  rnode_t ** b_tips = NULL;
  long ab_size = 0;

  /* go through all trees in tree file */
  while ((newick = getnextline(fp_tree)))
  {
    lineno++;

    /* skip empty lines or comments */
    if (is_emptyline(newick))
    {
      free(newick);
      continue;
    }

    rtree_t * rtree = bpp_parse_newick_string(newick); 
    if (!rtree)
      fatal("Cannot parse tree (%s:%ld)", opt_treefile, lineno);

    printf("Processing tree: %s\n", newick);

    /* reallocate space for A and B tips if insufficient */
    if (ab_size < rtree->tip_count)
    {
      if (a_tips) free(a_tips);
      if (b_tips) free(b_tips);
      a_tips = (rnode_t **)xmalloc((size_t)rtree->tip_count * sizeof(rnode_t *));
      b_tips = (rnode_t **)xmalloc((size_t)rtree->tip_count * sizeof(rnode_t *));
      ab_size = rtree->tip_count;
    }

    if (rtree->tip_count != sp_count)
      fatal("Mismatching number of species in MSA (%ld) and tree (%d)",sp_count,rtree->tip_count);

    /* link each tip node with its sequence
       TODO: Can be done more efficiently */
    for (i = 0; i < rtree->tip_count; ++i)
    {
      for (j = 0; j < sp_count; ++j)
      {
        if (!strcmp(rtree->nodes[i]->label,sp_labels[j]))
        {
          rtree->nodes[i]->seq_index = j;
          break;
        }
      }
      assert(j < sp_count);
    }

    
    for (i = 0; i < rtree->tip_count; ++i)
    {
      if (!strcmp(opt_outgroup, rtree->nodes[i]->label))
        break;
    }
    if (i == rtree->tip_count)
      fatal("Outgroup not found...");

    rnode_t * O = rtree->nodes[i];

    assert(O->parent);
    if (O->parent->parent)
      fatal("Outgroup must be a tip and direct descendent of root");


    double ** matrix;
    /* get clades */
    rnode_t * b;
    rnode_t * a;
    rnode_t * C;
    long total_nodes = rtree->tip_count+rtree->inner_count;
    matrix = (double **)xmalloc((size_t)total_nodes*sizeof(double *));
    for (i = 0; i < total_nodes; ++i)
    {
      matrix[i] = (double *)xmalloc((size_t)rtree->tip_count*sizeof(double));
      for (j = 0; j < rtree->tip_count; ++j)
        matrix[i][j] = -1;
    }

    double mgamma = 0;
    for (i = 0; i < rtree->tip_count+rtree->inner_count; ++i)
    {
      b = rtree->nodes[i];

      if (!b->parent) continue;
      if (b->parent && !b->parent->parent) continue;

      /* branch b is parental to node b */

      /* branch a is b's sibling */
      a = (b->parent->left == b) ? b->parent->right : b->parent->left;

      for (j = 0; j < rtree->tip_count; ++j)
      {
        /* skip outgroup */
        if (rtree->nodes[j] == O) continue;

        C = rtree->nodes[j];

        if (within_clade(C,a)) continue;
        if (within_clade(C,b)) continue;

        long a_tip_count = 0;
        long b_tip_count = 0;
        fill_tips(a,a_tips,&a_tip_count);
        fill_tips(b,b_tips,&b_tip_count);


        /* generate triplet */
        if (opt_debug)
        {
          char * astr = rtree_export_newick(a,NULL);
          char * bstr = rtree_export_newick(b,NULL);
          xdebug("Working on tree: ((%s,%s),%s);", astr,bstr,C->label);
          free(astr);
          free(bstr);
        }


        mgamma = median_A(condmsa,
                          a_tips,
                          a_tip_count,
                          b_tips,
                          b_tip_count,
                          C,
                          O,
                          vec,
                          abba_vec,
                          baba_vec,
                          bbaa_vec);
        matrix[i][j] = mgamma;
      
      }
    }

    printf("Matrix for tree on line %ld:\n", lineno);
    print_matrix(matrix,total_nodes,rtree->tip_count,O,rtree);
    for (i = 0; i < total_nodes; ++i)
      free(matrix[i]);
    free(matrix);


    char * pdf_filename;
    xasprintf(&pdf_filename, "tree.%ld.pdf", lineno);
    set_height(rtree->root);
    set_tau_uniformly(rtree);

    rtree_export_pdf(rtree, pdf_filename);
    free(pdf_filename);
    #if 0
    char * treestr = rtree_export_newick(rtree->root,NULL);
    printf("Newick: %s\n", treestr);
    free(treestr);
    #endif
    rtree_destroy(rtree,NULL);

    free(newick);
  }
  for (i = 0; i < sp_count; ++i)
    free(sp_labels[i]);
  free(sp_labels);
  if (a_tips) free(a_tips);
  if (b_tips) free(b_tips);

  list_clear(maplist, map_dealloc);
  free(maplist);

  fclose(fp_tree);
  msa_destroy(concat);
  msa_destroy(condmsa);
  free(vec);
  free(abba_vec);
  free(baba_vec);
  free(bbaa_vec);
}
