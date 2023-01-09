/*
    Copyright (C) 2016-2023 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

const unsigned int bpp_nt_normal[256] =
{
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, '-',   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, '?',
   0, 'A', 'B', 'C', 'D',   0,   0, 'G', 'H',   0,   0, 'K',   0, 'M', 'N', 'O',
   0,   0, 'R', 'S', 'T', 'T', 'V', 'W', 'X', 'Y',   0,   0,   0,   0,   0,   0,
   0, 'A', 'B', 'C', 'D',   0,   0, 'G', 'H',   0,   0, 'K',   0, 'M', 'N', 'O',
   0,   0, 'R', 'S', 'T', 'T', 'V', 'W', 'X', 'Y',   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
};


static void print_pretty_phylip_dna(FILE * fp,
                                    msa_t * msa,
                                    int pad,
                                    int every,
                                    unsigned int * weights)
{
  long i,j;
  if (!msa) return;

  fprintf(fp, "%d %d P\n", msa->count, msa->length);
  for (i = 0; i < msa->count; ++i)
  {
    fprintf(fp, "%-*s", pad, msa->label[i]);
    for (j = 0; j < (long)(msa->length); ++j)
    {
      /* note: prints an extra space before the sequence */
      if (j % every == 0)
        fprintf(fp, " ");
      fprintf(fp, "%c", (char)(bpp_nt_normal[(int)(msa->sequence[i][j])]));
    }
    fprintf(fp, "\n");
  }
  assert(weights);
  fprintf(fp, "%d", weights[0]);
  for (i = 1; i < msa->length; ++i)
    fprintf(fp, " %d", weights[i]);
  fprintf(fp, "\n");
}

static void print_pretty_phylip(FILE * fp,
                                msa_t * msa,
                                int pad,
                                int every,
                                unsigned int * weights)
{
  long i,j;
  if (!msa) return;

  if (msa->dtype == BPP_DATA_DNA)
  {
    print_pretty_phylip_dna(fp,msa,pad,every,weights);
    return;
  }

  fprintf(fp, "%d %d P\n", msa->count, msa->length);
  for (i = 0; i < msa->count; ++i)
  {
    fprintf(fp, "%-*s", pad, msa->label[i]);
    for (j = 0; j < (long)(msa->length); ++j)
    {
      /* note: prints an extra space before the sequence */
      if (j % every == 0)
        fprintf(fp, " ");
      fprintf(fp, "%c", msa->sequence[i][j]);
    }
    fprintf(fp, "\n");
  }
  assert(weights);
  fprintf(fp, "%d", weights[0]);
  for (i = 1; i < msa->length; ++i)
    fprintf(fp, " %d", weights[i]);
  fprintf(fp, "\n");
}

void msa_print_phylip(FILE * fp,
                      msa_t ** msa,
                      long count,
                      unsigned int ** weights)
{
  int every = 10;  /* separate data in sequence every 10 bases */
  int pad = 4;
  int maxlen = 0;
  long i,j,k;

  /* find length of longest sequence label */
  for (i = 0; i < count; ++i)
  {
    for (j = 0; j < msa[i]->count; ++j)
    {
      k = strlen(msa[i]->label[j]);
      if (k > maxlen)
        maxlen = k;
    }
  }

  for (i = 0; i < count; ++i)
  {
    print_pretty_phylip(fp, msa[i], maxlen+pad, every, weights[i]);
    fprintf(fp,"\n");
  }
}

void msa_count_ambiguous_sites(msa_t * msa, const unsigned int * map)
{
  int i,j;
  int amb;

  msa->amb_sites_count = 0;

  if (msa->dtype == BPP_DATA_AA) return;
  assert(msa->dtype == BPP_DATA_DNA);

  for (i = 0; i < msa->length; ++i)
  {
    amb = 0;
    for (j = 0; j < msa->count; ++j)
      amb |= map[(int)(msa->sequence[j][i])];

    if (amb)
      msa->amb_sites_count++;
  }
}
static int * mark_ambiguous_sites(msa_t * msa, const unsigned int * map)
{
  int i,j;
  int amb;
  int * ambvector = (int *)xcalloc(msa->length, sizeof(int));

  msa->amb_sites_count = 0;

  for (i = 0; i < msa->length; ++i)
  {
    amb = 0;
    for (j = 0; j < msa->count; ++j)
      amb |= map[(int)(msa->sequence[j][i])];

    if (amb)
    {
      ambvector[i] = 1;
      msa->amb_sites_count++;
    }
  }

  return ambvector;
}

static int remove_ambiguous(msa_t * msa, int * ambiguous)
{
  int i,j,k;
  int amb_count = 0;

  for (i = 0; i < msa->length; ++i)
    amb_count += ambiguous[i];

  /* if all sites contain ambigous characters exit with error */
  if (amb_count == msa->length)
    return 0;

  /* we will move all ambiguous sites to the right end of the alignment */

  i = 0;  j = msa->length-1;
  for (i = 0, j = msa->length-1; i < msa->length && j >= 0; ++i, --j)
  {
    /* find next ambiguous site from left to right */
    while (i < msa->length && !ambiguous[i])
      ++i;

    /* find next non-ambigous site from right to left */
    while (j >= 0 && ambiguous[j])
      --j;

    /* if we moved all sites containing ambiguities to the right, break */
    if (j < i) break;

    /* swap sites */
    for (k = 0; k < msa->count; ++k)
      SWAP(msa->sequence[k][i],msa->sequence[k][j]);

    /* swap mark */
    SWAP(ambiguous[j],ambiguous[i]);
  }

  /* all ambiguous sites should now be at the right end. We will place
     terminating zeroes to the new sequence lengths and update the length
     variable */
  msa->length -= amb_count;

  for (k = 0; k < msa->count; ++k)
    msa->sequence[k][msa->length] = 0;

  return 1;
}

int msa_remove_ambiguous(msa_t * msa)
{
  int * ambiguous;
  int rc;

  /* get a vector indicating which sites are ambiguous */
  ambiguous = mark_ambiguous_sites(msa,pll_map_amb);

  /* remove ambiugous sites from alignment */
  rc = remove_ambiguous(msa,ambiguous);

  free(ambiguous);

  return rc;
}

int msa_remove_missing_sequences(msa_t * msa)
{
  long i,j,k;
  long deleted = 0;
  const unsigned int * map;
  char ** seq;
  char ** label;

  map = (msa->dtype == BPP_DATA_DNA) ? pll_map_nt_missing : pll_map_aa_missing;

  /* go over the sequences and delete those that comprise of missing data */
  for (i = 0; i < msa->count; ++i)
  {
    for (j = 0; j < msa->length; ++j)
      if (!map[(int)(msa->sequence[i][j])])
        break;
    
    if (j == msa->length)
    {
      ++deleted;
      free(msa->sequence[i]);
      free(msa->label[i]);

      msa->sequence[i] = NULL;
      msa->label[i] = NULL;
    }
  }
  assert(deleted <= msa->count);

  if (msa->count == deleted)
    return -1;

  /* now remove those empty records in the msa structure */
  if (deleted)
  {
    seq   = (char **)xmalloc((size_t)(msa->count - deleted) * sizeof(char *));
    label = (char **)xmalloc((size_t)(msa->count - deleted) * sizeof(char *));

    /* populate new arrays */
    for (i = 0, k = 0; i < msa->count; ++i)
    {
      assert((msa->sequence[i] && msa->label[i]) ||
             (!msa->sequence[i] && !msa->label[i]));

      if (msa->sequence[i])
      {
        seq[k]   = msa->sequence[i];
        label[k] = msa->label[i];

        ++k;
      }
    }

    /* replace */
    free(msa->sequence);
    free(msa->label);
    msa->sequence = seq;
    msa->label = label;

    msa->count -= deleted;
  }
  return deleted;
}

void msa_destroy(msa_t * msa)
{
  int i;

  if (msa->label)
  {
    for (i = 0; i < msa->count; ++i)
      if (msa->label[i])
        free(msa->label[i]);
    free(msa->label);
  }

  if (msa->sequence)
  {
    for (i = 0; i < msa->count; ++i)
      if (msa->sequence[i])
        free(msa->sequence[i]);
    free(msa->sequence);
  }

  if (msa->freqs)
    free(msa->freqs);

  free(msa);
}
