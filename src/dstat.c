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

static double * abba_tbl = NULL;
static double * baba_tbl = NULL;
static long * patt_count = NULL;

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

  #if 0
  printf("Total number of sites: %d\n", msa->length);
  for (i = 0; i < msa->length; ++i)
  {
    long s = pats[i];
    printf(" pattern: %c%c%c%c code: %5ld abba: %f babba: %f  pats: %ld\n",
           msa->sequence[0][i],
           msa->sequence[1][i],
           msa->sequence[2][i],
           msa->sequence[3][i],
           s, abba_tbl[s], baba_tbl[s], patt_count[s]);
    abba += abba_tbl[s];
    baba += baba_tbl[s];
  }
  #endif

  printf("abba: %f\n", abba);
  printf("baba: %f\n", baba);

  free(pats);
  /* 
  decode_site(s);

  printf("Checking from precomputed table:\n");
  index = s[0] | (s[1] << 4) | (s[2] << 8) | (s[3] << 12);
  printf("abba: %f\n", abba_tbl[index]);
  printf("baba: %f\n", baba_tbl[index]);
  */
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

int debug_decode_site(unsigned int * s)
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

void cmd_dstat()
{
  long i,j,k,m;
  long msa_count;
  phylip_t * fd;
  msa_t ** msa_list;

  printf("Pre-computing table for site scores...\n");
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

  /* TODO: For now we only allow one alignment */
  assert(msa_count == 1);

  char ** taxa = split4(opt_dstat);

  printf("Tree: (((%s,%s),%s),%s);\n", taxa[0], taxa[1], taxa[2], taxa[3]);
  printf("Testing introgression between %s and %s, and between %s and %s\n",
         taxa[0], taxa[2], taxa[1], taxa[2]);


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
  for (i = 0; i < 4; ++i) free(labels[i]);
  free(labels);

  /* concatenate (possibly) multiple alignments and fill in missing data */
  msa_t * concat = phylip_concat(msa_list, msa_count);

  #if 0
  phylip_print(stdout, concat);
  #endif

  calculate_d(concat);

  for (i = 0; i < concat->count; ++i)
  {
    free(concat->label[i]);
    free(concat->sequence[i]);
  }
  free(concat->label);
  free(concat->sequence);
  free(concat);

  for (i = 0; i < 4; ++i)
    free(taxa[i]);
  free(taxa);

  for (i = 0; i < msa_count; ++i)
    msa_destroy(msa_list[i]);
  free(msa_list);

  if (abba_tbl)
    free(abba_tbl);
  if (baba_tbl)
    free(baba_tbl);
  if (patt_count)
    free(patt_count);

}
