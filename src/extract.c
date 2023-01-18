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

static char ** split(const char * s, const char * d, long * token_count)
{
  long i,k;
  long del_count = 0;
  char ** tokens = NULL;

  assert(strlen(d) == 1);

  /* compute number of commas in list of tips */
  for (i = 0; i < (long)strlen(s); ++i)
    if (s[i] == *d)
      ++del_count;
  
  tokens = (char **)xmalloc((size_t)(del_count+1) * sizeof(char *));

  k = 0;
  while (*s)
  {
    /* get next taxon */
    size_t token_len = strcspn(s,d);
    if (!token_len)
    {
      free(tokens);
      return NULL;
    }

    tokens[k++] = xstrndup(s, token_len);

    s += token_len;
    assert(*s == d[0] || *s == '\0');
    if (*s == d[0])
      ++s;
  }

  *token_count = del_count+1;

  if (*token_count == 1 && !tokens[0])
  {
    free(tokens);
    tokens = NULL;
  }

  return tokens;
}

static long ends_with(const char * x, const char * suffix)
{
  size_t xlen, slen;
  const char * p;

  xlen = strlen(x);
  slen = strlen(suffix);

  if (slen > xlen)
    return 0;

  p = x + xlen - slen;
  if (!strcmp(p,suffix))
    return 1;

  return 0;
}

static long starts_with(const char * x, const char * prefix)
{
  size_t xlen, plen;

  xlen = strlen(x);
  plen = strlen(prefix);

  if (plen > xlen)
    return 0;

  if (!strncmp(x,prefix,plen))
    return 1;

  return 0;
}

void cmd_extract()
{
  long i,j,k,m;
  long msa_count;
  long sp_count;
  long seq_count;
  long token_count;
  long index_size = 0;
  long * index = NULL;
  long seq_copy_count = 0;
  phylip_t * fp_in;
  char ** sp_tokens = NULL;
  char ** seq_tokens = NULL;
  msa_t ** msa_list;
  msa_t ** new_list;

  char ** tokens = split(opt_extract, ",", &token_count);
  if (!tokens)
    fatal("Cannot parse tokens");

  /* count number of specimens and sequences in list */
  sp_count = seq_count = 0;
  for (i = 0; i < token_count; ++i)
  {
    assert(tokens[i]);
    if (tokens[i][0] == '^')
      ++sp_count;
    else
      ++seq_count;
  }

  /* allocate arrays for storing speciments and sequences */
  if (sp_count)
    sp_tokens = (char **)xmalloc((size_t)sp_count * sizeof(char *));
  if (seq_count)
    seq_tokens = (char **)xmalloc((size_t)seq_count * sizeof(char *));

  /* separate specimen and sequences */
  sp_count = seq_count = 0;
  for (i = 0; i < token_count; ++i)
    if (tokens[i][0] == '^')
      sp_tokens[sp_count++] = tokens[i];
    else
      seq_tokens[seq_count++] = tokens[i];

  /* TODO: check for duplicates */

  #if 0
  /* print */
  printf("Specimens:\n");
  for (i = 0; i < sp_count; ++i)
    printf("%ld : %s\n", i, sp_tokens[i]);
  printf("Sequences:\n");
  for (i = 0; i < seq_count; ++i)
    printf("%ld : %s\n", i, seq_tokens[i]);
  #endif

  /* open phylip file */
  fp_in = phylip_open(opt_msafile, pll_map_fasta);
  if (!fp_in)
    fatal("Cannot open file %s", opt_msafile);

  /* read alignment */
  msa_list = phylip_parse_multisequential(fp_in, &msa_count);
  assert(msa_list);
  phylip_close(fp_in);

  /* filter out sequences */

  new_list = (msa_t **)xmalloc((size_t)msa_count*sizeof(msa_t *));

  for (m=0,i=0; i < msa_count; ++i)
  {
    /* create index array for marking sequences that will be copied */
    if (msa_list[i]->count > index_size)
    {
      if (index)
        free(index);
      index = (long *)xmalloc((size_t)msa_list[i]->count*sizeof(long));
      index_size = msa_list[i]->count;
    }
    memset(index,0,msa_list[i]->count*sizeof(long));
    seq_copy_count = 0;

    for (j = 0; j < msa_list[i]->count; ++j)
    {
      for (k = 0; k < sp_count && !index[j]; ++k)
        if (ends_with(msa_list[i]->label[j], sp_tokens[k]))
          break;
      if (k != sp_count)
        index[j] = 1; 

      for (k = 0; k < seq_count; ++k)
        if (starts_with(msa_list[i]->label[j], seq_tokens[k]))
          break;
      if (k != seq_count)
        index[j] = 1;
      
      if (index[j])
        seq_copy_count++;
    }

    if (!seq_copy_count) continue;

    msa_t * msa   = (msa_t *)xcalloc(1, sizeof(msa_t));
    msa->length   = msa_list[i]->length;
    msa->count    = seq_copy_count;
    msa->label    = (char **)xmalloc((size_t)seq_copy_count * sizeof(char *));
    msa->sequence =  (char **)xmalloc((size_t)seq_copy_count * sizeof(char *));

    /* copy */
    k = 0;
    for (j = 0; j < msa_list[i]->count; ++j)
    {
      if (index[j])
      {
        msa->label[k] = xstrdup(msa_list[i]->label[j]);
        msa->sequence[k++] = xstrdup(msa_list[i]->sequence[j]);
      }
    }
    new_list[m++] = msa;
  }

  FILE * fpout = opt_outfile ? xopen(opt_outfile,"w") : stdout;
  for (i = 0; i < m; ++i)
    phylip_print(fpout, new_list[i]);

  if (opt_outfile)
    fclose(fpout);

  /* dealloc */
  free(sp_tokens);
  free(seq_tokens);
  for (i = 0; i < token_count; ++i)
    free(tokens[i]);
  free(tokens);

  if (index) free(index);

  for (i = 0; i < msa_count; ++i)
    msa_destroy(msa_list[i]);
  free(msa_list);
  for (i = 0; i < m; ++i)
    msa_destroy(new_list[i]);
  free(new_list);
}
