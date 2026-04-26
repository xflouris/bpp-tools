/*
    Copyright (C) 2022-2026 Tomas Flouri

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

    Contact: Tomas Flouri <t.flouri@ucl.ac.uk>,
    Department of Genetics, Evolution and Environment,
    University College London, Gower Street, London WC1E 6BT, England
*/

#include "bpp-tools.h"

static int strptr_cmp(const void * a, const void * b)
{
  return strcmp(*(const char **)a, *(const char **)b);
}

/* Collect every label across every locus, sort lexicographically, and
   count distinct entries. If out_unique is non-NULL, it is set to a
   freshly-malloc'd array of pointers (into the original label storage)
   listing the unique labels in sorted order; the caller must free the
   array but not the strings. */
static long count_unique_labels(msa_t ** msa_list, long count,
                                char *** out_unique)
{
  long total = 0, i, j;
  for (i = 0; i < count; ++i)
    total += msa_list[i]->count;

  if (!total)
  {
    if (out_unique) *out_unique = NULL;
    return 0;
  }

  char ** all = (char **)xmalloc((size_t)total * sizeof(char *));
  long k = 0;
  for (i = 0; i < count; ++i)
    for (j = 0; j < msa_list[i]->count; ++j)
      all[k++] = msa_list[i]->label[j];

  qsort(all, (size_t)total, sizeof(char *), strptr_cmp);

  long unique = 1;
  for (i = 1; i < total; ++i)
    if (strcmp(all[i], all[i-1]) != 0)
      ++unique;

  if (out_unique)
  {
    char ** uniq = (char **)xmalloc((size_t)unique * sizeof(char *));
    long u = 0;
    uniq[u++] = all[0];
    for (i = 1; i < total; ++i)
      if (strcmp(all[i], all[i-1]) != 0)
        uniq[u++] = all[i];
    *out_unique = uniq;
  }

  free(all);
  return unique;
}

static long sum_pattern_weights(const msa_t * m)
{
  long s = 0, i;
  if (!m->pattern_weights) return 0;
  for (i = 0; i < m->length; ++i)
    s += m->pattern_weights[i];
  return s;
}

static const char * model_str(int compress_model)
{
  if (compress_model == COMPRESS_JC69)    return "JC69";
  if (compress_model == COMPRESS_GENERAL) return "GTR";
  return "-";
}

void cmd_info()
{
  long msa_count, i;
  phylip_t * fp_in;
  msa_t ** msa_list;

  if (!opt_msafile)
    fatal("Specify input alignment with --msa FILENAME");

  fp_in = phylip_open(opt_msafile, pll_map_fasta);
  if (!fp_in)
    fatal("Cannot open file %s", opt_msafile);

  msa_list = phylip_parse_multisequential(fp_in, &msa_count);
  assert(msa_list);
  phylip_close(fp_in);

  /* compression status across all loci */
  long n_compressed = 0;
  int  shared_model = -1;
  int  mixed_models = 0;
  for (i = 0; i < msa_count; ++i)
  {
    if (msa_list[i]->pattern_weights)
    {
      ++n_compressed;
      if (shared_model == -1)
        shared_model = msa_list[i]->compress_model;
      else if (shared_model != msa_list[i]->compress_model)
        mixed_models = 1;
    }
  }
  int any_compressed = (n_compressed > 0);
  int all_compressed = (n_compressed == msa_count);

  /* file / format header */
  printf("File:    %s\n", opt_msafile);
  if (!any_compressed)
    printf("Format:  uncompressed\n");
  else if (all_compressed && !mixed_models)
    printf("Format:  pattern-compressed (P %s)\n", model_str(shared_model));
  else if (all_compressed && mixed_models)
    printf("Format:  pattern-compressed (mixed JC69/GTR)\n");
  else
    printf("Format:  mixed (%ld compressed, %ld uncompressed)\n",
           n_compressed, msa_count - n_compressed);
  printf("\n");

  /* per-locus table */
  printf("Per-locus summary:\n");
  if (any_compressed)
    printf("  %5s  %5s  %8s  %5s  %12s\n",
           "idx", "seqs", "length", "model", "weights_sum");
  else
    printf("  %5s  %5s  %8s\n", "idx", "seqs", "length");

  int truncate = (msa_count > 10) && !opt_per_locus;

  for (i = 0; i < msa_count; ++i)
  {
    if (truncate && i == 5)
    {
      printf("  ...    (%ld loci omitted; pass --per-locus to show all)\n",
             msa_count - 10);
      i = msa_count - 6;  /* loop ++i lands on msa_count - 5 next iter */
      continue;
    }
    msa_t * m = msa_list[i];
    if (any_compressed)
    {
      if (m->pattern_weights)
        printf("  %5ld  %5d  %8d  %5s  %12ld\n",
               i + 1, m->count, m->length,
               model_str(m->compress_model),
               sum_pattern_weights(m));
      else
        printf("  %5ld  %5d  %8d  %5s  %12s\n",
               i + 1, m->count, m->length, "-", "-");
    }
    else
    {
      printf("  %5ld  %5d  %8d\n", i + 1, m->count, m->length);
    }
  }
  printf("\n");

  /* aggregate statistics */
  long total_seqs = 0;
  long total_len  = 0;
  long total_orig = 0;
  long min_len = (msa_count > 0) ? msa_list[0]->length : 0;
  long max_len = 0;
  for (i = 0; i < msa_count; ++i)
  {
    msa_t * m = msa_list[i];
    total_seqs += m->count;
    total_len  += m->length;
    if (m->length < min_len) min_len = m->length;
    if (m->length > max_len) max_len = m->length;
    total_orig += m->pattern_weights ? sum_pattern_weights(m) : m->length;
  }

  char ** unique_labels = NULL;
  long n_unique = count_unique_labels(msa_list, msa_count, &unique_labels);
  double mean_len = msa_count ? (double)total_len / (double)msa_count : 0.0;

  /* Each aggregate row uses a single space between the colon and the
     value, keeping a literal grep for e.g. "Number of alignments: 5"
     simple. Visual column alignment is intentionally not pursued here. */
  printf("Aggregate:\n");
  printf("  Number of alignments: %ld\n", msa_count);
  printf("  Total sequences: %ld\n", total_seqs);
  printf("  Unique sequence labels: %ld\n", n_unique);
  printf("  Per-locus length: mean=%.1f, min=%ld, max=%ld\n",
         mean_len, min_len, max_len);
  if (any_compressed)
  {
    printf("  Total patterns: %ld\n", total_len);
    printf("  Total original sites: %ld\n", total_orig);
    if (total_orig > 0)
      printf("  Compression ratio: %.3f (%ld / %ld)\n",
             (double)total_len / (double)total_orig,
             total_len, total_orig);
  }
  else
  {
    printf("  Total sites: %ld\n", total_len);
  }

  if (opt_show_labels && n_unique > 0)
  {
    printf("\nUnique sequence labels (%ld):\n", n_unique);
    for (i = 0; i < n_unique; ++i)
      printf("  %s\n", unique_labels[i]);
  }

  free(unique_labels);
  for (i = 0; i < msa_count; ++i)
    msa_destroy(msa_list[i]);
  free(msa_list);
}
