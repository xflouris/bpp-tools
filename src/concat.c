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

static char ** getuniqueseqlabels(msa_t ** msa_list,
                                  long msa_count,
                                  long max_seq_count,
                                  long * missing_data,
                                  long * count)
{
  long i,j,k,m;
  long all_labels_count = 0;
  long all_labels_size = max_seq_count;
  char ** all_labels = (char **)xmalloc((size_t)max_seq_count * sizeof(char *));
  char ** new_labels = (char **)xmalloc((size_t)max_seq_count * sizeof(char *));

  *missing_data = 0;

  for (i = 0; i < msa_count; ++i)
  {
    for (m=0,j=0; j < msa_list[i]->count; ++j)
    {
      for (k = 0; k < all_labels_count; ++k)
        if (!strcmp(msa_list[i]->label[j],all_labels[k]))
          break;

      /* if labeled already exists, move to next */
      if (k != all_labels_count)
        continue;

      /* new label */
      new_labels[m++] = xstrdup(msa_list[i]->label[j]);
    }
    
    if ((m && i) || msa_list[i]->count < max_seq_count)
    {
      /* missing data (sequences) found */
      *missing_data = 1;
    }

    /* new labels found */
    if (m)
    {
      /* check if there is enough space in all_labels */
      if (all_labels_size - all_labels_count < m)
      {
        /* reallocate */
        long allocsize = all_labels_count + m - (all_labels_size - all_labels_count);
        char ** tmp = (char **)xmalloc((size_t)allocsize*sizeof(char *));
        for (j = 0; j < all_labels_count; ++j)
          tmp[j] = all_labels[j];
        free(all_labels);
        all_labels = tmp;
        all_labels_size = allocsize;
      }

      for (j = 0; j < m; ++j)
        all_labels[all_labels_count++] = new_labels[j];
    }
  }
  free(new_labels);

  *count = all_labels_count;
  return all_labels;
}

msa_t * concatenate(msa_t ** msa_list, long msa_count)
{
  long i,j,k;
  long offset = 0;
  long max_seq_count = 0;
  long max_seq_size = 0;
  long missing_data = 0;
  long concat_seq_size = 0;
  long concat_seq_count = 0;
  msa_t * concat;
  char ** labels;

  /* Pattern-compression consistency check across loci. We allow:
      - all loci uncompressed (the classic path)
      - all loci compressed with the same compress_model (multi-locus
        compressed concatenation; patterns from each locus keep their own
        weights and are laid end-to-end in the concatenated MSA)
     We reject:
      - a mix of compressed and uncompressed (the semantics of the missing-
        data fill would not be well defined on the compressed side)
      - a mix of JC69 and GTR (the two re-encodings are incompatible) */
  int any_compressed = 0, all_compressed = 1;
  int shared_model = -1;
  for (i = 0; i < msa_count; ++i)
  {
    if (msa_list[i]->pattern_weights)
    {
      any_compressed = 1;
      if (shared_model == -1)
        shared_model = msa_list[i]->compress_model;
      else if (shared_model != msa_list[i]->compress_model)
        fatal("Cannot concatenate alignments with mixed compression models "
              "(JC69 + GTR)");
    }
    else
      all_compressed = 0;
  }
  if (any_compressed && !all_compressed)
    fatal("Cannot concatenate a mix of compressed and uncompressed "
          "alignments");

  /* single pattern-compressed locus: no concatenation needed; return a fresh
     deep copy that preserves sequences, labels, pattern_weights and
     compress_model so the caller can safely free msa_list independently */
  if (msa_count == 1 && msa_list[0]->pattern_weights)
  {
    msa_t * src = msa_list[0];
    concat = (msa_t *)xcalloc(1, sizeof(msa_t));
    concat->count = src->count;
    concat->length = src->length;
    concat->compress_model = src->compress_model;
    concat->sequence = (char **)xmalloc((size_t)src->count * sizeof(char *));
    concat->label = (char **)xmalloc((size_t)src->count * sizeof(char *));
    for (i = 0; i < src->count; ++i)
    {
      concat->sequence[i] = (char *)xmalloc((size_t)(src->length + 1));
      memcpy(concat->sequence[i], src->sequence[i], (size_t)src->length);
      concat->sequence[i][src->length] = 0;
      concat->label[i] = xstrdup(src->label[i]);
    }
    concat->pattern_weights =
      (unsigned int *)xmalloc((size_t)src->length * sizeof(unsigned int));
    memcpy(concat->pattern_weights, src->pattern_weights,
           (size_t)src->length * sizeof(unsigned int));
    return concat;
  }

  /* get length of concatenated alignment, max sequence count and length */
  for (i = 0; i < msa_count; ++i)
  {
    if (msa_list[i]->count > max_seq_count)
      max_seq_count = msa_list[i]->count;
    if (msa_list[i]->length > max_seq_size)
      max_seq_size = msa_list[i]->length;

    concat_seq_size += msa_list[i]->length;
  }

  /* create the list of sequence names */
  labels = getuniqueseqlabels(msa_list,
                              msa_count,
                              max_seq_count,
                              &missing_data,
                              &concat_seq_count);
  if (missing_data)
  {
    xwarn("Some alignments contain missing sequences");
    if (!opt_nachar)
    {
      xwarn("You can specify the missing data character using --nachar");
      xwarn("Using default symbol %c for missing data", opt_nachar_default);
      opt_nachar = (char *)xmalloc(sizeof(char)+1);
      opt_nachar[0] = opt_nachar_default;
      opt_nachar[1] = '\0';
    }
    else
    {
      xwarn("Using user-specified symbol %c for missing data", opt_nachar[0]);
    }
  }

  /* create concatenated alignment structure */
  concat = (msa_t *)xcalloc(1, sizeof(msa_t));
  concat->count = concat_seq_count;
  concat->length = concat_seq_size;
  concat->label = labels;
  concat->sequence = (char **)xmalloc((size_t)concat_seq_count*sizeof(char *));
  for (i = 0; i < concat_seq_count; ++i)
    concat->sequence[i] = (char *)xmalloc((size_t)(concat_seq_size+1) *
                                          sizeof(char));

  /* construct the concatenated alignment */
  for (i = 0; i < msa_count; ++i)
  {
    for (j = 0; j < concat_seq_count; ++j)
    {
      for (k = 0; k < msa_list[i]->count; ++k)
      {
        if (!strcmp(msa_list[i]->label[k], concat->label[j]))
          break;
      }
      
      if (k == msa_list[i]->count)
      {
        /* no sequence j found in current alignment, fill with missing data */
        memset(concat->sequence[j]+offset, opt_nachar[0], msa_list[i]->length);
      }
      else
      {
        /* sequence j found */
        memcpy(concat->sequence[j]+offset,
               msa_list[i]->sequence[k],
               msa_list[i]->length * sizeof(char));
      }
    }
    offset += msa_list[i]->length;
  }

  assert(offset == concat_seq_size);
  for (j = 0; j < concat_seq_count; ++j)
    concat->sequence[j][offset] = 0;

  /* multi-locus compressed: concatenate the per-locus weight vectors in the
     same order as the sequences were laid out */
  if (shared_model != -1)
  {
    long woff = 0;
    concat->compress_model = shared_model;
    concat->pattern_weights =
      (unsigned int *)xmalloc((size_t)concat_seq_size * sizeof(unsigned int));
    for (i = 0; i < msa_count; ++i)
    {
      memcpy(concat->pattern_weights + woff,
             msa_list[i]->pattern_weights,
             (size_t)msa_list[i]->length * sizeof(unsigned int));
      woff += msa_list[i]->length;
    }
    assert(woff == concat_seq_size);
  }

  return concat;
}

void cmd_concat()
{
  long i;
  long msa_count;
  char * outfile;
  char * partfile;
  phylip_t * fp_in;
  FILE * fp_out;
  msa_t * concat;
  msa_t ** msa_list;

  /* check that there is only one missing data character */
  if (opt_nachar)
  {
    if (strlen(opt_nachar) != 1)
      fatal("Option --nachar requires exactly one character");
  }

  /* open phylip file */
  fp_in = phylip_open(opt_msafile, pll_map_fasta);
  if (!fp_in)
    fatal("Cannot open file %s", opt_msafile);

  /* read alignment */
  msa_list = phylip_parse_multisequential(fp_in, &msa_count);
  assert(msa_list);
  phylip_close(fp_in);

  /* writer does not emit the `P MODEL` header / weights line, so a
     compressed input would produce a malformed output. */
  for (i = 0; i < msa_count; ++i)
    if (msa_list[i]->pattern_weights)
      fatal("--concat does not yet support pattern-compressed alignments");

  /* open output file */
  if (opt_outfile)
    outfile = xstrdup(opt_outfile);
  else
    xasprintf(&outfile, "%s.concat.txt", opt_msafile);

  fp_out = xopen(outfile,"w");
  
  /* concatenate */
  concat = concatenate(msa_list, msa_count);

  #if 1
  phylip_print(fp_out, concat);
  #endif

  /* create partition file */
  xasprintf(&partfile, "%s.part", opt_outfile);
  FILE * fp_part = xopen(partfile, "w");

  long offset;
  for (i = 0, offset = 1; i < msa_count; ++i)
  {
    fprintf(fp_part,
            "%ld-%ld, DNA, part%ld\n",
            offset,
            offset+msa_list[i]->length-1,
            i+1);
    offset += msa_list[i]->length;
  }
  assert(offset == concat->length+1);

  /* print concatenated alignment info */
  if (!opt_quiet)
  {
    printf("%ld alignments found in %s\n", msa_count, opt_msafile);
    printf("Number of sites in concatenated alignment: %d\n", concat->length);
    printf("Number of sequences in concatenated alignment: %d\n", concat->count);
    printf("Concatenated alignment stored in %s\n", outfile);
    printf("Partition information stored in %s\n", partfile);
  }
  
  /* dealloc */
  msa_destroy(concat);
  for (i = 0; i < msa_count; ++i)
    msa_destroy(msa_list[i]);
  free(msa_list);
  free(outfile);
  free(partfile);
  xclose(fp_out);
  xclose(fp_part);
}
