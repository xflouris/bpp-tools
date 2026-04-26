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
  long tag_count = 0;
  long label_count = 0;
  long species_count = 0;
  long index_size = 0;
  long newmsa_count = 0;
  long * tag_matches = NULL;
  long * label_matches = NULL;
  long * species_matches = NULL;
  long * index = NULL;
  list_t * maplist = NULL;
  phylip_t * fp_in;
  FILE * fp_out;
  char ** label_list = NULL;
  char ** species_list = NULL;
  char ** tag_list = NULL;
  msa_t ** msa_list;
  msa_t ** new_list;

  /* split labels and store in array */
  if (opt_label_list)
  {
    label_list = split(opt_label_list, ",", &label_count);
    if (!label_list)
      fatal("Cannot parse list of labels specified with --label_list");
    label_matches = (long *)xcalloc((size_t)label_count, sizeof(long));
  }

  /* split tags and store in array */
  if (opt_tag_list)
  {
    tag_list = split(opt_tag_list, ",", &tag_count);
    if (!tag_list)
      fatal("Cannot parse list of tags specified with --tag_list");
    tag_matches = (long *)xcalloc((size_t)tag_count, sizeof(long));
  }

  /* split species and store in array */
  if (opt_species_list)
  {
    species_list = split(opt_species_list, ",", &species_count);
    if (!species_list)
      fatal("Cannot parse list of species specified with --species_list");
    species_matches = (long *)xcalloc((size_t)species_count, sizeof(long));
  }

  if (!opt_label_list && !opt_tag_list && !opt_species_list)
    fatal("Please specify labels/tags/species to extract using:\n"
          "  --label_list\n"
          "  --tag_list\n"
          "  --tag_species");

  /* read map file if specified */
  if (opt_mapfile)
  {
    maplist = parse_mapfile(opt_mapfile);
    maplist_print(maplist);
  }

  if (opt_species_list && !opt_mapfile)
    fatal("Please specify a map file with --map when using --species_list");

  #if 0
  /* count number of specimens and sequences in list */
  sp_count = label_count = 0;
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
  #endif

  /* TODO: check for duplicates */

  /* print */
  if (opt_verbose)
  {
    if (opt_label_list)
    {
      printf("Labels:\n");
      for (i = 0; i < label_count; ++i)
        printf("%ld : %s\n", i, label_list[i]);
    }
    if (opt_tag_list)
    {
      printf("Tags:\n");
      for (i = 0; i < tag_count; ++i)
        printf("%ld : %s\n", i, tag_list[i]);
    }
    if (opt_species_list)
    {
      printf("Species:\n");
      maplist_print(maplist);
    }
  }

  /* check for msa file */
  if (!opt_msafile)
    fatal("Please specify an input PHYLIP file using --msa");

  /* open phylip file */
  fp_in = phylip_open(opt_msafile, pll_map_fasta);
  if (!fp_in)
    fatal("Cannot open file %s", opt_msafile);

  /* read alignment */
  msa_list = phylip_parse_multisequential(fp_in, &msa_count);
  assert(msa_list);
  phylip_close(fp_in);

  /* writer does not emit the `P MODEL` header / weights line, so a
     compressed input would produce a malformed output. Reject explicitly. */
  for (i = 0; i < msa_count; ++i)
    if (msa_list[i]->pattern_weights)
      fatal("--extract does not yet support pattern-compressed alignments");

  /* filter out sequences */

  new_list = (msa_t **)xmalloc((size_t)msa_count*sizeof(msa_t *));

  /* go through the alignments and check for label/tag/species match */
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

    /* zero-out index */
    memset(index,0,msa_list[i]->count*sizeof(long));
    newmsa_count = 0;

    for (j = 0; j < msa_list[i]->count; ++j)
    {
      for (k = 0; k < tag_count && !index[j]; ++k)
        if (ends_with_hat_str(msa_list[i]->label[j], tag_list[k]))
          break;
      if (k != tag_count)
      {
        index[j] = 1; 
        newmsa_count++;
        tag_matches[k]++;
        continue;
      }

      for (k = 0; k < label_count; ++k)
        if (starts_with(msa_list[i]->label[j], label_list[k]))
          break;
      if (k != label_count)
      {
        index[j] = 1;
        newmsa_count++;
        label_matches[k]++;
        continue;
      }
    }

    if (maplist)
    {
      for (j = 0; j < msa_list[i]->count; ++j)
      {
        list_item_t * li = maplist->head;
        while (li)
        {
          mapping_t * map = (mapping_t *)(li->data);
          if (ends_with_hat_str(msa_list[i]->label[j], map->individual))
            break;
          li = li->next;
        }
        if (li)
        {
          mapping_t * map = (mapping_t *)(li->data);
          for (k = 0; k < species_count; ++k)
          {
            if (!strcmp(map->species, species_list[k]))
            {
              species_matches[k]++;
              break;
            }
          }
          if (k && !index[j])
          {
            index[j] = 1;
            newmsa_count++;
          }
        }
      }
    }

    if (!newmsa_count) continue;

    /* allocate the alignment structure */
    msa_t * msa   = (msa_t *)xcalloc(1, sizeof(msa_t));
    msa->length   = msa_list[i]->length;
    msa->count    = newmsa_count;
    msa->label    = (char **)xmalloc((size_t)newmsa_count * sizeof(char *));
    msa->sequence = (char **)xmalloc((size_t)newmsa_count * sizeof(char *));

    /* copy the data */
    k = 0;
    for (j = 0; j < msa_list[i]->count; ++j)
    {
      if (index[j])
      {
        msa->label[k] = xstrdup(msa_list[i]->label[j]);
        msa->sequence[k++] = xstrdup(msa_list[i]->sequence[j]);
      }
    }

    /* add new alignment into a list */
    new_list[m++] = msa;
  }

  /* write alignments to file */
  fp_out = opt_outfile ? xopen(opt_outfile,"w") : stdout;
  for (i = 0; i < m; ++i)
    phylip_print(fp_out, new_list[i]);

  if (opt_outfile)
    xclose(fp_out);

  /* check for labels that were not found */
  if (opt_label_list)
  {
    long lcount = 0;
    for (i = 0; i < label_count; ++i)
    {
      if (label_matches[i])
        lcount++;
    }
    if (lcount != label_count)
    {
      xwarn("The following labels were not found in %s:", opt_msafile);
      for (i = 0; i < label_count; ++i)
        if (!label_matches[i])
          xwarn("  %s", label_list[i]);
    }
  }
  /* check for tags that were not found */
  if (opt_tag_list)
  {
    long tcount = 0;
    for (i = 0; i < tag_count; ++i)
    {
      if (tag_matches[i])
        tcount++;
    }
    if (tcount != tag_count)
    {
      xwarn("The following tags were not found in %s:", opt_msafile);
      for (i = 0; i < tag_count; ++i)
        if (!tag_matches[i])
          xwarn("  %s", tag_list[i]);
    }
  }
  /* check for species that were not found */
  if (opt_species_list)
  {
    long scount = 0;
    for (i = 0; i < species_count; ++i)
    {
      if (species_matches[i])
        scount++;
    }
    if (scount != species_count)
    {
      xwarn("The following species were not found in %s:", opt_msafile);
      for (i = 0; i < species_count; ++i)
        if (!species_matches[i])
          xwarn("  %s", species_list[i]);
    }
  }

  if (!opt_quiet)
  {
    long newmsa_total_seqcount = 0;
    for (i = 0; i < m; ++i)
      newmsa_total_seqcount += new_list[i]->count;
      
    printf("%ld sequences matched and copied to %ld new alignments\n", 
           newmsa_total_seqcount, m);

  }

  /* dealloc */
  if (label_list)
  {
    for (i = 0; i < label_count; ++i)
      free(label_list[i]);
    free(label_list);
  }

  if (tag_list)
  {
    for (i = 0; i < tag_count; ++i)
      free(tag_list[i]);
    free(tag_list);
  }

  if (species_list)
  {
    for (i = 0; i < species_count; ++i)
      free(species_list[i]);
    free(species_list);
    //map_dealloc(maplist); 
    list_clear(maplist, map_dealloc);
    free(maplist);
  }


  if (index) free(index);
  if (label_matches) free(label_matches);
  if (tag_matches) free(tag_matches);
  if (species_matches) free(species_matches);

  for (i = 0; i < msa_count; ++i)
    msa_destroy(msa_list[i]);
  free(msa_list);
  for (i = 0; i < m; ++i)
    msa_destroy(new_list[i]);
  free(new_list);
}
