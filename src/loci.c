/*
    Copyright (C) 2026 Tomas Flouri

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

/* Parse a comma-separated list of 1-indexed locus numbers (with optional
   hyphenated ranges) into a 0/1 mask of size n_loci.  Out-of-range or
   malformed input is fatal.  Caller owns the returned buffer.

   Spec format:  "1,2,5,8,20-100,7"
   - 1-indexed
   - whitespace allowed around tokens
   - duplicates and overlapping ranges silently merged (set semantics)
   - reversed ranges (5-3) are accepted and treated as (3-5)
*/
static unsigned char * parse_locus_csv_to_mask(const char * csv,
                                               long n_loci,
                                               const char * flag_name)
{
  unsigned char * mask = (unsigned char *)xcalloc((size_t)n_loci, 1);
  const char * p = csv;

  while (*p)
  {
    while (*p == ' ' || *p == '\t') ++p;
    if (!*p) break;

    char * end;
    long lo = strtol(p, &end, 10);
    if (end == p || lo <= 0)
      fatal("%s: malformed index near '%s'", flag_name, p);
    p = end;

    long hi = lo;
    while (*p == ' ' || *p == '\t') ++p;
    if (*p == '-')
    {
      ++p;
      while (*p == ' ' || *p == '\t') ++p;
      hi = strtol(p, &end, 10);
      if (end == p || hi <= 0)
        fatal("%s: malformed range in spec '%s'", flag_name, csv);
      p = end;
      if (hi < lo) { long t = lo; lo = hi; hi = t; }
    }

    if (lo > n_loci || hi > n_loci)
      fatal("%s: locus index %ld out of range (1..%ld)",
            flag_name, (hi > n_loci) ? hi : lo, n_loci);

    for (long i = lo; i <= hi; ++i)
      mask[i-1] = 1;

    while (*p == ' ' || *p == '\t') ++p;
    if (*p == ',') { ++p; continue; }
    if (*p)
      fatal("%s: unexpected character '%c' in spec", flag_name, *p);
  }

  return mask;
}

/* Write the selected loci to opt_outfile (or stdout when NULL), using the
   compressed writer for compressed MSAs and the plain writer otherwise.
   Multi-sequential blocks are separated by a single blank line, matching
   what phylip_parse_multisequential expects.

   keep == 1: emit MSAs where mask[i] is set.
   keep == 0: emit the complement. */
static void write_filtered_msa_list(msa_t ** msa_list,
                                    long count,
                                    const unsigned char * mask,
                                    int keep)
{
  long i;
  FILE * fp = stdout;

  if (opt_outfile)
    fp = xopen(opt_outfile, "w");

  int first = 1;
  for (i = 0; i < count; ++i)
  {
    int included = mask[i] ? 1 : 0;
    if (keep ? !included : included)
      continue;

    if (!first)
      fprintf(fp, "\n");
    first = 0;

    if (msa_list[i]->pattern_weights)
      phylip_print_compressed(fp, msa_list[i]);
    else
      phylip_print(fp, msa_list[i]);
  }

  if (opt_outfile)
    fclose(fp);
}

static void run_loci_filter(const char * csv,
                            int keep,
                            const char * flag_name)
{
  long i;

  if (!opt_msafile)
    fatal("Specify input alignment with --msa FILENAME");

  phylip_t * fd = phylip_open(opt_msafile, pll_map_fasta);
  if (!fd)
    fatal("Cannot open file %s", opt_msafile);

  long msa_count;
  msa_t ** msa_list = phylip_parse_multisequential(fd, &msa_count);
  assert(msa_list);
  phylip_close(fd);

  unsigned char * mask = parse_locus_csv_to_mask(csv, msa_count, flag_name);
  write_filtered_msa_list(msa_list, msa_count, mask, keep);

  free(mask);
  for (i = 0; i < msa_count; ++i)
    msa_destroy(msa_list[i]);
  free(msa_list);
}

void cmd_keep_loci(void)
{
  run_loci_filter(opt_keep_loci, 1, "--keep-loci");
}

void cmd_drop_loci(void)
{
  run_loci_filter(opt_drop_loci, 0, "--drop-loci");
}
