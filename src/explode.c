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

void cmd_explode()
{
  long i;
  long msa_count;
  char * filename;
  char * outfile;
  phylip_t * fp_in;
  FILE * fp_out;
  msa_t ** msa_list;

  /* open phylip file */
  fp_in = phylip_open(opt_msafile, pll_map_fasta);
  if (!fp_in)
    fatal("Cannot open file %s", opt_msafile);

  /* read alignment */
  msa_list = phylip_parse_multisequential(fp_in, &msa_count);
  assert(msa_list);
  phylip_close(fp_in);
  
  /* write separate files */
  outfile = opt_outfile ? xstrdup(opt_outfile) : xstrdup(opt_msafile);
  for (i = 0; i < msa_count; ++i)
  {
    xasprintf(&filename, "%s.%ld", outfile, i);
    fp_out = xopen(filename, "w");
    phylip_print(fp_out, msa_list[i]);

    msa_destroy(msa_list[i]);
    free(filename);
    fclose(fp_out);
  }

  free(outfile);
  free(msa_list);
}
