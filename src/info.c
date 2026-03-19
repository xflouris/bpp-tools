#include "bpp-tools.h"

void cmd_info()
{
  long msa_count;
  phylip_t * fp_in;
  msa_t ** msa_list;

  fp_in = phylip_open(opt_msafile, pll_map_fasta);
  if (!fp_in)
    fatal("Cannot open file %s", opt_msafile);

  /* read alignment */
  msa_list = phylip_parse_multisequential(fp_in, &msa_count);
  assert(msa_list);
  phylip_close(fp_in);

  printf("Number of alignments: %ld\n", msa_count);
}
