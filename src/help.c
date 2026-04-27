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

    Contact: Tomas Flouri <t.flouris@ucl.ac.uk>,
    Department of Genetics, Evolution and Environment,
    University College London, Gower Street, London WC1E 6BT, England
*/

#include "bpp-tools.h"

/*         0         1         2         3         4         5         6         7          */
/*         01234567890123456789012345678901234567890123456789012345678901234567890123456789 */

/* ---- top-level help ---------------------------------------------------- */

void cmd_help(void)
{
  fprintf(stdout,
          "Usage: %s [OPTIONS]\n", progname);

  fprintf(stdout,
          "\n"
          "General options:\n"
          "  --help                  Display help information\n"
          "  --help=COMMAND          Display detailed help for a specific command\n"
          "  --version               Display version information\n"
          "  --quiet                 Only output warnings and fatal errors to stderr\n"
          "\n"
          "Alignment information:\n"
          "  --info                  Display PHYLIP alignment information\n"
          "\n"
          "Alignment manipulation:\n"
          "  --concat                Concatenate multi-locus alignments\n"
          "  --explode               Split a multi-locus PHYLIP into per-locus files\n"
          "  --extract               Extract sequences by label/tag/species\n"
          "  --remove CSV            Remove sequences by label/tag/species\n"
          "  --keep-loci CSV         Keep selected loci from a multi-locus PHYLIP\n"
          "  --drop-loci CSV         Drop selected loci from a multi-locus PHYLIP\n"
          "  --compress              Pattern-compress an alignment (writes P [MODEL] PHYLIP)\n"
          "\n"
          "Introgression / hybridization analyses:\n"
          "  --dstat CSV             ABBA-BABA D-statistic for a quartet\n"
          "  --hyde CSV              HyDe hybridization detection for a quartet\n"
          "  --fbranch               f-branch statistic on a species tree\n"
          "\n"
          "Common input modifiers:\n"
          "  --msa FILENAME          Multi-locus PHYLIP alignment file\n"
          "  --map FILENAME          Individual-to-species mapping file\n"
          "  --treefile FILENAME     File containing one or more Newick trees\n"
          "  --outgroup LABEL        Outgroup tip label (fbranch)\n"
          "  --nachar CHAR           Character used for missing data (concat)\n"
          "                          (use `-` as FILENAME to read from stdin;\n"
          "                          at most one input flag per invocation)\n"
          "\n"
          "Output modifiers:\n"
          "  --out FILENAME          Output file (default varies by command;\n"
          "                          compress / keep-loci / drop-loci default\n"
          "                          to stdout. Use `-` as FILENAME to write\n"
          "                          to stdout explicitly. Not supported for\n"
          "                          --explode (where --out is a template).)\n"
          "\n"
          "Info modifiers:\n"
          "  --per-locus             Show every locus row (no truncation)\n"
          "  --show-labels           Append the full unique-label inventory\n"
          "\n"
          "Compression modifiers:\n"
          "  --model MODEL           Compression model for --compress: JC69 | GTR\n"
          "                          (case-insensitive; default is GTR)\n"
          "\n"
          "Filter lists (for --extract / --remove):\n"
          "  --label_list CSV        Comma-separated sequence labels to match\n"
          "  --tag_list CSV          Comma-separated sequence tags to match\n"
          "  --species_list CSV      Comma-separated species names to match\n"
          "                          (requires --map)\n"
          "\n"
          "Statistical / resampling modifiers:\n"
          "  --bscount INT           Number of bootstrap replicates (default: 1000)\n"
          "  --alpha REAL            Significance level for confidence intervals "
                                    "(default: 0.05)\n"
          "  --jackknife INT         Jackknife block size\n"
          "  --seed INT              RNG seed (default: time-based)\n"
          "\n"
          "Other modifiers:\n"
          "  --debug[=LEVEL]         Print debug information (level is optional)\n"
          "  --verbose[=LEVEL]       Increase verbosity (level is optional)\n"
          "  --all                   Run all 6 taxa permutations (dstat)\n"
          "\n"
          "Available commands (use --help=COMMAND for details):\n"
          "  info       concat     explode    extract    remove\n"
          "  compress   keep-loci  drop-loci  dstat      hyde       fbranch\n"
          "\n"
          "Examples:\n"
          "  %s --info --msa aln.phy\n"
          "  %s --dstat P1,P2,P3,O --msa aln.phy --map map.txt\n"
          "  %s --help=dstat\n",
          progname, progname, progname);
}

/* ---- per-command help -------------------------------------------------- */

static void help_info(void)
{
  fprintf(stdout,
    "Command: --info\n"
    "\n"
    "  Display a structured summary of a multi-locus PHYLIP alignment.\n"
    "\n"
    "  Output sections:\n"
    "    - File header: input path and format (uncompressed, pattern-\n"
    "      compressed, or mixed).\n"
    "    - Per-locus summary table: index, sequence count, length, and\n"
    "      `miss%%` (fraction of cells with characters '-', '?', 'N',\n"
    "      'X' — case-insensitive — i.e. the AMAS-style missing-data\n"
    "      rate). For pattern-compressed input two extra columns\n"
    "      appear before `miss%%`: `model` (JC69 or GTR) and\n"
    "      `weights_sum` (the original site count represented by that\n"
    "      locus). When the input has more than 10 loci the table is\n"
    "      truncated to first 5 + last 5 with an omission marker; pass\n"
    "      --per-locus to show every row.\n"
    "    - Aggregate section: total alignments, total sequences, unique\n"
    "      sequence labels, per-locus length statistics (mean / min /\n"
    "      max), the total site or pattern counts (and a compression\n"
    "      ratio when compressed), and the overall missing-cells rate\n"
    "      (weighted across patterns so the percentage matches the\n"
    "      unpacked alignment).\n"
    "\n"
    "  Both uncompressed PHYLIP and pattern-compressed PHYLIP (P JC69 /\n"
    "  P GTR) input are accepted, including mixed multi-locus files.\n"
    "\n"
    "Syntax:\n"
    "  bpp-tools --info --msa FILENAME [--per-locus] [--show-labels]\n"
    "\n"
    "Required arguments:\n"
    "  --msa FILENAME          Multi-locus PHYLIP alignment file.\n"
    "\n"
    "Optional arguments:\n"
    "  --per-locus             Print every locus row instead of the\n"
    "                          first-5/last-5 truncated table.\n"
    "  --show-labels           Append the full sorted list of unique\n"
    "                          sequence labels found across all loci.\n"
    "\n"
    "Examples:\n"
    "  bpp-tools --info --msa aln.phy\n"
    "  bpp-tools --info --msa large.phy --per-locus\n"
    "  bpp-tools --info --msa aln.phy --show-labels\n"
    "  bpp-tools --info --msa compressed.phy\n"
  );
}

static void help_concat(void)
{
  fprintf(stdout,
    "Command: --concat\n"
    "\n"
    "  Concatenate the loci of a multi-locus PHYLIP alignment into a single\n"
    "  alignment, filling missing species with a placeholder character.\n"
    "\n"
    "  Rejects pattern-compressed input (weight-vector concatenation is not\n"
    "  supported; `--compress` can be used on the uncompressed output).\n"
    "\n"
    "Syntax:\n"
    "  bpp-tools --concat --msa FILENAME [options]\n"
    "\n"
    "Required arguments:\n"
    "  --msa FILENAME          Multi-locus PHYLIP alignment file.\n"
    "\n"
    "Optional arguments:\n"
    "  --out FILENAME          Output file (default: <msa>.concat.txt).\n"
    "  --nachar CHAR           Missing-data placeholder (default: '?').\n"
    "\n"
    "Examples:\n"
    "  bpp-tools --concat --msa multilocus.phy --out combined.phy\n"
    "  bpp-tools --concat --msa multilocus.phy --nachar N\n"
  );
}

static void help_explode(void)
{
  fprintf(stdout,
    "Command: --explode\n"
    "\n"
    "  Split a multi-locus PHYLIP file into one file per locus.\n"
    "\n"
    "  The output filenames are <template>.1, <template>.2, ... where the\n"
    "  template defaults to the input file path or is taken from --out.\n"
    "  Rejects pattern-compressed input.\n"
    "\n"
    "Syntax:\n"
    "  bpp-tools --explode --msa FILENAME [options]\n"
    "\n"
    "Required arguments:\n"
    "  --msa FILENAME          Multi-locus PHYLIP alignment file.\n"
    "\n"
    "Optional arguments:\n"
    "  --out FILENAME          Output filename template (default: input file).\n"
    "\n"
    "Examples:\n"
    "  bpp-tools --explode --msa multilocus.phy\n"
    "  bpp-tools --explode --msa multilocus.phy --out locus\n"
  );
}

static void help_extract(void)
{
  fprintf(stdout,
    "Command: --extract\n"
    "\n"
    "  Extract sequences from a PHYLIP alignment that match one or more of\n"
    "  the given label/tag/species filters.\n"
    "\n"
    "  At least one of --label_list / --tag_list / --species_list must be\n"
    "  specified. Using --species_list also requires --map. Rejects\n"
    "  pattern-compressed input.\n"
    "\n"
    "Syntax:\n"
    "  bpp-tools --extract --msa FILENAME (filter) [options]\n"
    "\n"
    "Required arguments:\n"
    "  --msa FILENAME          PHYLIP alignment file.\n"
    "\n"
    "Filter options (at least one required):\n"
    "  --label_list CSV        Comma-separated sequence labels to match.\n"
    "  --tag_list CSV          Comma-separated sequence tags to match.\n"
    "  --species_list CSV      Comma-separated species names (requires --map).\n"
    "\n"
    "Optional arguments:\n"
    "  --map FILENAME          Individual-to-species mapping file.\n"
    "  --out FILENAME          Output file (default: printed to stdout/auto).\n"
    "\n"
    "Examples:\n"
    "  bpp-tools --extract --msa aln.phy --label_list Sp1,Sp2\n"
    "  bpp-tools --extract --msa aln.phy --species_list Adig,Agre "
    "--map map.txt\n"
  );
}

static void help_remove(void)
{
  fprintf(stdout,
    "Command: --remove CSV\n"
    "\n"
    "  Remove sequences matching a label/tag/species filter from a PHYLIP\n"
    "  alignment and write the remainder.\n"
    "\n"
    "  Rejects pattern-compressed input.\n"
    "\n"
    "Syntax:\n"
    "  bpp-tools --remove CSV --msa FILENAME [options]\n"
    "\n"
    "Required arguments:\n"
    "  --remove CSV            Comma-separated sequences/labels to remove.\n"
    "  --msa FILENAME          PHYLIP alignment file.\n"
    "\n"
    "Optional arguments:\n"
    "  --out FILENAME          Output file.\n"
    "\n"
    "Examples:\n"
    "  bpp-tools --remove Sp1,Sp2 --msa aln.phy --out filtered.phy\n"
  );
}

static void help_keep_loci(void)
{
  fprintf(stdout,
    "Command: --keep-loci CSV\n"
    "\n"
    "  Read a multi-locus PHYLIP alignment and write out a new multi-locus\n"
    "  alignment containing only the loci named in the spec.\n"
    "\n"
    "  The spec is a comma-separated list of 1-indexed locus numbers and\n"
    "  hyphenated ranges (e.g. `1,2,5,8,20-100,7`). Numbers may be unsorted.\n"
    "  Duplicates and overlapping ranges are silently merged. Out-of-range\n"
    "  indices are fatal. The output preserves the original locus order.\n"
    "\n"
    "  Pattern-compressed input is preserved per-block (the compressed\n"
    "  writer is used for any block whose `pattern_weights` is set).\n"
    "\n"
    "Syntax:\n"
    "  bpp-tools --keep-loci CSV --msa FILENAME [--out FILENAME]\n"
    "\n"
    "Required arguments:\n"
    "  --keep-loci CSV         Spec for which loci to keep.\n"
    "  --msa FILENAME          Multi-locus PHYLIP alignment file.\n"
    "\n"
    "Optional arguments:\n"
    "  --out FILENAME          Output file (default: stdout).\n"
    "\n"
    "Examples:\n"
    "  bpp-tools --keep-loci 1,3,5 --msa multilocus.phy --out subset.phy\n"
    "  bpp-tools --keep-loci 20-100 --msa multilocus.phy\n"
    "  bpp-tools --keep-loci 1,2,5,8,20-100,7 --msa multilocus.phy\n"
  );
}

static void help_drop_loci(void)
{
  fprintf(stdout,
    "Command: --drop-loci CSV\n"
    "\n"
    "  Read a multi-locus PHYLIP alignment and write out a new multi-locus\n"
    "  alignment with the loci named in the spec removed.\n"
    "\n"
    "  The spec format and parsing rules are identical to `--keep-loci`:\n"
    "  comma-separated 1-indexed locus numbers and hyphenated ranges, in\n"
    "  any order, with overlapping/duplicate ranges silently merged.\n"
    "  Out-of-range indices are fatal. The output preserves the original\n"
    "  locus order of the loci that survive.\n"
    "\n"
    "  This command is the complement of `--keep-loci`: running both with\n"
    "  the same spec on the same input produces non-overlapping outputs\n"
    "  whose locus counts sum to the input's locus count.\n"
    "\n"
    "Syntax:\n"
    "  bpp-tools --drop-loci CSV --msa FILENAME [--out FILENAME]\n"
    "\n"
    "Required arguments:\n"
    "  --drop-loci CSV         Spec for which loci to drop.\n"
    "  --msa FILENAME          Multi-locus PHYLIP alignment file.\n"
    "\n"
    "Optional arguments:\n"
    "  --out FILENAME          Output file (default: stdout).\n"
    "\n"
    "Examples:\n"
    "  bpp-tools --drop-loci 1,3,5 --msa multilocus.phy --out remainder.phy\n"
    "  bpp-tools --drop-loci 20-100 --msa multilocus.phy\n"
  );
}

static void help_compress(void)
{
  fprintf(stdout,
    "Command: --compress\n"
    "\n"
    "  Compress a PHYLIP alignment into pattern-compressed form.\n"
    "\n"
    "  Identical alignment columns are collapsed into a single pattern with\n"
    "  an attached integer weight. The output uses the header syntax\n"
    "  `<count> <length> P <MODEL>` and a trailing whitespace-separated\n"
    "  weights line. Multi-locus input is compressed locus-by-locus, with\n"
    "  each block separated by a blank line.\n"
    "\n"
    "  Two compression models are supported:\n"
    "    - GTR (default): byte-exact column match.\n"
    "    - JC69:          site-local allele renumbering before matching,\n"
    "                     giving more aggressive deduplication on\n"
    "                     well-behaved data.\n"
    "\n"
    "  Already-compressed input is rejected.\n"
    "\n"
    "Syntax:\n"
    "  bpp-tools --compress --msa FILENAME [options]\n"
    "\n"
    "Required arguments:\n"
    "  --msa FILENAME          PHYLIP alignment file to compress.\n"
    "\n"
    "Optional arguments:\n"
    "  --model MODEL           Compression model: JC69 | GTR "
                              "(case-insensitive;\n"
    "                          default: GTR).\n"
    "  --out FILENAME          Output file (default: stdout).\n"
    "\n"
    "Examples:\n"
    "  bpp-tools --compress --msa aln.phy --out aln.comp.phy\n"
    "  bpp-tools --compress --model JC69 --msa aln.phy --out aln.jc69.phy\n"
    "  bpp-tools --compress --msa aln.phy | head -1\n"
  );
}

static void help_dstat(void)
{
  fprintf(stdout,
    "Command: --dstat P1,P2,P3,O\n"
    "\n"
    "  Compute Patterson's D-statistic (ABBA-BABA test) for a quartet of\n"
    "  taxa.\n"
    "\n"
    "  D = (fABBA - fBABA) / (fABBA + fBABA)\n"
    "\n"
    "  Significant departure from 0 suggests gene flow between P3 and\n"
    "  either P1 (D<0) or P2 (D>0). Confidence intervals are produced via\n"
    "  both bootstrap (over loci and sites) and delete-one-locus\n"
    "  jackknife. GTR pattern-compressed input is accepted, including\n"
    "  multi-locus files where every block is compressed with the same\n"
    "  model: the bootstrap draws sites from the uncompressed alignment\n"
    "  by sampling patterns proportionally to their weights, giving the\n"
    "  same bootstrap distribution as a site-level bootstrap on the\n"
    "  unpacked data. JC69-compressed input is rejected, as are mixed\n"
    "  files (JC69+GTR or compressed+uncompressed).\n"
    "\n"
    "Syntax:\n"
    "  bpp-tools --dstat P1,P2,P3,O --msa FILENAME --map FILENAME [options]\n"
    "\n"
    "Required arguments:\n"
    "  --dstat CSV             Four species labels in order P1,P2,P3,O "
                              "where\n"
    "                          O is the outgroup.\n"
    "  --msa FILENAME          Multi-locus PHYLIP alignment file.\n"
    "  --map FILENAME          Individual-to-species mapping file.\n"
    "\n"
    "Optional arguments:\n"
    "  --bscount INT           Bootstrap replicates (default: 1000).\n"
    "  --alpha REAL            Significance level for CI (default: 0.05).\n"
    "  --seed INT              RNG seed for resampling.\n"
    "  --all                   Run all 6 permutations of the three "
                              "ingroup taxa.\n"
    "\n"
    "Examples:\n"
    "  bpp-tools --dstat Adig,Agre,Agra,Smic --msa aln.phy --map map.txt\n"
    "  bpp-tools --dstat P1,P2,P3,O --msa aln.phy --map map.txt --bscount "
    "10000\n"
  );
}

static void help_hyde(void)
{
  fprintf(stdout,
    "Command: --hyde P1,P2,P3,O\n"
    "\n"
    "  Estimate the admixture proportion gamma for a hybrid quartet using\n"
    "  the HyDe (Hybridization Detection) estimator.\n"
    "\n"
    "  gamma = (fBBAA - fBABA) / (fBBAA - 2*fBABA + fABBA)\n"
    "\n"
    "  A value of gamma near 0.5 is consistent with a first-generation\n"
    "  hybrid between P1 and P3; gamma near 0 or 1 indicates no or\n"
    "  complete introgression. GTR pattern-compressed input is accepted;\n"
    "  JC69-compressed input is rejected.\n"
    "\n"
    "Syntax:\n"
    "  bpp-tools --hyde P1,P2,P3,O --msa FILENAME --map FILENAME [options]\n"
    "\n"
    "Required arguments:\n"
    "  --hyde CSV              Four species labels: P1,P2,P3,O.\n"
    "  --msa FILENAME          Multi-locus PHYLIP alignment file.\n"
    "  --map FILENAME          Individual-to-species mapping file.\n"
    "\n"
    "Optional arguments:\n"
    "  --bscount INT           Bootstrap replicates (default: 1000).\n"
    "  --alpha REAL            Significance level (default: 0.05).\n"
    "  --seed INT              RNG seed.\n"
    "\n"
    "Examples:\n"
    "  bpp-tools --hyde Sp1,Sp2,Sp3,Out --msa aln.phy --map map.txt\n"
  );
}

static void help_fbranch(void)
{
  fprintf(stdout,
    "Command: --fbranch\n"
    "\n"
    "  Compute the f-branch statistic (Malinsky et al. 2018) that assigns\n"
    "  introgression signals to specific branches of a species tree.\n"
    "\n"
    "  For each internal branch b with sibling branch a, and each external\n"
    "  tip C (not in a or b and not the outgroup), the statistic\n"
    "\n"
    "    f_b(C) = median_A [ min_B [ gamma(A, B, C, O) ] ]\n"
    "\n"
    "  is reported, where A ranges over tips descending from a and B over\n"
    "  tips descending from b. Produces a matrix (branches x external\n"
    "  tips) and a PDF visualization of the tree with node indices.\n"
    "  GTR pattern-compressed input is accepted; JC69-compressed input\n"
    "  is rejected.\n"
    "\n"
    "Syntax:\n"
    "  bpp-tools --fbranch --treefile FILENAME --msa FILENAME\n"
    "                      --map FILENAME --outgroup LABEL\n"
    "\n"
    "Required arguments:\n"
    "  --treefile FILENAME     File containing one or more Newick trees.\n"
    "  --msa FILENAME          Multi-locus PHYLIP alignment file.\n"
    "  --map FILENAME          Individual-to-species mapping file.\n"
    "  --outgroup LABEL        Species label to use as outgroup.\n"
    "\n"
    "Examples:\n"
    "  bpp-tools --fbranch --treefile trees.txt --msa aln.phy \\\n"
    "            --map map.txt --outgroup Smic\n"
  );
}

/* ---- dispatcher -------------------------------------------------------- */

typedef struct
{
  const char * name;
  void (*func)(void);
} help_entry_t;

static const help_entry_t help_commands[] =
{
  { "info",      help_info     },
  { "concat",    help_concat   },
  { "explode",   help_explode  },
  { "extract",   help_extract  },
  { "remove",    help_remove   },
  { "keep-loci", help_keep_loci },
  { "drop-loci", help_drop_loci },
  { "compress",  help_compress },
  { "dstat",     help_dstat    },
  { "hyde",      help_hyde     },
  { "fbranch",   help_fbranch  },
  { NULL,        NULL          }
};

void cmd_help_command(const char * cmd)
{
  const help_entry_t * entry;

  for (entry = help_commands; entry->name; entry++)
  {
    if (!strcmp(cmd, entry->name))
    {
      entry->func();
      return;
    }
  }

  fprintf(stderr,
          "Error: Unknown command '%s'.\n\n"
          "Available commands:\n", cmd);

  for (entry = help_commands; entry->name; entry++)
    fprintf(stderr, "  %s\n", entry->name);

  fprintf(stderr,
          "\nUsage: bpp-tools --help=COMMAND\n"
          "   or: bpp-tools --help COMMAND\n");
}
