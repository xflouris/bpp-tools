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

#ifdef _MSC_VER
#include "getopt_win.h"
#endif

static char * progname;
static char progheader[80];
char * cmdline;

/* global error message buffer */
__THREAD int bpp_errno;
__THREAD char bpp_errmsg[200] = {0};

/* options */
long opt_arch;
long opt_explode;
long opt_help;
long opt_quiet;
long opt_seed;
long opt_version;
char * opt_msafile;
char * opt_outfile;
char * opt_dstat;
char * opt_extract;
char * opt_remove;

long mmx_present;
long sse_present;
long sse2_present;
long sse3_present;
long ssse3_present;
long sse41_present;
long sse42_present;
long popcnt_present;
long avx_present;
long avx2_present;
long altivec_present;

static struct option long_options[] =
{
  {"help",         no_argument,       0, 0 },  /*  0 */
  {"version",      no_argument,       0, 0 },  /*  1 */
  {"quiet",        no_argument,       0, 0 },  /*  2 */
  {"msa",          required_argument, 0, 0 },  /*  3 */
  {"dstat",        required_argument, 0, 0 },  /*  4 */
  {"out",          required_argument, 0, 0 },  /*  5 */
  {"explode",      no_argument,       0, 0 },  /*  6 */
  {"extract",      required_argument, 0, 0 },  /*  7 */
  {"remove",       required_argument, 0, 0 },  /*  8 */
  { 0, 0, 0, 0 }
};

void args_init(int argc, char ** argv)
{
  int option_index = 0;
  int c;
  
  /* set defaults */

  progname = argv[0];

  opt_arch = -1;
  opt_dstat = NULL;
  opt_explode = 0;
  opt_help = 0;
  opt_msafile = NULL;
  opt_outfile = NULL;
  opt_quiet = 0;
  opt_seed = -1;
  opt_version = 0;


  while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) == 0)
  {
    switch (option_index)
    {
      case 0:
        opt_help = 1;
        break;

      case 1:
        opt_version = 1;
        break;

      case 2:
        opt_quiet = 1;
        break;

      case 3:
        opt_msafile = xstrdup(optarg);
        break;

      case 4:
        opt_dstat = xstrdup(optarg);
        break;

      case 5:
        opt_outfile = xstrdup(optarg);
        break;

      case 6:
        opt_explode = 1;
        break;

      case 7:
        opt_extract = xstrdup(optarg);
        break;

      case 8:
        opt_remove = xstrdup(optarg);
        break;


      default:
        fatal("Internal error in option parsing");
    }
  }

  if (c != -1)
    exit(EXIT_FAILURE);

  int commands  = 0;

  /* check for number of independent commands selected */
  if (opt_version)
    commands++;
  if (opt_help)
    commands++;
  if (opt_dstat)
    commands++;
  if (opt_explode)
    commands++;
  if (opt_extract)
    commands++;
  if (opt_remove)
    commands++;

  /* if more than one independent command, fail */
  if (commands > 1)
    fatal("More than one command specified");

  #if 0
  /* if no command specified, turn on --help */
  if (!commands)
  {
    opt_help = 1;
    return;
  }
  #endif
}

static void dealloc_switches()
{
  if (opt_dstat) free(opt_dstat);
  if (opt_msafile) free(opt_msafile);
  if (opt_outfile) free(opt_outfile);
  if (opt_extract) free(opt_extract);
  if (opt_remove) free(opt_remove);
}

void cmd_none()
{
  if (!opt_quiet)
    fprintf(stderr,
            "For help, please enter: %s --help\n"
            "\n"
            "For further details, please see the manual by entering: man bpp-tools\n"
            "\n"
            "Example commands:\n"
            "\n"
            "bpp-tools --explode --msa FILENAME --output FILENAME\n"
            "bpp-tools --extract CSV --msa FILENAME --output FILENAME\n"
            "bpp-tools --remove CSV --msa FILENAME --output FILENAME\n"
            "bpp-tools --subsample CSV --msa FILENAME --output FILENAME\n"
            "bpp-tools --dstat CSV --msa FILENAME\n"
            "\n",
            progname);
}

void cmd_help()
{
  /*         0         1         2         3         4         5         6         7          */
  /*         01234567890123456789012345678901234567890123456789012345678901234567890123456789 */

  fprintf(stderr,
          "Usage: %s [OPTIONS]\n", progname);
  fprintf(stderr,
          "\n"
          "General options:\n"
          "  --help             display help information\n"
          "  --version          display version information\n"
          "  --quiet            only output warnings and fatal errors to stderr\n"
          "  --dstat taxa       run dstatistics\n"
          "\n"
         );

  /*         0         1         2         3         4         5         6         7          */
  /*         01234567890123456789012345678901234567890123456789012345678901234567890123456789 */
}

void getentirecommandline(int argc, char * argv[])
{
  int len = 0;
  int i;

  for (i = 0; i < argc; ++i)
    len += strlen(argv[i]);

  cmdline = (char *)xmalloc((size_t)(len + argc + 1));
  cmdline[0] = 0;

  for (i = 0; i < argc; ++i)
  {
    strcat(cmdline, argv[i]);
    strcat(cmdline, " ");
  }
}

void fillheader()
{
  snprintf(progheader, 80,
           "%s %s_%s, %1.fGB RAM, %ld cores",
           PROG_NAME, PROG_VERSION, PROG_ARCH,
           arch_get_memtotal() / 1024.0 / 1024.0 / 1024.0,
           arch_get_cores());
}

void show_header()
{
  fprintf(stdout, "%s\n", progheader);
  fprintf(stdout, "https://github.com/xflouris/bpp-tools\n");
  fprintf(stdout,"\n");
}

int main (int argc, char * argv[])
{
  fillheader();
  getentirecommandline(argc, argv);

  args_init(argc, argv);

  show_header();

  cpu_features_detect();
  cpu_features_show();
  if (!opt_version && !opt_help)
    cpu_setarch();

  if (opt_help)
  {
    cmd_help();
  }
  else if (opt_version)
  {
    ;
  }
  else if (opt_dstat)
  {
    cmd_dstat();
  }
  else if (opt_explode)
  {
    cmd_explode();
  }
  else if (opt_extract)
  {
    cmd_extract();
  }
  else if (opt_remove)
  {
    cmd_remove();
  }
  else
    cmd_none();

  dealloc_switches();
  free(cmdline);
  return (0);
}
