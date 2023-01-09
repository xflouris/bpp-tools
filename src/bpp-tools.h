/*
    Copyright (C) 2021-2023 Tomas Flouri

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

#include <assert.h>
#include <fcntl.h>
#include <search.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <sys/stat.h>
#include <stdint.h>
#include <inttypes.h>
#include <ctype.h>
#include <pthread.h>

#ifdef _MSC_VER
#include <pmmintrin.h>
#include <immintrin.h>
#else
#include <x86intrin.h>
#endif

#ifndef _MSC_VER
#include <getopt.h>
#endif

#ifndef _MSC_VER
#include <sys/time.h>
#include <unistd.h>
#endif

/* platform specific */

#if (defined(__BORLANDC__) || defined(_MSC_VER))
#define __THREAD __declspec(thread)
#else
#define __THREAD __thread
#endif

#ifdef _MSC_VER
#define PLL_ALIGN_HEADER(X) __declspec(align(X))
#define PLL_ALIGN_FOOTER(X)
#else
#define PLL_ALIGN_HEADER(X)
#define PLL_ALIGN_FOOTER(X) __attribute__((aligned(X)))
#endif

#ifndef _MSC_VER
#define xasprintf asprintf
#endif

#ifdef _MSC_VER
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif

/* constants */

#define PROG_NAME "bpp-tools"

#define PLL_STRING(x) #x
#define PLL_C2S(x) PLL_STRING(x)


#define VERSION_MAJOR 0
#define VERSION_MINOR 1
#define VERSION_PATCH 0

#define PROG_VERSION "v" PLL_C2S(VERSION_MAJOR) "." PLL_C2S(VERSION_MINOR) "." \
        PLL_C2S(VERSION_PATCH)

#ifdef __PPC__

#ifdef __LITTLE_ENDIAN__
#define PROG_CPU "ppc64le"
#else
#error "Big endian ppc64 CPUs not supported"
#endif

#else

#define PROG_CPU "x86_64"

#endif

#ifdef __APPLE__
#define PROG_OS "osx"
#include <sys/resource.h>
#include <sys/sysctl.h>
#endif

#ifdef __linux__
#define PROG_OS "linux"
#include <sys/resource.h>
#include <sys/sysinfo.h>
#endif

#ifdef _WIN32
#define PROG_OS "win"
#include <windows.h>
#include <psapi.h>
#endif

#define PROG_ARCH PROG_OS "_" PROG_CPU

#define BPP_FAILURE  0
#define BPP_SUCCESS  1

#define LINEALLOC 2048
#define ASCII_SIZE 256

#define BPP_DATA_DNA                    0
#define BPP_DATA_AA                     1

/* error codes */

#define ERROR_PHYLIP_SYNTAX            106
#define ERROR_PHYLIP_LONGSEQ           107
#define ERROR_PHYLIP_NONALIGNED        108
#define ERROR_PHYLIP_ILLEGALCHAR       109
#define ERROR_PHYLIP_UNPRINTABLECHAR   110
#define ERROR_PARSE_MORETHANEXPECTED   111
#define ERROR_PARSE_LESSTHANEXPECTED   112
#define ERROR_PARSE_INCORRECTFORMAT    113

/* libpll related definitions */

#define PLL_ALIGNMENT_CPU               8
#define PLL_ALIGNMENT_SSE              16
#define PLL_ALIGNMENT_AVX              32

#define PLL_ATTRIB_ARCH_CPU            0
#define PLL_ATTRIB_ARCH_SSE       (1 << 0)
#define PLL_ATTRIB_PATTERN_TIP    (1 << 4)

#define PLL_ATTRIB_ARCH_AVX       (1 << 1)
#define PLL_ATTRIB_ARCH_AVX2      (1 << 2)
#define PLL_ATTRIB_ARCH_AVX512    (1 << 3)
#define PLL_ATTRIB_ARCH_MASK         0xF

/* structures and data types */

typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;

typedef struct msa_s
{
  int count;
  int length;

  char ** sequence;
  char ** label;

  int amb_sites_count;
  int original_length;

  double * freqs;

  int dtype;
  int model;
  int original_index;

} msa_t;

typedef struct phylip_s
{
  FILE * fp;
  char * line;
  size_t line_size;
  size_t line_maxsize;
  char buffer[LINEALLOC];
  const unsigned int * chrstatus;
  long no;
  long filesize;
  long lineno;
  long stripped_count;
  long stripped[256];
} phylip_t;

typedef struct list_item_s
{
  void * data;
  struct list_item_s * next;
} list_item_t;

typedef struct list_s
{
  list_item_t * head;
  list_item_t * tail;
  long count;
} list_t;


typedef struct ht_item_s
{
  unsigned long key;
  void * value;
} ht_item_t;

typedef struct hashtable_s
{
  unsigned long table_size;
  unsigned long entries_count;
  list_t ** entries;
} hashtable_t;

typedef struct pair_s
{
  char * label;
  void * data;
} pair_t;

/* macros */

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifdef _MSC_VER
#define SWAP(x,y) do                                                  \
  {                                                                   \
    size_t s = MAX(sizeof(x),sizeof(y));                              \
    unsigned char * temp = (unsigned char *)malloc(s*sizeof(char));   \
    memcpy(temp,&y,s);                                                \
    memcpy(&y,&x,s);                                                  \
    memcpy(&x,temp,s);                                                \
    free(temp);                                                       \
  } while(0)
#else
#define SWAP(x,y) do { __typeof__ (x) _t = x; x = y; y = _t; } while(0)
#endif

#ifdef _MSC_VER
#define PLL_POPCOUNT pll_popcount
#define PLL_POPCOUNTL pll_popcount64
#define PLL_CTZ pll_ctz
#define xtruncate _chsize
#else
#define PLL_POPCOUNT __builtin_popcount
#define PLL_POPCOUNTL __builtin_popcountl
#define PLL_CTZ __builtin_ctz
#define xtruncate ftruncate
#endif

/* options */

extern long opt_arch;
extern long opt_explode;
extern long opt_help;
extern long opt_quiet;
extern long opt_seed;
extern long opt_version;
extern char * cmdline;
extern char * opt_msafile;
extern char * opt_outfile;
extern char * opt_dstat;
extern char * opt_extract;

/* common data */

extern __THREAD int bpp_errno;
extern __THREAD char bpp_errmsg[200];

extern const unsigned int pll_map_nt[256];
extern const unsigned int pll_map_nt_tcag[256];
extern const unsigned int pll_map_aa[256];
extern const unsigned int pll_map_fasta[256];
extern const unsigned int pll_map_amb[256];
extern const unsigned int pll_map_validjc69[16];
extern const unsigned int bpp_tolower_table[256];
extern const unsigned int pll_map_nt_missing[256];
extern const unsigned int pll_map_aa_missing[256];

extern long mmx_present;
extern long sse_present;
extern long sse2_present;
extern long sse3_present;
extern long ssse3_present;
extern long sse41_present;
extern long sse42_present;
extern long popcnt_present;
extern long avx_present;
extern long avx2_present;
extern long altivec_present;

/* functions in phylip.c */

phylip_t * phylip_open(const char * filename,
                       const unsigned int * map);

int phylip_rewind(phylip_t * fd);

void phylip_close(phylip_t * fd);

msa_t * phylip_parse_interleaved(phylip_t * fd);

msa_t * phylip_parse_sequential(phylip_t * fd);

msa_t ** phylip_parse_multisequential(phylip_t * fd, long * count);

void phylip_print(FILE * fp, const msa_t * msa);

/* functions in util.c */

#ifdef _MSC_VER
__declspec(noreturn) void fatal(const char * format, ...);
int xasprintf(char ** strp, const char * fmt, ...);
#else
void fatal(const char * format, ...) __attribute__ ((noreturn));
#endif
void progress_init(const char * prompt, unsigned long size);
void progress_update(unsigned long progress);
void progress_done(void);
void * xmalloc(size_t size);
void * xcalloc(size_t nmemb, size_t size);
void * xrealloc(void *ptr, size_t size);
char * xstrchrnul(char *s, int c);
char * xstrdup(const char * s);
char * xstrndup(const char * s, size_t len);
long getusec(void);
FILE * xopen(const char * filename, const char * mode);
void * pll_aligned_alloc(size_t size, size_t alignment);
void pll_aligned_free(void * ptr);
int xtolower(int c);

/* functions in arch.c */

uint64_t arch_get_memused(void);

uint64_t arch_get_memtotal(void);

long arch_get_cores(void);

/* functions in msa.c */

void msa_print_phylip(FILE * fp,
                      msa_t ** msa,
                      long count,
                      unsigned int ** weights);

void msa_destroy(msa_t * msa);

int msa_remove_ambiguous(msa_t * msa);

void msa_count_ambiguous_sites(msa_t * msa, const unsigned int * map);

int msa_remove_missing_sequences(msa_t * msa);

/* functions in dstat.c */

void cmd_dstat(void);

/* functions in explode.c */

void cmd_explode(void);

/* functions in hardware.c */

void cpu_features_show(void);

void cpu_features_detect(void);

void cpu_setarch(void);

#ifdef _MSC_VER
int pll_ctz(unsigned int x);
unsigned int pll_popcount(unsigned int x);
unsigned int pll_popcount64(unsigned long x);
#endif

/* functions in hash.c */

void * hashtable_find(hashtable_t * ht,
                      void * x,
                      unsigned long hash,
                      int (*cb_cmp)(void *, void *));

hashtable_t * hashtable_create(unsigned long items_count);

int hashtable_strcmp(void * x, void * y);

int hashtable_ptrcmp(void * x, void * y);

unsigned long hash_djb2a(char * s);

unsigned long hash_fnv(char * s);

int hashtable_insert(hashtable_t * ht,
                     void * x,
                     unsigned long hash,
                     int (*cb_cmp)(void *, void *));

void hashtable_insert_force(hashtable_t * ht,
                            void * x,
                            unsigned long hash);

void hashtable_destroy(hashtable_t * ht, void (*cb_dealloc)(void *));

int cb_cmp_pairlabel(void * a, void * b);

/* functions in list.c */

void list_append(list_t * list, void * data);

void list_prepend(list_t * list, void * data);

void list_clear(list_t * list, void (*cb_dealloc)(void *));

long list_reposition_tail(list_t * list, list_item_t * item);
long list_delitem(list_t * list, list_item_t * item, void (*cb_dealloc)(void *));

/* functions in extract.c */
void cmd_extract();
