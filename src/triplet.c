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

#define NT_CHARS 4

/* biallelic only */
static const unsigned int map_solid_and_biallelic[128] =
{
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  1,  0,  1,  0,  0,  0,  1,  0,  0,  0,  1,  0,  1,  0,  0,
  0,  0,  1,  1,  1,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,
  0,  1,  0,  1,  0,  0,  0,  1,  0,  0,  0,  1,  0,  1,  0,  0,
  0,  0,  1,  1,  1,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0
};

static const double allele_weights[128][4] =
{
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  1,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   1,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   1,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0, 0.5, 0.5},
  {  0,   0,   0,   0},
  {0.5, 0.5,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {0.5,   0, 0.5,   0},
  {  0, 0.5, 0.5,   0},
  {  0,   0,   0,   1},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {0.5,   0,   0, 0.5},
  {  0,   0,   0,   0},
  {  0, 0.5,   0, 0.5},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  1,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   1,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   1,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0, 0.5, 0.5},
  {  0,   0,   0,   0},
  {0.5, 0.5,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {0.5,   0, 0.5,   0},
  {  0, 0.5, 0.5,   0},
  {  0,   0,   0,   1},
  {  0,   0,   0,   0},
  {  0,   0,   0,   0},
  {0.5,   0,   0, 0.5},
  {  0,   0,   0,   0},
  {  0, 0.5,   0, 0.5}
};

static long getnumdigits(long n)
{
   if (n == 0) return 1;
   return (long)(floor(log10(labs(n)) + 1));
}

static int cb_cmp_double(const void* a, const void* b)
{
   double* x = (double*)a;
   double* y = (double*)b;

   if (*x > *y) return 1;
   if (*x < *y) return -1;

   return 0;
}

/* given the indices of four species (species), we condense the alignment (msa)
   by collecting all sequences from those species using the index structure
   and the maparray

   Params:
     msa      : the concatenated alignment
     maparray : the list of mappings (individual->species)
     index    : points to species element for given sequence in msa (0 to 3)
     species  : the 4 species to be used for d, points to maparray indices
     vec      : stores prob of each base for each site (size: 4x4*sites doubles)

*/
/* TODO: Get species tags */

static char invcharmap_acgt[16] =
{
  '-', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
  'T', 'W', 'Y', 'H', 'K', 'D', 'B', '-'
};

static long ends_with_hat_str(const char* x, const char* suffix)
{
   size_t xlen, slen;
   const char* p;
   const char* hat;

   xlen = strlen(x);
   slen = strlen(suffix);

   if (slen >= xlen)
      return 0;

   p = x + xlen - slen;
   hat = p - 1;
   if (!strcmp(p, suffix) && *hat == '^')
      return 1;

   return 0;
}

static time_t time_start;

static void timer_start()
{
   time_start = time(NULL);
}

static char * timer(char timestr[])
{
   time_t t;
   int h, m, s;

   t = time(NULL) - time_start;
   h = (int)t / 3600;
   m = (int)(t % 3600) / 60;
   s = (int)(t - (t / 60) * 60);
   if (h)
      sprintf(timestr, "%d:%02d:%02d", h, m, s);
   else
      sprintf(timestr, "%d:%02d", m, s);
   return(timestr);
}

static int get_ipop_for_sequences(int* ipop, msa_t* msa, list_t* maplist, char** species, int species_count)
{
   /* get ipop for sequences at locus into ipop */
   long seq, i;
   list_item_t* li;   /* item in list of sample names and species names */
   char* ind_name_in_seq, * sp;

   memset(ipop, -1, msa->count * sizeof(int));

   for (seq = 0; seq < msa->count; ++seq) {
      ind_name_in_seq = strchr(msa->label[seq], (int)'^');
      if (ind_name_in_seq == NULL) fatal("Use ^ to specify sample-label in seq %s.", msa->label[seq]);
      /* search in imap for individual name */
      li = maplist->head;
      while (li) {
         mapping_t* imap = (mapping_t*)(li->data);
         if (strcmp(ind_name_in_seq+1, imap->individual) == 0) {
            sp = imap->species;
            for (i = 0; i < species_count; i++)
               if (strcmp(sp, species[i]) == 0) break;
            assert(i < species_count);
            ipop[seq] = i;
         }
         li = li->next;
      }
      if (ipop[seq] == -1)
         fatal("sample %d in alignment, %s, not found in imap file", seq + 1, msa->label[seq]);
   }
   return 0;
}

void cmd_triplet()
{
   phylip_t* fd;
   msa_t** msa_list;
   long i, j, k, isp, locus, seq, site;
   long nloci, nspecies = 0, nsites_tot = 0;
   char* missing_species_sites, * missing_species_loci, ** species;
   char timestr[32] = "";
   double* species_site_vec;
   double mf[5] = {0}, sf[5] = { 0 };
   int* nsites_used;  /* nsites used differ among loci */
   int transform_triplet = 0, singleclean=1;

   /* change debug level */
   if(opt_debug) opt_debug = 3;
   timer_start();

   /* open phylip file */
   if (opt_debug)
      xdebug("Reading sequence alignments from %s\n", opt_msafile);
   fd = phylip_open(opt_msafile, pll_map_fasta);
   if (!fd)
      fatal("Cannot open file %s", opt_msafile);

   printf("Parsing sequence alignments %s .. %s\n", opt_msafile, timer(timestr));
   /* read alignment.  This takes both sequences and site patterns. */
   msa_list = phylip_parse_multisequential(fd, &nloci);
   assert(msa_list);
   phylip_close(fd);

   printf("Sequence alignments parsed %s .. %s\n", opt_msafile, timer(timestr));

   /* compress each locus in place */
   for (i = 0; i < nloci; ++i) {
      msa_t* m = msa_list[i];
      int len = m->length;
      unsigned int* w = NULL;
      if (msa_list[i]->pattern_weights == NULL) {
         w = compress_site_patterns(m->sequence, pll_map_nt, m->count, &len, COMPRESS_JC69, NULL);
         if (!w)
            fatal("compression failed for alignment %d", i + 1);
         m->pattern_weights = w;
         m->compress_model = COMPRESS_JC69;
         m->length = len;
      }
      if (opt_debug) {
         phylip_print_compressed(stdout, m);
         if (i < nloci - 1)
            fprintf(stdout, "\n");  /* blank line between multi-sequential blocks */
      }
   }
   if (!opt_debug) opt_debug = 1;

   nsites_used = (int*)xmalloc(nloci * sizeof(int));  

   if (!opt_mapfile)
      fatal("A map file needs to be specified with the --map option....");
   list_t* maplist = parse_mapfile(opt_mapfile);

   if (opt_debug)
      printf("compressed into site patterns.. %s\n", timer(timestr));

   /* collect distinct species names into species[] */
   list_item_t* li = maplist->head;
   char* sp;
   int nsamples = 0;
   while (li) {
      nsamples++;
      li = li->next;
   }
   species = (char**)xmalloc(nsamples * sizeof(char*));  /* nsamples > nspecies */
   li = maplist->head;
   nspecies = 0;
   while (li) {
      mapping_t* map = (mapping_t*)(li->data);
      sp = map->species;
      for (k = 0; k < nspecies; ++k)
         if (strcmp(sp, species[k]) == 0) break;
      if (k == nspecies)
         species[nspecies++] = xstrdup(sp);
      li = li->next;
   }
   if (nspecies < nsamples)
      species = (char**)xrealloc(species, nspecies * sizeof(char*));
   if (opt_debug) {
      printf("\n# of species found in the imap file: %2d\n", nspecies);
      for (k = 0; k < nspecies; ++k)
         printf("%d (%s)\n", k+1, species[k]);
   }
   for (locus = 0; locus < nloci; ++locus)
      nsites_tot += msa_list[locus]->length;

   missing_species_sites = (char*)xmalloc(nspecies * nsites_tot * sizeof(char));
   memset(missing_species_sites, 1, nspecies * nsites_tot * sizeof(char));
   missing_species_loci = (char*)xmalloc(nspecies * nloci * sizeof(char));
   memset(missing_species_loci, 1, nspecies * nloci * sizeof(char));
   species_site_vec = (double*)xmalloc(nspecies * nsites_tot * 4 * sizeof(double));
   memset(species_site_vec, 0, nspecies * nsites_tot * 4 * sizeof(double));

   /* generate species_site_vec[], and mf[] sf[] for a triplet */
   int* ipop = NULL;
   int locus_offset = 0;
   char b;
   double f[5] = {0}, s;

   int max_num_seqs = 0;
   for (locus = 0; locus < nloci; ++locus)
      if (max_num_seqs < msa_list[locus]->count)
         max_num_seqs = msa_list[locus]->count;
   ipop = (int*)xmalloc(max_num_seqs * sizeof(int));

   if (opt_debug)
      printf("\nConstructing species_site_vec for ACGT.. %s\n", timer(timestr));

   for (locus = 0; locus < nloci; locus_offset += msa_list[locus++]->length) {
      /* sp_labels = map_sequences_to_species(msa_list[locus], maplist, &index, &sp_seqcount, &sp_count); */
      get_ipop_for_sequences(ipop, msa_list[locus], maplist, species, nspecies);

      if (opt_debug > 1) {
         printf("\nipop for locus %ld (%ld seqs): ", locus + 1, msa_list[locus]->count);
         for (seq = 0; seq < msa_list[locus]->count; ++seq)
            printf(" %d", ipop[seq]);
      }
      for (site = 0; site < msa_list[locus]->length; ++site) {
         for (isp = 0; isp < nspecies; ++isp) {
            f[0] = f[1] = f[2] = f[3] = 0;  /* freqs for 4 bases at site from isp */
            /* averaging all sequences from the same species */
            for (seq = 0; seq < msa_list[locus]->count; ++seq) {
               if (isp != ipop[seq]) continue;
               b = msa_list[locus]->sequence[seq][site];
               if (map_solid_and_biallelic[b] == 0)
                  continue;
               missing_species_sites[isp * nsites_tot + locus_offset + site] = (char)0;
               missing_species_loci[isp * nloci + locus] = (char)0;
               for (i = 0; i < 4; i++) f[i] += allele_weights[b][i];
            }
            if (missing_species_sites[isp * nsites_tot + locus_offset + site])
               continue;
            s = f[0] + f[1] + f[2] + f[3];
            assert(s > 0);
            for (i = 0; i < 4; ++i) f[i] /= s;
            memcpy(species_site_vec + (isp * nsites_tot + locus_offset + site) * 4, f, 4 * sizeof(double));

            if (opt_debug > 2) {
               printf("\nspecies %ld (%s) locus %ld site %2ld ", isp + 1, species[isp], locus + 1, site + 1);
               printf(" [%d]", missing_species_sites[isp * nsites_tot + locus_offset + site]);
               printf(" ACGT %9.5f%9.5f%9.5f%9.5f", f[0], f[1], f[2], f[3]);
            }
         }
         if (opt_debug > 2) printf("\n");
      }
   }

   if (opt_debug)
      printf("\nspecies_site_vec for ACGT constructed.. %s", timer(timestr));

   int triplet[3] = { 0,1,2 }, itrip;
   int nloci_used, ipatt;
   double z, f01234[5] = {0};  /* frequencies for xxx xxy yxx xyx xyz for locus */

   for (itrip = 0; itrip < 1; itrip++) {
      if (opt_debug)
         printf("\ntriplet: %d (%s) %d (%s) %d (%s)\n", 
            triplet[0], species[triplet[0]], triplet[1], species[triplet[1]], triplet[2], species[triplet[2]]);
      locus_offset = 0;
      nloci_used = 0;
      mf[0] = mf[1] = mf[2] = mf[3] = 0;
      sf[0] = sf[1] = sf[2] = sf[3] = 0;
      memset(nsites_used, 0, nloci * sizeof(int));
      
      for (locus = 0; locus < nloci; locus_offset += msa_list[locus++]->length) {
         /* skip the locus if any species are missing */
         for (isp = 0; isp < 3; ++isp)
            if (missing_species_loci[triplet[isp] * nloci + locus]) break;
         if (isp < 3)
            continue;

         nsites_used[locus] = 0;
         f01234[0] = f01234[1] = f01234[2] = f01234[3] = 0;  
         for (site = 0; site < msa_list[locus]->length; ++site) {
            for (isp = 0; isp < 3; ++isp)
               if (missing_species_sites[triplet[isp] * nsites_tot + locus_offset + site]) break;
            if (isp < 3)
               continue;
            int weight = msa_list[locus]->pattern_weights[site];
            nsites_used[locus] += weight;

            if (singleclean) {
              int ipatt = 4;
              char z[3];
              for (j = 0; j < 3; j++)
                z[j] = msa_list[locus]->sequence[triplet[j]][site];
              if (z[0] != z[1] && z[0] != z[2] && z[1] != z[2]) ipatt = 4;
              else if (z[0] == z[1] && z[0] == z[2]) ipatt = 0;
              else if (z[0] == z[1]) ipatt = 1;
              else if (z[1] == z[2]) ipatt = 2;
              else if (z[0] == z[2]) ipatt = 3;
              f01234[ipatt] += weight;
            }
            else {
              double* sv[3];
              for (isp = 0; isp < 3; ++isp)
                sv[isp] = species_site_vec + (isp * nsites_tot + locus_offset + site) * 4;

              f[0] = f[1] = f[2] = f[3] = 0; /* frequencies xxx xxy yxx xyx xyz for site */
              for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                  for (k = 0; k < 4; k++) {
                    if (i != j && i != k && j != k) continue;
                    z = sv[0][i] * sv[1][j] * sv[2][k];
                    if (i == j && i == k)  f[0] += z;
                    else if (i == j)       f[1] += z;
                    else if (j == k)       f[2] += z;
                    else                   f[3] += z;
                  }
                }
              }
              for (i = 0; i < 4; i++) f01234[i] += weight * f[i];
            }
         }
         assert(nsites_used[locus]);
         for (i = 0; i < 4; i++) f01234[i] /= nsites_used[locus];
         f01234[4] = 1 - (f01234[0] + f01234[1] + f01234[2] + f01234[3]);

         for (i = 0; i < 5; i++) {
            double y = f01234[i], d;
            if (transform_triplet) {  /* this is necessary for transform only */
               if (y <= 0)
                  y = 0.5 / nsites_used[locus];
               else if (1 - y <= 0)
                  y = 1 - 0.5 / nsites_used[locus];
            }

            if (transform_triplet == 1)  /* logit: y=log(f) */
               y = log(y / (1 - y));
            else if (transform_triplet == 2) {
               y = log(y);
               /* y = (i >= 1 && i <= 4) ? log(y) : 0; */
            }

            d = y - mf[i];
            mf[i] += d / (nloci_used + 1.0);
            sf[i] += d * (y - mf[i]);
            /* sf[i] += d * d * nloci_used / (nloci_used + 1.0); this is equivalent */
         }

         nloci_used++;
         if (opt_debug > 2) {
            printf("\nlocus %3d pattern0123: %9.5f%9.5f%9.5f%9.5f", 
               locus + 1, f01234[0], f01234[1], f01234[2], f01234[3]);
         }
         if (nloci > 10000 && (locus + 1) % 1000 == 0)
            printf("\rcounting triplet site patterns %.0f%% done..", (locus+1.0)/nloci*100);
      }
      if (nloci_used > 1)
         for (i = 0; i < 5; i++)
            sf[i] = sqrt(sf[i] / (nloci_used - 1));
   }

   printf("\nnloci used: %d.. %s", nloci_used, timer(timestr));
   printf("\nobs-m: "); for (i = 0; i < 5; i++) printf("%12.8f", mf[i]);
   printf("\nobs-s: "); for (i = 0; i < 5; i++) printf("%12.8f", sf[i]);
   printf("\nobs-v: "); for (i = 0; i < 5; i++) printf("%12.8f", sf[i]* sf[i]);
   printf("\n1 = %.6f\n", sum(mf, 5));
   printf("\n");

   double nsites = 0;
   for (locus = 0; locus < nloci; locus++)
     nsites += (nsites_used[locus] - nsites) / (locus + 1.0);

   if (0)
     evaluate_models(mf, sf, nloci, (int)(nsites + 0.5), transform_triplet);
   else
     test_mean_var(mf, sf, nloci, (int)(nsites + 0.5), transform_triplet);

   for (i = 0; i < nloci; ++i)
      msa_destroy(msa_list[i]);
   free(msa_list);
   list_clear(maplist, map_dealloc);
   free(maplist);
   free(ipop);

   free(missing_species_sites);
   free(missing_species_loci);
   free(nsites_used);
   free(species_site_vec);
   for (k = 0; k < nspecies; ++k) free(species[k]);
   free(species);
}
