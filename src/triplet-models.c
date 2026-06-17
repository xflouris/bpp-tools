/* Functions here are from speciestree3s.c.
*  Early version is called speciestree3s-ZY2021.c
*/

#include "bpp-tools.h"
#define square(a) ((a)*(a))
#define rndexp(mean) (-(mean)*log(rndu(0)))

FILE* frub;
static double _mf[5], _sf[5];
static int _model = 0;
enum Models { MSC = 0, SISTER = 1, INFLOW = 2, OUTFLOW = 3, GHOST = 4, GHOSTC = 5 };
static char* modelstr[6] = { "MSC", "SISTER", "INFLOW", "OUTFLOW", "GHOST", "GHOSTC" };
static int _nloci, _length;
static int _transform_triplet;
/* 1 = logit: y = log(x/(1-x)); 2 = log: x = log(f) for x1 x2 x3 */
static char* transformstr[3] = { "no transform", "logit transform y = log x/(1-x)", "transform y = log(x123)"};

int ZhuYang2023(void);
void FromParaToRatiopara(double ratiopara[5], double tau0, double tau11s, double tau1s, double theta, double phi);
void p0124Fromb0b1(double p[5], double b[2]);
void generate_abd(double a[6], double para[3]);
void generate_p01234(double p[5], double para[3]);
void generate_v01234(double v[5], double para[3]);
void generate_p01234_msci(double p[5], double para[5]);
void generate_v01234_msci(double v[5], double para[5]);
double loglike_triplet(double ratiopara[5], int dimt);

#if(0)
int ZhuYang2023(void)
{
   char debug = 0;
   int nr = 1e0, ir, i, n[5], Gtree, GtreeE, nGtree[4];
   int Lalias[5], n2step = 0, nconcat = 0;         /* 2step2 and concat2 ignore ties */
   double mp01234[5] = { 0 }, p01234[5], b[2], Falias[5];
   double PG[4] = { 0 }, PGe[4] = { 0 }, PSGe = 0, Ptie = 0, z;
   double tau0 = 0.05, tau1 = 0.04, tau1t = 0.03, tau1s = 0.04, theta0 = 0.1, theta1 = 0.01, phi0 = 0.3;
   double mf0[5], sf0[5];

   //printf("tau0 tau1 theta0 theta1 = %.3f %.3f %.3f %.3f\n", tau0, tau1, theta0, theta1);
   printf("_model = MSC, INFLOW, OUTFLOW, GHOST, SISTER, GHOSTC\n");
   printf("nloci = %5d length = %5d\n\n", _nloci, _length);
   rnd_init(0);

   FILE* fout = fopen("output.txt", "w");

   for (ir = 0; ir < nr; ir++) {
      if (debug) { //Data Generation
         nGtree[0] = nGtree[1] = nGtree[2] = nGtree[3] = 0;
         for (i = 0; i < 5; i++) mf0[i] = sf0[i] = 0;
         for (locus = 0; locus < _nloci; locus++) {
            b[1] = rndexp(theta1 / 2);
            if (b[1] < tau0 - tau1) {
               Gtree = 1;
               b[1] += tau1;
               b[0] = tau0 - b[1] + rndexp(theta0 / 2);
            }
            else {
               Gtree = (int)(rndu(0) * 3) + 1;
               b[1] = tau0 + rndexp(theta0 / 6);
               b[0] = rndexp(theta0 / 2);
            }
            PG[Gtree]++;
            p0124Fromb0b1(p01234, b);
            if (Gtree == 2) { z = p01234[1]; p01234[1] = p01234[2]; p01234[2] = z; }
            else if (Gtree == 3) { z = p01234[1]; p01234[1] = p01234[3]; p01234[3] = z; }
            for (i = 0; i < 5; i++)  mp01234[i] += p01234[i];
            if (_length > 1) {
               MultiNomialAliasSetTable(5, p01234, Falias, Lalias);
               MultiNomialAlias(0, _length, 5, Falias, Lalias, n);
            }
            else {
               double r = rndu(0), s = p01234[0];
               n[0] = n[1] = n[2] = n[3] = n[4] = 0;
               for (k = 0; k < 5 - 1; k++)
                  if (r < s)  break;
                  else        s += p01234[k + 1];
               n[k] = 1;
            }
            for (i = 0; i < 5; i++) {
               double z = (double)n[i] / _length;
               if (n[i] == 0) z = 0.5 / _length;
               else if (n[i] == _length) z = 1 - 0.5 / _length;
               z = log(z / (1 - z));
               mf0[i] += z;
               sf0[i] += z * z;
            }
         }  /* for (locus) */
         for (i = 0; i < 5; i++) {
            _mf[i] = mf0[i] /= _nloci;
            _sf[i] = sf0[i] = sqrt(sf0[i] / _nloci - mf0[i] * mf0[i]);
         }
      }
      if (1) {
         //printf("\n%6d\n", ir);
         printf("f:");  for (i = 0; i < 4; i++) printf("%8.5f ", _mf[i]);  printf("\n");
         printf("s:");  for (i = 0; i < 4; i++) printf("%8.5f ", _sf[i]);  printf("\n");
      }
      evaluate_models(_mf, _sf, _nloci, _length, _transform_triplet);
   }
   return (0);
}
#endif

void FromParaToRatiopara(double ratiopara[5], double tau0, double tau1t, double tau1s, double theta, double phi)
{
/* first para is largest two, the other two tau's are expressed as ratios
*/
   if (_model == GHOST || _model == GHOSTC) { /* tau_o, tau0/tau_o, tau1s/tau0 */
      ratiopara[0] = tau1t;          /* tau_o */
      ratiopara[1] = tau0 / tau1t;
      ratiopara[2] = tau1s / tau0;
      ratiopara[3] = theta;
      ratiopara[4] = phi;
   }
   else {
      ratiopara[0] = tau0;
      ratiopara[1] = tau1t / tau1s;
      ratiopara[2] = tau1s / tau0;
      ratiopara[3] = theta;
      ratiopara[4] = phi;
   }
   return;
}

void p0124Fromb0b1(double p[5], double b[2])
{
/* This calculates p0,p1,p2,p4 for one locus, given b0 and b1.
   b0 is the internal branch length, and b1 is the terminal branch length,
   with b0 + b1 to be the age of the root.
*/
   double e1, e2, e3;
   if (b[0] < 0 || b[1] < 0) fatal("b0 b1 < 0");
   e1 = exp(-4. / 3 * b[1]);
   e2 = exp(-8. / 3 * (b[0] + b[1]));
   e3 = e1 * e2;
   e1 = e1 * e1;
   p[0] = (1. + 3 * e1 + 6 * e2 + 6 * e3) / 16;
   p[1] = (3. + 9 * e1 - 6 * e2 - 6 * e3) / 16;
   p[2] = p[3] = (3. - 3 * e1 + 6 * e2 - 6 * e3) / 16;
   p[4] = (6. - 6 * e1 - 12 * e2 + 12 * e3) / 16;
}

void generate_abd(double a[6], double para[3])
{
// from para: tau0, tau1, theta to a: a0, a1, b0, b1, d0
   double da, db, dd;//d? stands for denominator of ? 
   da = 3 + 4 * para[2];
   db = 3 + 2 * para[2];
   dd = 3 + 8 * para[2];
   double e0, e1, e02, e12, e04;//e1,2 stands for exp(-4 tau1,2/ 3)
   e0 = exp(-4 * para[0] / 3);
   e1 = exp(-4 * para[1] / 3);
   e02 = e0 * e0;
   e12 = e1 * e1;
   e04 = e02 * e02;
   a[0] = e02 / da; //a0
   a[1] = e12 / da; //a1
   a[2] = e0 / db; //b0
   a[3] = e1 / db; //b1
   a[4] = e04 / dd; //d0
   a[5] = exp(-2 * (para[0] - para[1]) / para[2]);
}

void generate_p01234(double p[5], double para[3])
{
   //from para: tau0, tau1, theta to p: p0, p1, p2, p4
   double r[6], a0, a1, a0b1;
   generate_abd(r, para);
   a0 = r[0];
   a1 = r[1];
   a0b1 = r[0] * r[3];
   p[0] = 1.0 / 16 * (1 + 18 * a0 + 54 * a0b1 + 9 * a1);
   p[1] = 3.0 / 16 * (1 - 6 * a0 - 18 * a0b1 + 9 * a1);
   p[2] = 3.0 / 16 * (1 + 6 * a0 - 18 * a0b1 - 3 * a1);
   p[3] = p[2];
   p[4] = 6.0 / 16 * (1 - 6 * a0 + 18 * a0b1 - 3 * a1);
}

void generate_v01234(double v[5], double para[3])
{
// from para: tau0, tau1, theta to v: v0, v1, v2, v4
   double a[6]; //a: a0, a1, b0, b1, d0
   generate_abd(a, para);
   double theta2, a02, a12, a02b1, a02b12, a0a1b1, b0d0, d1, d2, d3, d42, d4;

   theta2 = para[2] * para[2];
   a02 = a[0] * a[0];
   a12 = a[1] * a[1];
   a02b1 = a02 * a[3];
   a02b12 = a02b1 * a[3];
   a0a1b1 = a[0] * a[1] * a[3];
   b0d0 = a[2] * a[4];
   d1 = 3 + 8.0 * para[2];
   d2 = 1 + 2.0 * para[2];
   d3 = (9 + 8.0 * para[2]) * d1;
   d42 = (1 + 2.0 * para[2]) * (9 + 10.0 * para[2]);
   d4 = d42 * (3 + 4.0 * para[2]);
   v[0] = 27.0 / 16 * theta2 * (4.0 / d1 * a02 + 1.0 / d1 * a12 + 24.0 / d1 * a02b1
                                + 3.0 * (15.0 + 4 * para[2]) / d1 * a02b12 + 2.0 / d2 * a0a1b1
                                - 16.0 * para[2] * a[5] / d3 * a02 - 16.0 * para[2] * a[5] / d4 * b0d0);
   v[1] = 27.0 / 16 * theta2 * (4.0 / d1 * a02 + 9.0 / d1 * a12 + 24.0 / d1 * a02b1
                                + 3.0 * (15.0 + 4 * para[2]) / d1 * a02b12 - 6.0 / d2 * a0a1b1 -
                                16.0 * (3.0 + para[2]) * a[5] / d3 * a02 - 48.0 * (1.0 + para[2]) * a[5] / d4 * b0d0);
   v[2] = 27.0 / 16 * theta2 * (4.0 / d1 * a02 + 1.0 / d1 * a12 - 24.0 / d1 * a02b1
                                + 3.0 * (15.0 + 4 * para[2]) / d1 * a02b12 + 2.0 / d2 * a0a1b1 +
                                24.0 * (1.0 + 2.0 * para[2]) * a[5] / d3 * a02 + 8.0 * a[5] / d42 * b0d0);
   v[3] = v[2];
   v[4] = 27.0 / 4 * theta2 * (4.0 / d1 * a02 + 1.0 / d1 * a12 - 24.0 / d1 * a02b1
                               + 3.0 * (15.0 + 4 * para[2]) / d1 * a02b12 - 2.0 / d2 * a0a1b1
                               - 16.0 * para[2] * a[5] / d3 * a02 + 16.0 * para[2] * a[5] / d4 * b0d0);
}

void generate_p01234_msci(double p[5], double para[5])
{
// from para: tau0,tau1t, tau1s, theta phi to pa: p0a, p1a, p2a, p3a, p4a; Put tau_o in tau1t for _model GHOST
   double p1[5], p2[5], t0[3], phi = para[4];
   if (_model == MSC)	phi = 0;
   t0[0] = para[0];
   t0[1] = para[2];
   t0[2] = para[3];
   generate_p01234(p1, t0);
   if (_model == INFLOW || _model == SISTER || _model == MSC) {  /* inflow */
      t0[0] = para[0];
      t0[1] = para[1];
      t0[2] = para[3];
   }
   else if (_model == GHOST || _model == GHOSTC)
   {
      t0[0] = para[1];
      t0[1] = para[0];
      t0[2] = para[3];
   }
   else {    /* outflow */
      t0[0] = para[2];
      t0[1] = para[1];
      t0[2] = para[3];
   }
   generate_p01234(p2, t0);
   p[0] = (1 - phi) * p1[0] + phi * p2[0];
   p[1] = (1 - phi) * p1[1] + phi * p2[2];
   p[2] = (1 - phi) * p1[2] + phi * p2[1];
   p[3] = (1 - phi) * p1[2] + phi * p2[2];
   p[4] = (1 - phi) * p1[4] + phi * p2[4];

   if (_model == SISTER || _model == GHOSTC)
   {
      p[1] = (1 - phi) * p1[1] + phi * p2[1];
      p[2] = (1 - phi) * p1[2] + phi * p2[2];
   }
}

void generate_v01234_msci(double v[5], double para[5])
{
// from para: tau0, tau1t, tau1s, theta, phi to va: v0a, v1a, v2a, v3a, v4a; Put tau_o in tau1t for _model GHOST
   double p1[5], p2[5], v1[5], v2[5], para0[3], phi = para[4];
   if (_model == MSC)	phi = 0;
   para0[0] = para[0];
   para0[1] = para[2];
   para0[2] = para[3];
   generate_p01234(p1, para0);
   generate_v01234(v1, para0);
   if (_model == INFLOW || _model == SISTER || _model == MSC) { //tau0, tau1t
      para0[0] = para[0];
      para0[1] = para[1];
      para0[2] = para[3];
   }
   else if (_model == GHOST || _model == GHOSTC) { //tau_o, tau0
      para0[0] = para[1];
      para0[1] = para[0];
      para0[2] = para[3];
   }
   else { //tau1s, tau1t
      para0[0] = para[2];
      para0[1] = para[1];
      para0[2] = para[3];
   }
   generate_p01234(p2, para0);
   generate_v01234(v2, para0);
   v[0] = (1 - phi) * v1[0] + phi * v2[0] + phi * (1 - phi) * square(p1[0] - p2[0]);
   v[1] = (1 - phi) * v1[1] + phi * v2[2] + phi * (1 - phi) * square(p1[1] - p2[2]);
   v[2] = (1 - phi) * v1[2] + phi * v2[1] + phi * (1 - phi) * square(p1[2] - p2[1]);
   v[3] = (1 - phi) * v1[2] + phi * v2[2] + phi * (1 - phi) * square(p1[2] - p2[2]);
   v[4] = (1 - phi) * v1[4] + phi * v2[4] + phi * (1 - phi) * square(p1[4] - p2[4]);
   if (_model == SISTER || _model == GHOSTC) {
      v[1] = (1 - phi) * v1[1] + phi * v2[1] + phi * (1 - phi) * square(p1[1] - p2[1]);
      v[2] = (1 - phi) * v1[2] + phi * v2[2] + phi * (1 - phi) * square(p1[2] - p2[2]);
   }
}


int evaluate_models(double mf0[], double sf0[], int nloci0, int length0, int transform_triplet)
{
   int i, best_model = 0;
   double tau0 = 0.02, tau1s = 0.01, tau1t = 0.01, theta0 = 0.025, phi0 = 0.1;
   int dim = 5;
   double lnL[6], lnLt, e = 1e-9, ratiopara[5];
   double Rangepara[5][2] = { {0.0001,0.5}, {0.0001,0.999}, {0.0001,0.999}, {0.0001,0.2}, {0.001,0.999} }; // tau0, ratio1, ratio2, theta, phi
   double space0[1000];

   opt_debug = 1;
   if (opt_debug) frub = xopen("rub.txt", "w");

   _nloci = nloci0;  _length = length0;  _transform_triplet = transform_triplet;
   for (i = 0; i < 5; i++)	_mf[i] = mf0[i];
   for (i = 0; i < 5; i++)	_sf[i] = sf0[i];
   for (i = 0; i < 5; i++)	_sf[i] = square(_sf[i]);

   printf("\nmodel = 0 MSC, 1 SISTER, 2 INFLOW, 3 OUTFLOW, 4 GHOST, 5 GHOSTC\n");
   if(_transform_triplet) printf("\n%s used\n", transformstr[_transform_triplet]);

   for ( ; _model <= 5; _model++) {
      if (opt_debug) {
         fprintf(frub, "\n\n%-10s", modelstr[_model]);
         fprintf(frub, "\nnloci = %d  length = %d", _nloci, _length);
         fprintf(frub, "\nm: "); for (i = 0; i < 5; i++) fprintf(frub, "%12.8f", _mf[i]);
         fprintf(frub, "\nv: "); for (i = 0; i < 5; i++) fprintf(frub, "%12.8f", _sf[i]);
         fprintf(frub, "\n");
      }
      FromParaToRatiopara(ratiopara, tau0, tau1t, tau1s, theta0, phi0);
      ming2((opt_debug ? frub : NULL), &lnLt, loglike_triplet, NULL, ratiopara, Rangepara, space0, e, dim);
      printf("\n%-10s lnL = %12.3f ", modelstr[_model], -lnLt);
      if (_model == MSC)
         printf("tau0 = %8.4f tau1s = %8.4f theta = %8.4f", ratiopara[0], ratiopara[0] * ratiopara[2], ratiopara[3]);
      else if (_model == GHOST || _model == GHOSTC)
         printf("tauR = %8.4f tau1s = %8.4f tau_o = %8.4f theta = %8.4f phi = %8.4f", ratiopara[0] * ratiopara[2], ratiopara[0] * ratiopara[1] * ratiopara[2], ratiopara[0], ratiopara[3], ratiopara[4]);
      else
         printf("tau0 = %8.4f tau1t = %8.4f tau1s = %8.4f theta = %8.4f phi = %8.4f", ratiopara[0], ratiopara[0] * ratiopara[1] * ratiopara[2], ratiopara[0] * ratiopara[2], ratiopara[3], ratiopara[4]);

      lnL[_model] = -lnLt;
      /* model of gene flow is selected only if it is siginifantly better than M0. */
      if (_model && lnL[_model] > lnL[0] + 5.99 / 2) {
         printf(" **");
         if (best_model == 0 || lnL[best_model] < lnL[_model])
            best_model = _model;
      }
   }
   printf("\nSelected model: %s\n", modelstr[best_model]);

   if (opt_debug) fclose(frub);

   return (best_model);
}

int test_mean_var(double mf0[], double sf0[], int nloci0, int length0, int transform_triplet)
{
  int i;
  double tau0 = 0.02, tau1s = 0.01, tau1t = 0.01, theta0 = 0.025, phi0 = 0.1234567;
  double lnLt, ratiopara[5];
 
  opt_debug = 1;
  _nloci = nloci0;  _length = length0;  _transform_triplet = transform_triplet;
  assert(transform_triplet == 0);

  for (i = 0; i < 5; i++)	_mf[i] = mf0[i];
  for (i = 0; i < 5; i++)	_sf[i] = sf0[i];
  for (i = 0; i < 5; i++)	_sf[i] = square(_sf[i]);

  printf("\ninput model ID: 0: MSC, 1: SISTER, 2: INFLOW, 3: OUTFLOW, 4: GHOST, 5: GHOSTC? ");
  scanf("%d", &_model);

  printf("\ninput tau0, tau1t, tau1s, theta0, phi0?\n");
  scanf("%lf%lf%lf%lf%lf", &tau0, &tau1t, &tau1s, &theta0, &phi0);
  printf("\ntau0, tau1t, tau1s, theta, phi: %9.5f %9.5f %9.5f %9.5f %9.5f\n", 
         tau0, tau1t, tau1s, theta0, phi0);

  FromParaToRatiopara(ratiopara, tau0, tau1t, tau1s, theta0, phi0);
  lnLt = loglike_triplet(ratiopara, 5);
  exit(0);
}


double loglike_triplet(double ratiopara[5], int dimt)
{
   /* Ratiopara: tau0, tau1t/tau1s, tau1s/tau0, theta, phi, dimt is prepared for ming2
      model GHOST: tau_o, tau0/tau_o, tau1s/tau0
      Look at comments inside FromParaToRatiopara().
   */
   double lnL = 0, p[5], v[5], va[5], para[5];
   int m = _nloci, n = _length, j;

   if (_model == GHOST || _model == GHOSTC) {
      para[1] = ratiopara[0];            /* tau_o for ghost */
      para[0] = para[1] * ratiopara[1];  /* tau_r = tau_ABC */
      para[2] = para[0] * ratiopara[2];  /* tau_s = tau_AB */
      para[3] = ratiopara[3];            /* theta */
      para[4] = ratiopara[4];            /* phi */
      generate_p01234_msci(p, para);
      generate_v01234_msci(va, para);
   }
   else {
      para[0] = ratiopara[0];            /* tau_r = tau_ABC */
      para[2] = para[0] * ratiopara[2];  /* tau_s = tau_AB */
      para[1] = para[2] * ratiopara[1];  /* tau_h */
      para[3] = ratiopara[3];            /* theta */
      para[4] = ratiopara[4];            /* phi */
      generate_p01234_msci(p, para);
      generate_v01234_msci(va, para);
   }

   for (j = 0; j < 5; j++)
      v[j] = p[j] * (1.0 - p[j]) / n + va[j] * (n - 1.0) / n;

   if (_transform_triplet == 1)        /* logit: y = log(x/(1-x)) */
      for (j = 0; j <= 4; j++) {
         double ux = p[j], vx = v[j];
         double uu1uu = ux * ux * (1 - ux) * (1 - ux), u12u = (1 - 2 * ux) / uu1uu;
         p[j] = log(ux / (1 - ux)) - 0.5 * u12u * vx;
         v[j] = 1 / uu1uu * vx;
         if (1)
            v[j] += 0.5 * u12u * u12u * vx * vx;
      }
   else if (_transform_triplet == 2) { /* log: x = log(f) for x1 x2 x3*/
      p[0] = 0; v[0] = 1;
      for (j = 0; j <= 4; j++) {
         double ux = p[j], vx = v[j], vuu= vx/(ux*ux);
         p[j] = log(ux) - vuu/2;
         v[j] = vuu * (1 + vuu * vuu / 2);
      }
   }

   if (opt_debug) {
      printf("\nexp-p: "); for (j = 0; j < 5; j++) printf("%12.8f", p[j]);
      printf("\nexp-s: "); for (j = 0; j < 5; j++) printf("%12.8f", sqrt(v[j]));
      printf("\nexp-v: "); for (j = 0; j < 5; j++) printf("%12.8f", v[j]); printf("\n");
      printf("\n1 = %.6f\n", sum(p, 5));
   }
   if (_transform_triplet == 0) {
     double lambda = 5;
     for (j = 1; j <= 3; j++)
       lnL += square(_mf[j] - p[j]) / (p[j] * (1 - p[j])) * lambda
           + square(_sf[j]/v[j] - 1);
   }
   else {
     for (j = 1; j <= 3; j++)
       lnL += log(v[j]) + (_sf[j] + square(_mf[j] - p[j])) / v[j];
     lnL *= m / 2.0;
   }
   return (lnL);
}
