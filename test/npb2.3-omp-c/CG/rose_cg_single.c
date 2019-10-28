/*--------------------------------------------------------------------
  
  NAS Parallel Benchmarks 2.3 OpenMP C versions - CG
  This benchmark is an OpenMP C version of the NPB CG code.
  
  The OpenMP C versions are developed by RWCP and derived from the serial
  Fortran versions in "NPB 2.3-serial" developed by NAS.
  Permission to use, copy, distribute and modify this software for any
  purpose with or without fee is hereby granted.
  This software is provided "as is" without express or implied warranty.
  
  Send comments on the OpenMP C versions to pdp-openmp@rwcp.or.jp
  Information on OpenMP activities at RWCP is available at:
           http://pdplab.trc.rwcp.or.jp/pdperf/Omni/
  
  Information on NAS Parallel Benchmarks 2.3 is available at:
  
           http://www.nas.nasa.gov/NAS/NPB/
--------------------------------------------------------------------*/
/*--------------------------------------------------------------------
  Authors: M. Yarrow
           C. Kuszmaul
  OpenMP C version: S. Satoh
  
--------------------------------------------------------------------*/
/*
c---------------------------------------------------------------------
c  Note: please observe that in the routine conj_grad three 
c  implementations of the sparse matrix-vector multiply have
c  been supplied.  The default matrix-vector multiply is not
c  loop unrolled.  The alternate implementations are unrolled
c  to a depth of 2 and unrolled to a depth of 8.  Please
c  experiment with these to find the fastest for your particular
c  architecture.  If reporting timing results, any of these three may
c  be used without penalty.
c---------------------------------------------------------------------
*/
//#include "npb-C.h"
/*
  NAS Parallel Benchmarks 2.3 OpenMP C Versions
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(_OPENMP)
#include <omp.h>
#endif /* _OPENMP */
typedef int boolean;
typedef struct {
// 114 lv-analysis-out: bot
double real;
// 114 lv-analysis-in : bot
// 115 lv-analysis-out: bot
double imag;
// 115 lv-analysis-in : bot
}dcomplex;
#define TRUE    1
#define FALSE   0
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define pow2(a) ((a)*(a))
#define get_real(c) c.real
#define get_imag(c) c.imag
#define cadd(c,a,b) (c.real = a.real + b.real, c.imag = a.imag + b.imag)
#define csub(c,a,b) (c.real = a.real - b.real, c.imag = a.imag - b.imag)
#define cmul(c,a,b) (c.real = a.real * b.real - a.imag * b.imag, \
                     c.imag = a.real * b.imag + a.imag * b.real)
#define crmul(c,a,b) (c.real = a.real * b, c.imag = a.imag * b)
double randlc(double *x,double a);
extern void vranlc(int ,double *,double ,double *);
void timer_clear(int n);
void timer_start(int n);
void timer_stop(int n);
double timer_read(int n);
void c_print_results(char *name,char cclass,int n1,int n2,int n3,int niter,int nthreads,double t,double mops,char *optype,int passed_verification,char *npbversion,char *compiletime,char *cc,char *clink,char *c_lib,char *c_inc,char *cflags,char *clinkflags,char *rand);
//#include "npbparams.h"
/******************/
/* default values */
/******************/
#ifndef CLASS
#define CLASS 'B'
#endif
#if CLASS == 'S'
/* CLASS = S */
/*
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the classc of the NPB
c  in this directory. Do not modify it by hand.
*/
#define NA      1400
#define NONZER  7
#define NITER   15
#define SHIFT   10.0
#define RCOND   1.0e-1
#define CONVERTDOUBLE   FALSE
#endif
#if CLASS == 'W'
/* CLASS = W */
/*
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the classc of the NPB
c  in this directory. Do not modify it by hand.
*/
#define NA      7000
#define NONZER  8
#define NITER   15
#define SHIFT   12.0
#define RCOND   1.0e-1
#define CONVERTDOUBLE   FALSE
#endif
#if CLASS == 'A'
/* CLASS = A */
/*
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the classc of the NPB
c  in this directory. Do not modify it by hand.
*/
#define NA      14000
#define NONZER  11
#define NITER   15
#define SHIFT   20.0
#define RCOND   1.0e-1
#define CONVERTDOUBLE   FALSE
#endif
#if CLASS == 'B'
/* CLASS = B */
/*
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the classc of the NPB
c  in this directory. Do not modify it by hand.
*/
#define NA      75000
#define NONZER  13
#define NITER   75
#define SHIFT   60.0
#define RCOND   1.0e-1
#define CONVERTDOUBLE   FALSE
#endif
#if CLASS == 'C'
/* CLASS = C */
/*
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the classc of the NPB
c  in this directory. Do not modify it by hand.
*/
#define NA      150000
#define NONZER  15
#define NITER   75
#define SHIFT   110.0
#define RCOND   1.0e-1
#define CONVERTDOUBLE   FALSE
#endif
#define COMPILETIME "28 Oct 2014"
#define NPBVERSION "2.3"
#define CS1 "gcc"
#define CS2 "$(CC)"
#define CS3 "(none)"
#define CS4 "-I../common"
#define CS5 "-fopenmp -O2"
#define CS6 "-lm -fopenmp"
#define CS7 "randdp"
#define	NZ	NA*(NONZER+1)*(NONZER+1)+NA*(NONZER+2)
/* global variables */
/* common /partit_size/ */
// 116 lv-analysis-out: bot
static int naa;
// 116 lv-analysis-in : bot
// 117 lv-analysis-out: bot
static int nzz;
// 117 lv-analysis-in : bot
// 118 lv-analysis-out: bot
static int firstrow;
// 118 lv-analysis-in : bot
// 119 lv-analysis-out: bot
static int lastrow;
// 119 lv-analysis-in : bot
// 120 lv-analysis-out: bot
static int firstcol;
// 120 lv-analysis-in : bot
// 121 lv-analysis-out: bot
static int lastcol;
// 121 lv-analysis-in : bot
/* common /main_int_mem/ */
/* colidx[1:NZ] */
// 122 lv-analysis-out: bot
static int colidx[75000 * (13 + 1) * (13 + 1) + 75000 * (13 + 2) + 1];
// 122 lv-analysis-in : bot
/* rowstr[1:NA+1] */
// 123 lv-analysis-out: bot
static int rowstr[75000 + 1 + 1];
// 123 lv-analysis-in : bot
/* iv[1:2*NA+1] */
// 124 lv-analysis-out: bot
static int iv[2 * 75000 + 1 + 1];
// 124 lv-analysis-in : bot
/* arow[1:NZ] */
// 125 lv-analysis-out: bot
static int arow[75000 * (13 + 1) * (13 + 1) + 75000 * (13 + 2) + 1];
// 125 lv-analysis-in : bot
/* acol[1:NZ] */
// 126 lv-analysis-out: bot
static int acol[75000 * (13 + 1) * (13 + 1) + 75000 * (13 + 2) + 1];
// 126 lv-analysis-in : bot
/* common /main_flt_mem/ */
/* v[1:NA+1] */
// 127 lv-analysis-out: bot
static double v[75000 + 1 + 1];
// 127 lv-analysis-in : bot
/* aelt[1:NZ] */
// 128 lv-analysis-out: bot
static double aelt[75000 * (13 + 1) * (13 + 1) + 75000 * (13 + 2) + 1];
// 128 lv-analysis-in : bot
/* a[1:NZ] */
// 129 lv-analysis-out: bot
static double a[75000 * (13 + 1) * (13 + 1) + 75000 * (13 + 2) + 1];
// 129 lv-analysis-in : bot
/* x[1:NA+2] */
// 130 lv-analysis-out: bot
static double x[75000 + 2 + 1];
// 130 lv-analysis-in : bot
/* z[1:NA+2] */
// 131 lv-analysis-out: bot
static double z[75000 + 2 + 1];
// 131 lv-analysis-in : bot
/* p[1:NA+2] */
// 132 lv-analysis-out: bot
static double p[75000 + 2 + 1];
// 132 lv-analysis-in : bot
/* q[1:NA+2] */
// 133 lv-analysis-out: bot
static double q[75000 + 2 + 1];
// 133 lv-analysis-in : bot
/* r[1:NA+2] */
// 134 lv-analysis-out: bot
static double r[75000 + 2 + 1];
// 134 lv-analysis-in : bot
/* w[1:NA+2] */
// 135 lv-analysis-out: bot
static double w[75000 + 2 + 1];
// 135 lv-analysis-in : bot
/* common /urando/ */
// 136 lv-analysis-out: bot
static double amult;
// 136 lv-analysis-in : bot
// 137 lv-analysis-out: bot
static double tran;
// 137 lv-analysis-in : bot
/* function declarations */
static void conj_grad(int colidx[],int rowstr[],double x[],double z[],double a[],double p[],double q[],double r[],double w[],double *rnorm);
static void makea(int n,int nz,double a[],int colidx[],int rowstr[],int nonzer,int firstrow,int lastrow,int firstcol,int lastcol,double rcond,int arow[],int acol[],double aelt[],double v[],int iv[],double shift);
static void sparse(double a[],int colidx[],int rowstr[],int n,int arow[],int acol[],double aelt[],int firstrow,int lastrow,double x[],boolean mark[],int nzloc[],int nnza);
static void sprnvc(int n,int nz,double v[],int iv[],int nzloc[],int mark[]);
static int icnvrt(double x,int ipwr2);
static void vecset(int n,double v[],int iv[],int *nzv,int i,double val);
/*--------------------------------------------------------------------
      program cg
--------------------------------------------------------------------*/

int main(int argc,char **argv)
// 138 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_326,tmp_327,tmp_346}
{
// 141 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  int i;
// 141 lv-analysis-in : {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 142 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  int j;
// 142 lv-analysis-in : {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 143 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  int k;
// 143 lv-analysis-in : {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 144 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  int it;
// 144 lv-analysis-in : {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 145 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  int nthreads = 1;
// 145 lv-analysis-in : {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 146 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  double zeta;
// 146 lv-analysis-in : {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 147 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  double rnorm;
// 147 lv-analysis-in : {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 148 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  double norm_temp11;
// 148 lv-analysis-in : {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 149 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  double norm_temp12;
// 149 lv-analysis-in : {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 150 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  double t;
// 150 lv-analysis-in : {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 151 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  double mflops;
// 151 lv-analysis-in : {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 152 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  char cclass;
// 152 lv-analysis-in : {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 153 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  boolean verified;
// 153 lv-analysis-in : {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 154 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  double zeta_verify_value;
// 154 lv-analysis-in : {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 155 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  double epsilon;
// 155 lv-analysis-in : {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 156 lv-analysis-out: {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  firstrow = 1;
// 156 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 157 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  lastrow = 75000;
// 157 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 158 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  firstcol = 1;
// 158 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 159 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  lastcol = 75000;
// 159 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 160 lv-analysis-out: bot
  if (
// 161 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
75000 == 1400 && 13 == 7 && 75 == 15 && 60.0 == 10.0
// 161 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
) {
// 163 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
    cclass = 'S';
// 163 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 164 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
    zeta_verify_value = 8.5971775078648;
// 164 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  }
   else {
// 166 lv-analysis-out: bot
    if (
// 167 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
75000 == 7000 && 13 == 8 && 75 == 15 && 60.0 == 12.0
// 167 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
) {
// 169 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
      cclass = 'W';
// 169 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 170 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
      zeta_verify_value = 10.362595087124;
// 170 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
    }
     else {
// 172 lv-analysis-out: bot
      if (
// 173 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
75000 == 14000 && 13 == 11 && 75 == 15 && 60.0 == 20.0
// 173 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
) {
// 175 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
        cclass = 'A';
// 175 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 176 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
        zeta_verify_value = 17.130235054029;
// 176 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
      }
       else {
// 178 lv-analysis-out: bot
        if (
// 179 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
75000 == 75000 && 13 == 13 && 75 == 75 && 60.0 == 60.0
// 179 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
) {
// 181 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
          cclass = 'B';
// 181 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 182 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
          zeta_verify_value = 22.712745482631;
// 182 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
        }
         else {
// 184 lv-analysis-out: bot
          if (
// 185 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
75000 == 150000 && 13 == 15 && 75 == 75 && 60.0 == 110.0
// 185 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
) {
// 187 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
            cclass = 'C';
// 187 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 188 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
            zeta_verify_value = 28.973605592845;
// 188 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
          }
           else {
// 190 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
            cclass = 'U';
// 190 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
          }
// 184 lv-analysis-in : bot
        }
// 178 lv-analysis-in : bot
      }
// 172 lv-analysis-in : bot
    }
// 166 lv-analysis-in : bot
  }
// 160 lv-analysis-in : bot
// 191 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  printf("\n\n NAS Parallel Benchmarks 2.3 OpenMP C version - CG Benchmark\n");
// 191 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 193 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  printf(" Size: %10d\n",75000);
// 193 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 195 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
  printf(" Iterations: %5d\n",75);
// 195 lv-analysis-in : {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
// 197 lv-analysis-out: {tv_sec_72,tv_usec_73,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
  naa = 75000;
// 197 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
// 198 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
  nzz = 75000 * (13 + 1) * (13 + 1) + 75000 * (13 + 2);
// 198 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
/*--------------------------------------------------------------------
c  Initialize random number generator
c-------------------------------------------------------------------*/
// 199 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
  tran = 314159265.0;
// 199 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
// 200 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
  amult = 1220703125.0;
// 200 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
// 201 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
  zeta = randlc(&tran,amult);
// 201 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_326,tmp_327}
/*--------------------------------------------------------------------
c  
c-------------------------------------------------------------------*/
// 203 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_346}
  makea(naa,nzz,a,colidx,rowstr,13,firstrow,lastrow,firstcol,lastcol,1.0e-1,arow,acol,aelt,v,iv,60.0);
// 203 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_326,tmp_327,tmp_328,tmp_329,tmp_330,tmp_331,tmp_332,tmp_333,tmp_334,tmp_335,tmp_336,tmp_337,tmp_338,tmp_339,tmp_340,tmp_341,tmp_342,tmp_346}
/*---------------------------------------------------------------------
c  Note: as a result of the above call to makea:
c        values of j used in indexing rowstr go from 1 --> lastrow-firstrow+1
c        values of colidx which are col indexes go from firstcol --> lastcol
c        So:
c        Shift the col index vals from actual (firstcol --> lastcol ) 
c        to local, i.e., (1 --> lastcol-firstcol+1)
c---------------------------------------------------------------------*/
// 205 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
  
#pragma omp parallel private(it,i,j,k)
// 205 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
{
// 207 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
    
#pragma omp for nowait
// 207 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
// 208 lv-analysis-out: bot
    for (
// 209 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
j = 1
// 209 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
; 
// 210 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
j <= lastrow - firstrow + 1;
// 210 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
 j++) {
// 213 lv-analysis-out: bot
      for (
// 214 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
k = rowstr[j]
// 214 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,k_170,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
; 
// 215 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,k_170,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
k < rowstr[j + 1];
// 215 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,k_170,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
 k++) {
// 218 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,k_170,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
        colidx[k] = colidx[k] - firstcol + 1;
// 218 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,k_170,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
      }
// 213 lv-analysis-in : bot
    }
// 208 lv-analysis-in : bot
/*--------------------------------------------------------------------
c  set starting vector to (1, 1, .... 1)
c-------------------------------------------------------------------*/
// 219 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
    
#pragma omp for nowait
// 219 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
// 220 lv-analysis-out: bot
    for (
// 221 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
i = 1
// 221 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,i_168,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
; 
// 222 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,i_168,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
i <= 75000 + 1;
// 222 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,i_168,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
 i++) {
// 225 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,z_161,p_162,q_163,r_164,w_165,i_168,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
      x[i] = 1.0;
// 225 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,i_168,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
    }
// 220 lv-analysis-in : bot
// 226 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
    
#pragma omp single
// 226 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
// 227 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
    zeta = 0.0;
// 227 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
/*-------------------------------------------------------------------
c---->
c  Do one iteration untimed to init all code and data page tables
c---->                    (then reinit, start timing, to niter its)
c-------------------------------------------------------------------*/
// 228 lv-analysis-out: bot
    for (
// 229 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
it = 1
// 229 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
; 
// 230 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
it <= 1;
// 230 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
 it++) {
/*--------------------------------------------------------------------
c  The call to the conjugate gradient routine:
c-------------------------------------------------------------------*/
// 233 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
      conj_grad(colidx,rowstr,x,z,a,p,q,r,w,&rnorm);
// 233 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_326,tmp_327,tmp_328,tmp_329,tmp_330,tmp_331,tmp_332,tmp_333,tmp_334,tmp_335}
/*--------------------------------------------------------------------
c  zeta = shift + 1/(x.z)
c  So, first: (x.z)
c  Also, find norm of z
c  So, first: (z.z)
c-------------------------------------------------------------------*/
// 235 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
      
#pragma omp single
// 235 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
{
// 237 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
        norm_temp11 = 0.0;
// 237 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,cclass_179,zeta_verify_value_181,start_304}
// 238 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,cclass_179,zeta_verify_value_181,start_304}
        norm_temp12 = 0.0;
// 238 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
/* end single */
      }
// 239 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
      
#pragma omp for reduction(+:norm_temp11,norm_temp12)
// 239 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
// 240 lv-analysis-out: bot
      for (
// 241 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
j = 1
// 241 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
; 
// 242 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
j <= lastcol - firstcol + 1;
// 242 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
 j++) {
// 245 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
        norm_temp11 = norm_temp11 + x[j] * z[j];
// 245 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
// 246 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
        norm_temp12 = norm_temp12 + z[j] * z[j];
// 246 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
      }
// 240 lv-analysis-in : bot
// 247 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
      
#pragma omp single
// 247 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
// 248 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
      norm_temp12 = 1.0 / sqrt(norm_temp12);
// 248 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
/*--------------------------------------------------------------------
c  Normalize z to obtain x
c-------------------------------------------------------------------*/
// 249 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
      
#pragma omp for
// 249 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
// 250 lv-analysis-out: bot
      for (
// 251 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
j = 1
// 251 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
; 
// 252 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
j <= lastcol - firstcol + 1;
// 252 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
 j++) {
// 255 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
        x[j] = norm_temp12 * z[j];
// 255 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304}
      }
// 250 lv-analysis-in : bot
/* end of do one iteration untimed */
    }
// 228 lv-analysis-in : bot
/*--------------------------------------------------------------------
c  set starting vector to (1, 1, .... 1)
c-------------------------------------------------------------------*/
// 256 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181}
    
#pragma omp for nowait
// 256 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181}
// 257 lv-analysis-out: bot
    for (
// 258 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181}
i = 1
// 258 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,i_168,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181}
; 
// 259 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,i_168,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181}
i <= 75000 + 1;
// 259 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,i_168,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181}
 i++) {
// 262 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,z_161,p_162,q_163,r_164,w_165,i_168,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181}
      x[i] = 1.0;
// 262 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,i_168,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181}
    }
// 257 lv-analysis-in : bot
// 263 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181}
    
#pragma omp single
// 263 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181}
// 264 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181}
    zeta = 0.0;
// 264 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181}
/* end parallel */
  }
// 265 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181}
  timer_clear(1);
// 265 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181,tmp_326}
// 267 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181}
  timer_start(1);
// 267 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181,tmp_326}
/*--------------------------------------------------------------------
c---->
c  Main Iteration for inverse power method
c---->
c-------------------------------------------------------------------*/
// 269 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_346}
  
#pragma omp parallel private(it,i,j,k)
// 269 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_346}
{
// 271 lv-analysis-out: bot
    for (
// 272 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_346}
it = 1
// 272 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_346}
; 
// 273 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_346}
it <= 75;
// 273 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_346}
 it++) {
/*--------------------------------------------------------------------
c  The call to the conjugate gradient routine:
c-------------------------------------------------------------------*/
// 276 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
      conj_grad(colidx,rowstr,x,z,a,p,q,r,w,&rnorm);
// 276 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_326,tmp_327,tmp_328,tmp_329,tmp_330,tmp_331,tmp_332,tmp_333,tmp_334,tmp_335}
/*--------------------------------------------------------------------
c  zeta = shift + 1/(x.z)
c  So, first: (x.z)
c  Also, find norm of z
c  So, first: (z.z)
c-------------------------------------------------------------------*/
// 278 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_346}
      
#pragma omp single
// 278 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_346}
{
// 280 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_346}
        norm_temp11 = 0.0;
// 280 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,cclass_179,zeta_verify_value_181,start_304,tmp_346}
// 281 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,cclass_179,zeta_verify_value_181,start_304,tmp_346}
        norm_temp12 = 0.0;
// 281 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
/* end single */
      }
// 282 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
      
#pragma omp for reduction(+:norm_temp11,norm_temp12)
// 282 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
// 283 lv-analysis-out: bot
      for (
// 284 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
j = 1
// 284 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
; 
// 285 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
j <= lastcol - firstcol + 1;
// 285 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
 j++) {
// 288 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
        norm_temp11 = norm_temp11 + x[j] * z[j];
// 288 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
// 289 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
        norm_temp12 = norm_temp12 + z[j] * z[j];
// 289 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
      }
// 283 lv-analysis-in : bot
// 290 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,cclass_179,zeta_verify_value_181,start_304,tmp_346}
      
#pragma omp single
// 290 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,cclass_179,zeta_verify_value_181,start_304,tmp_346}
{
// 292 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,cclass_179,zeta_verify_value_181,start_304,tmp_346}
        norm_temp12 = 1.0 / sqrt(norm_temp12);
// 292 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
// 293 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,norm_temp11_175,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
        zeta = 60.0 + 1.0 / norm_temp11;
// 293 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
/* end single */
      }
// 294 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
      
#pragma omp master
// 294 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
{
// 296 lv-analysis-out: bot
        if (
// 297 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
it == 1
// 297 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
) {
// 299 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
          printf("   iteration           ||r||                 zeta\n");
// 299 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
        }
// 296 lv-analysis-in : bot
// 301 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
        printf("    %5d       %20.14e%20.13e\n",it,rnorm,zeta);
// 301 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
/* end master */
      }
/*--------------------------------------------------------------------
c  Normalize z to obtain x
c-------------------------------------------------------------------*/
// 303 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
      
#pragma omp for
// 303 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
// 304 lv-analysis-out: bot
      for (
// 305 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
j = 1
// 305 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
; 
// 306 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
j <= lastcol - firstcol + 1;
// 306 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
 j++) {
// 309 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
        x[j] = norm_temp12 * z[j];
// 309 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,j_169,it_171,nthreads_172,zeta_173,rnorm_174,norm_temp12_176,cclass_179,zeta_verify_value_181,start_304,tmp_346}
      }
// 304 lv-analysis-in : bot
/* end of main iter inv pow meth */
    }
// 271 lv-analysis-in : bot
#if defined(_OPENMP)
#endif /* _OPENMP */
/* end parallel */
  }
// 310 lv-analysis-out: {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,start_304,tmp_346}
  timer_stop(1);
// 310 lv-analysis-in : {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,start_304,tmp_326,tmp_346}
/*--------------------------------------------------------------------
c  End of timed section
c-------------------------------------------------------------------*/
// 312 lv-analysis-out: {nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,elapsed_305}
  t = timer_read(1);
// 312 lv-analysis-in : {nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,elapsed_305,tmp_326}
// 314 lv-analysis-out: {nthreads_172,zeta_173,t_177,cclass_179,zeta_verify_value_181,tmp_346}
  printf(" Benchmark completed\n");
// 314 lv-analysis-in : {nthreads_172,zeta_173,t_177,cclass_179,zeta_verify_value_181,tmp_346}
// 316 lv-analysis-out: {nthreads_172,zeta_173,t_177,cclass_179,zeta_verify_value_181,tmp_346}
  epsilon = 1.0e-10;
// 316 lv-analysis-in : {nthreads_172,zeta_173,t_177,cclass_179,zeta_verify_value_181,epsilon_182,tmp_346}
// 317 lv-analysis-out: bot
  if (
// 318 lv-analysis-out: {nthreads_172,zeta_173,t_177,cclass_179,zeta_verify_value_181,epsilon_182,tmp_346}
cclass != 'U'
// 318 lv-analysis-in : {nthreads_172,zeta_173,t_177,cclass_179,zeta_verify_value_181,epsilon_182,tmp_346}
) {
// 320 lv-analysis-out: {nthreads_172,zeta_173,t_177,cclass_179,zeta_verify_value_181,epsilon_182,tmp_346}
    double __temp0__ = zeta - zeta_verify_value;
// 320 lv-analysis-in : {nthreads_172,zeta_173,t_177,cclass_179,zeta_verify_value_181,epsilon_182,tmp_346}
// 321 lv-analysis-out: bot
    if (
// 322 lv-analysis-out: {nthreads_172,zeta_173,t_177,cclass_179,zeta_verify_value_181,epsilon_182,tmp_346}
fabs(__temp0__) <= epsilon
// 322 lv-analysis-in : {nthreads_172,zeta_173,t_177,cclass_179,zeta_verify_value_181,tmp_346}
) {
// 324 lv-analysis-out: {nthreads_172,zeta_173,t_177,cclass_179,zeta_verify_value_181,tmp_346}
      verified = 1;
// 324 lv-analysis-in : {nthreads_172,zeta_173,t_177,cclass_179,verified_180,zeta_verify_value_181,tmp_346}
// 325 lv-analysis-out: {nthreads_172,zeta_173,t_177,cclass_179,verified_180,zeta_verify_value_181,tmp_346}
      printf(" VERIFICATION SUCCESSFUL\n");
// 325 lv-analysis-in : {nthreads_172,zeta_173,t_177,cclass_179,verified_180,zeta_verify_value_181,tmp_346}
// 327 lv-analysis-out: {nthreads_172,zeta_173,t_177,cclass_179,verified_180,zeta_verify_value_181,tmp_346}
      printf(" Zeta is    %20.12e\n",zeta);
// 327 lv-analysis-in : {nthreads_172,zeta_173,t_177,cclass_179,verified_180,zeta_verify_value_181,tmp_346}
// 329 lv-analysis-out: {nthreads_172,zeta_173,t_177,cclass_179,verified_180,zeta_verify_value_181,tmp_346}
      double __temp1__ = zeta - zeta_verify_value;
// 329 lv-analysis-in : {nthreads_172,t_177,cclass_179,verified_180,__temp1___184,tmp_346}
// 330 lv-analysis-out: {nthreads_172,t_177,cclass_179,verified_180,__temp1___184,tmp_346}
      printf(" Error is   %20.12e\n",__temp1__);
// 330 lv-analysis-in : {nthreads_172,t_177,cclass_179,verified_180,tmp_346}
    }
     else {
// 333 lv-analysis-out: {nthreads_172,zeta_173,t_177,cclass_179,zeta_verify_value_181,tmp_346}
      verified = 0;
// 333 lv-analysis-in : {nthreads_172,zeta_173,t_177,cclass_179,verified_180,zeta_verify_value_181,tmp_346}
// 334 lv-analysis-out: {nthreads_172,zeta_173,t_177,cclass_179,verified_180,zeta_verify_value_181,tmp_346}
      printf(" VERIFICATION FAILED\n");
// 334 lv-analysis-in : {nthreads_172,zeta_173,t_177,cclass_179,verified_180,zeta_verify_value_181,tmp_346}
// 336 lv-analysis-out: {nthreads_172,zeta_173,t_177,cclass_179,verified_180,zeta_verify_value_181,tmp_346}
      printf(" Zeta                %20.12e\n",zeta);
// 336 lv-analysis-in : {nthreads_172,t_177,cclass_179,verified_180,zeta_verify_value_181,tmp_346}
// 338 lv-analysis-out: {nthreads_172,t_177,cclass_179,verified_180,zeta_verify_value_181,tmp_346}
      printf(" The correct zeta is %20.12e\n",zeta_verify_value);
// 338 lv-analysis-in : {nthreads_172,t_177,cclass_179,verified_180,tmp_346}
    }
// 321 lv-analysis-in : bot
  }
   else {
// 341 lv-analysis-out: {nthreads_172,t_177,cclass_179,tmp_346}
    verified = 0;
// 341 lv-analysis-in : {nthreads_172,t_177,cclass_179,verified_180,tmp_346}
// 342 lv-analysis-out: {nthreads_172,t_177,cclass_179,verified_180,tmp_346}
    printf(" Problem size unknown\n");
// 342 lv-analysis-in : {nthreads_172,t_177,cclass_179,verified_180,tmp_346}
// 344 lv-analysis-out: {nthreads_172,t_177,cclass_179,verified_180,tmp_346}
    printf(" NO VERIFICATION PERFORMED\n");
// 344 lv-analysis-in : {nthreads_172,t_177,cclass_179,verified_180,tmp_346}
  }
// 317 lv-analysis-in : bot
// 346 lv-analysis-out: bot
  if (
// 347 lv-analysis-out: {nthreads_172,t_177,cclass_179,verified_180,tmp_346}
t != 0.0
// 347 lv-analysis-in : {nthreads_172,t_177,cclass_179,verified_180,tmp_346}
) {
// 349 lv-analysis-out: {nthreads_172,t_177,cclass_179,verified_180,tmp_346}
    mflops = 2.0 * 75 * 75000 * (3.0 + (13 * (13 + 1)) + 25.0 * (5.0 + (13 * (13 + 1))) + 3.0) / t / 1000000.0;
// 349 lv-analysis-in : {nthreads_172,t_177,mflops_178,cclass_179,verified_180,tmp_346}
  }
   else {
// 351 lv-analysis-out: {nthreads_172,t_177,cclass_179,verified_180,tmp_346}
    mflops = 0.0;
// 351 lv-analysis-in : {nthreads_172,t_177,mflops_178,cclass_179,verified_180,tmp_346}
  }
// 346 lv-analysis-in : bot
// 352 lv-analysis-out: {nthreads_172,t_177,mflops_178,cclass_179,verified_180,tmp_346}
  c_print_results("CG",cclass,75000,0,0,75,nthreads,t,mflops,"          floating point",verified,"2.3","28 Oct 2014","gcc","$(CC)","(none)","-I../common","-fopenmp -O2","-lm -fopenmp","randdp");
// 352 lv-analysis-in : {tmp_326,tmp_327,tmp_328,tmp_329,tmp_330,tmp_331,tmp_332,tmp_333,tmp_334,tmp_335,tmp_336,tmp_337,tmp_338,tmp_339,tmp_340,tmp_341,tmp_342,tmp_343,tmp_344,tmp_345,tmp_346}
}
// 138 lv-analysis-in : {tv_sec_72,tv_usec_73,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_346}
/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void conj_grad(
/* colidx[1:nzz] */
int colidx[],
/* rowstr[1:naa+1] */
int rowstr[],
/* x[*] */
double x[],
/* z[*] */
double z[],
/* a[1:nzz] */
double a[],
/* p[*] */
double p[],
/* q[*] */
double q[],
/* r[*] */
double r[],
/* w[*] */
double w[],double *rnorm)
/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/
/*---------------------------------------------------------------------
c  Floaging point arrays here are named as in NPB1 spec discussion of 
c  CG algorithm
c---------------------------------------------------------------------*/
// 354 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_326,tmp_327,tmp_328,tmp_329,tmp_330,tmp_331,tmp_332,tmp_333,tmp_334,tmp_335}
{
// 357 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  static double d;
// 357 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 358 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  static double sum;
// 358 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 359 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  static double rho;
// 359 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 360 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  static double rho0;
// 360 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 361 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  static double alpha;
// 361 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 362 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  static double beta;
// 362 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 363 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  int i;
// 363 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 364 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  int j;
// 364 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 365 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  int k;
// 365 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 366 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  int cgit;
// 366 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 367 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  int cgitmax = 25;
// 367 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 368 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  
#pragma omp single nowait
// 368 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 369 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  rho = 0.0;
// 369 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
/*--------------------------------------------------------------------
c  Initialize the CG algorithm:
c-------------------------------------------------------------------*/
// 370 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  
#pragma omp for nowait
// 370 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 371 lv-analysis-out: bot
  for (
// 372 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j = 1
// 372 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
; 
// 373 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j <= naa + 1;
// 373 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
 j++) {
// 376 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,x_199,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    q[j] = 0.0;
// 376 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,x_199,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 377 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,x_199,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    z[j] = 0.0;
// 377 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,z_197,x_199,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 378 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,z_197,x_199,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    r[j] = x[j];
// 378 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,z_197,r_198,x_199,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 379 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,z_197,r_198,x_199,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    p[j] = r[j];
// 379 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,z_197,r_198,x_199,p_200,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 380 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,z_197,r_198,x_199,p_200,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    w[j] = 0.0;
// 380 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  }
// 371 lv-analysis-in : bot
/*--------------------------------------------------------------------
c  rho = r.r
c  Now, obtain the norm of r: First, sum squares of r elements locally...
c-------------------------------------------------------------------*/
// 381 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  
#pragma omp for reduction(+:rho)
// 381 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 382 lv-analysis-out: bot
  for (
// 383 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j = 1
// 383 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
; 
// 384 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j <= lastcol - firstcol + 1;
// 384 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
 j++) {
// 387 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    rho = rho + x[j] * x[j];
// 387 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,j_192,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  }
// 382 lv-analysis-in : bot
/*--------------------------------------------------------------------
c---->
c  The conj grad iteration loop
c---->
c-------------------------------------------------------------------*/
// 388 lv-analysis-out: bot
  for (
// 389 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
cgit = 1
// 389 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
; 
// 390 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
cgit <= cgitmax;
// 390 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
 cgit++) {
// 393 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    
#pragma omp single nowait
// 393 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
{
// 395 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
      rho0 = rho;
// 395 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 396 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
      d = 0.0;
// 396 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 397 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
      rho = 0.0;
// 397 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
/* end single */
    }
/*--------------------------------------------------------------------
c  q = A.p
c  The partition submatrix-vector multiply: use workspace w
c---------------------------------------------------------------------
C
C  NOTE: this version of the multiply is actually (slightly: maybe %5) 
C        faster on the sp2 on 16 nodes than is the unrolled-by-2 version 
C        below.   On the Cray t3d, the reverse is true, i.e., the 
C        unrolled-by-two version is some 10% faster.  
C        The unrolled-by-8 version below is significantly faster
C        on the Cray t3d - overall speed of code is 1.5 times faster.
*/
/* rolled version */
// 398 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    
#pragma omp for private(sum,k)
// 398 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 399 lv-analysis-out: bot
    for (
// 400 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j = 1
// 400 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
; 
// 401 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j <= lastrow - firstrow + 1;
// 401 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
 j++) {
// 404 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
      sum = 0.0;
// 404 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 405 lv-analysis-out: bot
      for (
// 406 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
k = rowstr[j]
// 406 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,rho_187,rho0_188,j_192,k_193,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
; 
// 407 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,rho_187,rho0_188,j_192,k_193,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
k < rowstr[j + 1];
// 407 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,rho_187,rho0_188,j_192,k_193,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
 k++) {
// 410 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,rho_187,rho0_188,j_192,k_193,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
        sum = sum + a[k] * p[colidx[k]];
// 410 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,rho_187,rho0_188,j_192,k_193,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
      }
// 405 lv-analysis-in : bot
// 411 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
      w[j] = sum;
// 411 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    }
// 399 lv-analysis-in : bot
/* unrolled-by-two version
#pragma omp for private(i,k)
        for (j = 1; j <= lastrow-firstrow+1; j++) {
	    int iresidue;
	    double sum1, sum2;
	    i = rowstr[j]; 
            iresidue = (rowstr[j+1]-i) % 2;
            sum1 = 0.0;
            sum2 = 0.0;
            if (iresidue == 1) sum1 = sum1 + a[i]*p[colidx[i]];
	    for (k = i+iresidue; k <= rowstr[j+1]-2; k += 2) {
		sum1 = sum1 + a[k]   * p[colidx[k]];
		sum2 = sum2 + a[k+1] * p[colidx[k+1]];
	    }
            w[j] = sum1 + sum2;
        }
*/
/* unrolled-by-8 version
#pragma omp for private(i,k,sum)
        for (j = 1; j <= lastrow-firstrow+1; j++) {
	    int iresidue;
            i = rowstr[j]; 
            iresidue = (rowstr[j+1]-i) % 8;
            sum = 0.0;
            for (k = i; k <= i+iresidue-1; k++) {
                sum = sum +  a[k] * p[colidx[k]];
            }
            for (k = i+iresidue; k <= rowstr[j+1]-8; k += 8) {
                sum = sum + a[k  ] * p[colidx[k  ]]
                          + a[k+1] * p[colidx[k+1]]
                          + a[k+2] * p[colidx[k+2]]
                          + a[k+3] * p[colidx[k+3]]
                          + a[k+4] * p[colidx[k+4]]
                          + a[k+5] * p[colidx[k+5]]
                          + a[k+6] * p[colidx[k+6]]
                          + a[k+7] * p[colidx[k+7]];
            }
            w[j] = sum;
        }
*/
// 412 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    
#pragma omp for
// 412 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 413 lv-analysis-out: bot
    for (
// 414 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j = 1
// 414 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
; 
// 415 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j <= lastcol - firstcol + 1;
// 415 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
 j++) {
// 418 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
      q[j] = w[j];
// 418 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    }
// 413 lv-analysis-in : bot
/*--------------------------------------------------------------------
c  Clear w for reuse...
c-------------------------------------------------------------------*/
// 419 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    
#pragma omp for nowait
// 419 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 420 lv-analysis-out: bot
    for (
// 421 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j = 1
// 421 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
; 
// 422 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j <= lastcol - firstcol + 1;
// 422 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
 j++) {
// 425 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
      w[j] = 0.0;
// 425 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    }
// 420 lv-analysis-in : bot
/*--------------------------------------------------------------------
c  Obtain p.q
c-------------------------------------------------------------------*/
// 426 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    
#pragma omp for reduction(+:d)
// 426 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 427 lv-analysis-out: bot
    for (
// 428 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j = 1
// 428 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
; 
// 429 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j <= lastcol - firstcol + 1;
// 429 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
 j++) {
// 432 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
      d = d + p[j] * q[j];
// 432 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    }
// 427 lv-analysis-in : bot
/*--------------------------------------------------------------------
c  Obtain alpha = rho / (p.q)
c-------------------------------------------------------------------*/
// 433 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    
#pragma omp single
// 433 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 434 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    alpha = rho0 / d;
// 434 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,alpha_189,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
/*--------------------------------------------------------------------
c  Save a temporary of rho
c-------------------------------------------------------------------*/
/*	rho0 = rho;*/
/*---------------------------------------------------------------------
c  Obtain z = z + alpha*p
c  and    r = r - alpha*q
c---------------------------------------------------------------------*/
// 435 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,alpha_189,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    
#pragma omp for
// 435 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,alpha_189,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 436 lv-analysis-out: bot
    for (
// 437 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,alpha_189,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j = 1
// 437 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,alpha_189,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
; 
// 438 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,alpha_189,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j <= lastcol - firstcol + 1;
// 438 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,alpha_189,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
 j++) {
// 441 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,alpha_189,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
      z[j] = z[j] + alpha * p[j];
// 441 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,alpha_189,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 442 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,alpha_189,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
      r[j] = r[j] - alpha * q[j];
// 442 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,alpha_189,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    }
// 436 lv-analysis-in : bot
/*---------------------------------------------------------------------
c  rho = r.r
c  Now, obtain the norm of r: First, sum squares of r elements locally...
c---------------------------------------------------------------------*/
// 443 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    
#pragma omp for reduction(+:rho)
// 443 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 444 lv-analysis-out: bot
    for (
// 445 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j = 1
// 445 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
; 
// 446 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j <= lastcol - firstcol + 1;
// 446 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
 j++) {
// 449 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
      rho = rho + r[j] * r[j];
// 449 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    }
// 444 lv-analysis-in : bot
/*--------------------------------------------------------------------
c  Obtain beta:
c-------------------------------------------------------------------*/
// 450 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    
#pragma omp single
// 450 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 451 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,rho0_188,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    beta = rho / rho0;
// 451 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,beta_190,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
/*--------------------------------------------------------------------
c  p = r + beta*p
c-------------------------------------------------------------------*/
// 452 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,beta_190,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    
#pragma omp for
// 452 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,beta_190,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 453 lv-analysis-out: bot
    for (
// 454 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,beta_190,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j = 1
// 454 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,beta_190,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
; 
// 455 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,beta_190,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j <= lastcol - firstcol + 1;
// 455 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,beta_190,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
 j++) {
// 458 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,beta_190,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
      p[j] = r[j] + beta * p[j];
// 458 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rho_187,beta_190,j_192,cgit_194,cgitmax_195,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    }
// 453 lv-analysis-in : bot
/* end of do cgit=1,cgitmax */
  }
// 388 lv-analysis-in : bot
/*---------------------------------------------------------------------
c  Compute residual norm explicitly:  ||r|| = ||x - A.z||
c  First, form A.z
c  The partition submatrix-vector multiply
c---------------------------------------------------------------------*/
// 459 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,z_197,r_198,x_199,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  
#pragma omp single nowait
// 459 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,z_197,r_198,x_199,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 460 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,z_197,r_198,x_199,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  sum = 0.0;
// 460 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,z_197,r_198,x_199,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 461 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,z_197,r_198,x_199,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  
#pragma omp for private(d, k)
// 461 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,z_197,r_198,x_199,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 462 lv-analysis-out: bot
  for (
// 463 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,z_197,r_198,x_199,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j = 1
// 463 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,j_192,z_197,r_198,x_199,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
; 
// 464 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,j_192,z_197,r_198,x_199,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
j <= lastrow - firstrow + 1;
// 464 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,j_192,z_197,r_198,x_199,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
 j++) {
// 467 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,j_192,z_197,r_198,x_199,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    d = 0.0;
// 467 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,j_192,z_197,r_198,x_199,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
// 468 lv-analysis-out: bot
    for (
// 469 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,j_192,z_197,r_198,x_199,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
k = rowstr[j]
// 469 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,j_192,k_193,z_197,r_198,x_199,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
; 
// 470 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,j_192,k_193,z_197,r_198,x_199,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
k <= rowstr[j + 1] - 1;
// 470 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,j_192,k_193,z_197,r_198,x_199,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
 k++) {
// 473 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,j_192,k_193,z_197,r_198,x_199,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
      d = d + a[k] * z[colidx[k]];
// 473 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,j_192,k_193,z_197,r_198,x_199,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    }
// 468 lv-analysis-in : bot
// 474 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,j_192,z_197,r_198,x_199,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
    w[j] = d;
// 474 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,j_192,z_197,r_198,x_199,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
  }
// 462 lv-analysis-in : bot
// 475 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,r_198,x_199,w_201,rnorm_205,start_304}
  
#pragma omp for
// 475 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,r_198,x_199,w_201,rnorm_205,start_304}
// 476 lv-analysis-out: bot
  for (
// 477 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,r_198,x_199,w_201,rnorm_205,start_304}
j = 1
// 477 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,j_192,r_198,x_199,w_201,rnorm_205,start_304}
; 
// 478 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,j_192,r_198,x_199,w_201,rnorm_205,start_304}
j <= lastcol - firstcol + 1;
// 478 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,j_192,r_198,x_199,w_201,rnorm_205,start_304}
 j++) {
// 481 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,j_192,x_199,w_201,rnorm_205,start_304}
    r[j] = w[j];
// 481 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,j_192,r_198,x_199,w_201,rnorm_205,start_304}
  }
// 476 lv-analysis-in : bot
/*--------------------------------------------------------------------
c  At this point, r contains A.z
c-------------------------------------------------------------------*/
// 482 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,r_198,x_199,rnorm_205,start_304}
  
#pragma omp for reduction(+:sum) private(d)
// 482 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,r_198,x_199,rnorm_205,start_304}
// 483 lv-analysis-out: bot
  for (
// 484 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,r_198,x_199,rnorm_205,start_304}
j = 1
// 484 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,j_192,r_198,x_199,rnorm_205,start_304}
; 
// 485 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,j_192,r_198,x_199,rnorm_205,start_304}
j <= lastcol - firstcol + 1;
// 485 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,j_192,r_198,x_199,rnorm_205,start_304}
 j++) {
// 488 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,j_192,r_198,x_199,rnorm_205,start_304}
    d = x[j] - r[j];
// 488 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,j_192,r_198,x_199,rnorm_205,start_304}
// 489 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,d_185,sum_186,j_192,r_198,x_199,rnorm_205,start_304}
    sum = sum + d * d;
// 489 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,sum_186,j_192,r_198,x_199,rnorm_205,start_304}
  }
// 483 lv-analysis-in : bot
// 490 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rnorm_205,start_304}
  
#pragma omp single
// 490 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rnorm_205,start_304}
{
// 492 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,rnorm_205,start_304}
     *rnorm = sqrt(sum);
// 492 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
/* end single */
  }
}
// 354 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,it_171,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,q_196,z_197,r_198,x_199,p_200,w_201,rowstr_202,a_203,colidx_204,rnorm_205,start_304}
/*---------------------------------------------------------------------
c       generate the test problem for benchmark 6
c       makea generates a sparse matrix with a
c       prescribed sparsity distribution
c
c       parameter    type        usage
c
c       input
c
c       n            i           number of cols/rows of matrix
c       nz           i           nonzeros as declared array size
c       rcond        r*8         condition number
c       shift        r*8         main diagonal shift
c
c       output
c
c       a            r*8         array for nonzeros
c       colidx       i           col indices
c       rowstr       i           row pointers
c
c       workspace
c
c       iv, arow, acol i
c       v, aelt        r*8
c---------------------------------------------------------------------*/

static void makea(int n,int nz,
/* a[1:nz] */
double a[],
/* colidx[1:nz] */
int colidx[],
/* rowstr[1:n+1] */
int rowstr[],int nonzer,int firstrow,int lastrow,int firstcol,int lastcol,double rcond,
/* arow[1:nz] */
int arow[],
/* acol[1:nz] */
int acol[],
/* aelt[1:nz] */
double aelt[],
/* v[1:n+1] */
double v[],
/* iv[1:2*n+1] */
int iv[],double shift)
// 493 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_326,tmp_327,tmp_328,tmp_329,tmp_330,tmp_331,tmp_332,tmp_333,tmp_334,tmp_335,tmp_336,tmp_337,tmp_338,tmp_339,tmp_340,tmp_341,tmp_342,tmp_346}
{
// 496 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  int i;
// 496 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 497 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  int nnza;
// 497 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 498 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  int iouter;
// 498 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 499 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  int ivelt;
// 499 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 500 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  int ivelt1;
// 500 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 501 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  int irow;
// 501 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 502 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  int nzv;
// 502 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
/*--------------------------------------------------------------------
c      nonzer is approximately  (int(sqrt(nnza /n)));
c-------------------------------------------------------------------*/
// 503 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  double size;
// 503 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 504 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  double ratio;
// 504 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 505 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  double scale;
// 505 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 506 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  int jcol;
// 506 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 507 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  size = 1.0;
// 507 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,size_213,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 508 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,size_213,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  double __temp2__ = 1.0 / ((double )n);
// 508 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,size_213,__temp2___217,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 509 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,size_213,__temp2___217,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  ratio = pow(rcond,__temp2__);
// 509 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,size_213,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 511 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  nnza = 0;
// 511 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
/*---------------------------------------------------------------------
c  Initialize colidx(n+1 .. 2n) to zero.
c  Used by sprnvc to mark nonzero positions
c---------------------------------------------------------------------*/
// 512 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  
#pragma omp parallel for
// 512 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 513 lv-analysis-out: bot
  for (
// 514 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
i = 1
// 514 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
; 
// 515 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
i <= n;
// 515 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
 i++) {
// 518 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,size_213,ratio_214,n_218,rcond_219,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
    colidx[n + i] = 0;
// 518 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  }
// 513 lv-analysis-in : bot
// 519 lv-analysis-out: bot
  for (
// 520 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
iouter = 1
// 520 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
; 
// 521 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
iouter <= n;
// 521 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
 iouter++) {
// 524 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304}
    nzv = nonzer;
// 524 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304}
// 525 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304}
    int *__temp3__ = &colidx[0];
// 525 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,__temp3___222,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304}
// 526 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,__temp3___222,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304}
    int *__temp4__ = &colidx[n];
// 526 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,__temp3___222,__temp4___223,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304}
// 527 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,__temp3___222,__temp4___223,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304}
    sprnvc(n,nzv,v,iv,__temp3__,__temp4__);
// 527 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_326,tmp_327,tmp_328,tmp_329,tmp_330,tmp_331}
// 529 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304}
    vecset(n,v,iv,&nzv,iouter,0.5);
// 529 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_326,tmp_327,tmp_328,tmp_329,tmp_330,tmp_331}
// 531 lv-analysis-out: bot
    for (
// 532 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
ivelt = 1
// 532 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
; 
// 533 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
ivelt <= nzv;
// 533 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
 ivelt++) {
// 536 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
      jcol = iv[ivelt];
// 536 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,nzv_212,size_213,ratio_214,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 537 lv-analysis-out: bot
      if (
// 538 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,nzv_212,size_213,ratio_214,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
jcol >= firstcol && jcol <= lastcol
// 538 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,nzv_212,size_213,ratio_214,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
) {
// 540 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,nzv_212,size_213,ratio_214,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
        scale = size * v[ivelt];
// 540 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 541 lv-analysis-out: bot
        for (
// 542 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
ivelt1 = 1
// 542 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
; 
// 543 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
ivelt1 <= nzv;
// 543 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
 ivelt1++) {
// 546 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
          irow = iv[ivelt1];
// 546 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 547 lv-analysis-out: bot
          if (
// 548 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
irow >= firstrow && irow <= lastrow
// 548 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
) {
// 550 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
            nnza = nnza + 1;
// 550 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 551 lv-analysis-out: bot
            if (
// 552 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
nnza > nz
// 552 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
) {
// 554 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
              printf("Space for matrix elements exceeded in makea\n");
// 554 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 556 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
              printf("nnza, nzmax = %d, %d\n",nnza,nz);
// 556 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 558 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
              printf("iouter = %d\n",iouter);
// 558 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 560 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
              exit(1);
// 560 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
            }
// 551 lv-analysis-in : bot
// 562 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
            acol[nnza] = jcol;
// 562 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 563 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,irow_211,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,shift_234,a_237,rowstr_238,start_304,tmp_346}
            arow[nnza] = irow;
// 563 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 564 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,shift_234,a_237,rowstr_238,start_304,tmp_346}
            aelt[nnza] = v[ivelt1] * scale;
// 564 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,ivelt_209,ivelt1_210,nzv_212,size_213,ratio_214,scale_215,jcol_216,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
          }
// 547 lv-analysis-in : bot
        }
// 541 lv-analysis-in : bot
      }
// 537 lv-analysis-in : bot
    }
// 531 lv-analysis-in : bot
// 565 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
    size = size * ratio;
// 565 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
  }
// 519 lv-analysis-in : bot
/*---------------------------------------------------------------------
c       ... add the identity * rcond to the generated matrix to bound
c           the smallest eigenvalue from below by rcond
c---------------------------------------------------------------------*/
// 566 lv-analysis-out: bot
  for (
// 567 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
i = firstrow
// 567 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
; 
// 568 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
i <= lastrow;
// 568 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
 i++) {
// 571 lv-analysis-out: bot
    if (
// 572 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
i >= firstcol && i <= lastcol
// 572 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
) {
// 574 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
      iouter = n + i;
// 574 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,iouter_208,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 575 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,iouter_208,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
      nnza = nnza + 1;
// 575 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,iouter_208,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 576 lv-analysis-out: bot
      if (
// 577 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,iouter_208,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
nnza > nz
// 577 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,iouter_208,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
) {
// 579 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,iouter_208,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
        printf("Space for matrix elements exceeded in makea\n");
// 579 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,iouter_208,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 581 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,iouter_208,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
        printf("nnza, nzmax = %d, %d\n",nnza,nz);
// 581 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,iouter_208,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 583 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,iouter_208,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
        printf("iouter = %d\n",iouter);
// 583 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 585 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
        exit(1);
// 585 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
      }
// 576 lv-analysis-in : bot
// 587 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,shift_234,a_237,rowstr_238,start_304,tmp_346}
      acol[nnza] = i;
// 587 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 588 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,shift_234,a_237,rowstr_238,start_304,tmp_346}
      arow[nnza] = i;
// 588 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,shift_234,a_237,rowstr_238,start_304,tmp_346}
// 589 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,shift_234,a_237,rowstr_238,start_304,tmp_346}
      aelt[nnza] = rcond - shift;
// 589 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,i_206,nnza_207,n_218,rcond_219,colidx_220,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
    }
// 571 lv-analysis-in : bot
  }
// 566 lv-analysis-in : bot
/*---------------------------------------------------------------------
c       ... make the sparse matrix from list of elements with duplicates
c           (v and iv are used as  workspace)
c---------------------------------------------------------------------*/
// 590 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,n_218,colidx_220,v_224,iv_225,firstrow_228,lastrow_229,acol_231,arow_232,aelt_233,a_237,rowstr_238,start_304}
  int *__temp5__ = &iv[0];
// 590 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,n_218,colidx_220,v_224,iv_225,firstrow_228,lastrow_229,acol_231,arow_232,aelt_233,__temp5___235,a_237,rowstr_238,start_304}
// 591 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,n_218,colidx_220,v_224,iv_225,firstrow_228,lastrow_229,acol_231,arow_232,aelt_233,__temp5___235,a_237,rowstr_238,start_304}
  int *__temp6__ = &iv[n];
// 591 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,n_218,colidx_220,v_224,firstrow_228,lastrow_229,acol_231,arow_232,aelt_233,__temp5___235,__temp6___236,a_237,rowstr_238,start_304}
// 592 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,n_218,colidx_220,v_224,firstrow_228,lastrow_229,acol_231,arow_232,aelt_233,__temp5___235,__temp6___236,a_237,rowstr_238,start_304}
  sparse(a,colidx,rowstr,n,arow,acol,aelt,firstrow,lastrow,v,__temp5__,__temp6__,nnza);
// 592 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_326,tmp_327,tmp_328,tmp_329,tmp_330,tmp_331,tmp_332,tmp_333,tmp_334,tmp_335,tmp_336,tmp_337,tmp_338}
}
// 493 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_346}
/*---------------------------------------------------
c       generate a sparse matrix from a list of
c       [col, row, element] tri
c---------------------------------------------------*/

static void sparse(
/* a[1:*] */
double a[],
/* colidx[1:*] */
int colidx[],
/* rowstr[1:*] */
int rowstr[],int n,
/* arow[1:*] */
int arow[],
/* acol[1:*] */
int acol[],
/* aelt[1:*] */
double aelt[],int firstrow,int lastrow,
/* x[1:n] */
double x[],
/* mark[1:n] */
boolean mark[],
/* nzloc[1:n] */
int nzloc[],int nnza)
/*---------------------------------------------------------------------
c       rows range from firstrow to lastrow
c       the rowstr pointers are defined for nrows = lastrow-firstrow+1 values
c---------------------------------------------------------------------*/
// 594 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,start_304,tmp_326,tmp_327,tmp_328,tmp_329,tmp_330,tmp_331,tmp_332,tmp_333,tmp_334,tmp_335,tmp_336,tmp_337,tmp_338}
{
// 597 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
  int nrows;
// 597 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
// 598 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
  int i;
// 598 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
// 599 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
  int j;
// 599 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
// 600 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
  int jajp1;
// 600 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
// 601 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
  int nza;
// 601 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
// 602 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
  int k;
// 602 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
// 603 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
  int nzrow;
// 603 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
// 604 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
  double xi;
// 604 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
/*--------------------------------------------------------------------
c    how many rows of result
c-------------------------------------------------------------------*/
// 605 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
  nrows = lastrow - firstrow + 1;
// 605 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
/*--------------------------------------------------------------------
c     ...count the number of triples in each row
c-------------------------------------------------------------------*/
// 606 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
  
#pragma omp parallel for
// 606 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
// 607 lv-analysis-out: bot
  for (
// 608 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
j = 1
// 608 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
; 
// 609 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
j <= n;
// 609 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
 j++) {
// 612 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,firstrow_248,n_249,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
    rowstr[j] = 0;
// 612 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,firstrow_248,n_249,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
// 613 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,firstrow_248,n_249,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
    mark[j] = 0;
// 613 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
  }
// 607 lv-analysis-in : bot
// 614 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
  rowstr[n + 1] = 0;
// 614 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
// 615 lv-analysis-out: bot
  for (
// 616 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
nza = 1
// 616 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,nza_243,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
; 
// 617 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,nza_243,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
nza <= nnza;
// 617 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,nza_243,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
 nza++) {
// 620 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,nza_243,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
    j = arow[nza] - firstrow + 1 + 1;
// 620 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
// 621 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
    rowstr[j] = rowstr[j] + 1;
// 621 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,nza_243,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
  }
// 615 lv-analysis-in : bot
// 622 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
  rowstr[1] = 1;
// 622 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
// 623 lv-analysis-out: bot
  for (
// 624 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
j = 2
// 624 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
; 
// 625 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
j <= nrows + 1;
// 625 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
 j++) {
// 628 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
    rowstr[j] = rowstr[j] + rowstr[j - 1];
// 628 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
  }
// 623 lv-analysis-in : bot
/*---------------------------------------------------------------------
c     ... rowstr(j) now is the location of the first nonzero
c           of row j of a
c---------------------------------------------------------------------*/
/*--------------------------------------------------------------------
c     ... do a bucket sort of the triples on the row index
c-------------------------------------------------------------------*/
// 629 lv-analysis-out: bot
  for (
// 630 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
nza = 1
// 630 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,nza_243,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
; 
// 631 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,nza_243,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
nza <= nnza;
// 631 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,nza_243,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
 nza++) {
// 634 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,nza_243,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,aelt_255,acol_257,x_258,nzloc_259,start_304}
    j = arow[nza] - firstrow + 1;
// 634 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,aelt_255,acol_257,x_258,nzloc_259,start_304}
// 635 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,aelt_255,acol_257,x_258,nzloc_259,start_304}
    k = rowstr[j];
// 635 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,k_244,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,aelt_255,acol_257,x_258,nzloc_259,start_304}
// 636 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,k_244,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,aelt_255,acol_257,x_258,nzloc_259,start_304}
    a[k] = aelt[nza];
// 636 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,k_244,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,acol_257,x_258,nzloc_259,start_304}
// 637 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,k_244,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,acol_257,x_258,nzloc_259,start_304}
    colidx[k] = acol[nza];
// 637 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
// 638 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
    rowstr[j] = rowstr[j] + 1;
// 638 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,nza_243,firstrow_248,n_249,rowstr_250,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
  }
// 629 lv-analysis-in : bot
/*--------------------------------------------------------------------
c       ... rowstr(j) now points to the first element of row j+1
c-------------------------------------------------------------------*/
// 639 lv-analysis-out: bot
  for (
// 640 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,n_249,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
j = nrows
// 640 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,n_249,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
; 
// 641 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,n_249,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
j >= 1;
// 641 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,n_249,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
 j--) {
// 644 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,n_249,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
    rowstr[j + 1] = rowstr[j];
// 644 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,n_249,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
  }
// 639 lv-analysis-in : bot
// 645 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,n_249,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
  rowstr[1] = 1;
// 645 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,n_249,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
/*--------------------------------------------------------------------
c       ... generate the actual output rows by adding elements
c-------------------------------------------------------------------*/
// 646 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,n_249,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
  nza = 0;
// 646 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,nza_243,n_249,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
// 647 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,nza_243,n_249,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
  
#pragma omp parallel for
// 647 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,nza_243,n_249,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
// 648 lv-analysis-out: bot
  for (
// 649 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,nza_243,n_249,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
i = 1
// 649 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,nza_243,n_249,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
; 
// 650 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,nza_243,n_249,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
i <= n;
// 650 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,nza_243,n_249,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
 i++) {
// 653 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,nza_243,n_249,rowstr_250,a_254,colidx_256,nzloc_259,start_304}
    x[i] = 0.0;
// 653 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,nza_243,n_249,rowstr_250,a_254,colidx_256,x_258,nzloc_259,start_304}
// 654 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,nza_243,n_249,rowstr_250,a_254,colidx_256,x_258,nzloc_259,start_304}
    mark[i] = 0;
// 654 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,nza_243,n_249,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
  }
// 648 lv-analysis-in : bot
// 655 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,nza_243,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
  jajp1 = rowstr[1];
// 655 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,jajp1_242,nza_243,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
// 656 lv-analysis-out: bot
  for (
// 657 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,jajp1_242,nza_243,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
j = 1
// 657 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,jajp1_242,nza_243,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
; 
// 658 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,jajp1_242,nza_243,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
j <= nrows;
// 658 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,jajp1_242,nza_243,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
 j++) {
// 661 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,jajp1_242,nza_243,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
    nzrow = 0;
// 661 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,jajp1_242,nza_243,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
/*--------------------------------------------------------------------
c          ...loop over the jth row of a
c-------------------------------------------------------------------*/
// 662 lv-analysis-out: bot
    for (
// 663 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,jajp1_242,nza_243,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
k = jajp1
// 663 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
; 
// 664 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
k < rowstr[j + 1];
// 664 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
 k++) {
// 667 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
      i = colidx[k];
// 667 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
// 668 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
      x[i] = x[i] + a[k];
// 668 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
// 669 lv-analysis-out: bot
      if (
// 670 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
mark[i] == 0 && x[i] != 0.0
// 670 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
) {
// 672 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,rowstr_250,a_254,colidx_256,x_258,start_304}
        mark[i] = 1;
// 672 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,start_304}
// 673 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,start_304}
        nzrow = nzrow + 1;
// 673 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,start_304}
// 674 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,start_304}
        nzloc[nzrow] = i;
// 674 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
      }
// 669 lv-analysis-in : bot
    }
// 662 lv-analysis-in : bot
/*--------------------------------------------------------------------
c          ... extract the nonzeros of this row
c-------------------------------------------------------------------*/
// 675 lv-analysis-out: bot
    for (
// 676 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
k = 1
// 676 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
; 
// 677 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
k <= nzrow;
// 677 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
 k++) {
// 680 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,k_244,nzrow_245,rowstr_250,a_254,colidx_256,x_258,nzloc_259,start_304}
      i = nzloc[k];
// 680 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,rowstr_250,a_254,colidx_256,x_258,nzloc_259,start_304}
// 681 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,rowstr_250,a_254,colidx_256,x_258,nzloc_259,start_304}
      mark[i] = 0;
// 681 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
// 682 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
      xi = x[i];
// 682 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,xi_246,rowstr_250,mark_251,a_254,colidx_256,nzloc_259,start_304}
// 683 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,xi_246,rowstr_250,mark_251,a_254,colidx_256,nzloc_259,start_304}
      x[i] = 0.0;
// 683 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,xi_246,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
// 684 lv-analysis-out: bot
      if (
// 685 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,xi_246,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
xi != 0.0
// 685 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,xi_246,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
) {
// 687 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,xi_246,rowstr_250,mark_251,x_258,nzloc_259,start_304}
        nza = nza + 1;
// 687 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,xi_246,rowstr_250,mark_251,x_258,nzloc_259,start_304}
// 688 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,xi_246,rowstr_250,mark_251,x_258,nzloc_259,start_304}
        a[nza] = xi;
// 688 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,x_258,nzloc_259,start_304}
// 689 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,i_240,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,x_258,nzloc_259,start_304}
        colidx[nza] = i;
// 689 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,k_244,nzrow_245,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
      }
// 684 lv-analysis-in : bot
    }
// 675 lv-analysis-in : bot
// 690 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,nza_243,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
    jajp1 = rowstr[j + 1];
// 690 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,jajp1_242,nza_243,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
// 691 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,jajp1_242,nza_243,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
    rowstr[j + 1] = nza + rowstr[1];
// 691 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nrows_239,j_241,jajp1_242,nza_243,rowstr_250,mark_251,a_254,colidx_256,x_258,nzloc_259,start_304}
  }
// 656 lv-analysis-in : bot
}
// 594 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,lastrow_247,firstrow_248,n_249,mark_251,nnza_252,arow_253,a_254,aelt_255,colidx_256,acol_257,x_258,nzloc_259,start_304}
/*---------------------------------------------------------------------
c       generate a sparse n-vector (v, iv)
c       having nzv nonzeros
c
c       mark(i) is set to 1 if position i is nonzero.
c       mark is all zero on entry and is reset to all zero before exit
c       this corrects a performance bug found by John G. Lewis, caused by
c       reinitialization of mark on every one of the n calls to sprnvc
---------------------------------------------------------------------*/

static void sprnvc(int n,int nz,
/* v[1:*] */
double v[],
/* iv[1:*] */
int iv[],
/* nzloc[1:n] */
int nzloc[],
/* mark[1:n] */
int mark[])
// 692 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_326,tmp_327,tmp_328,tmp_329,tmp_330,tmp_331}
{
// 695 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,n_267,nz_268,mark_269,nzloc_270,start_304}
  int nn1;
// 695 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,n_267,nz_268,mark_269,nzloc_270,start_304}
// 696 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,n_267,nz_268,mark_269,nzloc_270,start_304}
  int nzrow;
// 696 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,n_267,nz_268,mark_269,nzloc_270,start_304}
// 697 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,n_267,nz_268,mark_269,nzloc_270,start_304}
  int nzv;
// 697 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,n_267,nz_268,mark_269,nzloc_270,start_304}
// 698 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,n_267,nz_268,mark_269,nzloc_270,start_304}
  int ii;
// 698 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,n_267,nz_268,mark_269,nzloc_270,start_304}
// 699 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,n_267,nz_268,mark_269,nzloc_270,start_304}
  int i;
// 699 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,n_267,nz_268,mark_269,nzloc_270,start_304}
// 700 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,n_267,nz_268,mark_269,nzloc_270,start_304}
  double vecelt;
// 700 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
// 701 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
  double vecloc;
// 701 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
// 702 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
  nzv = 0;
// 702 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
// 703 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
  nzrow = 0;
// 703 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
// 704 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
  nn1 = 1;
// 704 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nn1_260,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
// 705 lv-analysis-out: bot
  do {
// 707 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nn1_260,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
    nn1 = 2 * nn1;
// 707 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nn1_260,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
  }while (
// 708 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nn1_260,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
nn1 < n
// 708 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nn1_260,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
);
// 705 lv-analysis-in : bot
/*--------------------------------------------------------------------
c    nn1 is the smallest power of two not less than n
c-------------------------------------------------------------------*/
// 709 lv-analysis-out: bot
  while(
// 710 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
nzv < nz
// 710 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
){
// 712 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
    vecelt = randlc(&tran,amult);
// 712 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_326,tmp_327}
/*--------------------------------------------------------------------
c   generate an integer between 1 and n in a portable manner
c-------------------------------------------------------------------*/
// 714 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
    vecloc = randlc(&tran,amult);
// 714 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_326,tmp_327}
// 716 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
    i = icnvrt(vecloc,nn1) + 1;
// 716 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,i_264,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
//continue;
// 717 lv-analysis-out: bot
    if (
// 718 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,i_264,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
i > n
// 718 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,i_264,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
) {
// 720 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzloc_270,start_304}
      break; 
// 720 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzloc_270,start_304}
    }
// 717 lv-analysis-in : bot
/*--------------------------------------------------------------------
c  was this integer generated already?
c-------------------------------------------------------------------*/
// 721 lv-analysis-out: bot
    if (
// 722 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,i_264,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
mark[i] == 0
// 722 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,i_264,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
) {
// 724 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,i_264,vecelt_265,n_267,nz_268,start_304}
      mark[i] = 1;
// 724 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,i_264,vecelt_265,n_267,nz_268,mark_269,start_304}
// 725 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,i_264,vecelt_265,n_267,nz_268,mark_269,start_304}
      nzrow = nzrow + 1;
// 725 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,i_264,vecelt_265,n_267,nz_268,mark_269,start_304}
// 726 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,i_264,vecelt_265,n_267,nz_268,mark_269,start_304}
      nzloc[nzrow] = i;
// 726 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,i_264,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
// 727 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,i_264,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
      nzv = nzv + 1;
// 727 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,i_264,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
// 728 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,i_264,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
      v[nzv] = vecelt;
// 728 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,i_264,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
// 729 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,i_264,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
      iv[nzv] = i;
// 729 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
    }
// 721 lv-analysis-in : bot
  }
// 709 lv-analysis-in : bot
// 730 lv-analysis-out: bot
  for (
// 731 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzloc_270,start_304}
ii = 1
// 731 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,ii_263,nzloc_270,start_304}
; 
// 732 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,ii_263,nzloc_270,start_304}
ii <= nzrow;
// 732 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,ii_263,nzloc_270,start_304}
 ii++) {
// 735 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,ii_263,nzloc_270,start_304}
    i = nzloc[ii];
// 735 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,ii_263,i_264,nzloc_270,start_304}
// 736 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,ii_263,i_264,nzloc_270,start_304}
    mark[i] = 0;
// 736 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,ii_263,nzloc_270,start_304}
  }
// 730 lv-analysis-in : bot
}
// 692 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,n_267,nz_268,mark_269,nzloc_270,start_304}
/*---------------------------------------------------------------------
* scale a double precision number x in (0,1) by a power of 2 and chop it
*---------------------------------------------------------------------*/

static int icnvrt(double x,int ipwr2)
// 737 lv-analysis-out: bot
{
// 740 lv-analysis-out: bot
  return (int )(ipwr2 * x);
// 740 lv-analysis-in : bot
}
// 737 lv-analysis-in : bot
/*--------------------------------------------------------------------
c       set ith element of sparse vector (v, iv) with
c       nzv nonzeros to val
c-------------------------------------------------------------------*/

static void vecset(int n,
/* v[1:*] */
double v[],
/* iv[1:*] */
int iv[],int *nzv,int i,double val)
// 741 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304,tmp_326,tmp_327,tmp_328,tmp_329,tmp_330,tmp_331}
{
// 744 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzv_277,iv_278,i_279,val_281,start_304}
  int k;
// 744 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzv_277,iv_278,i_279,val_281,start_304}
// 745 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzv_277,iv_278,i_279,val_281,start_304}
  boolean set;
// 745 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzv_277,iv_278,i_279,val_281,start_304}
// 746 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzv_277,iv_278,i_279,val_281,start_304}
  set = 0;
// 746 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,set_276,nzv_277,iv_278,i_279,val_281,start_304}
// 747 lv-analysis-out: bot
  for (
// 748 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,set_276,nzv_277,iv_278,i_279,val_281,start_304}
k = 1
// 748 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,k_275,set_276,nzv_277,iv_278,i_279,val_281,start_304}
; 
// 749 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,k_275,set_276,nzv_277,iv_278,i_279,val_281,start_304}
k <=  *nzv;
// 749 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,k_275,set_276,nzv_277,iv_278,i_279,val_281,start_304}
 k++) {
// 752 lv-analysis-out: bot
    if (
// 753 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,k_275,set_276,nzv_277,iv_278,i_279,val_281,start_304}
iv[k] == i
// 753 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,k_275,set_276,nzv_277,iv_278,i_279,val_281,start_304}
) {
// 755 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,k_275,nzv_277,iv_278,i_279,val_281,start_304}
      v[k] = val;
// 755 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,k_275,nzv_277,iv_278,i_279,val_281,start_304}
// 756 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,k_275,nzv_277,iv_278,i_279,val_281,start_304}
      set = 1;
// 756 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,k_275,set_276,nzv_277,iv_278,i_279,val_281,start_304}
    }
// 752 lv-analysis-in : bot
  }
// 747 lv-analysis-in : bot
// 757 lv-analysis-out: bot
  if (
// 758 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,set_276,nzv_277,i_279,val_281,start_304}
set == 0
// 758 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzv_277,i_279,val_281,start_304}
) {
// 760 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzv_277,i_279,val_281,start_304}
     *nzv =  *nzv + 1;
// 760 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzv_277,i_279,val_281,start_304}
// 761 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzv_277,i_279,val_281,start_304}
    v[ *nzv] = val;
// 761 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzv_277,i_279,start_304}
// 762 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzv_277,i_279,start_304}
    iv[ *nzv] = i;
// 762 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,start_304}
  }
// 757 lv-analysis-in : bot
}
// 741 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzv_277,iv_278,i_279,val_281,start_304}
/* cat ./common/c_print_results.c */
/*****************************************************************/
/******     C  _  P  R  I  N  T  _  R  E  S  U  L  T  S     ******/
/*****************************************************************/

void c_print_results(char *name,char cclass,int n1,int n2,int n3,int niter,int nthreads,double t,double mops,char *optype,int passed_verification,char *npbversion,char *compiletime,char *cc,char *clink,char *c_lib,char *c_inc,char *cflags,char *clinkflags,char *rand)
// 763 lv-analysis-out: {tmp_326,tmp_327,tmp_328,tmp_329,tmp_330,tmp_331,tmp_332,tmp_333,tmp_334,tmp_335,tmp_336,tmp_337,tmp_338,tmp_339,tmp_340,tmp_341,tmp_342,tmp_343,tmp_344,tmp_345,tmp_346}
{
// 766 lv-analysis-out: {name_283,cclass_284,n2_285,n3_286,n1_287,niter_288,nthreads_289,t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  char *evalue = "1000";
// 766 lv-analysis-in : {name_283,cclass_284,n2_285,n3_286,n1_287,niter_288,nthreads_289,t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
// 767 lv-analysis-out: {name_283,cclass_284,n2_285,n3_286,n1_287,niter_288,nthreads_289,t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  printf("\n\n %s Benchmark Completed\n",name);
// 767 lv-analysis-in : {cclass_284,n2_285,n3_286,n1_287,niter_288,nthreads_289,t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
// 769 lv-analysis-out: {cclass_284,n2_285,n3_286,n1_287,niter_288,nthreads_289,t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  printf(" Class           =                        %c\n",cclass);
// 769 lv-analysis-in : {n2_285,n3_286,n1_287,niter_288,nthreads_289,t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
// 771 lv-analysis-out: bot
  if (
// 772 lv-analysis-out: {n2_285,n3_286,n1_287,niter_288,nthreads_289,t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
n2 == 0 && n3 == 0
// 772 lv-analysis-in : {n2_285,n3_286,n1_287,niter_288,nthreads_289,t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
) {
/* as in IS */
// 774 lv-analysis-out: {n1_287,niter_288,nthreads_289,t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
    printf(" Size            =             %12d\n",n1);
// 774 lv-analysis-in : {niter_288,nthreads_289,t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  }
   else {
// 777 lv-analysis-out: {n2_285,n3_286,n1_287,niter_288,nthreads_289,t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
    printf(" Size            =              %3dx%3dx%3d\n",n1,n2,n3);
// 777 lv-analysis-in : {niter_288,nthreads_289,t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  }
// 771 lv-analysis-in : bot
// 779 lv-analysis-out: {niter_288,nthreads_289,t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  printf(" Iterations      =             %12d\n",niter);
// 779 lv-analysis-in : {nthreads_289,t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
// 781 lv-analysis-out: {nthreads_289,t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  printf(" Threads         =             %12d\n",nthreads);
// 781 lv-analysis-in : {t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
// 783 lv-analysis-out: {t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  printf(" Time in seconds =             %12.2f\n",t);
// 783 lv-analysis-in : {mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
// 785 lv-analysis-out: {mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  printf(" Mop/s total     =             %12.2f\n",mops);
// 785 lv-analysis-in : {optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
// 787 lv-analysis-out: {optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  printf(" Operation type  = %24s\n",optype);
// 787 lv-analysis-in : {passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
// 789 lv-analysis-out: bot
  if (
// 790 lv-analysis-out: {passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
passed_verification
// 790 lv-analysis-in : {npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
) {
// 792 lv-analysis-out: {npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
    printf(" Verification    =               SUCCESSFUL\n");
// 792 lv-analysis-in : {npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  }
   else {
// 795 lv-analysis-out: {npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
    printf(" Verification    =             UNSUCCESSFUL\n");
// 795 lv-analysis-in : {npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  }
// 789 lv-analysis-in : bot
// 797 lv-analysis-out: {npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  printf(" Version         =             %12s\n",npbversion);
// 797 lv-analysis-in : {compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
// 799 lv-analysis-out: {compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  printf(" Compile date    =             %12s\n",compiletime);
// 799 lv-analysis-in : {cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
// 801 lv-analysis-out: {cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  printf("\n Compile options:\n");
// 801 lv-analysis-in : {cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
// 803 lv-analysis-out: {cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  printf("    CC           = %s\n",cc);
// 803 lv-analysis-in : {clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
// 805 lv-analysis-out: {clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  printf("    CLINK        = %s\n",clink);
// 805 lv-analysis-in : {c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
// 807 lv-analysis-out: {c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  printf("    C_LIB        = %s\n",c_lib);
// 807 lv-analysis-in : {c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
// 809 lv-analysis-out: {c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
  printf("    C_INC        = %s\n",c_inc);
// 809 lv-analysis-in : {cflags_300,clinkflags_301,rand_302,tmp_346}
// 811 lv-analysis-out: {cflags_300,clinkflags_301,rand_302,tmp_346}
  printf("    CFLAGS       = %s\n",cflags);
// 811 lv-analysis-in : {clinkflags_301,rand_302,tmp_346}
// 813 lv-analysis-out: {clinkflags_301,rand_302,tmp_346}
  printf("    CLINKFLAGS   = %s\n",clinkflags);
// 813 lv-analysis-in : {rand_302,tmp_346}
// 815 lv-analysis-out: {rand_302,tmp_346}
  printf("    RAND         = %s\n",rand);
// 815 lv-analysis-in : {tmp_346}
#ifdef SMP
#endif
/*    printf( "\n\n" );
    printf( " Please send the results of this run to:\n\n" );
    printf( " NPB Development Team\n" );
    printf( " Internet: npb@nas.nasa.gov\n \n" );
    printf( " If email is not available, send this to:\n\n" );
    printf( " MS T27A-1\n" );
    printf( " NASA Ames Research Center\n" );
    printf( " Moffett Field, CA  94035-1000\n\n" );
    printf( " Fax: 415-604-3957\n\n" );*/
}
// 763 lv-analysis-in : {name_283,cclass_284,n2_285,n3_286,n1_287,niter_288,nthreads_289,t_290,mops_291,optype_292,passed_verification_293,npbversion_294,compiletime_295,cc_296,clink_297,c_lib_298,c_inc_299,cflags_300,clinkflags_301,rand_302,tmp_346}
/*
cat ./common/c_timers.c
*/
/*
#include "wtime.h"
#if defined(IBM)
#define wtime wtime
#elif defined(CRAY)
#define wtime WTIME
#else
#define wtime wtime_
#endif
*/
/*  Prototype  */
void wtime(double *t);
/*****************************************************************/
/******         E  L  A  P  S  E  D  _  T  I  M  E          ******/
/*****************************************************************/

double elapsed_time()
// 817 lv-analysis-out: {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,start_304,n_310,tmp_346}
{
// 820 lv-analysis-out: {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,start_304,n_310,tmp_346}
  double t;
// 820 lv-analysis-in : {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310,tmp_346}
// 821 lv-analysis-out: {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310,tmp_346}
  wtime(&t);
// 821 lv-analysis-in : {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310,tmp_326,tmp_346}
// 823 lv-analysis-out: {nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310}
  return t;
// 823 lv-analysis-in : {nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,start_304,n_310}
}
// 817 lv-analysis-in : {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,start_304,n_310,tmp_346}
// 824 lv-analysis-out: bot
double start[64];
// 824 lv-analysis-in : bot
// 825 lv-analysis-out: bot
double elapsed[64];
// 825 lv-analysis-in : bot
/*****************************************************************/
/******            T  I  M  E  R  _  C  L  E  A  R          ******/
/*****************************************************************/

void timer_clear(int n)
// 826 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181,tmp_326}
{
// 829 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181,n_306}
  elapsed[n] = 0.0;
// 829 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181}
}
// 826 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181,n_306}
/*****************************************************************/
/******            T  I  M  E  R  _  S  T  A  R  T          ******/
/*****************************************************************/

void timer_start(int n)
// 830 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181,tmp_326}
{
// 833 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181,n_307}
  start[n] = elapsed_time();
// 833 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181,start_304}
}
// 830 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,a_159,x_160,z_161,p_162,q_163,r_164,w_165,nthreads_172,zeta_173,rnorm_174,cclass_179,zeta_verify_value_181,n_307}
/*****************************************************************/
/******            T  I  M  E  R  _  S  T  O  P             ******/
/*****************************************************************/

void timer_stop(int n)
// 834 lv-analysis-out: {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,start_304,tmp_326,tmp_346}
{
// 837 lv-analysis-out: {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,start_304,n_310,tmp_346}
  double t;
// 837 lv-analysis-in : {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,start_304,n_310,tmp_346}
// 838 lv-analysis-out: {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,start_304,n_310,tmp_346}
  double now;
// 838 lv-analysis-in : {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,start_304,n_310,tmp_346}
// 839 lv-analysis-out: {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,start_304,n_310,tmp_346}
  now = elapsed_time();
// 839 lv-analysis-in : {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,start_304,n_310,tmp_346}
// 841 lv-analysis-out: {nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,start_304,now_309,n_310}
  t = now - start[n];
// 841 lv-analysis-in : {nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_308,n_310}
// 842 lv-analysis-out: {nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_308,n_310}
  elapsed[n] += t;
// 842 lv-analysis-in : {nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,elapsed_305}
}
// 834 lv-analysis-in : {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,start_304,n_310,tmp_346}
/*****************************************************************/
/******            T  I  M  E  R  _  R  E  A  D             ******/
/*****************************************************************/

double timer_read(int n)
// 843 lv-analysis-out: {nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,elapsed_305,tmp_326}
{
// 846 lv-analysis-out: {nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,elapsed_305,n_311}
  return elapsed[n];
// 846 lv-analysis-in : {nthreads_172,zeta_173,cclass_179,zeta_verify_value_181}
}
// 843 lv-analysis-in : {nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,elapsed_305,n_311}

void wtime(double *t)
// 847 lv-analysis-out: {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310,tmp_326,tmp_346}
{
// 850 lv-analysis-out: {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310,t_314,tmp_346}
  static int sec = - 1;
// 850 lv-analysis-in : {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310,sec_312,t_314,tmp_346}
// 851 lv-analysis-out: {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310,sec_312,t_314,tmp_346}
  struct timeval tv;
// 851 lv-analysis-in : {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310,sec_312,tv_313,t_314,tmp_346}
// 852 lv-analysis-out: {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310,sec_312,tv_313,t_314,tmp_346}
  gettimeofday(&tv,((void *)0));
// 852 lv-analysis-in : {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310,sec_312,tv_313,t_314,tmp_346}
//  gettimeofday(&tv, (struct timezone *)0);
// 854 lv-analysis-out: bot
  if (
// 855 lv-analysis-out: {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310,sec_312,tv_313,t_314}
sec < 0
// 855 lv-analysis-in : {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310,sec_312,tv_313,t_314}
) {
// 857 lv-analysis-out: {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310,tv_313,t_314}
    sec = tv . tv_sec;
// 857 lv-analysis-in : {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310,sec_312,tv_313,t_314}
  }
// 854 lv-analysis-in : bot
// 858 lv-analysis-out: {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310,sec_312,tv_313,t_314}
   *t = (tv . tv_sec - sec) + 1.0e-6 * tv . tv_usec;
// 858 lv-analysis-in : {nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310}
}
// 847 lv-analysis-in : {tv_sec_72,tv_usec_73,nthreads_172,zeta_173,cclass_179,zeta_verify_value_181,t_303,start_304,n_310,t_314,tmp_346}
// common/c_randdp.c
/*
*/
#if defined(USE_POW)
#define r23 pow(0.5, 23.0)
#define r46 (r23*r23)
#define t23 pow(2.0, 23.0)
#define t46 (t23*t23)
#else
#define r23 (0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5)
#define r46 (r23*r23)
#define t23 (2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0)
#define t46 (t23*t23)
#endif
/*c---------------------------------------------------------------------
c---------------------------------------------------------------------*/

double randlc(double *x,double a)
// 859 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,tmp_326,tmp_327}
{
/*c---------------------------------------------------------------------
c---------------------------------------------------------------------*/
/*c---------------------------------------------------------------------
c
c   This routine returns a uniform pseudorandom double precision number in the
c   range (0, 1) by using the linear congruential generator
c
c   x_{k+1} = a x_k  (mod 2^46)
c
c   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
c   before repeating.  The argument A is the same as 'a' in the above formula,
c   and X is the same as x_0.  A and X must be odd double precision integers
c   in the range (1, 2^46).  The returned value RANDLC is normalized to be
c   between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
c   the new seed x_1, so that subsequent calls to RANDLC using the same
c   arguments will generate a continuous sequence.
c
c   This routine should produce the same results on any computer with at least
c   48 mantissa bits in double precision floating point data.  On 64 bit
c   systems, double precision should be disabled.
c
c   David H. Bailey     October 26, 1990
c
c---------------------------------------------------------------------*/
// 862 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a_324,x_325}
  double t1;
// 862 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a_324,x_325}
// 863 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a_324,x_325}
  double t2;
// 863 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a_324,x_325}
// 864 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a_324,x_325}
  double t3;
// 864 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a_324,x_325}
// 865 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a_324,x_325}
  double t4;
// 865 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a_324,x_325}
// 866 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a_324,x_325}
  double a1;
// 866 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a_324,x_325}
// 867 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a_324,x_325}
  double a2;
// 867 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a_324,x_325}
// 868 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a_324,x_325}
  double x1;
// 868 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a_324,x_325}
// 869 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a_324,x_325}
  double x2;
// 869 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,x2_322,a_324,x_325}
// 870 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,x2_322,a_324,x_325}
  double z;
// 870 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,x2_322,a_324,x_325}
/*c---------------------------------------------------------------------
c   Break A into two parts such that A = 2^23 * A1 + A2.
c---------------------------------------------------------------------*/
// 871 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,x2_322,a_324,x_325}
  t1 = 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * a;
// 871 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,t1_315,x2_322,a_324,x_325}
// 872 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,t1_315,x2_322,a_324,x_325}
  a1 = ((int )t1);
// 872 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,t1_315,a1_319,x2_322,a_324,x_325}
// 873 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,t1_315,a1_319,x2_322,a_324,x_325}
  a2 = a - 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * a1;
// 873 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,t1_315,a1_319,a2_320,x2_322,x_325}
/*c---------------------------------------------------------------------
c   Break X into two parts such that X = 2^23 * X1 + X2, compute
c   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
c   X = 2^23 * Z + A2 * X2  (mod 2^46).
c---------------------------------------------------------------------*/
// 874 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,t1_315,a1_319,a2_320,x2_322,x_325}
  t1 = 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 *  *x;
// 874 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,t1_315,a1_319,a2_320,x2_322,x_325}
// 875 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,t1_315,a1_319,a2_320,x2_322,x_325}
  x1 = ((int )t1);
// 875 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a1_319,a2_320,x1_321,x2_322,x_325}
// 876 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a1_319,a2_320,x1_321,x2_322,x_325}
  x2 =  *x - 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * x1;
// 876 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a1_319,a2_320,x1_321,x2_322,x_325}
// 877 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a1_319,a2_320,x1_321,x2_322,x_325}
  t1 = a1 * x2 + a2 * x1;
// 877 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,t1_315,a2_320,x2_322,x_325}
// 878 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,t1_315,a2_320,x2_322,x_325}
  t2 = ((int )(0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * t1));
// 878 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,t1_315,t2_316,a2_320,x2_322,x_325}
// 879 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,t1_315,t2_316,a2_320,x2_322,x_325}
  z = t1 - 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * t2;
// 879 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a2_320,x2_322,z_323,x_325}
// 880 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a2_320,x2_322,z_323,x_325}
  t3 = 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * z + a2 * x2;
// 880 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,t3_317,x_325}
// 881 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,t3_317,x_325}
  t4 = ((int )(0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * (0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5) * t3));
// 881 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,t3_317,t4_318,x_325}
// 882 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,t3_317,t4_318,x_325}
   *x = t3 - 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * (2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0) * t4;
// 882 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,x_325}
// 883 lv-analysis-out: {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,x_325}
  return 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * (0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5) *  *x;
// 883 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304}
}
// 859 lv-analysis-in : {tv_sec_72,tv_usec_73,naa_146,nzz_147,firstrow_148,lastrow_149,firstcol_150,lastcol_151,colidx_152,rowstr_153,iv_154,arow_155,acol_156,v_157,aelt_158,a_159,x_160,z_161,p_162,q_163,r_164,w_165,amult_166,tran_167,nthreads_172,rnorm_174,cclass_179,zeta_verify_value_181,nnza_207,iouter_208,nzv_212,size_213,ratio_214,n_218,rcond_219,colidx_220,nonzer_221,v_224,iv_225,firstcol_226,lastcol_227,firstrow_228,lastrow_229,nz_230,acol_231,arow_232,aelt_233,shift_234,a_237,rowstr_238,nzrow_261,nzv_262,vecelt_265,n_267,nz_268,mark_269,nzloc_270,start_304,a_324,x_325}
