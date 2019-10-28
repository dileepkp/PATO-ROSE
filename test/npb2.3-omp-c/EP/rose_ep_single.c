/*--------------------------------------------------------------------
  
  NAS Parallel Benchmarks 2.3 OpenMP C versions - EP
  This benchmark is an OpenMP C version of the NPB EP code.
  
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
  Author: P. O. Frederickson 
          D. H. Bailey
          A. C. Woo
  OpenMP C version: S. Satoh
  
--------------------------------------------------------------------*/
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
void vranlc(int n,double *x_seed,double a,double *y);
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
#define CLASS    'S'
#define M       24
#define CONVERTDOUBLE   FALSE
#endif
#if CLASS == 'W'
/* CLASS = W */
/*
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the classc of the NPB
c  in this directory. Do not modify it by hand.
*/
#define CLASS    'W'
#define M       25
#define CONVERTDOUBLE   FALSE
#endif
#if CLASS == 'A'
/* CLASS = A */
/*
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the classc of the NPB
c  in this directory. Do not modify it by hand.
*/
#define CLASS    'A'
#define M       28
#define CONVERTDOUBLE   FALSE
#endif
#if CLASS == 'B'
/* CLASS = B */
/*
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the classc of the NPB
c  in this directory. Do not modify it by hand.
*/
#define CLASS    'B'
#define M       30
#define CONVERTDOUBLE   FALSE
#endif
#if CLASS == 'C'
/* CLASS = C */
/*
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the classc of the NPB
c  in this directory. Do not modify it by hand.
*/
#define CLASS    'C'
#define M       32
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
/* parameters */
#define	MK		16
#define	MM		(M - MK)
#define	NN		(1 << MM)
#define	NK		(1 << MK)
#define	NQ		10
#define EPSILON		1.0e-8
#define	A		1220703125.0
#define	S		271828183.0
#define	TIMERS_ENABLED	FALSE
/* global variables */
/* common /storage/ */
// 116 lv-analysis-out: bot
static double x[2 * (1 << 16)];
// 116 lv-analysis-in : bot
// 117 lv-analysis-out: bot
#pragma omp threadprivate(x)
// 117 lv-analysis-in : bot
// 118 lv-analysis-out: bot
static double q[10];
// 118 lv-analysis-in : bot
/*--------------------------------------------------------------------
      program EMBAR
c-------------------------------------------------------------------*/
/*
c   This is the serial version of the APP Benchmark 1,
c   the "embarassingly parallel" benchmark.
c
c   M is the Log_2 of the number of complex pairs of uniform (0, 1) random
c   numbers.  MK is the Log_2 of the size of each batch of uniform random
c   numbers.  MK can be set for convenience on a given system, since it does
c   not affect the results.
*/

int main(int argc,char **argv)
// 119 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_279,tmp_280,tmp_299}
{
// 122 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  double Mops;
// 122 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 123 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  double t1;
// 123 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 124 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  double t2;
// 124 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 125 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  double t3;
// 125 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 126 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  double t4;
// 126 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 127 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  double x1;
// 127 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 128 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  double x2;
// 128 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 129 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  double sx;
// 129 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 130 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  double sy;
// 130 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 131 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  double tm;
// 131 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 132 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  double an;
// 132 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 133 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  double tt;
// 133 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 134 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  double gc;
// 134 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 135 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  double dum[3] = {(1.0), (1.0), (1.0)};
// 135 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 136 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  int np;
// 136 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 137 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  int ierr;
// 137 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 138 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  int node;
// 138 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 139 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  int no_nodes;
// 139 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 140 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  int i;
// 140 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 141 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  int ik;
// 141 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 142 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  int kk;
// 142 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 143 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  int l;
// 143 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 144 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  int k;
// 144 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 145 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  int nit;
// 145 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 146 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  int ierrcode;
// 146 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 147 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  int no_large_nodes;
// 147 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 148 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  int np_add;
// 148 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 149 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  int k_offset;
// 149 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 150 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  int j;
// 150 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 151 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  int nthreads = 1;
// 151 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 152 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  boolean verified;
// 152 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
/* character*13 */
// 153 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  char size[13 + 1];
// 153 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,size_179,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
/*
c   Because the size of the problem is too large to store in a 32-bit
c   integer for some classes, we put it into a string (for printing).
c   Have to strip off the decimal point put in there by the floating
c   point print statement (internal file)
*/
// 154 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,size_179,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  printf("\n\n NAS Parallel Benchmarks 2.3 OpenMP C version - EP Benchmark\n");
// 154 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,size_179,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 156 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,size_179,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  double __temp0__ = (double )(30 + 1);
// 156 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,size_179,__temp0___180,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 157 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,size_179,__temp0___180,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  double __temp1__ = pow(2.0,__temp0__);
// 157 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,size_179,__temp1___181,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 158 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,size_179,__temp1___181,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  sprintf(size,"%12.0f",__temp1__);
// 158 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,size_179,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 160 lv-analysis-out: bot
  for (
// 161 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,size_179,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
j = 13
// 161 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,j_176,nthreads_177,size_179,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
; 
// 162 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,j_176,nthreads_177,size_179,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
j >= 1;
// 162 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,j_176,nthreads_177,size_179,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
 j--) {
// 165 lv-analysis-out: bot
    if (
// 166 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,j_176,nthreads_177,size_179,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
size[j] == '.'
// 166 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,j_176,nthreads_177,size_179,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
) {
// 168 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,j_176,nthreads_177,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
      size[j] = 32;
// 168 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,j_176,nthreads_177,size_179,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
    }
// 165 lv-analysis-in : bot
  }
// 160 lv-analysis-in : bot
// 169 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,size_179,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  printf(" Number of random numbers generated: %13s\n",size);
// 169 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 171 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
  verified = 0;
// 171 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
/*
c   Compute the number of "batches" of random number pairs generated 
c   per processor. Adjust if the number of processors does not evenly 
c   divide the total number
*/
// 172 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
  np = 1 << 30 - 16;
// 172 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
/*
c   Call the random number generator functions and initialize
c   the x-array to reduce the effects of paging on the timings.
c   Also, call all mathematical functions that are used. Make
c   sure these initializations cannot be eliminated as dead code.
*/
// 173 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
  double *__temp2__ = &dum[0];
// 173 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,__temp2___182,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
// 174 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,__temp2___182,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
  double __temp3__ = dum[1];
// 174 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,__temp2___182,__temp3___183,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
// 175 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,__temp2___182,__temp3___183,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
  double *__temp4__ = &dum[2];
// 175 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,__temp2___182,__temp3___183,__temp4___184,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
// 176 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,__temp2___182,__temp3___183,__temp4___184,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
  vranlc(0,__temp2__,__temp3__,__temp4__);
// 176 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_279,tmp_280,tmp_281,tmp_282}
// 178 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
  double *__temp5__ = &dum[1];
// 178 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
// 179 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
  double __temp6__ = dum[2];
// 179 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
// 180 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
  dum[0] = randlc(__temp5__,__temp6__);
// 180 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
// 181 lv-analysis-out: bot
  for (
// 182 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
i = 0
// 182 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
; 
// 183 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
i < 2 * (1 << 16);
// 183 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
 i++) {
// 186 lv-analysis-out: {tv_sec_72,tv_usec_73,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
    x[i] = - 1.0e99;
// 186 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
  }
// 181 lv-analysis-in : bot
// 187 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
  double __temp7__ = 1.0 > 1.0?1.0 : 1.0;
// 187 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,__temp7___187,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
// 188 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,__temp7___187,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
  double __temp8__ = fabs(__temp7__);
// 188 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,__temp8___188,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
// 189 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,__temp8___188,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
  double __temp9__ = sqrt(__temp8__);
// 189 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,__temp9___189,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
// 190 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,__temp9___189,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
  Mops = log(__temp9__);
// 190 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_299}
// 192 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200}
  timer_clear(1);
// 192 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_279}
// 194 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200}
  timer_clear(2);
// 194 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_279}
// 196 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200}
  timer_clear(3);
// 196 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_279}
// 198 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200}
  timer_start(1);
// 198 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_279}
// 200 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
  vranlc(0,&t1,1220703125.0,x);
// 200 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_279,tmp_280,tmp_281,tmp_282}
/*   Compute AN = A ^ (2 * NK) (mod 2^46). */
// 202 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  t1 = 1220703125.0;
// 202 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
// 203 lv-analysis-out: bot
  for (
// 204 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
i = 1
// 204 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
; 
// 205 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
i <= 16 + 1;
// 205 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
 i++) {
// 208 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
    t2 = randlc(&t1,t1);
// 208 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_279,tmp_280}
  }
// 203 lv-analysis-in : bot
// 210 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,dum_161,np_162,nthreads_177,verified_178,start_242,tmp_299}
  an = t1;
// 210 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,an_158,dum_161,np_162,nthreads_177,verified_178,start_242,tmp_299}
// 211 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,an_158,dum_161,np_162,nthreads_177,verified_178,start_242,tmp_299}
  tt = 271828183.0;
// 211 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,an_158,dum_161,np_162,nthreads_177,verified_178,start_242,tmp_299}
// 212 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,an_158,dum_161,np_162,nthreads_177,verified_178,start_242,tmp_299}
  gc = 0.0;
// 212 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,an_158,gc_160,dum_161,np_162,nthreads_177,verified_178,start_242,tmp_299}
// 213 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,an_158,gc_160,dum_161,np_162,nthreads_177,verified_178,start_242,tmp_299}
  sx = 0.0;
// 213 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,an_158,gc_160,dum_161,np_162,nthreads_177,verified_178,start_242,tmp_299}
// 214 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,an_158,gc_160,dum_161,np_162,nthreads_177,verified_178,start_242,tmp_299}
  sy = 0.0;
// 214 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,nthreads_177,verified_178,start_242,tmp_299}
// 215 lv-analysis-out: bot
  for (
// 216 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,nthreads_177,verified_178,start_242,tmp_299}
i = 0
// 216 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,nthreads_177,verified_178,start_242,tmp_299}
; 
// 217 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,nthreads_177,verified_178,start_242,tmp_299}
i <= 10 - 1;
// 217 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,nthreads_177,verified_178,start_242,tmp_299}
 i++) {
// 220 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,nthreads_177,verified_178,start_242,tmp_299}
    q[i] = 0.0;
// 220 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,nthreads_177,verified_178,start_242,tmp_299}
  }
// 215 lv-analysis-in : bot
/*
c   Each instance of this loop may be performed independently. We compute
c   the k offsets separately to take into account the fact that some nodes
c   have more numbers to generate than others
*/
// 221 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,nthreads_177,verified_178,start_242,tmp_299}
  k_offset = - 1;
// 221 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,start_242,tmp_299}
// 222 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,start_242,tmp_299}
  
#pragma omp parallel copyin(x)
// 222 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,start_242,tmp_299}
{
// 224 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,start_242,tmp_299}
    double t1;
// 224 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,start_242,tmp_299}
// 225 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,start_242,tmp_299}
    double t2;
// 225 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,start_242,tmp_299}
// 226 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,start_242,tmp_299}
    double t3;
// 226 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,start_242,tmp_299}
// 227 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,start_242,tmp_299}
    double t4;
// 227 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,start_242,tmp_299}
// 228 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,start_242,tmp_299}
    double x1;
// 228 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,start_242,tmp_299}
// 229 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,start_242,tmp_299}
    double x2;
// 229 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,start_242,tmp_299}
// 230 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,start_242,tmp_299}
    int kk;
// 230 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,start_242,tmp_299}
// 231 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,start_242,tmp_299}
    int i;
// 231 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,start_242,tmp_299}
// 232 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,start_242,tmp_299}
    int ik;
// 232 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,ik_198,start_242,tmp_299}
// 233 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,ik_198,start_242,tmp_299}
    int l;
// 233 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,ik_198,start_242,tmp_299}
/* private copy of q[0:NQ-1] */
// 234 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,ik_198,start_242,tmp_299}
    double qq[10];
// 234 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,ik_198,qq_200,start_242,tmp_299}
// 235 lv-analysis-out: bot
    for (
// 236 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,ik_198,qq_200,start_242,tmp_299}
i = 0
// 236 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
; 
// 237 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
i < 10;
// 237 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
 i++) {
// 240 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,start_242,tmp_299}
      qq[i] = 0.0;
// 240 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
    }
// 235 lv-analysis-in : bot
// 241 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,ik_198,qq_200,start_242,tmp_299}
    
#pragma omp for reduction(+:sx,sy) schedule(static)
// 241 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,ik_198,qq_200,start_242,tmp_299}
// 242 lv-analysis-out: bot
    for (
// 243 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,ik_198,qq_200,start_242,tmp_299}
k = 1
// 243 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,ik_198,qq_200,start_242,tmp_299}
; 
// 244 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,ik_198,qq_200,start_242,tmp_299}
k <= np;
// 244 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,ik_198,qq_200,start_242,tmp_299}
 k++) {
// 247 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,ik_198,qq_200,start_242}
      kk = k_offset + k;
// 247 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,kk_196,ik_198,qq_200,start_242}
// 248 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,kk_196,ik_198,qq_200,start_242}
      t1 = 271828183.0;
// 248 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,kk_196,ik_198,qq_200,start_242}
// 249 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,kk_196,ik_198,qq_200,start_242}
      t2 = an;
// 249 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,kk_196,ik_198,qq_200,start_242}
/*      Find starting seed t1 for this kk. */
// 250 lv-analysis-out: bot
      for (
// 251 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,kk_196,ik_198,qq_200,start_242}
i = 1
// 251 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,kk_196,i_197,ik_198,qq_200,start_242}
; 
// 252 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,kk_196,i_197,ik_198,qq_200,start_242}
i <= 100;
// 252 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,kk_196,i_197,ik_198,qq_200,start_242}
 i++) {
// 255 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,kk_196,i_197,qq_200,start_242}
        ik = kk / 2;
// 255 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,kk_196,i_197,ik_198,qq_200,start_242}
// 256 lv-analysis-out: bot
        if (
// 257 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,kk_196,i_197,ik_198,qq_200,start_242}
2 * ik != kk
// 257 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
) {
// 259 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
          t3 = randlc(&t1,t2);
// 259 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_279,tmp_280}
        }
// 256 lv-analysis-in : bot
// 261 lv-analysis-out: bot
        if (
// 262 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
ik == 0
// 262 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
) {
// 264 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
          break; 
// 264 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
        }
// 261 lv-analysis-in : bot
// 265 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
        t3 = randlc(&t2,t2);
// 265 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_279,tmp_280}
// 267 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
        kk = ik;
// 267 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,kk_196,i_197,ik_198,qq_200,start_242}
      }
// 250 lv-analysis-in : bot
/*      Compute uniform pseudorandom numbers. */
// 268 lv-analysis-out: bot
      if (
// 269 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
0 == 1
// 269 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
) {
// 271 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200}
        timer_start(3);
// 271 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_279}
      }
// 268 lv-analysis-in : bot
// 273 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
      int __temp10__ = 2 * (1 << 16);
// 273 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,__temp10___201,start_242}
// 274 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,__temp10___201,start_242}
      double *__temp11__ = x - 1;
// 274 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,__temp10___201,__temp11___202,start_242}
// 275 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,__temp10___201,__temp11___202,start_242}
      vranlc(__temp10__,&t1,1220703125.0,__temp11__);
// 275 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_279,tmp_280,tmp_281,tmp_282}
// 277 lv-analysis-out: bot
      if (
// 278 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
0 == 1
// 278 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
) {
// 280 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
        timer_stop(3);
// 280 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_279,tmp_299}
      }
// 277 lv-analysis-in : bot
/*
c       Compute Gaussian deviates by acceptance-rejection method and 
c       tally counts in concentric square annuli.  This loop is not 
c       vectorizable.
*/
// 282 lv-analysis-out: bot
      if (
// 283 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
0 == 1
// 283 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
) {
// 285 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200}
        timer_start(2);
// 285 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_279}
      }
// 282 lv-analysis-in : bot
// 287 lv-analysis-out: bot
      for (
// 288 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,ik_198,qq_200,start_242,tmp_299}
i = 0
// 288 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
; 
// 289 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
i < 1 << 16;
// 289 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
 i++) {
// 292 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
        x1 = 2.0 * x[2 * i] - 1.0;
// 292 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t2_191,x1_194,i_197,ik_198,qq_200,start_242,tmp_299}
// 293 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t2_191,x1_194,i_197,ik_198,qq_200,start_242,tmp_299}
        x2 = 2.0 * x[2 * i + 1] - 1.0;
// 293 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t2_191,x1_194,x2_195,i_197,ik_198,qq_200,start_242,tmp_299}
// 294 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t2_191,x1_194,x2_195,i_197,ik_198,qq_200,start_242,tmp_299}
        t1 = x1 * x1 + x2 * x2;
// 294 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,x1_194,x2_195,i_197,ik_198,qq_200,start_242,tmp_299}
// 295 lv-analysis-out: bot
        if (
// 296 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,x1_194,x2_195,i_197,ik_198,qq_200,start_242,tmp_299}
t1 <= 1.0
// 296 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,x1_194,x2_195,i_197,ik_198,qq_200,start_242,tmp_299}
) {
// 298 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,x1_194,x2_195,i_197,ik_198,start_242,tmp_299}
          double __temp12__ = - 2.0 * log(t1) / t1;
// 298 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,x1_194,x2_195,i_197,ik_198,__temp12___203,start_242,tmp_299}
// 299 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,x1_194,x2_195,i_197,ik_198,__temp12___203,start_242,tmp_299}
          t2 = sqrt(__temp12__);
// 299 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,x1_194,x2_195,i_197,ik_198,start_242,tmp_299}
/* Xi */
// 301 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,x1_194,x2_195,i_197,ik_198,start_242,tmp_299}
          t3 = x1 * t2;
// 301 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,t3_192,x2_195,i_197,ik_198,start_242,tmp_299}
/* Yi */
// 302 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,t3_192,x2_195,i_197,ik_198,start_242,tmp_299}
          t4 = x2 * t2;
// 302 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,t3_192,t4_193,i_197,ik_198,start_242,tmp_299}
// 303 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,t3_192,t4_193,i_197,ik_198,start_242,tmp_299}
          l = ((fabs(t3) > fabs(t4)?fabs(t3) : fabs(t4)));
// 303 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,t3_192,t4_193,i_197,ik_198,l_199,start_242,tmp_299}
/* counts */
// 304 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,t3_192,t4_193,i_197,ik_198,l_199,start_242,tmp_299}
          qq[l] += 1.0;
// 304 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,t3_192,t4_193,i_197,ik_198,qq_200,start_242,tmp_299}
/* sum of Xi */
// 305 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,t3_192,t4_193,i_197,ik_198,qq_200,start_242,tmp_299}
          sx = sx + t3;
// 305 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,t4_193,i_197,ik_198,qq_200,start_242,tmp_299}
/* sum of Yi */
// 306 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,t4_193,i_197,ik_198,qq_200,start_242,tmp_299}
          sy = sy + t4;
// 306 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
        }
// 295 lv-analysis-in : bot
      }
// 287 lv-analysis-in : bot
// 307 lv-analysis-out: bot
      if (
// 308 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
0 == 1
// 308 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
) {
// 310 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
        timer_stop(2);
// 310 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_279,tmp_299}
      }
// 307 lv-analysis-in : bot
    }
// 242 lv-analysis-in : bot
// 312 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,ik_198,qq_200,start_242,tmp_299}
    
#pragma omp critical
// 312 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,ik_198,qq_200,start_242,tmp_299}
{
// 314 lv-analysis-out: bot
      for (
// 315 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,ik_198,qq_200,start_242,tmp_299}
i = 0
// 315 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
; 
// 316 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
i <= 10 - 1;
// 316 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
 i++) {
// 319 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
        q[i] += qq[i];
// 319 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
      }
// 314 lv-analysis-in : bot
    }
#if defined(_OPENMP)
#endif /* _OPENMP */    
/* end of parallel region */
  }
// 320 lv-analysis-out: bot
  for (
// 321 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
i = 0
// 321 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
; 
// 322 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
i <= 10 - 1;
// 322 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
 i++) {
// 325 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
    gc = gc + q[i];
// 325 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  }
// 320 lv-analysis-in : bot
// 326 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
  timer_stop(1);
// 326 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_279,tmp_299}
// 328 lv-analysis-out: {q_147,sx_155,sy_156,gc_160,nthreads_177,verified_178,elapsed_243}
  tm = timer_read(1);
// 328 lv-analysis-in : {q_147,sx_155,sy_156,gc_160,nthreads_177,verified_178,elapsed_243,tmp_279}
// 330 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nthreads_177,verified_178,tmp_299}
  nit = 0;
// 330 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
// 331 lv-analysis-out: bot
  if (
// 332 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
30 == 24
// 332 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
) {
// 334 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
    double __temp13__ = (sx - - 3.247834652034740e3) / sx;
// 334 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
// 335 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
    double __temp24__ = (sy - - 6.958407078382297e3) / sy;
// 335 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
// 336 lv-analysis-out: bot
    if (
// 337 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
fabs(__temp13__) <= 1.0e-8 && fabs(__temp24__) <= 1.0e-8
// 337 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
) {
// 339 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,tmp_299}
      verified = 1;
// 339 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
    }
// 336 lv-analysis-in : bot
  }
   else {
// 341 lv-analysis-out: bot
    if (
// 342 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
30 == 25
// 342 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
) {
// 344 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
      double __temp14__ = (sx - - 2.863319731645753e3) / sx;
// 344 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
// 345 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
      double __temp25__ = (sy - - 6.320053679109499e3) / sy;
// 345 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
// 346 lv-analysis-out: bot
      if (
// 347 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
fabs(__temp14__) <= 1.0e-8 && fabs(__temp25__) <= 1.0e-8
// 347 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
) {
// 349 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,tmp_299}
        verified = 1;
// 349 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
      }
// 346 lv-analysis-in : bot
    }
     else {
// 351 lv-analysis-out: bot
      if (
// 352 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
30 == 28
// 352 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
) {{
// 355 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
          double __temp15__ = (sx - - 4.295875165629892e3) / sx;
// 355 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
// 356 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
          double __temp26__ = (sy - - 1.580732573678431e4) / sy;
// 356 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
// 357 lv-analysis-out: bot
          if (
// 358 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
fabs(__temp15__) <= 1.0e-8 && fabs(__temp26__) <= 1.0e-8
// 358 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
) {
// 360 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,tmp_299}
            verified = 1;
// 360 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
          }
// 357 lv-analysis-in : bot
// 361 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
          printf("Debug: 231, sx is:%f, sy is:%f\n",sx,sy);
// 361 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
        }
      }
       else {
// 364 lv-analysis-out: bot
        if (
// 365 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
30 == 30
// 365 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
) {
// 367 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
          double __temp16__ = (sx - 4.033815542441498e4) / sx;
// 367 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
// 368 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
          double __temp27__ = (sy - - 2.660669192809235e4) / sy;
// 368 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
// 369 lv-analysis-out: bot
          if (
// 370 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
fabs(__temp16__) <= 1.0e-8 && fabs(__temp27__) <= 1.0e-8
// 370 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
) {
// 372 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,tmp_299}
            verified = 1;
// 372 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
          }
// 369 lv-analysis-in : bot
        }
         else {
// 374 lv-analysis-out: bot
          if (
// 375 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
30 == 32
// 375 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
) {
// 377 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
            double __temp17__ = (sx - 4.764367927995374e4) / sx;
// 377 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
// 378 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
            double __temp28__ = (sy - - 8.084072988043731e4) / sy;
// 378 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
// 379 lv-analysis-out: bot
            if (
// 380 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
fabs(__temp17__) <= 1.0e-8 && fabs(__temp28__) <= 1.0e-8
// 380 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
) {
// 382 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,tmp_299}
              verified = 1;
// 382 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
            }
// 379 lv-analysis-in : bot
          }
// 374 lv-analysis-in : bot
        }
// 364 lv-analysis-in : bot
      }
// 351 lv-analysis-in : bot
    }
// 341 lv-analysis-in : bot
  }
// 331 lv-analysis-in : bot
// 383 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
  double __temp18__ = (double )(30 + 1);
// 383 lv-analysis-in : {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
// 384 lv-analysis-out: {q_147,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
  Mops = pow(2.0,__temp18__) / tm / 1000000.0;
// 384 lv-analysis-in : {q_147,Mops_148,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
// 385 lv-analysis-out: {q_147,Mops_148,sx_155,sy_156,tm_157,gc_160,nit_171,nthreads_177,verified_178,tmp_299}
  printf("EP Benchmark Results: \nCPU Time = %10.4f\nN = 2^%5d\nNo. Gaussian Pairs = %15.0f\nSums = %25.15e %25.15e\nCounts:\n",tm,30,gc,sx,sy);
// 385 lv-analysis-in : {q_147,Mops_148,tm_157,nit_171,nthreads_177,verified_178,tmp_299}
// 387 lv-analysis-out: bot
  for (
// 388 lv-analysis-out: {q_147,Mops_148,tm_157,nit_171,nthreads_177,verified_178,tmp_299}
i = 0
// 388 lv-analysis-in : {q_147,Mops_148,tm_157,i_166,nit_171,nthreads_177,verified_178,tmp_299}
; 
// 389 lv-analysis-out: {q_147,Mops_148,tm_157,i_166,nit_171,nthreads_177,verified_178,tmp_299}
i <= 10 - 1;
// 389 lv-analysis-in : {q_147,Mops_148,tm_157,i_166,nit_171,nthreads_177,verified_178,tmp_299}
 i++) {
// 392 lv-analysis-out: {q_147,Mops_148,tm_157,i_166,nit_171,nthreads_177,verified_178,tmp_299}
    double __temp19__ = q[i];
// 392 lv-analysis-in : {q_147,Mops_148,tm_157,i_166,nit_171,nthreads_177,verified_178,__temp19___215,tmp_299}
// 393 lv-analysis-out: {q_147,Mops_148,tm_157,i_166,nit_171,nthreads_177,verified_178,__temp19___215,tmp_299}
    printf("%3d %15.0f\n",i,__temp19__);
// 393 lv-analysis-in : {q_147,Mops_148,tm_157,i_166,nit_171,nthreads_177,verified_178,tmp_299}
  }
// 387 lv-analysis-in : bot
// 395 lv-analysis-out: {Mops_148,tm_157,nit_171,nthreads_177,verified_178,tmp_299}
  int __temp20__ = 30 + 1;
// 395 lv-analysis-in : {Mops_148,tm_157,nit_171,nthreads_177,verified_178,__temp20___216,tmp_299}
// 396 lv-analysis-out: {Mops_148,tm_157,nit_171,nthreads_177,verified_178,__temp20___216,tmp_299}
  c_print_results("EP",'B',__temp20__,0,0,nit,nthreads,tm,Mops,"Random numbers generated",verified,"2.3","28 Oct 2014","gcc","$(CC)","(none)","-I../common","-fopenmp -O2","-lm -fopenmp","randdp");
// 396 lv-analysis-in : {tmp_279,tmp_280,tmp_281,tmp_282,tmp_283,tmp_284,tmp_285,tmp_286,tmp_287,tmp_288,tmp_289,tmp_290,tmp_291,tmp_292,tmp_293,tmp_294,tmp_295,tmp_296,tmp_297,tmp_298,tmp_299}
// 398 lv-analysis-out: bot
  if (
// 399 lv-analysis-out: {tmp_299}
0 == 1
// 399 lv-analysis-in : {tmp_299}
) {
// 401 lv-analysis-out: {tmp_299}
    double __temp21__ = timer_read(1);
// 401 lv-analysis-in : {__temp21___217,tmp_299}
// 402 lv-analysis-out: {__temp21___217,tmp_299}
    printf("Total time:     %f",__temp21__);
// 402 lv-analysis-in : {tmp_299}
// 404 lv-analysis-out: {tmp_299}
    double __temp22__ = timer_read(2);
// 404 lv-analysis-in : {__temp22___218,tmp_299}
// 405 lv-analysis-out: {__temp22___218,tmp_299}
    printf("Gaussian pairs: %f",__temp22__);
// 405 lv-analysis-in : {tmp_299}
// 407 lv-analysis-out: {tmp_299}
    double __temp23__ = timer_read(3);
// 407 lv-analysis-in : {__temp23___219,tmp_299}
// 408 lv-analysis-out: {__temp23___219,tmp_299}
    printf("Random numbers: %f",__temp23__);
// 408 lv-analysis-in : {tmp_299}
  }
// 398 lv-analysis-in : bot
}
// 119 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_299}
/* cat ./common/c_print_results.c */
/*****************************************************************/
/******     C  _  P  R  I  N  T  _  R  E  S  U  L  T  S     ******/
/*****************************************************************/

void c_print_results(char *name,char cclass,int n1,int n2,int n3,int niter,int nthreads,double t,double mops,char *optype,int passed_verification,char *npbversion,char *compiletime,char *cc,char *clink,char *c_lib,char *c_inc,char *cflags,char *clinkflags,char *rand)
// 410 lv-analysis-out: {tmp_279,tmp_280,tmp_281,tmp_282,tmp_283,tmp_284,tmp_285,tmp_286,tmp_287,tmp_288,tmp_289,tmp_290,tmp_291,tmp_292,tmp_293,tmp_294,tmp_295,tmp_296,tmp_297,tmp_298,tmp_299}
{
// 413 lv-analysis-out: {name_221,cclass_222,n2_223,n3_224,n1_225,niter_226,nthreads_227,t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  char *evalue = "1000";
// 413 lv-analysis-in : {name_221,cclass_222,n2_223,n3_224,n1_225,niter_226,nthreads_227,t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
// 414 lv-analysis-out: {name_221,cclass_222,n2_223,n3_224,n1_225,niter_226,nthreads_227,t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  printf("\n\n %s Benchmark Completed\n",name);
// 414 lv-analysis-in : {cclass_222,n2_223,n3_224,n1_225,niter_226,nthreads_227,t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
// 416 lv-analysis-out: {cclass_222,n2_223,n3_224,n1_225,niter_226,nthreads_227,t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  printf(" Class           =                        %c\n",cclass);
// 416 lv-analysis-in : {n2_223,n3_224,n1_225,niter_226,nthreads_227,t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
// 418 lv-analysis-out: bot
  if (
// 419 lv-analysis-out: {n2_223,n3_224,n1_225,niter_226,nthreads_227,t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
n2 == 0 && n3 == 0
// 419 lv-analysis-in : {n2_223,n3_224,n1_225,niter_226,nthreads_227,t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
) {
/* as in IS */
// 421 lv-analysis-out: {n1_225,niter_226,nthreads_227,t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
    printf(" Size            =             %12d\n",n1);
// 421 lv-analysis-in : {niter_226,nthreads_227,t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  }
   else {
// 424 lv-analysis-out: {n2_223,n3_224,n1_225,niter_226,nthreads_227,t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
    printf(" Size            =              %3dx%3dx%3d\n",n1,n2,n3);
// 424 lv-analysis-in : {niter_226,nthreads_227,t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  }
// 418 lv-analysis-in : bot
// 426 lv-analysis-out: {niter_226,nthreads_227,t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  printf(" Iterations      =             %12d\n",niter);
// 426 lv-analysis-in : {nthreads_227,t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
// 428 lv-analysis-out: {nthreads_227,t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  printf(" Threads         =             %12d\n",nthreads);
// 428 lv-analysis-in : {t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
// 430 lv-analysis-out: {t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  printf(" Time in seconds =             %12.2f\n",t);
// 430 lv-analysis-in : {mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
// 432 lv-analysis-out: {mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  printf(" Mop/s total     =             %12.2f\n",mops);
// 432 lv-analysis-in : {optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
// 434 lv-analysis-out: {optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  printf(" Operation type  = %24s\n",optype);
// 434 lv-analysis-in : {passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
// 436 lv-analysis-out: bot
  if (
// 437 lv-analysis-out: {passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
passed_verification
// 437 lv-analysis-in : {npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
) {
// 439 lv-analysis-out: {npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
    printf(" Verification    =               SUCCESSFUL\n");
// 439 lv-analysis-in : {npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  }
   else {
// 442 lv-analysis-out: {npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
    printf(" Verification    =             UNSUCCESSFUL\n");
// 442 lv-analysis-in : {npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  }
// 436 lv-analysis-in : bot
// 444 lv-analysis-out: {npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  printf(" Version         =             %12s\n",npbversion);
// 444 lv-analysis-in : {compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
// 446 lv-analysis-out: {compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  printf(" Compile date    =             %12s\n",compiletime);
// 446 lv-analysis-in : {cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
// 448 lv-analysis-out: {cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  printf("\n Compile options:\n");
// 448 lv-analysis-in : {cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
// 450 lv-analysis-out: {cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  printf("    CC           = %s\n",cc);
// 450 lv-analysis-in : {clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
// 452 lv-analysis-out: {clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  printf("    CLINK        = %s\n",clink);
// 452 lv-analysis-in : {c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
// 454 lv-analysis-out: {c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  printf("    C_LIB        = %s\n",c_lib);
// 454 lv-analysis-in : {c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
// 456 lv-analysis-out: {c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
  printf("    C_INC        = %s\n",c_inc);
// 456 lv-analysis-in : {cflags_238,clinkflags_239,rand_240,tmp_299}
// 458 lv-analysis-out: {cflags_238,clinkflags_239,rand_240,tmp_299}
  printf("    CFLAGS       = %s\n",cflags);
// 458 lv-analysis-in : {clinkflags_239,rand_240,tmp_299}
// 460 lv-analysis-out: {clinkflags_239,rand_240,tmp_299}
  printf("    CLINKFLAGS   = %s\n",clinkflags);
// 460 lv-analysis-in : {rand_240,tmp_299}
// 462 lv-analysis-out: {rand_240,tmp_299}
  printf("    RAND         = %s\n",rand);
// 462 lv-analysis-in : {tmp_299}
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
// 410 lv-analysis-in : {name_221,cclass_222,n2_223,n3_224,n1_225,niter_226,nthreads_227,t_228,mops_229,optype_230,passed_verification_231,npbversion_232,compiletime_233,cc_234,clink_235,c_lib_236,c_inc_237,cflags_238,clinkflags_239,rand_240,tmp_299}
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
// 464 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,n_248,tmp_299}
{
// 467 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,n_248,tmp_299}
  double t;
// 467 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248,tmp_299}
// 468 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248,tmp_299}
  wtime(&t);
// 468 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248,tmp_279,tmp_299}
// 470 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248}
  return t;
// 470 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,n_248}
}
// 464 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,n_248,tmp_299}
// 471 lv-analysis-out: bot
double start[64];
// 471 lv-analysis-in : bot
// 472 lv-analysis-out: bot
double elapsed[64];
// 472 lv-analysis-in : bot
/*****************************************************************/
/******            T  I  M  E  R  _  C  L  E  A  R          ******/
/*****************************************************************/

void timer_clear(int n)
// 473 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_279}
{
// 476 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,n_244}
  elapsed[n] = 0.0;
// 476 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200}
}
// 473 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,n_244}
/*****************************************************************/
/******            T  I  M  E  R  _  S  T  A  R  T          ******/
/*****************************************************************/

void timer_start(int n)
// 477 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,tmp_279}
{
// 480 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,n_245}
  start[n] = elapsed_time();
// 480 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
}
// 477 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,n_245}
/*****************************************************************/
/******            T  I  M  E  R  _  S  T  O  P             ******/
/*****************************************************************/

void timer_stop(int n)
// 481 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_279,tmp_299}
{
// 484 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,n_248,tmp_299}
  double t;
// 484 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,n_248,tmp_299}
// 485 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,n_248,tmp_299}
  double now;
// 485 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,n_248,tmp_299}
// 486 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,n_248,tmp_299}
  now = elapsed_time();
// 486 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,n_248,tmp_299}
// 488 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,now_247,n_248}
  t = now - start[n];
// 488 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t_246,n_248}
// 489 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t_246,n_248}
  elapsed[n] += t;
// 489 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,elapsed_243}
}
// 481 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,n_248,tmp_299}
/*****************************************************************/
/******            T  I  M  E  R  _  R  E  A  D             ******/
/*****************************************************************/

double timer_read(int n)
// 490 lv-analysis-out: {q_147,sx_155,sy_156,gc_160,nthreads_177,verified_178,elapsed_243,tmp_279}
{
// 493 lv-analysis-out: {q_147,sx_155,sy_156,gc_160,nthreads_177,verified_178,elapsed_243,n_249}
  return elapsed[n];
// 493 lv-analysis-in : {q_147,sx_155,sy_156,gc_160,nthreads_177,verified_178}
}
// 490 lv-analysis-in : {q_147,sx_155,sy_156,gc_160,nthreads_177,verified_178,elapsed_243,n_249}

void wtime(double *t)
// 494 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248,tmp_279,tmp_299}
{
// 497 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248,t_252,tmp_299}
  static int sec = - 1;
// 497 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248,sec_250,t_252,tmp_299}
// 498 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248,sec_250,t_252,tmp_299}
  struct timeval tv;
// 498 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248,sec_250,tv_251,t_252,tmp_299}
// 499 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248,sec_250,tv_251,t_252,tmp_299}
  gettimeofday(&tv,((void *)0));
// 499 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248,sec_250,tv_251,t_252,tmp_299}
//  gettimeofday(&tv, (struct timezone *)0);
// 501 lv-analysis-out: bot
  if (
// 502 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248,sec_250,tv_251,t_252}
sec < 0
// 502 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248,sec_250,tv_251,t_252}
) {
// 504 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248,tv_251,t_252}
    sec = tv . tv_sec;
// 504 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248,sec_250,tv_251,t_252}
  }
// 501 lv-analysis-in : bot
// 505 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248,sec_250,tv_251,t_252}
   *t = (tv . tv_sec - sec) + 1.0e-6 * tv . tv_usec;
// 505 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248}
}
// 494 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,t_241,start_242,n_248,t_252,tmp_299}
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
// 506 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_279,tmp_280}
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
// 509 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_262,x_263}
  double t1;
// 509 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_262,x_263}
// 510 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_262,x_263}
  double t2;
// 510 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_262,x_263}
// 511 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_262,x_263}
  double t3;
// 511 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_262,x_263}
// 512 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_262,x_263}
  double t4;
// 512 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_262,x_263}
// 513 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_262,x_263}
  double a1;
// 513 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_262,x_263}
// 514 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_262,x_263}
  double a2;
// 514 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_262,x_263}
// 515 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_262,x_263}
  double x1;
// 515 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_262,x_263}
// 516 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_262,x_263}
  double x2;
// 516 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x2_260,a_262,x_263}
// 517 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x2_260,a_262,x_263}
  double z;
// 517 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x2_260,a_262,x_263}
/*c---------------------------------------------------------------------
c   Break A into two parts such that A = 2^23 * A1 + A2.
c---------------------------------------------------------------------*/
// 518 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x2_260,a_262,x_263}
  t1 = 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * a;
// 518 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t1_253,x2_260,a_262,x_263}
// 519 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t1_253,x2_260,a_262,x_263}
  a1 = ((int )t1);
// 519 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t1_253,a1_257,x2_260,a_262,x_263}
// 520 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t1_253,a1_257,x2_260,a_262,x_263}
  a2 = a - 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * a1;
// 520 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t1_253,a1_257,a2_258,x2_260,x_263}
/*c---------------------------------------------------------------------
c   Break X into two parts such that X = 2^23 * X1 + X2, compute
c   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
c   X = 2^23 * Z + A2 * X2  (mod 2^46).
c---------------------------------------------------------------------*/
// 521 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t1_253,a1_257,a2_258,x2_260,x_263}
  t1 = 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 *  *x;
// 521 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t1_253,a1_257,a2_258,x2_260,x_263}
// 522 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t1_253,a1_257,a2_258,x2_260,x_263}
  x1 = ((int )t1);
// 522 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a1_257,a2_258,x1_259,x2_260,x_263}
// 523 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a1_257,a2_258,x1_259,x2_260,x_263}
  x2 =  *x - 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * x1;
// 523 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a1_257,a2_258,x1_259,x2_260,x_263}
// 524 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a1_257,a2_258,x1_259,x2_260,x_263}
  t1 = a1 * x2 + a2 * x1;
// 524 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t1_253,a2_258,x2_260,x_263}
// 525 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t1_253,a2_258,x2_260,x_263}
  t2 = ((int )(0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * t1));
// 525 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t1_253,t2_254,a2_258,x2_260,x_263}
// 526 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t1_253,t2_254,a2_258,x2_260,x_263}
  z = t1 - 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * t2;
// 526 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a2_258,x2_260,z_261,x_263}
// 527 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a2_258,x2_260,z_261,x_263}
  t3 = 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * z + a2 * x2;
// 527 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t3_255,x_263}
// 528 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t3_255,x_263}
  t4 = ((int )(0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * (0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5) * t3));
// 528 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t3_255,t4_256,x_263}
// 529 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,t3_255,t4_256,x_263}
   *x = t3 - 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * (2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0) * t4;
// 529 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_263}
// 530 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_263}
  return 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * (0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5) *  *x;
// 530 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
}
// 506 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_262,x_263}
/*c---------------------------------------------------------------------
c---------------------------------------------------------------------*/

void vranlc(int n,double *x_seed,double a,double *y)
// 531 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,tmp_279,tmp_280,tmp_281,tmp_282}
{
/* void vranlc (int n, double *x_seed, double a, double y[]) { */
/*c---------------------------------------------------------------------
c---------------------------------------------------------------------*/
/*c---------------------------------------------------------------------
c
c   This routine generates N uniform pseudorandom double precision numbers in
c   the range (0, 1) by using the linear congruential generator
c
c   x_{k+1} = a x_k  (mod 2^46)
c
c   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
c   before repeating.  The argument A is the same as 'a' in the above formula,
c   and X is the same as x_0.  A and X must be odd double precision integers
c   in the range (1, 2^46).  The N results are placed in Y and are normalized
c   to be between 0 and 1.  X is updated to contain the new seed, so that
c   subsequent calls to VRANLC using the same arguments will generate a
c   continuous sequence.  If N is zero, only initialization is performed, and
c   the variables X, A and Y are ignored.
c
c   This routine is the standard version designed for scalar or RISC systems.
c   However, it should produce the same results on any single processor
c   computer with at least 48 mantissa bits in double precision floating point
c   data.  On 64 bit systems, double precision should be disabled.
c
c---------------------------------------------------------------------*/
// 534 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_275,x_seed_276,n_277,y_278}
  int i;
// 534 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_275,x_seed_276,n_277,y_278}
// 535 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_275,x_seed_276,n_277,y_278}
  double x;
// 535 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
// 536 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
  double t1;
// 536 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
// 537 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
  double t2;
// 537 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
// 538 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
  double t3;
// 538 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
// 539 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
  double t4;
// 539 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
// 540 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
  double a1;
// 540 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
// 541 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
  double a2;
// 541 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
// 542 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
  double x1;
// 542 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
// 543 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
  double x2;
// 543 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
// 544 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
  double z;
// 544 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
/*c---------------------------------------------------------------------
c   Break A into two parts such that A = 2^23 * A1 + A2.
c---------------------------------------------------------------------*/
// 545 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a_275,x_seed_276,n_277,y_278}
  t1 = 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * a;
// 545 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,t1_266,a_275,x_seed_276,n_277,y_278}
// 546 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,t1_266,a_275,x_seed_276,n_277,y_278}
  a1 = ((int )t1);
// 546 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a1_270,a_275,x_seed_276,n_277,y_278}
// 547 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a1_270,a_275,x_seed_276,n_277,y_278}
  a2 = a - 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * a1;
// 547 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a1_270,a2_271,x_seed_276,n_277,y_278}
// 548 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a1_270,a2_271,x_seed_276,n_277,y_278}
  x =  *x_seed;
// 548 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a1_270,a2_271,x_seed_276,n_277,y_278}
/*c---------------------------------------------------------------------
c   Generate N results.   This loop is not vectorizable.
c---------------------------------------------------------------------*/
// 549 lv-analysis-out: bot
  for (
// 550 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,a1_270,a2_271,x_seed_276,n_277,y_278}
i = 1
// 550 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,x_265,a1_270,a2_271,x_seed_276,n_277,y_278}
; 
// 551 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,x_265,a1_270,a2_271,x_seed_276,n_277,y_278}
i <= n;
// 551 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,x_265,a1_270,a2_271,x_seed_276,n_277,y_278}
 i++) {
/*c---------------------------------------------------------------------
c   Break X into two parts such that X = 2^23 * X1 + X2, compute
c   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
c   X = 2^23 * Z + A2 * X2  (mod 2^46).
c---------------------------------------------------------------------*/
// 554 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,x_265,a1_270,a2_271,x_seed_276,n_277,y_278}
    t1 = 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * x;
// 554 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,x_265,t1_266,a1_270,a2_271,x_seed_276,n_277,y_278}
// 555 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,x_265,t1_266,a1_270,a2_271,x_seed_276,n_277,y_278}
    x1 = ((int )t1);
// 555 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,x_265,a1_270,a2_271,x1_272,x_seed_276,n_277,y_278}
// 556 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,x_265,a1_270,a2_271,x1_272,x_seed_276,n_277,y_278}
    x2 = x - 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * x1;
// 556 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,a1_270,a2_271,x1_272,x2_273,x_seed_276,n_277,y_278}
// 557 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,a1_270,a2_271,x1_272,x2_273,x_seed_276,n_277,y_278}
    t1 = a1 * x2 + a2 * x1;
// 557 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,t1_266,a1_270,a2_271,x2_273,x_seed_276,n_277,y_278}
// 558 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,t1_266,a1_270,a2_271,x2_273,x_seed_276,n_277,y_278}
    t2 = ((int )(0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * t1));
// 558 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,t1_266,t2_267,a1_270,a2_271,x2_273,x_seed_276,n_277,y_278}
// 559 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,t1_266,t2_267,a1_270,a2_271,x2_273,x_seed_276,n_277,y_278}
    z = t1 - 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * t2;
// 559 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,a1_270,a2_271,x2_273,z_274,x_seed_276,n_277,y_278}
// 560 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,a1_270,a2_271,x2_273,z_274,x_seed_276,n_277,y_278}
    t3 = 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * z + a2 * x2;
// 560 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,t3_268,a1_270,a2_271,x_seed_276,n_277,y_278}
// 561 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,t3_268,a1_270,a2_271,x_seed_276,n_277,y_278}
    t4 = ((int )(0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * (0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5) * t3));
// 561 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,t3_268,t4_269,a1_270,a2_271,x_seed_276,n_277,y_278}
// 562 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,t3_268,t4_269,a1_270,a2_271,x_seed_276,n_277,y_278}
    x = t3 - 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * (2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0 * 2.0) * t4;
// 562 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,x_265,a1_270,a2_271,x_seed_276,n_277,y_278}
// 563 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,x_265,a1_270,a2_271,x_seed_276,n_277,y_278}
    y[i] = 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * (0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5) * x;
// 563 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,i_264,x_265,a1_270,a2_271,x_seed_276,n_277,y_278}
  }
// 549 lv-analysis-in : bot
// 564 lv-analysis-out: {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,x_265,x_seed_276}
   *x_seed = x;
// 564 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242}
}
// 531 lv-analysis-in : {tv_sec_72,tv_usec_73,x_146,q_147,t1_149,sx_155,sy_156,an_158,gc_160,dum_161,np_162,i_166,k_170,k_offset_175,nthreads_177,verified_178,t1_190,t2_191,i_197,ik_198,qq_200,start_242,a_275,x_seed_276,n_277,y_278}
