/*--------------------------------------------------------------------
  
  NAS Parallel Benchmarks 2.3 OpenMP C versions - IS
  This benchmark is an OpenMP C version of the NPB IS code.
  
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
  Author: M. Yarrow
  OpenMP C version: S. Satoh
  
--------------------------------------------------------------------*/
//#include "npb-C.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#if defined(_OPENMP)
#include <omp.h>
#endif /* _OPENMP */
typedef int boolean;
typedef struct {
// 118 lv-analysis-out: bot
double real;
// 118 lv-analysis-in : bot
// 119 lv-analysis-out: bot
double imag;
// 119 lv-analysis-in : bot
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
extern double randlc(double *,double );
extern void vranlc(int ,double *,double ,double *);
void timer_clear(int n);
void timer_start(int n);
void timer_stop(int n);
double timer_read(int n);
void c_print_results(char *name,char cclass,int n1,int n2,int n3,int niter,int nthreads,double t,double mops,char *optype,int passed_verification,char *npbversion,char *compiletime,char *cc,char *clink,char *c_lib,char *c_inc,char *cflags,char *clinkflags,char *rand);
//#include "npbparams.h"
/*
   This file is generated automatically by the setparams utility.
   It sets the number of processors and the classc of the NPB
   in this directory. Do not modify it by hand.   */
#define COMPILETIME "07 Mar 2013"
#define NPBVERSION "2.3"
#define CC "identityTranslator "
#define CFLAGS "-rose:openmp:lowering "
#define CLINK "$(CC)"
#define CLINKFLAGS "-lm"
#define C_LIB "/export/tmp.liao6/workspace/thrifty/build64..."
#define C_INC "-I../common"
#include <stdlib.h>
#include <stdio.h>
#if defined(_OPENMP)
#include <omp.h>
#endif /* _OPENMP */
/*****************************************************************/
/* For serial IS, buckets are not really req'd to solve NPB1 IS  */
/* spec, but their use on some machines improves performance, on */
/* other machines the use of buckets compromises performance,    */
/* probably because it is extra computation which is not req'd.  */
/* (Note: Mechanism not understood, probably cache related)      */
/* Example:  SP2-66MhzWN:  50% speedup with buckets              */
/* Example:  SGI Indy5000: 50% slowdown with buckets             */
/* Example:  SGI O2000:   400% slowdown with buckets (Wow!)      */
/*****************************************************************/
/* #define USE_BUCKETS  */
/* buckets are not used in the OpenMP C version */
#define CLASS 'A'
/******************/
/* default values */
/******************/
#ifndef CLASS
#define CLASS 'A'
#endif
/*************/
/*  CLASS S  */
/*************/
#if CLASS == 'S'
#define  TOTAL_KEYS_LOG_2    16
#define  MAX_KEY_LOG_2       11
#define  NUM_BUCKETS_LOG_2   9
#endif
/*************/
/*  CLASS W  */
/*************/
#if CLASS == 'W'
#define  TOTAL_KEYS_LOG_2    20
#define  MAX_KEY_LOG_2       16
#define  NUM_BUCKETS_LOG_2   10
#endif
/*************/
/*  CLASS A  */
/*************/
#if CLASS == 'A'
#define  TOTAL_KEYS_LOG_2    23
#define  MAX_KEY_LOG_2       19
#define  NUM_BUCKETS_LOG_2   10
#endif
/*************/
/*  CLASS B  */
/*************/
#if CLASS == 'B'
#define  TOTAL_KEYS_LOG_2    25
#define  MAX_KEY_LOG_2       21
#define  NUM_BUCKETS_LOG_2   10
#endif
/*************/
/*  CLASS C  */
/*************/
#if CLASS == 'C'
#define  TOTAL_KEYS_LOG_2    27
#define  MAX_KEY_LOG_2       23
#define  NUM_BUCKETS_LOG_2   10
#endif
#define  TOTAL_KEYS          (1 << TOTAL_KEYS_LOG_2)
#define  MAX_KEY             (1 << MAX_KEY_LOG_2)
#define  NUM_BUCKETS         (1 << NUM_BUCKETS_LOG_2)
#define  NUM_KEYS            TOTAL_KEYS
#define  SIZE_OF_BUFFERS     NUM_KEYS  
#define  MAX_ITERATIONS      10
#define  TEST_ARRAY_SIZE     5
/*************************************/
/* Typedef: if necessary, change the */
/* size of int here by changing the  */
/* int type to, say, long            */
/*************************************/
typedef int INT_TYPE;
/********************/
/* Some global info */
/********************/
/* used by full_verify to get */
// 120 lv-analysis-out: bot
INT_TYPE *key_buff_ptr_global;
// 120 lv-analysis-in : bot
/* copies of rank info        */
// 121 lv-analysis-out: bot
int passed_verification;
// 121 lv-analysis-in : bot
/************************************/
/* These are the three main arrays. */
/* See SIZE_OF_BUFFERS def above    */
/************************************/
// 122 lv-analysis-out: bot
INT_TYPE key_array[1 << 23];
// 122 lv-analysis-in : bot
// 123 lv-analysis-out: bot
INT_TYPE key_buff1[1 << 23];
// 123 lv-analysis-in : bot
// 124 lv-analysis-out: bot
INT_TYPE key_buff2[1 << 23];
// 124 lv-analysis-in : bot
// 125 lv-analysis-out: bot
INT_TYPE partial_verify_vals[5];
// 125 lv-analysis-in : bot
#ifdef USE_BUCKETS
#endif
/**********************/
/* Partial verif info */
/**********************/
// 126 lv-analysis-out: bot
INT_TYPE test_index_array[5];
// 126 lv-analysis-in : bot
// 127 lv-analysis-out: bot
INT_TYPE test_rank_array[5];
// 127 lv-analysis-in : bot
// 128 lv-analysis-out: bot
INT_TYPE S_test_index_array[5] = {(48427), (17148), (23627), (62548), (4431)};
// 128 lv-analysis-in : bot
// 129 lv-analysis-out: bot
INT_TYPE S_test_rank_array[5] = {(0), (18), (346), (64917), (65463)};
// 129 lv-analysis-in : bot
// 130 lv-analysis-out: bot
INT_TYPE W_test_index_array[5] = {(357773), (934767), (875723), (898999), (404505)};
// 130 lv-analysis-in : bot
// 131 lv-analysis-out: bot
INT_TYPE W_test_rank_array[5] = {(1249), (11698), (1039987), (1043896), (1048018)};
// 131 lv-analysis-in : bot
// 132 lv-analysis-out: bot
INT_TYPE A_test_index_array[5] = {(2112377), (662041), (5336171), (3642833), (4250760)};
// 132 lv-analysis-in : bot
// 133 lv-analysis-out: bot
INT_TYPE A_test_rank_array[5] = {(104), (17523), (123928), (8288932), (8388264)};
// 133 lv-analysis-in : bot
// 134 lv-analysis-out: bot
INT_TYPE B_test_index_array[5] = {(41869), (812306), (5102857), (18232239), (26860214)};
// 134 lv-analysis-in : bot
// 135 lv-analysis-out: bot
INT_TYPE B_test_rank_array[5] = {(33422937), (10244), (59149), (33135281), (99)};
// 135 lv-analysis-in : bot
// 136 lv-analysis-out: bot
INT_TYPE C_test_index_array[5] = {(44172927), (72999161), (74326391), (129606274), (21736814)};
// 136 lv-analysis-in : bot
// 137 lv-analysis-out: bot
INT_TYPE C_test_rank_array[5] = {(61147), (882988), (266290), (133997595), (133525895)};
// 137 lv-analysis-in : bot
/***********************/
/* function prototypes */
/***********************/
static double randlc2(double *X,double *A);
void full_verify();
/*
 *    FUNCTION RANDLC (X, A)
 *
 *  This routine returns a uniform pseudorandom double precision number in the
 *  range (0, 1) by using the linear congruential generator
 *
 *  x_{k+1} = a x_k  (mod 2^46)
 *
 *  where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
 *  before repeating.  The argument A is the same as 'a' in the above formula,
 *  and X is the same as x_0.  A and X must be odd double precision integers
 *  in the range (1, 2^46).  The returned value RANDLC is normalized to be
 *  between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
 *  the new seed x_1, so that subsequent calls to RANDLC using the same
 *  arguments will generate a continuous sequence.
 *
 *  This routine should produce the same results on any computer with at least
 *  48 mantissa bits in double precision floating point data.  On Cray systems,
 *  double precision should be disabled.
 *
 *  David H. Bailey     October 26, 1990
 *
 *     IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 *     SAVE KS, R23, R46, T23, T46
 *     DATA KS/0/
 *
 *  If this is the first call to RANDLC, compute R23 = 2 ^ -23, R46 = 2 ^ -46,
 *  T23 = 2 ^ 23, and T46 = 2 ^ 46.  These are computed in loops, rather than
 *  by merely using the ** operator, in order to insure that the results are
 *  exact on all systems.  This code assumes that 0.5D0 is represented exactly.
 */
/*****************************************************************/
/*************           R  A  N  D  L  C             ************/
/*************                                        ************/
/*************    portable random number generator    ************/
/*****************************************************************/

static double randlc2(double *X,double *A)
// 138 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239,tmp_250,tmp_251}
{
// 141 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  static int KS = 0;
// 141 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 142 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  static double R23;
// 142 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 143 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  static double R46;
// 143 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 144 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  static double T23;
// 144 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 145 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  static double T46;
// 145 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 146 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  double T1;
// 146 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 147 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  double T2;
// 147 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 148 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  double T3;
// 148 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 149 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  double T4;
// 149 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 150 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  double A1;
// 150 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 151 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  double A2;
// 151 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 152 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  double X1;
// 152 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 153 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  double X2;
// 153 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 154 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  double Z;
// 154 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 155 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  int i;
// 155 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 156 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  int j;
// 156 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 157 lv-analysis-out: bot
  if (
// 158 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,KS_168,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
KS == 0
// 158 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
) {
// 160 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
    R23 = 1.0;
// 160 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 161 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
    R46 = 1.0;
// 161 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 162 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
    T23 = 1.0;
// 162 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 163 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
    T46 = 1.0;
// 163 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 164 lv-analysis-out: bot
    for (
// 165 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
i = 1
// 165 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,i_182,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
; 
// 166 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,i_182,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
i <= 23;
// 166 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,i_182,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
 i++) {
// 169 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,i_182,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
      R23 = 0.50 * R23;
// 169 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,i_182,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 170 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,i_182,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
      T23 = 2.0 * T23;
// 170 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,i_182,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
    }
// 164 lv-analysis-in : bot
// 171 lv-analysis-out: bot
    for (
// 172 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
i = 1
// 172 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,i_182,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
; 
// 173 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,i_182,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
i <= 46;
// 173 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,i_182,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
 i++) {
// 176 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,i_182,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
      R46 = 0.50 * R46;
// 176 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,i_182,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 177 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,i_182,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
      T46 = 2.0 * T46;
// 177 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,i_182,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
    }
// 171 lv-analysis-in : bot
// 178 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
    KS = 1;
// 178 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  }
// 157 lv-analysis-in : bot
/*  Break A into two parts such that A = 2^23 * A1 + A2 and set X = N.  */
// 179 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  T1 = R23 *  *A;
// 179 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 180 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  j = T1;
// 180 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,j_183,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 181 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,j_183,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  A1 = j;
// 181 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A1_177,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 182 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A1_177,A2_178,X2_180,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  A2 =  *A - T23 * A1;
// 182 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A1_177,A2_178,X2_180,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
/*  Break X into two parts such that X = 2^23 * X1 + X2, compute
    Z = A1 * X2 + A2 * X1  (mod 2^23), and then
    X = 2^23 * Z + A2 * X2  (mod 2^46).                            */
// 183 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A1_177,A2_178,X2_180,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  T1 = R23 *  *X;
// 183 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A1_177,A2_178,X2_180,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 184 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A1_177,A2_178,X2_180,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  j = T1;
// 184 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,A1_177,A2_178,X2_180,j_183,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 185 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,A1_177,A2_178,X2_180,j_183,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  X1 = j;
// 185 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,A1_177,A2_178,X1_179,X2_180,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 186 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,A1_177,A2_178,X1_179,X2_180,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  X2 =  *X - T23 * X1;
// 186 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,A1_177,A2_178,X1_179,X2_180,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 187 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,A1_177,A2_178,X1_179,X2_180,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  T1 = A1 * X2 + A2 * X1;
// 187 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 188 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R23_169,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  j = (R23 * T1);
// 188 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,j_183,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 189 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R46_170,T23_171,T46_172,T1_173,A2_178,X2_180,j_183,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  T2 = j;
// 189 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R46_170,T23_171,T46_172,T1_173,T2_174,A2_178,X2_180,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 190 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R46_170,T23_171,T46_172,T1_173,T2_174,A2_178,X2_180,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  Z = T1 - T23 * T2;
// 190 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R46_170,T23_171,T46_172,A2_178,X2_180,Z_181,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 191 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R46_170,T23_171,T46_172,A2_178,X2_180,Z_181,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  T3 = T23 * Z + A2 * X2;
// 191 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R46_170,T46_172,T3_175,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 192 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R46_170,T46_172,T3_175,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  j = (R46 * T3);
// 192 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R46_170,T46_172,T3_175,j_183,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 193 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R46_170,T46_172,T3_175,j_183,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  T4 = j;
// 193 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R46_170,T46_172,T3_175,T4_176,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 194 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R46_170,T46_172,T3_175,T4_176,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
   *X = T3 - T46 * T4;
// 194 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R46_170,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 195 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,R46_170,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  return R46 *  *X;
// 195 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
}
// 138 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,A_184,X_185,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
/*****************************************************************/
/*************      C  R  E  A  T  E  _  S  E  Q      ************/
/*****************************************************************/

void create_seq(double seed,double a)
// 196 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_250,tmp_251}
{
// 199 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,seed_190,a_191,iteration_209,nthreads_211,start_239}
  double x;
// 199 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 200 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,seed_190,a_191,iteration_209,nthreads_211,start_239}
  int i;
// 200 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 201 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,seed_190,a_191,iteration_209,nthreads_211,start_239}
  int j;
// 201 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 202 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,seed_190,a_191,iteration_209,nthreads_211,start_239}
  int k;
// 202 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 203 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,seed_190,a_191,iteration_209,nthreads_211,start_239}
  k = (1 << 19) / 4;
// 203 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 204 lv-analysis-out: bot
  for (
// 205 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
i = 0
// 205 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
; 
// 206 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
i < 1 << 23;
// 206 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
 i++) {
// 209 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
    x = randlc2(&seed,&a);
// 209 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239,tmp_250,tmp_251}
// 211 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
    x += randlc2(&seed,&a);
// 211 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 212 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
    x += randlc2(&seed,&a);
// 212 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 213 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
    x += randlc2(&seed,&a);
// 213 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,x_186,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
// 214 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,x_186,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
    key_array[i] = (k * x);
// 214 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_187,k_189,seed_190,a_191,iteration_209,nthreads_211,start_239}
  }
// 204 lv-analysis-in : bot
}
// 196 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,seed_190,a_191,iteration_209,nthreads_211,start_239}
/*****************************************************************/
/*************    F  U  L  L  _  V  E  R  I  F  Y     ************/
/*****************************************************************/

void full_verify()
// 215 lv-analysis-out: {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,timecounter_212,tmp_270}
{
// 218 lv-analysis-out: {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,timecounter_212,tmp_270}
  INT_TYPE i;
// 218 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,timecounter_212,tmp_270}
// 219 lv-analysis-out: {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,timecounter_212,tmp_270}
  INT_TYPE j;
// 219 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,timecounter_212,tmp_270}
// 220 lv-analysis-out: {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,timecounter_212,tmp_270}
  INT_TYPE k;
// 220 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,timecounter_212,tmp_270}
// 221 lv-analysis-out: {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,timecounter_212,tmp_270}
  INT_TYPE m;
// 221 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,timecounter_212,tmp_270}
// 222 lv-analysis-out: {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,timecounter_212,tmp_270}
  INT_TYPE unique_keys;
// 222 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,timecounter_212,tmp_270}
/*  Now, finally, sort the keys:  */
// 223 lv-analysis-out: bot
  for (
// 224 lv-analysis-out: {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,timecounter_212,tmp_270}
i = 0
// 224 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,i_192,nthreads_211,timecounter_212,tmp_270}
; 
// 225 lv-analysis-out: {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,i_192,nthreads_211,timecounter_212,tmp_270}
i < 1 << 23;
// 225 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,i_192,nthreads_211,timecounter_212,tmp_270}
 i++) {
// 228 lv-analysis-out: {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,i_192,nthreads_211,timecounter_212,tmp_270}
    key_array[--key_buff_ptr_global[key_buff2[i]]] = key_buff2[i];
// 228 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,i_192,nthreads_211,timecounter_212,tmp_270}
  }
// 223 lv-analysis-in : bot
/*  Confirm keys correctly sorted: count incorrectly sorted keys, if any */
// 229 lv-analysis-out: {passed_verification_151,key_array_152,nthreads_211,timecounter_212,tmp_270}
  j = 0;
// 229 lv-analysis-in : {passed_verification_151,key_array_152,j_193,nthreads_211,timecounter_212,tmp_270}
// 230 lv-analysis-out: bot
  for (
// 231 lv-analysis-out: {passed_verification_151,key_array_152,j_193,nthreads_211,timecounter_212,tmp_270}
i = 1
// 231 lv-analysis-in : {passed_verification_151,key_array_152,i_192,j_193,nthreads_211,timecounter_212,tmp_270}
; 
// 232 lv-analysis-out: {passed_verification_151,key_array_152,i_192,j_193,nthreads_211,timecounter_212,tmp_270}
i < 1 << 23;
// 232 lv-analysis-in : {passed_verification_151,key_array_152,i_192,j_193,nthreads_211,timecounter_212,tmp_270}
 i++) {
// 235 lv-analysis-out: bot
    if (
// 236 lv-analysis-out: {passed_verification_151,key_array_152,i_192,j_193,nthreads_211,timecounter_212,tmp_270}
key_array[i - 1] > key_array[i]
// 236 lv-analysis-in : {passed_verification_151,key_array_152,i_192,j_193,nthreads_211,timecounter_212,tmp_270}
) {
// 238 lv-analysis-out: {passed_verification_151,key_array_152,i_192,j_193,nthreads_211,timecounter_212,tmp_270}
      j++;
// 238 lv-analysis-in : {passed_verification_151,key_array_152,i_192,j_193,nthreads_211,timecounter_212,tmp_270}
    }
// 235 lv-analysis-in : bot
  }
// 230 lv-analysis-in : bot
// 239 lv-analysis-out: bot
  if (
// 240 lv-analysis-out: {passed_verification_151,j_193,nthreads_211,timecounter_212,tmp_270}
j != 0
// 240 lv-analysis-in : {passed_verification_151,j_193,nthreads_211,timecounter_212,tmp_270}
) {
// 242 lv-analysis-out: {passed_verification_151,j_193,nthreads_211,timecounter_212,tmp_270}
    printf("Full_verify: number of keys out of sort: %d\n",j);
// 242 lv-analysis-in : {passed_verification_151,nthreads_211,timecounter_212,tmp_270}
  }
   else {
// 245 lv-analysis-out: {passed_verification_151,nthreads_211,timecounter_212}
    passed_verification++;
// 245 lv-analysis-in : {passed_verification_151,nthreads_211,timecounter_212}
  }
// 239 lv-analysis-in : bot
}
// 215 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,timecounter_212,tmp_270}
/*****************************************************************/
/*************             R  A  N  K             ****************/
/*****************************************************************/

void rank(int iteration)
// 246 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_250,tmp_270}
{
// 249 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  INT_TYPE i;
// 249 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
// 250 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  INT_TYPE j;
// 250 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
// 251 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  INT_TYPE k;
// 251 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
// 252 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  INT_TYPE l;
// 252 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
// 253 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  INT_TYPE m;
// 253 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
// 254 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  INT_TYPE shift = 19 - 10;
// 254 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
// 255 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  INT_TYPE key;
// 255 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
// 256 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  INT_TYPE min_key_val;
// 256 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
// 257 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  INT_TYPE max_key_val;
// 257 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
// 258 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  INT_TYPE prv_buff1[1 << 19];
// 258 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
// 259 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  
#pragma omp master
// 259 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
{
// 261 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
    key_array[iteration] = iteration;
// 261 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
// 262 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
    key_array[iteration + 10] = (1 << 19) - iteration;
// 262 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
/*  Determine where the partial verify test keys are, load into  */
/*  top of array bucket_size                                     */
// 263 lv-analysis-out: bot
    for (
// 264 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i = 0
// 264 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
; 
// 265 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i < 5;
// 265 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
 i++) {
// 268 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
      partial_verify_vals[i] = key_array[test_index_array[i]];
// 268 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
    }
// 263 lv-analysis-in : bot
/*  Clear the work array */
// 269 lv-analysis-out: bot
    for (
// 270 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i = 0
// 270 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
; 
// 271 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i < 1 << 19;
// 271 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
 i++) {
// 274 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
      key_buff1[i] = 0;
// 274 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
    }
// 269 lv-analysis-in : bot
  }
// 275 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  
#pragma omp barrier
// 275 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
// 276 lv-analysis-out: bot
  for (
// 277 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i = 0
// 277 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
; 
// 278 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i < 1 << 19;
// 278 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
 i++) {
// 281 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
    prv_buff1[i] = 0;
// 281 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  }
// 276 lv-analysis-in : bot
/*  Copy keys into work array; keys in key_array will be reused each iter. */
// 282 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  
#pragma omp for nowait
// 282 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
// 283 lv-analysis-out: bot
  for (
// 284 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i = 0
// 284 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
; 
// 285 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i < 1 << 23;
// 285 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
 i++) {
// 288 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
    key_buff2[i] = key_array[i];
// 288 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
/*  Ranking of all keys occurs in this section:                 */
/*  In this section, the keys themselves are used as their 
    own indexes to determine how many of each there are: their
    individual population                                       */
/* Now they have individual key   */
// 289 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
    prv_buff1[key_buff2[i]]++;
// 289 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  }
// 283 lv-analysis-in : bot
/* population                     */
// 290 lv-analysis-out: bot
  for (
// 291 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i = 0
// 291 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
; 
// 292 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i < (1 << 19) - 1;
// 292 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
 i++) {
// 295 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
    prv_buff1[i + 1] += prv_buff1[i];
// 295 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  }
// 290 lv-analysis-in : bot
// 296 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  
#pragma omp critical
// 296 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
{
// 298 lv-analysis-out: bot
    for (
// 299 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i = 0
// 299 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
; 
// 300 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i < 1 << 19;
// 300 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
 i++) {
// 303 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
      key_buff1[i] += prv_buff1[i];
// 303 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,prv_buff1_206,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
    }
// 298 lv-analysis-in : bot
  }
/*  To obtain ranks of each key, successively add the individual key
    population, not forgetting to add m, the total of lesser keys,
    to the first key population                                          */
// 304 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  
#pragma omp barrier
// 304 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
// 305 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
  
#pragma omp master
// 305 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
{
/* This is the partial verify test section */
/* Observe that test_rank_array vals are   */
/* shifted differently for different cases */
// 307 lv-analysis-out: bot
    for (
// 308 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i = 0
// 308 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
; 
// 309 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i < 5;
// 309 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
 i++) {
/* test vals were put here */
// 312 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
      k = partial_verify_vals[i];
// 312 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
// 313 lv-analysis-out: bot
      if (
// 314 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
0 <= k && k <= (1 << 23) - 1
// 314 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
) {
        switch(
// 317 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
'A'
// 317 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
){
          case 'S':
{
// 321 lv-analysis-out: bot
            if (
// 322 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i <= 2
// 322 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
) {
// 324 lv-analysis-out: bot
              if (
// 325 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
key_buff1[k - 1] != test_rank_array[i] + iteration
// 325 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
) {
// 327 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                printf("Failed partial verification: iteration %d, test key %d\n",iteration,i);
// 327 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
               else {
// 330 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                passed_verification++;
// 330 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
// 324 lv-analysis-in : bot
            }
             else {
// 332 lv-analysis-out: bot
              if (
// 333 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
key_buff1[k - 1] != test_rank_array[i] - iteration
// 333 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
) {
// 335 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                printf("Failed partial verification: iteration %d, test key %d\n",iteration,i);
// 335 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
               else {
// 338 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                passed_verification++;
// 338 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
// 332 lv-analysis-in : bot
            }
// 321 lv-analysis-in : bot
// 339 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
            break; 
// 339 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
          }
// 319 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
          case 'W':
{
// 342 lv-analysis-out: bot
            if (
// 343 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i < 2
// 343 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
) {
// 345 lv-analysis-out: bot
              if (
// 346 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
key_buff1[k - 1] != test_rank_array[i] + (iteration - 2)
// 346 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
) {
// 348 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                printf("Failed partial verification: iteration %d, test key %d\n",iteration,i);
// 348 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
               else {
// 351 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                passed_verification++;
// 351 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
// 345 lv-analysis-in : bot
            }
             else {
// 353 lv-analysis-out: bot
              if (
// 354 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
key_buff1[k - 1] != test_rank_array[i] - iteration
// 354 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
) {
// 356 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                printf("Failed partial verification: iteration %d, test key %d\n",iteration,i);
// 356 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
               else {
// 359 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                passed_verification++;
// 359 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
// 353 lv-analysis-in : bot
            }
// 342 lv-analysis-in : bot
// 360 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
            break; 
// 360 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
          }
// 340 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
          case 'A':
{
// 363 lv-analysis-out: bot
            if (
// 364 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i <= 2
// 364 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
) {
// 366 lv-analysis-out: bot
              if (
// 367 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
key_buff1[k - 1] != test_rank_array[i] + (iteration - 1)
// 367 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
) {
// 369 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                printf("Failed partial verification: iteration %d, test key %d\n",iteration,i);
// 369 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
               else {
// 372 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                passed_verification++;
// 372 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
// 366 lv-analysis-in : bot
            }
             else {
// 374 lv-analysis-out: bot
              if (
// 375 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
key_buff1[k - 1] != test_rank_array[i] - (iteration - 1)
// 375 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
) {
// 377 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                printf("Failed partial verification: iteration %d, test key %d\n",iteration,i);
// 377 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
               else {
// 380 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                passed_verification++;
// 380 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
// 374 lv-analysis-in : bot
            }
// 363 lv-analysis-in : bot
// 381 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
            break; 
// 381 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
          }
// 361 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
          case 'B':
{
// 384 lv-analysis-out: bot
            if (
// 385 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i == 1 || i == 2 || i == 4
// 385 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
) {
// 387 lv-analysis-out: bot
              if (
// 388 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
key_buff1[k - 1] != test_rank_array[i] + iteration
// 388 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
) {
// 390 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                printf("Failed partial verification: iteration %d, test key %d\n",iteration,i);
// 390 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
               else {
// 393 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                passed_verification++;
// 393 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
// 387 lv-analysis-in : bot
            }
             else {
// 395 lv-analysis-out: bot
              if (
// 396 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
key_buff1[k - 1] != test_rank_array[i] - iteration
// 396 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
) {
// 398 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                printf("Failed partial verification: iteration %d, test key %d\n",iteration,i);
// 398 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
               else {
// 401 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                passed_verification++;
// 401 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
// 395 lv-analysis-in : bot
            }
// 384 lv-analysis-in : bot
// 402 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
            break; 
// 402 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
          }
// 382 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
          case 'C':
{
// 405 lv-analysis-out: bot
            if (
// 406 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
i <= 2
// 406 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
) {
// 408 lv-analysis-out: bot
              if (
// 409 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
key_buff1[k - 1] != test_rank_array[i] + iteration
// 409 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
) {
// 411 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                printf("Failed partial verification: iteration %d, test key %d\n",iteration,i);
// 411 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
               else {
// 414 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                passed_verification++;
// 414 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
// 408 lv-analysis-in : bot
            }
             else {
// 416 lv-analysis-out: bot
              if (
// 417 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
key_buff1[k - 1] != test_rank_array[i] - iteration
// 417 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
) {
// 419 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                printf("Failed partial verification: iteration %d, test key %d\n",iteration,i);
// 419 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
               else {
// 422 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
                passed_verification++;
// 422 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
              }
// 416 lv-analysis-in : bot
            }
// 405 lv-analysis-in : bot
// 423 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
            break; 
// 423 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
          }
// 403 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,i_197,k_199,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
        }
      }
// 313 lv-analysis-in : bot
    }
// 307 lv-analysis-in : bot
/*  Make copies of rank info for use by full_verify: these variables
    in rank are local; making them global slows down the code, probably
    since they cannot be made register by compiler                        */
// 424 lv-analysis-out: bot
    if (
// 425 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239}
iteration == 10
// 425 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239}
) {
// 427 lv-analysis-out: {tv_sec_72,tv_usec_73,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239}
      key_buff_ptr_global = key_buff1;
// 427 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239}
    }
// 424 lv-analysis-in : bot
/* end master */
  }
}
// 246 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_207,iteration_209,nthreads_211,start_239,tmp_270}
/*****************************************************************/
/*************             M  A  I  N             ****************/
/*****************************************************************/

int main(int argc,char **argv)
// 428 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,start_239,tmp_250,tmp_251,tmp_270}
{
// 431 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,start_239,tmp_270}
  int i;
// 431 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,start_239,tmp_270}
// 432 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,start_239,tmp_270}
  int iteration;
// 432 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,iteration_209,start_239,tmp_270}
// 433 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,iteration_209,start_239,tmp_270}
  int itemp;
// 433 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,iteration_209,start_239,tmp_270}
// 434 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,iteration_209,start_239,tmp_270}
  int nthreads = 1;
// 434 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,iteration_209,nthreads_211,start_239,tmp_270}
// 435 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,iteration_209,nthreads_211,start_239,tmp_270}
  double timecounter;
// 435 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,iteration_209,nthreads_211,start_239,tmp_270}
// 436 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,iteration_209,nthreads_211,start_239,tmp_270}
  double maxtime;
// 436 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,iteration_209,nthreads_211,start_239,tmp_270}
/*  Initialize the verification arrays if a valid class */
// 437 lv-analysis-out: bot
  for (
// 438 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,iteration_209,nthreads_211,start_239,tmp_270}
i = 0
// 438 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
; 
// 439 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
i < 5;
// 439 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
 i++) {
    switch(
// 443 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
'A'
// 443 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
){
      case 'S':
{
// 447 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
        test_index_array[i] = S_test_index_array[i];
// 447 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
// 448 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
        test_rank_array[i] = S_test_rank_array[i];
// 448 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
// 449 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
        break; 
// 449 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
      }
// 445 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
      case 'A':
{
// 452 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
        test_index_array[i] = A_test_index_array[i];
// 452 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
// 453 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
        test_rank_array[i] = A_test_rank_array[i];
// 453 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
// 454 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
        break; 
// 454 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
      }
// 450 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
      case 'W':
{
// 457 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
        test_index_array[i] = W_test_index_array[i];
// 457 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
// 458 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
        test_rank_array[i] = W_test_rank_array[i];
// 458 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
// 459 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
        break; 
// 459 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
      }
// 455 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
      case 'B':
{
// 462 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
        test_index_array[i] = B_test_index_array[i];
// 462 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
// 463 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
        test_rank_array[i] = B_test_rank_array[i];
// 463 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
// 464 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
        break; 
// 464 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
      }
// 460 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
      case 'C':
{
// 467 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
        test_index_array[i] = C_test_index_array[i];
// 467 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
// 468 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
        test_rank_array[i] = C_test_rank_array[i];
// 468 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
// 469 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
        break; 
// 469 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
      }
// 465 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,i_208,iteration_209,nthreads_211,start_239,tmp_270}
    }
  }
// 437 lv-analysis-in : bot
// 470 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
  ;
// 470 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
/*  Printout initial NPB info */
// 471 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
  printf("\n\n NAS Parallel Benchmarks 2.3 OpenMP C version - IS Benchmark\n\n");
// 471 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
// 473 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
  int __temp0__ = 1 << 23;
// 473 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,__temp0___214,start_239,tmp_270}
// 474 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,__temp0___214,start_239,tmp_270}
  printf(" Size:  %d  (class %c)\n",__temp0__,'A');
// 474 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
// 476 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
  printf(" Iterations:   %d\n",10);
// 476 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
/*  Initialize timer  */
// 478 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239}
  timer_clear(0);
// 478 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_250}
/*  Generate random number sequence and subsequent keys on all procs */
/* Random number gen seed */
// 480 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239}
  create_seq(314159265.00,1220703125.00);
// 480 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_250,tmp_251}
/* Random number gen mult */
/*  Do one interation for free (i.e., untimed) to guarantee initialization of  
    all data and code pages and respective tables */
// 482 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
  
#pragma omp parallel
// 482 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
// 483 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
  rank(1);
// 483 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_250,tmp_270}
/*  Start verification counter */
// 485 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,nthreads_211,tmp_270}
  passed_verification = 0;
// 485 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,nthreads_211,tmp_270}
// 486 lv-analysis-out: bot
  if (
// 487 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,nthreads_211,tmp_270}
'A' != 'S'
// 487 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,nthreads_211,tmp_270}
) {
// 489 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,nthreads_211,tmp_270}
    printf("\n   iteration\n");
// 489 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,nthreads_211,tmp_270}
  }
// 486 lv-analysis-in : bot
/*  Start timer  */
// 491 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,nthreads_211}
  timer_start(0);
// 491 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,nthreads_211,tmp_250}
/*  This is the main iteration */
// 493 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,nthreads_211,start_239,tmp_270}
  
#pragma omp parallel private(iteration)
// 493 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,nthreads_211,start_239,tmp_270}
// 494 lv-analysis-out: bot
  for (
// 495 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,nthreads_211,start_239,tmp_270}
iteration = 1
// 495 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
; 
// 496 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
iteration <= 10;
// 496 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
 iteration++) {
// 499 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
    
#pragma omp master
// 499 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
// 500 lv-analysis-out: bot
    if (
// 501 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
'A' != 'S'
// 501 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
) {
// 503 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
      printf("        %d\n",iteration);
// 503 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
    }
// 500 lv-analysis-in : bot
// 505 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_270}
    rank(iteration);
// 505 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_250,tmp_270}
#if defined(_OPENMP)	
#endif /* _OPENMP */	
  }
// 494 lv-analysis-in : bot
/*  End of timing, obtain maximum time of all processors */
// 507 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,start_239,tmp_270}
  timer_stop(0);
// 507 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,start_239,tmp_250,tmp_270}
// 509 lv-analysis-out: {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,elapsed_240}
  timecounter = timer_read(0);
// 509 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,elapsed_240,tmp_250}
/*  This tests that keys are in sequence: sorting of last ranked key seq
    occurs here, but is an untimed operation                             */
// 511 lv-analysis-out: {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,timecounter_212,tmp_270}
  full_verify();
// 511 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,timecounter_212,tmp_270}
/*  The final printout  */
// 513 lv-analysis-out: bot
  if (
// 514 lv-analysis-out: {passed_verification_151,nthreads_211,timecounter_212,tmp_270}
passed_verification != 5 * 10 + 1
// 514 lv-analysis-in : {passed_verification_151,nthreads_211,timecounter_212,tmp_270}
) {
// 516 lv-analysis-out: {nthreads_211,timecounter_212,tmp_270}
    passed_verification = 0;
// 516 lv-analysis-in : {passed_verification_151,nthreads_211,timecounter_212,tmp_270}
  }
// 513 lv-analysis-in : bot
// 517 lv-analysis-out: {passed_verification_151,nthreads_211,timecounter_212,tmp_270}
  int __temp1__ = 1 << 23;
// 517 lv-analysis-in : {passed_verification_151,nthreads_211,timecounter_212,__temp1___215,tmp_270}
// 518 lv-analysis-out: {passed_verification_151,nthreads_211,timecounter_212,__temp1___215,tmp_270}
  double __temp2__ = ((double )(10 * (1 << 23))) / timecounter / 1000000.;
// 518 lv-analysis-in : {passed_verification_151,nthreads_211,timecounter_212,__temp1___215,__temp2___216,tmp_270}
// 519 lv-analysis-out: {passed_verification_151,nthreads_211,timecounter_212,__temp1___215,__temp2___216,tmp_270}
  c_print_results("IS",'A',__temp1__,0,0,10,nthreads,timecounter,__temp2__,"keys ranked",passed_verification,"2.3","07 Mar 2013","identityTranslator ","$(CC)","/export/tmp.liao6/workspace/thrifty/build64...","-I../common","-rose:openmp:lowering ","-lm","randlc2");
// 519 lv-analysis-in : {tmp_250,tmp_251,tmp_252,tmp_253,tmp_254,tmp_255,tmp_256,tmp_257,tmp_258,tmp_259,tmp_260,tmp_261,tmp_262,tmp_263,tmp_264,tmp_265,tmp_266,tmp_267,tmp_268,tmp_269,tmp_270}
// 521 lv-analysis-out: {}
  return 0;
// 521 lv-analysis-in : {}
/**************************/
/*  E N D  P R O G R A M  */
}
// 428 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,S_test_index_array_158,S_test_rank_array_159,W_test_index_array_160,W_test_rank_array_161,A_test_index_array_162,A_test_rank_array_163,B_test_index_array_164,B_test_rank_array_165,C_test_index_array_166,C_test_rank_array_167,start_239,tmp_270}
/**************************/
/* cat ./common/c_print_results.c */
/*****************************************************************/
/******     C  _  P  R  I  N  T  _  R  E  S  U  L  T  S     ******/
/*****************************************************************/

void c_print_results(char *name,char cclass,int n1,int n2,int n3,int niter,int nthreads,double t,double mops,char *optype,int passed_verification,char *npbversion,char *compiletime,char *cc,char *clink,char *c_lib,char *c_inc,char *cflags,char *clinkflags,char *rand)
// 522 lv-analysis-out: {tmp_250,tmp_251,tmp_252,tmp_253,tmp_254,tmp_255,tmp_256,tmp_257,tmp_258,tmp_259,tmp_260,tmp_261,tmp_262,tmp_263,tmp_264,tmp_265,tmp_266,tmp_267,tmp_268,tmp_269,tmp_270}
{
// 525 lv-analysis-out: {name_218,cclass_219,n2_220,n3_221,n1_222,niter_223,nthreads_224,t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  char *evalue = "1000";
// 525 lv-analysis-in : {name_218,cclass_219,n2_220,n3_221,n1_222,niter_223,nthreads_224,t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
// 526 lv-analysis-out: {name_218,cclass_219,n2_220,n3_221,n1_222,niter_223,nthreads_224,t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  printf("\n\n %s Benchmark Completed\n",name);
// 526 lv-analysis-in : {cclass_219,n2_220,n3_221,n1_222,niter_223,nthreads_224,t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
// 528 lv-analysis-out: {cclass_219,n2_220,n3_221,n1_222,niter_223,nthreads_224,t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  printf(" Class           =                        %c\n",cclass);
// 528 lv-analysis-in : {n2_220,n3_221,n1_222,niter_223,nthreads_224,t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
// 530 lv-analysis-out: bot
  if (
// 531 lv-analysis-out: {n2_220,n3_221,n1_222,niter_223,nthreads_224,t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
n2 == 0 && n3 == 0
// 531 lv-analysis-in : {n2_220,n3_221,n1_222,niter_223,nthreads_224,t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
) {
/* as in IS */
// 533 lv-analysis-out: {n1_222,niter_223,nthreads_224,t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
    printf(" Size            =             %12d\n",n1);
// 533 lv-analysis-in : {niter_223,nthreads_224,t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  }
   else {
// 536 lv-analysis-out: {n2_220,n3_221,n1_222,niter_223,nthreads_224,t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
    printf(" Size            =              %3dx%3dx%3d\n",n1,n2,n3);
// 536 lv-analysis-in : {niter_223,nthreads_224,t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  }
// 530 lv-analysis-in : bot
// 538 lv-analysis-out: {niter_223,nthreads_224,t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  printf(" Iterations      =             %12d\n",niter);
// 538 lv-analysis-in : {nthreads_224,t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
// 540 lv-analysis-out: {nthreads_224,t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  printf(" Threads         =             %12d\n",nthreads);
// 540 lv-analysis-in : {t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
// 542 lv-analysis-out: {t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  printf(" Time in seconds =             %12.2f\n",t);
// 542 lv-analysis-in : {mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
// 544 lv-analysis-out: {mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  printf(" Mop/s total     =             %12.2f\n",mops);
// 544 lv-analysis-in : {optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
// 546 lv-analysis-out: {optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  printf(" Operation type  = %24s\n",optype);
// 546 lv-analysis-in : {passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
// 548 lv-analysis-out: bot
  if (
// 549 lv-analysis-out: {passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
passed_verification
// 549 lv-analysis-in : {npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
) {
// 551 lv-analysis-out: {npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
    printf(" Verification    =               SUCCESSFUL\n");
// 551 lv-analysis-in : {npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  }
   else {
// 554 lv-analysis-out: {npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
    printf(" Verification    =             UNSUCCESSFUL\n");
// 554 lv-analysis-in : {npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  }
// 548 lv-analysis-in : bot
// 556 lv-analysis-out: {npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  printf(" Version         =             %12s\n",npbversion);
// 556 lv-analysis-in : {compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
// 558 lv-analysis-out: {compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  printf(" Compile date    =             %12s\n",compiletime);
// 558 lv-analysis-in : {cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
// 560 lv-analysis-out: {cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  printf("\n Compile options:\n");
// 560 lv-analysis-in : {cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
// 562 lv-analysis-out: {cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  printf("    CC           = %s\n",cc);
// 562 lv-analysis-in : {clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
// 564 lv-analysis-out: {clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  printf("    CLINK        = %s\n",clink);
// 564 lv-analysis-in : {c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
// 566 lv-analysis-out: {c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  printf("    C_LIB        = %s\n",c_lib);
// 566 lv-analysis-in : {c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
// 568 lv-analysis-out: {c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
  printf("    C_INC        = %s\n",c_inc);
// 568 lv-analysis-in : {cflags_235,clinkflags_236,rand_237,tmp_270}
// 570 lv-analysis-out: {cflags_235,clinkflags_236,rand_237,tmp_270}
  printf("    CFLAGS       = %s\n",cflags);
// 570 lv-analysis-in : {clinkflags_236,rand_237,tmp_270}
// 572 lv-analysis-out: {clinkflags_236,rand_237,tmp_270}
  printf("    CLINKFLAGS   = %s\n",clinkflags);
// 572 lv-analysis-in : {rand_237,tmp_270}
// 574 lv-analysis-out: {rand_237,tmp_270}
  printf("    RAND         = %s\n",rand);
// 574 lv-analysis-in : {tmp_270}
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
// 522 lv-analysis-in : {name_218,cclass_219,n2_220,n3_221,n1_222,niter_223,nthreads_224,t_225,mops_226,optype_227,passed_verification_228,npbversion_229,compiletime_230,cc_231,clink_232,c_lib_233,c_inc_234,cflags_235,clinkflags_236,rand_237,tmp_270}
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
// 576 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,start_239,n_245,tmp_270}
{
// 579 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,start_239,n_245,tmp_270}
  double t;
// 579 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245,tmp_270}
// 580 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245,tmp_270}
  wtime(&t);
// 580 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245,tmp_250,tmp_270}
// 582 lv-analysis-out: {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245}
  return t;
// 582 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,start_239,n_245}
}
// 576 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,start_239,n_245,tmp_270}
// 583 lv-analysis-out: bot
double start[64];
// 583 lv-analysis-in : bot
// 584 lv-analysis-out: bot
double elapsed[64];
// 584 lv-analysis-in : bot
/*****************************************************************/
/******            T  I  M  E  R  _  C  L  E  A  R          ******/
/*****************************************************************/

void timer_clear(int n)
// 585 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,tmp_250}
{
// 588 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,n_241}
  elapsed[n] = 0.0;
// 588 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239}
}
// 585 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,iteration_209,nthreads_211,start_239,n_241}
/*****************************************************************/
/******            T  I  M  E  R  _  S  T  A  R  T          ******/
/*****************************************************************/

void timer_start(int n)
// 589 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,nthreads_211,tmp_250}
{
// 592 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,nthreads_211,n_242}
  start[n] = elapsed_time();
// 592 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,nthreads_211,start_239}
}
// 589 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff1_153,key_buff2_154,partial_verify_vals_155,test_index_array_156,test_rank_array_157,nthreads_211,n_242}
/*****************************************************************/
/******            T  I  M  E  R  _  S  T  O  P             ******/
/*****************************************************************/

void timer_stop(int n)
// 593 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,start_239,tmp_250,tmp_270}
{
// 596 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,start_239,n_245,tmp_270}
  double t;
// 596 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,start_239,n_245,tmp_270}
// 597 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,start_239,n_245,tmp_270}
  double now;
// 597 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,start_239,n_245,tmp_270}
// 598 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,start_239,n_245,tmp_270}
  now = elapsed_time();
// 598 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,start_239,n_245,tmp_270}
// 600 lv-analysis-out: {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,start_239,now_244,n_245}
  t = now - start[n];
// 600 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_243,n_245}
// 601 lv-analysis-out: {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_243,n_245}
  elapsed[n] += t;
// 601 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,elapsed_240}
}
// 593 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,start_239,n_245,tmp_270}
/*****************************************************************/
/******            T  I  M  E  R  _  R  E  A  D             ******/
/*****************************************************************/

double timer_read(int n)
// 602 lv-analysis-out: {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,elapsed_240,tmp_250}
{
// 605 lv-analysis-out: {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,elapsed_240,n_246}
  return elapsed[n];
// 605 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211}
}
// 602 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,elapsed_240,n_246}

void wtime(double *t)
// 606 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245,tmp_250,tmp_270}
{
// 609 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245,t_249,tmp_270}
  static int sec = - 1;
// 609 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245,sec_247,t_249,tmp_270}
// 610 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245,sec_247,t_249,tmp_270}
  struct timeval tv;
// 610 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245,sec_247,tv_248,t_249,tmp_270}
// 611 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245,sec_247,tv_248,t_249,tmp_270}
  gettimeofday(&tv,((void *)0));
// 611 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245,sec_247,tv_248,t_249,tmp_270}
//  gettimeofday(&tv, (struct timezone *)0);
// 613 lv-analysis-out: bot
  if (
// 614 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245,sec_247,tv_248,t_249}
sec < 0
// 614 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245,sec_247,tv_248,t_249}
) {
// 616 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245,tv_248,t_249}
    sec = tv . tv_sec;
// 616 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245,sec_247,tv_248,t_249}
  }
// 613 lv-analysis-in : bot
// 617 lv-analysis-out: {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245,sec_247,tv_248,t_249}
   *t = (tv . tv_sec - sec) + 1.0e-6 * tv . tv_usec;
// 617 lv-analysis-in : {key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245}
}
// 606 lv-analysis-in : {tv_sec_72,tv_usec_73,key_buff_ptr_global_150,passed_verification_151,key_array_152,key_buff2_154,nthreads_211,t_238,start_239,n_245,t_249,tmp_270}
