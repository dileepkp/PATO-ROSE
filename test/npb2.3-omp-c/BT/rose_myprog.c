#define MSIZE 200
// 0 lv-analysis-out: bot
int n;
// 0 lv-analysis-in : {}
// 1 lv-analysis-out: bot
int m;
// 1 lv-analysis-in : {}
// 2 lv-analysis-out: bot
int mits;
// 2 lv-analysis-in : {}
// 3 lv-analysis-out: bot
double tol;
// 3 lv-analysis-in : {}
// 4 lv-analysis-out: bot
double relax = 1.0;
// 4 lv-analysis-in : {}
// 5 lv-analysis-out: bot
double alpha = 0.0543;
// 5 lv-analysis-in : {}
// 6 lv-analysis-out: bot
double u[200][200];
// 6 lv-analysis-in : {}
// 7 lv-analysis-out: bot
double f[200][200];
// 7 lv-analysis-in : {}
// 8 lv-analysis-out: bot
double uold[200][200];
// 8 lv-analysis-in : {}
// 9 lv-analysis-out: bot
double dx;
// 9 lv-analysis-in : {}
// 10 lv-analysis-out: bot
double dy;
// 10 lv-analysis-in : {}

void initialize()
// 11 lv-analysis-out: {n_0,m_1,alpha_5}
{
// 14 lv-analysis-out: {n_0,m_1,alpha_5}
  int i;
// 14 lv-analysis-in : {n_0,m_1,alpha_5}
// 15 lv-analysis-out: {n_0,m_1,alpha_5}
  int j;
// 15 lv-analysis-in : {n_0,m_1,alpha_5}
// 16 lv-analysis-out: {n_0,m_1,alpha_5}
  int xx;
// 16 lv-analysis-in : {n_0,m_1,alpha_5}
// 17 lv-analysis-out: {n_0,m_1,alpha_5}
  int yy;
// 17 lv-analysis-in : {n_0,m_1,alpha_5}
//  double PI = 3.1415926;
// -->dx@112:2
// 18 lv-analysis-out: {n_0,m_1,alpha_5}
  dx = 2.0 / (n - 1);
// 18 lv-analysis-in : {n_0,m_1,alpha_5,dx_9}
//-->dy@113:2
// 19 lv-analysis-out: {n_0,m_1,alpha_5,dx_9}
  dy = 2.0 / (m - 1);
// 19 lv-analysis-in : {n_0,m_1,alpha_5,dx_9,dy_10}
/* Initialize initial condition and RHS */
//#pragma omp parallel for private(i,j,xx,yy)
// 20 lv-analysis-out: bot
  for (
// 21 lv-analysis-out: {n_0,m_1,alpha_5,dx_9,dy_10}
i = 0
// 21 lv-analysis-in : {n_0,m_1,alpha_5,dx_9,dy_10,i_11}
; 
// 22 lv-analysis-out: {n_0,m_1,alpha_5,dx_9,dy_10,i_11}
i < n;
// 22 lv-analysis-in : {n_0,m_1,alpha_5,dx_9,dy_10,i_11}
 i++) 
// 24 lv-analysis-out: bot
    for (
// 25 lv-analysis-out: {n_0,m_1,alpha_5,dx_9,dy_10,i_11}
j = 0
// 25 lv-analysis-in : {n_0,m_1,alpha_5,dx_9,dy_10,i_11,j_12}
; 
// 26 lv-analysis-out: {n_0,m_1,alpha_5,dx_9,dy_10,i_11,j_12}
j < m;
// 26 lv-analysis-in : {n_0,m_1,alpha_5,dx_9,dy_10,i_11,j_12}
 j++) {
/* -1 < x < 1 */
// 29 lv-analysis-out: {n_0,m_1,alpha_5,dx_9,dy_10,i_11,j_12}
      xx = ((int )(- 1.0 + dx * (i - 1)));
// 29 lv-analysis-in : {n_0,m_1,alpha_5,dx_9,dy_10,i_11,j_12,xx_13}
/* -1 < y < 1 */
// 30 lv-analysis-out: {n_0,m_1,alpha_5,dx_9,dy_10,i_11,j_12,xx_13}
      yy = ((int )(- 1.0 + dy * (j - 1)));
// 30 lv-analysis-in : {n_0,m_1,alpha_5,dx_9,dy_10,i_11,j_12,xx_13,yy_14}
// 31 lv-analysis-out: {n_0,m_1,alpha_5,dx_9,dy_10,i_11,j_12,xx_13,yy_14}
      u[i][j] = 0.0;
// 31 lv-analysis-in : {n_0,m_1,alpha_5,dx_9,dy_10,i_11,j_12,xx_13,yy_14}
// 32 lv-analysis-out: {n_0,m_1,alpha_5,dx_9,dy_10,i_11,j_12,xx_13,yy_14}
      f[i][j] = - 1.0 * alpha * (1.0 - (xx * xx)) * (1.0 - (yy * yy)) - 2.0 * (1.0 - (xx * xx)) - 2.0 * (1.0 - (yy * yy));
// 32 lv-analysis-in : {n_0,m_1,alpha_5,dx_9,dy_10,i_11,j_12}
    }
// 24 lv-analysis-in : {}
// 20 lv-analysis-in : {}
}
// 11 lv-analysis-in : {n_0,m_1,alpha_5}
