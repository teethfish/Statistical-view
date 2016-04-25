/************************************************************************
*
*
*
*
*
*
*************************************************************************/

#include "main.h"
#include "dataproc.h"

double *fux;
double *fuy;
double *fuz;
double *fvx;
double *fvy;
double *fvz;
double *fwx;
double *fwy;
double *fwz;
double *wx;
double *wy;
double *wz;
double *dissipation_3d;
min_max_struct *g_u;
min_max_struct *g_v;
min_max_struct *g_w;
min_max_struct *g_wx;
min_max_struct *g_wy;
min_max_struct *g_wz;

void malloc_dataproc(void)
{
  // Init arrays for velocity gradient and vorticity
  fux = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  fuy = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  fuz = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  fvx = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  fvy = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  fvz = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  fwx = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  fwy = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  fwz = (double*) malloc(dom.Gcc.s3 * sizeof(double));

  wx = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  wy = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  wz = (double*) malloc(dom.Gcc.s3 * sizeof(double));

  dissipation_3d = (double*) malloc(dom.Gcc.s3 * sizeof(double));

  // min and max value on each plane
  g_u = (min_max_struct*) malloc(layer * sizeof(min_max_struct));
  g_v = (min_max_struct*) malloc(layer * sizeof(min_max_struct));
  g_w = (min_max_struct*) malloc(layer * sizeof(min_max_struct));
  g_wx = (min_max_struct*) malloc(layer * sizeof(min_max_struct));
  g_wy = (min_max_struct*) malloc(layer * sizeof(min_max_struct));
  g_wz = (min_max_struct*) malloc(layer * sizeof(min_max_struct));


  for (int k = 0; k < layer; k++) {
    g_u[k].min = 0.0;
    g_u[k].max = 0.0;
    g_v[k].min = 0.0;
    g_v[k].max = 0.0;
    g_w[k].min = 0.0;
    g_w[k].max = 0.0;
    g_wx[k].min = 0.0;
    g_wx[k].max = 0.0;
    g_wy[k].min = 0.0;
    g_wy[k].max = 0.0;
    g_wz[k].min = 0.0;
    g_wz[k].max = 0.0;
  }

}

void free_dataproc(void)
{
  free(fux);
  free(fuy);
  free(fuz);
  free(fvx);
  free(fvy);
  free(fvz);
  free(fwx);
  free(fwy);
  free(fwz);

  free(wx);
  free(wy);
  free(wz);

  free(dissipation_3d);

  free(g_u);
  free(g_v);
  free(g_w);
  free(g_wx);
  free(g_wy);
  free(g_wz);
}

void calculate_gradient(void)
{
  gradient_x(fux, uf, dom.dx, dom.Gcc.in, dom.Gcc.jn, dom.Gcc.kn);
  gradient_x(fvx, vf, dom.dx, dom.Gcc.in, dom.Gcc.jn, dom.Gcc.kn);
  gradient_x(fwx, wf, dom.dx, dom.Gcc.in, dom.Gcc.jn, dom.Gcc.kn);
  gradient_y(fuy, uf, dom.dy, dom.Gcc.in, dom.Gcc.jn, dom.Gcc.kn);
  gradient_y(fvy, vf, dom.dy, dom.Gcc.in, dom.Gcc.jn, dom.Gcc.kn);
  gradient_y(fwy, wf, dom.dy, dom.Gcc.in, dom.Gcc.jn, dom.Gcc.kn);
  gradient_z(fuz, uf, dom.dz, dom.Gcc.in, dom.Gcc.jn, dom.Gcc.kn);
  gradient_z(fvz, vf, dom.dz, dom.Gcc.in, dom.Gcc.jn, dom.Gcc.kn);
  gradient_z(fwz, wf, dom.dz, dom.Gcc.in, dom.Gcc.jn, dom.Gcc.kn);

}

void vorticity(void)
{
  sub(dom.Gcc.s3, fwy, fvz, wx);
  sub(dom.Gcc.s3, fuz, fwx, wy);
  sub(dom.Gcc.s3, fvx, fuy, wz);
}

void get_dissipation(double nu)
{
  double *tmp1;
  tmp1 = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  double *tmp2;
  tmp2 = (double*) malloc(dom.Gcc.s3 * sizeof(double));

  // 2(dudx^2 + dvdy^2 + dwdz^2)
  multiply(dom.Gcc.s3, fux, fux, tmp1);
  multiply(dom.Gcc.s3, fvy, fvy, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, fwz, fwz, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  for (int i = 0; i < dom.Gcc.s3; i++) {
    tmp2[i] = 2.0;
  }
  multiply(dom.Gcc.s3, tmp1, tmp2, tmp1);
 
  // dudy^2 + dvdx^2 + 2dudy*dvdx
  multiply(dom.Gcc.s3, fuy, fuy, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, fvx, fvx, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, fuy, fvx, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);
  multiply(dom.Gcc.s3, fuy, fvx, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  // dudz^2 + dwdx^2 + 2*dudz*dwdx
  multiply(dom.Gcc.s3, fuz, fuz, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, fwx, fwx, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, fuz, fwx, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);
  multiply(dom.Gcc.s3, fuz, fwx, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  // dvdz^2 + dwdy^2 + 2dvdz*dwdy
  multiply(dom.Gcc.s3, fvz, fvz, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, fwy, fwy, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, fvz, fwy, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);
  multiply(dom.Gcc.s3, fvz, fwy, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  for (int i = 0; i < dom.Gcc.s3; i++) {
    dissipation_3d[i] = tmp1[i]*nu;
  }

  free(tmp1);
  free(tmp2);
}









