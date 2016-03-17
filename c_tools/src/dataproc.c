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
min_max_struct *g_u;
min_max_struct *g_v;
min_max_struct *g_w;

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

  for (int i = 0; i < dom.Gcc.s3; i++) {
    fux[i] = 0.0;
    fuy[i] = 0.0;
  }

  // min and max value on each plane
  g_u = (min_max_struct*) malloc(layer * sizeof(min_max_struct));
  g_v = (min_max_struct*) malloc(layer * sizeof(min_max_struct));
  g_w = (min_max_struct*) malloc(layer * sizeof(min_max_struct));

  for (int k = 0; k < layer; k++) {
    g_u[k].min = 0.0;
    g_u[k].max = 0.0;
    g_v[k].min = 0.0;
    g_v[k].max = 0.0;
    g_w[k].min = 0.0;
    g_w[k].max = 0.0;
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

  free(g_u);
  free(g_v);
  free(g_w);
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









