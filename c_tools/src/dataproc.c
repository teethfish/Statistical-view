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

double *u_deficit;
double *v_deficit;
double *w_deficit;
double *k_cross_vel;
double *dissipation_cross_vel;

double *mean_fux;
double *mean_fuy;
double *mean_fuz;
double *mean_fvx;
double *mean_fvy;
double *mean_fvz;
double *mean_fwx;
double *mean_fwy;
double *mean_fwz;

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

  mean_fux = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  mean_fuy = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  mean_fuz = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  mean_fvx = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  mean_fvy = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  mean_fvz = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  mean_fwx = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  mean_fwy = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  mean_fwz = (double*) malloc(dom.Gcc.s3 * sizeof(double));

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

  for (int i = 0; i < dom.Gcc.s3; i++) {
    mean_fux[i] = 0.0;
    mean_fuy[i] = 0.0;
    mean_fuz[i] = 0.0;
    mean_fvx[i] = 0.0;
    mean_fvy[i] = 0.0;
    mean_fvz[i] = 0.0;
    mean_fwx[i] = 0.0;
    mean_fwy[i] = 0.0;
    mean_fwz[i] = 0.0;
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

  free(mean_fux);
  free(mean_fuy);
  free(mean_fuz);
  free(mean_fvx);
  free(mean_fvy);
  free(mean_fvz);
  free(mean_fwx);
  free(mean_fwy);
  free(mean_fwz);

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

void get_fluctuation_dissipation(double *dissipation_3d)
{
  for (int i = 0; i < dom.Gcc.s3; i++) {
    fux[i] -= mean_fux[i];
    fuy[i] -= mean_fuy[i];
    fuz[i] -= mean_fuz[i];
    fvx[i] -= mean_fvx[i];
    fvy[i] -= mean_fvy[i];
    fvz[i] -= mean_fvz[i];
    fwx[i] -= mean_fwx[i];
    fwy[i] -= mean_fwy[i];
    fwz[i] -= mean_fwz[i];
  }
  get_dissipation(dissipation_3d, nu, fux, fuy, fuz, fvx, fvy, fvz, fwx, fwy, fwz);
}
    
void get_mean_dissipation(void)
{
  double *mean_dissipation_3d;
  mean_dissipation_3d = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  get_dissipation(mean_dissipation_3d, nu, mean_fux, mean_fuy, mean_fuz, mean_fvx, mean_fvy, mean_fvz, mean_fwx, mean_fwy, mean_fwz);
  
  double *mean_dissipation_2d;
  mean_dissipation_2d = (double*) malloc(layer * dom.Gcc.s2 * sizeof(double));
  extract_surface_z(zpos, layer, mean_dissipation_3d, mean_dissipation_3d, mean_dissipation_3d, mean_dissipation_2d, mean_dissipation_2d, mean_dissipation_2d);

  double *mean_dissipation_cross_vel;
  mean_dissipation_cross_vel = (double*) malloc(layer * dom.Gcc.s1 * sizeof(double));  
  for (int k = 0; k < layer; k++) {
    for (int i = 0; i < dom.Gcc.s1; i++) {
      mean_dissipation_cross_vel[i+k*dom.Gcc.s1] = 0.0;
    }
  }
  extract_line_data_from_2d(mean_dissipation_2d, mean_dissipation_cross_vel, 1.0);

  record_2d_wake_analysis_vel_cross("mean_dissipation_cross_plane", mean_dissipation_cross_vel);

  free(mean_dissipation_3d);
  free(mean_dissipation_2d);
  free(mean_dissipation_cross_vel);
} 

void get_dissipation(double *dissipation_3d, double nu, double *dudx, double *dudy, double *dudz, double *dvdx, double *dvdy, double *dvdz, double *dwdx, double *dwdy, double *dwdz)
{
  double *tmp1;
  tmp1 = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  double *tmp2;
  tmp2 = (double*) malloc(dom.Gcc.s3 * sizeof(double));

  // 2(dudx^2 + dvdy^2 + dwdz^2)
  multiply(dom.Gcc.s3, dudx, dudx, tmp1);
  multiply(dom.Gcc.s3, dvdy, dvdy, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, dwdz, dwdz, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  for (int i = 0; i < dom.Gcc.s3; i++) {
    tmp2[i] = 2.0;
  }
  multiply(dom.Gcc.s3, tmp1, tmp2, tmp1);
 
  // dudy^2 + dvdx^2 + 2dudy*dvdx
  multiply(dom.Gcc.s3, dudy, dudy, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, dvdx, dvdx, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, dudy, dvdx, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);
  multiply(dom.Gcc.s3, dudy, dvdx, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  // dudz^2 + dwdx^2 + 2*dudz*dwdx
  multiply(dom.Gcc.s3, dudz, dudz, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, dwdx, dwdx, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, dudz, dwdx, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);
  multiply(dom.Gcc.s3, dudz, dwdx, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  // dvdz^2 + dwdy^2 + 2dvdz*dwdy
  multiply(dom.Gcc.s3, dvdz, dvdz, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, dwdy, dwdy, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, dvdz, dwdy, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);
  multiply(dom.Gcc.s3, dvdz, dwdy, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  for (int i = 0; i < dom.Gcc.s3; i++) {
    dissipation_3d[i] = tmp1[i]*nu;
  }

  free(tmp1);
  free(tmp2);
}

void calculate_mean_vel_gradient(void)
{
  double tmp = 1./nFiles;
  for (int i = 0; i < dom.Gcc.s3; i++) {
    mean_fux[i] += fux[i]*tmp; 
    mean_fuy[i] += fuy[i]*tmp;
    mean_fuz[i] += fuz[i]*tmp;
    mean_fvx[i] += fvx[i]*tmp;
    mean_fvy[i] += fvy[i]*tmp;
    mean_fvz[i] += fvz[i]*tmp;
    mean_fwx[i] += fwx[i]*tmp;
    mean_fwy[i] += fwy[i]*tmp;
    mean_fwz[i] += fwz[i]*tmp;
  }
}  


void malloc_2d_wake_analysis(void)
{
  u_deficit = (double*) malloc(layer * sizeof(double));
  v_deficit = (double*) malloc(layer * sizeof(double));
  w_deficit = (double*) malloc(layer * sizeof(double));
  k_cross_vel = (double*) malloc(layer * dom.Gcc.s1 * sizeof(double));
  dissipation_cross_vel = (double*) malloc(layer * dom.Gcc.s1 * sizeof(double));

  // init the arrays
  for (int k = 0; k < layer; k++) {
    u_deficit[k] = 0.0;
    v_deficit[k] = 0.0;
    w_deficit[k] = 0.0;
    for (int i = 0; i < dom.Gcc.s1; i++) {
      k_cross_vel[i+k*dom.Gcc.s1] = 0.0;
      dissipation_cross_vel[i+k*dom.Gcc.s1] = 0.0;
    }
  }

}

void free_2d_wake_analysis(void)
{
  free(u_deficit);
  free(v_deficit);
  free(w_deficit);

  free(k_cross_vel);
  free(dissipation_cross_vel);
}


void record_2d_wake_analysis_vel_deficit(char *name, double *fx, double *fy, double *fz)
{
  // create the file for velocity_deficit
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/dataproc/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "w");
  if(rec == NULL) {
    fprintf(stderr, "Could not open file %s\n", name);
    exit(EXIT_FAILURE);
  }

  fprintf(rec, "%-14s", "layer");
  fprintf(rec, "%-10s", "u_deficit");
  fprintf(rec, "%-10s", "v_deficit");
  fprintf(rec, "%-10s", "w_deficit");

  fprintf(rec, "\n");
  for (int k = 0; k < layer; k++) {
    fprintf(rec, "%-15f", zpos[k]);
    fprintf(rec, "%-15f", fx[k]);
    fprintf(rec, "%-15f", fy[k]);
    fprintf(rec, "%-15f", fz[k]);
    fprintf(rec, "\n");
  }
}

void record_2d_wake_analysis_vel_cross(char *name, double *f)
{
  // create the file for velocity_deficit
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/dataproc/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "w");
  if(rec == NULL) {
    fprintf(stderr, "Could not open file %s\n", name);
    exit(EXIT_FAILURE);
  }
  for (int k = 0; k < layer; k++) {
    fprintf(rec, "%s_%-9d", "layer",k);
  }
  fprintf(rec, "\n");
  for (int i = 0; i < dom.Gcc.s1; i++) {
    for (int k = 0; k < layer; k++) {
      fprintf(rec, "%-15f", f[i+k*dom.Gcc.s1]);
    }
    fprintf(rec, "\n");
  }
}

