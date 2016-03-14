/************************************************************************
*
*
*
*
*
*
*************************************************************************/

#include "dataproc_2d.h"

int layer;
double *u_mean;
double *v_mean;
double *w_mean;
double *u_rms;
double *v_rms;
double *w_rms;
double *u_2d;
double *v_2d;
double *w_2d;
double *dissipation;

void malloc_2d(void)
{
  int size = layer*dom.Gcc.s2;
  u_2d = (double*) malloc(size * sizeof(double));
  v_2d = (double*) malloc(size * sizeof(double));
  w_2d = (double*) malloc(size * sizeof(double));

  u_mean = (double*) malloc(layer * sizeof(double));
  v_mean = (double*) malloc(layer * sizeof(double));
  w_mean = (double*) malloc(layer * sizeof(double));

  u_rms = (double*) malloc(layer * sizeof(double));
  v_rms = (double*) malloc(layer * sizeof(double));
  w_rms = (double*) malloc(layer * sizeof(double));

  dissipation = (double*) malloc(layer * sizeof(double)); 
}

void free_2d(void)
{
  free(u_mean);
  free(v_mean);
  free(w_mean);

  free(u_rms);
  free(v_rms);
  free(w_rms);
  
  free(dissipation);

  free(u_2d);
  free(v_2d);
  free(w_2d);
}  

void extract_surface_z(double *zpos, int layer, double *uf, double *vf, double *wf, double *u_2d, double *v_2d, double *w_2d)
{  
  int L; 
  int H;
  double ddz = 1./dom.dz;
  double zpos_c = 0.0;
  for(int k = 0; k < layer; k++){
    // find the nearest lower cell-center index
    L = floor((zpos[k] - dom.xe - dom.dz/2.0)*ddz);
    zpos_c = zpos[k] - (dom.dz/2.0 + L*dom.dz);
    // find the nearest higher cell-center index
    H = L + 1;
    if(H > dom.zn)H = L;
    if(L < 1)L = 1;
   
    //interpolate 3d data to 2d plane
    for(int i = 0; i < dom.Gcc.in; i++){
      for(int j = 0; j < dom.Gcc.jn; j++){
	u_2d[i+j*dom.Gcc.s1] = uf[i+j*dom.Gcc.s1+L*dom.Gcc.s2] + (uf[i+j*dom.Gcc.s1+H*dom.Gcc.s2] - uf[i+j*dom.Gcc.s1+L*dom.Gcc.s2])*zpos_c/dom.dz 
	v_2d[i+j*dom.Gcc.s1] = vf[i+j*dom.Gcc.s1+L*dom.Gcc.s2] + (vf[i+j*dom.Gcc.s1+H*dom.Gcc.s2] - vf[i+j*dom.Gcc.s1+L*dom.Gcc.s2])*zpos_c/dom.dz 
	w_2d[i+j*dom.Gcc.s1] = wf[i+j*dom.Gcc.s1+L*dom.Gcc.s2] + (wf[i+j*dom.Gcc.s1+H*dom.Gcc.s2] - wf[i+j*dom.Gcc.s1+L*dom.Gcc.s2])*zpos_c/dom.dz 
      }
    }
  }
}

void mean_vel_surface(int layer, double *u_2d, double *v_2d, double *w_2d)
{
  double *tmp;
  tmp = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  for(int k = 0; k < layer; k++){
    u_mean[k] = mean(dom.Gcc.s1, u_2d + k*dom.Gcc.s1);
    v_mean[k] = mean(dom.Gcc.s1, v_2d + k*dom.Gcc.s1);
    w_mean[k] = mean(dom.Gcc.s1, w_2d + k*dom.Gcc.s1);
    
    multiply(dom.Gcc.s1, u_2d + k*dom.Gcc.s1, u_2d + k*dom.Gcc.s1, tmp);
    u_rms[k] = mean(dom.Gcc.s1, tmp);
    
    multiply(dom.Gcc.s1, v_2d + k*dom.Gcc.s1, v_2d + k*dom.Gcc.s1, tmp);
    v_rms[k] = mean(dom.Gcc.s1, tmp);
    
    multiply(dom.Gcc.s1, w_2d + k*dom.Gcc.s1, w_2d + k*dom.Gcc.s1, tmp);
    w_rms[k] = mean(dom.Gcc.s1, tmp);

  free(tmp);
}

double calculate_gradient(void)
{
  gradient_x(uf, dx, fux);
  gradient_x(vf, dx, fvx);
  gradient_x(wf, dx, fwx);
  gradient_y(uf, dy, fuy);
  gradient_y(vf, dy, fvy);
  gradient_y(wf, dy, fwy);
  gradient_z(uf, dz, fuz);
  gradient_z(vf, dz, fvz);
  gradient_z(wf, dz, fwz);
}

void vorticity(void)
{
  // Init arrays for vorticity
  wx = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  wy = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  wz = (double*) malloc(dom.Gcc.s3 * sizeof(double));

  sub(dom.Gcc.s3, fwy, fvz, wx);
  sub(dom.Gcc.s3, fuz, fwx, wy);
  sub(dom.Gcc.s3, fvx, fuy, wz);
}

double dissipation_3d(double nu)
{

  double *tmp1;
  tmp1 = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  double *tmp2;
  tmp2 = (double*) malloc(dom.Gcc.s3 * sizeof(double));

  multiply(dom.Gcc.s3, fux, fux, tmp1);
  multiply(dom.Gcc.s3, fvx, fvx, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, fwx, fwx, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, fuy, fuy, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, fvy, fvy, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, fwy, fwy, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, fuz, fuz, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, fvz, fvz, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  multiply(dom.Gcc.s3, fwz, fwz, tmp2);
  add(dom.Gcc.s3, tmp1, tmp2, tmp1);

  extract_surface_z(double *zpos, int layer, double *uf, double *vf, double *wf, double *u_2d, double *v_2d, double *w_2d)  
 





