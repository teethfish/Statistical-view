/************************************************************************
*
*
*
*
*
*
*************************************************************************/

#include "main.h"
#include "dataproc_2d.h"

/*****FUNCTIONS********/
void malloc_2d(void);

void free_2d(void);

void extract_surface_z(double *zpos, int layer, double *uf, double *vf, double *wf, double *u_2d, double *v_2d, double *w_2d);

void mean_vel_surface(int layer, double *u_2d, double *v_2d, double *w_2d);

void dissipation_2d(double nu);

void record_2d_init(char *name);

void record_2d(char *name);

/*** VARIABLES *******/
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
  //value on planes
  int size = layer*dom.Gcc.s2;
  u_2d = (double*) malloc(size * sizeof(double));
  v_2d = (double*) malloc(size * sizeof(double));
  w_2d = (double*) malloc(size * sizeof(double));

  //1d array
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

void analyze_2d(char *name)
{
  malloc_2d();
  // interpolate velocity from 3d to 2d
  extract_surface_z(zpos, layer, uf, vf, wf, u_2d, v_2d, w_2d);

  // calculate mean and rms velocity on each plane
  mean_vel_surface(layer, u_2d, v_2d, w_2d);

  dissipation_2d(nu);

  record_2d(name);

  free_2d();
} 

void extract_surface_z(double *zpos, int layer, double *uf, double *vf, double *wf, double *u_2d, double *v_2d, double *w_2d)
{  
  int L; 
  int H;
  double ddz = 1./dom.dz;
  double zpos_c = 0.0;
  for (int k = 0; k < layer; k++) {
    // find the nearest lower cell-center index
    L = floor((zpos[k] - dom.zs - dom.dz*0.5)*ddz) + 1;
    zpos_c = zpos[k] - dom.zs - (dom.dz/2.0 + (L - 1)*dom.dz);
    zpos_c = zpos_c*ddz;
    // find the nearest higher cell-center index
    H = L + 1;
    if (H > dom.zn) H = L;
    if (L < 1) L = 1;

    //interpolate 3d data to 2d plane
    for (int i = 0; i < dom.Gcc.in; i++) {
      for (int j = 0; j < dom.Gcc.jn; j++) {
	u_2d[i+j*dom.Gcc.s1+k*dom.Gcc.s2] = uf[i+j*dom.Gcc.s1+L*dom.Gcc.s2] + (uf[i+j*dom.Gcc.s1+H*dom.Gcc.s2] - uf[i+j*dom.Gcc.s1+L*dom.Gcc.s2])*zpos_c;
	v_2d[i+j*dom.Gcc.s1+k*dom.Gcc.s2] = vf[i+j*dom.Gcc.s1+L*dom.Gcc.s2] + (vf[i+j*dom.Gcc.s1+H*dom.Gcc.s2] - vf[i+j*dom.Gcc.s1+L*dom.Gcc.s2])*zpos_c;
 	w_2d[i+j*dom.Gcc.s1+k*dom.Gcc.s2] = wf[i+j*dom.Gcc.s1+L*dom.Gcc.s2] + (wf[i+j*dom.Gcc.s1+H*dom.Gcc.s2] - wf[i+j*dom.Gcc.s1+L*dom.Gcc.s2])*zpos_c;      
      }
    }
  }
}

void mean_vel_surface(int layer, double *u_2d, double *v_2d, double *w_2d)
{
  double *tmp;
  tmp = (double*) malloc(dom.Gcc.s2 * sizeof(double));
 
  for (int k = 0; k < layer; k++) {
    u_mean[k] = mean(dom.Gcc.s2, u_2d + k*dom.Gcc.s2);
    v_mean[k] = mean(dom.Gcc.s2, v_2d + k*dom.Gcc.s2);
    w_mean[k] = mean(dom.Gcc.s2, w_2d + k*dom.Gcc.s2);
   
    for (int i = 0; i < dom.Gcc.s2; i++) {
      tmp[i] = u_mean[k];
    }
    //calculate (u - u_mean)*(u - u_mean)
    sub(dom.Gcc.s2, u_2d + k*dom.Gcc.s2, tmp, tmp); 
    multiply(dom.Gcc.s2, tmp, tmp, tmp);
    u_rms[k] = mean(dom.Gcc.s2, tmp);
    u_rms[k] = sqrt(u_rms[k]); 

    for (int i = 0; i < dom.Gcc.s2; i++) {
      tmp[i] = v_mean[k];
    } 
    sub(dom.Gcc.s2, v_2d + k*dom.Gcc.s2, tmp, tmp);
    multiply(dom.Gcc.s2, tmp, tmp, tmp);
    v_rms[k] = mean(dom.Gcc.s2, tmp);
    v_rms[k] = sqrt(v_rms[k]);  

    for (int i = 0; i < dom.Gcc.s2; i++) {
      tmp[i] = w_mean[k];
    } 
    sub(dom.Gcc.s2, w_2d + k*dom.Gcc.s2, tmp, tmp);
    multiply(dom.Gcc.s2, tmp, tmp, tmp);
    w_rms[k] = mean(dom.Gcc.s2, tmp);
    w_rms[k] = sqrt(w_rms[k]);
  }

  free(tmp);
}

void dissipation_2d(double nu)
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

  double *tmp;
  tmp = (double*) malloc(layer*dom.Gcc.s2 * sizeof(double));

  extract_surface_z(zpos, layer, tmp1, tmp1, tmp1, tmp, tmp, tmp);

  for (int k = 0; k < layer; k++) {
    dissipation[k] = mean(dom.Gcc.s2, tmp+k*dom.Gcc.s2);
    dissipation[k] = nu * dissipation[k];
    //printf("dissipation on plane %d is %f\n", k, dissipation[k]);
  }

  free(tmp1);
  free(tmp2);
  free(tmp);
}   

void record_2d_init(char *name)
{
  // create the file for each surface
  for (int k = 0; k < layer; k++) {
    char path[FILE_NAME_SIZE] = "";
    sprintf(path, "%s/record/%s_%d", ROOT_DIR, name, k);
    FILE *rec = fopen(path, "w");
    if(rec == NULL) {
      fprintf(stderr, "Could not open file %s_%d\n", name, k);
      exit(EXIT_FAILURE);
    }
   
    fprintf(rec, "%-15s", "time");
    fprintf(rec, "%s_%-3d", "dissipation",k);
    fprintf(rec, "%s_%-8d", "u_mean",k);
    fprintf(rec, "%s_%-8d", "v_mean",k);
    fprintf(rec, "%s_%-8d", "w_mean",k);
    fprintf(rec, "%s_%-9d", "u_rms",k);
    fprintf(rec, "%s_%-9d", "v_rms",k);
    fprintf(rec, "%s_%-9d", "w_rms",k);

    // close the file
    fclose(rec);
  }
}

void record_2d(char *name)
{
  // open the file
  for (int k = 0; k < layer; k++) {
    char path[FILE_NAME_SIZE] = "";
    sprintf(path, "%s/record/%s_%d", ROOT_DIR, name, k);
    FILE *rec = fopen(path, "r+");
    if(rec == NULL) {
      record_2d_init(name);
      rec = fopen(path, "r+");
    }
   
    // move to the end of the file
    fseek(rec, 0, SEEK_END);

    fprintf(rec, "\n");
    fprintf(rec, "%-15d", tt);
    fprintf(rec, "%-15f", dissipation[k]);
    fprintf(rec, "%-15f", u_mean[k]);
    fprintf(rec, "%-15f", v_mean[k]);
    fprintf(rec, "%-15f", w_mean[k]);
    fprintf(rec, "%-15f", u_rms[k]);
    fprintf(rec, "%-15f", v_rms[k]);
    fprintf(rec, "%-15f", w_rms[k]);
   
    // close the file
    fclose(rec);
  }
}


