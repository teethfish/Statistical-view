/************************************************************************
*
*
*
*
*
*
*************************************************************************/

#include "main.h"
#include "3d.h"

/**** FUNCTIONS ****/
void dissipation_3d(double nu);

void mean_velocity(void);

void rms_velocity(void);

void record_3d_init(char *name);

void record_3d(char *name);

/**** VARIABLES ****/
/*double *fux;
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
*/
double dissipation;
double u_mean;
double v_mean;
double w_mean;
double u_rms;
double v_rms;
double w_rms;

/*void malloc_3d(void)
{
  // Init arrays for vorticity
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

}*/

/*void free_3d(void)
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
}*/
  
void analyze_3d(char *name)
{
  dissipation_3d(nu);
  printf("dissipation is %f\n", dissipation);

  //calculate mean velocity in the entire domain
  u_mean = mean(dom.Gcc.s3, uf);
  v_mean = mean(dom.Gcc.s3, vf);
  w_mean = mean(dom.Gcc.s3, wf);

  //calculate rms velocity in the entire domain
  rms_velocity();

  //write output file
  record_3d(name);
}
    
/*void calculate_gradient(void)
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
}*/

/*void vorticity(void)
{
  sub(dom.Gcc.s3, fwy, fvz, wx);
  sub(dom.Gcc.s3, fuz, fwx, wy);
  sub(dom.Gcc.s3, fvx, fuy, wz);
}*/

void dissipation_3d(double nu)
{
  dissipation = 0.0;
  
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

  dissipation = mean(dom.Gcc.s3, tmp1);
  dissipation *= nu;

  free(tmp1);
  free(tmp2);
}
 
void mean_velocity(void)
{
  u_mean = 0.0;
  v_mean = 0.0;
  w_mean = 0.0;

  u_mean = mean(dom.Gcc.s3, uf);
  v_mean = mean(dom.Gcc.s3, vf);
  w_mean = mean(dom.Gcc.s3, wf);
}

void rms_velocity(void)
{
  u_rms = 0.0;
  v_rms = 0.0;
  w_rms = 0.0;
  
  double *tmp;
  tmp = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  
  multiply(dom.Gcc.s3, uf, uf, tmp);
  u_rms = mean(dom.Gcc.s3, tmp);
  u_rms = sqrt(u_rms);
  
  multiply(dom.Gcc.s3, vf, vf, tmp);
  v_rms = mean(dom.Gcc.s3, tmp);
  v_rms = sqrt(v_rms);

  multiply(dom.Gcc.s3, wf, wf, tmp);
  w_rms = mean(dom.Gcc.s3, tmp);
  w_rms = sqrt(w_rms);

  free(tmp);
}


void record_3d_init(char *name)
{
  // create the file
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/record/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "w");
  if(rec == NULL) {
    fprintf(stderr, "Could not open file %s\n", name);
    exit(EXIT_FAILURE);
  }

 fprintf(rec, "%-15s", "time");
 fprintf(rec, "%-15s", "dissipation");
 fprintf(rec, "%-15s", "u_mean");
 fprintf(rec, "%-15s", "v_mean");
 fprintf(rec, "%-15s", "w_mean");
 fprintf(rec, "%-15s", "u_rms");
 fprintf(rec, "%-15s", "v_rms");
 fprintf(rec, "%-15s", "w_rms");

  // close the file
  fclose(rec);
}

void record_3d(char *name)
{
  // open the file
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/record/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "r+");
  if(rec == NULL) {
    record_3d_init(name);
    rec = fopen(path, "r+");
  }

  // move to the end of the file
  fseek(rec, 0, SEEK_END);

  fprintf(rec, "\n");
  fprintf(rec, "%-15d", tt);
  fprintf(rec, "%-15f", dissipation);
  fprintf(rec, "%-15f", u_mean);
  fprintf(rec, "%-15f", v_mean);
  fprintf(rec, "%-15f", w_mean);
  // only one particle
  fprintf(rec, "%-15f", u_rms);
  fprintf(rec, "%-15f", v_rms);
  fprintf(rec, "%-15f", w_rms);

  // close the file
  fclose(rec);
}
