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

void compute_pdf_3d(double *f, int M, min_max_struct *g_f, double *r_f, double *pdf_f, double divider);

void phase_change(void);

void store_min_max_3d(double *f, min_max_struct *g_f);

void record_3d_scalar(char *name, int N, double *r, double *f);
/**** VARIABLES ****/
double dissipation;
double u_mean;
double v_mean;
double w_mean;
double u_rms;
double v_rms;
double w_rms;
 
void analyze_3d(char *name)
{
  dissipation_3d(nu);
  printf("dissipation is %f\n", dissipation);

  //calculate mean velocity in the entire domain
  phase_change();

  u_mean = mean(dom.Gcc.s3, uf);
  v_mean = mean(dom.Gcc.s3, vf);
  w_mean = mean(dom.Gcc.s3, wf);

  layer = 1;
  store_min_max_3d(uf, g_u);
  store_min_max_3d(vf, g_v);
  store_min_max_3d(wf, g_w);
  store_min_max_3d(wx, g_wx);
  store_min_max_3d(wy, g_wy);
  store_min_max_3d(wz, g_wz);

  //calculate rms velocity in the entire domain
  rms_velocity();

  //write output file
  record_3d(name);
}

void analyze_pdf_3d(int M, int Ns, int Ne)
{
  layer = 1;
  double *r_u;
  double *r_v;
  double *r_w;
  r_u = (double*) malloc(M *layer * sizeof(double));
  r_v = (double*) malloc(M *layer * sizeof(double));
  r_w = (double*) malloc(M *layer * sizeof(double));
  double *r_wx;
  double *r_wy;
  double *r_wz;
  r_wx = (double*) malloc(M *layer * sizeof(double));
  r_wy = (double*) malloc(M *layer * sizeof(double));
  r_wz = (double*) malloc(M *layer * sizeof(double));

  double *pdf_u;
  double *pdf_v;
  double *pdf_w;
  pdf_u = (double*) malloc(M * layer * sizeof(double));
  pdf_v = (double*) malloc(M * layer * sizeof(double));
  pdf_w = (double*) malloc(M * layer * sizeof(double));
  double *pdf_wx;
  double *pdf_wy;
  double *pdf_wz;
  pdf_wx = (double*) malloc(M * layer * sizeof(double));
  pdf_wy = (double*) malloc(M * layer * sizeof(double));
  pdf_wz = (double*) malloc(M * layer * sizeof(double));

  //initialize the pdf array to be zeros
  for (int i = 0; i < M*layer; i++) {
    pdf_u[i] = 0.0;
    pdf_v[i] = 0.0;
    pdf_w[i] = 0.0;
    pdf_wx[i] = 0.0;
    pdf_wy[i] = 0.0;
    pdf_wz[i] = 0.0;
  }

  for (int k = 0; k < layer; k++) {
    printf("g_u.min is %f\n", g_u[k].min);
    printf("g_u.max is %f\n", g_u[k].max);
    printf("g_v.min is %f\n", g_v[k].min);
    printf("g_v.max is %f\n", g_v[k].max);
    printf("g_w.min is %f\n", g_w[k].min);
    printf("g_w.max is %f\n", g_w[k].max);
    printf("g_wx.min is %f\n", g_wx[k].min);
    printf("g_wx.max is %f\n", g_wx[k].max);
    printf("g_wy.min is %f\n", g_wy[k].min);
    printf("g_wy.max is %f\n", g_wy[k].max);
    printf("g_wz.min is %f\n", g_wz[k].min);
    printf("g_wz.max is %f\n", g_wz[k].max);    
  }

  double divider = 1.0/(Ne-Ns);
  for (int i = Ns; i < Ne; i++) {
    tt = i;
    cgns_fill_flow();

    calculate_gradient();
    vorticity();

    phase_change();
     
    compute_pdf_3d(uf, M, g_u, r_u, pdf_u, divider);
    compute_pdf_3d(vf, M, g_v, r_v, pdf_v, divider);
    compute_pdf_3d(wf, M, g_w, r_w, pdf_w, divider);
    compute_pdf_3d(wx, M, g_wx, r_wx, pdf_wx, divider);
    compute_pdf_3d(wy, M, g_wy, r_wy, pdf_wy, divider);
    compute_pdf_3d(wz, M, g_wz, r_wz, pdf_wz, divider);
  }
  printf("finished time iteration for pdf\n");

  record_3d_scalar("pdf_u.dat", M, r_u, pdf_u);
  record_3d_scalar("pdf_v.dat", M, r_v, pdf_v);
  record_3d_scalar("pdf_w.dat", M, r_w, pdf_w);
  record_3d_scalar("pdf_wx.dat", M, r_wx, pdf_wx);
  record_3d_scalar("pdf_wy.dat", M, r_wy, pdf_wy);
  record_3d_scalar("pdf_wz.dat", M, r_wz, pdf_wz);

  free(r_u);
  free(r_v);
  free(r_w);
  free(pdf_u);
  free(pdf_v);
  free(pdf_w);
  free(r_wx);
  free(r_wy);
  free(r_wz);
  free(pdf_wx);
  free(pdf_wy);
  free(pdf_wz);
}

void compute_pdf_3d(double *f, int M, min_max_struct *g_f, double *r_f, double *pdf_f, double divider)
{
  double *tmp;
  tmp = (double*) malloc(M * sizeof(double));

  get_pdf(f, dom.Gcc.s3, g_f[0].min, g_f[0].max, M, r_f, tmp);
  for (int i = 0; i < M; i++) {
    pdf_f[i] += divider*tmp[i];
  }
}


void phase_change(void)
{
  int tmp;
  //int fluid_cell = 0;
  for (int i = 0; i < dom.Gcc.s3; i++) {
    tmp = phase[i] < 0;
    //fluid_cell += tmp;
    uf[i] = tmp*uf[i];
    vf[i] = tmp*vf[i];
    wf[i] = tmp*wf[i];
    wx[i] = tmp*wx[i];
    wy[i] = tmp*wx[i];
    wz[i] = tmp*wz[i];
  }
 //printf("fluid cell is %d\n", fluid_cell);
 //return fluid_cell;
}  

void store_min_max_3d(double *f, min_max_struct *g_f)
{
  double tmp = 0.0;
  tmp = get_min(f, dom.Gcc.s3);
  if(tmp < g_f[0].min) g_f[0].min = tmp;
  tmp = get_max(f, dom.Gcc.s3);
  if(tmp > g_f[0].max) g_f[0].max = tmp;
}
    
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
  sprintf(path, "%s/dataproc/%s", ROOT_DIR, name);
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

void record_3d_scalar(char *name, int N, double *r, double *f)
{
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/dataproc/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "w");
  if(rec == NULL) {
    fprintf(stderr, "Could not open file %s\n", name);
    exit(EXIT_FAILURE);
  }

  fprintf(rec, "%-15s", "pdf_x");
  fprintf(rec, "%-15s", "pdf_y");
  fprintf(rec, "\n");

  for (int i = 0; i < N; i++) {
    fprintf(rec, "%-15f", r[i]);
    fprintf(rec, "%-15f", f[i]);
    fprintf(rec, "\n");
  }
}

void record_3d(char *name)
{
  // open the file
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/dataproc/%s", ROOT_DIR, name);
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
