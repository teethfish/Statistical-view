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
void dissipation_in_cage(void);
/*
 * calculate the averaged dissipation within a shell near each particle
 */

void mean_velocity(void);

void rms_velocity(void);

void record_3d_init(char *name);

void record_3d(char *name);

void compute_pdf_3d(double *f, int N, int M, min_max_struct *g_f, double *r_f, double *pdf_f, double divider);
/*
 * f is the data array
 * N is the length of the data array
 * M is the number of bins for pdf
 * g_f stores the min and max number of f
 * r_f stores the returned bin information
 * pdf_f stores the pdf of f
 * divider stores the percentage of the dataset f
*/

int phase_change(void);
/*
 * omit the u,v,w and wx,wy,wz within a particle
 * returns the number of fluid cells, which is the meaningfull length of u v w wx wy wz
 */

void store_min_max_3d(double *f, min_max_struct *g_f);

void record_3d_scalar(char *name, int N, double *r, double *f);

void record_evolution_3d_init(char *name, int N);

void record_evolution_3d(char *name, int N, double *f);

/**** VARIABLES ****/
double dissipation_overall;
double u_mean;
double v_mean;
double w_mean;
double u_rms;
double v_rms;
double w_rms;
int fluid_cell; 
void analyze_3d(char *name)
{
  dissipation_in_cage();
  dissipation_overall = 0.0;

  //dissipation_overall = mean(dom.Gcc.s3, dissipation_3d);
  //printf("dissipation is %f\n", dissipation_overall);

  /*//calculate mean velocity in the entire domain
  fluid_cell = phase_change();

  u_mean = mean(fluid_cell, uf);
  v_mean = mean(fluid_cell, vf);
  w_mean = mean(fluid_cell, wf);

  layer = 1;
  store_min_max_3d(uf, g_u);
  store_min_max_3d(vf, g_v);
  store_min_max_3d(wf, g_w);
  store_min_max_3d(wx, g_wx);
  store_min_max_3d(wy, g_wy);
  store_min_max_3d(wz, g_wz);

  //calculate rms velocity in the entire domain
  rms_velocity();
  */
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

    fluid_cell = phase_change();
  
    compute_pdf_3d(uf, fluid_cell, M, g_u, r_u, pdf_u, divider);
    compute_pdf_3d(vf, fluid_cell, M, g_v, r_v, pdf_v, divider);
    compute_pdf_3d(wf, fluid_cell, M, g_w, r_w, pdf_w, divider);
    compute_pdf_3d(wx, fluid_cell, M, g_wx, r_wx, pdf_wx, divider);
    compute_pdf_3d(wy, fluid_cell, M, g_wy, r_wy, pdf_wy, divider);
    compute_pdf_3d(wz, fluid_cell, M, g_wz, r_wz, pdf_wz, divider);
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

void compute_pdf_3d(double *f, int N, int M, min_max_struct *g_f, double *r_f, double *pdf_f, double divider)
{
  double *tmp;
  tmp = (double*) malloc(M * sizeof(double));

  get_pdf(f, N, g_f[0].min, g_f[0].max, M, r_f, tmp);
  for (int i = 0; i < M; i++) {
    pdf_f[i] += divider*tmp[i];
  }
}


int phase_change(void)
{
  int tmp;
  int j = 0;
  for (int i = 0; i < dom.Gcc.s3; i++) {
    tmp = phase[i] < 0;
    if(tmp == 1){ 
      uf[j] = uf[i];
      vf[j] = vf[i];
      wf[j] = wf[i];
      wx[j] = wx[i];
      wy[j] = wx[i];
      wz[j] = wz[i];
      j = j+tmp;
    }
  }
 return j;
}  

void store_min_max_3d(double *f, min_max_struct *g_f)
{
  double tmp = 0.0;
  tmp = get_min(f, dom.Gcc.s3);
  if(tmp < g_f[0].min) g_f[0].min = tmp;
  tmp = get_max(f, dom.Gcc.s3);
  if(tmp > g_f[0].max) g_f[0].max = tmp;
}

void dissipation_in_cage(void)
{

  // calculate dissipation in the cage for each particle 
  double *dissipation_each_part;
  dissipation_each_part = (double*) malloc(nparts * sizeof(double));
  double *ncells;
  ncells = (double *) malloc(nparts * sizeof(double));
 
  double x, y, z, a, d;
  double cx1, cx2, cy1, cy2, cz1, cz2;

  for (int i = 0; i < nparts; i++) {
    dissipation_each_part[i] = 0.0;
    ncells[i] = 0;
  }

  for (int i = 0; i < dom.Gcc.in; i++) {
    x = dom.xs + (i + 0.5)*dom.dx;
    for (int j = 0; j < dom.Gcc.jn; j++) {
      y = dom.ys + (j + 0.5)*dom.dy;
      for (int k = 0; k < dom.Gcc.kn; k++) {
        z = dom.zs + (k + 0.5)*dom.dz; 
        for (int p = 0; p < nparts; p++) {
          a = parts[p].r;
          d = 1.15*parts[p].r;

          cx1 = parts[p].x;
					if(parts[p].x - d < dom.xs) cx1 = parts[p].x + dom.xl;
          if(parts[p].x + d > dom.xe) cx1 = parts[p].x - dom.xl;
          cx2 = parts[p].x;
       
          cy1 = parts[p].y;
          if(parts[p].y - d < dom.ys) cy1 = parts[p].y + dom.yl;
          if(parts[p].y + d > dom.ye) cy1 = parts[p].y - dom.yl;
          cy2 = parts[p].y;      
    
          cz1 = parts[p].z;
          if(parts[p].z - d < dom.zs) cz1 = parts[p].z + dom.zl;
          if(parts[p].z + d > dom.ze) cz1 = parts[p].z - dom.zl;
          cz2 = parts[p].z;

          int check_1 = ((x - cx1)*(x-cx1) + (y - cy1)*(y-cy1) + (z - cz1)*(z-cz1) - d*d < 0) && ((x - cx1)*(x-cx1) + (y - cy1)*(y-cy1) + (z - cz1)*(z-cz1) - a*a > 0);
          int check_2 = ((x - cx1)*(x-cx1) + (y - cy1)*(y-cy1) + (z - cz2)*(z-cz2) - d*d < 0) && ((x - cx1)*(x-cx1) + (y - cy1)*(y-cy1) + (z - cz2)*(z-cz2) - a*a > 0);
          int check_3 = ((x - cx1)*(x-cx1) + (y - cy2)*(y-cy2) + (z - cz1)*(z-cz1) - d*d < 0) && ((x - cx1)*(x-cx1) + (y - cy2)*(y-cy2) + (z - cz1)*(z-cz1) - a*a > 0);
          int check_4 = ((x - cx1)*(x-cx1) + (y - cy2)*(y-cy2) + (z - cz2)*(z-cz2) - d*d < 0) && ((x - cx1)*(x-cx1) + (y - cy2)*(y-cy2) + (z - cz2)*(z-cz2) - a*a > 0);
          int check_5 = ((x - cx2)*(x-cx2) + (y - cy1)*(y-cy1) + (z - cz1)*(z-cz1) - d*d < 0) && ((x - cx2)*(x-cx2) + (y - cy1)*(y-cy1) + (z - cz1)*(z-cz1) - a*a > 0);
          int check_6 = ((x - cx2)*(x-cx2) + (y - cy1)*(y-cy1) + (z - cz2)*(z-cz2) - d*d < 0) && ((x - cx2)*(x-cx2) + (y - cy1)*(y-cy1) + (z - cz2)*(z-cz2) - a*a > 0);
          int check_7 = ((x - cx2)*(x-cx2) + (y - cy2)*(y-cy2) + (z - cz1)*(z-cz1) - d*d < 0) && ((x - cx2)*(x-cx2) + (y - cy2)*(y-cy2) + (z - cz1)*(z-cz1) - a*a > 0);
          int check_8 = ((x - cx2)*(x-cx2) + (y - cy2)*(y-cy2) + (z - cz2)*(z-cz2) - d*d < 0) && ((x - cx2)*(x-cx2) + (y - cy2)*(y-cy2) + (z - cz2)*(z-cz2) - a*a > 0);	
          int check = (check_1 + check_2 + check_3 + check_4 + check_5 + check_6 + check_7 + check_8)!= 0;
          dissipation_each_part[p] += dissipation_3d[i + j*dom.Gcc.s1 + k*dom.Gcc.s2]*check;
          ncells[p] += 1.0*check; 
        }
      }
    }
  } 
   
  double *dissipation_total;
  dissipation_total = (double*) malloc(1 * sizeof(double));
  for (int i = 0; i < nparts; i++) {
    dissipation_each_part[i] *= nu;
    dissipation_each_part[i] = dissipation_each_part[i] / ncells[i];
    dissipation_total[0] += dissipation_each_part[i]/(double) nparts;
  }
  record_evolution_3d("dissipation_each_particle", nparts, dissipation_each_part);
  record_evolution_3d("dissipation_total_particle", 1, dissipation_total);  
  record_evolution_3d("ncells_each_particle", nparts, ncells); 
 
  free(dissipation_each_part);
  free(dissipation_total);
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
 
  double mean_value = mean(fluid_cell, uf);
  for (int i = 0; i < fluid_cell; i++) {
    tmp[i] = uf[i] - mean_value;
  } 
  multiply(fluid_cell, tmp, tmp, tmp);
  u_rms = mean(fluid_cell, tmp);
  u_rms = sqrt(u_rms);
 
  mean_value = mean(fluid_cell, vf);
  for (int i = 0; i < fluid_cell; i++) {
    tmp[i] = vf[i] - mean_value;
  }
  multiply(fluid_cell, tmp, tmp, tmp);
  v_rms = mean(fluid_cell, tmp);
  v_rms = sqrt(v_rms);

  mean_value = mean(fluid_cell, wf);
  for (int i = 0; i < fluid_cell; i++) {
    tmp[i] = wf[i] - mean_value;
  }
  multiply(fluid_cell, tmp, tmp, tmp);
  w_rms = mean(fluid_cell, tmp);
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

void record_evolution_3d_init(char *name, int N)
{
  // create the file for each surface
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/dataproc/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "w");
  if(rec == NULL) {
    fprintf(stderr, "Could not open file %s\n", name);
    exit(EXIT_FAILURE);
  }
  
  fprintf(rec, "%-15s", "time");
  for (int i = 0; i < N; i++) {
    fprintf(rec, "%s_%-10d ", "part",i);
  } 

  // close the file
  fclose(rec);
}  

void record_evolution_3d(char *name, int N, double *f)
{
  // open the file
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/dataproc/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "r+");
  if(rec == NULL) {
    record_evolution_3d_init(name, N);
    rec = fopen(path, "r+");
  }

  // move to the end of the file
  fseek(rec, 0, SEEK_END);

  fprintf(rec, "\n");
  fprintf(rec, "%-15d", tt);
  for (int i = 0; i < N; i++) {
    fprintf(rec, "%-15f ", f[i]);
  }
  // close the file
  fclose(rec);
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
  fprintf(rec, "%-15f", dissipation_overall);
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
