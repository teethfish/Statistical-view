/************************************************************************
*
*
*
*
*
*
*************************************************************************/

#ifndef _DATAPROC_H
#define _DATAPROC_H

#include "main.h"
#include "cgns_reader.h"
#include "tool.h"

/**** VARIABLES ****/
//min and max
typedef struct min_max_struct {
  double min;
  double max;
} min_max_struct;

extern double *fux;
extern double *fuy;
extern double *fuz;
extern double *fvx;
extern double *fvy;
extern double *fvz;
extern double *fwx;
extern double *fwy;
extern double *fwz;
extern double *wx;
extern double *wy;
extern double *wz;
extern double *dissipation_3d;

extern min_max_struct *g_u;
extern min_max_struct *g_v;
extern min_max_struct *g_w;
extern min_max_struct *g_wx;
extern min_max_struct *g_wy;
extern min_max_struct *g_wz;

extern double *u_deficit;
extern double *v_deficit;
extern double *w_deficit;
extern double *k_cross_vel;
extern double *dissipation_cross_vel;

extern double *mean_fux;
extern double *mean_fuy;
extern double *mean_fuz;
extern double *mean_fvx;
extern double *mean_fvy;
extern double *mean_fvz;
extern double *mean_fwx;
extern double *mean_fwy;
extern double *mean_fwz;

/**** FUNCTIONS ****/
void malloc_dataproc(void);

void free_dataproc(void);

void calculate_gradient(void);

void vorticity(void);

void get_dissipation(double *dissipation_3d, double nu, double *dudx, double *dudy, double *dudz, double *dvdx, double *dvdy, double *dvdz, double *dwdx, double *dwdy, double *dwdz);

void get_fluctuation_dissipation(double *dissipation_3d);

void get_mean_dissipation(void);

void calculate_mean_vel_gradient(void);

void malloc_2d_wake_analysis(void);

void free_2d_wake_analysis(void);

void record_2d_wake_analysis_vel_deficit(char *name, double *fx, double *fy, double *fz);

void record_2d_wake_analysis_vel_cross(char *name, double *f);

#endif
