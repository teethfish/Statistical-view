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

extern min_max_struct *g_u;
extern min_max_struct *g_v;
extern min_max_struct *g_w;
extern min_max_struct *g_wx;
extern min_max_struct *g_wy;
extern min_max_struct *g_wz;

/**** FUNCTIONS ****/
void malloc_dataproc(void);

void free_dataproc(void);

void calculate_gradient(void);

void vorticity(void);

#endif
