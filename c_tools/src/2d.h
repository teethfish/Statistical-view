/************************************************************************
*
*
*
*
*
*
*************************************************************************/
#ifndef _2D_H
#define _2D_H

#include "main.h"
#include "cgns_reader.h"
#include "dataproc.h"
#include "tool.h"
#include <cgnslib.h>
#include <dirent.h>

/**** VARIABLES ****/
//min and max
/*typedef struct min_max_struct {
  double min;
  double max;
} min_max_struct;

extern min_max_struct *g_u;
extern min_max_struct *g_v;
extern min_max_struct *g_w;
*/

/**** FUNCTIONS ****/

void analyze_2d(char *name);

void analyze_pdf_2d(int M, int Ns, int Ne);

void extract_surface_z(double *zpos, int layer, double *uf, double *vf, double *wf, double *u_2d, double *v_2d, double *w_2d);

void extract_line_data_from_2d(double *f, double *f_cross_vel, double avg_factor);
/*
 *  * extract data on a line from a 2d surface
 *  * f is the array to store 2d data, length layer*in*jn
 *  * f_cross_vel store the line data from a 2d plane, length in*layer
 */

void extract_point_data_from_2d(double *f, double *f_deficit, double avg_factor);
/*
 *  * extract data on a point from 2d surface
 *  * f is the array to store 2d data, length layer*in*jn
 *  * f_deficit is the array to store point data from a 2d plane, length layer
 */    


#endif
