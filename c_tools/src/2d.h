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

#endif
