/************************************************************************
*
*
*
*
*
*
*************************************************************************/
#ifndef _3D_H
#define _3D_H

#include "main.h"
#include "cgns_reader.h"
#include "dataproc.h"
#include "tool.h"
#include <cgnslib.h>
#include <dirent.h>

/**** VARIABLES ****/
/*extern double *fux;
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
*/
/**** FUNCTIONS ****/

void analyze_3d(char *name);
void analyze_pdf_3d(int M, int Ns, int Ne);

#endif




