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

/**** FUNCTIONS ****/

void analyze_3d(char *name);
void analyze_pdf_3d(int M, int Ns, int Ne);

void variable_in_part_surround_cage(double *f, double rs, double re, char *name_each, char *name_total);
/*
 *  * calculate the averaged scarlar within a shell near each particle
 *  * rs is the inner shell radius normalized by particle radius
 *  * re is the outer shell radius normalized by particle radius
 *  * name_each will record the mean value for each particle into the file
 *  * name_total will record the mean of all particles into the file
 */



#endif




