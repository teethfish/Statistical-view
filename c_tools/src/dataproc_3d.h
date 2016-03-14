/************************************************************************
*
*
*
*
*
*
*************************************************************************/

#include "main.h"
#include "cgns_reader.h"
#include "tool.h"
#include <cgnslib.h>
#include <dirent.h>

/**** VARIABLES ****/
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

/**** FUNCTIONS ****/
void malloc_3d(void);

void free_3d(void);

void calculate_gradient(void);

void vorticity(void);

void analyze_3d(char *name);






