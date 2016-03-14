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
extern double dissipation;
extern double u_mean;
extern double v_mean;
extern double w_mean;
extern double u_rms;
extern double v_rms;
extern double w_rms;
*/
/**** FUNCTIONS ****/
void malloc_3d(void);
void free_3d(void);

void analyze_3d(char *name);
void calculate_gradient(void);
void vorticity(void);

void dissipation_3d(double nu);
void mean_velocity(void);

void rms_velocity(void);
void record_3d_init(char *name);
void record_3d(char *name);






