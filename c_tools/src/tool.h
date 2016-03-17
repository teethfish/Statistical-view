/************************************************************************
*
*
*
*
*
*
*************************************************************************/
#ifndef _TOOL_H
#define _TOOL_H

#include <stdlib.h>
#include <math.h>
#include <string.h>

/**** VARIABLES ****/

/**** FUNCTIONS ****/
void gradient_x(double *gradx, double *f, double dx,int in, int jn, int kn);

void gradient_y(double *grady, double *f, double dy,int in, int jn, int kn);

void gradient_z(double *gradz, double *f, double dz,int in, int jn, int kn);

// basic functions
double mean(int N, double *f);

void add(int N, double *fx, double *fy, double *f_add);

void sub(int N, double *fx, double *fy, double *f_sub);

void multiply(int N, double *fx, double *fy, double *f_multi);

void divide(int N, double *fx, double *fy, double *f_divide);

void copy(double *f1, double *f2, int N);


// find min and max value for a 1d array
//inline double min2(double x, double y);

double get_min(double *a, int n);

//inline double max2(double x, double y);

double get_max(double *a, int n);

//advanced function:pdf and correlation
void get_pdf(double *f, int N, double min_f, double max_f, int M, double *r, double *pdf);
/****
  * f input array
  * N number of datapoints
  * M number of bins
  * r range of datapoints in
  * pdf the output density function
***/

void correlation(int N, int lag, double *fx, double *fy, double *cxy);
/****
  * N is the total length for signal fx and fy
  * lag is the maximum lag that supposed to calculated, usually it takes as N/2
  * cxy is the cross-correlation array, with the first term is the co-variance of fx and fy; 
  *******LAG=0 WILL GIVE THE AUTO-VARIANCE OF FX AND FY
***/

#endif

