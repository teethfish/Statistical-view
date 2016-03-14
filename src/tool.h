/************************************************************************
*
*
*
*
*
*
*************************************************************************/

/**** VARIABLES ****/

/**** FUNCTIONS ****/
void gradient_x(double *gradx, double *f, double dx,int in, int jn, int kn);

void gradient_y(double *grady, double *f, double dy,int in, int jn, int kn);

void gradient_z(double *gradz, double *f, double dz,int in, int jn, int kn);

//void vorticity(double *wx, double *wy, double *wz);

double mean(int N, double *f);
void add(int N, double *fx, double *fy, double *f_add);
void sub(int N, double *fx, double *fy, double *f_sub);
void multiply(int N, double *fx, double *fy, double *f_multi);
void divide(int N, double *fx, double *fy, double *f_divide);

