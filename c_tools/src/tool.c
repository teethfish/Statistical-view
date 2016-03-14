/**********************************************************************
*
*
*
***********************************************************************/

void gradient_x(double *gradx, double *f, double dx, int in, int jn, int kn)
{
  int s1 = in;
  int s2 = in * jn;
  for(int i = 1; i < in-1; i++){
    for(int j = 0; j < jn; j++){
      for(int k = 0; k < kn; k++){
	 gradx[i+j*s1+k*s2] = (f[i+j*s1+k*s2+1] - f[i+j*s1+k*s2-1])/dx/2.0; 	
	}
    }
  }
  for(int j = 0; j < jn; j++){
    for(int k = 0; k < kn; k++){
      gradx[j*s1+k*s2] = (f[1+j*s1+k*s2]-f[j*s1+k*s2])/dx;
      gradx[in-1 + j*s1+k*s2] = (f[in-1 + j*s1+k*s2] - f[in-2+j*s1+k*s2])/dx;    
    }
  }
}

void gradient_y(double *grady, double *f, double dy, int in, int jn, int kn)
{
  int s1 = in;
  int s2 = in*jn;
  for(int i = 0; i < in; i++){
    for(int j = 1; j < jn - 1; j++){
      for(int k = 0; k < kn; k++){
	grady[i+j*s1+k*s2] = (f[i+(j+1)*s1+k*s2] - f[i+(j-1)*s1+k*s2])/dy/2.0;
      }
    }
  }
 for(int i = 0; i < in; i++){
    for(int k = 0; k < kn; k++){
      grady[i+k*s2] = (f[i+s1+k*s2] - f[i+k*s2])/dy;    
      grady[i+(jn-1)*s1+k*s2] = (f[i+(jn-1)*s1+k*s2] - f[i+(jn-2)*s1+k*s2])/dy;
    }
  }
} 

void gradient_z(double *gradz, double *f, double dz, int in, int jn, int kn)
{
  int s1 = in;
  int s2 = in*jn;
  for(int i = 0; i < in; i++){
    for(int j = 0; j < jn; j++){
      for(int k = 1; k < kn-1; k++){
	gradz[i+j*s1+k*s2] = (f[i+j*s1+(k+1)*s2] - f[i+j*s1+(k-1)*s2])/dz/2.0;
      }
    }
  }

  for(int i = 0; i < in; i++){
    for(int j = 0; j < jn; j++){
       gradz[i+j*s1] = (f[i+j*s1+s2] - f[i+j*s1])/dz;
      gradz[i+j*s1+(kn-1)*s2] = (f[i+j*s1+(kn-1)*s2] - f[i+j*s1+(kn-2)*s2])/dz;
    }
  }
}
/*
void vorticity(double *wx, double *wy, double *wz)
{
  // Init arrays for vorticity
  wx = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  wy = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  wz = (double*) malloc(dom.Gcc.s3 * sizeof(double));
   
  double *fuy;
  double *fuz;
  double *fvx;
  double *fvz;
  double *fwx;
  double *fwy;
 
  fuy = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  fuz = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  fvx = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  fvz = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  fwx = (double*) malloc(dom.Gcc.s3 * sizeof(double));
  fwy = (double*) malloc(dom.Gcc.s3 * sizeof(double));    

  gradient_y(vf, dy, fuy);
  gradient_z(uf, dx, fuz);
  gradient_x(vf, dx, fvx);
  gradient_z(vf, dz, fvz);
  gradient_x(wf, dx, fwx);
  gradient_y(wf, dy, fwy);

  sub(dom.Gcc.s3, fwy, fvz, wx);
  sub(dom.Gcc.s3, fuz, fwx, wy);
  sub(dom.Gcc.s3, fvx, fuy, wz);

  free(fuy);
  free(fuz);
  free(fvx);
  free(fvz);
  free(fwx);
  free(fwy);
}
*/
double mean(int N, double *f)
{
  double mean_f = 0.0;
  for(int i = 0; i < N; i++){
    mean_f = mean_f + f[i];
  }
  mean_f = mean_f/N;
  return mean_f;
}


  
void add(int N, double *fx, double *fy, double *f_add)
{
  for(int i = 0; i < N; i++){
    f_add[i] = fx[i] + fy[i];
  }
}

void sub(int N, double *fx, double *fy, double *f_sub)
{
  for(int i = 0; i < N; i++){
    f_sub[i] = fx[i] - fy[i];
  }
}

void multiply(int N, double *fx, double *fy, double *f_multi)
{
  for(int i = 0; i < N; i++){
    f_multi[i] = fx[i]*fy[i];
  }
}

void divide(int N, double *fx, double *fy, double *f_divide)
{
  for(int i = 0; i < N; i++){
    f_divide[i] = fx[i]/fy[i];
  }
}

    
 
