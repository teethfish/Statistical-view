#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define FILE_NAME_SIZE 20
#define INPUT_DIR '/Users/yayun/Documents/MATLAB/dataprocessing/seeder'

void main() 
{
  // This function is used to give a random field. Each particle can move freely and it will check the distance between another particle
  // bias is the ratio of how  much it can move 
  
  // seed the random number generator
  srand(time(NULL));

  int Nx = 5; 
  int Ny = 4;
  int Nz = 4;
  int NP = Nx*Ny*Nz; //Tototal number of particles
  double a = 0.5; //radius of particles
  double bias = 0.2;
  double xs = 0.0;
  double xe = 5;
  double ys = 0.0;
  double ye = 5;
  double zs = 0.0;
  double ze = 5;

  //int times = 110000;//total pertubation times;
  char filename[FILE_NAME_SIZE];  //define the name of output file
  int w_t = 100000; //start writing the results from w_t
  int n_w = 100; //write n_w data files
  int d_w = 100; //write results every d_w steps
  int times = w_t + (n_w - 1)*d_w;//total pertubation times;
    
  printf("Running bluebottle seeder for %d particles...\n\n", Nx*Ny*Nz);
  fflush(stdout);
  int fail = 0;

  double xl = fabs(xe - xs);
  double yl = fabs(ye - ys);
  double zl = fabs(ze - zs);
  double dx = xl/Nx; //dx in the distance between centers of two nearby particles in x direction
  double dy = yl/Ny;
  double dz = zl/Nz;

  double d_pair = 100.0 * a;
  double d_x = 0.0;
  double d_y = 0.0;
  double d_z = 0.0;
  double d_min = 100 * a;

  if(dx < 2*a)
  {
  printf(" Too many particles in x direction\n");
  fail = !fail;
  }
  if(dy < 2*a)
  {
  printf("Too many particles in y direction\n");
  fail = !fail; 
  }
  if(dz < 2*a)
  {
  printf("Too many particles in z direction\n");
  fail = !fail; 
  } 
   if(fail) {
    printf("...bluebottle seeder done.\n\n");
    exit(EXIT_FAILURE);
  }

  double *acceptance;  //acceptance rate 
  acceptance = (double*)calloc(times, sizeof(double));
  for (int t = 0; t < times; t++) acceptance[t] = 0; 
  double *x;
  x = (double*)calloc(NP, sizeof(double));
  double *y;
  y = (double*)calloc(NP, sizeof(double));
  double *z;
  z = (double*)calloc(NP, sizeof(double));
  double *x_new;
  x_new = (double*)calloc(NP, sizeof(double));
  double *y_new;
  y_new = (double*)calloc(NP, sizeof(double));
  double *z_new;
  z_new = (double*)calloc(NP, sizeof(double));
  int *move;
  move = (int*)calloc(NP, sizeof(int));
  // Set the initial regular domain
  int np = 0;
  for (int k = 0; k < Nz; k++)
  { 
  	for (int j = 0; j < Ny; j++)
  	{  
  		for (int i = 0; i < Nx; i++)
  		{
          x[np] = xs + (2*i + 1)*dx/2; 
          y[np] = ys + (2*j + 1)*dy/2;
          z[np] = zs + (2*k + 1)*dz/2; 
          np = np + 1;
		  }
	  }
  }
 
  for (int t = 0; t < times; t++)
  {
    for (int i = 0; i < NP; i++)//all particles move once
    {
      x_new[i] = x[i] + bias*(-1 + 2*((double) rand() /(double) (RAND_MAX)))*a*2.0;
      y_new[i] = y[i] + bias*(-1 + 2*((double) rand() /(double) (RAND_MAX)))*a*2.0;
      z_new[i] = z[i] + bias*(-1 + 2*((double) rand() /(double) (RAND_MAX)))*a*2.0;     
      move[i] = 0;
      if (x_new[i] > xe) x_new[i] = x_new[i] - xl;
      if (x_new[i] < xs) x_new[i] = x_new[i] + xl;
      if (y_new[i] > ye) y_new[i] = y_new[i] - yl;
      if (y_new[i] < ys) y_new[i] = y_new[i] + yl;
      if (z_new[i] > ze) z_new[i] = z_new[i] - zl;
      if (z_new[i] < zs) z_new[i] = z_new[i] + zl;      
    }  

    // calculate the distance between particles
    for (int j = 0; j < NP; j++)  
    {
      d_min = 100*a*2;
      for (int k = 0; k < NP; k++)
      {
        //printf("y= %f\n",y[k]);         
        d_x = x_new[j] - x[k];
        d_y = y_new[j] - y[k];
        d_z = z_new[j] - z[k];
        double tmp = fabs(d_y);
        if (fabs(d_x) > xl/2.0) d_x = xl - fabs(d_x);
        if (fabs(d_y) > yl/2.0) d_y = yl - fabs(d_y);
        if (fabs(d_z) > zl/2.0) d_z = zl - fabs(d_z);
        if (k == j)d_pair = 100*a*2;        
        else d_pair = sqrt(d_x*d_x + d_y*d_y + d_z*d_z);       
        if(d_pair <= d_min) d_min = d_pair;
      }
      if(d_min > 2*a)
      {
        move[j] = 1;
        acceptance[t] = acceptance[t] + 1.0/NP;
      }  
    }
    //update the new position
    for (int i = 0; i < NP; i++)
    {
      if (move[i] == 1)
      {
        x[i] = x_new[i];
        y[i] = y_new[i];
        z[i] = z_new[i];
      }
    }
    //write the results
    for(int p = 0; p <= n_w; p++)
    {
      if(t == w_t + p * d_w)
      {
      printf("Writing part_seeder.input...");
      printf("time = %d\n",t);
      printf("acceptance = %f\n",acceptance[t-1]); 
      fflush(stdout);

      sprintf(filename, "myrandom_%d.txt", p);
      //char filename[FILE_NAME_SIZE] = (filename,t);
      FILE *ofile = fopen(filename, "w");
      if(ofile == NULL) 
      {
       fprintf(stderr, "Could not open file %s\n", filename);
       exit(EXIT_FAILURE);
      }
      for(int i = 0; i < NP; i++) 
      {
       fprintf(ofile, "%f %f %f\n", x[i], y[i], z[i]);
      }
      // close the file
      fclose(ofile);
      }  
   }
  }  
  for( int t = 1; t < times; t++)
  {
    acceptance[0] = acceptance[0] + acceptance[t];
  }
  acceptance[0] = acceptance[0]/times;
  printf("acceptance = %f\n",acceptance[0]); 

  printf("done.\n");
  printf("\n...bluebottle seeder done.\n\n");
  fflush(stdout);

  free(x);
  free(y);
  free(z);			
}  			
  			
  			
  			
  			
  			
  			
  			
  			
  			
  			
  			
  			
 
