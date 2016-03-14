#include "main.h"

// Define global variables declared in header file
int dev_start;
double tStart;
double tEnd;
int tt;

int main(int argc, char *argv[]) 
{
  if(--argc > 0){
	tStart = atof(argv[1]);
	tEnd = atof(argv[2]);
  }
  printf("start time is %f\n", tStart);
  printf("end time is %f\n", tEnd);

  // Read input file; read tStart and tEnd
  //main_read_input();

  // Read and sort output directory for finding files within our time limits
  //init_part_files();
  init_flow_files();
  // Create output directory
  create_output_dir();
  // Get number of particles
  //cgns_read_nparts();

  // Initialize domain and flow arrays
  domain_init();

  // Initialize partstruct and flow vars
  //parts_init();
  flow_init();

  // Messy hack for taking advantage of CUDA_VISIBLE_DEVICES in SLURM
  //dev_start = read_devices();

  // Allocate device memory
  //cuda_dev_malloc();
  //cuda_dom_push(); 
  //cuda_part_push();
  //cuda_part_pull();

  analyze_3d("3d_data.dat");

  /*// Calculate tetrad stats
  for (tt = 0; tt < nFiles; tt++) {
    // Read in new data and push to device
    cgns_fill_flow();
    analyze_3d("3d_data.dat");  
  }*/

  printf("finish time iteration!\n");

  // write phaseAvereaged to file
  //write_part_data();
   
  // Free and exit
  free_flow_vars();
  //free_part_vars();
  //cuda_dev_free();
  return EXIT_SUCCESS;
}


