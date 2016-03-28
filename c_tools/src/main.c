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

  // Read input file; read tStart and tEnd
  main_read_input();
  
  // Read and sort output directory for finding files within our time limits
  init_part_files();
  init_flow_files();
  // Create output directory
  //create_output_dir();
  // Get number of particles
  nparts = cgns_read_nparts();
  // Initialize domain and flow arrays
  domain_init();
  // Initialize partstruct and flow vars
  parts_init();
  flow_init();

  //before time iteration, malloc memory
  malloc_dataproc();
  // Generate statistics
  for (tt = 0; tt < nFiles; tt++) {
    // Read in new data and push to device
    cgns_fill_flow();
    cgns_fill_parts();
 
    calculate_gradient();
    vorticity();
    analyze_3d("3d_data.dat");  
    //analyze_2d("2d_data");
    printf("data processing %d is finished!\n", tt);
  }
  // calculate the pdf scalar in one plane
  //analyze_pdf_2d(500, 0, nFiles);
  //analyze_pdf_3d(200, 500, nFiles);

  //free memory
  free_dataproc();  
  printf("finish time iteration!\n");
 
  // Free and exit
  free_flow_vars();
  //free_part_vars();
  return EXIT_SUCCESS;
}


