#include <stdlib.h>
#include <complex.h>
#include <stdio.h>

#include "../caobs_main/refvap.h"

int main(int argc, char *argv[]){
  
  double nu = 22000000000.0; // 22 GHz
  //double nu = 2.20000010E+10;
  double t = 293.0; // Roughly 20 Celcius
  double pvap = 2317.0;
  double pdry = 101325.0 - pvap;

  double complex nvap = refvap(nu, t, pdry, pvap);
  fprintf(stdout, "%.8E\n", nu);
  fprintf(stdout, "%f\n", t);
  fprintf(stdout, "%f\n", pdry);
  fprintf(stdout, "%f\n", pvap);
  fprintf(stdout, "(%f, %.8E)\n", creal(nvap), cimag(nvap));
  
  return EXIT_SUCCESS;
}
