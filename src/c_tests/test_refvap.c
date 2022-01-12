#include <stdlib.h>
#include <complex.h>
#include <stdio.h>

#include "../caobs_main/refvap.h"

int main(int argc, char *argv[]){
  
  float nu = 22000000000.0; // 22 GHz
  float t = 293.0; // Roughly 20 Celcius
  float pvap = 2317.0;
  float pdry = 101325.0 - pvap;

  float complex nvap = refvap(nu, t, pdry, pvap);
  fprintf(stdout, "%.8E\n", nu);
  fprintf(stdout, "%f\n", t);
  fprintf(stdout, "%f\n", pdry);
  fprintf(stdout, "%f\n", pvap);
  fprintf(stdout, "(%f, %.8E)\n", creal(nvap), cimag(nvap));
  
  return EXIT_SUCCESS;
}
