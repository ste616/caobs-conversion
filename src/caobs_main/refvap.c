#include <complex.h>
#include <math.h>

#include "refvap.h"

// REAL in Fortran is 4-bytes float 
// REAL (KIND=8) in Fortran is 8-bytes float
// Assume that C/C++ double is 8 bytes and float is 4 bytes

float complex refvap(const float nu,
		     const float t,
		     const float pdry,
		     const float pvap){

  // Table of the microwave water lines.
  double nu0[NL] = {22.235080,  67.813960,  119.995940, 183.310117, 321.225644, 325.152919,
		    336.187000, 380.197372, 390.134508, 437.346667, 439.150812, 443.018295,
		    448.001075, 470.888947, 474.689127,	488.491133, 503.568532, 504.482692,
		    556.936002, 620.700807, 658.006500, 752.033227, 841.073593, 859.865000,
		    899.407000, 902.555000, 906.205524, 916.171582, 970.315022, 987.926764};
  
  float b1[NL] = {0.1090,  0.0011,  0.0007,  2.3000, 0.0464,  1.5400,  0.0010,  11.9000, 
		  0.0044,  0.0637,  0.9210,  0.1940, 10.6000, 0.3300,  1.2800,  0.2530,  
		  0.0374,  0.0125,  510.000, 5.0900, 0.2740,  250.000, 0.0130,  0.1330,  
		  0.0550,  0.0380,  0.1830,  8.5600, 9.1600,  138.000};
  
  float b2[NL] = {2.143, 8.730, 8.347, 0.653, 6.156, 1.515, 9.802, 1.018, 
		  7.318, 5.015, 3.561, 5.015, 1.370, 3.561, 2.342, 2.814, 
		  6.693, 6.693, 0.114, 2.150, 7.767, 0.336, 8.113, 7.989, 
		  7.845, 8.360, 5.039, 1.369, 1.842, 0.178};

  float b3[NL] = {27.84E-3, 27.60E-3, 27.00E-3, 28.35E-3, 21.40E-3, 27.00E-3, 26.50E-3, 27.60E-3,
		  19.00E-3, 13.70E-3, 16.40E-3, 14.40E-3, 23.80E-3, 18.20E-3, 19.80E-3, 24.90E-3,
		  11.50E-3, 11.90E-3, 30.00E-3, 22.30E-3, 30.00E-3, 28.60E-3, 14.10E-3, 28.60E-3,
		  28.60E-3, 26.40E-3, 23.40E-3, 25.30E-3, 24.00E-3, 28.60E-3};


  //Convert to the units of Liebe.
  float theta = 300.0 / t;
  float e     = 1.0E-3 * pvap;
  float p     = 1.0E-3 * pdry;
  float f     = 1.0E-9 * nu;

  double nr = 2.39 * e * theta + 41.6 * e * theta * theta + 6.47E-6*pow(f,2.05)*e*pow(theta,2.4);
  double ni = (0.915*1.40E-6*p + 5.41E-5*e*theta*theta*theta)* f*e*pow(theta,2.5);


  // Sum the contributions of the lines.
  for (int i = 0; i < NL; i++){
    double s     = b1[i] * e  * pow(theta,3.5) * exp (b2[i] * (1 - theta));
    double gamma = b3[i] * (p * pow(theta,0.8) + 4.80 * e * theta);
    double x     = (nu0[i] - f) * (nu0[i] - f) + gamma * gamma;
    double y     = (nu0[i] + f) * (nu0[i] + f) + gamma * gamma;
    double z     = (nu0[i] + gamma * gamma / nu0[i]);
    nr    = nr + s * ((z - f) / x + (z + f) / y - 2 / nu0[i]);
    ni    = ni + s * ((1/x+1/y)*gamma*f/nu0[i]);
  }


  // Return the result
  complex float result = CMPLX(nr, ni);
  
  return result;
}
