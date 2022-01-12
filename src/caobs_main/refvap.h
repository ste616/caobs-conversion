#ifndef REFVAP_H
#define REFVAP_H

#include <complex.h>

#define NL 30

float complex refvap(const float nu,
		     const float t,
		     const float pdry,
		     const float pvap);

#endif
