#ifndef __DISTRIB_LOH__
#define __DISTRIB_LOH__
#include <gsl/gsl_vector.h>

void distrib_LOH_init(void);
double distrib_LOH(void* ignore, gsl_vector* x);
void distrib_LOH_free(void);

#endif
