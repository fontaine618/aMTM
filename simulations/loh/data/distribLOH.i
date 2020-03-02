%module distribLOH

%{
#include "distribLOH.h"
%}

%init %{
  distrib_LOH_init();
%}

%callback("%s_cb");
double distrib_LOH(void* ignore, gsl_vector* x);
%nocallback;
