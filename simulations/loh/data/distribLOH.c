#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_matrix.h>
#include "distribLOH.h"

#define BOUND 1e10
static const int NR = 40; /*no. of rows*/
static gsl_matrix* LOH;
static const double eps = 1e-10;

/*
load data
*/
void distrib_LOH_init() {
  /*load input data*/
  LOH = gsl_matrix_alloc(NR, 2);
  FILE* fdata = fopen("BarrettsLOH.dat", "r");
  if(!fdata)
    exit(1);
  gsl_matrix_fscanf(fdata, LOH);
  fclose(fdata);
}

static double log_invLogitDeriv(double y) {
  double ey = exp(y);
  return log( (ey * (1.0 + ey) - ey * ey) / ((1.0+ey) * (1.0+ey)) );
}

static double log_jacobian(double eta_1, double pi1_1, double pi2_1, double omega) {
  return( log_invLogitDeriv(eta_1) + log_invLogitDeriv(pi1_1) + log_invLogitDeriv(pi2_1));
}

static double ldbb(int x, int size, double prob, double omega) {
  if(omega < eps)
    return( log(gsl_ran_binomial_pdf(x, prob, size)) );

  double theta = 1.0 / omega;

  return( gsl_sf_lnchoose(size, x)
	  - gsl_sf_lnbeta(theta*(1-prob), theta*prob)
	  + gsl_sf_lnbeta((double) size - (double) x + theta*(1-prob), (double) x + theta*prob) );
}

static double my_lf(int x, int n, double eta, double pi1, double pi2, double omega) {
  double ans = eta * gsl_ran_binomial_pdf(x, pi1, n) +
    (1.0 - eta) * exp(ldbb(x, n, pi2, omega));
  return log(ans);
}

static double invlogit(double a) {
  return exp(a) / (1.0 + exp(a));
}

static double logit(double a) {
  return log(a / (1.0 - a));
}

/*parameters' unnormalized posterior distribution function
x: parameters vector*/
double distrib_LOH(void* ignore, gsl_vector* x) {
  double* theta = x->data;
  double eta = invlogit(theta[0]);
  double pi1 = invlogit(theta[1]);
  double pi2 = invlogit(theta[2]);
  double omega = exp(theta[3]) / (2.0 * ( 1.0 + exp(theta[3]) ) );
  if(!(isfinite(eta) && isfinite(pi1) && isfinite(pi2) && isfinite(omega)))
    return log(0.0); /*out of domain*/
  for(int i=0; i<3; i++) {
    if(fabs(theta[i]) > BOUND)
      return log(0.0); /*out of domain*/
  }
  if((eta < eps)  || (pi1 < eps) || (pi2 < eps) || (theta[3] < -30.0) || (theta[3] > 30.0) ||
     (eta > 1.0-eps) || (pi1 > 1.0-eps) || (pi2 > 1.0-eps)) {
    return log(0.0); /*out of domain*/
  }
  double ans  = 0.0;
  for(int i=0; i<NR && isfinite(ans); i++)
    ans += my_lf(gsl_matrix_get(LOH, i, 0), gsl_matrix_get(LOH, i, 1),
		 eta, pi1, pi2, omega);
  ans += log_jacobian(theta[0], theta[1], theta[2], theta[3]);
  return ans;
}

/* release memory */
void distrib_LOH_free(void) {
  gsl_matrix_free(LOH);
}
