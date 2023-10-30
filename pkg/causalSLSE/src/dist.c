#include <R.h>
#include <Rmath.h>

void F77_SUB(fpnorm)(double* p, double* x, double* mu, double* sigma,
		    int* lower_tail, int* log_p)
{
  *p = pnorm5(*x, *mu, *sigma, *lower_tail, *log_p);
}

