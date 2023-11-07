#include <R.h>
#include <Rmath.h>

void F77_SUB(fpnorm)(double* p, double* x, double* mu, double* sigma)
{
  *p = pnorm5(*x, *mu, *sigma, 1, 0);
}

void F77_SUB(fpt)(double* p, double* x, double* df)
{
  *p = pt(*x, *df, 1, 0);
}

void F77_SUB(fpf)(double* p, double* x, double* df1, double* df2)
{
  *p = pf(*x, *df1, *df2, 1, 0);
}

void F77_SUB(sortf)(double* x, int* n)
{
  R_rsort(x, *n);
}


