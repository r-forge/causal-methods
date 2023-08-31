#ifndef R_SLSE_H
#define R_SLSE_H

#include <R_ext/RS.h>


void F77_SUB(selic)(double *y0, double *y1, double *x0, double *x1,
		    int *n0, int *n1, int *p, double *tol, double *k0,
		    int *nk0, int *mnk0, int *tnk0, double *k1,
		    int *nk1, int *mnk1, int *tnk1, double *pval0,
		    double *pval1, double *spval, int *npval, double *bic, double *aic,
		    int *w0bic, int *w1bic, int *w0aic, int *w1aic);
		
#endif		     
