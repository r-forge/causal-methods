#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern void F77_NAME(pvalb)(void *, void *, void *, void *, void *,
			    void *, void *, void *, void *, void *, void *);

extern void F77_NAME(pvalf)(void *, void *, void *, void *, void *,
			    void *, void *, void *, void *, void *, void *);

extern void F77_NAME(selmodel)(void *, void *, void *, void *, void *,
			       void *, void *, void *, void *, void *,
			       void *, void *, void *, void *, void *,
			       void *, void *, void *, void *, void *,
			       void *, void *, void *, void *, void *,
			       void *, void *, void *, void *, void *,
			       void *, void *);

static const R_FortranMethodDef fortranMethods[] = {
  {"pvalb", (DL_FUNC) &F77_NAME(pvalb), 11},
  {"pvalf", (DL_FUNC) &F77_NAME(pvalf), 11},
  {"selmodel", (DL_FUNC) &F77_NAME(selmodel), 32},  
  {NULL, NULL, 0}
};

void R_init_causalSLSE(DllInfo *dll)
   {
     R_registerRoutines(dll,
			NULL, NULL, 
			fortranMethods, NULL);
     R_useDynamicSymbols(dll, FALSE);
   };





