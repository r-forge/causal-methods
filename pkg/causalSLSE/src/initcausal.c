#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern void F77_NAME(pvalb)(void *, void *, void *, void *, void *,
			    void *, void *, void *, void *, void *, void *);

extern void F77_NAME(pvalf)(void *, void *, void *, void *, void *,
			    void *, void *, void *, void *, void *, void *);

extern void F77_NAME(selcmodel)(void *, void *, void *, void *, void *,
			       void *, void *, void *, void *, void *,
			       void *, void *, void *, void *, void *,
			       void *, void *, void *, void *, void *,
			       void *, void *, void *, void *, void *,
			       void *, void *, void *, void *, void *,
			       void *, void *);

extern void F77_NAME(selmodel)(void *, void *, void *, void *, void *,
			       void *, void *, void *, void *, void *,
			       void *, void *, void *, void *, void *,
			       void *, void *, void *, void *, void *);

extern void F77_NAME(cvfct)(void *, void *, void *, void *, void *,
			    void *, void *);

extern void F77_NAME(llr)(void *, void *, void *, void *, void *,
			  void *, void *, void *, void *);

extern void F77_NAME(gridcv)(void *, void *, void *, void *, void *,
			     void *, void *, void *, void *, void *);

extern void F77_NAME(brentcv)(void *, void *, void *, void *, void *,
			      void *, void *, void *, void *, void *,
			      void *);

extern void F77_NAME(findnn)(void *, void *, void *, void *, void *,
			     void *, void *, void *, void *, void *,
			     void *, void *);

static const R_FortranMethodDef fortranMethods[] = {
  {"pvalb", (DL_FUNC) &F77_NAME(pvalb), 11},
  {"pvalf", (DL_FUNC) &F77_NAME(pvalf), 11},
  {"selcmodel", (DL_FUNC) &F77_NAME(selcmodel), 32},
  {"selmodel", (DL_FUNC) &F77_NAME(selmodel), 20},
  {"cvfct", (DL_FUNC) &F77_NAME(cvfct), 7},
  {"llr", (DL_FUNC) &F77_NAME(llr), 9},
  {"gridcv", (DL_FUNC) &F77_NAME(gridcv), 10},
  {"brentcv", (DL_FUNC) &F77_NAME(brentcv), 11},
  {"findnn", (DL_FUNC) &F77_NAME(findnn), 12},          
  {NULL, NULL, 0}
};

void R_init_causalSLSE(DllInfo *dll)
   {
     R_registerRoutines(dll,
			NULL, NULL, 
			fortranMethods, NULL);
     R_useDynamicSymbols(dll, FALSE);
   }





