#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern void F77_NAME(selic)(void *, void *, void *, void *, void *,
			    void *, void *, void *, void *, void *,
			    void *, void *, void *, void *, void *,
			    void *, void *, void *, void *, void *,
			    void *, void *, void *, void *, void *, void *);

extern void F77_NAME(myls)(void *, void *, void *, void *, void *,
			   void *, void *, void *, void *, void *, void *);

extern void F77_NAME(mypnorm)(void *, void *, void *, void *, void *);


static const R_FortranMethodDef fortranMethods[] = {
  {"selic", (DL_FUNC) &F77_NAME(selic), 26},
  {"myls", (DL_FUNC) &F77_NAME(myls), 11},
  {"mypnorm", (DL_FUNC) &F77_NAME(mypnorm), 5},
  {NULL, NULL, 0}
};

void R_init_causalSLSE(DllInfo *dll)
   {
     R_registerRoutines(dll,
			NULL, NULL, 
			fortranMethods, NULL);
     R_useDynamicSymbols(dll, FALSE);
   };




