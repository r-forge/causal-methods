#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include "slse.h"

static const R_FortranMethodDef fortranMethods[] = {
  {"selic", (DL_FUNC) &F77_SUB(selic), 26},
  {NULL, NULL, 0}
};

void R_init_causalSLSE(DllInfo *dll)
   {
     R_registerRoutines(dll,
			NULL, NULL, 
			fortranMethods, NULL);
     R_useDynamicSymbols(dll, FALSE);
   }

