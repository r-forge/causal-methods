#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include "tlse.h"

static const R_FortranMethodDef fortranMethods[] = {
  {"selic", (DL_FUNC) &F77_SUB(selic), 25},
  {NULL, NULL, 0}
};

void R_init_causalTLSE(DllInfo *dll)
   {
     R_registerRoutines(dll,
			NULL, NULL, 
			fortranMethods, NULL);
     R_useDynamicSymbols(dll, FALSE);
   }

