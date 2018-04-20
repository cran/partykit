
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "rfweights.h"

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}
#define REGCALL(name) R_RegisterCCallable("partykit", #name, (DL_FUNC) &name)

static const R_CallMethodDef callMethods[] = {
    CALLDEF(R_rfweights, 4),    
    {NULL, NULL, 0}
};

void attribute_visible R_init_partykit
(
    DllInfo *dll
) {

    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
    REGCALL(R_rfweights);
}
