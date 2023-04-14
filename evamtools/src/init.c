#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP C_grad_loop_j(SEXP, SEXP, SEXP);
extern SEXP C_kronvec(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"C_grad_loop_j", (DL_FUNC) &C_grad_loop_j, 3},
    {"C_kronvec",     (DL_FUNC) &C_kronvec,     5},
    {NULL, NULL, 0}
};

void R_init_evamtools(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
