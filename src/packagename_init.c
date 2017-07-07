#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 *    Check these declarations against the C/Fortran source code.
 *    */

/* .Call calls */
extern SEXP FUNLDA_ebmme_cpp_binned(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FUNLDA_newtissue(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FUNLDA_predictlogsum(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"FUNLDA_ebmme_cpp_binned", (DL_FUNC) &FUNLDA_ebmme_cpp_binned, 12},
    {"FUNLDA_newtissue",        (DL_FUNC) &FUNLDA_newtissue,         7},
    {"FUNLDA_predictlogsum",    (DL_FUNC) &FUNLDA_predictlogsum,     6},
    {NULL, NULL, 0}
};

void R_init_FUNLDA(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

