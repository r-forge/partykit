#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP Ginv(SEXP, SEXP, SEXP);
SEXP Rinv_lower_Amos_bound(SEXP, SEXP);
SEXP Rinv_upper_Amos_bound(SEXP, SEXP);
SEXP A_GCF(SEXP, SEXP);
SEXP Aprime_GCF(SEXP, SEXP, SEXP);
SEXP solve_kappa_Newton_Fourier(SEXP, SEXP);
SEXP circ_density(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP circ_score(SEXP, SEXP, SEXP, SEXP, SEXP);

static R_CallMethodDef callMethods[] = {
  {"Ginv", (DL_FUNC) &Ginv, 3},
  {"Rinv_lower_Amos_bound", (DL_FUNC) &Rinv_lower_Amos_bound, 2},
  {"Rinv_upper_Amos_bound", (DL_FUNC) &Rinv_upper_Amos_bound, 2},
  {"A_GCF", (DL_FUNC) &A_GCF, 2},
  {"Aprime_GCF", (DL_FUNC) &Aprime_GCF, 3},
  {"solve_kappa_Newton_Fourier", (DL_FUNC) &solve_kappa_Newton_Fourier, 2},
  {"circ_density", (DL_FUNC) &circ_density, 6},
  {"circ_score", (DL_FUNC) &circ_score, 5},
  {NULL, NULL, 0}
};

void R_init_sourcetools(DllInfo* info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

