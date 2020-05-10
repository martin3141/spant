// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_spant(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

// #include <R_ext/RS.h>
// #include <stdlib.h> // for NULL
// #include <R_ext/Rdynload.h>
// 
// /* FIXME: 
//    Check these declarations against the C/Fortran source code.
// */
// 
// /* .C calls */
// extern void indx(void *, void *, void *, void *, void *);
// extern void matMaxs(void *, void *, void *, void *, void *);
// 
// /* .Fortran calls */
// extern void F77_NAME(hfti)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
// extern void F77_NAME(nnls)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
// extern void F77_NAME(pnnls)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
// extern void F77_NAME(svdrs)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
// 
// static const R_CMethodDef CEntries[] = {
//     {"indx",    (DL_FUNC) &indx,    5},
//     {"matMaxs", (DL_FUNC) &matMaxs, 5},
//     {NULL, NULL, 0}
// };
// 
// static const R_FortranMethodDef FortranEntries[] = {
//     {"hfti",  (DL_FUNC) &F77_NAME(hfti),  13},
//     {"nnls",  (DL_FUNC) &F77_NAME(nnls),  11},
//     {"pnnls", (DL_FUNC) &F77_NAME(pnnls), 12},
//     {"svdrs", (DL_FUNC) &F77_NAME(svdrs),  9},
//     {NULL, NULL, 0}
// };
// 
// void R_init_lsei(DllInfo *dll)
// {
//     R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
//     R_useDynamicSymbols(dll, FALSE);
// }


