
#include "rfweights.h"

SEXP R_rfweights (SEXP fdata, SEXP fnewdata, SEXP weights) {

    SEXP ans;
    double *dans, *dweights;
    int *id, *ind;
    int Ntree = LENGTH(fdata), Ndata, Nnewdata;
    int OOB = LENGTH(fnewdata) == 0;
    
    Ndata = LENGTH(VECTOR_ELT(fdata, 0));
    if (OOB) {
        Nnewdata = Ndata;
        fnewdata = fdata;
    } else {
        Nnewdata = LENGTH(VECTOR_ELT(fnewdata, 0));
    }
    
    PROTECT(ans = allocMatrix(REALSXP, Ndata, Nnewdata));
    dans = REAL(ans);
    
    for (int i = 0; i < Ndata * Nnewdata; i++)
        dans[i] = 0.0;
        
    for (int b = 0; b < Ntree; b++) {
        dweights = REAL(VECTOR_ELT(weights, b));
        ind = INTEGER(VECTOR_ELT(fnewdata, b));
        id = INTEGER(VECTOR_ELT(fdata, b));
        for (int j = 0; j < Nnewdata; j++) {
            if (OOB & (dweights[j] > 0)) continue;
            for (int i = 0; i < Ndata; i++) {
                if (id[i] == ind[j]) 
                    dans[j * Ndata + i] += dweights[i];
            }
        }
    }
    
    UNPROTECT(1);
    return(ans);
}
