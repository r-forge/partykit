#include <Rdefines.h>
#include <omp.h>
#include <R.h>
#include <Rmath.h>



double _Ginv(double rho, double alpha, double beta) {

    double res;
    double sigma = rho*rho;
    double s1    = 1. - sigma;
    res = rho * (alpha + sqrt(alpha*alpha * sigma + beta*beta * s1)) / s1;
    return res;

}

/* Ginv R/C interface function, implementation of the
   movMF:::Ginv function. Uses _Ginv internal function. 
   Copyright Statement: Modified version of function distributed under GPL-2
   in R-package 'movMF'(version 0.2-2) by Kurt Hornik and Bettina Gruen*/
SEXP Ginv(SEXP rho, SEXP alpha, SEXP beta) {

    int i, n = length(rho);
    double *rhoptr   = REAL(rho);
    double *alphaptr = REAL(alpha);
    double *betaptr  = REAL(beta);

    /* Return vector */
    SEXP rval = PROTECT(allocVector(REALSXP,n));
    double *rvalptr  = REAL(rval);

    for ( i = 0; i < n; i++ ) {
        rvalptr[i] = _Ginv(rhoptr[i], alphaptr[i], betaptr[i]);
    }

    UNPROTECT(1);
    return rval;

}



/* _Rinv_lower_Amos_bound internal function, implementation of the equation
   from movMF:::Rinv_lower_Amos_bound. Uses internal _Ginv function.
   See Rinv_lower_Amos_bound function for R/C interaction.
   Copyright Statement: Modified version of function distributed under GPL-2
   in R-package 'movMF'(version 0.2-2) by Kurt Hornik and Bettina Gruen*/
double _Rinv_lower_Amos_bound(double rho, double nu) {

    double tmp1, tmp2, res;
    tmp1 = _Ginv(rho, nu, nu + 2.);
    if ( nu == 0. ) {
        tmp2 = _Ginv(rho, 0.5, 0.8660254);
    } else {
        tmp2 = _Ginv(rho, nu + 0.5, sqrt((nu + 0.5) * (nu + 1.5)));
    }
    if ( tmp1 > tmp2 ) { res = tmp1; } else { res = tmp2; }
    return res;

}

/* Rinv_lower_Amos_bound R/C interface function, implementation of the
   movMF:::Rinv_lower_Amos_bound function.
   Uses _Rinv_lower_Amos_bound internal function.
   Copyright Statement: Modified version of function distributed under GPL-2
   in R-package 'movMF'(version 0.2-2) by Kurt Hornik and Bettina Gruen*/
SEXP Rinv_lower_Amos_bound(SEXP rho, SEXP nu) {

    int i, n = length(rho);
    double *rhoptr = REAL(rho);
    double *nuptr  = REAL(nu);

    /* Return vector */
    SEXP rval = PROTECT(allocVector(REALSXP,n));
    double *rvalptr  = REAL(rval);

    for ( i = 0; i < n; i++ ) {
        rvalptr[i] = _Rinv_lower_Amos_bound(rhoptr[i], nuptr[i]);
    }

    UNPROTECT(1);
    return rval;

}

/* _Rinv_upper_Amos_bound is the internal function. Would not really
   be required by Rinv_upper_Amos_bound, but I am using it in other
   functions. Thus, it is nicer to have a separate function here.
   For the R/C implementation check Rinv_upper_Amos_bound. */
double _Rinv_upper_Amos_bound(double rho, double nu) {

    double res;
    if ( nu == 0 ) {
        res = _Ginv(rho, 0.5, 1.5);
    } else {
        res = _Ginv(rho, nu + 0.5, nu + 1.5);
    }
    return res;

}
/* Rinv_upper_Amos_bound R/C interface function, implementation of the
   movMF:::Rinv_upper_Amos_bound function. Uses internal _Ginv function.
   Copyright Statement: Modified version of function distributed under GPL-2
   in R-package 'movMF'(version 0.2-2) by Kurt Hornik and Bettina Gruen*/
SEXP Rinv_upper_Amos_bound(SEXP rho, SEXP nu) {

    int i, n = length(rho);
    double *rhoptr = REAL(rho);
    double *nuptr  = REAL(nu);

    /* Return vector */
    SEXP rval = PROTECT(allocVector(REALSXP,n));
    double *rvalptr  = REAL(rval);

    for ( i = 0; i < n; i++ ) {
        rvalptr[i] = _Rinv_upper_Amos_bound(rhoptr[i], nuptr[i]);
    }

    UNPROTECT(1);
    return rval;

}


/* Literally a copy of the movMF C function C_mycfP
   Compute
     I_{\nu}(z) / I_{\nu - 1}(z)
   via the Gauss continued fraction, implemented using Eqn 3.2'
   in Gautschi and Slavik (1978), modified to use a_0 = z/(2\nu)
   instead of a_0 = 1.
   See also "Forward Series Recurrence Algorithm" in
   http://dlmf.nist.gov/3.10#iii.
   Copyright Statement: Modified version of function distributed under GPL-2
   in R-package 'movMF'(version 0.2-2) by Kurt Hornik and Bettina Gruen*/
double _mycfP(double z, double nu, double y, double tol) {

    double xG, rho, s, t, u, v;
    double res, p;
    int i, k;

    xG = z * z / 4.;
    p = res = z / (2. * nu);
    rho = 0.;
    u = nu * (nu - 1.);
    v = 2. * nu;
    k = 1;
    while ( fabs(p) > tol * s ) {
        u += v;
        v += 2.;
        t = xG * (1. + rho);
        rho = - t / (u + t);
        p   *= rho;
        res += p;
        //REprintf("k: %d p: %g s: %g\n", k, p, s);
        k++;
    }
    return res;

}




/* Internal function of the movMF:::A function.
   Solely for method "PCF" (see R manual) and with fixed tolerance
   to 1e-6 corresponding to the R default.
   Uses internal _A_PCF function to compute the results.
   Copyright Statement: Modified version of function distributed under GPL-2
   in R-package 'movMF'(version 0.2-2) by Kurt Hornik and Bettina Gruen*/
double _A_PCF(double kappa, double d, double tol) {

    double res;
    /* The one "d" most often used in the Newton kappa aproximation scheme */
    if ( d == 2. ) {
        if ( kappa > tol ) {
            res = _mycfP(kappa, 1., 1., tol);
        } else {
            double kappa3 = kappa*kappa*kappa;
            double kappa5 = kappa3*kappa*kappa;
            res = kappa * 0.5 - kappa3 * 0.0625 + 2. * kappa5 / 192;
        }
    } else {
        if ( kappa > tol ) {
            res = _mycfP(kappa, d * 0.5, 1., tol);
        } else {
            double kappa3 = kappa*kappa*kappa;
            double kappa5 = kappa3*kappa*kappa;
            double d2 = d*d;
            res = kappa / d - kappa3 / (d2 * (d + 2.)) + 
                2 * kappa5 / (d2*d * (d + 2.) * (d + 4.));
        }
    }

    return res;

}


/* R/C implementation of the movMF:::A function.
   Solely for method "PCF" (see R manual) and with fixed tolerance
   to 1e-06 corresponding to the R default.
   Uses internal _A_PCF function to compute the results.
   Copyright Statement: Modified version of function distributed under GPL-2
   in R-package 'movMF'(version 0.2-2) by Kurt Hornik and Bettina Gruen*/
SEXP A_PCF(SEXP kappa, SEXP d) {
    
    int i, n = length(kappa), nd = length(d);
    double tol = 1e-06;

    double *kappaptr = REAL(kappa);
    double *dptr     = REAL(d);

    SEXP A = PROTECT(allocVector(REALSXP,n));
    double *Aptr = REAL(A);

    /* "d" is allowd to have length 1 or the same length
     * as the input vector kappa. Exceptions are not handled! */
    if ( nd == 1 ) {
        for ( i = 0; i < n; i++ ) {
            Aptr[i] = _A_PCF(kappaptr[i], dptr[0], tol);
        }
    } else {
        for ( i = 0; i < n; i++ ) {
            Aptr[i] = _A_PCF(kappaptr[i], dptr[i], tol);
        }
    }

    UNPROTECT(1);
    return A;

}


/* Implementation of movMF:::Aprime. This is the internal function,
   for R/C interaction see Aprime_PCF. Implementation currently only
   for method = "PCF" (see R manual). 
   Copyright Statement: Modified version of function distributed under GPL-2
   in R-package 'movMF'(version 0.2-2) by Kurt Hornik and Bettina Gruen*/
double _Aprime_PCF(double kappa, double d, double A, double tol) {

    double res, tmp;
    if ( kappa > tol ) {
        if ( isnan(A) ) {
            tmp = _A_PCF(kappa, d, tol);
        } else {
            tmp = A;
        }
        /* Most common value in Newton iteration */
        if ( d == 2. ) {
            res = 1. - tmp*tmp - tmp / kappa;
        } else {
            res = 1. - tmp*tmp - tmp * (d - 1.) / kappa;
        }
    } else {
        /* Most common value in Newton iteration */
        double kappa2 = kappa * kappa;
        if ( d == 2. ) {
            res = 0.5 - 0.1875 * kappa2 + 0.05208333 * kappa2 * kappa2;
        } else {
            res = 1. / d - 3. / (d*d * (d + 2.)) * kappa2 + 
                  10. / (d*d*d * (d + 2.) * (d + 4.)) * kappa2 * kappa2;
        }
    }

    return res;

}

/* R/C implementation of the movMF:::Aprime function.
   Solely for method "PCF" (see R manual) and with fixed tolerance
   to 1e-06 corresponding to the R default.
   Uses internal _Aprime_PCF function to compute the results.
   Copyright Statement: Modified version of function distributed under GPL-2
   in R-package 'movMF'(version 0.2-2) by Kurt Hornik and Bettina Gruen*/
SEXP Aprime_PCF(SEXP kappa, SEXP d, SEXP A) {
    
    int i, n = length(kappa), nd = length(d), nA = length(A);
    double tol = 1e-06;
    double tmp;

    double di, Ai;

    double *kappaptr = REAL(kappa);
    double *dptr     = REAL(d);
    double *Aptr     = REAL(A);

    SEXP Aprime = PROTECT(allocVector(REALSXP,n));
    double *Aprimeptr = REAL(Aprime);

    /* "d" is allowd to have length 1 or the same length
     * as the input vector kappa. Exceptions are not handled! */
    for ( i = 0; i < n; i++ ) {
        if ( nd == 1 ) { di = dptr[0]; } else { di = dptr[i]; }
        if ( nA == 1 ) { Ai = Aptr[0]; } else { Ai = Aptr[i]; }

        Aprimeptr[i] = _Aprime_PCF(kappaptr[i], di, Ai, tol);
    }

    UNPROTECT(1);
    return Aprime;

}


/* C implementation of the Newton-Fourier kappa estimation.
   Tolerance is fixed to 1e-06, the R default. */

/*
       lower <- movMF:::Rinv_lower_Amos_bound(r, 0)
       upper <- movMF:::Rinv_upper_Amos_bound(r, 0)
       iter <- 1L
       while (iter <= maxiter) {
           A <- movMF:::A(lower, 2)
           Aprime <- movMF:::Aprime(lower, 2, A = A)
           lower <- lower - (A - r)/Aprime
           A <- movMF:::A(upper, 2)
           upper <- upper - (A - r)/Aprime
           if ((upper - lower) < tol * (lower + upper)) {
               if ((upper - lower) < -tol * (lower + upper))
                   stop("no convergence")
               break
           }
           iter <- iter + 1L
       }
       return((lower + upper)/2)
*/
    
   /*Copyright Statement: Modified version of function distributed under GPL-2
   in R-package 'movMF'(version 0.2-2) by Kurt Hornik and Bettina Gruen*/
SEXP solve_kappa_Newton_Fourier(SEXP r, SEXP maxiter) {

    int iter, i, n = length(r);
    int *maxiterptr = INTEGER(maxiter);

    double lower, upper, A, Aprime;
    double tol = 1e-6;
    double *rptr = REAL(r);

    /* Result vector */
    SEXP res = PROTECT(allocVector(REALSXP,n));
    double *resptr = REAL(res);

    i=0;
    for ( i = 0; i < n; i++ ) {
        lower = _Rinv_lower_Amos_bound(rptr[i], 0.);
        upper = _Rinv_upper_Amos_bound(rptr[i], 0.);

        /* Start iteration */
        iter = 1;
        while ( iter <= maxiterptr[0] ) {

            A      = _A_PCF(lower, 2., tol);
            Aprime = _Aprime_PCF(lower, 2., A, tol);
            lower  = lower - (A - rptr[i]) / Aprime;
            A      = _A_PCF(upper, 2., tol);
            upper  = upper - (A - rptr[i]) / Aprime;

            if ( (upper-lower) < (tol * (lower+upper)) ) {
                if ( (upper-lower) < (-1. * tol * (lower+upper)) ) {
                    Rprintf("DAMN IT, no convergence!\n");
                }
                break;
            }
            /* Increase increment */
            iter += 1;
        }

        /* Store the result */
        resptr[i] = (lower+upper) / 2.;
    }

    UNPROTECT(1);
    return res;

}



/* Von Mises density function.
   @param eta is a vector with two elements, the first one mu on the link
          scale, the second one kappa on the link scale.
   @param weights is a vector of weights.
   @param dlog. Integer, if 0 densities will be evaluated, else log-densities.
   @param sum. Integer, if 0 no sum will be computed, the densities will
          be returned element-wise. Else the sum is returned. This makes
          only sense with dlog = 1 (which then results in the log-likelihood).
   @return Returns a vector of length(y) with the densities or log-densities
          if sum = 0, else one single element with the sum. */
SEXP circ_density(SEXP y, SEXP eta, SEXP weights, SEXP dlog, SEXP sum, SEXP ncores) {

    int i, n;
    int *ncoresptr = INTEGER(ncores);
    int *dlogptr   = INTEGER(dlog);
    int *sumptr    = INTEGER(sum);

    double *yptr   = REAL(y);
    double *etaptr = REAL(eta);
    double *wptr   = REAL(weights);

    if ( sumptr[0] == 0 ) { n = length(y); } else { n = 1; }
    SEXP rval = PROTECT(allocVector(REALSXP,n));
    double *rvalptr = REAL(rval);

    /* Convert parameters from link-scale to parameter scale */
    double mu    = 2. * atan(etaptr[0]);
    double kappa = exp(etaptr[1]);
    double BE    = -1.837877 - log(bessel_i(kappa, 0, 2));
    //double BE    = -log(6.283185 * bessel_i(kappa, 0, 2));

    /* Multithreading via omp */
    //Rprintf("[C] Using %d cores now ...\n", ncoresptr[0]);
    omp_set_num_threads(ncoresptr[0]);

    if ( sumptr[0] == 0 ) {
        #pragma omp parallel
        {
            #pragma omp for
            for ( i = 0; i < length(y); i++ ) {
                rvalptr[i] = BE + kappa * (cos(yptr[i] - mu) - 1.);
                if ( dlogptr[0] == 0 ) { rvalptr[i] = exp(rvalptr[i]); }
            }
        }
    } else {
        double tmp, res = 0.;
        #pragma omp parallel private(tmp)
        {
            #pragma omp for reduction(+:res)
            for ( i = 0; i < length(y); i++ ) {
                tmp = BE + kappa * (cos(yptr[i] - mu) - 1.);
                if ( dlogptr[0] == 0 ) { res += exp(tmp); } else { res += tmp; }
            }
        }
        /* Store sum */
        rvalptr[0] = res;
    }

    /* Close open threads */
    omp_set_num_threads(1);

    /* Unprotect the one variable rval we have crated and return result */
    UNPROTECT(1);
    return rval;

}


/* Von Mises score function.
   @param eta is a vector with two elements, the first one mu on the link
          scale, the second one kappa on the link scale.
   @param weights is a vector of weights.
   @param sum. Integer, if 0 no sum will be computed per parameter.
   @return Returns a matrix of size length(y)*2 with element wise scores
          or a matrix of size 1*2 if sum = 1. */
SEXP circ_score(SEXP y, SEXP eta, SEXP weights, SEXP sum, SEXP ncores) {


    int i, n;
    int *sumptr    = INTEGER(sum);
    int *ncoresptr = INTEGER(ncores);

    double *yptr   = REAL(y);
    double *etaptr = REAL(eta);
    double *wptr   = REAL(weights);

    /* Setting up results matrix. Either a 1x2 matrix to take up
     * the sums for both parameters (mu and kappa) if sum = 1, else
     * a matrix of size n x 2 to take up the scores for each observation
     * for both parameters mu/kappa. */
    if ( sumptr[0] == 0 ) { n = length(y); } else { n = 1; }
    SEXP rval = PROTECT(allocMatrix(REALSXP,n,2));
    double *rvalptr = REAL(rval);

    /* Convert parameters from link-scale to parameter scale */
    double mu    = 2. * atan(etaptr[0]);
    double kappa = exp(etaptr[1]);
    double tmp;

    /* Initialize number of parallel threads if requested */
    //Rprintf("[C] Using %d cores now ...\n", ncoresptr[0]);
    omp_set_num_threads(ncoresptr[0]);

    /* Compute fixed part of score for mu once */
    double M1 = 2. * kappa / (pow(tan(mu * 0.5),2.) + 1.);
    /* Same for score for kappa: compute bessel once (constant) */
    double BE = bessel_i(kappa, 1, 2) / bessel_i(kappa, 0, 2);

    /* If sum == 0 (do not compute sum of scores) store element-wise
     * scores on rvalptr */
    if ( sumptr[0] == 0 ) {
        #pragma omp parallel
        {
            #pragma omp for
            for ( i = 0; i < length(y); i++ ) {
                rvalptr[i]   = sin(yptr[i] - mu) * M1;  // mu
                rvalptr[n+i] = kappa * ( cos(yptr[i] - mu) - BE ); // kappa
            }
        }
    } else {
        /* Using omp parallel to sum up element-wise scores
         * for both parameters mu/kappa. Store elements on smu
         * and skappa, omp takes care of sum up (reduction statement)
         * smu/skappa from all threads. */
        double smu = 0., skappa = 0.;
        #pragma omp parallel
        {
            #pragma omp for reduction(+:smu,skappa)
            for ( i = 0; i < length(y); i++ ) {
                smu    += sin(yptr[i] - mu) * M1; // mu
                skappa += kappa * ( cos(yptr[i] - mu) - BE ); // kappa
            }
        }
        /* Store score sums for mu and kappa */
        rvalptr[0] = smu;
        rvalptr[1] = skappa;
    }

    /* Close open threads */
    omp_set_num_threads(1);

    /* Unprotect the one variable rval we have crated and return result */
    UNPROTECT(1);
    return rval;

}














////* ============ MULTITHREADING GAUSSIAN DISTRIBUTION ======================= */
///
////* C-function to compute the density of a pargaussian (multithreading 
///   gaussian distribution).
///   @param y: observations where to evaluate the density
///   @param mu: location parameter of the gaussian distribution
///   @param sigma: scale parameter of the gaussian distribution
///   @return Returns real array of length(y) with the densities.
///   Created: 2017-02-24, Reto Stauffer.
///*/
///SEXP bamlss_pargaussian_density(SEXP y, SEXP mu, SEXP sigma, SEXP logarithm, SEXP nthreads )
///{
///   int i, n = length(y);
///
///   int *logptr = INTEGER(logarithm);
///   int *nthreadsptr = INTEGER(nthreads);
///   double *yptr     = REAL(y);
///   double *muptr    = REAL(mu);
///   double *sigmaptr = REAL(sigma);
///
///   SEXP rval = PROTECT(allocVector(REALSXP,n));
///   double *rvalptr = REAL(rval);
///
///   double sig2, ymu;
///
///   omp_set_num_threads( nthreadsptr[0] );
///   #pragma omp parallel
///   {
///      #pragma omp for firstprivate(sig2,ymu)
///      for(i = 0; i < n; i++) {
///         sig2 = sigmaptr[i] * sigmaptr[i];
///         ymu  = yptr[i] - muptr[i];
///         if ( logptr[0] == 0 ) {
///            rvalptr[i] = 1./sqrt(6.283185307179*sig2) * exp( -ymu*ymu/2./sig2 );
///         } else {
///            rvalptr[i] = -log(sqrt(6.283185307179*sig2)) - ymu*ymu/2./sig2;
///         }
///      }
///   }
///   omp_set_num_threads( 1 );
///   UNPROTECT(1);
///   return rval;
///}

