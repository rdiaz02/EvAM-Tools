## Files MHN__*.R: MHN__UtilityFunctions.R, MHN__RegularizedOptimization.R,
## MHN__ModelConstruction.R, MHN__Likelihood.R, MHN__InlineFunctions.R,,
## MHN__ExampleApplications.R

##   Files obtained from https://github.com/RudiSchill/MHN
##   Commit 49a8cc0 from 2018-08-16
##   (we have added the "MHN__" and made minor modifications to conform to usage
##   within an R package).

##   License: MIT, which can be combined with software under the AGPL 3 as used
##   by the rest of this project. (Updated on 2022-04-04).

##   Author of code: from commit history, most likely Rudolf Schill.
  
##   Authors of paper/project: Schill, R., Solbrig, S., Wettig, T., & Spang, R.
  
##   Paper: Schill, R., Solbrig, S., Wettig, T., & Spang, R. (2020). Modelling
##   cancer progression using Mutual Hazard Networks. Bioinformatics, 36(1),
##   241–249. http://dx.doi.org/10.1093/bioinformatics/btz513


# TW: inline functions for operations on 2^n-dimensional vectors
# useful docs: [1] http://adv-r.had.co.nz/C-interface.html
#              [2] https://github.com/cjgeyer/mat

## require("inline")

#--------------------------------------------------------------------------------

## # multiplies a Kronecker product with a vector
## code_kronvec <- "
##   #include <R_ext/BLAS.h>
##   double *ptheta = REAL(Theta);
##   int i          = asInteger(i_);
##   double *px     = REAL(x);
##   int diag       = asLogical(diag_);
##   int transp     = asLogical(transp_);

##   int n       = length(Theta);
##   int nx      = length(x);
##   int nxhalf  = nx/2;
##   double mOne = -1;
##   double zero = 0;

##   SEXP out     = PROTECT(allocVector(REALSXP, nx)); 
##   double *pout = REAL(out);
##   double *ptmp = malloc(nx*sizeof(double));
##   double *px1,*px2;
##   double theta;
##   int incx     = 1;
##   int incx2    = 2;
##   int j;

##   // the argument of the function must not be modified
##   F77_CALL(dcopy)(&nx,px,&incx,ptmp,&incx);

##   for (j=1; j<=n; j++) {

##     // matrix shuffle

##     F77_CALL(dcopy)(&nxhalf,ptmp,&incx2,pout,&incx);
##     F77_CALL(dcopy)(&nxhalf,ptmp+1,&incx2,pout+nxhalf,&incx);

##     theta = ptheta[j-1];
##     px1   = pout;
##     px2   = pout + nxhalf;

##     // Kronecker product (daxpby is slightly slower than dcopy + dscal)

##     if (j == i) {
##       if (!transp) {
##         F77_CALL(dcopy)(&nxhalf,px1,&incx,px2,&incx);
##         F77_CALL(dscal)(&nxhalf,&theta,px2,&incx);
##         if (diag) {
##           F77_CALL(dcopy)(&nxhalf,px2,&incx,px1,&incx);
##           F77_CALL(dscal)(&nxhalf,&mOne,px1,&incx);
##         } else {
##           F77_CALL(dscal)(&nxhalf,&zero,px1,&incx);
##         }
##       } else {
##         F77_CALL(dcopy)(&nxhalf,px2,&incx,px1,&incx);
##         F77_CALL(dscal)(&nxhalf,&theta,px1,&incx);
##         if (diag) {
##           F77_CALL(dcopy)(&nxhalf,px1,&incx,px2,&incx);
##           F77_CALL(dscal)(&nxhalf,&mOne,px2,&incx);
##         } else {
##           F77_CALL(dscal)(&nxhalf,&zero,px2,&incx);
##         }
##       }
##     } else {
##       F77_CALL(dscal)(&nxhalf,&theta,px2,&incx);
##     }  

##     F77_CALL(dcopy)(&nx,pout,&incx,ptmp,&incx);
##   }

##   free(ptmp);
##   UNPROTECT(1);

##   return out;
## "
## kronvec <- cfunction(signature(Theta="numeric", i_="integer", x="numeric", diag_="logical", transp_="logical"),
##                      body=code_kronvec,
##                      language=c("C"))


kronvec <- function(Theta, i_, x, diag_, transp_)
    .Call(C_kronvec, Theta, i_, x, diag_, transp_,
          PACKAGE = "evamtools")


## #--------------------------------------------------------------------------------
## # loop over j in Grad
## code_grad_loop_j <- "
##   #include <R_ext/BLAS.h>
##   int i        = asInteger(i_);
##   int n        = asInteger(n_);
##   double *pr   = REAL(r);
##   int nx       = length(r);
##   int nxhalf   = nx/2;
##   SEXP G       = PROTECT(allocVector(REALSXP, n));
##   double *pG   = REAL(G);
##   double *pone = malloc(nxhalf*sizeof(double));
##   double *ptmp = malloc(nx*sizeof(double));
##   int incx     = 1;
##   int incx2    = 2;
##   int j;

##   #pragma omp parallel for shared(pone)
##   for (j=0; j<nxhalf; j++) pone[j] = 1;

##   for (j=0; j<n; j++) {
##     // matrix shuffle
##     F77_CALL(dcopy)(&nxhalf,pr,&incx2,ptmp,&incx);
##     F77_CALL(dcopy)(&nxhalf,pr+1,&incx2,ptmp+nxhalf,&incx);
##     F77_CALL(dcopy)(&nx,ptmp,&incx,pr,&incx);

##     // sums as dot products with 1
##     pG[j] = F77_CALL(ddot)(&nxhalf,pr+nxhalf,&incx,pone,&incx);
##     if (j == (i-1)) pG[j] += F77_CALL(ddot)(&nxhalf,pr,&incx,pone,&incx);
##   }

##   UNPROTECT(1);
##   free(pone);
##   free(ptmp);

##   return G;
## "
## grad_loop_j <- cfunction(signature(i_="integer", n_="integer", r="numeric"),
##                          body=code_grad_loop_j,
##                          language=c("C"),
##                          cppargs=c("-fopenmp -ftree-vectorize"))


grad_loop_j <- function(i_, n_, r) .Call(C_grad_loop_j, i_, n_, r,
                                         PACKAGE = "evamtools")


## RDU: for lack of a better place

.onUnload <- function (libpath) {
  library.dynam.unload("evamtools", libpath)
}
