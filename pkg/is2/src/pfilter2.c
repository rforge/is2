// -*- C++ -*-
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


static R_INLINE SEXP makearray (int rank, int *dim) {
  int nprotect = 0;
  int *dimp, k;
  double *xp;
  SEXP dimx, x;
  PROTECT(dimx = NEW_INTEGER(rank)); nprotect++;
  dimp = INTEGER(dimx); 
  for (k = 0; k < rank; k++) dimp[k] = dim[k];
  PROTECT(x = allocArray(REALSXP,dimx)); nprotect++;
  xp = REAL(x);
  for (k = 0; k < length(x); k++) xp[k] = NA_REAL;
  UNPROTECT(nprotect);
  return x;
}
static R_INLINE SEXP matchnames (SEXP x, SEXP names, const char *where) {
  int nprotect = 0;
  int n = length(names);
  int *idx, k;
  SEXP index, nm;
  PROTECT(nm = AS_CHARACTER(names)); nprotect++;
  PROTECT(index = match(x,names,0)); nprotect++;
  idx = INTEGER(index);
  for (k = 0; k < n; k++) {
    if (idx[k]==0) 
      error("variable '%s' not found among the %s",
	    CHARACTER_DATA(STRING_ELT(nm,k)),
	    where);
    idx[k] -= 1;
  }
  UNPROTECT(nprotect);
  return index;
}
static R_INLINE void setrownames (SEXP x, SEXP names, int n) {
  int nprotect = 0;
  SEXP dimnms, nm;
  PROTECT(nm = AS_CHARACTER(names)); nprotect++;
  PROTECT(dimnms = allocVector(VECSXP,n)); nprotect++;
  SET_ELEMENT(dimnms,0,nm);	// set row names
  SET_DIMNAMES(x,dimnms);
  UNPROTECT(nprotect);
}



void renormalize(double *w, int n, int log) {
  int i;

  if (log) {
    double max=w[0];
    for (i=1; i<n; i++) {
      if (w[i]>max) max = w[i];
    }
    for (i=0; i<n; i++) {
      w[i] = exp(w[i]-max);
    }
  }

  double sum=0;
  for (i=0; i<n; i++) {
    sum += w[i];
  }
  for (i=0; i<n; i++) {
    w[i] /= sum;
  }
}


static void nosort_resamp (int nw, double *w, int np, int *p, int offset) 
{
  int i, j;
  double du, u;

  for (j = 1; j < nw; j++) w[j] += w[j-1];

  if (w[nw-1] <= 0.0)
    error("non-positive sum of weights");

  du = w[nw-1] / ((double) np);
  u = runif(-du,0);

  for (i = 0, j = 0; j < np; j++) {   
    u += du;
    while (u > w[i]) i++;
    p[j] = i;
      
    
  }
 
if (offset)      // add offset if needed
    for (j = 0; j < np; j++) p[j] += offset;

}

SEXP systematic_resampling (SEXP weights)
{
  int nprotect = 0;
  int n;
  SEXP perm;

  GetRNGstate();
  n = LENGTH(weights);
  PROTECT(perm = NEW_INTEGER(n)); nprotect++;
  nosort_resamp(n,REAL(weights),n,INTEGER(perm),1);
  UNPROTECT(nprotect);
  PutRNGstate();
  return(perm);
}



// examines weights for filtering failure
// computes log likelihood and effective sample size
// computes (if desired) prediction mean, prediction variance, filtering mean.
// it is assumed that ncol(x) == ncol(params).
// weights are used in filtering mean computation.
// if length(weights) == 1, an unweighted average is computed.
// returns all of the above in a named list
SEXP pfilter2_computations (SEXP x, SEXP params, SEXP Np,
			   SEXP rw, SEXP rw_sd,
			   SEXP predmean, SEXP predvar,
			   SEXP filtmean, SEXP onepar,
			   SEXP weights, SEXP tol)
{
  int nprotect = 0;
  SEXP pm = R_NilValue, pv = R_NilValue, fm = R_NilValue;
  SEXP rw_names, ess, fail, loglik;
  SEXP newstates = R_NilValue, newparams = R_NilValue;
  SEXP retval, retvalnames;
  double *xpm = 0, *xpv = 0, *xfm = 0, *xw = 0, *xx = 0, *xp = 0, *xpw=0;
  int *xpa=0;
  SEXP dimX, dimP, newdim, Xnames, Pnames, pindex;
  SEXP pw=R_NilValue,pa=R_NilValue, psample=R_NilValue;
  int *dim, *pidx, lv, np;
  int nvars, npars = 0, nrw = 0, nreps, offset, nlost;
  int do_rw, do_pm, do_pv, do_fm, do_par_resamp, all_fail = 0;
  double sum, sumsq, vsq, ws, w, toler;
  int j, k;

  PROTECT(dimX = GET_DIM(x)); nprotect++;
  dim = INTEGER(dimX);
  nvars = dim[0]; nreps = dim[1];
  xx = REAL(x);
  PROTECT(Xnames = GET_ROWNAMES(GET_DIMNAMES(x))); nprotect++;

  PROTECT(dimP = GET_DIM(params)); nprotect++;
  dim = INTEGER(dimP);
  npars = dim[0];
  if (nreps != dim[1])
    error("'states' and 'params' do not agree in second dimension");
  PROTECT(Pnames = GET_ROWNAMES(GET_DIMNAMES(params))); nprotect++;

  np = INTEGER(AS_INTEGER(Np))[0]; // number of particles to resample

  PROTECT(rw_names = GET_NAMES(rw_sd)); nprotect++; // names of parameters undergoing random walk

  do_rw = *(LOGICAL(AS_LOGICAL(rw))); // do random walk in parameters?
  do_pm = *(LOGICAL(AS_LOGICAL(predmean))); // calculate prediction means?
  do_pv = *(LOGICAL(AS_LOGICAL(predvar)));  // calculate prediction variances?
  do_fm = *(LOGICAL(AS_LOGICAL(filtmean))); // calculate filtering means?
  do_par_resamp = *(LOGICAL(AS_LOGICAL(onepar))); // are all cols of 'params' the same?
  do_par_resamp = !do_par_resamp || do_rw || (np != nreps); // should we do parameter resampling?

  PROTECT(ess = NEW_NUMERIC(1)); nprotect++; // effective sample size
  PROTECT(loglik = NEW_NUMERIC(1)); nprotect++; // log likelihood
  PROTECT(fail = NEW_LOGICAL(1)); nprotect++;	// particle failure?

  xw = REAL(weights); 
  toler = *(REAL(tol));		// failure tolerance
  
    
  // check the weights and compute sum and sum of squares
  for (k = 0, w = 0, ws = 0, nlost = 0; k < nreps; k++) {
    
    if (xw[k] >= 0) {	
     
      w += xw[k];
      ws += xw[k]*xw[k];
    } else {			// this particle is lost
      xw[k] = 0;
      nlost++;
    }
  }
  if (nlost >= nreps) all_fail = 1; // all particles are lost
  if (all_fail) {
    *(REAL(loglik)) = log(toler); // minimum log-likelihood
    *(REAL(ess)) = 0;		  // zero effective sample size
  } else {
    *(REAL(loglik)) = log(w/((double) nreps)); // mean of weights is likelihood
    *(REAL(ess)) = w*w/ws;	// effective sample size
  }
  *(LOGICAL(fail)) = all_fail;

  if (do_rw) {
    // indices of parameters undergoing random walk
    PROTECT(pindex = matchnames(Pnames,rw_names,"parameters")); nprotect++; 
    xp = REAL(params);
    pidx = INTEGER(pindex);
    nrw = LENGTH(rw_names);
    lv = nvars+nrw;
  } else {
    pidx = NULL;
    lv = nvars;
  }

  if (do_pm || do_pv) {
    PROTECT(pm = NEW_NUMERIC(lv)); nprotect++;
    xpm = REAL(pm);
  }

  if (do_pv) {
    PROTECT(pv = NEW_NUMERIC(lv)); nprotect++;
    xpv = REAL(pv);
  }

  if (do_fm) {
    if (do_rw) {
      PROTECT(fm = NEW_NUMERIC(nvars+npars)); nprotect++;
    } else {
      PROTECT(fm = NEW_NUMERIC(nvars)); nprotect++;
    }
    xfm = REAL(fm);
  }
  
  PROTECT(pa = NEW_INTEGER(np)); nprotect++;
  xpa = INTEGER(pa);
  
  
  
  for (j = 0; j < nvars; j++) {	// state variables

    // compute prediction mean
    if (do_pm || do_pv) {
      for (k = 0, sum = 0; k < nreps; k++) sum += xx[j+k*nvars];
      sum /= ((double) nreps);
      xpm[j] = sum;
    }

    // compute prediction variance
    if (do_pv) {	
      for (k = 0, sumsq = 0; k < nreps; k++) {
	vsq = xx[j+k*nvars]-sum;
	sumsq += vsq*vsq;
      }
      xpv[j] = sumsq / ((double) (nreps - 1));
      
    }

    //  compute filter mean
    if (do_fm) {
      if (all_fail) {		// unweighted average
	for (k = 0, ws = 0; k < nreps; k++) ws += xx[j+k*nvars]; 
	xfm[j] = ws/((double) nreps);
      } else { 			// weighted average
	for (k = 0, ws = 0; k < nreps; k++) ws += xx[j+k*nvars]*xw[k]; 
	xfm[j] = ws/w;
      }
    }

  }

  // compute means and variances for parameters (if needed)
  if (do_rw) {
    for (j = 0; j < nrw; j++) {
      offset = pidx[j];		// position of the parameter

      if (do_pm || do_pv) {
	for (k = 0, sum = 0; k < nreps; k++) sum += xp[offset+k*npars];
	sum /= ((double) nreps);
	xpm[nvars+j] = sum;
      }

      if (do_pv) {
	for (k = 0, sumsq = 0; k < nreps; k++) {
	  vsq = xp[offset+k*npars]-sum;
	  sumsq += vsq*vsq;
	}
	xpv[nvars+j] = sumsq / ((double) (nreps - 1));
      }

    }

    if (do_fm) {
      for (j = 0; j < npars; j++) {
	if (all_fail) {		// unweighted average
	  for (k = 0, ws = 0; k < nreps; k++) ws += xp[j+k*npars];
	  xfm[nvars+j] = ws/((double) nreps);
	} else {		// weighted average
	  for (k = 0, ws = 0; k < nreps; k++) ws += xp[j+k*npars]*xw[k];
	  xfm[nvars+j] = ws/w;
	}
      }
    }
  }

  GetRNGstate();

  if (!all_fail) { // resample the particles unless we have filtering failure
    int xdim[2];
    //int sample[np];
    double *ss = 0, *st = 0, *ps = 0, *pt = 0;

    // create storage for new states
    xdim[0] = nvars; xdim[1] = np;
    PROTECT(newstates = makearray(2,xdim)); nprotect++;
    setrownames(newstates,Xnames,2);
    ss = REAL(x);
    st = REAL(newstates);

    // create storage for new parameters
    if (do_par_resamp) {
      xdim[0] = npars; xdim[1] = np;
      PROTECT(newparams = makearray(2,xdim)); nprotect++;
      setrownames(newparams,Pnames,2);
      ps = REAL(params);
      pt = REAL(newparams);
    }
    
    PROTECT(pw = NEW_NUMERIC(nreps)); nprotect++;
    xpw = REAL(pw);
    for (k = 0; k < nreps; k++)
      xpw[k]=REAL(weights)[k];
    nosort_resamp(nreps,REAL(weights),np,xpa,0);
    for (k = 0; k < np; k++) { // copy the particles
      for (j = 0, xx = ss+nvars*xpa[k]; j < nvars; j++, st++, xx++) 
	*st = *xx;
      
          
	        
      if (do_par_resamp) {
	for (j = 0, xp = ps+npars*xpa[k]; j < npars; j++, pt++, xp++){
    *pt = *xp;
   
	} 
	  
      }
    }

  } else { // don't resample: just drop 3rd dimension in x prior to return
    
    PROTECT(newdim = NEW_INTEGER(2)); nprotect++;
    dim = INTEGER(newdim);
    dim[0] = nvars; dim[1] = nreps;
    SET_DIM(x,newdim);
    setrownames(x,Xnames,2);

  }
    
  if (do_rw) { // if random walk, adjust prediction variance and move particles
    xx = REAL(rw_sd);
    xp = (all_fail || !do_par_resamp) ? REAL(params) : REAL(newparams);
    nreps = (all_fail) ? nreps : np;

    for (j = 0; j < nrw; j++) {
      offset = pidx[j];
      vsq = xx[j];
      if (do_pv) {
	xpv[nvars+j] += vsq*vsq;
      }
      for (k = 0; k < nreps; k++)
	xp[offset+k*npars] += rnorm(0,vsq);
    }
  }
  
  renormalize(xpw,nreps,0);
  PutRNGstate();

  PROTECT(retval = NEW_LIST(10)); nprotect++;
  PROTECT(retvalnames = NEW_CHARACTER(10)); nprotect++;
  SET_STRING_ELT(retvalnames,0,mkChar("fail"));
  SET_STRING_ELT(retvalnames,1,mkChar("loglik"));
  SET_STRING_ELT(retvalnames,2,mkChar("ess"));
  SET_STRING_ELT(retvalnames,3,mkChar("states"));
  SET_STRING_ELT(retvalnames,4,mkChar("params"));
  SET_STRING_ELT(retvalnames,5,mkChar("pm"));
  SET_STRING_ELT(retvalnames,6,mkChar("pv"));
  SET_STRING_ELT(retvalnames,7,mkChar("fm"));
  SET_STRING_ELT(retvalnames,8,mkChar("weight"));
  SET_STRING_ELT(retvalnames,9,mkChar("pa"));
  
  SET_NAMES(retval,retvalnames);

  SET_ELEMENT(retval,0,fail);
  SET_ELEMENT(retval,1,loglik);
  SET_ELEMENT(retval,2,ess);
  
  if (all_fail) {
    SET_ELEMENT(retval,3,x);
  } else {
    SET_ELEMENT(retval,3,newstates);
  }

  if (all_fail || !do_par_resamp) {
    SET_ELEMENT(retval,4,params);
  } else {
    SET_ELEMENT(retval,4,newparams);
  }

  if (do_pm) {
    SET_ELEMENT(retval,5,pm);
  }
  if (do_pv) {
    SET_ELEMENT(retval,6,pv);
  }
  if (do_fm) {
    SET_ELEMENT(retval,7,fm);
  }
  SET_ELEMENT(retval,8,pw);
  SET_ELEMENT(retval,9,pa);
  UNPROTECT(nprotect);
  return(retval);
}
