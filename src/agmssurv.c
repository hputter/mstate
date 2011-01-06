/*  agmssurv.c	*/
/*
** Modeled after agsurv2.c from the survival library by Terry Therneau
** Called from msfit
**
**  Input
**    n - no of subjects
**    p - no of vars in xmat
**    var - NEW: 1 if variances are to be calculated, 0 otherwise
**	  method - CHANGED: method for handling ties, 1=breslow, 2=efron
**    H - NEW: number of strata in coxph call
**    K - no of curves (= no of transitions)
**    nt - NEW: no of unique time points joint from all strata
**    y[n,3] - matrix containing (start, stop, event)
**    score[n] - vector of weights exp(beta^T Z_i)
**    xmat[n,p] - data matrix that generated the Cox fit
**    varcov[p,p] - inverse Fisher information for estimation of beta
**    strata[H+1] - CHANGED: strata[0] = 0, strata[1:H]= last obs strata 1,2, etc in y
**    kstrata[K] - NEW: kstrata[k] = stratum in newdata[k,]
**    unt[nt] - NEW: vector of unique time points joint from all strata
**    newx[K,p] - matrix with covariates Z for each individual
**    newrisk[K] - vector containing exp(\beta^\top Z*) from newdata
**
** Output
**    Haz[nt*K] - NEW: OUTPUT, sort of what used to be surv in agsurv2, dimension has changed
**    varHaz[nt*K*(K+1)/2] - NEW: OUTPUT, like varh in agsurv2, variance and covariance of estimated Hazard at unt
**
**  Work
**    d[3*p] - WORKING, contains d, a, a2
**    work[nt*K*(p+1)] - NEW: WORKING
**
**  Input must be sorted by (event before censor) within stop time within strata,
*/
#include <math.h>
#include <R.h>

void agmssurv(
			int		*sn, /* n, no of subjects */
			int		*sp, /* p, no of vars in xmat */
			int     *svar, /* var, NEW: 1 if variances are to be calculated, 0 otherwise */
			int		*smethod, /* method (for handling ties), CHANGED: 1=breslow, 2=efron */
			int		*sH, /* H, NEW: no of strata in coxph call */
			int		*sK, /* K, no of curves (= no of transitions) */
			int		*snt, /* nt, NEW: no of unique time points joint from all strata */
			double	*y, /* matrix containing (start, stop, event) */
			double	*score, /* vector of weights exp(beta^T Z_i) */
			double	*xmat, /* data matrix that generated the Cox fit */
			double	*varcov, /* inverse Fisher information for estimation of beta */
			int		*strata, /* CHANGED: strata[0] = 0, strata[1:H]= last obs strata 1,2, etc in y */
			int		*kstrata, /* NEW: kstrata[k] = stratum in newdata[k,] */
			double	*unt, /* NEW: vector of unique time points joint from all strata */
			double	*newx, /* matrix with covariates Z for each individual */
			double	*newrisk, /* vector containing exp(\beta^\top Z*) from newdata */
			double	*Haz, /* NEW: OUTPUT, sort of what used to be surv in agsurv2, dimension has changed */
			double	*varHaz, /* NEW: OUTPUT, like varh in agsurv2, variance and covariance of estimated Hazard at unt */
			double	*d, /* WORKING, contains d, a, a2 */
			double	*work /* NEW, WORKING */
			)
{
    int i,j,k,l;
    double hazard, varhaz;
    double *start, *stop, *event;
    int n, p, var, H, K, nt, method;
    double *a, *a2, *tmp, *eta;
    int k1, k2, k12, K12;
    double *covar, *imatinv, *covar2;
    int idx, thestrat, person;
    double time, weight=0, denom, e_denom;
    double crisk, deaths;
    double temp, downwt, d2;

/*	Rprintf("Entering agmssurv ...\n");
	R_FlushConsole(); */

    n = *sn;  p = *sp; var = *svar, H = *sH, K = *sK, nt = *snt; method = *smethod;
    start = y;
    stop = y+n;
    event = y+n+n;
    a = d+p;
    a2 = a+p;
    tmp = work;
    eta = work + nt*K;

/*	Rprintf("n = %d, p = %d, var = %d, H = %d, K = %d, nt = %d, method = %d\n\n",n,p,var,H,K,nt,method);
	R_FlushConsole(); */

    /*
    **  Set up vectors
    */
    covar = xmat; /* dmatrix(xmat, n, p) */
    imatinv = varcov; /* dmatrix(varcov, p, p) */
    covar2 = newx; /* dmatrix(newx, K, p) */

    for (k=0; k<K; k++) {
		/*
		** For-loop to store Haz, and two components to calculate varHaz:
		** - tmp: first element of (12) in de Wreede et al. (2009)
		** - eta: part second element of (12) in de Wreede et al. (2009)
		**        that comes before I^{-1}
		*/
		idx = 0;
		thestrat = kstrata[k];
//		Rprintf("\n\nk = %d, thestrat = %d\n\n\n",k,thestrat);
		crisk = newrisk[k];
		hazard = 0;
		varhaz = 0;
		for (j=0; j<p; j++) d[j] = 0;
		for (person=strata[thestrat-1]; person<strata[thestrat];) {
//			Rprintf("person = %d, event[person] = %d\n",person,event[person]);
			if (event[person]==0) person++;
			else {
				/*
				** compute the mean and denominator over the risk set
				*/
				denom = 0;
				e_denom = 0;
				for (j=0; j<p; j++) {
					a[j] = 0;
					a2[j] = 0;
				}
				time = stop[person];
				deaths = 0;
//				Rprintf("time = %6.4f\n",time);
				/* propogate forward last hazard, varhaz, d, to Haz, tmp, eta */
				while (unt[idx] < time) {
					Haz[nt*k + idx] = hazard;
					if (var==1) {
						tmp[nt*k + idx] = varhaz;
						for (j=0; j<p; j++)
							eta[nt*p*k + p*idx + j] = d[j];
//						Rprintf("idx=%d, unt[idx]=%6.4f, time=%6.4f\n",idx,unt[idx],time);
					}
					idx++;
				}
				for (i=person; i<strata[thestrat]; i++) {
					if (start[i] < time) {
						weight = score[i]/crisk;
						denom += weight;
						for (j=0; j<p; j++) {
							a[j] += weight*(covar[n*j+i]- covar2[K*j+k]);
							/* check */
//							Rprintf("k=%d, i=%d, j=%d, covar[n*j+i]=%6.4f, covar2[K*j+k]=%6.4f\n",k,i,j,covar[n*j+i],covar2[K*j+k]);
						}
					}
					if (stop[i]==time && event[i]==1) {
						deaths += 1;
						e_denom += weight;
						for (j=0; j<p; j++) {
							a2[j] += weight*(covar[n*j+i]- covar2[K*j+k]);
							/* check */
//							Rprintf("k=%d, i=%d, j=%d, covar[n*j+i]=%6.4f, covar2[K*j+k]=%6.4f\n",k,i,j,covar[n*j+i],covar2[K*j+k]);
						}
					}
				}
				/*
				** Add results of all events for this time point
				*/
				/* Since modeled after survfit, which calls agsurv2.c, use of names
				is close to Therneau and Grambsch. Here is what the terms used in the program
				correspond to in terms of the notation of de Wreede et al. (2009).
				(All hats in betahat are suppressed.)
				denom = S_q^(0)(t)/exp(beta^T Z*); Hazard += exp(beta^T Z*)/S_q^(0)(t) dN_q(t)
				varHaz += (exp(beta^T Z*)/S_q^(0)(t))^2 dN_q(t)
				a = (S_q^(1)(t) - Z* S_q^(0)(t))/exp(beta^T Z*)
				d += a/denom^2 = (S_q^(1)(t)/S_q^(0)(t) - Z*) exp(beta^T Z*) dN_q(t)/S_q^(0)(t) */
				temp = 0;
				for (i=person; i<strata[thestrat] && stop[i]==time; i++) {
					if (event[i]==1) {
						downwt = temp/deaths;
						if (method==2) { /* efron */
							d2 = (denom - downwt*e_denom);
							hazard += 1/d2;
						}
						else hazard += 1/denom; /* breslow */
						if (var==1) {
							if (method==2) { /* efron */
								d2 = (denom - downwt*e_denom);
								varhaz += 1/(d2*d2);
								for (j=0; j<p; j++)
									d[j] += (a[j] - downwt*a2[j])/(d2*d2);
							}
							else { /* breslow */
								varhaz += 1/(denom*denom);
								for (j=0; j<p; j++)
									d[j] += a[j]/(denom*denom);
							}
						}
						temp++;
					}
					person++;
				}
			}
		}
		/* After running through all individuals in the stratum, store results in Haz, tmp, eta */
		for (i=idx; i<nt; i++) {
			Haz[nt*k + i] = hazard;
			if (var==1) {
				tmp[nt*k + i] = varhaz;
				for (j=0; j<p; j++)
					eta[nt*p*k + p*i + j] = d[j];
			}
		}
		/* Check */
/*		for (i=0; i<nt; i++) {
			Rprintf("k=%d, i=%d, Haz[nt*k+i]=%6.4f, tmp[nt*k+i]=%6.4f, eta[nt*p*k+p*i+j]=:\n",k,i,Haz[nt*k+i],tmp[nt*k+i]);
			for (j=0; j<p; j++)
				Rprintf("  j=%d: %6.4f",j,eta[nt*p*k+p*i+j]);
			Rprintf("\n");
		} */
	}
	if (var==1) {
//		Rprintf("\n\n\n\n");
		k12 = 0; K12 = K*(K+1)/2;
		for (k1=0; k1<K; k1++) {
			for (k2=k1; k2<K; k2++) {
				for (i=0; i<nt; i++) {
					varHaz[K12*i+k12] = 0.0;
					if (kstrata[k1]==kstrata[k2]) varHaz[K12*i+k12] = tmp[nt*k1+i];
//					Rprintf("k=%d, i=%d, k1=%d, k2=%d, k12=%d, kstrata[k1]=%d, kstrata[k2]=%d\n",k,i,k1,k2,k12,kstrata[k1],kstrata[k2]);
//					Rprintf("\ttmp[nt*k1+i]=%6.4f, varHaz[K12*i+k12]=%6.4f\n",tmp[nt*k1+i],varHaz[K12*i+k12]);
					for (j=0; j<p; j++)
						for (l=0; l<p; l++)
							varHaz[K12*i+k12] += eta[nt*p*k1+p*i+j]*eta[nt*p*k2+p*i+l]*imatinv[p*j+l];
//					Rprintf("k=%d, i=%d, k1=%d, k2=%d, k12=%d, varHaz[K12*i+k12]=%6.4f\n",k,i,k1,k2,k12,varHaz[K12*i+k12]);
				}
				++k12;
			}
		}
	}
}
