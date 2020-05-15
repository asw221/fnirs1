/*
 *  Marginal LFD.c
 *  LGCP
 *
 *  Created by Timothy D. Johnson on 04/16/18.
 *  Copyright 2015 University of Michigan. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "randgen.h"
#include "cholesky.h"
#include <time.h>
#include "fNIRS.h"

extern int MAX_ITER;
extern int BURN_IN;
extern sDLM *dlmStruc;
extern int maxP;
extern FILE *flog;

inline double marg_ll_t(double res,double prior_var,double data_var,double *basis,const int dim) 
{
    double V;
    double dotProd(double *x,double *y,const int dim);
    
    V = prior_var*dotProd(basis,basis,dim)+data_var;
    return log(V) + res*res/V;
}

double marg_ll(double *res,double prior_var,double *data_prec,double *V,const int dim,const int N)
{
    double ll = 0;
    
    for (int t=0;t<N;t++)
        ll += marg_ll_t(res[t],prior_var,1./data_prec[t],&(V[t*dim]),dim);
    
    return -0.5*ll;
}

inline void marg_res(double *res,double *Y,double *Xb,double *delta,const int N,const int P)
{
    double Wd;
    
    for (int t=0;t<P;t++)
        res[t] = Y[t] - Xb[t];
    
    for (int t=P;t<N;t++) {
        Wd = 0;
        for (int i=0;i<P;i++)
            Wd +=  Y[t-i-1]*delta[t*P+i];
        res[t] = Y[t] - Xb[t] - Wd;
    }
}

void remove_knot2(REP *rep,double *V,double *knots,const int sdegree,int remove)
{
    int i,j,k;
    double *data;
    void bsplinebasis(int,int,int,double,double *,int,double *);
    
    k = 0;
    for (i=0;i<rep->nKnots;i++) {
        if (i != remove) {
            knots[k] = rep->knots[i];
            k++;
        }
    }
    
    data = (double *)calloc(rep->dim_V[0],sizeof(double));
    for (i=0;i<rep->dim_V[0];i++)
        data[i] = i;
 
    for (i=0;i<rep->dim_V[0];i++) {
        bsplinebasis(0,0,sdegree,data[i],knots,rep->nKnots-1,&V[i*(rep->dim_V[1]-1)]);
    }
    
    free(data);
}

void remove_basis(POP *pop,REP *rep,int sdegree,unsigned long *seed)
{
    double *res,*V,*knots;
    void calH(double *H,double *A,double *m,const int nrow,const int ncol,const int P); 
    void adddoubles(double *C,double *A,double *B,const int dim,const int sign);
    void copydouble(double *B,double *A,const int dim);
    
   if (rep->nKnots-2*sdegree > 0) {
        double vareta = 1./rep->preceta;
    
        V = (double *)calloc(rep->dim_V[1]*rep->dim_V[0],sizeof(double)); 
        res = (double *)calloc(rep->dim_V[0],sizeof(double));
    
        calH(V,rep->V,rep->delta,(const int)rep->dim_V[0],(const int)rep->dim_V[1],(const int)rep->dim_W[1]);
        adddoubles(V,rep->V,V,(const int)(rep->dim_V[0]*rep->dim_V[1]),-1);
 
        marg_res(res,rep->Y,rep->Xbeta,rep->delta,(const int)rep->dim_X[0],(const int)rep->dim_W[1]);
    
        double old_ll = marg_ll(res,vareta,rep->d_Y,V,(const int)rep->dim_V[1],(const int)rep->dim_V[0]);
        free(V);
    
        int remove = runiform_n(rep->nKnots-2*sdegree,seed) + sdegree;
        V = (double *)calloc((rep->dim_V[1]-1)*rep->dim_V[0],sizeof(double));
        knots = (double *)calloc(rep->nKnots-1,sizeof(double));
        remove_knot2(rep,V,knots,(const int)sdegree,remove);
    
        double new_ll = marg_ll(res,vareta,rep->d_Y,V,(const int)(rep->dim_V[1]-1),(const int)rep->dim_V[0]);

        double ll_ratio = new_ll - old_ll;
        int nKnots = rep->nKnots-2*sdegree;
        double lprior_ratio = log((double)nKnots/rgamma((double)pop->knots,1,seed)) + log((double)rep->dim_X[0]);
        double lprop_ratio = log((double)nKnots/(double)rep->dim_X[0]);
        
        if (log(kiss(seed)) < ll_ratio + lprior_ratio + lprop_ratio) {
            (rep->dim_V[1])--;
            (rep->nKnots)--;

            rep->eta = (double *)realloc(rep->eta,sizeof(double)*(rep->dim_V[1]));
            rep->knots = (double *)realloc(rep->knots,sizeof(double)*(rep->nKnots));
            copydouble(rep->knots,knots,(const int)rep->nKnots);
    
            rep->V = (double *)realloc(rep->V,rep->dim_V[0]*rep->dim_V[1]*sizeof(double));
            copydouble(rep->V,V,(const int)(rep->dim_V[0]*rep->dim_V[1]));  
        }
        free(knots);
        free(res);
        free(V);
   }
}

void add_knot2(REP *rep,double *V,double *knots,int sdegree,double position)
{
    int ii,i,j,k,flag,pos,nKnots;
    double *data;
    void bsplinebasis(int,int,int,double,double *,int,double *);
    
    flag = 1;
    for (i=0;i<rep->nKnots;i++) {
        if (rep->knots[i] < position) {
            knots[i] = rep->knots[i];
        }
        else
            break;
    }    
    knots[i] = position;
    for (j=i;j<rep->nKnots;j++)
        knots[j+1] = rep->knots[j];
     
    data = (double *)calloc(rep->dim_V[0],sizeof(double));
    for (i=0;i<rep->dim_V[0];i++)
        data[i] = i;
    
    nKnots = rep->nKnots+1;
    for (i=0;i<rep->dim_V[0];i++)
        bsplinebasis(0,0,sdegree,data[i],knots,nKnots,&V[i*(rep->dim_V[1]+1)]);
    
    free(data);
}

void add_basis(POP *pop,REP *rep,int sdegree,unsigned long *seed)
{
    int MAX_KNOTS = 100;
    double *res,*V,*knots;
    void calH(double *H,double *A,double *m,const int nrow,const int ncol,const int P); 
    void adddoubles(double *C,double *A,double *B,const int dim,const int sign);
    void copydouble(double *B,double *A,const int dim);
   
   if (rep->nKnots-2*sdegree < MAX_KNOTS) {
        double vareta = 1./rep->preceta;
    
        V = (double *)calloc(rep->dim_V[1]*rep->dim_V[0],sizeof(double)); 
        res = (double *)calloc(rep->dim_V[0],sizeof(double));
    
        calH(V,rep->V,rep->delta,(const int)rep->dim_V[0],(const int)rep->dim_V[1],(const int)rep->dim_W[1]);
        adddoubles(V,rep->V,V,(const int)(rep->dim_V[0]*rep->dim_V[1]),-1);
 
        marg_res(res,rep->Y,rep->Xbeta,rep->delta,(const int)rep->dim_X[0],(const int)rep->dim_W[1]);
    
        double old_ll = marg_ll(res,vareta,rep->d_Y,V,(const int)rep->dim_V[1],(const int)rep->dim_V[0]);
        free(V);
    
        double  position = runif_atob(seed,0,rep->dim_V[0]-1);
        for (int i = sdegree;i<rep->nKnots-sdegree;i++){
            if (fabs(position-rep->knots[i]) < 0.01) {
                position = runif_atob(seed,0,rep->dim_V[0]-1);
                i = sdegree;
            }
        }

        V = (double *)calloc((rep->dim_V[1]+1)*rep->dim_V[0],sizeof(double));
        knots = (double *)calloc(rep->nKnots+1,sizeof(double));
        add_knot2(rep,V,knots,sdegree,position);
     
        double new_ll = marg_ll(res,vareta,rep->d_Y,V,(const int)(rep->dim_V[1]+1),(const int)rep->dim_V[0]);

        double ll_ratio = new_ll - old_ll;
        int nKnots = rep->nKnots-2*sdegree + 1;
        double lprior_ratio = log(rgamma((double)pop->knots,1,seed)/(double)nKnots) - log((double)rep->dim_X[0]);
        double lprop_ratio = log((double)rep->dim_X[0]/(double)nKnots);
        
        if (log(kiss(seed)) < ll_ratio + lprior_ratio + lprop_ratio) {
            (rep->dim_V[1])++;
            (rep->nKnots)++;

            rep->eta = (double *)realloc(rep->eta,sizeof(double)*(rep->dim_V[1]));
            rep->knots = (double *)realloc(rep->knots,sizeof(double)*(rep->nKnots));
            copydouble(rep->knots,knots,(const int)rep->nKnots);
    
            rep->V = (double *)realloc(rep->V,rep->dim_V[0]*rep->dim_V[1]*sizeof(double));
            copydouble(rep->V,V,(const int)(rep->dim_V[0]*rep->dim_V[1]));  
        }
        free(knots);
        free(res);
        free(V);
   }
}

