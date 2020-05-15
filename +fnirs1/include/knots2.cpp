//
//  knots.cpp
//  
//
//  Created by Timothy Johnson on 10/14/16.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "randgen.h"
#include "cholesky.h"
#include <time.h>
#include "fNIRS.h"

extern sDLM *dlmStruc;

void draw_knot_locations(REP *rep,int sdegree,int *flag,unsigned long *seed)
{
    int i,j,k,l,fl2=0;
    static double old_like,new_like;
    double like_ratio,prior_ratio,*data,*knots;
    double *V,*knots2,tmp,max,*eta,*eta2;
    void bsplinebasis(int j,int k,int SplineOrder,double x,double *knots,int Nknots,double *b);
    void calculate_marginal_residuals(REP *rep,int P);
    void calculate_residuals(REP *rep,int P);
    void calAx(double *Ax,double *A,double *x,const int nrow,const int ncol);   
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);    
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calApVinvA(double *ApA,double **A,double *Vinv,const int nrow,const int ncol);
    void dlm_forward_filter2(REP *rep,sDLM *dlmStruc);  

    knots = (double *)calloc(rep->nKnots,sizeof(double));
    knots2 = (double *)calloc(rep->nKnots,sizeof(double));
    
    double *Qinv = (double *)calloc(rep->dim_V[0],sizeof(double));

    // calc current likelihood
    if (*flag == 0) {
        *flag = 1;
        dlm_forward_filter2(rep,dlmStruc);
        calWdelta(rep->Wm,rep->W,rep->m,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
        calculate_marginal_residuals(rep,rep->P);
        for (i=0;i<rep->P;i++)
            Qinv[i] = rep->d_Y[i];
        for (i=rep->P;i<rep->dim_V[0];i++)
            Qinv[i] = 1./dlmStruc[i].Q;

        old_like = 0;
        for (i=0;i<rep->dim_X[0];i++)
            old_like += log(Qinv[i])+Qinv[i]*rep->residuals[i]*rep->residuals[i];

    }
    data = (double *)calloc(rep->dim_V[0],sizeof(double));
    for (i=0;i<rep->dim_V[0];i++)
        data[i] = i;
    
    for (i=0;i<rep->nKnots;i++)
        knots[i] = rep->knots[i];
    
    V = (double *)calloc(rep->dim_V[0]*rep->dim_V[1],sizeof(double));
    eta = (double *)calloc(rep->dim_V[1],sizeof(double));
    eta2 = (double *)calloc(rep->dim_V[1],sizeof(double));
    for (i=0;i<rep->dim_V[1];i++)
        eta[i] = rep->eta[i];
    
           // attempt change of each interior knot in order
    j = runiform_n(rep->nKnots-sdegree*2,seed) + sdegree;    

        knots[j] = rnorm(knots[j],rep->prop_sd[3],seed);  // propose knot location
        for (i = sdegree;i<rep->nKnots-sdegree;i++){
            if (!(i==j)) {
                if (fabs(knots[j]-knots[i]) < 0.01) {
                    free(V);
                    free(eta);
                    free(eta2);
                    free(knots);
                    free(knots2);
                    free(data);
                    return;
                }
            } 
        }

        (rep->attempt[3])++;

        if (knots[j] > 0 && knots[j] < rep->dim_V[0]-1) {  // ensure knot is still interior knot
            
            for (i=0;i<rep->nKnots;i++)
                knots2[i] = knots[i];
            for (i=0;i<rep->dim_V[1];i++)
                eta2[i] = eta[i];
            
            for (i=0;i<rep->nKnots-sdegree-1;i++) {  // reorder knots so their locations are sequential
                for (k=i+1;k<rep->nKnots-sdegree;k++) {
                    if (knots2[k] < knots2[i]) {
                        tmp = knots2[k];
                        knots2[k] = knots2[i];
                        knots2[i] = tmp; 
                        tmp = eta2[k];
                        eta2[k] = eta2[i];
                        eta2[i] = tmp;                     
                    }
                }
            }
            
            for (i=0;i<rep->dim_V[0];i++)
                for (k=0;k<rep->dim_V[1];k++)
                    V[i*rep->dim_V[1]+k] = 0;

            for (i=0;i<rep->dim_V[0];i++) {  // create new bspline basis
                bsplinebasis(0,0,sdegree,data[i],knots2,rep->nKnots,&V[i*rep->dim_V[1]]);
            }
 
           // calculate Veta, calculate residuals
    
            calAx(rep->Veta,V,eta2,(const int)rep->dim_V[0],(const int)(rep->dim_V[1]));
            calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],rep->P);
            dlm_forward_filter2(rep,dlmStruc); 
            calWdelta(rep->Wm,rep->W,rep->m,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
            int P = rep->dim_W[1];
          
            calculate_marginal_residuals(rep,rep->P);

            // calculate new log likelihood
            for (i=0;i<rep->P;i++)
                Qinv[i] = rep->d_Y[i];
            for (i=rep->P;i<rep->dim_V[0];i++)
                Qinv[i] = 1./dlmStruc[i].Q;
            new_like = 0;
            for (i=0;i<rep->dim_X[0];i++)
                new_like += log(Qinv[i]) + Qinv[i]*rep->residuals[i]*rep->residuals[i];
           
            like_ratio = -0.5*(new_like-old_like);
           
            prior_ratio = 0;//log(prb_new/prb_old);
            
  //          printf("%lf %lf %lf\n",kiss(seed),like_ratio,prior_ratio);
            if (log(kiss(seed)) < like_ratio + prior_ratio) { // accept new location
                old_like = new_like;
                (rep->accept[3])++;
                rep->knots[j] = knots[j];
//                rep->eta[j] = eta[j];
                old_like = new_like;
                fl2 = 1;
              }
            else {   // reject new location
                knots[j] = rep->knots[j];
//                eta[j] = rep->eta[j];
            }
        }
        else { // reject new location, not an interior knot
            knots[j] = rep->knots[j];
//            eta[j] = rep->eta[j];
         }
//    }
    
    //  save basis, Veta, residuals
    if (fl2) {
        for (i=0;i<rep->nKnots-sdegree-1;i++) {
            for (k=i+1;k<rep->nKnots-sdegree;k++) {
                if (rep->knots[k] < rep->knots[i]) {
                    tmp = rep->knots[k];
                    rep->knots[k] = rep->knots[i];
                    rep->knots[i] = tmp;
                    tmp = rep->eta[k];
                    rep->eta[k] = rep->eta[i];
                    rep->eta[i] = tmp;
                }
            }
        }
//printf("nKnots - sdegree = %d dim_V[1] = %d\n",rep->nKnots-sdegree,rep->dim_V[1]);
        for (i=0;i<rep->dim_V[0];i++)
            for (k=0;k<rep->dim_V[1];k++)
                rep->V[i*rep->dim_V[1]+k] = 0;

        for (i=0;i<rep->dim_V[0];i++) {  // create new bspline basis
            bsplinebasis(0,0,sdegree,data[i],rep->knots,rep->nKnots,&rep->V[i*rep->dim_V[1]]);
        }
 
    }
    else {
    calAx(rep->Veta,rep->V,rep->eta,(const int)rep->dim_V[0],(const int)(rep->dim_V[1]));

    calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],rep->P);
//    dlm_forward_filter2(rep,dlmStruc); 
//    calWdelta(rep->Wm,rep->W,rep->m,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
//    calWdelta(rep->Wdelta,rep->W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
//    calculate_residuals(rep,rep->P);
    //calculate_marginal_residuals(rep,rep->P);
    }
    // free allocated scratch memory
    free(Qinv);
    free(V);
    free(eta);
    free(eta2);   
    free(knots);
    free(knots2);
    free(data);
}

int calculate_death_rates(double *death_rate,double *Death_rate,double Birth_rate,double prior_rate,REP *rep,POP *pop,int sdegree,double full_likelihood) {
    int i,j,k;
    double *partial_likelihood;
    double max;
    double *knots,*eta,*V;
    double *Veta,*W,*Wm;
    
    void calculate_marginal_residuals(REP *,int);
    void remove_knot(REP *rep,double *V,double *eta,double *knots,const int sdegree,int remove);
    void calAx(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void knot_death_rate(double *,REP *,double *,double,double,double,int);
    void dlm_forward_filter2(REP *rep,sDLM *dlmStruc); 
 
    partial_likelihood = (double *)calloc(rep->nKnots-2*sdegree,sizeof(double));
     
    double *Qinv = (double *)calloc(rep->dim_V[0],sizeof(double));
    double *tmpVeta = rep->Veta;
   
    if (rep->nKnots-2*sdegree > 0) {
        knots = (double *)calloc(rep->nKnots-1,sizeof(double));
        eta = (double *)calloc(rep->dim_V[1]-1,sizeof(double));
        V = (double *)calloc(rep->dim_V[0]*(rep->dim_V[1]-1),sizeof(double));
        Veta = (double *)calloc(rep->dim_V[0],sizeof(double));
        W = (double *)calloc(rep->dim_W[0]*rep->dim_W[1],sizeof(double));
        Wm = (double *)calloc(rep->dim_W[0],sizeof(double)); 
//        for (i=0;i<rep->dim_V[0];i++)
//            V[i] = (double *)calloc(rep->dim_V[1]-1,sizeof(double));
        
        for (i=0;i<rep->nKnots-2*sdegree;i++) {
            
            for (j=0;j<rep->dim_V[0];j++)
                for (k=0;k<rep->dim_V[1]-1;k++)
                    V[j*(rep->dim_V[1]-1)+k] = 0;
            
               /* remove each knot in turn and calculate partial likelihoods */
            remove_knot(rep,V,eta,knots,4,i);
            
            /* calculate Veta */
            
            calAx(Veta,V,eta,(const int)rep->dim_V[0],(const int)(rep->dim_V[1]-1));
            rep->Veta = Veta;
            /* update W and Wdelta */
            calW(rep->W,rep->Y,rep->Xbeta,Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],rep->P);
  
            dlm_forward_filter2(rep,dlmStruc);
            calWdelta(rep->Wm,rep->W,rep->m,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
            
            calculate_marginal_residuals(rep,rep->P);
            
            for (j=0;j<rep->P;j++)
                Qinv[j] = rep->d_Y[j];
            for (j=rep->P;j<rep->dim_V[0];j++)
                Qinv[j] = 1./dlmStruc[j].Q;
            partial_likelihood[i] = 0;
            for (j=0;j<rep->dim_X[0];j++) {
                partial_likelihood[i] += log(Qinv[j]) + Qinv[j]*rep->residuals[j]*rep->residuals[j];
            }
            partial_likelihood[i] *= -0.5;
        }
        
        free(knots);
        free(eta);
//        for (i=0;i<rep->dim_V[0];i++)
//            free(V[i]);
        free(V);
        free(Veta);
        free(W);
        free(Wm);
    }
    free(Qinv);
    rep->Veta = tmpVeta;
    
    /* CALCULATE DEATH RATE FOR EACH COMPONENT */
    int dflag = 0;
    
    if (rep->nKnots-2*sdegree > 0) {
        dflag = 1;
        knot_death_rate(death_rate,rep,partial_likelihood,full_likelihood,Birth_rate,
                                     prior_rate,sdegree);
    }
    if (dflag) {
        *Death_rate = death_rate[0];
        for (i=1;i<rep->nKnots-2*sdegree;i++) {
            max = (*Death_rate > death_rate[i]) ? *Death_rate:death_rate[i];
            *Death_rate = log(exp(*Death_rate-max) + exp(death_rate[i]-max)) + max;
        }
        
        for (i=0;i<rep->nKnots-2*sdegree;i++) {
            death_rate[i] -= *Death_rate;
            death_rate[i] = exp(death_rate[i]);
            
        }
        for (i=1;i<rep->nKnots-2*sdegree;i++)
            death_rate[i] += death_rate[i-1];
        *Death_rate = exp(*Death_rate);
    }
    else
        *Death_rate = 0;

    free(partial_likelihood);

    return dflag;
}

void birth(REP *rep,POP *pop,double *full_likelihood,unsigned long *seed) {
    int i,j,pos_idx,sdegree=4;
    double position,full_ldens,new_eta;
    double *teta,*tknots,*VV;
    
    int add_knot(REP *rep,double *V,double *knots,double *eta,int sdegree,double new_eta,double position);
    void calculate_marginal_residuals(REP *rep,int P);
    void calAx(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void dlm_forward_filter2(REP *,sDLM *);
    
    double *Qinv = (double *)calloc(rep->dim_V[0],sizeof(double));

    position = runif_atob(seed,0,rep->dim_V[0]-1);
    for (i = sdegree;i<rep->nKnots-sdegree;i++){
        if (fabs(position-rep->knots[i]) < 0.1) {
                    return;
         }
    }
     /* simulate eta from its full cond */
 
    VV = (double *)calloc(rep->dim_V[0]*(rep->dim_V[1]+1),sizeof(double));
//    for (i=0;i<rep->dim_V[0];i++)
//        VV[i] = (double *)calloc(rep->dim_V[1]+1,sizeof(double));
    
    teta = (double *)calloc(rep->dim_V[1]+1,sizeof(double));
    tknots = (double *)calloc(rep->nKnots+1,sizeof(double));
    
 //   free(rep->eta);
    pos_idx = add_knot(rep,VV,tknots,teta,4,new_eta,position);
    
    teta[pos_idx] = rnorm(0,1/sqrt(rep->preceta),seed);
  
    (rep->dim_V[1])++;
    (rep->nKnots)++;

    rep->eta = (double *)realloc(rep->eta,sizeof(double)*(rep->dim_V[1]));
//    rep->eta = (double *)calloc(rep->dim_V[1]+1,sizeof(double));
    for (i=0;i<rep->dim_V[1];i++)
        rep->eta[i] = teta[i];
    free(teta);
    
//    free(rep->knots);
//    rep->knots = (double *)calloc(rep->nKnots+1,sizeof(double));
    rep->knots = (double *)realloc(rep->knots,sizeof(double)*(rep->nKnots));
    for (i=0;i<rep->nKnots;i++)
        rep->knots[i] = tknots[i];
    free(tknots);
    
    
    
//    for (i=0;i<rep->dim_V[0];i++)
//        free(rep->V[i]);
    //      free(rep->V);
    
    //      rep->V = (double **)calloc(rep->dim_V[0],sizeof(double *));
//    for (i=0;i<rep->dim_V[0];i++)
        rep->V = (double *)realloc(rep->V,rep->dim_V[0]*rep->dim_V[1]*sizeof(double));
//    rep->V[i] = (double *)calloc(rep->dim_V[1],sizeof(double));
    
    //           FILE *fouta; fouta = fopen("VVn.out","w");
    for (i=0;i<rep->dim_V[0];i++) {
        for (j=0;j<rep->dim_V[1];j++) {
            rep->V[i*rep->dim_V[1]+j] = VV[i*rep->dim_V[1]+j];
            //                   fprintf(fouta,"%lf ",VV[i][j]);
        }
        //              fprintf(fouta,"\n");
    }
    //         fclose(fouta);
    
 //   for (i=0;i<rep->dim_V[0];i++)
 //       free(VV[i]);
    free(VV);
    
//    for (i=0;i<rep->dim_V[0];i++)
//        free(rep->J[i]);
//    free(rep->J);
    
//    rep->J = (double **)calloc(rep->dim_V[0],sizeof(double *));
 //   for (i=0;i<rep->dim_V[0];i++)
 //       rep->J[i] = (double *)realloc(rep->J[i],rep->dim_V[1]*sizeof(double));
//    rep->J[i] = (double *)calloc(rep->dim_V[1],sizeof(double));
    
    /* calculate Veta */
    
    calAx(rep->Veta,rep->V,rep->eta,(const int)rep->dim_V[0],(const int)rep->dim_V[1]);
    
    /* update W and Wdelta */
    calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],rep->P);
    dlm_forward_filter2(rep,dlmStruc);
    calWdelta(rep->Wm,rep->W,rep->m,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
    
    calculate_marginal_residuals(rep,rep->P);
    
    for (i=0;i<rep->P;i++)
        Qinv[i] = rep->d_Y[i];
    for (i=rep->P;i<rep->dim_V[0];i++)
        Qinv[i] = 1./dlmStruc[i].Q;
   
    *full_likelihood = 0;
    for (i=0;i<rep->dim_X[0];i++)
        *full_likelihood += log(Qinv[i]) + Qinv[i]*rep->residuals[i]*rep->residuals[i];
    *full_likelihood *= -0.5;
    
    free(Qinv);
}

void death(double *death_rate,REP *rep,POP *pop,double *full_likelihood,int sdegree,unsigned long *seed) {
    int i,j,remove;
    double *eta,*tknots;
    
    void remove_knot(REP *rep,double *V,double *eta,double *knots,const int sdegree,int remove);
    void calculate_marginal_residuals(REP *rep,int P);
    void calAx(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void dlm_forward_filter2(REP*,sDLM*);

    double *Qinv = (double *)calloc(rep->dim_V[0],sizeof(double));

    remove = (int)rmultinomial(death_rate,(long)(rep->nKnots-2*sdegree),seed);
//    for (i=0;i<rep->dim_V[0];i++)
//        free(rep->V[i]);
//    free(rep->V);
 
//    rep->V = (double **)calloc(rep->dim_V[0],sizeof(double *));
 //   for (i=0;i<rep->dim_V[0];i++)
        rep->V = (double *)realloc(rep->V,rep->dim_V[0]*(rep->dim_V[1]-1)*sizeof(double));
//    rep->V[i] = (double *)calloc(rep->dim_V[1]-1,sizeof(double));
    
//    for (i=0;i<rep->dim_V[0];i++)
//        free(rep->J[i]);
//    free(rep->J);
    
//    rep->J = (double **)calloc(rep->dim_V[0],sizeof(double *));
//    for (i=0;i<rep->dim_V[0];i++)
//        rep->J[i] = (double *)calloc(rep->dim_V[1]-1,sizeof(double));
    
    eta = (double *)calloc(rep->dim_V[1]-1,sizeof(double));
    
    tknots = (double *)calloc(rep->nKnots-1,sizeof(double));
    remove_knot(rep,rep->V,eta,tknots,4,remove);
    (rep->nKnots)--;
    (rep->dim_V[1])--;
    
//    free(rep->eta);
//    rep->eta = (double *)calloc(rep->dim_V[1],sizeof(double));
    rep->eta = (double *)realloc(rep->eta,rep->dim_V[1]*sizeof(double));
    for (i=0;i<rep->dim_V[1];i++)
        rep->eta[i] = eta[i];
    free(eta);
    
//    free(rep->knots);
//    rep->knots = (double *)calloc(rep->nKnots,sizeof(double));
    rep->knots = (double *)realloc(rep->knots,rep->nKnots*sizeof(double));
    for (i=0;i<rep->nKnots;i++)
        rep->knots[i] = tknots[i];
    free(tknots);
    
    /* calculate Veta */
    
    calAx(rep->Veta,rep->V,rep->eta,(const int)rep->dim_V[0],(const int)rep->dim_V[1]);
    
    /* update W and Wdelta */
    calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],rep->P);
    dlm_forward_filter2(rep,dlmStruc);
    calWdelta(rep->Wm,rep->W,rep->m,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
    
    calculate_marginal_residuals(rep,rep->P);
    
    for (i=0;i<rep->P;i++)
        Qinv[i] = rep->d_Y[i];
    for (i=rep->P;i<rep->dim_V[0];i++)
        Qinv[i] = 1./dlmStruc[i].Q;
    *full_likelihood = 0;
    for (i=0;i<rep->dim_X[0];i++)
        *full_likelihood += log(Qinv[i]) + Qinv[i]*rep->residuals[i]*rep->residuals[i];
    *full_likelihood *= -0.5;
    
    free(Qinv);
}

int knot_birth_death(REP *rep,POP *pop,const int sdegree,int iter,unsigned long *seed){
    
    int i,j,k,aaa,dflag,max_knots = 100;
    double Birth_rate,prior_rate,T,full_likelihood,max,S;
    double Death_rate,*death_rate,position,*partial_likelihood;
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calAx(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calculate_residuals(REP *rep,int P);
    void dlm_forward_filter2(REP*,sDLM*);
    
    void calculate_marginal_residuals(REP *rep,int P);
   
    double *Qinv = (double *)calloc(rep->dim_V[0],sizeof(double));

    aaa = 0;
    
    // set simulation time

    prior_rate = pop->knots;//rgamma(pop->knots,1,seed); // prior on number of knots given prior_rate is Poisson(prior_rate).  This results in a neg binomial prior.
    Birth_rate = prior_rate;
//    T= 1./Birth_rate;
    T = 1;
//    calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],rep->P);
    dlm_forward_filter2(rep,dlmStruc);
    calWdelta(rep->Wm,rep->W,rep->m,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
    calculate_marginal_residuals(rep,rep->P);
    
    // calculate likelihood;
    for (i=0;i<rep->P;i++)
        Qinv[i] = rep->d_Y[i];
    for (i=rep->P;i<rep->dim_V[0];i++)
        Qinv[i] = 1./dlmStruc[i].Q;
    full_likelihood = 0;
    for (i=0;i<rep->dim_X[0];i++)
        full_likelihood += log(Qinv[i]) + Qinv[i]*rep->residuals[i]*rep->residuals[i];
    full_likelihood *= -0.5;
    S = 0.0;

    while (1) {
        aaa++;
  
        death_rate = (double *)calloc(rep->nKnots-2*sdegree,sizeof(double));

        dflag = calculate_death_rates(death_rate,&Death_rate,Birth_rate,prior_rate,rep,pop,sdegree,full_likelihood);
        
         /* SIMULATE TIME, S, TO NEXT JUMP */
        if (!dflag) Death_rate = 0;
 //                        printf("s %lu %lu %lu\n",seed[0],seed[1],seed[2]);
       S += rexp(Birth_rate + Death_rate,seed);
// printf("%lf %lf %lf\n",S,T,Death_rate);
        if (S > T)  break;
        if (aaa > 1000) break;
        /* SIMULATE JUMP TYPE (BIRTH OR DEATH) */
        if (!dflag)
            max = 1.1;
        else if (rep->nKnots-sdegree >= max_knots)
            max = -0.1;
        else
            max = Birth_rate/(Birth_rate + Death_rate);

 //                         printf("kiss %lu %lu %lu\n",seed[0],seed[1],seed[2]);
      if (kiss(seed) < max) {  /* a birth occurs */

           birth(rep,pop,&full_likelihood,seed);

            /* COMPUTE Partial LIKELIHOODS */
            /* first need to calculate the full cond of a new eta given all others */
            
        }
        
        else { /* a death occurs */
            /* remove a knot */
            if (rep->nKnots-2*sdegree > 1) 
                death(death_rate,rep,pop,&full_likelihood,4,seed);
        }
        if (dflag) {
            free(death_rate);
        }
    }
//    printf("aaa = %d\n",aaa);
    if (dflag)
        free(death_rate);
        
    calAx(rep->Veta,rep->V,rep->eta,(const int)rep->dim_V[0],(const int)rep->dim_V[1]);

    calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],rep->P);
    dlm_forward_filter2(rep,dlmStruc);
    calWdelta(rep->Wm,rep->W,rep->m,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
    calWdelta(rep->Wdelta,rep->W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
    

 //   calculate_marginal_residuals(rep,rep->P);
 //   calculate_residuals(rep,rep->P);
    free(Qinv);
    return(aaa);
}


void knot_death_rate(double *death_rate,REP *rep,double *partial_likelihood,double full_likelihood,double Birth_rate,double prior_rate,
                        int sdegree)
{
    int i,nKnots;
    
    nKnots = rep->nKnots-2*sdegree;
    
    if (nKnots > 0) {
        for (i=0;i<nKnots;i++) {
            death_rate[i] = log(Birth_rate/prior_rate) + partial_likelihood[i] - full_likelihood;
//            printf("%g %g %g\n",partial_likelihood[i],full_likelihood,death_rate[i]);
        }
    }
    else {
        if (rep->nKnots==0)
            printf("No knots available\n");
        else {
            death_rate = (double *)calloc(nKnots,sizeof(double));
            death_rate[0] = -1e300;
        }
    }
}

void remove_knot(REP *rep,double *V,double *eta,double *knots,const int sdegree,int remove)
{
    int i,j,k;
    double *data;
    void bsplinebasis(int,int,int,double,double *,int,double *);
    
    k = 0;
    for (i=0;i<rep->nKnots;i++) {
        if (i != remove+sdegree) {
            knots[k] = rep->knots[i];
            k++;
        }
    }
    k = 0;
    for (j=0;j<rep->dim_V[1];j++) {
        if (j != remove+sdegree-2) {
            eta[k] = rep->eta[j];
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

int add_knot(REP *rep,double *V,double *knots,double *eta,int sdegree,double new_eta,double position)
{
    int ii,i,j,k,flag,pos,nKnots;
    double *data;
    void bsplinebasis(int,int,int,double,double *,int,double *);
    
    flag = 1;
    for (i=0;i<rep->nKnots;i++) {
        if (rep->knots[i] < position) {
            knots[i] = rep->knots[i];
//            printf("i = %d %lf %lf\n",i,rep->knots[i],rep->knots[i+1]);fflush(stdout);
        }
        else
            break;
    }
    k = i-2;
//    printf("k = %d i = %d position = %lf\n",k,i,position);
    knots[i] = position;
    for (j=i;j<rep->nKnots;j++)
        knots[j+1] = rep->knots[j];
    
    for (i=0;i<k;i++)
        eta[i] = rep->eta[i];
    eta[k] = 0;//new_eta;
    pos = k;
    k++;
    for (i=k;i<rep->dim_V[1]+1;i++)
        eta[i] = rep->eta[i-1];

    data = (double *)calloc(rep->dim_V[0],sizeof(double));
    for (i=0;i<rep->dim_V[0];i++)
        data[i] = i;
    
    nKnots = rep->nKnots+1;
    for (i=0;i<rep->dim_V[0];i++)
        bsplinebasis(0,0,sdegree,data[i],knots,nKnots,&V[i*(rep->dim_V[1]+1)]);
    
    free(data);
    return pos;
}

