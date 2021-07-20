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

double *pos_prob;
extern int max_knots;
extern double xi;

double *calc_prob_dist_location(double *location,int K) 
{
    double *p,denom=0;
    p = (double *)calloc(K+1,sizeof(double));
    for (int i=0;i<K+1;i++) {
        p[i] = (location[i] - location[i-1]);
//        p[i] = 1./(location[i] - location[i-1]);
    }
    for (int i=1;i<K+1;i++)
        p[i] += p[i-1];
    for (int i=0;i<K+1;i++)
        p[i] /= p[K];
    return(p);
}

double *get_prob_birth(double *location,int K)
{
    double *P,*probs;
    double *newlocation = (double *)calloc(K-1,sizeof(double));
    probs = (double *)calloc(K,sizeof(double));
    
    
//    printf("\nK = %d\n",K);
    for (int j=0;j<K-1;j++) {
        int l=0;
        for (int i=0;i<K-1;i++) {
            if (i == j) 
                l=1;
            newlocation[i] = location[i+l];
        }
        P = calc_prob_dist_location(newlocation,K-2);
//        for (int r=0;r<K-1;r++)
//        printf("%3d %lf \n",r,P[r]);
//    printf("\n");fflush(NULL);exit(0);
        if (j==0)
            probs[j] = P[j];
        else
            probs[j] = P[j]-P[j-1];        
        free(P);
    }
    free(newlocation);
    return(probs);
}

void Ax_Spline(double *Ax,double *A,double *x,const int nrow,const int ncol,const int sdegree,double *knots,const int Nknots) 
{
    int start,end;

    int i = 0;
    for (int j=0;j<Nknots;j++) {
        while (i < knots[j]) {
            if (i > nrow-1)
                break;
            start = (0 > j-sdegree) ? 0:(j-sdegree); 
            end = (ncol < j+1) ? ncol:(j+1);
            Ax[i] = 0;
            for (int k=start;k<end;k++)
                Ax[i] += A[i*ncol+k]*x[k];
            i++; 
        }
    }
}
void draw_knot_locations(REP *rep,int sdegree,unsigned long *seed)
{
    int i,j,k;
    double old_like,new_like;
    double like_ratio,prior_ratio,*knots;
    double *knots2,tmp,*eta,*eta2;
    void bsplinebasis(int j,int k,int SplineOrder,double x,double *knots,int Nknots,double *b);
    void calculate_residualsA(double *,double *,double *,double *,double *,int,int);
    void calAx(double *Ax,double *A,double *x,const int nrow,const int ncol);   
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);    
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calApVinvA(double *ApA,double **A,double *Vinv,const int nrow,const int ncol);
 
    knots = (double *)calloc(rep->nKnots,sizeof(double));
    knots2 = (double *)calloc(rep->nKnots,sizeof(double));
    
    // calc current likelihood
 //   calculate_residuals(rep->residuals,rep->Y,rep->Xbeta,rep->Veta,rep->Wdelta,xi,rep->P,rep->dim_X[0]);

    old_like = 0;
    for (i=0;i<rep->dim_X[0];i++)
        old_like += rep->d_Y[i]*rep->residuals[i]*rep->residuals[i];
   
     
    for (i=0;i<rep->nKnots;i++)
        knots[i] = rep->knots[i];
    
    eta = (double *)calloc(rep->dim_V[1],sizeof(double));
    eta2 = (double *)calloc(rep->dim_V[1],sizeof(double));
    for (i=0;i<rep->dim_V[1];i++)
        eta[i] = rep->eta[i];
    
           // attempt change of a random knot location
    j = runiform_n(rep->nKnots-sdegree*2,seed) + sdegree;    
    
    knots[j] = rnorm(knots[j],rep->prop_sd[3],seed);//runif_atob(seed,knots[j-1],knots[j+1]);  // propose knot location
    for (i = sdegree;i<rep->nKnots-sdegree;i++){
        if (!(i==j)) {
            if (fabs(knots[j]-knots[i]) < 1e-100) {
                free(eta);
                free(eta2);
                free(knots);
                free(knots2);
                return;
            }
        } 
    }
        (rep->attempt[3])++;

    int lower,upper;
/*    int T1 = 9550/10;
    int T2 = 9650/10;
    if (rep->knots[j] < T1) {
        lower = 0;
        upper = T1;
    }
    else if (rep->knots[j] > T2) {
        lower = T2;
        upper = rep->dim_V[0]-1;
    }
    else {
        lower = T1;
        upper = T2;
    }*/
    lower = 0;
    upper = rep->dim_V[0]-1;
    
    if (knots[j] > lower && knots[j] < upper) {  // ensure knot is still interior knot
//        if ((knots[j] > knots[j-1]) && (knots[j] < knots[j+1])) {  // ensure knot is still interior knot
            
        for (i=0;i<rep->nKnots;i++)
            knots2[i] = knots[i];
        for (i=0;i<rep->dim_V[1];i++)
            eta2[i] = eta[i];
            
        for (i=sdegree;i<rep->nKnots-sdegree-1;i++) {  // reorder knots so their locations are sequential
            for (k=i+1;k<rep->nKnots-sdegree;k++) {
                if (knots2[k] < knots2[i]) {
                    tmp = knots2[k];
                    knots2[k] = knots2[i];
                    knots2[i] = tmp; 
                    tmp = eta2[k-2];
                    eta2[k-2] = eta2[i-2];
                    eta2[i-2] = tmp;                     
                }
            }
        }
 //           for (i=0;i<rep->dim_V[0]*rep->dim_V[1];i++)
 //                   rep->V[i] = 0;
        double *V = (double *)calloc(rep->dim_V[0]*rep->dim_V[1],sizeof(double));
        for (i=0;i<rep->dim_V[0]*rep->dim_V[1];i++)
            V[i] = rep->V[i];
                
        for (i=0;i<rep->dim_V[0];i++)   // create new bspline basis
            if ((rep->data[i] >= knots2[j-2]) &&  (rep->data[i] <= knots2[j+2]))
                bsplinebasis(0,0,sdegree,rep->data[i],knots2,rep->nKnots,&V[i*rep->dim_V[1]]);
            

           // calculate Veta, calculate residuals
        double *Veta = (double *)calloc(rep->dim_V[0],sizeof(double));
        double *W = (double *)calloc(rep->dim_W[0]*rep->dim_W[1],sizeof(double));
        double *Wdelta = (double *)calloc(rep->dim_W[0],sizeof(double));
        double *res = (double *)calloc(rep->dim_X[0],sizeof(double));
        Ax_Spline(Veta,V,eta2,(const int)rep->dim_V[0],(const int)rep->dim_V[1],4,knots2,(const int)rep->nKnots);
        calW(W,rep->Y,rep->Xbeta,Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],rep->P);
        calWdelta(Wdelta,W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
 
        calculate_residualsA(res,rep->Y,rep->Xbeta,Veta,Wdelta,rep->P,rep->dim_X[0]);
            // calculate new log likelihood
        new_like = 0;
        for (i=0;i<rep->dim_X[0];i++)
            new_like += rep->d_Y[i]*res[i]*res[i];
           
        like_ratio = -0.5*(new_like-old_like);
           
        prior_ratio = 0;//log(prb_new/prb_old);
            
  //          printf("%lf %lf %lf\n",kiss(seed),like_ratio,prior_ratio);
        rep->alpha_MH[3] = exp(like_ratio + prior_ratio);
        if (rep->alpha_MH[3] > 1) rep->alpha_MH[3] = 1.;
           
        if (log(kiss(seed)) < like_ratio + prior_ratio) { // accept new location
            (rep->accept[3])++;
            free(rep->eta);
            rep->eta = eta2;
            free(rep->knots);
            rep->knots = knots2;
            free(rep->V);
            rep->V = V;
            free(rep->Veta);
            rep->Veta = Veta;
            free(rep->W);
            rep->W = W;
            free(rep->Wdelta);
            rep->Wdelta = Wdelta;
            free(rep->residuals);
            rep->residuals = res;
        }
        else {
            free(V);
            free(eta2);
            free(knots2);
            free(Veta);
            free(W);
            free(Wdelta);
            free(res);
        }
    }
    else {
        rep->alpha_MH[3] = 0;
        free(knots2);
        free(eta2);
    }
    // free allocated scratch memory
    
    free(eta); 
    free(knots);
}

int calculate_death_rates(long double *death_rate,double *Death_rate,double Birth_rate,double prior_rate,REP *rep,int sdegree,long double full_likelihood) {
    int i,j,k;
    long double *partial_likelihood;
    long double max;
    double *knots,*eta,*V;
    double *Veta,*W,*Wdelta,*res;
    long double DR;
    
    void remove_knot(REP *rep,double *V,double *eta,double *knots,const int sdegree,int remove);
    void calAx(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void knot_death_rate(long double *,REP *,long double *,long double,double,double,int);
    void calculate_residualsA(double *,double *,double *,double *,double *,int,int);
    partial_likelihood = (long double *)calloc(rep->nKnots-2*sdegree,sizeof(long double));
    
    if (rep->nKnots-2*sdegree > 0) {
        res = (double *)calloc(rep->dim_X[0],sizeof(double)); 
        knots = (double *)calloc(rep->nKnots-1,sizeof(double));
        eta = (double *)calloc(rep->dim_V[1]-1,sizeof(double));
        V = (double *)calloc(rep->dim_V[0]*(rep->dim_V[1]-1),sizeof(double));
        Veta = (double *)calloc(rep->dim_V[0],sizeof(double));
        W = (double *)calloc(rep->dim_W[0]*rep->dim_W[1],sizeof(double));
        Wdelta = (double *)calloc(rep->dim_W[0],sizeof(double)); 
//        for (i=0;i<rep->dim_V[0];i++)
//            V[i] = (double *)calloc(rep->dim_V[1]-1,sizeof(double));
        
       for (i=0;i<rep->nKnots-2*sdegree;i++) {
            
            for (j=0;j<rep->dim_V[0]*(rep->dim_V[1]-1);j++)
                V[j] = 0;
            for (j=0;j<rep->dim_V[0];j++) {
                int kk = 0;
                for (k=0;k<rep->dim_V[1];k++) {
                    if (k != i+sdegree-2) {
//                        if (k >= rep->dim_V[1]) {
//                            printf("276 %d %d %d\n",k,i,rep->dim_V[1]);fflush(NULL);}
                        V[j*(rep->dim_V[1]-1)+kk] = rep->V[j*rep->dim_V[1]+k];
                        kk++;
                    }
                }  
            }      
           
               /* remove each knot in turn and calculate partial likelihoods */
//         printf("AAA\n");fflush(NULL);
             remove_knot(rep,V,eta,knots,4,i);
//         printf("BBB\n");fflush(NULL);
             /* calculate Veta */
            
            Ax_Spline(Veta,V,eta,(const int)rep->dim_V[0],(const int)(rep->dim_V[1]-1),4,knots,(const int)(rep->nKnots-1));
//            calAx(Veta,V,eta,(const int)rep->dim_V[0],(const int)(rep->dim_V[1]-1));
            
            int P = rep->dim_W[1];
            /* update W and Wdelta */
            calW(W,rep->Y,rep->Xbeta,Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],P);
            calWdelta(Wdelta,W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
            
            calculate_residualsA(res,rep->Y,rep->Xbeta,Veta,Wdelta,P,rep->dim_X[0]);

            partial_likelihood[i] = (long double)rep->d_Y[0]*(long double)res[0]*(long double)res[0];
            for (j=1;j<rep->dim_X[0];j++)
                partial_likelihood[i] += (long double)rep->d_Y[j]*(long double)res[j]*(long double)res[j];
            partial_likelihood[i] *= -0.5;
        }

        free(knots);
        free(eta);
        free(V);
        free(Veta);
        free(W);
        free(Wdelta);
        free(res);
    }
    
    /* CALCULATE DEATH RATE FOR EACH COMPONENT */
    int dflag = 0;

    if (rep->nKnots-2*sdegree > 0) {
        dflag = 1;
        knot_death_rate(death_rate,rep,partial_likelihood,full_likelihood,Birth_rate,prior_rate,sdegree);
    }
    if (dflag) {
        DR = death_rate[0];
        for (i=1;i<rep->nKnots-2*sdegree;i++) {
            max = (DR > death_rate[i]) ? DR:death_rate[i];
            if (DR == max)
                DR = log1pl(expl(death_rate[i]-max)) + max;
            else
                DR = log1pl(expl(DR-max)) + max;
            death_rate[i] = DR;
           /* the above computes    *Death_rate = log(exp(*Death_rate-max) + exp(death_rate[i]-max)) + max;*/
        }
       
        for (i=0;i<rep->nKnots-2*sdegree;i++) 
            death_rate[i] = expl(death_rate[i]-DR);
        *Death_rate = (double)expl(DR);
    }
    else
        *Death_rate = 0;

    free(partial_likelihood);

    return dflag;
}

void birth(REP *rep,long double *full_likelihood,int iter,unsigned long *seed) {
    int i,j,pos_idx,sdegree=4;
    double position,new_eta;
    double *teta,*tknots,*VV;
    
    int add_knot(REP *rep,double *V,double *knots,double *eta,int sdegree,double new_eta,double position);
    void calculate_residualsA(double *,double *,double *,double *,double *,int,int);
    void calAx(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);

/*    double *pos_prob = calc_prob_dist_location(&rep->knots[sdegree],rep->nKnots-2*sdegree);
    if (rep->nKnots-2*sdegree) 
        j = (int)rmultinomial(pos_prob,rep->nKnots-2*sdegree,seed);  
    else {
        j = 0;
    }
    position = runif_atob(seed,rep->knots[sdegree-1+j],rep->knots[sdegree+j]);  
*/
/*    int T1 = 9550/10;
    int T2 = 9650/10;
    double region;
    if (kiss(seed) < 0.5)
        position = runif_atob(seed,T1,T2);
    else {
        if (kiss(seed) < 0.5)
            position = runif_atob(seed,0,T1-1);
        else 
            position = runif_atob(seed,T2+1,rep->dim_V[0]-1);
    }  */      
    position = runif_atob(seed,0,rep->dim_V[0]-1);
    for (i = sdegree;i<rep->nKnots-sdegree;i++){
        if (fabs(position-rep->knots[i]) < 1e-100) {
/*            if (rep->nKnots-2*sdegree) 
                j = (int)rmultinomial(pos_prob,rep->nKnots-2*sdegree,seed);  
            else
                j = 0;
            position = runif_atob(seed,rep->knots[sdegree-1+j],rep->knots[sdegree+j]);  
*/
            position = runif_atob(seed,0,rep->dim_V[0]-1);
            i = sdegree;
          //  return;
        }
    }

//    free(pos_prob);
     /* simulate eta from its full cond */

    for (pos_idx=0;pos_idx<rep->nKnots;pos_idx++) 
        if (position < rep->knots[pos_idx]) 
            break;
        
    pos_idx = pos_idx - 2;
//printf("pos_idx = %d ",pos_idx);     
    VV = (double *)calloc(rep->dim_V[0]*(rep->dim_V[1]+1),sizeof(double));
    for (int j=0;j<rep->dim_V[0];j++) {
        int kk = 0;
        for (int k=0;k<rep->dim_V[1]+1;k++) {
            if (k != pos_idx) {
 //               if (kk >= rep->dim_V[1]) {
 //                   printf("405 %d %d \n",kk,rep->dim_V[1]);fflush(NULL);}
                VV[j*(rep->dim_V[1]+1)+k] = rep->V[j*rep->dim_V[1]+kk];
                kk++;
            }
        }  
    }      
    
    teta = (double *)calloc(rep->dim_V[1]+1,sizeof(double));
    tknots = (double *)calloc(rep->nKnots+1,sizeof(double));
    
    pos_idx = add_knot(rep,VV,tknots,teta,4,new_eta,position);
    teta[pos_idx] = rnorm(0,1./sqrt(rep->preceta),seed);
  
    (rep->dim_V[1])++;
    (rep->nKnots)++;

//    rep->eta = (double *)realloc(rep->eta,sizeof(double)*(rep->dim_V[1]));
    free(rep->eta);
    rep->eta = teta;
 /*   for (i=0;i<rep->dim_V[1];i++)
        rep->eta[i] = teta[i];
    free(teta);*/
    
    free(rep->knots);
    rep->knots = tknots;
/*    rep->knots = (double *)realloc(rep->knots,sizeof(double)*(rep->nKnots));
    for (i=0;i<rep->nKnots;i++)
        rep->knots[i] = tknots[i];
    free(tknots);*/
    
    free(rep->V);
    rep->V = VV;
    
/*    rep->V = (double *)realloc(rep->V,rep->dim_V[0]*rep->dim_V[1]*sizeof(double));
    
    for (i=0;i<rep->dim_V[0];i++) {
        for (j=0;j<rep->dim_V[1];j++) {
            rep->V[i*rep->dim_V[1]+j] = VV[i*rep->dim_V[1]+j];
        }
    }
    
    free(VV);*/
    
    /* calculate Veta */
    
    Ax_Spline(rep->Veta,rep->V,rep->eta,(const int)rep->dim_V[0],(const int)rep->dim_V[1],4,rep->knots,(const int)rep->nKnots);
//    calAx(rep->Veta,rep->V,rep->eta,(const int)rep->dim_V[0],(const int)rep->dim_V[1]);
    
    int P = rep->dim_W[1];
    /* update W and Wdelta */
    calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],P);
    calWdelta(rep->Wdelta,rep->W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
    
    calculate_residualsA(rep->residuals,rep->Y,rep->Xbeta,rep->Veta,rep->Wdelta,rep->P,rep->dim_X[0]);
    
    
    *full_likelihood = (long double)rep->d_Y[0]*(long double)rep->residuals[0]*(long double)rep->residuals[0];
    for (i=1;i<rep->dim_X[0];i++)
        *full_likelihood += (long double)rep->d_Y[i]*(long double)rep->residuals[i]*(long double)rep->residuals[i];
    *full_likelihood *= (long double)-0.5;
    
    
}

void death(long double *death_rate,REP *rep,long double *full_likelihood,int sdegree,unsigned long *seed) {
    int i,remove;
    double *eta,*tknots;
    
    void remove_knot(REP *rep,double *V,double *eta,double *knots,const int sdegree,int remove);
    void calculate_residualsA(double *,double *,double *,double *,double *,int,int);
    void calAx(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);

    remove = (int)rmultinomial2(death_rate,(long)(rep->nKnots-2*sdegree),seed);
//    for (i=0;i<rep->dim_V[0];i++)
//        free(rep->V[i]);
//    free(rep->V);
 
//        rep->V = (double *)realloc(rep->V,rep->dim_V[0]*(rep->dim_V[1]-1)*sizeof(double));
       double *V = (double *)calloc(rep->dim_V[0]*(rep->dim_V[1]-1),sizeof(double));

        
//    for (i=0;i<rep->dim_V[0];i++)
//        free(rep->J[i]);
//    free(rep->J);
    
//    rep->J = (double **)calloc(rep->dim_V[0],sizeof(double *));
//    for (i=0;i<rep->dim_V[0];i++)
//        rep->J[i] = (double *)calloc(rep->dim_V[1]-1,sizeof(double));
    
    for (int j=0;j<rep->dim_V[0];j++) {
        int kk = 0;
        for (int k=0;k<rep->dim_V[1];k++) {
            if (k != remove+sdegree-2) {
                V[j*(rep->dim_V[1]-1)+kk] = rep->V[j*rep->dim_V[1]+k];
                kk++;
            }
        }  
    }      
    eta = (double *)calloc(rep->dim_V[1]-1,sizeof(double));
    
    tknots = (double *)calloc(rep->nKnots-1,sizeof(double));

    remove_knot(rep,V,eta,tknots,4,remove);
    free(rep->V);
    rep->V = V;
    (rep->nKnots)--;
    (rep->dim_V[1])--;
    
    free(rep->eta);
    rep->eta = eta;
    
    free(rep->knots);
    rep->knots = tknots;
    
    /* calculate Veta */
   
    Ax_Spline(rep->Veta,rep->V,rep->eta,(const int)rep->dim_V[0],(const int)rep->dim_V[1],4,rep->knots,(const int)rep->nKnots);
//    calAx(rep->Veta,rep->V,rep->eta,(const int)rep->dim_V[0],(const int)rep->dim_V[1]);
    
    int P = rep->dim_W[1];
    /* update W and Wdelta */
    calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],P);
    calWdelta(rep->Wdelta,rep->W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
    
    
    calculate_residualsA(rep->residuals,rep->Y,rep->Xbeta,rep->Veta,rep->Wdelta,rep->P,rep->dim_X[0]);
    
    *full_likelihood = (long double)rep->d_Y[0]*(long double)rep->residuals[0]*(long double)rep->residuals[0];
    for (i=1;i<rep->dim_X[0];i++)
        *full_likelihood += (long double)rep->d_Y[i]*(long double)rep->residuals[i]*(long double)rep->residuals[i];
    *full_likelihood *= -0.5;
}

int knot_birth_death(REP *rep,int pknots,const int sdegree,int iter,unsigned long *seed){
    
    int i,aaa,dflag;
    double Birth_rate,T,max,S;
    double Death_rate,prior_rate;
    long double *death_rate,full_likelihood;

    aaa = 0;
    
    // set simulation time

//    prior_rate = rep->nKnots-2*sdegree;  // gives a uniform prior on the number of knots (1 - MAX NUMBER)
                                         // not the prior rate, but gives correct death rate in knot_death_rate
//    prior_rate = rgamma(pknots,1,seed);  // gives a negative binomial prior on the number of knots
    prior_rate = pknots;
    Birth_rate = 1.;
    T = 1.;

    // calculate likelihood;
    full_likelihood = (long double)rep->d_Y[0]*(long double)rep->residuals[0]*(long double)rep->residuals[0];
    for (i=1;i<rep->dim_X[0];i++)
        full_likelihood += (long double)rep->d_Y[i]*(long double)rep->residuals[i]*(long double)rep->residuals[i];
    full_likelihood *= -0.5;
    S = 0.0;

    while (1) {
        aaa++;

        if (rep->nKnots-2*sdegree > 1) {
            death_rate = (long double *)calloc(rep->nKnots-2*sdegree,sizeof(long double));
            
            dflag = calculate_death_rates(death_rate,&Death_rate,Birth_rate,prior_rate,rep,sdegree,full_likelihood);
        }
        else {
            Death_rate = 0;
            dflag = 0;
        }
         /* SIMULATE TIME, S, TO NEXT JUMP */

        S += rexp(Birth_rate + Death_rate,seed);
        if (S > T)  break;
        if (aaa > 1000) break;
        /* SIMULATE JUMP TYPE (BIRTH OR DEATH) */
        
        max = Birth_rate/(Birth_rate + Death_rate);
        if (kiss(seed) < max) {  /* a birth occurs */
            /* add a knot */
            if (rep->nKnots-2*sdegree < max_knots)
                birth(rep,&full_likelihood,iter,seed);
       }       
        else { /* a death occurs */
            /* remove a knot */
            death(death_rate,rep,&full_likelihood,sdegree,seed);
        }
        if (dflag) {
            free(death_rate);
        }
//        free(pos_prob);
    }
//    printf("aaa = %d knots = %d\n",aaa,rep->nKnots);
    if (dflag)
        free(death_rate);
//  free(pos_prob);
    return(aaa);
}


void knot_death_rate(long double *death_rate,REP *rep,long double *partial_likelihood,long double full_likelihood,double Birth_rate,double prior_rate,int sdegree)
{
    int i,nKnots;
    long double pi = 3.141592653589793116; 
    double beta =  1.;
    double alpha = beta*prior_rate;
    double N,min = 1e100;
    double *birth_probs,birth_prob;
    nKnots = rep->nKnots-2*sdegree; // # knots
    N = rep->knots[nKnots+sdegree]; // length of TS
    long double twopi = 2.*pi;
    long double lratio;
    
 //   long double PR = (long double)(prior_rate);
    long double PR = (long double)(nKnots+alpha-1)/(long double)(beta+1);
    if (nKnots > 0) {
        lratio = logl((long double)Birth_rate/PR) - full_likelihood ;
//        lratio = logl((long double)Birth_rate/(long double)prior_rate) - full_likelihood ;
        for (i=0;i<nKnots;i++) {
           death_rate[i] =  lratio + partial_likelihood[i];// + 0.5*logl((long double)rep->preceta/1e-6) - 0.5*((long double)rep->preceta-1e-6)*(long double)rep->eta[i+2]*(long double)rep->eta[i+2];
        }
    }
    else {
        if (rep->nKnots==0)
            printf("No knots available\n");
        else {
            death_rate = (long double *)calloc(nKnots,sizeof(double));
            death_rate[0] = -1e300;
        }
    }
}

void remove_knot(REP *rep,double *V,double *eta,double *knots,const int sdegree,int remove)
{
    int i,j,k;
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

    for (i=0;i<rep->dim_V[0];i++) {
        if ((rep->data[i] >= rep->knots[remove+sdegree-2]) && (rep->data[i] <= rep->knots[remove+sdegree+2]))
            bsplinebasis(0,0,sdegree,rep->data[i],knots,rep->nKnots-1,&V[i*(rep->dim_V[1]-1)]);
    }
}

int add_knot(REP *rep,double *V,double *knots,double *eta,int sdegree,double new_eta,double position)
{
    int i,j,k,pos,nKnots;
    void bsplinebasis(int,int,int,double,double *,int,double *);
    
    for (i=0;i<rep->nKnots;i++) {
        if (rep->knots[i] < position) {
            knots[i] = rep->knots[i];
        }
        else
            break;
    }
    k = i-2;
    knots[i] = position;
    int pos2 = i;
    for (j=i;j<rep->nKnots;j++)
        knots[j+1] = rep->knots[j];
    
    for (i=0;i<k;i++)
        eta[i] = rep->eta[i];
    eta[k] = 0;
    pos = k;
    k++;
    for (i=k;i<rep->dim_V[1]+1;i++)
        eta[i] = rep->eta[i-1];

     nKnots = rep->nKnots+1;
    
    for (i=0;i<rep->dim_V[0];i++)
        if ((rep->data[i] >= knots[pos2-2]) && (rep->data[i] <= knots[pos2+2])) 
            bsplinebasis(0,0,sdegree,rep->data[i],knots,nKnots,&V[i*(rep->dim_V[1]+1)]);
    
    return pos;
}

