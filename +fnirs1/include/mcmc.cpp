/*
 *  mcmc.c
 *  LGCP
 *
 *  Created by Timothy D. Johnson on 11/10/15.
 *  Copyright 2015 University of Michigan. All rights reserved.
 *
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "cholesky.h"
#include "fNIRS.h"
// #include "printToMatlab.h"
#include "randgen.h"


extern int MAX_ITER;
extern int BURN_IN;
extern sDLM *dlmStruc;
extern int maxP;
extern FILE *flog;
double md1,md2;
int mP;

void adjust_acceptance(double x,double *X,double rate)
{
	double y;
	
	y = 1. + 2.*(x-rate)*(x-rate)*(x-rate);
	if (y < .8)
		y = .8;
	if (y > 1.2)
		y = 1.2;
	*X *= y;
}

void adjust_acceptance2(double x,double *X,double rate)
{
    double y;
    
    y = 1. + 10.*(x-rate)*(x-rate)*(x-rate);
    if (y < .8)
        y = .8;
    if (y > 1.2)
        y = 1.2;
    *X *= y;
}


void mcmc(POP *pop,unsigned long *seed) {
    int i,j,k,iter,isub,irep;
    double **loglik;
    int *tableP;
    int nrow_loglik;
    SUB *sub;
    REP *rep;
    FILE *fout;
    
    void draw_knot_locations(REP *rep,int sdegree,int *flag,unsigned long *seed);
    int knot_birth_death(REP *rep,POP *pop,const int sdegree,int iter,unsigned long *seed);
 
//    void mcmc_delta1(REP *rep,int iter,dlm_tst_flag flag,unsigned long *seed);
//    void mcmc_delta2(REP *rep,int iter,dlm_tst_flag flag,unsigned long *seed);
    void DLMtst(REP *rep,int iter,dlm_tst_flag flag,unsigned long *seed);
    void DLM(REP *rep,unsigned long *seed);
    void draw_precYstart(REP *rep,unsigned long *seed);
    void draw_preceta(REP *rep,unsigned long *seed);
    void draw_reprec(POP *pop,SUB *sub,unsigned long *seed);
 
/*    void calculate_residuals(REP *rep,int P);
    void calculate_res1(REP *rep,int P);
    void calculate_res2(REP *rep,int P);
    void calculate_res3(REP *rep);
    void calculate_res4(REP *rep,int P);*/
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    
    void draw_beta_eta(POP *pop,SUB *sub,REP *rep,unsigned long *seed);
//    void draw_beta_eta(POP *pop,SUB *sub,REP *rep,unsigned long *seed,int iter,int isub);
    void draw_pop_beta(POP *pop,SUB *sub,unsigned long *seed);
    
    void DIC(POP *pop,REP *rep,unsigned long *seed);
   
    pop->ED = 0;
       
    int sdegree = 4;
    int aaa = 0;

    for (iter=0;iter<=MAX_ITER;iter++) {
        if (!(iter%100)) { 
            printf("%d",iter);
            fflush(stdout);
        }
        if (!(iter%20) && (iter%100)) {
            printf(".");
            fflush(stdout);
        }
        
        for (isub=0;isub<pop->N_SUBS;isub++) {
            sub = &(pop->sub[isub]);
            
            for (irep=0;irep<sub->N_REPS;irep++) {
                rep = &(sub->rep[irep]);
                
               if (iter > 0)
                   aaa = knot_birth_death(rep,pop,sdegree,iter,seed);
 
                fprintf(rep->fout_nknots,"%d ",rep->nKnots-2*sdegree);fflush(rep->fout_nknots);
                
                 for (i=0;i<rep->nKnots-sdegree*2;i++) {
                    int flag = 0;
                    draw_knot_locations(rep,sdegree,&flag,seed);
                }
 
                if (rep->dim_W[1] > 0) {
                    DLMtst(rep,iter,fdelta1,seed);
                    DLMtst(rep,iter,fdelta2,seed);
                    DLMtst(rep,iter,fP,seed);
                    DLM(rep,seed);
                }
                draw_beta_eta(pop,sub,rep,seed);
                
                draw_preceta(rep,seed);
                draw_precYstart(rep,seed);
            }   
        }
     
        if (pop->GRP) {
            draw_reprec(pop,pop->sub,seed);
            draw_pop_beta(pop,pop->sub,seed);
        }
        
        if (!(iter%50)) {
            for (isub=0;isub<pop->N_SUBS;isub++) {
                sub = &(pop->sub[isub]);
                
                for (irep=0;irep<sub->N_REPS;irep++) {
                    rep = &(sub->rep[irep]);
 
                    if ((iter <= BURN_IN) && (iter > 1)) {
                        double rt;

                        rt = (double)rep->accept[0]/(double)rep->attempt[0];
                        adjust_acceptance2(rt,&(rep->prop_sd[0]),0.35);
                        rep->accept[0] = rep->attempt[0] = 0;
                        if (!(iter%100))
                            fprintf(flog,"\t\t rt = %lf proposal = %lf\n",rt,rep->prop_sd[0]);

                        rt = (double)rep->accept[1]/(double)rep->attempt[1];
                        adjust_acceptance2(rt,&(rep->prop_sd[1]),0.35);
                        pop->sub[isub].rep[irep].accept[1] = pop->sub[isub].rep[irep].attempt[1] = 0;
                        if (!(iter%100))
                            fprintf(flog,"\t\t rt = %lf proposal = %lf\n",rt,pop->sub[isub].rep[irep].prop_sd[1]);
   
                        rt = (double)rep->accept[3]/(double)rep->attempt[3];
                        adjust_acceptance2(rt,&(rep->prop_sd[3]),0.35);
                        rep->accept[3] = rep->attempt[3] = 0;
                        if (!(iter%100))
                            fprintf(flog,"\t\t rt = %lf proposal = %lf\n",rt,rep->prop_sd[3]);
                    }
                }
            }
        }
        
        if (!(iter%100)) {
            for (isub=0;isub<pop->N_SUBS;isub++) {
                sub = &(pop->sub[isub]);
                
                for (irep=0;irep<sub->N_REPS;irep++) {
                    rep = &(sub->rep[irep]);
                    
                    fprintf(flog,"iter = %6d\t Sub %d, Rep %d \t int knots = %d",iter,isub,irep,rep->nKnots-2*sdegree);
                    fprintf(flog,"\t precY0 = %10.6lf \t preceta = %10.6lf\n",rep->d_Y[rep->dim_X[0]-1],rep->preceta);
                    fprintf(flog,"\t df_delta1 = %10.6lf \t df_delta2 = %10.6lf \t P = %3d\n",rep->df_delta1,rep->df_delta2,rep->P);
                    
                }
            }
 
            for (k=0;k<pop->Ns;k++) {
                fprintf(flog,"Stimulus %d\n",k);
                if (pop->GRP) {
                    fprintf(flog,"\t %g ",pop->re_prec[k]);
                    fprintf(flog,"\n\n");
                }
                for (i=0;i<pop->N_SUBS;i++) {
                    sub = &(pop->sub[i]);
                    for (j=0;j<pop->Nb;j++) {
                        fprintf(flog,"\t %10.6lf \t ",sub->beta[k*pop->Nb+j]);
                    }
                    fprintf(flog,"\n");
                }
                fprintf(flog,"\n");
                if (pop->GRP) {
                    for (i=0;i<pop->Ncov*pop->Nb;i++)
                        fprintf(flog,"\t %10.6lf \t ",pop->beta[k*pop->Ncov*pop->Nb+i]);
                    fprintf(flog,"\n\n");
                }
            }
        }
        fflush(flog);
        
        if ((iter>BURN_IN)) {
            for (isub=0;isub<pop->N_SUBS;isub++) {
                sub = &(pop->sub[isub]);
                
                for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
                    rep = &(sub->rep[irep]);

                    for (i=rep->P;i<rep->dim_W[0];i++) {
                        for (j=0;j<rep->dim_W[1];j++) {
                            rep->mdelta[i*rep->dim_W[1]+j] += rep->delta[i*rep->dim_W[1]+j];
                            rep->mdelta2[i*rep->dim_W[1]+j] += rep->delta[i*rep->dim_W[1]+j]*rep->delta[i*rep->dim_W[1]+j];
                        }
                    }
 
                    for (i=0;i<rep->dim_X[0];i++) {
                        rep->md_Y[i] += 1./sqrt(rep->d_Y[i]);
                        rep->sd_Y[i] += 1./rep->d_Y[i];
                    }
                    
                    for (i=0;i<rep->dim_V[1];i++)
                        fprintf(rep->fout_eta,"%lf ",rep->eta[i]);
                    fprintf(rep->fout_eta,"\n");
                   
                    DIC(pop,rep,seed);
 
                    for (i=0;i<rep->dim_X[0];i++) {
                        rep->mean_fit[i] += rep->Y[i] - rep->residuals[i];
                        rep->std_res[i]  += rep->residuals[i]*sqrt(rep->d_Y[i]);
                        rep->mean_res[i] += rep->residuals[i];
                        rep->mean_d_Y[i] += rep->d_Y[i];
                    }
                                       
                     
                    for (i=0;i<rep->dim_X[0];i++)
                        rep->mXbeta[i] += rep->Xbeta[i];
                    for (i=0;i<rep->dim_V[0];i++)
                        rep->mVeta[i] += rep->Veta[i];
                    for (i=0;i<rep->dim_W[0];i++)
                        rep->mWdelta[i] += rep->Wdelta[i];
                    
                    fprintf(rep->fout_dlm,"%lf %lf %d\n",rep->df_delta1,rep->df_delta2,rep->P);
                    
//                    fprintf(rep->fout_nknots,"%d ",rep->nKnots-2*sdegree);

                    for (i=sdegree;i<rep->nKnots-sdegree;i++)
                        fprintf(rep->fout_knots,"%lf ",rep->knots[i]);
                    fprintf(rep->fout_knots,"\n");  
                }
                
                for (i=0;i<pop->Nb*(pop->Ns);i++)
                    fprintf(sub->fout_beta,"%lf ",sub->beta[i]);
                fprintf(sub->fout_beta,"\n");
            }

            for (i=0;i<pop->Ncov*pop->Nb*pop->Ns;i++)
                fprintf(pop->fout_beta,"%lf ",pop->beta[i]);
            fprintf(pop->fout_beta,"\n");
            for (i=0;i<pop->Ns;i++)
                fprintf(pop->fout_reprec,"%lf ",pop->re_prec[i]);
            fprintf(pop->fout_reprec,"\n");
            fflush(NULL);      
        }
    }
    
    fclose(pop->fout_beta);
    fclose(pop->fout_reprec);

    for (isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
            fclose(rep->fout_dlm);
            fclose(rep->fout_eta);
            fclose(rep->fout_nknots);
            fclose(rep->fout_knots);
        }
        fclose(sub->fout_beta);
    }

    for (isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
            for (i=rep->P;i<rep->dim_W[0];i++) {
                for (j=0;j<rep->dim_W[1];j++) {
                    rep->mdelta[i*rep->dim_W[1]+j] /= (double)(MAX_ITER-BURN_IN);
                    fprintf(rep->fout_delta,"%lf ",rep->mdelta[i*rep->dim_W[1]+j]);
                }
                fprintf(rep->fout_delta,"\n");                  
            }                
            fclose(rep->fout_delta);
        }
    }
    
   for (isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
            for (i=0;i<rep->dim_X[0];i++) {
                    rep->md_Y[i] /= (double)(MAX_ITER-BURN_IN);
                    fprintf(rep->fout_stdev,"%lf ",rep->md_Y[i]);
            }                
            fprintf(rep->fout_stdev,"\n");                  
            fclose(rep->fout_stdev);
        }
    }     
      
    for (isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
            for (i=0;i<rep->dim_W[0];i++)
                fprintf(rep->fout_wdelta,"%lf ",rep->mWdelta[i]/(double)((MAX_ITER-BURN_IN)));
            fclose(rep->fout_wdelta);
        }
    }
    
    for (isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
            for (i=0;i<rep->dim_V[0];i++)
                fprintf(rep->fout_veta,"%lf ",rep->mVeta[i]/(double)((MAX_ITER-BURN_IN)));
            fclose(rep->fout_veta);
        }
    }

    for (isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
            for (i=0;i<rep->dim_X[0];i++)
                fprintf(rep->fout_Xbeta,"%lf ",rep->mXbeta[i]/(double)((MAX_ITER-BURN_IN)));
            fclose(rep->fout_Xbeta);
        }
    }

    pop->ED /= ((MAX_ITER-BURN_IN));
    for (isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
            for (i=0;i<rep->dim_X[0];i++) {
                rep->std_res[i] /= ((MAX_ITER-BURN_IN));
                rep->mean_res[i] /= ((MAX_ITER-BURN_IN));
                rep->mean_fit[i] /= ((MAX_ITER-BURN_IN));
                rep->mean_d_Y[i] /= ((MAX_ITER-BURN_IN));
            }
       }
    }
    double DE=0;
    double mlogsqrt2pi = -0.9189385332046727;
    for (isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
            for (i=maxP;i<rep->dim_X[0];i++)
                DE += mlogsqrt2pi + sqrt(rep->mean_d_Y[i]) - 0.5*(rep->mean_res[i])*(rep->mean_res[i]*rep->mean_d_Y[i]);
        }
    }
    DE *= -2;
    fprintf(flog,"DIC = %lf pD = %lf\n",2*pop->ED - DE,pop->ED - DE);
    fout = fopen("./log/DIC.log","w");
    fprintf(fout, "DIC = %15.3lf, pD = %15.3lf\n",2*pop->ED - DE,pop->ED - DE);
    fclose(fout);
    
    for (isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
            for (i=0;i<rep->dim_X[0];i++)
                fprintf(rep->fout_res,"%lf ",rep->mean_res[i]);
            fclose(rep->fout_res);
        }
    }

    for (isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
            for (i=0;i<rep->dim_X[0];i++)
                fprintf(rep->fout_stdres,"%lf ",rep->std_res[i]);
            fclose(rep->fout_stdres);
        }
    }

    for (isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
           for (i=0;i<rep->dim_X[0];i++)
                fprintf(rep->fout_fit,"%lf ",rep->mean_fit[i]);
            fclose(rep->fout_fit);
        }
    }
    fflush(NULL);   
}

void draw_precYstart(REP *rep,unsigned long *seed) {
    double ALPHA=0.001,BETA=0.001;
    double a=0, b=0;
    
    a = (double)rep->P/2. + ALPHA;
    
    b = 0;
    for (int i=0;i<rep->P;i++)
        b += (rep->Y[i]-rep->Veta[i])*(rep->Y[i]-rep->Veta[i]);
    
    rep->precYstart = rgamma(a,b,seed);
    for (int i=0;i<rep->P;i++)
        rep->d_Y[i] = rep->precYstart;

}

void calH(double *H,double *A,double *m,const int nrow,const int ncol,const int P) {
    int t,j,k;
    double *tmp;
    
    for (j=0;j<nrow;j++) 
        for (k=0;k<ncol;k++)
            H[j*ncol+k] = 0;
    
    for (t=P;t<nrow;t++) 
        for (j=0;j<P;j++) 
            for (k=0;k<ncol;k++) 
                H[t*ncol+k] += A[(t-j-1)*ncol+k]*m[t*P+j];
}

void calApVinvX(double *ApX,double *A,double *X,double *Vinv,const int nrow,const int ncol) {
    int i,j;
    
    for (i=0;i<ncol;i++)
        ApX[i] = 0;
        
    for (j=0;j<nrow;j++) 
        for (i=0;i<ncol;i++)
            ApX[i] += A[j*ncol+i]*Vinv[j]*X[j];  
}

void calApVinvA(double *ApA,double *A,double *Vinv,const int nrow,const int ncol) {
// only for diagonal matrices, V
    int i,j,k;
        
    for (i=0;i<ncol;i++)
        for (j=0;j<ncol;j++)
            ApA[i*ncol+j] = 0;
   
     for (k=0;k<nrow;k++) {
        double V = Vinv[k];
        double *tmp = &A[k*ncol];
        for (i=0;i<ncol;i++) {
            double *tmpout = &ApA[i*ncol];
            for (j=0;j<ncol;j++)
                tmpout[j] += V*tmp[i]*tmp[j];
        }
    }
}

void calAx(double *Ax,double *A,double *x,const int nrow,const int ncol) {
    int i,j;
    
    for (i=0;i<nrow;i++) {
        Ax[i] = 0;
        for (j=0;j<ncol;j++)
            Ax[i] += A[i*ncol+j]*x[j];
    }
    
}

void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol) {
    int i,j;
    
     for (i=0;i<nrow;i++) {
        Ax[i] = 0;
        if (i>= ncol) {
            for (j=0;j<ncol;j++)
                Ax[i] += A[i*ncol+j]*x[i*ncol+j];
            }
    }
    
}

void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P) {
    int t,i,j;
        
    for (i=0;i<nrow;i++) {
        int min = (i < ncol) ? i:ncol;
        for (j=0;j<min;j++) {
            W[i*ncol+j] = Y[i-j-1] - Xb[i-j-1] - Ve[i-j-1];
       }
    }
}

void calW2(double *W,double *Y,const int nrow,const int ncol,const int P) {
    int t,i,j;
        
    for (i=0;i<nrow;i++) {
        int min = (i < ncol) ? i:ncol;
        for (j=0;j<min;j++) {
            W[i*ncol+j] = Y[i-j-1];
       }
    }
}

void calculate_residuals(REP *rep,int P) {
    int i,j;
    
    double tmp = 0;
    for (i=0;i<P;i++)
        rep->residuals[i] = rep->Y[i] - rep->Veta[i] - rep->Xbeta[i];
    for (i=P;i<rep->dim_X[0];i++)
        rep->residuals[i] = rep->Y[i] - (rep->Veta[i] + rep->Wdelta[i] + rep->Xbeta[i]);
    
}

void calculate_marginal_residuals(REP *rep,int P) {
    int i,j;
    
    double tmp = 0;
    for (i=0;i<P;i++)
        rep->residuals[i] = rep->Y[i] - rep->Veta[i] - rep->Xbeta[i];
    for (i=P;i<rep->dim_X[0];i++)
        rep->residuals[i] = rep->Y[i] - (rep->Veta[i] + rep->Wm[i] + rep->Xbeta[i]);
    
}

/*void calculate_res1(REP *rep,int P) {
    int i,j;
    
    for (i=0;i<P;i++)
        rep->residuals1[i] = 0;
    for (i=P;i<rep->dim_X[0];i++) {
        rep->residuals1[i] = rep->Y[i] - rep->Wdelta[i] - rep->Veta[i];
    }
    
}

void calculate_res2(REP *rep,int P) {

    for (int i=0;i<P;i++)
        rep->residuals2[i] = 0;
    for (int i=P;i<rep->dim_X[0];i++) {
        rep->residuals2[i] = rep->Y[i] - rep->Wdelta[i] - rep->Xbeta[i] - rep->Veta[i];
    }
}
*/
void calculate_res3(REP *rep) {
    for (int i=0;i<rep->dim_X[0];i++) {
        rep->residuals3[i] = rep->Y[i] -  rep->Xbeta[i] - rep->Veta[i];
    }
}
/*
void calculate_res4(REP *rep,int P) {
    for (int i=0;i<P;i++) {
        rep->residuals4[i] = rep->Y[i];
    }
   for (int i=P;i<rep->dim_V[0];i++)
        rep->residuals4[i] = rep->Y[i] - rep->Xbeta[i] - rep->Wdelta[i];
}
*/
void calculate_res5(REP *rep,int P) {
    for (int i=0;i<P;i++) {
        rep->residuals5[i] = rep->Y[i];
    }
    for (int i=P;i<rep->dim_V[0];i++)
        rep->residuals5[i] = rep->Y[i] - rep->Wdelta[i];
}

void draw_preceta(REP *rep,unsigned long *seed) {
    double ALPHA=0, BETA=0;
    double a, b;
    ALPHA = 1.0;
    BETA =  1.0;
    b = 0;
    
    for (int i=0;i<rep->dim_V[1];i++)
            b += rep->eta[i]*rep->eta[i];
    
    a = 0.5*rep->dim_V[1] + ALPHA;
    
    b = 0.5*b + BETA;
    
    double tmp = rgamma(a,b,seed);
    if (tmp > 0.000001)
        rep->preceta = tmp;
}

void draw_reprec(POP *pop,SUB *sub,unsigned long *seed) {
    int N,M;
    double ALPHA, BETA,tmp;
    double *Xb;

    ALPHA = 1.1;//10;
    BETA  = 1.1;//10;
        
    Xb = (double *)calloc(pop->Ns*pop->Nb,sizeof(double));
    for (int is=0;is<pop->Ns;is++) {
        tmp = 0;
        for (int i=0;i<pop->N_SUBS;i++) {
            M = sub[i].dim_X[0];
            N = sub[i].dim_X[1];
            calAx(Xb,sub[i].X,pop->beta,(const int)M,(const int)N);
            for (int j=0;j<M;j++)
                tmp += (sub[i].beta[j]-Xb[j])*(sub[i].beta[j]-Xb[j]);
        }
        double tmp2 = rgamma(ALPHA + 0.5*(double)pop->N_SUBS*pop->Nb,BETA + 0.5*tmp,seed);
        if (tmp2 > 0.01)
            pop->re_prec[is] = tmp2;
//        printf("%lf\n",pop->re_prec[is]);
    }
    free(Xb);
}

/*void draw_reprec(POP *pop,SUB *sub,unsigned long *seed) {
    int N,M;
    double ALPHA, BETA,*tmp;
    double *Xb;
    ALPHA = 1.1;
    BETA  = 1.1;
    
    M = pop->Ns*pop->Nb;
    N = pop->Ncov*M;
    Xb = (double *)calloc(M,sizeof(double));
    tmp = (double *)calloc(M,sizeof(double));
    

    for (int isub=0;isub<pop->N_SUBS;isub++) {
        calAx(Xb,sub[isub].X,pop->beta,(const int)M,(const int)N);
        for (int i=0;i<M;i++)
            tmp[i] += (sub[isub].beta[i]-Xb[i])*(sub[isub].beta[i]-Xb[i]);
    }
    for (int i=0;i<M;i++)
        pop->re_prec[i] = rgamma(ALPHA + 0.5*(double)pop->N_SUBS,BETA + 0.5*tmp[i],seed);
    
    free(Xb);
    free(tmp);

}*/

void draw_beta_eta(POP *pop,SUB *sub,REP *rep,unsigned long *seed) {
    int i,j,ncol;
    double *M,*YY,*V,*VpV,*J,*X,*betaeta;
          
    ncol = (rep->dim_X[1] + rep->dim_V[1]);
    
    betaeta = (double *)calloc(ncol,sizeof(double));

    V = (double *)calloc(rep->dim_V[0]*ncol,sizeof(double));
    
    J = (double *)calloc(rep->dim_V[0]*ncol,sizeof(double));
   
    X = (double *)calloc(rep->dim_V[0]*ncol,sizeof(double));

    for (i=0;i<rep->dim_V[0];i++)
        for (j=0;j<rep->dim_X[1];j++)
            X[i*ncol+j] = rep->X[i*rep->dim_X[1]+j];
     for (i=0;i<rep->dim_V[0];i++)
        for (j=rep->dim_X[1];j<ncol;j++)
            X[i*ncol+j] = rep->V[i*rep->dim_V[1]+j-rep->dim_X[1]];
    int cnt = 0;
    for (i=0;i<rep->dim_X[1];i++)
        betaeta[i] = sub->beta[i];
    for (i=rep->dim_X[1];i<ncol;i++)
        betaeta[i] = rep->eta[i-rep->dim_X[1]];
 
     calH(J,X,rep->delta,(const int)rep->dim_V[0],(const int)ncol,(const int)rep->P);
  
    /* calculate H(X,delta)beta */
    /* calculate J(V,delta)eta */
    double *Jbe;
    Jbe = (double *)calloc(rep->dim_V[0],sizeof(double));
    calAx(Jbe,J,betaeta,(const int)rep->dim_V[0],(const int)ncol);
    
    free(betaeta);
    /* calculate Y - Wdelta - J(V,delta)eta - H(X,delta)beta*/
    
    calculate_res5(rep,rep->P);

    for (i=0;i<rep->dim_V[0];i++)
        rep->residuals5[i] -= Jbe[i];
       
    free(Jbe);
   
    M = (double *)calloc(ncol,sizeof(double));
    VpV = (double *)calloc(ncol*ncol,sizeof(double));   
  
    /* Add X and -J(V,delta)*/
    
    for (i=0;i<rep->dim_V[0];i++)
        for (j=0;j<ncol;j++)
            V[i*ncol+j] = X[i*ncol+j] - J[i*ncol+j];
            
    calApVinvX(M,V,rep->residuals5,rep->d_Y,rep->dim_V[0],ncol);
   
    /* calculate  (X - J)'(X - J) */
    
    calApVinvA(VpV,V,rep->d_Y,(const int)rep->dim_V[0],(const int)ncol);
   
    for (i=rep->dim_X[1];i<ncol;i++)
        VpV[i*ncol+i] += rep->preceta;

    if (pop->non_parm || pop->GRP) {
        int Nrow = rep->dim_X[1];
        int Ncol = pop->Ncov*Nrow;
        double *Xb = (double *)calloc(Nrow,sizeof(double));
        calAx(Xb,sub->X,pop->beta,(const int)Nrow,(const int)Ncol);
         for (int is=0;is<pop->Ns;is++) {
            for (int i=0;i<pop->Nb;i++) {
                    M[i+is*pop->Nb] += pop->re_prec[is]*Xb[i+is*pop->Nb];
                    VpV[(i+is*pop->Nb)*ncol+i+is*pop->Nb] += pop->re_prec[is];
            }
        }
        free(Xb);
   }
   
    int err = cholesky_decomp2vec(VpV,ncol);
    
    double *mean;
    mean = (double *)calloc(ncol,sizeof(double));
   
   if (err) {  // err = 1 means P is SPD
        err = forward_substitution2vec(VpV,M,ncol);
        err = cholesky_backsub2vec(VpV,M,ncol);
        err = rmvnorm3vec(mean,VpV,ncol,M,seed,1);
    }
    else {
        printf("error in draw_beta_eta, precision is not SPD\n");
        exit(0);
    }

    /*  Will need to change this so that each rep has its own beta, for cases with multiple reps per subject. And r.e. about sub->beta and r.e. about pop->beta. CODE changes required */

    for (i=0;i<rep->dim_X[1];i++)
        sub->beta[i] = mean[i];
    for (i=rep->dim_X[1];i<ncol;i++)
        rep->eta[i-rep->dim_X[1]] = mean[i];
        
    /* calculate Veta and Xbeta */
    
    calAx(rep->Veta,rep->V,rep->eta,(const int)rep->dim_V[0],(const int)rep->dim_V[1]);
    calAx(rep->Xbeta,rep->X,sub->beta,(const int)rep->dim_X[0],(const int)rep->dim_X[1]);

    /* update W and Wdelta */
    calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],rep->P);
    calWdelta(rep->Wdelta,rep->W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);

    free(J);
    free(X);
    free(VpV);
    free(M);
    free(mean);
    free(V);
  }


void draw_pop_beta(POP *pop,SUB *sub,unsigned long *seed) {
    double *P,*M,*tmpM,*tmpP,*reprec;
    int N = pop->Ncov*pop->Ns*pop->Nb;
    int SB = pop->Ns*pop->Nb;
 
    M = (double *)calloc(N,sizeof(double));
    tmpM = (double *)calloc(N,sizeof(double));
    P = (double *)calloc(N*N,sizeof(double));
    tmpP = (double *)calloc(N*N,sizeof(double));
    reprec = (double *)calloc(SB,sizeof(double));
 
    for (int i=0;i<pop->Ns;i++)
        for (int j=0;j<pop->Nb;j++)
            reprec[i*pop->Nb+j] = pop->re_prec[i];

    for (int isub=0;isub<pop->N_SUBS;isub++) {
        calApVinvX(tmpM,sub[isub].X,sub[isub].beta,reprec,(const int)SB,(const int)N);
        for (int i=0;i<N;i++)
            M[i] += tmpM[i];
    }

    for (int isub=0;isub<pop->N_SUBS;isub++) {
        calApVinvA(tmpP,sub[isub].X,reprec,(const int)SB,(const int)N);
        for (int i=0;i<N;i++)
            for (int j=0;j<N;j++)
                P[i*N+j] += tmpP[i*N+j];
    }
    
/*    for (int is=0;is<pop->Ns;is++) {
        for (int isub=0;isub<pop->N_SUBS;isub++) {
            for (int j=0;j<pop->Nb;j++) {
                P[(j+is*pop->Nb)*N + (j+is*pop->Nb)] += pop->re_prec[j+is*pop->Nb];
            }
        }
    }*/

    for (int i=0;i<N;i++)
        P[i*N+i] += 0.1;
    int err = cholesky_decomp2vec(P,N);
    
    if (err) {  // err = 1 means P is SPD
        err = forward_substitution2vec(P,M,N);
        err = cholesky_backsub2vec(P,M,N);
        err = rmvnorm3vec(pop->beta,P,N,M,seed,1);
    }
    else {
        printf("error in draw_pop_beta, precision is not SPD\n");
        exit(0);
    }
    free(reprec);
    free(P);
    free(M);
    free(tmpM);
    free(tmpP);
}



void DIC(POP *pop,REP *rep,unsigned long *seed) {
    int i;
    double tmp;
    double mlogsqrt2pi = -0.9189385332046727;
    
//      calWdelta(rep->Wm,rep->W,rep->m,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
//     calculate_marginal_residuals(rep,rep->P);
    calculate_residuals(rep,rep->P);
    tmp = 0;
    for (i=maxP;i<rep->dim_X[0];i++)
            tmp += mlogsqrt2pi + sqrt(rep->d_Y[i]) - 0.5*rep->d_Y[i]*(rep->residuals[i])*(rep->residuals[i]);

    pop->ED += -2*tmp;
}
