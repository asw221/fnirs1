/*
 *  mcmc.c
 *  LGCP
 *
 *  Created by Timothy D. Johnson on 11/10/15.
 *  Copyright 2015 University of Michigan. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <omp.h>
#include "randgen.h"
#include "cholesky.h"
#include <time.h>
#include "fNIRS.h"

extern int DEBUG;
extern int MAX_ITER;
extern int BURN_IN;
extern int maxP;
extern int nthd_sys;
extern FILE *flog;
extern int PPP;

void adjust_MCMC_proposal_dist(double target,double alpha,double *eps,int iter)
// target is the global acceptance probability target
// alpha is the MCMC proposal probability
// gamma controls how much movement aroung mu
// eps0 is initial sd
// *eps is update sd
// iter is the iteration number
{
    static double kappa = 0.75;
    
    double tmp = pow((double)(iter),-kappa);
    double tmp2 = log(*eps);
    
    tmp2 = tmp2 - tmp*(target-alpha);
    *eps = exp(tmp2);
}

//void update_prop_dist(REP *rep,FILE *fwhy,int iter)
void update_prop_dist(REP *rep,int iter)
{
    if ((iter <= BURN_IN)) {
 
        adjust_MCMC_proposal_dist(0.4,rep->alpha_MH[0],&rep->prop_sd[0],iter);
        adjust_MCMC_proposal_dist(0.4,rep->alpha_MH[1],&rep->prop_sd[1],iter);
                
        if (iter == BURN_IN)
            rep->attempt[0] = rep->accept[0] = 0;          
    }
//    fprintf(fwhy,"%lf %lf ",sub->prop_sd[0],sub->prop_sd[1]);
//    fprintf(fwhy,"\n");
}
 

void mcmc(POP *pop,unsigned long *seed) {
    int i,j,iter,isub,irep;
    unsigned long seed2[3];
    SUB *sub;
    REP *rep;
 
 
    void draw_knot_locations(REP *rep,int sdegree,unsigned long *seed);
    int knot_birth_death(REP *rep,int pknots,const int sdegree,int iter,unsigned long *seed);
 
    void DLM(REP *rep,int iter,unsigned long *seed);
    void draw_precYstart(REP *rep,unsigned long *seed);
    void draw_precY(REP *rep,unsigned long *seed);
    void draw_reprec(POP *pop,SUB *sub,unsigned long *seed);
    void draw_re_rep_prec(POP *pop,SUB *sub,unsigned long *seed);
 
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    
    void draw_beta_eta(POP *pop,SUB *sub,unsigned long *seed);
    void draw_pop_beta(POP *pop,SUB *sub,unsigned long *seed);
    void draw_sub_beta(POP *pop,SUB *sub,unsigned long *seed);    
    double DLMloglik(REP *rep,double sfreq,int,unsigned long *seed);
  
   
    pop->ED = 0;
       
    int sdegree = 4;

    int NN = 0;
    for (i=0;i<pop->N_SUBS;i++) {
        for (j=0;j<pop->sub[i].N_REPS;j++) {
            NN++;
        }
    }
    
    int *subidx = (int *)calloc(NN,sizeof(int));
    REP **tmprep = (REP **)calloc(NN,sizeof(REP *));
    SUB **tmpsub = (SUB **)calloc(pop->N_SUBS,sizeof(SUB *));
    for (i=0;i<pop->N_SUBS;i++) {
        sub = &(pop->sub[i]);
        tmpsub[i] = sub;
    }
    int ij=0;
    for (i=0;i<pop->N_SUBS;i++) {
        sub = &(pop->sub[i]);
        for (j=0;j<pop->sub[i].N_REPS;j++) {
            rep = &(sub->rep[j]);            
            tmprep[ij] = rep;
            subidx[ij] = i;
            ij++;
        }
    }
 
/*    int nthd;
    if (pop->N_SUBS > 1) 
        nthd = nthd_sys;
    else
        nthd = 1;
    omp_set_dynamic(0);
    fprintf(flog,"%d\n",nthd);
    omp_set_num_threads(1);
    */
    int *v = (int *)malloc(NN*sizeof(int));
    
    for (iter=0;iter<=MAX_ITER;iter++) {
        if (!(iter%100)) { 
            printf("%d",iter);
            fflush(stdout);
        }
        if (!(iter%20) && (iter%100)) {
            printf(".");
            fflush(stdout);
        }
       

 /*      for (int i=0;i<NN;i++)
            v[i] = i;
        permute_sample_int(v,NN,seed);
*/
        for (int i=0;i<3;i++)
            seed2[i] = seed[i];
            
 //       #pragma omp parallel for firstprivate(seed2) private(j)
//printf("iter = %d\n",iter);fflush(NULL);
        for (i=0;i<pop->N_SUBS;i++) {
            for (j=0;j<pop->sub[i].N_REPS;j++) {
//printf("0 %s \n",pop->sub[i].rep[j].dataname);fflush(NULL);
                if (pop->sub[i].rep[j].dim_X[0] != pop->sub[i].rep[j].dim_V[0])
                    printf("%s X[0] = %d V[0] = %d\n",pop->sub[i].rep[j].dataname,pop->sub[i].rep[j].dim_X[0],pop->sub[i].rep[j].dim_V[0]);
                if (pop->sub[i].rep[j].dim_X[0] != pop->sub[i].rep[j].dim_W[0])
                    printf("%s X[0] = %d W[0] = %d\n",pop->sub[i].rep[j].dataname,pop->sub[i].rep[j].dim_X[0],pop->sub[i].rep[j].dim_W[0]);
                fflush(NULL);     
                knot_birth_death(&(pop->sub[i].rep[j]),pop->knots,sdegree,iter,seed);
            }
        }  
//printf("\n");fflush(NULL);              
//        #pragma omp parallel for firstprivate(seed2) private(j)
        for (i=0;i<pop->N_SUBS;i++) {
            for (j=0;j<pop->sub[i].N_REPS;j++) {
//printf("A %s \n",pop->sub[i].rep[j].dataname);fflush(NULL);
                int imax = pop->sub[i].rep[j].nKnots-sdegree*2;
                imax = (int)ceil((double)imax/4.);
                for (int k=0;k<imax;k++)
                    draw_knot_locations(&(pop->sub[i].rep[j]),sdegree,seed);
            }
        }
//        #pragma omp parallel for firstprivate(seed2) private(j)
//printf("A\n");fflush(NULL);

        for (i=0;i<pop->N_SUBS;i++) {
// printf("B %s \n",pop->sub[i].rep[0].dataname);fflush(NULL);
 //           for (j=0;j<pop->sub[i].N_REPS;j++) {
                draw_beta_eta(pop,&(pop->sub[i]),seed);
 //           }
        }
//printf("C\n");fflush(NULL);
  
//        #pragma omp parallel for firstprivate(seed2) private(j)
        for (i=0;i<pop->N_SUBS;i++) {
            for (j=0;j<pop->sub[i].N_REPS;j++) {
//printf("D %s \n",pop->sub[i].rep[j].dataname);fflush(NULL);
                draw_precYstart(&(pop->sub[i].rep[j]),seed);
            }
        }
       
 //       #pragma omp parallel for firstprivate(seed2) private(j)
        for (i=0;i<pop->N_SUBS;i++) {
            for (j=0;j<pop->sub[i].N_REPS;j++) {
//printf("E %s \n",pop->sub[i].rep[j].dataname);fflush(NULL);
                DLM(&(pop->sub[i].rep[j]),iter,seed);
                if (iter > 0)
                    update_prop_dist(&(pop->sub[i].rep[j]),iter);
            }
        }
        
        if (pop->GRP) {
//printf("F\n");fflush(NULL);

            draw_re_rep_prec(pop,pop->sub,seed);
//printf("D\n");fflush(NULL);
            
            int isub2;
            if (!pop->No_replicates) {
//printf("D1\n");fflush(NULL);
//                #pragma omp parallel for
                for (isub2=0;isub2<pop->N_SUBS;isub2++)
                    draw_sub_beta(pop,tmpsub[isub2],seed);
 
                draw_reprec(pop,pop->sub,seed);
            } 
//printf("E\n");fflush(NULL);
            draw_pop_beta(pop,pop->sub,seed);
//printf("F\n");fflush(NULL);
        }
        
/*        if (!(iter%50) && (iter > 0)) {
            for (isub=0;isub<pop->N_SUBS;isub++) {
                sub = &(pop->sub[isub]);
                
                for (irep=0;irep<sub->N_REPS;irep++) {
                    rep = &(sub->rep[irep]);
 
                    if ((iter <= BURN_IN) && (iter > 1)) {
                        double rt;
                     
                        rt = (double)rep->accept[0]/(double)rep->attempt[0];
                        adjust_acceptance2(rt,&(rep->prop_sd[0]),0.35);
                        rep->accept[0] = rep->attempt[0] = 0;

                        rt = (double)rep->accept[1]/(double)rep->attempt[1];
                        adjust_acceptance2(rt,&(rep->prop_sd[1]),0.35);
                        rep->accept[1] = rep->attempt[1] = 0;
                        
                        rt = (double)rep->accept[3]/(double)rep->attempt[3];
                        adjust_acceptance2(rt,&(rep->prop_sd[3]),0.35);
                        rep->accept[3] = rep->attempt[3] = 0;                     

                    }
                }
            }
        }*/
        
        if (!(iter%10)) {
            for (isub=0;isub<pop->N_SUBS;isub++) {
                sub = &(pop->sub[isub]);
                
                for (irep=0;irep<sub->N_REPS;irep++) {
                    rep = &(sub->rep[irep]);
                    
                    fprintf(flog,"iter = %6d\t Sub %s, Rep %d \t int knots = %d",iter,rep->dataname,irep,rep->nKnots-2*sdegree);
                    fprintf(flog,"\t beta[1] = %10.6lf \t %.6lf %.6lf %d\n",rep->beta[1],rep->df_delta1,rep->df_delta2,rep->P);
       //             fprintf(flog,"%.6lf %.6lf %d\n",rep->df_delta1,rep->df_delta2,rep->P);
                }
            }

            if (pop->GRP) {
                fprintf(flog,"Iter %6d Pop Parms:\n",iter);
                for (i=0;i<pop->Ncov;i++)
                    fprintf(flog,"%30s = %10.6lf\n",pop->covnames[i],pop->beta[i]);
                fprintf(flog,"\n\n");
            }
        }
        fflush(flog);
        
        if ((iter>BURN_IN)) {
            for (isub=0;isub<pop->N_SUBS;isub++) {
                sub = &(pop->sub[isub]);
                
                for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
                    rep = &(sub->rep[irep]);
 
                    for (i=0;i<rep->dim_X[0];i++) {
                        rep->sd_Y[i] += rep->d_Y[i];
                    }
                    
//                    for (i=0;i<rep->dim_V[1];i++)
//                        fprintf(rep->fout_eta,"%lf ",rep->eta[i]);
//                    fprintf(rep->fout_eta,"\n");
                   
 
                    for (i=0;i<rep->dim_X[0];i++) {
                        rep->mean_fit[i] += rep->Y[i] - rep->residuals[i];
                        rep->mean_res[i] += rep->residuals[i];
                        rep->mean_d_Y[i] += rep->d_Y[i];
                    }
                                       
                     
                    for (i=0;i<rep->dim_X[0];i++)
                        rep->mXbeta[i] += rep->Xbeta[i];
                    for (i=0;i<rep->dim_V[0];i++)
                        rep->mVeta[i] += rep->Veta[i];
                    for (i=0;i<rep->dim_W[0];i++)
                        rep->mWdelta[i] += rep->Wdelta[i];
                }
            }
        }          
         if ((iter>0)) {
            for (isub=0;isub<pop->N_SUBS;isub++) {
                sub = &(pop->sub[isub]);
                
                for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
                    rep = &(sub->rep[irep]);
                  
                    if (DEBUG) {
                        fprintf(rep->fout_dlm,"%lf %lf %d\n",rep->df_delta1,rep->df_delta2,rep->P);
                    
                         fprintf(rep->fout_nknots,"%d ",rep->nKnots-2*sdegree);

                        for (i=sdegree;i<rep->nKnots-sdegree;i++)
                            fprintf(rep->fout_knots,"%lf ",rep->knots[i]);
                        fprintf(rep->fout_knots,"\n");  
                    }
                    
                    for (i=0;i<pop->Nb*(pop->Ns);i++)
                        fprintf(rep->fout_beta,"%lf ",rep->beta[i]);
                    fprintf(rep->fout_beta,"\n");  

                }
                
//                for (i=0;i<pop->Nb*(pop->Ns);i++)
//                    fprintf(sub->fout_beta,"%lf ",sub->beta[i]);
//                fprintf(sub->fout_beta,"\n");
            }

            for (i=0;i<pop->Ncov;i++)
                fprintf(pop->fout_beta,"%lf ",pop->beta[i]);
            fprintf(pop->fout_beta,"\n");
            
 //           for (i=0;i<pop->Nb*pop->Ns;i++)
 //               fprintf(pop->fout_reprec,"%lf ",pop->re_prec[i]);
 //           fprintf(pop->fout_reprec,"\n");

            int len = pop->Nb*pop->Ns*pop->sub[0].N_REPS;
            len *= len;
            for (i=0;i<len;i++)
                fprintf(pop->fout_re_rep_prec,"%lf ",pop->re_rep_prec[i]);
            fprintf(pop->fout_re_rep_prec,"\n");

            fflush(NULL);      
        }
    }
    free(v);
    free(tmprep);
    free(tmpsub);
    free(subidx);
    
    fclose(pop->fout_beta);
//    fclose(pop->fout_reprec);
    fclose(pop->fout_re_rep_prec);

    for (isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
            if (DEBUG) {
                fclose(rep->fout_dlm);
                fclose(rep->fout_eta);
                fclose(rep->fout_nknots);
                fclose(rep->fout_knots);
            }
            fclose(rep->fout_beta);
         }
//        fclose(sub->fout_beta);
    }

 /*   for (isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
            for (i=0;i<rep->dim_W[0];i++) {
                for (j=0;j<rep->dim_W[1];j++) {
                    rep->mdelta[i*rep->dim_W[1]+j] /= (double)(MAX_ITER-BURN_IN);
                    fprintf(rep->fout_delta,"%lf ",rep->mdelta[i*rep->dim_W[1]+j]);
                }
                fprintf(rep->fout_delta,"\n");                  
            }                
            fclose(rep->fout_delta);
        }
    }*/
    
   for (isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
            for (i=0;i<rep->dim_X[0];i++) {
                    rep->sd_Y[i] =  sqrt((double)(MAX_ITER-BURN_IN)/rep->sd_Y[i]);
                    fprintf(rep->fout_stdev,"%lf ",rep->sd_Y[i]);
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
                rep->mean_res[i] /= ((MAX_ITER-BURN_IN));
                rep->mean_fit[i] /= ((MAX_ITER-BURN_IN));
                rep->mean_d_Y[i] /= ((MAX_ITER-BURN_IN));
                rep->std_res[i] = rep->mean_res[i]*sqrt(rep->d_Y[i]);
            }
       }
    }
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
    double ALPHA=3,BETA=2;
    double a=0, b=0;
    
    a = (double)rep->P/2. + ALPHA;
    
    b = 0;
    for (int i=0;i<rep->P;i++)
        b += (rep->Y[i]-rep->Veta[i])*(rep->Y[i]-rep->Veta[i]);
    b = BETA + 0.5*b;
    rep->precYstart = rgamma(a,b,seed);
    for (int i=0;i<rep->P;i++)
        rep->d_Y[i] = rep->precYstart;

}

void draw_precY(REP *rep,unsigned long *seed)
{
    int N,P;
    double ALPHA=0.01,BETA=0.01;
    double a,b,factor;
    void calculate_residuals(REP *rep,int P);
    
    N = rep->dim_X[0];
    P = rep->P;
     
    calculate_residuals(rep,P);

    a = (double)(N-P)/2. + ALPHA;
    
    b = 0;
    for (int i=P;i<N;i++) {
        if (rep->dY[i])
            b += rep->residuals[i]*rep->residuals[i];
        else
            b += rep->residuals[i]*rep->residuals[i]*rep->LY;
    }        
    b = BETA + 0.5*b;
    double tmp = rgamma(a,b,seed);
    
    for (int i=P;i<N;i++) {
        if (rep->dY[i])
            factor = 1;
        else
            factor = rep->LY;

        rep->d_Y[i] = factor*tmp;
    }
}


void calH(double *H,double *A,double *m,const int nrow,const int ncol,const int P) {
    int t,j,k;
    
    for (j=0;j<nrow;j++) 
        for (k=0;k<ncol;k++)
            H[j*ncol+k] = 0;
    
    for (t=P;t<nrow;t++) 
        for (j=0;j<P;j++) 
            for (k=0;k<ncol;k++) 
                H[t*ncol+k] += A[(t-j-1)*ncol+k]*m[t*P+j];
}

void calApVinvX(double *ApX,double *A,double *X,double *Vinv,const int nrow,const int ncol) {
    int i;
    double *At;
    double *Vx;
    void transpose(double *At,double *A,const int nrow,const int ncol);
    void Ax(double *B,double *A,double *x,const int nrow,const int ncol);
    
    At = (double *)calloc(nrow*ncol,sizeof(double));
    
    transpose(At,A,nrow,ncol);
    
    Vx = (double *)malloc(nrow*sizeof(double));
 
    for (i=0;i<nrow;i++)    
        Vx[i] = Vinv[i]*X[i];
        
    for (i=0;i<ncol;i++)
        ApX[i] = 0;
   
    Ax(ApX,At,Vx,ncol,nrow);
                      
    free(At); 
    free(Vx);
}

void calApVinvA(double *ApA,double *A,double *Vinv,const int nrow,const int ncol) {
// only for diagonal matrices, V
    int i,j;
    double *At,*VA;
    void transpose(double *At,double *A,const int nrow,const int ncol);
    void AxB(double *C,double *A,double *B,const int nrowA,const int mid,const int ncolB);
        
    At = (double *)calloc(nrow*ncol,sizeof(double));
    VA = (double *)calloc(nrow*ncol,sizeof(double));
    
    for (i=0;i<nrow;i++)
        for (j=0;j<ncol;j++)
            VA[i*ncol+j] = Vinv[i]*A[i*ncol+j];
            
    transpose(At,A,nrow,ncol);
        
    for (i=0;i<ncol*ncol;i++)
        ApA[i] = 0;

    AxB(ApA,At,VA,ncol,nrow,ncol);

        
    free(At);
    free(VA);

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
    int i,j;
        
    for (i=0;i<nrow;i++) {
        int min = (i < ncol) ? i:ncol;
        for (j=0;j<min;j++) {
            W[i*ncol+j] = Y[i-j-1] - Xb[i-j-1] - Ve[i-j-1];
       }
    }
}

void calW2(double *W,double *Y,const int nrow,const int ncol,const int P) {
    int i,j;
        
    for (i=0;i<nrow;i++) {
        int min = (i < ncol) ? i:ncol;
        for (j=0;j<min;j++) {
            W[i*ncol+j] = Y[i-j-1];
       }
    }
}

void calculate_residuals(REP *rep,int P) {
    int i;
    
    for (i=0;i<P;i++)
        rep->residuals[i] = rep->Y[i] - rep->Veta[i] - rep->Xbeta[i];
    for (i=P;i<rep->dim_X[0];i++)
        rep->residuals[i] = rep->Y[i] - (rep->Veta[i] + rep->Wdelta[i] + rep->Xbeta[i]);
    
}

void calculate_residualsA(double *res,double *Y,double *Xb,double *Ve,double *Wd,int P,int length) {
    int i;
    
    for (i=0;i<P;i++)
        res[i] = Y[i] - Ve[i] - Xb[i];
    for (i=P;i<length;i++)
        res[i] = Y[i] - Ve[i] - Wd[i] - Xb[i];
    
}

void calculate_marginal_residuals(REP *rep,int P) {
    int i;
    
    for (i=0;i<P;i++)
        rep->residuals[i] = rep->Y[i] - rep->Veta[i] - rep->Xbeta[i];
    for (i=P;i<rep->dim_X[0];i++)
        rep->residuals[i] = rep->Y[i] - (rep->Veta[i] + rep->Wm[i] + rep->Xbeta[i]);
    
}

void calculate_res3(REP *rep) {
    for (int i=0;i<rep->dim_X[0];i++) {
        rep->residuals3[i] = rep->Y[i] -  rep->Xbeta[i] - rep->Veta[i];
    }
}

void calculate_res5(REP *rep,int P) {
    for (int i=0;i<P;i++) {
        rep->residuals5[i] = rep->Y[i];
    }
    for (int i=P;i<rep->dim_V[0];i++)
        rep->residuals5[i] = rep->Y[i] - rep->Wdelta[i];
}

inline void marg_res2(double *res,double *Y,double *delta,const int N,const int P)
{
    double Wd;
    
    for (int t=0;t<P;t++)
        res[t] = Y[t];
    
    for (int t=P;t<N;t++) {
        Wd = 0;
        for (int i=0;i<P;i++)
            Wd +=  Y[t-i-1]*delta[t*P+i];
        res[t] = Y[t] - Wd;
    }
}

void draw_preceta(REP *rep,unsigned long *seed) {
    double ALPHA=0, BETA=0;
    double a, b;
    ALPHA = 0.1;
    BETA =  0.1;
    b = 0;
    
    for (int i=0;i<rep->dim_V[1];i++)
            b += rep->eta[i]*rep->eta[i];
    
    a = 0.5*rep->dim_V[1] + ALPHA;
    
    b = 0.5*b + BETA;
    
    double tmp = rgamma(a,b,seed);
    if ((.01 > tmp) && (tmp < 10.))
        rep->preceta = tmp;
}


void draw_reprec(POP *pop,SUB *sub,unsigned long *seed) {
    int N,M;
    double ALPHA, BETA,*rss;
    double *Xb;
    ALPHA = 1;//3;
    BETA  = 1;//2;
    
    M = pop->Ns*pop->Nb;
    N = pop->Ncov;
    Xb = (double *)calloc(M,sizeof(double));
    rss = (double *)calloc(M,sizeof(double));
    
    for (int isub=0;isub<pop->N_SUBS;isub++) {
        calAx(Xb,sub[isub].X,pop->beta,(const int)M,(const int)N);
        for (int i=0;i<M;i++)
            rss[i] += (sub[isub].beta[i]-Xb[i])*(sub[isub].beta[i]-Xb[i]);
    }
    for (int i=0;i<M;i++)   
        pop->re_prec[i] =  rgamma(ALPHA + 0.5*(double)pop->N_SUBS,BETA + 0.5*rss[i],seed);    
    free(Xb);
    free(rss);
}


void draw_re_rep_prec(POP *pop,SUB *sub,unsigned long *seed) {
    int N,M;
    double ALPHA, BETA,*Sinv,*Xb;
    REP *rep;
    
    int ntasks = sub[0].N_REPS;
    M = pop->Nb*pop->Ns*ntasks;
    N = pop->Ncov;
    Xb = (double *)calloc(M,sizeof(double));
       
    Sinv = (double *)calloc(M*M,sizeof(double));
    for (int i=0;i<M;i++)
        Sinv[i+i*M] = 1.;
        

    for (int i=0;i<pop->N_SUBS;i++) {
        if (pop->No_replicates) {
            calAx(Xb,sub[i].X,pop->beta,(const int)M,(const int)N);
            double *beta = (double *)calloc(M,sizeof(double));
            int cnt = 0;
            for (int t=0;t<ntasks;t++) {
                if (t>0)
                    cnt += pop->Nb*pop->Ns;
                for (int j=0;j<pop->Nb*pop->Ns;j++)
                    beta[j+t*cnt] = sub[i].rep[t].beta[j];
            }        
            for (int j=0;j<M;j++)
                for (int k=0;k<M;k++) 
                    Sinv[j*M + k] += (beta[j]-Xb[j])*(beta[k]-Xb[k]);  
            free(beta); 
        }
/*        else {
            for (int k=0;k<sub[i].N_REPS;k++) {
                rep = &(sub[i].rep[k]);
                rss += (rep->beta[j+is*pop->Nb]-sub[i].beta[j+is*pop->Nb])*(rep->beta[j+is*pop->Nb]-sub[i].beta[j+is*pop->Nb]);
                cnt++;
            }*/
    }
    rwishart(pop->re_rep_prec,Sinv,M,pop->N_SUBS+M,seed);
    
    free(Sinv);
    free(Xb);
}

void draw_beta_eta(POP *pop,SUB *sub,unsigned long *seed) {
    int i,j,ncol,*nc;
    double *M,*V,*VpV,*J,*X,*betaeta,*mean,*prec;
    int ntasks = sub->N_REPS;
    REP *rep;
    
    nc = (int *)calloc(ntasks,sizeof(int));
    ncol = 0;
    for (int t=0;t<ntasks;t++) {
        rep = &(sub->rep[t]);
        nc[t] = rep->dim_X[1]+rep->dim_V[1];
        ncol += nc[t]; 
    }
    mean = (double *)calloc(ncol,sizeof(double));
    prec = (double *)calloc(ncol*ncol,sizeof(double));
    
    int cnt = 0;
//   printf("0\n");fflush(NULL);
    for (int t=0;t<ntasks;t++) {
        rep = &(sub->rep[t]);
//      printf("t = %d\n",t);fflush(NULL);
        betaeta = (double *)calloc(nc[t],sizeof(double));
        
        V = (double *)calloc(rep->dim_V[0]*nc[t],sizeof(double));
    
        J = (double *)calloc(rep->dim_V[0]*nc[t],sizeof(double));
   
        X = (double *)calloc(rep->dim_V[0]*nc[t],sizeof(double));

        for (i=0;i<rep->dim_V[0];i++)
            for (j=0;j<rep->dim_X[1];j++)
                X[i*nc[t]+j] = rep->X[i*rep->dim_X[1]+j];
        for (i=0;i<rep->dim_V[0];i++)
            for (j=rep->dim_X[1];j<nc[t];j++)
                X[i*nc[t]+j] = rep->V[i*rep->dim_V[1]+j-rep->dim_X[1]];
  
        for (i=0;i<rep->dim_X[1];i++)
            betaeta[i] = rep->beta[i];
        for (i=rep->dim_X[1];i<nc[t];i++)
            betaeta[i] = rep->eta[i-rep->dim_X[1]];
 
        calH(J,X,rep->delta,(const int)rep->dim_V[0],(const int)nc[t],(const int)rep->P);
  
        /* calculate J(V,delta)eta */
        double *Jbe;
        Jbe = (double *)calloc(rep->dim_V[0],sizeof(double));
        calAx(Jbe,J,betaeta,(const int)rep->dim_V[0],(const int)nc[t]);
    
        free(betaeta);
        // calculate Y - Wdelta - J(V,delta)eta
    
        calculate_res5(rep,rep->P);

        for (i=0;i<rep->dim_V[0];i++)
            rep->residuals5[i] -= Jbe[i];

        free(Jbe);
    
        M = (double *)calloc(nc[t],sizeof(double));
        VpV = (double *)calloc(nc[t]*nc[t],sizeof(double));   
  
        /* Add X and -J(V,delta)*/
    
        for (i=0;i<rep->dim_V[0];i++)
            for (j=0;j<nc[t];j++)
                V[i*nc[t]+j] = X[i*nc[t]+j] - J[i*nc[t]+j];
            
        calApVinvX(M,V,rep->residuals5,rep->d_Y,rep->dim_V[0],nc[t]);

        /* calculate  (X - J)'(X - J) */
    
        calApVinvA(VpV,V,rep->d_Y,(const int)rep->dim_V[0],(const int)nc[t]);
        free(V);
        free(J);
        free(X);
   
        for (i=rep->dim_X[1];i<nc[t];i++)
            VpV[i*nc[t]+i] += rep->preceta;
            
        // now combine M and VpV from the tasks
        if (t>0)
            cnt += nc[t-1];
        for (i=0;i<nc[t];i++)
            mean[i+cnt] = M[i];
        for (i=0;i<nc[t];i++)
            for (j=0;j<nc[t];j++) 
                prec[(i+cnt)*ncol+(j+cnt)] = VpV[i*nc[t]+j];


        free(M);
        free(VpV);
    }
//   printf("1\n");fflush(NULL);
    if (pop->non_parm || pop->GRP) {
        if (pop->No_replicates) {
            int nrow = pop->Nb*pop->Ns;
            int Nrow = pop->Nb*pop->Ns*ntasks;
            int Ncol = pop->Ncov;
            double *Xb = (double *)calloc(Nrow,sizeof(double));
            calAx(Xb,sub->X,pop->beta,(const int)Nrow,(const int)Ncol);
            double *tmp = (double *)calloc(Nrow,sizeof(double));
            calAx(tmp,pop->re_rep_prec,Xb,(const int)Nrow,(const int)Nrow);
            cnt = 0;
            Ncol = Nrow;
            for (int t=0;t<ntasks;t++) {
                if (t>0)
                    cnt += nc[t-1];
                for (i=0;i<nrow;i++)
                    mean[i+cnt] += tmp[i+t*nrow];
//                printf("t = %d a\n",t);fflush(NULL);   
                int cnt2 = 0;
                for (int t2=0;t2<ntasks;t2++) {
                    if (t2>0)
                        cnt2 += nc[t2-1];
                    for (i=0;i<nrow;i++) 
                        for (j=0;j<nrow;j++)
                            prec[(i+cnt)*ncol + (j+cnt2)] += pop->re_rep_prec[(i+t*nrow)*Nrow + (j+t2*nrow)]; 
          
 //                       printf("%s %d %d %d %d \t %d %d %d %d\n",sub->rep[t].dataname,ncol,ncol,(i+cnt),j+cnt2,Nrow,Nrow,(i+t*nrow),j+t2*nrow);
                }
//                printf("t = %d b\n",t);fflush(NULL);   
            }
//            printf("\n");fflush(NULL);
/*            double *dest = (double *)calloc(nrow*nrow,sizeof(double));
            double *dest2 = (double *)calloc(ncol*ncol,sizeof(double));
            for (int r=0;r<ntasks;r++) {
                for (int c=0;c<ntaskts;c++) {
                    start_row = r*nrow;
                    start_col = c*nrow;
                    get_submatrix(dest,pop->re_rep_prec,(const int)Nrow,(const int)Nrow,(const int)start_row,(const int)start_col,(const int)nrow,(const int)nrow);
                    start_row = 
                    start_col =
                    put_submatrix(dest2,dest,(const int)ncol,(const int)ncol,(const int)start_row,(const int)start_col,(const int)nrow,(const int)nrow);
                }
            }
            adddoubles(prec,prec,dest2,ncol*ncol,1);

            free(dest);
            free(dest2);*/
/*            for (int is=0;is<pop->Ns;is++) {
                for (int i=0;i<pop->Nb;i++) {
                    M[i+is*pop->Nb] += pop->re_rep_prec[is*pop->Nb+i]*Xb[i+is*pop->Nb];
                    VpV[(i+is*pop->Nb)*ncol+i+is*pop->Nb] += pop->re_rep_prec[is*pop->Nb+i];
                }
            }*/
            free(Xb);
            free(tmp);
        }
        else {
            for (int is=0;is<pop->Ns;is++) {
                for (int i=0;i<pop->Nb;i++) {
                    M[i+is*pop->Nb] += pop->re_rep_prec[is*pop->Nb+i]*sub->beta[i+is*pop->Nb];
                    VpV[(i+is*pop->Nb)*ncol+i+is*pop->Nb] += pop->re_rep_prec[is*pop->Nb+i];
                }
            }
        }
    }
//printf("2\n");fflush(NULL);
//    double *MM  = (double *)calloc(ncol,sizeof(double));
//    V = (double *)calloc(ncol*ncol,sizeof(double));
//    for (int i=0;i<ncol;i++) MM[i] = mean[i];
//    for (int i=0;i<ncol*ncol;i++)
//        V[i] = prec[i];
    int err = cholesky_decomp2vec(prec,ncol);
//printf("3\n");fflush(NULL);

    double *mean2;
    mean2 = (double *)calloc(ncol,sizeof(double));
   
   if (err) {  // err = 1 means P is SPD
        err = forward_substitution2vec(prec,mean,ncol);
        err = cholesky_backsub2vec(prec,mean,ncol);
        err = rmvnorm3vec(mean2,prec,ncol,mean,seed,1);
    }
    else {
        printf("error in draw_beta_eta, precision is not SPD, Data name = %s\n",rep->dataname);
        fprintf(flog,"error in draw_beta_eta, precision is not SPD, Data name = %s\n",rep->dataname);
        fprintf(flog,"prec_eta = %lf\n",rep->preceta);
        fprintf(flog,"Precision matrix dim = %d x %d\n",ncol,ncol);
        fprintf(flog,"Precision matrix data\n");
 /*       for (int i=0;i<ncol*ncol;i++)
            fprintf(flog,"%lf ",V[i]);
        fprintf(flog,"\n \n V matrix %d %d\n",rep->dim_V[0],rep->dim_V[1]);
        for (i=0;i<rep->dim_V[0]*rep->dim_V[1];i++)
            fprintf(flog,"%lf ",rep->V[i]);*/
        fprintf(flog,"\n");fflush(NULL);
        fprintf(flog,"re_rep_prec = \n");
        for (int i=0;i<pop->Ns*pop->Nb*ntasks;i++)
            fprintf(flog,"%lf ",pop->re_rep_prec[i]);
        fprintf(flog,"\n ");
   /*     fprintf(flog,"d_Y = ");
        for (int i=0;i<rep->dim_V[0];i++)
            fprintf(flog,"%lf ",rep->d_Y[i]);
        fprintf(flog,"\n");fflush(NULL);
        fprintf(flog,"residuals = ");
        for (int i=0;i<rep->dim_V[0];i++)
            fprintf(flog,"%lf ",rep->residuals5[i]);
        fprintf(flog,"\n");fflush(NULL);*/
        exit(0);
    }
//printf("4\n");fflush(NULL);

    cnt = 0;
    for (int t=0;t<ntasks;t++) {
        rep = &(sub->rep[t]);
        if (t>0)
            cnt += nc[t-1];

        for (i=0;i<rep->dim_X[1];i++)
            rep->beta[i] = mean2[i+t*cnt];
    
        for (i=rep->dim_X[1];i<nc[t];i++)
           rep->eta[i-rep->dim_X[1]] = mean2[i+t*cnt];
        
        /* calculate Veta and Xbeta */
        calAx(rep->Veta,rep->V,rep->eta,(const int)rep->dim_V[0],(const int)rep->dim_V[1]);
        calAx(rep->Xbeta,rep->X,rep->beta,(const int)rep->dim_X[0],(const int)rep->dim_X[1]);

        /* update W and Wdelta */
        calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],rep->P);
        calWdelta(rep->Wdelta,rep->W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
    }
 //   printf("5\n");fflush(NULL);
//    free(VpV);
    free(mean2);
    free(prec);
    free(mean);
//    free(V);
//    free(MM);
    free(nc);
}

void draw_sub_beta(POP *pop,SUB *sub,unsigned long *seed) {

    int N = pop->Nb*pop->Ns;
    double *M = (double *)calloc(N,sizeof(double));
    double *P = (double *)calloc(N*N,sizeof(double));
    double *beta = (double *)calloc(N,sizeof(double));
    void calAx(double *b,double *A,double *x,const int nrow,const int ncol);
       
    
    for (int i=0;i<N;i++)
        P[i*N+i] = sub->N_REPS*pop->re_rep_prec[i] + pop->re_prec[i];
        
    for (int i=0;i<N;i++)
        for (int j=0;j<sub->N_REPS;j++)
            M[i] += pop->re_rep_prec[i]*sub->rep[j].beta[i];
            
    calAx(beta,sub->X,pop->beta,(const int)sub->dim_X[0],(const int)sub->dim_X[1]);        

    for (int i=0;i<N;i++)
        M[i] += pop->re_prec[i]*beta[i];

    int err = cholesky_decomp2vec(P,N);
    
    if (!err) {  // err = 1 means P is SPD
        printf("error in draw_sub_beta, precision is not SPD\n");
        exit(0);
    }
    err = forward_substitution2vec(P,M,N);
    err = cholesky_backsub2vec(P,M,N);
    err = rmvnorm3vec(sub->beta,P,N,M,seed,1);
    
    free(M);
    free(P);
    free(beta);
}

void draw_pop_beta(POP *pop,SUB *sub,unsigned long *seed) {
    double *P,*M,*tmpM,*tmpP;
    int N = pop->Ncov;
    int SB = pop->Nb*pop->Ns;
    double dotProd(double *x,double *y,const int dim);
    void Ax(double *B,double *A,double *x,const int nrow,const int ncol);
    void transpose(double *At,double *A,const int nrow,const int ncol);
    void AxB(double *C,double *A,double *B,const int nrowA,const int mid,const int ncolB);

  
    M = (double *)calloc(N,sizeof(double));
    tmpM = (double *)calloc(N,sizeof(double));
    P = (double *)calloc(N*N,sizeof(double));
    tmpP = (double *)calloc(N*N,sizeof(double));
     
    for (int isub=0;isub<pop->N_SUBS;isub++) {
        if (pop->No_replicates) { 
            int SB2 = SB*sub[isub].N_REPS;
            double *beta = (double *)calloc(SB2,sizeof(double));
            int cnt = 0;
            for (int t=0;t<sub[isub].N_REPS;t++) {
                if (t>0)
                    cnt += SB;
                for (int i=0;i<SB;i++)
                    beta[i+t*cnt] = sub[isub].rep[t].beta[i];
            }  
            double *tmp = (double *)calloc(SB2,sizeof(double));
            Ax(tmp,pop->re_rep_prec,beta,(const int)SB2,(const int)SB2);
            double *At = (double *)calloc(SB2*N,sizeof(double));
            transpose(At,sub[isub].X,(const int)SB2,(const int)N);
            Ax(tmpM,At,tmp,(const int)N,(const int)SB2);
//            printf("x\n");fflush(NULL);
            free(tmp);
//            printf("y\n");fflush(NULL);
            double *tmp2 = (double *)calloc(SB2*N,sizeof(double));
//            printf("Z\n");fflush(NULL);
            AxB(tmp2,pop->re_rep_prec,sub[isub].X,(const int)SB2,(const int)SB2,(const int)N);
            for (int i=0;i<N*N;i++)
                tmpP[i] = 0;
            AxB(tmpP,At,tmp2,(const int)N,(const int)SB2,(const int)N);
//            printf("xx\n");fflush(NULL);
            free(tmp2);
//            printf("yy\n");fflush(NULL);
            free(At);
//            printf("zz\n");fflush(NULL);
            free(beta);
 //           printf("xyz\n");fflush(NULL);
 //           calApVinvX(tmpM,sub[isub].X,sub[isub].rep[0].beta,pop->re_rep_prec,(const int)SB,(const int)N);
//            calApVinvA(tmpP,sub[isub].X,pop->re_rep_prec,(const int)SB,(const int)N);
        }
        else {
            calApVinvX(tmpM,sub[isub].X,sub[isub].beta,pop->re_prec,(const int)SB,(const int)N);
            calApVinvA(tmpP,sub[isub].X,pop->re_prec,(const int)SB,(const int)N);
        }
        
        for (int i=0;i<N;i++)
            M[i] += tmpM[i];

        for (int i=0;i<N*N;i++)
            P[i] += tmpP[i];
    }
    
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

    free(P);
    free(M);
    free(tmpM);
    free(tmpP);
}

