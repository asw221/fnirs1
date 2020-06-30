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
#include <omp.h>
#include "randgen.h"
#include "cholesky.h"
#include <time.h>
#include "fNIRS.h"

#include <iostream>

extern int MAX_ITER;
extern int BURN_IN;
extern int maxP;
extern int nthd_sys;
extern FILE *flog;
extern int PPP;

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
    int nrow_loglik;
    SUB *sub;
    REP *rep;
    FILE *fout;
    clock_t start,end;
    void draw_knot_locations(REP *rep,int sdegree,int *flag,unsigned long *seed);
    int knot_birth_death(REP *rep,POP *pop,const int sdegree,int iter,unsigned long *seed);
 
//    void mcmc_delta1(REP *rep,int iter,dlm_tst_flag flag,unsigned long *seed);
//    void mcmc_delta2(REP *rep,int iter,dlm_tst_flag flag,unsigned long *seed);
    void DLM(REP *rep,int iter,unsigned long *seed);
    void draw_precYstart(REP *rep,unsigned long *seed);
    void draw_preceta(REP *rep,unsigned long *seed);
    void draw_reprec(POP *pop,SUB *sub,unsigned long *seed);
    void draw_re_rep_prec(POP *pop,SUB *sub,unsigned long *seed);
    void draw_dY(REP *rep,unsigned long *seed);
    void draw_prob(REP *rep,unsigned long *seed);
    void draw_precY(REP *rep,unsigned long *seed);
    void draw_hrf_mean(REP *rep,unsigned long *seed);
 
/*    void calculate_residuals(REP *rep,int P);
    void calculate_res1(REP *rep,int P);
    void calculate_res2(REP *rep,int P);
    void calculate_res3(REP *rep);f
    void calculate_res4(REP *rep,int P);*/
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    
    void draw_beta_eta(POP *pop,SUB *sub,REP *rep,unsigned long *seed);
//    void draw_beta_eta(POP *pop,SUB *sub,REP *rep,unsigned long *seed,int iter,int isub);
    void draw_pop_beta(POP *pop,SUB *sub,unsigned long *seed);
    void draw_sub_beta(POP *pop,SUB *sub,unsigned long *seed);    

    void add_basis(POP *pop,REP *rep,int sdegree,unsigned long *seed);
    void remove_basis(POP *pop,REP *rep,int sdegree,unsigned long *seed);
    
    void DIC(POP *pop,REP *rep,unsigned long *seed);
   
    pop->ED = 0;
       
    int sdegree = 4;
    int aaa = 0;

    int NN = 0;//pop->N_SUBS*pop->sub[0].N_REPS;
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
 
    int nthd = (int)ceil(
			 (pop->N_SUBS * pop->No_replicates) /
			 (int)ceil((pop->N_SUBS * pop->No_replicates) / (double)nthd_sys)
			 );
    if (pop->N_SUBS == 1)
        nthd = 1;

    // std::cout << "nthd = " << nthd << std::endl;
    // std::cout << "nthd_sys = " << nthd_sys << std::endl;
    
    omp_set_num_threads(nthd);

    int *v = (int *)malloc(NN*sizeof(int));

//     // Andrew (4):
//     int imax, flag = 0;
// #pragma omp parallel prifate(j, k, imax, flag)
//     {
// #pragma omp barrier
    for (iter=0;iter<=MAX_ITER;iter++) {
        if (!(iter%100)) { 
            printf("%d",iter);
            fflush(stdout);
        }
        if (!(iter%20) && (iter%100)) {
            printf(".");
            fflush(stdout);
        }
       

	// Andrew, comment (+2):
        // int imax,flag = 0;
	imax = 0;
	flag = 0;
        for (int i=0;i<NN;i++)
            v[i] = i;
        permute_sample_int(v,NN,seed);
      
#pragma omp parallel for private(j,k,imax,flag)
        for (i=0;i<NN;i++) {
            if (iter == 0) {

	      std::cout << "i = " << i << ",  " << "iThread = "
			<< ( ::omp_get_thread_num() )
			<< std::endl;
	      
//                draw_beta_eta(pop,tmpsub[subidx[v[i]]],tmprep[v[i]],seed);
                for (j=0;j<50;j++) {//printf("j = %d v[i] = %d\n",j,v[i]);
                    tmprep[v[i]]->P = 0;
//                  start = clock();
                    aaa = knot_birth_death(tmprep[v[i]],pop,sdegree,iter,seed);
//                  end = clock();
//                  if (!(iter%100)) printf("birth-death time = \t%lf\t tries = %d\n",(double)(end-start)/CLOCKS_PER_SEC,aaa);
            
//                  #pragma omp parallel for private(i,j,imax)
//                  for (j=0;j<NN;j++) {
//                    start = clock();
                    
                    imax = tmprep[v[i]]->nKnots-sdegree*2;
                    for (k=0;k<imax;k++) 
                        draw_knot_locations(tmprep[v[i]],sdegree,&flag,seed);
                    
                    draw_beta_eta(pop,tmpsub[subidx[v[i]]],tmprep[v[i]],seed);
                    draw_precY(tmprep[v[i]],seed);
                    draw_preceta(tmprep[v[i]],seed);
//                  end = clock();
//                if (!(iter%100)) printf("draw knots time = \t%lf\t knots = %d\n",(double)(end-start)/CLOCKS_PER_SEC,imax);
                   tmprep[v[i]]->P = PPP ;
                }
            }
            else {
//                if (iter < BURN_IN) {

          
//                  start = clock();
                    aaa = knot_birth_death(tmprep[v[i]],pop,sdegree,iter,seed);
//                  end = clock();
//                  if (!(iter%100)) printf("birth-death time = \t%lf\t tries = %d\n",(double)(end-start)/CLOCKS_PER_SEC,aaa);
            
//                  #pragma omp parallel for private(i,j,imax)
//                  for (j=0;j<NN;j++) {
//                    start = clock();
          
                    imax = tmprep[v[i]]->nKnots-sdegree*2;
                    for (k=0;k<imax;k++) 
                    draw_knot_locations(tmprep[v[i]],sdegree,&flag,seed);
                    
//                  end = clock();
//                if (!(iter%100)) printf("draw knots time = \t%lf\t knots = %d\n",(double)(end-start)/CLOCKS_PER_SEC,imax);
//                }
            }     
 //        #pragma omp parallel for private(i)
//        for (i=0;i<NN;i++) {
            start = clock();
            draw_beta_eta(pop,tmpsub[subidx[v[i]]],tmprep[v[i]],seed);
            draw_preceta(tmprep[v[i]],seed);
            draw_precYstart(tmprep[v[i]],seed);
            end = clock();
 //           if (!(iter%100)) printf("other time = \t\t%lf\n",(double)(end-start)/CLOCKS_PER_SEC);
 
            start = clock();
            DLM(tmprep[v[i]],iter,seed);
            end = clock();
 //           if (!(iter%100)) printf("DLM time = \t\t%lf\t iter = %d\n\n",(double)(end-start)/CLOCKS_PER_SEC,iter);
       
//            draw_dY(tmprep[i],seed);
//            draw_precY(tmprep[i],seed);        
        }
	// #pragma omp barrier  // -- Andrew
         
        if (pop->GRP) {
            draw_re_rep_prec(pop,pop->sub,seed);
            
            int isub2;
            if (!pop->No_replicates) {
                #pragma omp parallel for
                for (isub2=0;isub2<pop->N_SUBS;isub2++)
                    draw_sub_beta(pop,tmpsub[isub2],seed);
		// #pragma omp barrier  // -- Andrew
 
                draw_reprec(pop,pop->sub,seed);
            } 
            draw_pop_beta(pop,pop->sub,seed);
        }
        
        if (!(iter%50) && (iter > 0)) {
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
                        rep->accept[1] = rep->attempt[1] = 0;
                        if (!(iter%100))
                           fprintf(flog,"\t\t rt = %lf proposal = %lf\n",rt,pop->sub[isub].rep[irep].prop_sd[1]);
                        
                        rt = (double)rep->accept[3]/(double)rep->attempt[3];
                        adjust_acceptance2(rt,&(rep->prop_sd[3]),0.35);
                        rep->accept[3] = rep->attempt[3] = 0;
                        if (!(iter%100))
                            fprintf(flog,"\t\t rt = %lf proposal = %lf\n",rt,rep->prop_sd[3]);

//                        rt = (double)rep->accept[2]/(double)rep->attempt[2];
 //                       fprintf(flog,"rt = %lf %d\n",rt,rep->accept[2]);
//                        adjust_acceptance2(rt,&(rep->prop_sd[2]),0.35);
//                        rep->accept[2] = rep->attempt[2] = 0;
//                            if (!(iter%100))
//                                fprintf(flog,"\t\t rt = %lf proposal = %lf mean = %lf\n\n",rt,rep->prop_sd[2],rep->mhrf);
                      

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
                    fprintf(flog,"\t precY0 = %10.6lf \t preceta = %10.6lf\n",rep->precYstart,rep->preceta);
                    fprintf(flog,"\t df_delta1 = %10.6lf \t df_delta2 = %10.6lf \t P = %3d\n",rep->df_delta1,rep->df_delta2,rep->P);
 //                   fprintf(flog,"\t cnt = %d, pi = %lf\n",rep->dim_X[0]-rep->d_Ycnt,rep->pi);
                    
                }
            }
            
            for (k=0;k<pop->Ns;k++) {
                fprintf(flog,"Condition %d\n",k);
                if (pop->GRP) {
                    for (int i=0;i<pop->Nb;i++)
                        fprintf(flog,"\t %g \t %g",pop->re_prec[k*pop->Nb+i],sqrt(1/pop->re_rep_prec[k*pop->Nb+i]));
                    fprintf(flog,"\n\n");
                }
                for (i=0;i<pop->N_SUBS;i++) {
                    sub = &(pop->sub[i]);
                    for (int irep=0;irep<sub->N_REPS;irep++) {
                        for (j=0;j<pop->Nb;j++) {
                            fprintf(flog,"\tsub %d:\t %10.6lf ",i,sub->rep[irep].beta[k+j*pop->Ns]);
                        }
                       // fprintf(flog,"\n");
                    }
                    fprintf(flog,"\n");
                }
                fprintf(flog,"\n");
                if (pop->GRP) {
                    if (!pop->No_replicates) {
                        double meansub = 0;
                        for (j=0;j<pop->N_SUBS;j++) {
                            for (i=0;i<pop->Nb;i++) {
                                fprintf(flog,"\t %10.6lf ",pop->sub[j].beta[k+i*pop->Ns]);
                                if (i == 0)
                                    meansub += pop->sub[j].beta[k+i*pop->Ns];
                             }
                             fprintf(flog,"\n");
                        }
                        meansub /= pop->N_SUBS;
                        fprintf(flog,"\n\n");
                        fprintf(flog,"\t %10.6lf \n",meansub);
                    }
//                    for (i=0;i<pop->Ncov*pop->Nb;i++)
//                        fprintf(flog,"\t %10.6lf \t ",pop->beta[k*pop->Ncov*pop->Nb+i]);
//                    for (i=0;i<pop->Nc;i++)
//                        fprintf(flog,"\t %10.6lf \t ",pop->beta[k*pop->Nc+i]);
//                     fprintf(flog,"\n\n");
                }
            }
            if (pop->GRP) {
                fprintf(flog,"Pop Parms:\n");
                for (i=0;i<pop->Ncov;i++)
                    fprintf(flog,"%30s = %10.6lf\n",pop->covnames[i],pop->beta[i]);
                fprintf(flog,"\n\n");
            }
        }
  //      fflush(flog);
        
        if ((iter>BURN_IN)) {
            for (isub=0;isub<pop->N_SUBS;isub++) {
                sub = &(pop->sub[isub]);
                
                for (irep=0;irep<pop->sub[isub].N_REPS;irep++) {
                    rep = &(sub->rep[irep]);

                    for (i=rep->P;i<rep->dim_W[0];i++) {
                        for (j=0;j<rep->dim_W[1];j++) {
                            rep->mdelta[i*maxP+j] += rep->delta[i*rep->dim_W[1]+j];
                            rep->mdelta2[i*rep->dim_W[1]+j] += rep->delta[i*rep->dim_W[1]+j]*rep->delta[i*rep->dim_W[1]+j];
                        }
                        for (j=rep->dim_W[1];j<maxP;j++)
                            rep->mdelta[i*maxP+j] = 0;
                    }
 /*                    FILE *ftmp = fopen("veta.dat","w");
                    for (i=0;i<rep->dim_V[0];i++)
                        fprintf(ftmp,"%lf ",rep->Veta[i]);
                    fprintf(ftmp,"\n");
                    fclose(ftmp);
                    ftmp = fopen("Wdelta.dat","w");
                    for (i=0;i<rep->dim_W[0];i++)
                        fprintf(ftmp,"%lf ",rep->Wdelta[i]);
                    fprintf(ftmp,"\n");
                    fclose(ftmp);                   
                    for (i=0;i<rep->dim_X[1];i++)
                        printf("%lf\n",sub->beta[i]);
                    exit(0);*/

                    for (i=0;i<rep->dim_X[0];i++) {
  //                      rep->md_Y[i] += 1./sqrt(rep->d_Y[i]);
                        rep->sd_Y[i] += rep->d_Y[i];
                    }
                    
                    for (i=0;i<rep->dim_V[1];i++)
                        fprintf(rep->fout_eta,"%lf ",rep->eta[i]);
                    fprintf(rep->fout_eta,"\n");
                   
 //                   DIC(pop,rep,seed);
 
                    for (i=0;i<rep->dim_X[0];i++) {
                        rep->mean_fit[i] += rep->Y[i] - rep->residuals[i];
   //                     rep->std_res[i]  += rep->residuals[i]*sqrt(rep->d_Y[i]);
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
                    
                    fprintf(rep->fout_nknots,"%d ",rep->nKnots-2*sdegree);

                    for (i=sdegree;i<rep->nKnots-sdegree;i++)
                        fprintf(rep->fout_knots,"%lf ",rep->knots[i]);
                    fprintf(rep->fout_knots,"\n");  

                    for (i=0;i<pop->Nb*(pop->Ns);i++)
                        fprintf(rep->fout_beta,"%lf ",rep->beta[i]);
                    fprintf(rep->fout_beta,"\n");  

                }
                
                for (i=0;i<pop->Nb*(pop->Ns);i++)
                    fprintf(sub->fout_beta,"%lf ",sub->beta[i]);
                fprintf(sub->fout_beta,"\n");
            }

            for (i=0;i<pop->Ncov;i++)
                fprintf(pop->fout_beta,"%lf ",pop->beta[i]);
            fprintf(pop->fout_beta,"\n");
            
            for (i=0;i<pop->Nb*pop->Ns;i++)
                fprintf(pop->fout_reprec,"%lf ",pop->re_prec[i]);
            fprintf(pop->fout_reprec,"\n");

            for (i=0;i<pop->Nb*pop->Ns;i++)
                fprintf(pop->fout_re_rep_prec,"%lf ",pop->re_rep_prec[i]);
            fprintf(pop->fout_re_rep_prec,"\n");


            fflush(NULL);      
        }
    }
    }  // -- Andrew (close omp parallel section)
    free(v);
    free(tmprep);
    free(tmpsub);
    free(subidx);
    
    fclose(pop->fout_beta);
    fclose(pop->fout_reprec);
    fclose(pop->fout_re_rep_prec);

    for (isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
            fclose(rep->fout_dlm);
            fclose(rep->fout_eta);
            fclose(rep->fout_beta);
            fclose(rep->fout_nknots);
            fclose(rep->fout_knots);
        }
        fclose(sub->fout_beta);
    }

    for (isub=0;isub<pop->N_SUBS;isub++) {
        sub = &(pop->sub[isub]);
        for (irep=0;irep<sub->N_REPS;irep++) {
            rep = &(sub->rep[irep]);
            for (i=0;i<rep->dim_W[0];i++) {
                for (j=0;j<maxP;j++) {
                    rep->mdelta[i*maxP+j] /= (double)(MAX_ITER-BURN_IN);
                    fprintf(rep->fout_delta,"%lf ",rep->mdelta[i*maxP+j]);
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
 //                   rep->md_Y[i] /= (double)(MAX_ITER-BURN_IN);
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
                rep->std_res[i] = rep->mean_res[i]/rep->sd_Y[i];
            }
       }
    }
/*    double DE=0;
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
*/    
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
    
    b = BETA + 0.5*b;
    
    rep->precYstart = rgamma(a,b,seed);
    for (int i=0;i<rep->P;i++)
        rep->d_Y[i] = rep->precYstart;

}

void draw_dY(REP *rep,unsigned long *seed) {
    double ratio,prob;
    
    void zerooutdouble(double *A,const int dim);
    void Ax(double *b,double *A,double *x,const int nrow,const int ncol);
    double dotProd(double *x,double *y,const int dim);
    
    prob = rep->pi;
    int P = rep->P;
    int N = rep->dim_X[0];
    double phi;
    
    rep->d_Ycnt = 0;
    for (int t=0;t<N;t++) {
        if (rep->dY[t]) {
            phi = rep->d_Y[t];
            break;
        }
    }
    for (int t=P;t<N;t++) {
        ratio = rep->residuals[t]*rep->residuals[t];//*rep->d_Y[t];
        ratio = (rep->LY*(1.-prob)/prob) *exp(0.5*phi*(1.-rep->LY)*ratio);
        ratio = -log1p(ratio);
        if (log(kiss(seed)) < ratio) {
            rep->dY[t] = 1;
            (rep->d_Ycnt)++;
        }
        else
            rep->dY[t] = 0;
    }
}

void draw_prob(REP *rep,unsigned long *seed) {
    double ALPHA = 95;
    double BETA =  5;
    
    int N = rep->dim_X[0];
    int P = rep->P;
    int cnt = rep->d_Ycnt;
    rep->pi = rbeta(ALPHA+cnt,BETA+N-P-cnt,seed);
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
            
/*    for (i=0;i<ncol;i++) 
        for (j=0;j<nrow;j++)
            ApX[i] += At[i*nrow+j]*Vinv[j]*X[j]; 
 */           
    free(At); 
    free(Vx);
}

void calApVinvA(double *ApA,double *A,double *Vinv,const int nrow,const int ncol) {
// only for diagonal matrices, V
    int i,j,k;
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
/*     for (k=0;k<nrow;k++) {
        double V = Vinv[k];
        double *tmp = &A[k*ncol];
        for (i=0;i<ncol;i++) {
            double *tmpout = &ApA[i*ncol];
            for (j=0;j<ncol;j++)
                tmpout[j] += V*tmp[i]*tmp[j];
        }
    }*/
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
    if (tmp < 10)
        rep->preceta = tmp;
}


void draw_reprec(POP *pop,SUB *sub,unsigned long *seed) {
    int N,M;
    double ALPHA, BETA,*rss;
    double *Xb;
    ALPHA = 1;//3;
    BETA  = 1;//2;
    
//    M = pop->Ns*pop->Nb;
//    N = pop->Ncov*M;
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
    double ALPHA, BETA,rss,*Xb;
    REP *rep;
    
    M = pop->Nb*pop->Ns;
    N = pop->Ncov;
    Xb = (double *)calloc(M,sizeof(double));
    ALPHA = 1;//-0.5;//3;
    BETA  = 1;//0;//2;
        
    for (int is=0;is<pop->Ns;is++) {
        for (int j=0;j<pop->Nb;j++) {
            rss = 0;
            int cnt = 0;
            for (int i=0;i<pop->N_SUBS;i++) {
                if (pop->No_replicates) {
                    calAx(Xb,sub[i].X,pop->beta,(const int)M,(const int)N);
                    rep = &(sub[i].rep[0]);
                    rss += (rep->beta[j+is*pop->Nb]-Xb[j+is*pop->Nb])*(rep->beta[j+is*pop->Nb]-Xb[j+is*pop->Nb]);
                    cnt++;
                }
                else {
                    for (int k=0;k<sub[i].N_REPS;k++) {
                        rep = &(sub[i].rep[k]);
                        rss += (rep->beta[j+is*pop->Nb]-sub[i].beta[j+is*pop->Nb])*(rep->beta[j+is*pop->Nb]-sub[i].beta[j+is*pop->Nb]);
                        cnt++;
                    }
                }
            }
            double tmp = rgamma(ALPHA + 0.5*(double)cnt,BETA + 0.5*rss,seed);
//            printf("cnt = %d, rss = %lf, tmp = %lf\n",cnt,rss,tmp);fflush(NULL);
            if (tmp < 1000.)
                pop->re_rep_prec[is*pop->Nb+j] = tmp;
        }
    }
    free(Xb);
}

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
        betaeta[i] = rep->beta[i];
    for (i=rep->dim_X[1];i<ncol;i++)
        betaeta[i] = rep->eta[i-rep->dim_X[1]];
 
     calH(J,X,rep->delta,(const int)rep->dim_V[0],(const int)ncol,(const int)rep->P);
  
    /* calculate J(V,delta)eta */
    double *Jbe;
    Jbe = (double *)calloc(rep->dim_V[0],sizeof(double));
    calAx(Jbe,J,betaeta,(const int)rep->dim_V[0],(const int)ncol);
    
    free(betaeta);
    // calculate Y - Wdelta - J(V,delta)eta
    
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
        if (pop->No_replicates) {
            int Nrow = pop->Nb*pop->Ns;
            int Ncol = pop->Ncov;
            double *Xb = (double *)calloc(Nrow,sizeof(double));
            calAx(Xb,sub->X,pop->beta,(const int)Nrow,(const int)Ncol);
            for (int is=0;is<pop->Ns;is++) {
                for (int i=0;i<pop->Nb;i++) {
                    M[i+is*pop->Nb] += pop->re_rep_prec[is*pop->Nb+i]*Xb[i+is*pop->Nb];
                    VpV[(i+is*pop->Nb)*ncol+i+is*pop->Nb] += pop->re_rep_prec[is*pop->Nb+i];
                }
            }
            free(Xb);
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

    for (i=0;i<rep->dim_X[1];i++)
        rep->beta[i] = mean[i];
    
    for (i=rep->dim_X[1];i<ncol;i++)
        rep->eta[i-rep->dim_X[1]] = mean[i];
        
    /* calculate Veta and Xbeta */
    
    calAx(rep->Veta,rep->V,rep->eta,(const int)rep->dim_V[0],(const int)rep->dim_V[1]);
    calAx(rep->Xbeta,rep->X,rep->beta,(const int)rep->dim_X[0],(const int)rep->dim_X[1]);

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
 
    M = (double *)calloc(N,sizeof(double));
    tmpM = (double *)calloc(N,sizeof(double));
    P = (double *)calloc(N*N,sizeof(double));
    tmpP = (double *)calloc(N*N,sizeof(double));
     
    for (int isub=0;isub<pop->N_SUBS;isub++) {
        if (pop->No_replicates) { 
            calApVinvX(tmpM,sub[isub].X,sub[isub].rep[0].beta,pop->re_rep_prec,(const int)SB,(const int)N);
            calApVinvA(tmpP,sub[isub].X,pop->re_rep_prec,(const int)SB,(const int)N);
        }
        else {
            calApVinvX(tmpM,sub[isub].X,sub[isub].beta,pop->re_prec,(const int)SB,(const int)N);
            calApVinvA(tmpP,sub[isub].X,pop->re_prec,(const int)SB,(const int)N);
        }
        
        for (int i=0;i<N;i++)
            M[i] += tmpM[i];

        for (int i=0;i<N;i++)
            for (int j=0;j<N;j++)
                P[i*N+j] += tmpP[i*N+j];
    }
    
//    for (int i=0;i<N;i++)
//        P[i*N+i] += 0.0001;//0.0001;

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
