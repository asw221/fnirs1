#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fNIRS.h"

extern FILE *flog;

double dens_gamma(double x,double alpha,double beta) {
    double ll;
    
    ll = alpha*log(beta) - lgamma(alpha) + (alpha-1)*log(x) - beta*x;
        
    return ll;
}

double tden(double x,double mean,double var,double df) {
    double tmp,value = 0;
    double pi = 3.14159;
    
    tmp = (x-mean)*(x-mean)/var;
    value = lgamma(0.5*(df+1)) - lgamma(0.5*df) - 0.5*log(df*pi*var) - 0.5*(df+1)*log1p(tmp/df);
    return value;
}

int compdbl(const void *c1,const void *c2) // compares two doubles, used in qsort
{
 double a,b;
 a = *(double *)c1;
 b = *(double *)c2;
 if (a < b) return -1;
 return 1;
}

void compute_mean_sd(double *mean,double *sd,const double *x,const int len) {
    double N;
    double tmp;
    
    *mean = 0;
    *sd = 0;
    
    for (int i=0;i<len;i++) {
        tmp = x[i];
        *mean += tmp;
        *sd += tmp*tmp;
    }
    N = (double)len;
    *mean /= N;
    *sd = sqrt((N/(N-1))*(*sd/N - *mean* *mean));
}

void standardize_data(double *x,const int N) {
    double mean,sd;
    
    compute_mean_sd(&mean,&sd,(const double *)x,N);
    for (int i=0;i<N;i++)
//        x[i] = (x[i] - mean)/0.000001986;
        x[i] = (x[i] - mean)/sd;
}

void compute_quantile(double *quantile,const double prob,const double *x,const int len) {
    //expects x to be sorted
    //this method is the same as the default in R version 3.4.3
    
    int fh;
    double h,N;
              
    N = (double)len;
      
    h = (N - 1)*prob + 1;
    fh = (int)floor(h);

    *quantile = x[fh-1] + (h - fh)*(x[fh] - x[fh-1]);  
}

void compute_credible_interval(double *interval,const double prob,double *x,const int len) {
    double plow,phigh;
  
    qsort(x,(size_t)len,sizeof(double *),compdbl);
    
    plow = (1 - prob)/2;
    compute_quantile(&interval[0],(const double)plow,(const double *)x,len);
    phigh = 1 - (1-prob)/2;
    compute_quantile(&interval[1],(const double)phigh,(const double *)x,len);
}

void compute_stats(POP *pop,const double cred_int,const int Niter) {
    int N;
    double *X,*Xt,mean,sd;
    double interval[2];
    FILE *fparmest,*fout;

    fparmest = fopen("./log/Parameter_Estimates.log","w");
 
    if (pop->GRP) {
        fout = fopen("./log/pop_beta.log","r");

        N = 0;
        X = (double *)calloc(Niter*pop->Ncov,sizeof(double));
        while (fscanf(fout,"%lf ",&X[N]) != EOF)
            N++;
        fclose(fout);
        
        Xt = (double *)calloc(Niter*pop->Ncov,sizeof(double));
        for (int i=0;i<Niter;i++)
            for (int j=0;j<pop->Ncov;j++)
                Xt[j*Niter + i]  = X[i*pop->Ncov + j];
        free(X);
       
 //       N /= Niter*pop->Ns;
 //       N = pop->Nb*pop->Nc;
 //       printf("N = %d, Niter = %d pop->Ns = %d\n",N,Niter,pop->Ns);
        for (int i=0;i<pop->Ncov;i++) {
            compute_mean_sd(&mean,&sd,(const double *)&Xt[i*Niter],(const int)Niter);
            compute_credible_interval(interval,cred_int,&Xt[i*Niter],(const int)Niter);
            if (i == 0)
                fprintf(fparmest,"Pop Level Parameters: \n");
            fprintf(fparmest,"\t%30s mean = %.3lf\tsd = %.3lf\n",pop->covnames[i],mean,sd);
            fprintf(fparmest,"\t%30s 95%% Cred.Int. = (%.3lf, %.3lf)*\n\n"," ",interval[0],interval[1]);
        }
        free(Xt);
    }

    fprintf(fparmest,"\n\n\nParameter Summary for Subjects\n");

    char *CC = (char *)calloc(300,sizeof(char));
    char *S = (char *)calloc(400,sizeof(char));
    
    X = (double *)calloc((Niter)*pop->Nb*pop->Ns,sizeof(double));
    Xt = (double *)calloc((Niter)*pop->Nb*pop->Ns,sizeof(double));
    for (int isub=0;isub<pop->N_SUBS;isub++) {
        for (int irep=0;irep<pop->sub[isub].N_REPS;irep++) {
            
            CC = strcpy(CC,pop->sub[isub].rep[irep].dataname);
            CC = strtok(CC,".");
            S = strcpy(S,"./log/");
            S = strcat(S,CC);
            S = strcat(S,"_beta.log");
            fout = fopen(S,"r");
//            printf("%s\n",S);
            N = 0;
            while (fscanf(fout,"%lf ",&X[N]) != EOF)
                N++;
            fclose(fout);
            for (int i=0;i<(Niter);i++)
                for (int j=0;j<pop->Nb*pop->Ns;j++)
                    Xt[j*Niter + i]  = X[i*pop->Nb*pop->Ns + j];

//    NEED TO FIX THIS
            fprintf(fparmest,"\n\n\t%s\n",CC);
            int ii = 0;
            for (int istim = 0;istim<pop->Ns;istim++) {
                fprintf(fparmest,"\n\t\tCond = %d\n\n",istim);
                
                for (int j=0;j<pop->Nb;j++) {
                    ii = Niter*(j*pop->Ns + istim);
                    compute_mean_sd(&mean,&sd,(const double *)&Xt[ii],(const int)Niter);
                    compute_credible_interval(interval,cred_int,&Xt[ii],(const int)Niter);
 //                   ii += Niter;
                    
                    if (j==0)
                        fprintf(fparmest,"\n\t\t\tHRF\n\n");
                    else if (j==1)
                        fprintf(fparmest,"\n\t\t\tHRF Temporal Derivative\n\n");
                    else
                        fprintf(fparmest,"\n\t\t\tHRF Dispersion Derivative\n\n");
           
                    fprintf(fparmest,"\t\t\t\tmean = %.3lf\tsd = %.3lf\n",mean,sd);
  
                    fprintf(fparmest,"\t\t\t\t95%% Cred.Int. = (%.3lf, %.3lf)*\n",interval[0],interval[1]);
                }
            }
        }
        if (!pop->No_replicates) {
            S = strcpy(S,"./log/sub_");
            S = strcat(S,CC);
            S = strcat(S,"_beta.log");
            fout = fopen(S,"r");

            N = 0;
            while (fscanf(fout,"%lf ",&X[N]) != EOF)
                N++;
            fclose(fout);
            for (int i=0;i<Niter;i++)
                for (int j=0;j<pop->Nb*pop->Ns;j++)
                    Xt[j*Niter + i]  = X[i*pop->Nb*pop->Ns + j];

            N /= Niter*pop->Ns;
            
            fprintf(fparmest,"\n\n\t%s, Subject Level Results\n",CC);
            for (int istim = 0;istim<pop->Ns;istim++) {
                fprintf(fparmest,"\n\t\tCond. = %d\n\n",istim);
                
                for (int j=0;j<N;j++) {
                    if (!j)
                        fprintf(fparmest,"\n\t\t\tHRF\n\n");
                    else
                        fprintf(fparmest,"\n\t\t\tHRF Temporal Derivative\n\n");
                    
                   compute_mean_sd(&mean,&sd,(const double *)&Xt[(istim*N+j)*Niter],(const int)Niter);
                   compute_credible_interval(interval,cred_int,&Xt[(istim*N+j)*Niter],(const int)Niter);
 
                   fprintf(fparmest,"\t\t\t\tHb0:\tmean = %.3lf\tsd = %.3lf\n",mean,sd);

                   fprintf(fparmest,"\t\t\t\t    \t95%% Cred.Int. = (%.3lf, %.3lf)*\n\n",interval[0],interval[1]);
                }   
            }
        }
    }
    fprintf(fparmest,"\n\n95%% Cred.Int. based on the 0.025 and 0.975 quantiles.\nQuantiles estimated using the default method in R version 3.4.3.\n");
    fclose(fparmest);
  
    free(Xt);
    free(X);
    free(CC);
    free(S);
}

void center_covars(double *X,const int nrow,const int ncol) 
{
    if (ncol > 1) {
        for (int j=0;j<ncol;j++) {
            double mean = 0;
            int cnt = 0;

            for (int i=0;i<nrow;i++) {
                mean += X[i*ncol+j];
                if ((X[i*ncol+j] == 1) || (X[i*ncol+j] == 0))
                    cnt++;                     
            }    
            
            if (cnt != nrow) {
                mean /= (double)(nrow);
            
                for (int i=0;i<nrow;i++)
                    X[i*ncol+j] -= mean;
            }
        }
    }
}



