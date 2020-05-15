/*
 *  hrf.cpp
 *  fNIRS
 *
 *  Created by Timothy Johnson on 11/08/15.
 *  Copyright 2015 University of Michigan. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/usr/local/include/fftw3.h"
#include "fNIRS.h"
#include "randgen.h"

double *DCT_basis(int N,double T,int period,int *K)
{
    int k;
    double *X;
    double pi = 3.14159265359;
    
    *K = (int)floor(2*T/(double)period);
    X = (double *)calloc(N*(*K+1),sizeof(double));
//    for (int i=0;i<N;i++)
//        X[i] = (double *)calloc(*K+1,sizeof(double));

    for (int i=0;i<N;i++) {
        for (k=0;k<=*K;k++) {
            X[i*(*K+1)+k] = cos(pi/(double)N *((double)i + 0.5)*(double)k);
        }
    }
    return X;
}

//HRF = function(T,a1=6,b1=1,a2=16,b2=1,c=1/6) {
//    HRF = b1^a1/gamma(a1) *T^(a1-1) * exp(-b1*T) - c*b2^a2/gamma(a2) *T^(a2-1) * exp(-b2*T)
//    return(HRF)
//}

void proj(double *out,double *v,double *u,int length)
{
    double num=0,den=0;
    
    for (int i=0;i<length;i++) {
        num += v[i]*u[i];
        den += u[i]*u[i];
    }
    num = num/den;
    for (int i=0;i<length;i++)
        out[i] = num*u[i];
}

double **canonical_HRF(int T,double freq,int *dim_HRF,double m,double var,int type)
{
    int num;
    double a1 = m*m/var;
    double b1 = m/var;
    double a2 = 16;
    double b2 = 1;
    double c = 1./6.;
    double x,**HRF;
    double const1,const2;
    
    const1 = pow(b1,a1)/tgamma(a1);
    const2 = c*pow(b2,a2)/tgamma(a2);
    
    num = T*freq;
    dim_HRF[0] = num;
    switch (type) {
        case 0:
            dim_HRF[1] = 1;
            HRF = (double **)calloc(num,sizeof(double *));
            for (int i=0;i<num;i++)
                HRF[i] = (double *)calloc(dim_HRF[1],sizeof(double));
            break;
        case 1:
            dim_HRF[1] = 2;
            HRF = (double **)calloc(num,sizeof(double *));
            for (int i=0;i<num;i++)
                HRF[i] =(double *)calloc(dim_HRF[1],sizeof(double));
            break;
        case 2: default:
            dim_HRF[1] = 3;
            HRF = (double **)calloc(num,sizeof(double *));
            for (int i=0;i<num;i++)
                HRF[i] =(double *)calloc(dim_HRF[1],sizeof(double));
            break;
    }
    for (int i=0;i<num;i++) {
        x = (double)i/(double)freq;
        switch (type) {
            case 0:
                HRF[i][0] = const1*pow(x,a1-1)*exp(-b1*x) - const2*pow(x,a2-1)*exp(-b2*x);
                break;
            case 1:
                HRF[i][0] = const1*pow(x,a1-1)*exp(-b1*x) - const2*pow(x,a2-1)*exp(-b2*x);
                HRF[i][1] = const1*pow(x,a1-1)*exp(-b1*x)*( (a1-1)/x - b1  ) - const2*pow(b2,a2)*pow(x,a2-1)*exp(-b2*x)*((a2-1)/x - b2);
                break;
            case 2: default:
                HRF[i][0] = const1*pow(x,a1-1)*exp(-b1*x) - const2*pow(x,a2-1)*exp(-b2*x);
                HRF[i][1] = const1*pow(x,a1-1)*exp(-b1*x)*( (a1-1)/x - b1  ) - const2*pow(b2,a2)*pow(x,a2-1)*exp(-b2*x)*((a2-1)/x - b2);
                HRF[i][2] = pow(x,a1-1)/tgamma(a1)*exp(-b1*x)*(a1*pow(b1,a1-1) - x*pow(b1,a1) );
                break;
        }
    }
//    HRF[0][1] = 0;
    double *out = (double *)calloc(num,sizeof(double));
    double *v = (double *)calloc(num,sizeof(double));
    double *u = (double *)calloc(num,sizeof(double));
    
    switch (type) {
        case 0:
            double tmp;
            for (int i=0;i<num;i++)
                tmp += HRF[i][0]*HRF[i][0];
            tmp = sqrt(tmp);
            for (int i=0;i<num;i++)
                HRF[i][0] /= tmp;
            break;
        case 1:
            for (int i=0;i<num;i++) {
                u[i] = HRF[i][0];
                v[i] = HRF[i][1];
            }
            proj(out,v,u,num);
            for (int i=0;i<num;i++) {
                HRF[i][1] = v[i] - out[i];
                v[i] = HRF[i][2];
            }
           
            
            tmp = 0;
            for (int i=0;i<num;i++)
                tmp += HRF[i][0]*HRF[i][0];
            tmp = sqrt(tmp);
            for (int i=0;i<num;i++)
                HRF[i][0] /= tmp;
            tmp = 0;
            for (int i=0;i<num;i++)
                tmp += HRF[i][1]*HRF[i][1];
            tmp = sqrt(tmp);
            for (int i=0;i<num;i++)
                HRF[i][1] /= tmp;
             break;
        case 2: default:
            for (int i=0;i<num;i++) {
                u[i] = HRF[i][0];
                v[i] = HRF[i][1];
            }
            proj(out,v,u,num);
            for (int i=0;i<num;i++) {
                HRF[i][1] = v[i] - out[i];
                v[i] = HRF[i][2];
            }
            proj(out,v,u,num);
            for (int i=0;i<num;i++)
                HRF[i][2] = v[i] - out[i];
            for (int i=0;i<num;i++)
                u[i] = HRF[i][1];
            proj(out,v,u,num);
            for (int i=0;i<num;i++)
                HRF[i][2] -= out[i];
            
            tmp = 0;
            for (int i=0;i<num;i++)
                tmp += HRF[i][0]*HRF[i][0];
            tmp = sqrt(tmp);
            for (int i=0;i<num;i++)
                HRF[i][0] /= tmp;
            tmp = 0;
            for (int i=0;i<num;i++)
                tmp += HRF[i][1]*HRF[i][1];
            tmp = sqrt(tmp);
            for (int i=0;i<num;i++)
                HRF[i][1] /= tmp;
            tmp = 0;
            for (int i=0;i<num;i++)
                tmp += HRF[i][2]*HRF[i][2];
            tmp = sqrt(tmp);
            for (int i=0;i<num;i++)
                HRF[i][2] /= tmp;
            break;
    }
    
    
    free(out);
    free(v);
    free(u);
    
    return HRF;
}

double *convolve(double **design,double **hrf,int *dim_design,int *dim_hrf)
{
    fftw_plan  setup_fwd1;
    fftw_plan  setup_fwd2;
    fftw_plan  setup_inv;
    fftw_complex *XYZ1,*XYZ2;
    int i=0,m=0,g=0;
    double *X,*A;
    
    //	printf("%d\n",grid);
    
    g = (int)ceil(log((double)dim_design[0])/log(2.0));
    m = (int)pow(2,g);

    XYZ1 = (fftw_complex *)fftw_malloc(m*sizeof(fftw_complex));
    XYZ2 = (fftw_complex *)fftw_malloc(m*sizeof(fftw_complex));

/*    if (!fftw_init_threads()) {
        exit(0);
    }
    int nthreads = 4;
    fftw_plan_with_nthreads(nthreads);*/

    setup_fwd1 = fftw_plan_dft_1d(m,XYZ1,XYZ1,FFTW_FORWARD,FFTW_MEASURE);
    setup_fwd2 = fftw_plan_dft_1d(m,XYZ2,XYZ2,FFTW_FORWARD,FFTW_MEASURE);
    setup_inv = fftw_plan_dft_1d(m,XYZ2,XYZ2,FFTW_BACKWARD,FFTW_MEASURE);

    // setup to perform fft
    
    X = (double *)calloc(dim_design[0]*dim_design[1]*dim_hrf[1],sizeof(double));
//    for (i=0;i<dim_design[0];i++)
//        X[i] = (double *)calloc(dim_design[1]*dim_hrf[1],sizeof(double));

    for (int i=0;i<dim_design[0];i++)
        X[i*dim_design[1]*dim_hrf[1] + 0] = 1;

    for (int k=0;k<dim_hrf[1];k++) {
        for (int j=0;j<dim_design[1];j++) {
            for (i=0;i<m;i++)
                XYZ1[i][0] = XYZ1[i][1] = 0;
            for (i=0;i<dim_design[0];i++)
                XYZ1[i][0] = design[i][j];
            
            fftw_execute(setup_fwd1);
            
        
            for (i=0;i<m;i++)
                XYZ2[i][0] = XYZ2[i][1] = 0;
            
            for (i=0;i<dim_hrf[0];i++)
                XYZ2[i][0] = hrf[i][k];
        
            
            fftw_execute(setup_fwd2);
            
            // multiply
            
            double tmpR,tmpI;
            for (i=0;i<m;i++) {
                tmpR = XYZ1[i][0]*XYZ2[i][0] - XYZ1[i][1]*XYZ2[i][1];
                tmpI = XYZ1[i][0]*XYZ2[i][1] + XYZ1[i][1]*XYZ2[i][0];
                XYZ2[i][0] = tmpR;
                XYZ2[i][1] = tmpI;
            }
            
            fftw_execute(setup_inv);
            
            for (i=0;i<dim_design[0];i++)
                X[i*dim_design[1]*dim_hrf[1] + k + j*(dim_hrf[1])] = XYZ2[i][0]/(double)m;
        
        }
    }
    fftw_free(XYZ1);
    fftw_free(XYZ2);
    fftw_destroy_plan(setup_fwd1);
    fftw_destroy_plan(setup_fwd2);
    fftw_destroy_plan(setup_inv);
    fftw_cleanup();

    return X;
}


double dnorm(double *res,double *prec,int dim)
{
    void Ax(double *b,double *A,double *x,const int nrow,const int ncol); 
    double dotProd(double *x,double *y,const int dim);
    double *b,y;
      
    b = (double *)calloc(dim,sizeof(double));
       
    Ax(b,prec,res,(const int)dim,(const int)dim);
    y = -0.5*dotProd(res,b,(const int)dim);
    
    free(b);
    return y;
}

void draw_hrf_mean(REP *rep,unsigned long *seed) {
    int *dim_HRF;
    double **HRF;
    double beta = 1,old_prior,new_prior;
  
    void calculate_residuals(REP *rep,int P);
    void Ax(double *b,double *A,double *x,const int nrow,const int ncol); 
    void calculate_residuals(REP *);
    void calAx(double *Ax,double *A,double *x,const int nrow,const int ncol);
    double *subsample_design(double *mat,int len,int ncol,int *sublen,double true_freq,double samp_freq);
    void copydouble(double *B,double *A,const int dim);

   // draw proposal

    double new_m = rnorm(rep->mhrf,rep->prop_sd[2],seed);
    double v = new_m;

    rep->attempt[2] += 1;
    if ((new_m > 0) && (new_m < 12)) {
 
    old_prior = (6-1)*log(rep->mhrf) - 1*rep->mhrf;
    new_prior = (6-1)*log(new_m) - 1*new_m;

    double old_loglik = 0;
    // calculate old log likelihood
    dim_HRF = (int *)calloc(2,sizeof(double));
    for (int i=0;i<rep->dim_X[0];i++) {
        double res = rep->residuals[i];
        old_loglik += dnorm(&res,&(rep->d_Y[i]),1);
    }
    // compute new HRF
//    printf("new_m = %lf %lf\n",new_m,v);
    HRF = canonical_HRF(35,50,dim_HRF,new_m,v,0);
//    printf("dim_HRF = %d %d\n",dim_HRF[0],dim_HRF[1]);
    // convolve
    double *X;
        
    X = convolve(rep->design,HRF,rep->dim_design,dim_HRF);

/*FILE *ftmp;
ftmp = fopen("./log/tmp.log","w");
for (int i=0;i<dim_HRF[0];i++) {
    for (int j=0;j<dim_HRF[1];j++)
        fprintf(ftmp,"%lf ",HRF[i][j]);
    fprintf(ftmp,"\n");
}
fclose(ftmp);exit(0);*/
    // subsample design
            
    int subN;  
//    if (50 > 5) {
        double *XXX;     
        XXX = subsample_design(X,rep->dim_design[0],rep->dim_design[1]*dim_HRF[1],&subN,50,5); 
        free(X);   
//    }
//    else
//        XXX = X;
    // rescale design matrix
            
    double max = 0;
    double tmp = 0;

    for (int j=0;j<rep->dim_X[1];j++) { // normalize
        tmp = 0;
        for (int i=0;i<rep->dim_X[0];i++)
            tmp = (tmp > XXX[i*rep->dim_X[1]+j]) ? tmp:XXX[i*rep->dim_X[1]+j];
        max = (tmp > max) ? tmp:max;
    }
    for (int j=0;j<rep->dim_X[1];j++) { // normalize
        for (int i=0;i<rep->dim_X[0];i++)
           XXX[i*rep->dim_X[1]+j] /= max;
    }
    
    calAx(rep->Xbeta,XXX,rep->beta,(const int)rep->dim_X[0],(const int)rep->dim_X[1]);

    // calculate_residuals;
    calculate_residuals(rep,rep->P);
 
    // calculate new log likelihood
    double new_loglik = 0;
    for (int i=0;i<rep->dim_X[0];i++) {
        double res = rep->residuals[i];//printf("res = %lf\n",res);fflush(NULL);
        new_loglik += dnorm(&res,&(rep->d_Y[i]),1);
    }
    // accept reject
    double mh_ratio = (new_loglik - old_loglik) + (new_prior - old_prior);
     
//    np = old_prior = 0;
   
//    printf("%lf %lf %lf\n",new_loglik,old_loglik, mh_ratio);fflush(stdout);
    if (log(kiss(seed)) < mh_ratio) {
        rep->accept[2] += 1;
 //       printf("here\n");
        rep->mhrf = new_m;
        copydouble(rep->X,XXX,(const int)rep->dim_X[0]*rep->dim_X[1]);
    }
    else {
        calAx(rep->Xbeta,rep->X,rep->beta,(const int)rep->dim_X[0],(const int)rep->dim_X[1]);

        // calculate_residuals;
        calculate_residuals(rep,rep->P);
    }
    free(XXX);
    for (int i=0;i<dim_HRF[0];i++)
        free(HRF[i]);
    free(HRF);
    free(dim_HRF);
    } 
}

