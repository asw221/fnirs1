/*
 *  hrf.cpp
 *  fNIRS
 *
 *  Created by Timothy Johnson on 11/08/15.
 *  Copyright 2015 University of Michigan. All rights reserved.
 *
 */

#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// #include "/usr/local/include/fftw3.h"


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

double **canonical_HRF(int T,double freq,int *dim_HRF,int type)
{
    int num;
    double a1 = 6;
    double b1 = 1;
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
                HRF[i] =(double *)calloc(dim_HRF[1],sizeof(double));
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
    HRF[0][1] = 0;
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
    
    return X;
}
