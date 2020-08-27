#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "fNIRS.h"
#include "cholesky.h"
#include "randgen.h"

extern int maxP;
extern sDLM *dlmStruc;

double dotProd(double *x,double *y,const int dim)
{
    double dp = 0;
    
    for (int i=0;i<dim;i++)
        dp += x[i]*y[i];
    
    return dp;
}

void Ax(double *B,double *A,double *x,const int nrow,const int ncol)
{
    for (int i=0;i<nrow;i++)
        B[i] = dotProd(&A[i*ncol],x,ncol);
}

void cA(double *C,const double c,double *A,const int dim)
{
    for (int i=0;i<dim;i++)
        C[i] = c*A[i];
}

void copydouble(double *B,double *A,const int dim)
{
    for (int i=0;i<dim;i++)
        B[i] = A[i];
}

void transpose(double *At,double *A,const int nrow,const int ncol)
{
    for (int i=0;i<nrow;i++)
        for (int j=0;j<ncol;j++)
            At[j*nrow+i] = A[i*ncol+j];
}

void AxB(double *C,double *A,double *B,const int nrowA,const int mid,const int ncolB)
{
    for (int k=0,rB=0;k<mid;k++,rB+=ncolB) {
        //int rB = k*ncolB;
        for (int i=0,nrow=0,rA=0;i<nrowA;i++,nrow+=ncolB,rA+=mid) {
            //int rA += mid;
            for (int j=0;j<ncolB;j++)
                C[nrow+j] += A[rA+k]*B[rB+j];
        }
    }
}
 
void addtodouble(double *A,double *B,const int dim,const int sign)
{
    for (int i=0;i<dim;i++)
        A[i] += sign*B[i];
}

void adddoubles(double *C,double *A,double *B,const int dim,const int sign)
{
    for (int i=0;i<dim;i++)
        C[i] = A[i] + sign*B[i];
}

double dlm_forward_filter_draw(double *res,int N,sDLM *dlmStruc,const int Pmax, const int P, const double beta,const double delta,double *W)
{
    double n0,S0;
    double Q,*A,*CC,*R;
    double f,e;
    double tden(double x,double mean,double var,double df);
    
    A = (double *)calloc(P,sizeof(double));
    CC = (double *)calloc(P*P,sizeof(double));
    R = (double *)calloc(P*P,sizeof(double));
      
    for (int i=0;i<P*P;i++)
            dlmStruc[P-1].C[i] = 0;
    for (int i=0;i<P;i++)
            dlmStruc[P-1].C[i*P+i] = 1;
    
    
    S0 = 1;
    n0 = 1;            
    S0 = S0/n0;
    double n,d;
    
    n = n0;
    d = n0*S0;
    dlmStruc[P-1].S = d/n;

    for (int i=0;i<P;i++)
        dlmStruc[P-1].m[i] = 0;
    
    double ll = 0;
    for (int t=P;t<N;t++) {
//        for (int i=0;i<P;i++)
//            dlmStruc[t].a[i] = dlmStruc[t-1].m[i];
             
        f = 0;
        for (int i=0;i<P;i++)
            f += W[t*P+i]*dlmStruc[t-1].m[i];
//            f += W[t*P+i]*dlmStruc[t].a[i];
    
        
        for (int i=0;i<P;i++)
            for (int j=0;j<P;j++) 
                R[i*P+j] = dlmStruc[t-1].C[i*P+j]/delta;
 
        Q = dlmStruc[t-1].S;
        for (int i=0;i<P;i++)
            for (int j=0;j<P;j++)
                Q += W[t*P+i]*R[i*P+j]*W[t*P+j];
        
        for (int i=0;i<P;i++) {
            A[i] = 0;
            for (int j=0;j<P;j++)
                A[i] += R[i*P+j]*W[t*P+j];
            A[i] /= Q; 
        }                
                           
        e = res[t] - f;        
        
        if (t >= Pmax)
            ll += tden(res[t],f,Q,beta*n);

        n = n*beta + 1;
        d = beta*d + dlmStruc[t-1].S*e*e/Q;
        dlmStruc[t].S = d/n;
                    
        double tmp = dlmStruc[t].S/dlmStruc[t-1].S;
        for (int i=0;i<P;i++) {
            for (int j=0;j<P;j++) {
                dlmStruc[t].C[i*P+j] = tmp*(R[i*P+j] - A[i]*A[j]*Q);
            }
        }
                
        for (int i=0;i<P;i++)
            dlmStruc[t].m[i] = dlmStruc[t-1].m[i] + A[i]*e;
 //           dlmStruc[t].m[i] = dlmStruc[t].a[i] + A[i]*e;
 
        transpose(CC,dlmStruc[t].C,(const int)P,(const int)P);
        addtodouble(CC,dlmStruc[t].C,(const int)P*P,1);
        cA(dlmStruc[t].C,(const double)(1./2.),CC,(const int)P*P);
 
    }
    free(CC);
    free(A);
    free(R);
    return ll;
}

void dlm_forward_filter(REP *rep,sDLM *dlmStruc)
{
    int P;
    double n0,S0,*CC,*R;
    double Q,*A;
    double f,e;
 
    P = rep->P;  
    A = (double *)calloc(P,sizeof(double));
    CC = (double *)calloc(P*P,sizeof(double));
    R = (double *)calloc(P*P,sizeof(double));
        
    for (int i=0;i<P*P;i++)
            dlmStruc[P-1].C[i] = 0;
    for (int i=0;i<P;i++)
            dlmStruc[P-1].C[i*P+i] = 1;
 
    S0 = 1;
    n0 = 1;            
    S0 = S0/n0;
    
    dlmStruc[P-1].n = n0;
    dlmStruc[P-1].d = n0*S0;
    dlmStruc[P-1].S = S0;

    for (int i=0;i<P;i++)
        dlmStruc[P-1].m[i] = 0;
    
    
    for (int t=P;t<rep->dim_X[0];t++) {
//        for (int i=0;i<P;i++)
//            dlmStruc[t].a[i] = dlmStruc[t-1].m[i];
             
        f = 0;
        for (int i=0;i<P;i++)
            f += rep->W[t*P+i]*dlmStruc[t-1].m[i];
//            f += rep->W[t*P+i]*dlmStruc[t].a[i];
    
        
        for (int i=0;i<P;i++)
            for (int j=0;j<P;j++) 
                R[i*P+j] = dlmStruc[t-1].C[i*P+j]/rep->df_delta1;
 
        Q = dlmStruc[t-1].S;
        for (int i=0;i<P;i++)
            for (int j=0;j<P;j++)
                Q += rep->W[t*P+i]*R[i*P+j]*rep->W[t*P+j];
                
        for (int i=0;i<P;i++) {
            A[i] = 0;
            for (int j=0;j<P;j++)
                A[i] += R[i*P+j]*rep->W[t*P+j];
            A[i] /= Q; 
        }                
                           
        e = rep->residuals3[t] - f;        
        
        dlmStruc[t].n = dlmStruc[t-1].n*rep->df_delta2 + 1;
        dlmStruc[t].d = rep->df_delta2*dlmStruc[t-1].d + dlmStruc[t-1].S*e*e/Q;
        if (isnan(dlmStruc[t].d)) {
            printf("\n forward filter error\n");
            printf("delta2 = %lf, d[t-1] = %lf, S[t-1] = %lf, e = %lf, Q = %lf\n",
                     rep->df_delta2,dlmStruc[t-1].d,dlmStruc[t-1].S,e,Q);
            printf("t = %d, res[t] = %lf f = %lf\n",t,rep->residuals3[t],f);
            printf("W = \n");
            for (int i=0;i<P;i++)
                printf("%lf ",rep->W[t*P+i]);
            printf("m = \n");
            for (int i=0;i<P;i++)
                printf("%lf ",dlmStruc[t-1].m[i]);
            fflush(NULL);
             exit(1);
        }
        dlmStruc[t].S = dlmStruc[t].d/dlmStruc[t].n;
                    
        double tmp = dlmStruc[t].S/dlmStruc[t-1].S;
        for (int i=0;i<P;i++) {
            for (int j=0;j<P;j++) {
                dlmStruc[t].C[i*P+j] = tmp*(R[i*P+j] - A[i]*A[j]*Q);
            }
        }
                
        for (int i=0;i<P;i++)
            dlmStruc[t].m[i] = dlmStruc[t-1].m[i] + A[i]*e;
//           dlmStruc[t].m[i] = dlmStruc[t].a[i] + A[i]*e;

        transpose(CC,dlmStruc[t].C,(const int)P,(const int)P);
        addtodouble(CC,dlmStruc[t].C,(const int)P*P,1);
        cA(dlmStruc[t].C,(const double)(1./2.),CC,(const int)P*P);
 
    }
    free(CC);
    free(A);
    free(R);
}

void dlm_backward_sampling(REP *rep,sDLM *dlmStruc,const int P,unsigned long *seed)
{    
    double *m,*Var,*V2,rg;
       
    m = (double *)calloc(P,sizeof(double));
    Var = (double *)calloc(P*P,sizeof(double));   
    V2 = (double *)calloc(P*P,sizeof(double));   

    for (int t=rep->dim_X[0]-2;t>=P;t--) {
        
        for (int i=0;i<P;i++)
            m[i] = dlmStruc[t].m[i] + rep->df_delta1*(rep->delta[(t+1)*P + i] - dlmStruc[t].m[i]);
//            m[i] = dlmStruc[t].m[i] + rep->df_delta1*(rep->delta[(t+1)*P + i] - dlmStruc[t+1].a[i]);


        // compute mean
                
        double tmp = 1.-rep->df_delta1;
        for (int i=0;i<P*P;i++)
            Var[i] = dlmStruc[t].C[i]*tmp;
                
        for (int i=0;i<P*P;i++)
            V2[i] = Var[i];   
        // draw new delta[t]
        int err = cholesky_decomp2vec(Var,P);
        if (err) {  // err = 1 means C is SPD
            err = rmvtvec(&(rep->delta[t*P]),Var,P,dlmStruc[t].n,m,seed);
        int flag = 0;
        for (int i=0;i<P;i++) {
            if (isnan(rep->delta[t*P+i]))
                flag = 1;
        }
        if (flag) {
            printf("P = %d, df = %lf\n",P,dlmStruc[t].n);
            printf("delta  = \n");
            for (int i=0;i<P;i++)
                printf("%lf ",rep->delta[t*P+i]);
            printf("\n m = \n");
            for (int i=0;i<P;i++)
                printf("%lf ",dlmStruc[t].m[i]);
            printf("\n\n");
            for (int i=0;i<P*P;i++)
                printf("%lf \n",Var[i]);
           fflush(NULL);
           exit(0);
        }
        }
        else {
            printf("error in dlm_backward_sampling, var is not SPD\n");
            printf("P = %d %lf %lf %s\n",P,rep->df_delta1,rep->df_delta2,rep->dataname);
            printf("Var = \n");
            for (int i=0;i<P*P;i++)
                printf("%lf ",V2[i]);
            fflush(NULL);
            printf("C = \n");
            for (int i=0;i<P*P;i++)
               printf("%lf ",dlmStruc[t].C[i]);
            printf("1 - dfdelta1 = %lf\n",1-rep->df_delta1);
            fflush(NULL);
            exit(0);
        }
        
        // draw new prec[t]
        rep->d_Y[t] = rep->df_delta2*rep->d_Y[t+1];
        rg = rgamma(0.5*(1-rep->df_delta2)*dlmStruc[t].n,0.5*dlmStruc[t].d,seed);
        rep->d_Y[t] += rg;
        
//     if (isnan(rep->d_Y[t])) {
/*     if (rep->d_Y[t] > 1e4) {
            printf("what2 %lf %lf\n",0.5*(1-rep->df_delta2)*dlmStruc[t].n,0.5*dlmStruc[t].d);fflush(NULL);
            printf("rgamma = %lf\n",rg);}
     
     int flag = 0;
    
     if (isnan(rep->d_Y[t])) 
            flag = 1;
     if (flag) {
             printf("\n%d %lf\n",t,rep->d_Y[t]);
             printf("%lf, %lf, %lf %lf\n",rep->df_delta2,dlmStruc[t].n,dlmStruc[t].d,rep->d_Y[t+1]);
             printf("alpha = %lf beta = %lf\n",0.5*(1-rep->df_delta2)*dlmStruc[t].n,0.5*dlmStruc[t].d);
          printf("%s DLM %lf %lf %d\n",rep->dataname,rep->df_delta1,rep->df_delta2,rep->P);fflush(NULL);
        exit(0);
     }*/
     
    }
    free(m);
    free(Var);
    free(V2);
}

double DLMloglik(REP *rep,double sfreq,int NITER,unsigned long *seed) {
    int P,Pmax,prop_P;
    double new_loglik,max_loglik;
    double beta,delta;
    double prop_delta;
    double prop_beta,factor = 5.;
    double minb,mind;
    double maxb,maxd;
    sDLM *dlmStruc;

    double *W,*xbeta,*veta;
    void standardize_data(double *x,const int len);       
    void kernel_reg(double *f,double *Y,double *X,double lambda,int N);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    double runif_atob(unsigned long *,double,double);
      void calculate_res3(REP *rep);

    int N = rep->N;
    double *f = (double *)calloc(N,sizeof(double));
    double *data = (double *)calloc(N,sizeof(double));
    double *Y = (double *)calloc(N,sizeof(double));
    
     for (int i=0;i<N;i++) {
        Y[i] = rep->Y[i];
        data[i] = i; 
    }          
    standardize_data(Y,(const int)N);
    kernel_reg(f,Y,data,factor*sfreq,N);
    if (NITER < 100000) {
    for (int i=0;i<N;i++) { 
        f[i] = Y[i] - f[i];
    }
    }
    else{   
        calculate_res3(rep);
        for (int i=0;i<N;i++) { 
            f[i] = rep->residuals3[i];
        }
    }
/*    FILE *fkernel;
    char *S = (char *)calloc(300,sizeof(char));
    S = strcpy(S,"kern_smooth_");
    S = strcat(S,rep->dataname);
    S = strcat(S,".dat");
    fkernel = fopen(S,"w");
    for (int i=0;i<N;i++)
        fprintf(fkernel,"%lf ",f[i]);
    fclose(fkernel);
    free(S); */
    max_loglik = -1e100;
    Pmax = 20;
    W = (double *)calloc(N*Pmax,sizeof(double));
    xbeta = (double *)calloc(N*Pmax,sizeof(double));
    veta = (double *)calloc(N*Pmax,sizeof(double)); 
    
    int dimsize = Pmax;
    dlmStruc = (sDLM *)calloc(N,sizeof(sDLM));
    for (int k=0;k<N;k++)
        dlmStruc[k].m = (double *)calloc(dimsize,sizeof(double));
    for (int k=0;k<N;k++)
        dlmStruc[k].C = (double *)calloc(dimsize*dimsize,sizeof(double));

    mind = 0.90;maxd=1;
    minb = 0.80;maxb=1;
    for (int prop_P=1;prop_P<=Pmax;prop_P++) {
        calW(W,f,xbeta,veta,(const int)N,(const int)prop_P,prop_P);
        for (prop_delta=mind;prop_delta<maxd;prop_delta+=0.005) {
            for (prop_beta=minb;prop_beta<maxb;prop_beta+=0.005) {
                new_loglik = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,
                  (const int)prop_P,(const double)prop_beta,(const double)prop_delta,W);
                if (new_loglik > max_loglik) {
                    P = prop_P;
                    beta = prop_beta;
                    delta = prop_delta;
                    max_loglik = new_loglik;
                }
            } 
        } 
     }
/*    mind = 0.90;maxd=0.99999;
    minb = 0.70;maxb=0.99999;
    for (int prop_P=1;prop_P<=Pmax;prop_P++) {
        calW(W,f,xbeta,veta,(const int)N,(const int)prop_P,prop_P);
        for (int i=0;i<NITER;i++) {
  //      prop_P = runiform_n(Pmax,seed)+1;
  //      calW(W,f,xbeta,veta,(const int)N,(const int)prop_P,prop_P);
            prop_delta = runif_atob(seed,mind,maxd);
            prop_beta = runif_atob(seed,minb,maxb);
            new_loglik = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,
                  (const int)prop_P,(const double)prop_beta,(const double)prop_delta,W);
            if (new_loglik > max_loglik) {
                P = prop_P;
                beta = prop_beta;
                delta = prop_delta;
                max_loglik = new_loglik;
            }
        }   
     }*/
/*     for (int prop_P=1;prop_P<=Pmax;prop_P++) {
        calW(W,f,xbeta,veta,(const int)N,(const int)prop_P,prop_P);
        minb = (0.70 > (beta - 0.075)) ? 0.70:(beta - 0.075);
        maxb = (0.99999 < (beta + 0.075)) ? 0.99999:(beta + 0.075);
        mind = (0.90 > (delta - 0.025)) ? 0.90:(delta - 0.025);
        maxd = (0.99999 < (delta + 0.025)) ? 0.99999:(delta + 0.025);
        for (int i=0;i<200;i++) {
            prop_delta = runif_atob(seed,mind,maxd);
            prop_beta = runif_atob(seed,minb,maxb);
            new_loglik = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)prop_beta,(const double)prop_delta,W);
            if (new_loglik > max_loglik) {
                rep->P = P = prop_P;
                beta = prop_beta;
                delta = prop_delta;
                max_loglik = new_loglik;
            }
        }
        }
    for (int prop_P=1;prop_P<=Pmax;prop_P++) {
        calW(W,f,xbeta,veta,(const int)N,(const int)prop_P,prop_P);
        minb = (0.70 > (beta - 0.0375)) ? 0.70:(beta - 0.0375);
        maxb = (0.99999 < (beta + 0.0375)) ? 0.99999:(beta + 0.0375);
        mind = (0.90 > (delta - 0.0125)) ? 0.90:(delta - 0.0125);
        maxd = (0.99999 < (delta + 0.0125)) ? 0.99999:(delta + 0.0125);
        for (int i=0;i<200;i++) {
            prop_delta = runif_atob(seed,mind,maxd);
            prop_beta = runif_atob(seed,minb,maxb);
            new_loglik = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)prop_beta,(const double)prop_delta,W);
            if (new_loglik > max_loglik) {
                rep->P = P = prop_P;
                beta = prop_beta;
                delta = prop_delta;
                max_loglik = new_loglik;
            }
        }
    }*/
  
    rep->P = P;
    rep->df_delta1 = delta;
    rep->df_delta2 = beta;
    free(f);
    free(data);
    free(Y);
    free(W);
    free(xbeta);
    free(veta);
    for (int k=0;k<N;k++)
        free(dlmStruc[k].m);
    for (int k=0;k<N;k++)
        free(dlmStruc[k].C);
    free(dlmStruc);
    return (max_loglik);
}

double DLMloglik2(REP *rep,double sfreq,unsigned long *seed) {
    int P,Pmax;
    double new_loglik,max_loglik;
    double beta,delta,u,l;
    double prop_delta,fl,fu;
    double prop_beta,factor = 5.;
    double minb,mind;
    double maxb,maxd;
    double maxbeta,maxdelta;
    sDLM *dlmStruc;

    double *W,*xbeta,*veta;
    void standardize_data(double *x,const int len);       
    void kernel_reg(double *f,double *Y,double *X,double lambda,int N);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    double runif_atob(unsigned long *,double,double);
  
    int N = rep->N;
    double *f = (double *)calloc(N,sizeof(double));
    double *data = (double *)calloc(N,sizeof(double));
    double *Y = (double *)calloc(N,sizeof(double));
    
     for (int i=0;i<N;i++) {
        Y[i] = rep->Y[i];
        data[i] = i; 
    }          
    standardize_data(Y,(const int)N);
    kernel_reg(f,Y,data,factor*sfreq,N);
    FILE *fkernel;
    char *S = (char *)calloc(300,sizeof(char));
//    char *sss = (char *)calloc(20,sizeof(char));
    S = strcpy(S,"kern_smooth_");
    S = strcat(S,rep->dataname);
    S = strcat(S,".dat");
    fkernel = fopen(S,"w");
    for (int i=0;i<N;i++)
        fprintf(fkernel,"%lf ",Y[i] - f[i]);
    fclose(fkernel);
    for (int i=0;i<N;i++) { 
        f[i] = Y[i] - f[i];
    }
    free(S); 
    max_loglik = -1e100;
    Pmax = 20;
    W = (double *)calloc(N*Pmax,sizeof(double));
    xbeta = (double *)calloc(N*Pmax,sizeof(double));
    veta = (double *)calloc(N*Pmax,sizeof(double)); 
    
    int dimsize = Pmax;
    dlmStruc = (sDLM *)calloc(N,sizeof(sDLM));
    for (int k=0;k<N;k++)
        dlmStruc[k].m = (double *)calloc(dimsize,sizeof(double));
    for (int k=0;k<N;k++)
        dlmStruc[k].C = (double *)calloc(dimsize*dimsize,sizeof(double));

 //   printf("delta beta P\n");fflush(NULL);
    for (int prop_P=1;prop_P<=Pmax;prop_P++) {
        calW(W,f,xbeta,veta,(const int)N,(const int)prop_P,prop_P);
        
        minb = 0.50;
        maxb = 0.99999;
        l = mind = 0.50;
        u = maxd = 0.99999;
        beta = minb;
        fl = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)beta,(const double)mind,W);
        fu = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)beta,(const double)maxd,W);
        while (fabs(fu-fl) > 1e-6) {
            if (fl <= fu) {
                l = (l+u)/2.;
            }
            else {
                u = (l+u)/2.;
            }
            fl = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)beta,(const double)l,W);
            fu = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)beta,(const double)u,W);
        }
        delta = (u+l)/2.;
        fl = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)minb,(const double)delta,W);
        fu = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)maxb,(const double)delta,W);
        l = minb;
        u = maxb;
        while (fabs(fu-fl) > 1e-6) {
            if (fl <= fu) {
                l = (l+u)/2.;
            }
            else {
                u = (l+u)/2.;
            }
            fl = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)l,(const double)delta,W);
            fu = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)u,(const double)delta,W);
        }
        beta = (u+l)/2.;
        new_loglik = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)beta,(const double)delta,W);
        if (new_loglik >= max_loglik) {
            max_loglik = new_loglik;
            maxbeta = beta;
            maxdelta = delta;
            P = prop_P;
        }
 //        printf("P = %d max_delta = %lf max_beta = %lf %.10g\n",prop_P,delta,beta,new_loglik);

     }

    rep->P = P;
    rep->df_delta1 = maxdelta;
    rep->df_delta2 = maxbeta;
    free(f);
    free(data);
    free(Y);
    free(W);
    free(xbeta);
    free(veta);
    for (int k=0;k<N;k++)
        free(dlmStruc[k].m);
    for (int k=0;k<N;k++)
        free(dlmStruc[k].C);
    free(dlmStruc);
    return (max_loglik);
}

double DLMloglik3(REP *rep,double sfreq,unsigned long *seed) {
    int P,Pmax;
    double new_loglik,max_loglik;
    double beta,delta,u,l;
    double prop_delta,fl,fu;
    double prop_beta,factor = 5.;
    double minb,mind;
    double maxb,maxd;
    double old_beta,old_delta;
    double maxbeta,maxdelta;
    sDLM *dlmStruc;

    double *W,*xbeta,*veta;
    void standardize_data(double *x,const int len);       
    void kernel_reg(double *f,double *Y,double *X,double lambda,int N);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    double runif_atob(unsigned long *,double,double);
  
    int N = rep->N;
    double *f = (double *)calloc(N,sizeof(double));
    double *data = (double *)calloc(N,sizeof(double));
    double *Y = (double *)calloc(N,sizeof(double));
    
     for (int i=0;i<N;i++) {
        Y[i] = rep->Y[i];
        data[i] = i; 
    }          
    standardize_data(Y,(const int)N);
    kernel_reg(f,Y,data,factor*sfreq,N);
/*    FILE *fkernel;
    char *S = (char *)calloc(300,sizeof(char));
    S = strcpy(S,"kern_smooth_");
    S = strcat(S,rep->dataname);
    S = strcat(S,".dat");
    fkernel = fopen(S,"w");
    for (int i=0;i<N;i++)
        fprintf(fkernel,"%lf ",Y[i] - f[i]);
    fclose(fkernel);
    free(S); */
    for (int i=0;i<N;i++) { 
        f[i] = Y[i] - f[i];
    }
    max_loglik = -1e100;
    Pmax = 20;
    W = (double *)calloc(N*Pmax,sizeof(double));
    xbeta = (double *)calloc(N*Pmax,sizeof(double));
    veta = (double *)calloc(N*Pmax,sizeof(double)); 
    
    int dimsize = Pmax;
    dlmStruc = (sDLM *)calloc(N,sizeof(sDLM));
    for (int k=0;k<N;k++)
        dlmStruc[k].m = (double *)calloc(dimsize,sizeof(double));
    for (int k=0;k<N;k++)
        dlmStruc[k].C = (double *)calloc(dimsize*dimsize,sizeof(double));

 //   printf("delta beta P\n");fflush(NULL);
    for (int prop_P=1;prop_P<=Pmax;prop_P++) {
        calW(W,f,xbeta,veta,(const int)N,(const int)prop_P,prop_P);
  //      printf("prop_P = %d\n",prop_P);fflush(NULL);
        mind = 0.90;
        maxd = 0.99999;
        minb = 0.70;
        maxb = 0.99999;
        beta = (maxb-minb)/2.;
        delta = (maxd-mind)/2.;
        old_delta = old_beta = 0;
        while ((fabs(beta - old_beta) > 1e-4) && (fabs(delta - old_delta) > 1e-4)) {
 //       printf("%lf %lf\n",fabs(beta - old_beta),fabs(delta - old_delta));fflush(NULL);
            l = mind;
            u = maxd;
            fl = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)beta,(const double)mind,W);
            fu = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)beta,(const double)maxd,W);
            int counter = 0;
            while (fabs(fu-fl) > 1e-4) {
//                printf("beta %lf\n",fabs(fu-fl));fflush(NULL);
                if (fl < fu) {
                    l = (l+u)/2.;
                }
                else {
                    u = (l+u)/2.;
                }
                fl = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)beta,(const double)l,W);
                fu = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)beta,(const double)u,W);
                counter++;
                if (counter > 20)
                    break;
            }
            old_delta = delta;
            delta = (u+l)/2.;
            fl = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)minb,(const double)delta,W);
            fu = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)maxb,(const double)delta,W);
            l = minb;
            u = maxb;
            counter = 0;
            while (fabs(fu-fl) > 1e-4) {
//               printf("delta %lf\n",fabs(fu-fl));fflush(NULL);
              if (fl < fu) {
                    l = (l+u)/2.;
                }
                else {
                    u = (l+u)/2.;
                }
                fl = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)l,(const double)delta,W);
                fu = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)u,(const double)delta,W);
                counter++;
                if (counter > 20)
                    break;
            }
            old_beta = beta;
            beta = (u+l)/2.;
        }
        new_loglik = dlm_forward_filter_draw(f,N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)beta,(const double)delta,W);
        if (new_loglik > max_loglik) {
            max_loglik = new_loglik;
            maxbeta = beta;
            maxdelta = delta;
            P = prop_P;
        }
 //        printf("P = %d max_delta = %lf max_beta = %lf %.10g\n",prop_P,delta,beta,new_loglik);
    }

    rep->P = P;
    rep->df_delta1 = maxdelta;
    rep->df_delta2 = maxbeta;
    free(f);
    free(data);
    free(Y);
    free(W);
    free(xbeta);
    free(veta);
    for (int k=0;k<N;k++)
        free(dlmStruc[k].m);
    for (int k=0;k<N;k++)
        free(dlmStruc[k].C);
    free(dlmStruc);
    return (max_loglik);
}

void DLMtst(REP *rep,sDLM *dlmStruc,dlm_tst_flag flag,unsigned long *seed) {
    double old_loglik,new_loglik;
    double *W;
    void calculate_res3(REP *rep);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    double tden(double x,double mean,double var,double df);
    
    int P = rep->P;
    W = rep->W;

    calculate_res3(rep);

    int prop_P;
    double new_logprop,old_logprop;

    if (flag == fP) {
        double *CDF1,*PDF1;
        CDF1 = (double *)calloc((maxP-1+1),sizeof(double));
        PDF1 = (double *)calloc((maxP-1+1),sizeof(double));
        double lfact = 0;
        for (int i=1;i<=maxP;i++) {
            lfact += log((double)i);
            PDF1[i-1] =  exp(-P + log(pow((double)P,(double)i)) - lfact);
            CDF1[i-1] = PDF1[i-1];
        }
        for (int i=2;i<=maxP;i++)
            CDF1[i-1] += CDF1[i-2];
        
        for (int i=1;i<=maxP;i++)
            PDF1[i-1] /= CDF1[maxP-1];
        for (int i=1;i<=maxP;i++)
            CDF1[i-1] /= CDF1[maxP-1];
        
        prop_P = trunc_rpois(CDF1,1,maxP,seed);
        new_logprop = log(PDF1[prop_P-1]);
 
        lfact = 0;
        for (int i=1;i<=maxP;i++) {
            lfact += log((double)i);
            PDF1[i-1] =  exp(-prop_P + log(pow((double)prop_P,(double)i)) - lfact);
            CDF1[i-1] = PDF1[i-1];
        }
        for (int i=2;i<=maxP;i++)
            CDF1[i-1] += CDF1[i-2];
        
        for (int i=1;i<=maxP;i++)
            PDF1[i-1] /= CDF1[maxP-1];
        for (int i=1;i<=maxP;i++)
            CDF1[i-1] /= CDF1[maxP-1];
 
        old_logprop = log(PDF1[P-1]);
  
        free(PDF1);
        free(CDF1);
    }
    else {
        prop_P = P;
    }


    int Pmax = (P > prop_P) ? P:prop_P;
    
    // Calculate old log likelihood
    
    old_loglik = dlm_forward_filter_draw(rep->residuals3,rep->dim_X[0],dlmStruc,(const int)Pmax,(const int)P,(const double)rep->df_delta2,(const double)rep->df_delta1,W);
//    old_loglik = dlm_forward_filter_draw(rep,dlmStruc,(const int)Pmax,(const int)P,(const double)rep->df_delta2,(const double)rep->df_delta1,W);
    
    double prop_delta,prop_beta;
    double low,high;
    
    double range = rep->prop_sd[0];
    low = rep->df_delta1 - range;
    low = (low < 0.75) ? 0.75:low;
    high = rep->df_delta1 + range;
    high = (high > 0.999999) ? 0.999999:high;
   
    prop_delta = runif_atob(seed,low,high);
    if (flag == fdelta1)
        new_logprop = -log(high-low);
 
    low = prop_delta - range;
    low = (low < 0.75) ? 0.75:low;
    high = prop_delta + range;
    high = (high > 0.999999) ? 0.999999:high;
    if (flag == fdelta1)
        old_logprop = -log(high-low);
    
    range = rep->prop_sd[1];
    low = rep->df_delta2 - range;
    low = (low < 0.50) ? 0.50:low;
    high = rep->df_delta2 + range;
    high = (high > 0.999999) ? 0.999999:high;
    
    prop_beta = runif_atob(seed,low,high);
    if (flag == fdelta2)
        new_logprop = -log(high-low);
 
    low = prop_beta - range;
    low = (low < 0.50) ? 0.50:low;
    high = prop_beta + range;
    high = (high > 0.999999) ? 0.999999:high;
    if (flag == fdelta2)
        old_logprop = -log(high-low);
     
    if (flag == fdelta1) {
        prop_beta = rep->df_delta2;
    }   
    else if (flag == fdelta2)
        prop_delta = rep->df_delta1;
    else if (flag == fP) {
        prop_delta = rep->df_delta1;
        prop_beta = rep->df_delta2;
    }
    
    if (prop_P != P) {
        calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)prop_P,P);
        W = rep->W;
    }
    
    // Calculate new log likelihood
     
    new_loglik = dlm_forward_filter_draw(rep->residuals3,rep->dim_X[0],dlmStruc,(const int)Pmax,(const int)prop_P,(const double)prop_beta,(const double)prop_delta,W);
//    new_loglik = dlm_forward_filter_draw(rep,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)prop_beta,(const double)prop_delta,W);
  
    double new_log_prior,old_log_prior;
    switch (flag) {
        case fdelta1:
//           new_log_prior = (0.999*(0.2*rep->dim_X[0]) - 1)*log(prop_delta) + (0.001*(0.2*rep->dim_X[0]) - 1)*log(1.-prop_delta);
//            old_log_prior = (0.999*(0.2*rep->dim_X[0]) - 1)*log(rep->df_delta1) + (0.001*(0.2*rep->dim_X[0]) - 1)*log(1.-rep->df_delta1);
            new_log_prior = old_log_prior = 0;
            break;        
        case fdelta2: 
//            new_log_prior = (0.8*(0.2*rep->dim_X[0]) - 1)*log(prop_beta) + (0.2*(0.2*rep->dim_X[0]) - 1)*log(1.-prop_beta);
//            old_log_prior = (0.8*(0.2*rep->dim_X[0]) - 1)*log(rep->df_delta2) + (0.2*(0.2*rep->dim_X[0]) - 1)*log(1.-rep->df_delta2);
              new_log_prior = old_log_prior = 0;           
            break;
        case fP: default:
            new_log_prior = old_log_prior = 0;
            break;
    }
    if (log(kiss(seed)) < ((new_loglik - old_loglik) + (new_log_prior - old_log_prior) - (new_logprop - old_logprop))) {
        rep->df_delta1 = prop_delta;
        rep->df_delta2 = prop_beta;
        rep->P = prop_P;
        rep->dim_W[1] = rep->P;
        if (flag != fP)
            (rep->accept[flag])++;
    }
    if (flag != fP)
        (rep->attempt[flag])++;

    calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],rep->P);
}

void DLMtst2(REP *rep,sDLM *dlmStruc,dlm_tst_flag flag,unsigned long *seed) {
    double old_loglik,new_loglik;
    double *W;
    sDLM *tmp_dlmStruc;    
    void calculate_res3(REP *rep);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    double tden(double x,double mean,double var,double df);
    
    int P = rep->P;
    W = rep->W;

    calculate_res3(rep);

    int prop_P;
    double new_logprop,old_logprop;

    if (flag == fP) {
        if (P == 1) {
            if (kiss(seed) < 0.5) {
                prop_P = P+1;
                old_logprop = log(1./3.);
            }
            else {
                prop_P = P+2;
                old_logprop = log(.25);
            }
            new_logprop = log(0.5);
        }
        else if (P==2) {
            if (kiss(seed) < 1./3.) {
                prop_P = P-1;
                old_logprop = log(0.5);
                new_logprop = log(1./3.);
            }
            else {
                if (kiss(seed) < 0.5) {
                    prop_P = P+1;
                    old_logprop = log(0.25);
                    new_logprop = log(1./3.);
                }
                else {
                    prop_P = P+2;
                    old_logprop = log(0.25);
                    new_logprop = log(1./3.);
                }
            }
        }
        else if (P == maxP-1) {
            if (kiss(seed) < 0.5) {
                prop_P = P-1;
                old_logprop = log(1./3.);
            }
            else {
                prop_P = P-2;
                old_logprop = log(.25);
            }
            new_logprop = log(0.5);
        }
        else if (P == maxP-2) {
            if (kiss(seed) < 1./3.) {
                prop_P = P+1;
                old_logprop = log(0.5);
                new_logprop = log(1./3.);
            }
            else {
                if (kiss(seed) < 0.5) {
                    prop_P = P-1;
                    old_logprop = log(0.25);
                    new_logprop = log(1./3.);
                }
                else {
                    prop_P = P-2;
                    old_logprop = log(0.25);
                    new_logprop = log(1./3.);
                }
            }
        }
        else {
            if (kiss(seed) < 0.5) { //drop down
                if (kiss(seed) < 0.5) 
                    prop_P = P-2;
                else 
                    prop_P = P-1;
                
                new_logprop = log(0.25);
                if (prop_P > 2) 
                    old_logprop = log(0.25);
                else if (prop_P == 1)
                    old_logprop = log(0.5);
                else //(prop_P == 2)
                    old_logprop = log(1./3.);
            }
            else { // pop up
                if (kiss(seed) < 0.5) 
                    prop_P = P+1;
                else 
                    prop_P = P+2;
               
                new_logprop = log(0.25);
                if (prop_P < maxP-2)
                    old_logprop = log(0.25);
                else if (prop_P == maxP-1)
                    old_logprop = log(0.5);
                else // (prop_P == maxP-2)
                    old_logprop = log(1./3.);
            }
        }
    }
    else {
        prop_P = P;
    }


    int Pmax = (P > prop_P) ? P:prop_P;
    
    // Calculate old log likelihood
    
    old_loglik = dlm_forward_filter_draw(rep->residuals3,rep->N,dlmStruc,(const int)Pmax,(const int)P,(const double)rep->df_delta2,(const double)rep->df_delta1,W);
    
    double prop_delta,prop_beta;
    double low,high;
    
    double range = rep->prop_sd[0];
    low = rep->df_delta1 - range;
    low = (low < 0.75) ? 0.75:low;
    high = rep->df_delta1 + range;
    high = (high > 0.999999) ? 0.999999:high;
   
    prop_delta = runif_atob(seed,low,high);
    if (flag == fdelta1)
        new_logprop = -log(high-low);
 
    low = prop_delta - range;
    low = (low < 0.75) ? 0.75:low;
    high = prop_delta + range;
    high = (high > 0.999999) ? 0.999999:high;
    if (flag == fdelta1)
        old_logprop = -log(high-low);
    
    range = rep->prop_sd[1];
    low = rep->df_delta2 - range;
    low = (low < 0.50) ? 0.50:low;
    high = rep->df_delta2 + range;
    high = (high > 0.999999) ? 0.999999:high;
    
    prop_beta = runif_atob(seed,low,high);
    if (flag == fdelta2)
        new_logprop = -log(high-low);
 
    low = prop_beta - range;
    low = (low < 0.50) ? 0.50:low;
    high = prop_beta + range;
    high = (high > 0.999999) ? 0.999999:high;
    if (flag == fdelta2)
        old_logprop = -log(high-low);
     
    if (flag == fdelta1) {
        prop_beta = rep->df_delta2;
    }   
    else if (flag == fdelta2)
        prop_delta = rep->df_delta1;
    else if (flag == fP) {
        prop_delta = rep->df_delta1;
        prop_beta = rep->df_delta2;
    }
    
 //   if (prop_P != P) {
 //       calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)prop_P,P);
 //       W = rep->W;
 //   }
    
    // Calculate new log likelihood
    if (flag == fP) {
        W = (double *)calloc(rep->dim_W[0]*prop_P,sizeof(double));
        calW(W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)prop_P,P);
 //       W = rep->W;
// create new dlmStruc and call with it
        int dimsize = prop_P;
        tmp_dlmStruc = (sDLM *)calloc(rep->dim_W[0],sizeof(sDLM));
        for (int k=0;k<rep->dim_W[0];k++)
            tmp_dlmStruc[k].m = (double *)calloc(dimsize,sizeof(double));
        for (int k=0;k<rep->dim_W[0];k++)
            tmp_dlmStruc[k].C = (double *)calloc(dimsize*dimsize,sizeof(double));

        new_loglik = dlm_forward_filter_draw(rep->residuals3,rep->N,tmp_dlmStruc,(const int)Pmax,(const int)prop_P,(const double)prop_beta,(const double)prop_delta,W);
    }
    else {
         new_loglik = dlm_forward_filter_draw(rep->residuals3,rep->N,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)prop_beta,(const double)prop_delta,rep->W);
    }
     
 //   new_loglik = dlm_forward_filter_draw(rep,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)prop_beta,(const double)prop_delta,W);
  
    double new_log_prior,old_log_prior;
    switch (flag) {
        case fdelta1:
//           new_log_prior = (0.999*(0.2*rep->dim_X[0]) - 1)*log(prop_delta) + (0.001*(0.2*rep->dim_X[0]) - 1)*log(1.-prop_delta);
//            old_log_prior = (0.999*(0.2*rep->dim_X[0]) - 1)*log(rep->df_delta1) + (0.001*(0.2*rep->dim_X[0]) - 1)*log(1.-rep->df_delta1);
            new_log_prior = old_log_prior = 0;
            break;        
        case fdelta2: 
//            new_log_prior = (0.8*(0.2*rep->dim_X[0]) - 1)*log(prop_beta) + (0.2*(0.2*rep->dim_X[0]) - 1)*log(1.-prop_beta);
//            old_log_prior = (0.8*(0.2*rep->dim_X[0]) - 1)*log(rep->df_delta2) + (0.2*(0.2*rep->dim_X[0]) - 1)*log(1.-rep->df_delta2);
              new_log_prior = old_log_prior = 0;           
            break;
        case fP: default:
            new_log_prior = old_log_prior = 0;
            break;
    }
    if (log(kiss(seed)) < ((new_loglik - old_loglik) + (new_log_prior - old_log_prior) - (new_logprop - old_logprop))) {
        rep->df_delta1 = prop_delta;
        rep->df_delta2 = prop_beta;
        rep->P = prop_P;
        rep->dim_W[1] = rep->P;
        if (flag != fP)
            (rep->accept[flag])++;
        if (flag == fP) {
            free(rep->W);
            rep->W = W;
            // free rep->dlmStuct and point it to new one
            int len = rep->dim_X[0];
            for (int k=0;k<len;k++)
                free(rep->dlmStruc[k].m);
            for (int k=0;k<len;k++)
                free(rep->dlmStruc[k].C);
            free(rep->dlmStruc);
            rep->dlmStruc = tmp_dlmStruc;
        }
    }
    else {
        if (flag == fP) {
            free(W);
            // free up tmp dlmStruc
            int len = rep->dim_X[0];
            for (int k=0;k<len;k++)
                free(tmp_dlmStruc[k].m);
            for (int k=0;k<len;k++)
                free(tmp_dlmStruc[k].C);
            free(tmp_dlmStruc);
        }
    }
    if (flag != fP)
        (rep->attempt[flag])++;

//    calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],rep->P);
}

void sampleDLM(REP *rep,sDLM *dlmStruc,unsigned long *seed) {
    double *Var;
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calculate_residuals(REP *rep,int P);
    void calculate_res3(REP *);
  
    int P = rep->P;
    calculate_res3(rep);

    // Forward_filter
  
    dlm_forward_filter(rep,dlmStruc);

    Var = (double *)calloc(P*P,sizeof(double));
    int t = rep->dim_X[0]-1; 

    for (int i=0;i<P*P;i++) 
        Var[i] = dlmStruc[t].C[i];
        
    int err = cholesky_decomp2vec(Var,P);
    if (err) {  // err = 1 means P is SPD
        err = rmvtvec(&(rep->delta[t*P]),Var,P,dlmStruc[t].n,dlmStruc[t].m,seed);
    }
    else {
        printf("error in DLM, var is not SPD\n");
        exit(0);
    }
    rep->d_Y[t] = rgamma(dlmStruc[t].n/2.,dlmStruc[t].d/2.,seed);     

    free(Var);
    
    // Backward sampling
    
    dlm_backward_sampling(rep,dlmStruc,(const int)P,seed);
    
    calWdelta(rep->Wdelta,rep->W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
    calculate_residuals(rep,P);
}

void DLM(REP *rep,int iter,unsigned long *seed) {
//    DLMtst2(rep,rep->dlmStruc,fdelta1,seed);
//    DLMtst2(rep,rep->dlmStruc,fdelta2,seed);
//    DLMtst2(rep,rep->dlmStruc,fP,seed);
    
//    free(rep->delta);
//        rep->delta = (double *)calloc(rep->dim_W[0]*rep->P,sizeof(double));        

    sampleDLM(rep,rep->dlmStruc,seed);
}

