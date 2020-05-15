#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fNIRS.h"
#include "cholesky.h"
#include "randgen.h"

extern int maxP;
extern sDLM *dlmStruc;

double dlm_forward_filter_draw(REP *rep,sDLM *dlmStruc,const int Pmax, const int P, const double beta,const double delta,double *W)
{
    double n0,S0;
    double Q,*A;
    double f,e;
    double tden(double x,double mean,double var,double df);
    
    A = (double *)calloc(P,sizeof(double));
        
    for (int i=0;i<P*P;i++)
            dlmStruc[P-1].C[i] = 0;
    for (int i=0;i<P;i++)
            dlmStruc[P-1].C[i*P+i] = 1;
    S0 = 1;
    n0 = 1;            
    S0 = S0/n0;
    
    dlmStruc[P-1].n = n0;
    dlmStruc[P-1].d = dlmStruc[P-1].n*S0;
    dlmStruc[P-1].S = dlmStruc[P-1].d/dlmStruc[P-1].n;

    for (int i=0;i<P;i++)
        dlmStruc[P-1].m[i] = 0;
    
    double ll = 0;
    for (int t=P;t<rep->dim_X[0];t++) {
        for (int i=0;i<P;i++)
            dlmStruc[t].a[i] = dlmStruc[t-1].m[i];
             
        f = 0;
        for (int i=0;i<P;i++)
            f += W[t*P+i]*dlmStruc[t].a[i];
    
        
        for (int i=0;i<P;i++)
            for (int j=0;j<P;j++) 
                dlmStruc[t].R[i*P+j] = dlmStruc[t-1].C[i*P+j]/delta;
 
        Q = dlmStruc[t-1].S;
        for (int i=0;i<P;i++)
            for (int j=0;j<P;j++)
                Q += W[t*P+i]*dlmStruc[t].R[i*P+j]*W[t*P+j];
        
        for (int i=0;i<P;i++) {
            A[i] = 0;
            for (int j=0;j<P;j++)
                A[i] += dlmStruc[t].R[i*P+j]*W[t*P+j];
            A[i] /= Q; 
        }                
                           
        e = rep->residuals3[t] - f;        
        
        if (t >= Pmax)
            ll += tden(rep->residuals3[t],f,Q,dlmStruc[t].n);

        dlmStruc[t].n = dlmStruc[t].n*beta + 1;
        dlmStruc[t].d = beta*(dlmStruc[t].d) + dlmStruc[t-1].S*e*e/Q;
        dlmStruc[t].S = dlmStruc[t].d/dlmStruc[t].n;
                    
        double tmp = dlmStruc[t].S/dlmStruc[t-1].S;
        for (int i=0;i<P;i++) {
            for (int j=0;j<P;j++) {
                dlmStruc[t].C[i*P+j] = tmp*(dlmStruc[t].R[i*P+j] - A[i]*A[j]*Q);
            }
        }
                
        for (int i=0;i<P;i++)
            dlmStruc[t].m[i] = dlmStruc[t].a[i] + A[i]*e;
 
    }
    free(A);
    return ll;
}

void dlm_forward_filter(REP *rep,sDLM *dlmStruc)
{
    int P;
    double n0,S0;
    double Q,*A;
    double f,e;
 
    P = rep->P;  
    A = (double *)calloc(P,sizeof(double));
        
    for (int i=0;i<P*P;i++)
            dlmStruc[P-1].C[i] = 0;
    for (int i=0;i<P;i++)
            dlmStruc[P-1].C[i*P+i] = 1;
    S0 = 1;
    n0 = 1;            
    S0 = S0/n0;
    
    dlmStruc[P-1].n = n0;
    dlmStruc[P-1].d = n0*S0;
    dlmStruc[P-1].S = dlmStruc[P-1].d/dlmStruc[P-1].n;

    for (int i=0;i<P;i++)
        dlmStruc[P-1].m[i] = 0;
    
    
    for (int t=P;t<rep->dim_X[0];t++) {
        for (int i=0;i<P;i++)
            dlmStruc[t].a[i] = dlmStruc[t-1].m[i];
             
        f = 0;
        for (int i=0;i<P;i++)
            f += rep->W[t*P+i]*dlmStruc[t].a[i];
    
        
        for (int i=0;i<P;i++)
            for (int j=0;j<P;j++) 
                dlmStruc[t].R[i*P+j] = dlmStruc[t-1].C[i*P+j]/rep->df_delta1;
 
        Q = dlmStruc[t-1].S;
        for (int i=0;i<P;i++)
            for (int j=0;j<P;j++)
                Q += rep->W[t*P+i]*dlmStruc[t].R[i*P+j]*rep->W[t*P+j];
        
        for (int i=0;i<P;i++) {
            A[i] = 0;
            for (int j=0;j<P;j++)
                A[i] += dlmStruc[t].R[i*P+j]*rep->W[t*P+j];
            A[i] /= Q; 
        }                
                           
        e = rep->residuals3[t] - f;        
        
        dlmStruc[t].n = dlmStruc[t-1].n*rep->df_delta2 + 1;
        dlmStruc[t].d = rep->df_delta2*dlmStruc[t-1].d + dlmStruc[t-1].S*e*e/Q;
        dlmStruc[t].S = dlmStruc[t].d/dlmStruc[t].n;
                    
        double tmp = dlmStruc[t].S/dlmStruc[t-1].S;
        for (int i=0;i<P;i++) {
            for (int j=0;j<P;j++) {
                dlmStruc[t].C[i*P+j] = tmp*(dlmStruc[t].R[i*P+j] - A[i]*A[j]*Q);
            }
        }
                
        for (int i=0;i<P;i++)
            dlmStruc[t].m[i] = dlmStruc[t].a[i] + A[i]*e;
 
    }
    free(A);
}

void dlm_backward_sampling_delta1(REP *rep,sDLM *dlmStruc,const int P,unsigned long *seed)
{
 //   int n;
    double *m,*Var;
       
    m = (double *)calloc(P,sizeof(double));
    Var = (double *)calloc(P*P,sizeof(double));   

    for (int t=rep->dim_X[0]-2;t>=P;t--) {
        
        for (int i=0;i<P;i++)
            m[i] = dlmStruc[t].m[i] + rep->df_delta1*(rep->delta[(t+1)*P + i] - dlmStruc[t+1].a[i]);


        // compute mean
                
        double tmp = 1-rep->df_delta1;
        for (int i=0;i<P*P;i++)
            Var[i] = dlmStruc[t].C[i]*tmp;
                
               
        // draw new delta[t]
        int err = cholesky_decomp2vec(Var,P);
        if (err) {  // err = 1 means C is SPD
            err = rmvtvec(&(rep->delta[t*P]),Var,P,dlmStruc[t].n,m,seed);
        }
        else {
           printf("error in dlm_backward_filter, var is not SPD\n");
           exit(0);
        }
    }
    free(m);
    free(Var);
}

void dlm_backward_sampling(REP *rep,sDLM *dlmStruc,const int P,unsigned long *seed)
{
    double *m,*Var;
       
    m = (double *)calloc(P,sizeof(double));
    Var = (double *)calloc(P*P,sizeof(double));   

    for (int t=rep->dim_X[0]-2;t>=P;t--) {
        
        for (int i=0;i<P;i++)
            m[i] = dlmStruc[t].m[i] + rep->df_delta1*(rep->delta[(t+1)*P + i] - dlmStruc[t+1].a[i]);


        // compute mean
                
        double tmp = 1-rep->df_delta1;
        for (int i=0;i<P*P;i++)
            Var[i] = dlmStruc[t].C[i]*tmp;
                
               
        // draw new delta[t]
        int err = cholesky_decomp2vec(Var,P);
        if (err) {  // err = 1 means C is SPD
            err = rmvtvec(&(rep->delta[t*P]),Var,P,dlmStruc[t].n,m,seed);
        }
        else {
           printf("error in dlm_backward_filter, var is not SPD\n");
           exit(0);
        }
        
        // draw new prec[t]
        rep->d_Y[t] = rep->df_delta2*rep->d_Y[t+1];
        rep->d_Y[t] += rgamma(0.5*(1-rep->df_delta2)*dlmStruc[t].n,0.5*dlmStruc[t].d,seed);
    }
    free(m);
    free(Var);
}

double dlm_backward_sampling_delta2(REP *rep,sDLM *dlmStruc,const int P)
{
    double ll; 
    double dens_gamma(double,double,double);  
     
    for (int t=P;t<rep->dim_X[0];t++) {
        if (t == P)
            ll = dens_gamma(rep->d_Y[t],0.5*dlmStruc[t-1].n,0.5*dlmStruc[t-1].d);
        else
            ll += dens_gamma(rep->d_Y[t],0.5*dlmStruc[t-1].n,0.5*dlmStruc[t-1].d);
    }
/*    for (int t=rep->dim_X[0]-1;t>=P;t--) {
        if (t == rep->dim_X[0]-1)
            ll = dens_gamma(rep->d_Y[t],0.5*dlmStruc[t].n,0.5*dlmStruc[t].d);
        else
            ll += dens_gamma(rep->d_Y[t],0.5*(1-rep->df_delta2)*dlmStruc[t].n,0.5*dlmStruc[t].d);
    }*/
    
}

void mcmc_delta1(REP *rep,int iter,dlm_tst_flag flag,unsigned long *seed) {
    double n,d;
    double old_loglik,new_loglik;
    double f,e,S_old,*W;
    void calculate_res3(REP *rep);
    void calculate_residuals(REP *rep,int P);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    double tden(double x,double mean,double var,double df);
 //   double log_likelihood(REP *rep);
    void DLMdelta1(REP *rep,unsigned long *seed);        
   
    int P = rep->P;
    W = rep->W;

    calculate_res3(rep);
    calculate_residuals(rep,rep->P);

    int prop_P;
    double new_log_dens,old_log_dens;
    int sign;
    double u;

    new_log_dens = old_log_dens = 0;
   
    // Calculate old log likelihood
    
 //   old_loglik = log_likelihood(rep);

    double low,high,range;
    
   
    range = rep->prop_sd[0];
    low = rep->df_delta1 - range;
    low = (low < 0.8) ? 0.8:low;
    high = rep->df_delta1 + range;
    high = (high > 0.999999) ? 0.999999:high;
 
        
    double save_delta1 = rep->df_delta1;   
    rep->df_delta1 = runif_atob(seed,low,high);
    new_log_dens = 1/(high-low);
 
    low = rep->df_delta1 - range;
    low = (low < 0.8) ? 0.8:low;
    high = rep->df_delta1 + range;
    high = (high > 0.999999) ? 0.999999:high;

    old_log_dens = 1/(high-low);
 
            
    double *delta = (double *)calloc(rep->dim_W[0]*maxP,sizeof(double));
    for (int i=0;i<rep->dim_W[0]*maxP;i++)
        delta[i] = rep->delta[i];
        
    // Calculate new log likelihood
    
    DLMdelta1(rep,seed);
    calWdelta(rep->Wdelta,rep->W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
    calculate_residuals(rep,P);
    calculate_res3(rep);
    
//    new_loglik = log_likelihood(rep);

    if (log(kiss(seed)) < ((new_loglik - old_loglik) - (new_log_dens - old_log_dens))) {
        (rep->accept[flag])++;
    }
    else {
        rep->df_delta1 = save_delta1; 
        for (int i=0;i<rep->dim_W[0]*maxP;i++)
            rep->delta[i] = delta[i];
        calWdelta(rep->Wdelta,rep->W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
        calculate_residuals(rep,P);
        calculate_res3(rep);
    }
    (rep->attempt[flag])++;
    free(delta);       
}

void mcmc_delta2(REP *rep,int iter,dlm_tst_flag flag,unsigned long *seed) {
    double n,d;
    double old_loglik,new_loglik;
    double f,e,S_old;
    void calculate_residuals(REP *rep,int P);
    void calculate_res3(REP *rep);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    double tden(double x,double mean,double var,double df);
    double log_likelihood(REP *rep);
    double DLMdelta2(REP *rep);        
   

    calculate_res3(rep);
//    calculate_residuals(rep,rep->P);

    double new_log_prop,old_log_prop;

   
    // Calculate old log likelihood
    
    old_loglik = DLMdelta2(rep);

    double low,high,range;
    
    range = rep->prop_sd[1];
    low = rep->df_delta2 - range;
    low = (low < 0.1) ? 0.1:low;
    high = rep->df_delta2 + range;
    high = (high > 0.999999) ? 0.999999:high;
 
        
    double save_delta2 = rep->df_delta2;   
    rep->df_delta2 = runif_atob(seed,low,high);
    new_log_prop = log(1./(high-low));
 
    low = rep->df_delta2 - range;
    low = (low < 0.1) ? 0.1:low;
    high = rep->df_delta2 + range;
    high = (high > 0.999999) ? 0.999999:high;

    old_log_prop = log(1./(high-low));
 
            
/*    double *d_Y = (double *)calloc(rep->dim_X[0],sizeof(double));
    for (int i=0;i<rep->dim_X[0];i++)
        d_Y[i] = rep->d_Y[i];
*/        
    // Calculate new log likelihood
    
    new_loglik = DLMdelta2(rep);
    
//    new_loglik = log_likelihood_delta2(rep);

    if (log(kiss(seed)) < ((new_loglik - old_loglik) - (new_log_prop - old_log_prop))) {
        (rep->accept[flag])++;
    }
    else {
        rep->df_delta2 = save_delta2; 
//        for (int i=0;i<rep->dim_X[0];i++)
//            rep->d_Y[i] = d_Y[i];
    }
    (rep->attempt[flag])++;
//    free(d_Y);       
}


void DLMtst(REP *rep,int iter,dlm_tst_flag flag,unsigned long *seed) {
    double old_loglik,new_loglik;
    double f,e,S_old,*W;
    void calculate_res3(REP *rep);
    void calculate_residuals(REP *rep,int P);
    void calW(double *W,double *Y,double *Xb,double *Ve,const int nrow,const int ncol,const int P);
    double tden(double x,double mean,double var,double df);
    double log_likelihood(REP *rep);
      
    int P = rep->P;
    W = rep->W;

    calculate_res3(rep);

    int prop_P;
    double new_log_dens,old_log_dens;
    int sign;
    double u;
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
        new_log_dens = log(PDF1[prop_P-1]);
 
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
 
        old_log_dens = log(PDF1[P-1]);
  
        free(PDF1);
        free(CDF1);
    }
    else {
        new_log_dens = old_log_dens = 0;
        prop_P = P;
    }


    int Pmax = (P > prop_P) ? P:prop_P;
    
    // Calculate old log likelihood
    
    old_loglik = dlm_forward_filter_draw(rep,dlmStruc,(const int)Pmax,(const int)P,(const double)rep->df_delta2,(const double)rep->df_delta1,W);

    double prop_delta,prop_beta;
    double low,high;
    double range = rep->prop_sd[0];
    low = rep->df_delta1 - range;
    low = (low < 0.8) ? 0.8:low;
    high = rep->df_delta1 + range;
    high = (high > 0.99999) ? 0.99999:high;
    
    prop_delta = runif_atob(seed,low,high);
    
    range = rep->prop_sd[1];
    low = rep->df_delta2 - range;
    low = (low < 0.7) ? 0.7:low;
    high = rep->df_delta2 + range;
    high = (high > 0.99999) ? 0.99999:high;
    
    prop_beta = runif_atob(seed,low,high);
     
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
    rep->df_delta2 = prop_beta;
    new_loglik = dlm_forward_filter_draw(rep,dlmStruc,(const int)Pmax,(const int)prop_P,(const double)prop_beta,(const double)prop_delta,W);

    if (log(kiss(seed)) < ((new_loglik - old_loglik) - (new_log_dens - old_log_dens))) {
        rep->df_delta1 = prop_delta;
        rep->df_delta2 = prop_beta;
        rep->P = prop_P;
        rep->dim_W[1] = rep->P;
        if (flag < 2)
            (rep->accept[flag])++;
    }
    if (flag < 2)
        (rep->attempt[flag])++;
        
    calW(rep->W,rep->Y,rep->Xbeta,rep->Veta,(const int)rep->dim_W[0],(const int)rep->dim_W[1],rep->P);
}

void DLM(REP *rep,unsigned long *seed) {
    double *Var;
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calculate_residuals(REP *rep,int P);
    void calculate_res3(REP *);
  
    int P = rep->P;
 //   calculate_res3(rep,P);

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

void DLMdelta1(REP *rep,unsigned long *seed) {
    double *Var;
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calculate_residuals(REP *rep,int P);
    void calculate_res3(REP *);
  
    int P = rep->P;
 //   calculate_res3(rep,P);

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

    free(Var);
    
    // Backward sampling
    
    dlm_backward_sampling_delta1(rep,dlmStruc,(const int)P,seed);
    
//    calWdelta(rep->Wdelta,rep->W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
//    calculate_residuals(rep,P);
}

double DLMdelta2(REP *rep) {
    double *Var;
    void calWdelta(double *Ax,double *A,double *x,const int nrow,const int ncol);
    void calculate_residuals(REP *rep,int P);
    void calculate_res3(REP *);
  
 //   int P = rep->P;
 //   calculate_res3(rep,P);

    // Forward_filter
  
    dlm_forward_filter(rep,dlmStruc);

 //   Var = (double *)calloc(P*P,sizeof(double));
 //   int t = rep->dim_X[0]-1; 

/*    for (int i=0;i<P*P;i++) 
        Var[i] = dlmStruc[t].C[i];
        
    int err = cholesky_decomp2vec(Var,P);
    if (err) {  // err = 1 means P is SPD
        err = rmvtvec(&(rep->delta[t*P]),Var,P,n,dlmStruc[t].m,seed);
    }
    else {
        printf("error in DLM, var is not SPD\n");
        exit(0);
    }*/
 //   rep->d_Y[t] = rgamma(dlmStruc[t].n/2.,dlmStruc[t].d/2.,seed);     

 //   free(Var);
    
    // Backward sampling
    
    double ll = dlm_backward_sampling_delta2(rep,dlmStruc,(const int)rep->P);
    return ll;
//    calWdelta(rep->Wdelta,rep->W,rep->delta,(const int)rep->dim_W[0],(const int)rep->dim_W[1]);
//    calculate_residuals(rep,P);
}

