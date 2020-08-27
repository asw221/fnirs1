#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void kernel_reg(double *f,double *Y,double *X,double lambda,int N)
{
    // tricube kernel
    int lambda_int = (int)floor(lambda);
    double lambda_cubed = lambda*lambda*lambda;
    double factor = 70./81.;
    
    double *K = (double *)calloc(lambda_int+1,sizeof(double));
    for (int i=0;i<lambda_int+1;i++) {
        double tmp = 1 - i*i*i/lambda_cubed;
        K[i] = factor * tmp*tmp*tmp;
    }

    int lower,upper;
    double denom;
    for (int t=0;t<N;t++) {
        lower = (0 >= (X[t] - lambda_int)) ? 0:(X[t] - lambda_int);
        upper = ((X[t] + lambda_int) <= (N-1)) ? (X[t] + lambda_int):(N-1);
        
        denom = 0;
        for (int i=lower;i<=upper;i++)
            denom += K[abs(i-(int)X[t])];
        
        f[t] = 0;
        for (int i=lower;i<=upper;i++)
            f[t] += K[abs(i-(int)X[t])]*Y[i];
        f[t] /= denom;
    }
}
