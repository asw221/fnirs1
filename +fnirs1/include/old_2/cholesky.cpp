#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "randgen.h"
#include "cholesky.h"

/***************************************************************************/
/*                                                                         */
/* cholesky_decomp                                                         */
/*                                                                         */
/* This routine takes a symmetric positive definite matrix and performs    */
/* a cholesky decomposition.  It replaces the lower triangular part of     */
/* A with G.                                                               */
/*                                                                         */
/* The cholesky decomposition, decomposes A into A = GG'.  This version of */
/* the algorithm is an outer product (rank-1) update of the principal      */
/* submatrices of A.                                                       */
/*                                                                         */
/* See "Matrix Computations" by Golub and Van Loan.  2nd edition, Johns    */
/* Hopkins University Press, 1989. p. 423ff                                */
/*                                                                         */
/* Input  A - a SPD matrix                                                 */
/*        num_col - size of A                                              */
/*                                                                         */
/* Returns 1 upon success, 0 if A not SPD                                  */
/***************************************************************************/
int cholesky_decomp(double **A, int num_col)
{
  int k,i,j;

    for (k=0;k<num_col;k++) {
        if (A[k][k] <= 0) {
            printf("cholesky_decomp: error matrix A not SPD %g %d\n",A[k][k],k);
            return 0;
        }
        
        A[k][k] = sqrt(A[k][k]);
        
        for (j=k+1;j<num_col;j++)
            A[j][k] /= A[k][k];
        
        for (j=k+1;j<num_col;j++)
            for (i=j;i<num_col;i++)
                A[i][j] = A[i][j] - A[i][k] * A[j][k];
    }
    return 1;
    
/* this is the GAXPY update */
/*    double tmp;
    for (j=0;j<num_col;j++) {
        if (j > 0) {
            for (i=j;i<num_col;i++) {
                tmp = 0;
                for (k=0;k<j;k++)
                    tmp += A[i][k]*A[j][k];
                A[i][j] -= tmp;
            }
            if (A[j][j] <=0) {
                printf("cholesky_decomp: error matrix A not SPD %lf %d\n",A[k][k],k);
                return 0;
            }
        }
        tmp = sqrt(A[j][j]);
        for (i=j;i<num_col;i++)
            A[i][j] /= tmp;
    }
    return 1;*/
}

int cholesky_decomp2(double **A, int num_col)
{
    // replaces upper triangular part of A with G' where A = GG'
    int k,i,j;
    
    for (k=0;k<num_col;k++) {
        double *tmp2 = A[k];
        if (tmp2[k] <= 0) {
            printf("cholesky_decomp2: error matrix A not SPD %g %d\n",A[k][k],k);
            return 0;
        }
        
        tmp2[k] = sqrt(tmp2[k]);
        double tmp = 1./tmp2[k];
        for (j=k+1;j<num_col;j++)
            tmp2[j] *= tmp;//A[k][k];
        
        for (j=k+1;j<num_col;j++) {
            double *tmp3 = A[j];
            for (i=j;i<num_col;i++)
                tmp3[i] -= tmp2[i] * tmp2[j];
        }
    }
    return 1;
}

int cholesky_decomp2vec(double *A, int num_col)
{
    // replaces upper triangular part of A with G' where A = GG'
    int k,i,j;
    
    for (k=0;k<num_col;k++) {
        double *tmp2 = &A[k*num_col];
        if (tmp2[k] <= 0) {
            printf("cholesky_decomp2vec: error matrix A not SPD %g %d\n",A[k+num_col*k],k);
            return 0;
        }
        
        tmp2[k] = sqrt(tmp2[k]);
        double tmp = 1./tmp2[k];
        for (j=k+1;j<num_col;j++)
            tmp2[j] *= tmp;//A[k][k];
        
        for (j=k+1;j<num_col;j++) {
            double *tmp3 = &A[j*num_col];
            for (i=j;i<num_col;i++)
                tmp3[i] -= tmp2[i] * tmp2[j];
        }
    }
    return 1;
}

void cholesky_invert2(int len,double **G)
{
    /* takes G' from GG' = A and computes A inverse */
    int i,j,k;
    double **INV;

    int forward_substitution2(double **U,double *b,int n);
    int cholesky_backsub2(double **G,double *b,int n);

    INV = (double **)calloc(len,sizeof(double *));
    for (i=0;i<len;i++)
        INV[i] = (double *)calloc(len,sizeof(double));
    for (i=0;i<len;i++)
        INV[i][i] = 1;
    
    for (k=0;k<len;k++) {
        forward_substitution2(G,INV[k],len);
        cholesky_backsub2(G,INV[k],len);
    }
    
    for (i=0;i<len;i++)
        for (j=0;j<len;j++)
            G[i][j] = INV[i][j];
 
    for (i=0;i<len;i++)
        free(INV[i]);
    free(INV);
}


void cholesky_invert(int len,double **G)
{
/* takes G from GG' = A and computes A inverse */
  int i,j,k;
  double temp,**INV;

  INV = (double **)calloc(len,sizeof(double *));
  for (i=0;i<len;i++)
    INV[i] = (double *)calloc(len,sizeof(double));
  for (i=0;i<len;i++)
    INV[i][i] = 1;
  
  for (k=0;k<len;k++) {
    INV[k][0] /= G[0][0];
    for (i=1;i<len;i++) {
      temp = 0.0;
      for (j=0;j<i;j++)
	temp += G[i][j]*INV[k][j];
      INV[k][i] = (INV[k][i] - temp)/G[i][i];
    }
    INV[k][len-1] /= G[len-1][len-1];
    for (i=len-2;i>=0;i--) {
      temp = 0.0;
      for (j=i+1;j<len;j++)
	temp += G[j][i]*INV[k][j];
      INV[k][i] = (INV[k][i] - temp)/G[i][i];
    }
  }
  for (i=0;i<len;i++)
    for (j=0;j<len;j++)
      G[i][j] = INV[i][j];
  for (i=0;i<len;i++)
    free(INV[i]);
  free(INV);
}

double Determinant(double **A,int dim)
{
	double det;

	switch (dim) {
		case 2:
			det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
			break;
		case 3:
			det = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] -
				  A[2][0]*A[1][1]*A[0][2] - A[2][1]*A[1][2]*A[0][0] - A[2][2]*A[1][0]*A[0][1]; 
			break;
		default:
			det = 1;
			break;
	}
	return det;
}

int MatrixInverse(double **A,int dim)
{
	int i,j,flag;
	double det,tmp,**B;
	
	switch (dim) {
	case 2:
		det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
		if (fabs(det) > 1e-10) {
			det = 1.0/det;
			tmp = A[0][0];
			A[0][0] = A[1][1]*det;
			A[1][1] = tmp*det;
			A[0][1] = -A[0][1]*det;
			A[1][0] = -A[1][0]*det;
			flag = 1;
		}
		else
			flag = 0;
		break;
	case 3:
		det = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] -
			  A[2][0]*A[1][1]*A[0][2] - A[2][1]*A[1][2]*A[0][0] - A[2][2]*A[1][0]*A[0][1]; 
		if (fabs(det) > 1e-10) {
			B = (double **)calloc(dim,sizeof(double *));
			for (i=0;i<dim;i++)
				B[i] = (double *)calloc(dim,sizeof(double));
			det = 1.0/det;
			B[0][0] = (A[1][1]*A[2][2] - A[2][1]*A[1][2])*det;
			B[0][1] = -(A[0][1]*A[2][2] - A[2][1]*A[0][2])*det;
			B[0][2] = (A[0][1]*A[1][2] - A[1][1]*A[0][2])*det;
			B[1][0] = -(A[1][0]*A[2][2] - A[2][0]*A[1][2])*det;
			B[1][1] = (A[0][0]*A[2][2] - A[2][0]*A[0][2])*det;
			B[1][2] = -(A[0][0]*A[1][2] - A[1][0]*A[0][2])*det;
			B[2][0] = (A[1][0]*A[2][1] - A[2][0]*A[1][1])*det;
			B[2][1] = -(A[0][0]*A[2][1] - A[2][0]*A[0][1])*det;
			B[2][2] = (A[0][0]*A[1][1] - A[1][0]*A[0][1])*det;
			for (i=0;i<dim;i++)
				for (j=0;j<dim;j++)
					A[i][j] = B[i][j];
			for (i=0;i<dim;i++)
				free(B[i]);
			free(B);
			flag = 1;
		}
		else 
			flag = 0;
		break;
	default:
		flag = 0;
		break;
	}
	return flag;
}

int banded_cholesky_decomp(double **A,int num_col,int p)
{
  int i,j,k,lk,lambda;
  double tmp;
  for (j=0;j<num_col;j++) {
    lk = (0 > j-p+1) ? 0:j-p+1;
    for (k=lk;k<j;k++) {
      lambda = ((lambda=k+p) < num_col) ? lambda : num_col;  
      for (i=j;i<lambda;i++) {
	tmp = A[j][k]*A[i][k];
	A[i][j] -= tmp;
      }
    }
    lambda = ((lambda=j+p) < num_col) ? lambda: num_col;
    if (A[j][j] <= 0)
      return 0;
    tmp = sqrt(A[j][j]);
    for (i=j;i<lambda;i++)
      A[i][j] /= tmp;
  }
  return 1;
}

void banded_cholesky_invert(int len,double **H,int p)
{
/* takes G from GG' = A and computes A inverse  where A is a banded matrix with bandwidth p*/
  int i,j,k,lj;
  double temp,**INV;

  INV = (double **)calloc(len,sizeof(double *)); 
  for (i=0;i<len;i++)
    INV[i] = (double *)calloc(len,sizeof(double));
  for (i=0;i<len;i++)
    INV[i][i] = 1;
  
  for (k=0;k<len;k++) {
    INV[k][0] /= H[0][0];
    for (i=1;i<len;i++) {
      temp = 0.0;
      lj = (0>(lj=i-p+1)) ? 0:lj;
      for (j=lj;j<i;j++)
	temp += H[i][j]*INV[k][j];
      INV[k][i] = (INV[k][i] - temp)/H[i][i];
    }
    INV[k][len-1] /= H[len-1][len-1];
    for (i=len-2;i>=0;i--) {
      temp = 0.0;
      lj = (len<(lj=i+p)) ? len:lj;
      for (j=i+1;j<lj;j++)
	temp += H[j][i]*INV[k][j];
      INV[k][i] = (INV[k][i] - temp)/H[i][i];
    }
  }
  for (i=0;i<len;i++)
    for (j=0;j<len;j++)
      H[i][j] = INV[i][j];
  for (i=0;i<len;i++)
    free(INV[i]);
  free(INV);
}

int rmvnorm(double *result,double **A,int size_A,double *mean,unsigned long *seed,int flag)
{
	int i,j;
	double *runiv;
    // A is upper triangular
	runiv = (double *)calloc(size_A,sizeof(double));
	j = 1;
	if (!flag)
		j = cholesky_decomp2(A,size_A);
	if (j) {
		for (i=0;i<size_A;i++) 
			runiv[i] = snorm(seed);
		for (i=0;i<size_A;i++) {
			result[i] = 0;
			for (j=0;j<=i;j++)
				result[i] += A[j][i]*runiv[j];  // uses A'
			result[i] += mean[i];
		}
		free(runiv);
		return 1;
	}
	else {
		free(runiv);
		return 0;
	}
}

int rmvnorm2(double *result,double **A,int size_A,double *mean,unsigned long *seed,int flag)
{
    // A is precision matrix
	int i,j;
	double *runiv;

	runiv = (double *)calloc(size_A,sizeof(double));
	j = 1;
	if (!flag)
		j = cholesky_decomp(A,size_A);
	if (j) {
		for (i=0;i<size_A;i++) 
			runiv[i] = snorm(seed);
		cholesky_backsub(A,runiv,size_A);
		for (i=0;i<size_A;i++) 
			result[i] = mean[i] + runiv[i];
		free(runiv);
		return 1;
	}
	else {
		free(runiv);
		return 0;
	}
}

int rmvnorm3(double *result,double **A,int size_A,double *mean,unsigned long *seed,int flag)
{
    // A is precision matrix
    int i,j;
    double *runiv;
    
    runiv = (double *)calloc(size_A,sizeof(double));
    j = 1;
    if (!flag)
        j = cholesky_decomp2(A,size_A);
    if (j) {
        for (i=0;i<size_A;i++)
            runiv[i] = snorm(seed);
        cholesky_backsub2(A,runiv,size_A);
        for (i=0;i<size_A;i++)
            result[i] = mean[i] + runiv[i];
        free(runiv);
        return 1;
    }
    else {
        free(runiv);
        return 0;
    }
}

int rmvnorm3vec(double *result,double *A,int size_A,double *mean,unsigned long *seed,int flag)
{
    // A is precision matrix
    int i,j;
    double *runiv;
    
    runiv = (double *)calloc(size_A,sizeof(double));
    j = 1;
    if (!flag)
        j = cholesky_decomp2vec(A,size_A);
    if (j) {
        for (i=0;i<size_A;i++)
            runiv[i] = snorm(seed);
        cholesky_backsub2vec(A,runiv,size_A);
        for (i=0;i<size_A;i++)
            result[i] = mean[i] + runiv[i];
        free(runiv);
        return 1;
    }
    else {
        free(runiv);
        return 0;
    }
}


int rmvt(double *result,double **A,int size_A,double df,double *mean,unsigned long *seed) 
{
    //  A is upper triangular matrix from Cholesky decomp
    double x;
    double *ran_norm;
    ran_norm = (double *)calloc(size_A,sizeof(double));
    for (int i=0;i<size_A;i++)
        ran_norm[i] = snorm(seed);
    x = rgamma(0.5*df,0.5,seed);
    
    for (int i=0;i<size_A;i++) {
        result[i] = 0;
	    for (int j=0;j<=i;j++)
		    result[i] += A[j][i]*ran_norm[j];  // uses A'
		result[i] = mean[i] + result[i]*sqrt(df/x);
	}

	free(ran_norm);
	return 1;
}

int rmvtvec(double *result,double *A,int size_A,double df,double *mean,unsigned long *seed)
{
    //  A is upper triangular matrix from Cholesky decomp
    double x;
    double *ran_norm,alpha;
    double beta = 0.5;
    ran_norm = (double *)calloc(size_A,sizeof(double));
    for (int i=0;i<size_A;i++)
        ran_norm[i] = snorm(seed);
    alpha = 0.5*df;
    x = rgamma(alpha,beta,seed);
    if (isnan(x))
        printf("in rmvtvec, rgamma returned nan\n");
    while (x < 0) {
        printf("in rmvtvec, rgamma returned a negative number %.20lf alpha = %.20lf beta = %.20lf %lu %lu %lu \n",x,alpha,beta,seed[0],seed[1],seed[2]);
        x = rgamma(alpha,beta,seed);
     }   
    for (int i=0;i<size_A;i++) {
        result[i] = 0;
        for (int j=0;j<=i;j++)
            result[i] += A[j*size_A + i]*ran_norm[j];  // uses A'
        result[i] = mean[i] + result[i]*sqrt(df/x);
    }
    
    free(ran_norm);
    return 1;
}

int rwishart(double **result,double **S,int size_S,int df,unsigned long *seed)
{
	int i,j,k,flag;
	double *x,*zero;
	
	if ((double)df < (double)size_S)
		return 0;
	
	flag = 0;
	x = (double *)calloc(size_S,sizeof(double));
	zero = (double *)calloc(size_S,sizeof(double));
	
	for (i=0;i<size_S;i++) 
		for (j=0;j<size_S;j++) 
			result[i][j] = 0;
	for (k=0;k<df;k++) {
		if (k > 0)
			flag = 1;
		if (rmvnorm(x,S,size_S,zero,seed,flag)) {
			for (i=0;i<size_S;i++) 
				for (j=0;j<size_S;j++) 
					result[i][j] += x[i]*x[j];
		}
		else {
			free(zero);
			free(x);
			return 0;
		}
	}
	free(zero);
	free(x);
	return 1;
}

/*double xAx(double **A,double *x,int dim)
{
	int i,j;
	double ans,*tmp;

	tmp = (double *)calloc(dim,sizeof(double));
	for (i=0;i<dim;i++)
		for (j=0;j<dim;j++)
			tmp[i] += A[i][j]*x[j];
	ans = 0;
	for (i=0;i<dim;i++)
		ans += x[i]*tmp[i];
	
	free(tmp);
	return ans;
}*/

/***************************************************************************/
/*                                                                         */
/* forward_substitution                                                    */
/*                                                                         */
/* Given a lower triangular matrix L and a vector b, forward_substition    */
/* solves the equation Lx = b for x and replaces b with x.                 */
/*                                                                         */
/* Input  L - lower triangular matrix (n by n)                             */
/*        b - a n-dim vector                                               */
/*        n - number of columns in L                                       */
/*                                                                         */
/* Returns 1 on success                                                    */
/*                                                                         */
/*                                                                         */
/* backward_substitution                                                   */
/*                                                                         */
/* Given an upper triangular matrix U and a vector b, backward_substitution*/
/* solves the equation Ux = b for x and replaces b with x.                 */
/*                                                                         */
/* Input  U - upper triangular matrix (n by n)                             */
/*        b - a n-dim vector                                               */
/*        n - number of columns in U                                       */
/*                                                                         */
/* Returns 1 on success                                                    */
/*                                                                         */
/*                                                                         */
/* cholesky_backsub                                                        */
/*                                                                         */
/* Given an cholesky lower triangular matrix G and a vector b,             */
/* cholesky_backsub solves the equation G'x = b for x and replaces b       */
/* with x.                                                                 */
/*                                                                         */
/* Input  G - lower triangular cholesky matrix (n by n)                    */
/*        b - a n-dim vector                                               */
/*        n - number of columns in G                                       */
/*                                                                         */
/* note: use forward_substitution if you want to solve Gx = b.             */
/*                                                                         */
/* Returns 1 on success                                                    */
/*                                                                         */
/***************************************************************************/

int forward_substitution(double **L,double *b,int n)
{
  int i,j;
  double temp;

    b[0] /= L[0][0];
    
    for (i=1;i<n;i++)
    {
        temp = 0.0;
        
        for (j=0;j<i;j++)
            temp += L[i][j]*b[j];
        
        b[i] = (b[i] - temp)/L[i][i];
    }
    
    return 1;
}

int forward_substitution2(double **U,double *b,int n)
{
    int i,j;
    double temp;
    
    b[0] /= U[0][0];
    
    for (i=1;i<n;i++)
    {
        temp = 0.0;
        
        for (j=0;j<i;j++)
            temp += U[j][i]*b[j];
        
        b[i] = (b[i] - temp)/U[i][i];
    }
    
    return 1;
}

int forward_substitution2vec(double *U,double *b,int n)
{
    int i,j;
    double temp;
    
    b[0] /= U[0];
    
    for (i=1;i<n;i++)
    {
        temp = 0.0;
        
        for (j=0;j<i;j++)
            temp += U[j*n+i]*b[j];
        
        b[i] = (b[i] - temp)/U[i+n*i];
    }
    
    return 1;
}

int backward_substitution(double **U,double *b,int n)
{
  int i,j;
  double temp;

  b[n-1] /= U[n-1][n-1];

  for (i=n-2;i>=0;i--)
    {
      temp = 0.0;

      for (j=i+1;j<n;j++)
	temp += U[i][j]*b[j];

      b[i] = (b[i] - temp)/U[i][i];
    }

  return 1;
}


int cholesky_backsub(double **G,double *b,int n)
{
    int i,j;
    double temp;
    
    b[n-1] /= G[n-1][n-1];
    
    for (i=n-2;i>=0;i--)
    {
        temp = 0.0;
        
        for (j=i+1;j<n;j++)
            temp += G[j][i]*b[j];
        
        b[i] = (b[i] - temp)/G[i][i];
    }
    
    return 1;
}

int cholesky_backsub2(double **G,double *b,int n)
{
    int i,j;
    double temp;
    
    b[n-1] /= G[n-1][n-1];
    
    for (i=n-2;i>=0;i--)
    {
        temp = 0.0;
        
        for (j=i+1;j<n;j++)
            temp += G[i][j]*b[j];
        
        b[i] = (b[i] - temp)/G[i][i];
    }
    
    return 1;
}

int cholesky_backsub2vec(double *G,double *b,int n)
{
    int i,j;
    double temp;
    
    b[n-1] /= G[(n-1)+(n-1)*n];
    
    for (i=n-2;i>=0;i--)
    {
        temp = 0.0;
        
        for (j=i+1;j<n;j++)
            temp += G[i*n+j]*b[j];
        
        b[i] = (b[i] - temp)/G[i*n+i];
    }
    
    return 1;
}


