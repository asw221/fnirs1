int rwishart(double **r,double **,int,int,unsigned long *);
int rmvnorm(double *result,double **A,int size_A,double *mean,unsigned long *seed,int flag);
int rmvnorm2(double *result,double **A,int size_A,double *mean,unsigned long *seed,int flag);
int rmvnorm3(double *result,double **A,int size_A,double *mean,unsigned long *seed,int flag);
int rmvnorm3vec(double *result,double *A,int size_A,double *mean,unsigned long *seed,int flag);
int cholesky_decomp(double **A, int num_col);
int cholesky_decomp2(double **A, int num_col);
int cholesky_decomp2vec(double *A, int num_col);
void cholesky_invert(int len,double **G);
void cholesky_invert2(int len,double **G);
double Determinant(double **A,int dim);
int MatrixInverse(double **A,int dim);
int banded_cholesky_decomp(double **A,int num_col,int p);
void banded_cholesky_invert(int len,double **H,int p);
int forward_substitution(double **L,double *b,int n);
int forward_substitution2(double **L,double *b,int n);
int forward_substitution2vec(double *L,double *b,int n);
int backward_substitution(double **U,double *b,int n);
int cholesky_backsub(double **G,double *b,int n);
int cholesky_backsub2(double **G,double *b,int n);
int cholesky_backsub2vec(double *G,double *b,int n);
int rmvt(double *result,double **A,int size_A,double df,double *mean,unsigned long *seed);
int rmvtvec(double *result,double *A,int size_A,double df,double *mean,unsigned long *seed);






