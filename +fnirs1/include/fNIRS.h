enum dlm_tst_flag {fdelta1,fdelta2,fP};


typedef struct dlm{
    double S;
    double n;
    double d;
    double *m;
    double *C;
} sDLM;

typedef struct replication{
    char *dataname;
    char *designname;
    double **design;
    int N;
    int    dim_design[2];
    double mhrf;
    double K;
    double *pi_Y;
    double pi_eta;
    double pi_delta;
    int *dim_X;
    int *dim_W;
    int *dim_V;
    int P;
    double df_delta1;
    double df_delta2;
    int     *dY;
    double *d_Y;
    int     d_Ycnt;
    double pi;
    double *md_Y;
    double *sd_Y;
    double *Y;
    double *precY;
    double LY;
    double precYstart;
    double *eta;
    double preceta;
    double *delta;
    double precdelta;
    double *beta;
    double *Xbeta;
    double *mXbeta;
    double *Wdelta;
    double *mWdelta;
    double *Wm;
    double *X;
    double *W;
    double *V;
    double *Veta;
    double *mVeta;
    sDLM *dlmStruc;
    
    double *mdelta;
    double *mdelta2; 
    double *std_res;
    double *mean_res;
    double *mean_d_Y;
    double *mean_fit;
    double md1;
    double md2;
    double *residuals;
    double *residuals1;
    double *residuals2;
    double *residuals3;
    double *residuals4;
    double *residuals5;
        
    double *prop_sd;
    int *attempt;
    int *accept;
    int maxsteps;
    int nKnots;
    double *knots;
    
    FILE *fout_beta;
    FILE *fout_eta;
    FILE *fout_dlm;
    FILE *fout_res;
    FILE *fout_stdres;
    FILE *fout_delta;
    FILE *fout_nknots;
    FILE *fout_knots;
    FILE *fout_stdev;
    FILE *fout_wdelta;
    FILE *fout_veta;
    FILE *fout_fit;
    FILE *fout_X;
    FILE *fout_Y;
    FILE *fout_HRF;
    FILE *fout_Xbeta;
} REP;


typedef struct subdata{
    int N_REPS;
    double freq;
    double subsampled_freq;
    REP *rep;
    double *beta;  // subject level effects  (Nb*Ns)
    int *dim_X;
    double *X; // (Nb*Ns)x(Ncov*Nb*Ns)  covariate information for each subject
    
    FILE *fout_beta;
} SUB;


typedef struct popdata{
    int GRP;
    int N_SUBS;
    int No_replicates;
    int N_RE;  // number of random effects
    int N_Fixed; // number of Fixed effects
    int Ncov; // number of covariate
    int Nb;     // number of basis in the HRF
    int Ns;      // number of stimuli (or events)
    int non_parm; // canonical HRF = 0, non-parm B-SPlINE = 1
    int knots; // number of bspline bases
    double ED;
    SUB *sub;
    double *beta;  // population level effects, (Ncov*Nb*Ns)
    double *re_rep_prec; //replicate level random effects
    double *re_prec; // sub level random effects, (Ns)x(Ns) diagonal
    char **covnames;
    FILE *fout_reprec;
    FILE *fout_re_rep_prec;
    FILE *fout_beta;
} POP;
