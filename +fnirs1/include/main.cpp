#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "cholesky.h"
#include "randgen.h"
#include "fNIRS.h"

double M;
int PREC = 64;
sDLM *dlmStruc;
int maxP = 50;
int maxDLM;
int MAX_ITER = 7500;
int BURN_IN = 2500;
int nthd_sys = 1;
FILE *flog,*fseed;
int PPP;

int main (int argc, const char * argv[]) {
// command line arguments
//  1 - setup file name

    char *S;
	unsigned long *seed;
    POP *pop;
    FILE *fnthread;
    
    void delete_REP_Struc(REP *rep);
    void delete_SUB_Struc(SUB *sub);
    void delete_POP_Struc(POP *pop);
    void load_config_info(POP *,const char *,unsigned long *seed);
    void load_data_structs(POP *,int);
    void mcmc(POP *pop,unsigned long *seed);
    void compute_stats(POP *pop,const double cred_int,const int Niter);   
    
	M = exp(-PREC*log(2.0));  /* PREC should be set to the compiler precisions 32 or 64 bit */

    seed = (unsigned long *)calloc(3,sizeof(unsigned long));
    pop = (POP *)calloc(1,sizeof(POP));
    PPP = 5;  // starting AR degree 
    
 
    int iss = system("mkdir -p log"); 
    flog = fopen("./log/output.log","w");
 
 
    S = (char *)calloc(100,sizeof(char));
    iss = system("touch mnt_xyz");
    iss = system("uname -s >> mnt_xyz");
    fnthread = fopen("mnt_xyz","r");
    iss = fscanf(fnthread,"%s\n",S);
    if (!(iss = strcmp(S,"Darwin")))
        iss = system("sysctl -n hw.ncpu >> mnt_xyz");
    else if (!(iss = strcmp(S,"Linux")))
        iss = system("nproc --all >> mnt_xyz");
    else {
        fprintf(flog,"Kernel %s not supported",S);
        printf("Kernel %s not supported",S);
        exit(0);
    }
    rewind(fnthread);
    iss = fscanf(fnthread,"%s\n %d",S,&nthd_sys);
    iss = system("rm mnt_xyz");
    fclose(fnthread);
    

    fprintf(flog,"Kernel: %s, Number of Threads = %3d\n",S,nthd_sys);
    free(S);
    
    load_config_info(pop,argv[1],seed);
 
    load_data_structs(pop,PPP);
/*** CALL MCMC ***/
    fprintf(flog,"ENTERING MCMC\n");fflush(NULL);
    mcmc(pop,seed);
    fprintf(flog,"EXITING MCMC\n");fflush(NULL);
/*****************/

    compute_stats(pop,0.95,(const int)(MAX_ITER-BURN_IN));    
    for (int isub=0;isub<pop->N_SUBS;isub++) {
        SUB *sub = &(pop->sub[isub]);
        for (int irep=0;irep<sub->N_REPS;irep++) {
            REP *rep = &(sub->rep[irep]);
            delete_REP_Struc(rep);
        }
        delete_SUB_Struc(sub);
    }
    delete_POP_Struc(pop);
    free(pop);
 /*   for (int k=0;k<maxDLM;k++)
        free(dlmStruc[k].m);
    for (int k=0;k<maxDLM;k++)
        free(dlmStruc[k].C);
    free(dlmStruc);*/
 
    rewind(fseed);
	for (int i=0;i<3;i++)
		fprintf(fseed,"%lu ",seed[i]);
	fclose(fseed);
	free(seed);

    fclose(flog);
    
    printf("\n\nPlease see file Parameter_Estimates.log for parameter estimates and summaries\n");
    
    return 0;
}

