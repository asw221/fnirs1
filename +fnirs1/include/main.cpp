
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
int maxP = 100;
int MAX_ITER = 7500;
int BURN_IN = 2500;
FILE *flog,*fseed;

int main (int argc, const char * argv[]) {
  // command line arguments
  //  1 - setup file name

  int PPP;
  unsigned long *seed;
  POP *pop;
    
  void load_config_info(POP *,const char *,unsigned long *seed);
  void load_data_structs(POP *,int);
  void mcmc(POP *pop,unsigned long *seed);
  void compute_stats(POP *pop,const double cred_int,const int Niter);   
    
  M = exp(-PREC*log(2.0));  /* PREC should be set to the compiler precisions 32 or 64 bit */

  seed = (unsigned long *)calloc(3,sizeof(unsigned long));
  pop = (POP *)calloc(1,sizeof(POP));
  PPP = 15;  // starting AR degree 
    
  int iss = system("mkdir -p log"); 
  flog = fopen("./log/output.log","w");

  load_config_info(pop,argv[1],seed);
 
  load_data_structs(pop,PPP);

  /*** CALL MCMC ***/
  fprintf(flog,"ENTERING MCMC\n");fflush(NULL);
  mcmc(pop,seed);
  fprintf(flog,"EXITING MCMC\n");fflush(NULL);
  /*****************/

  compute_stats(pop,0.95,(const int)(MAX_ITER-BURN_IN));    

  rewind(fseed);
  for (int i=0;i<3;i++)
    fprintf(fseed,"%lu ",seed[i]);
  fclose(fseed);
  free(seed);

  fclose(flog);
    
  printf("\n\nPlease see file Parameter_Estimates.log for parameter estimates and summaries\n");
    
  return 0;
}

