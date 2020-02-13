
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "cholesky.h"
#include "fNIRS.h"
#include "MatlabDataArray.hpp"
#include "matlabError.h"
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "printToMatlab.h"
#include "randgen.h"



/*
 * MATLAB interface to Dr. Johnson's DLM programs. Not intended to be called
 * directly by end-users. Instead see the function/M-file fnirs1.dlm
*/


double M;
int PREC = 64;
sDLM *dlmStruc;
int maxP = 100;
int MAX_ITER = 7500;
int BURN_IN = 2500;
FILE *flog,*fseed;



void load_config_info(POP *,const char *,unsigned long *seed);
void load_data_structs(POP *,int);
void mcmc(POP *pop,unsigned long *seed);
void compute_stats(POP *pop,const double cred_int,const int Niter);



class MexFunction : public matlab::mex::Function 
{
public:
  void operator() (
    matlab::mex::ArgumentList output, 
    matlab::mex::ArgumentList input
  ) {
    matlab::data::CharArray setupFileChar = input[0];
    
    int PPP;
    unsigned long *seed;
    POP *pop;
    
    
    M = exp(-PREC*log(2.0));  /* PREC should be set to the compiler precisions 32 or 64 bit */

    seed = (unsigned long *)calloc(3,sizeof(unsigned long));
    pop = (POP *)calloc(1,sizeof(POP));
    PPP = 15;  // starting AR degree 
    
    int iss = system("mkdir -p log"); 
    flog = fopen("./log/output.log","w");

    // load_config_info(pop,argv[1],seed);
    try {
      load_config_info(pop, setupFileChar.toAscii().c_str(), seed);
    }
    catch (...) {
      matlabError("Error loading model configuration information. Double-check setup file");
    }

    try {
      load_data_structs(pop,PPP);
    }
    catch (...) {
      matlabError("Error loading data.");
    }

    /*** CALL MCMC ***/
    try {
      printToMatlab("Starting MCMC...\n");
      fprintf(flog,"ENTERING MCMC\n");fflush(NULL);
      mcmc(pop,seed);
      printToMatlab("Ending MCMC\n");
      fprintf(flog,"EXITING MCMC\n");fflush(NULL);
    }
    catch (...) {
      matlabError("Error within MCMC. Check output logs");
    }
    /*****************/

    try {
      compute_stats(pop,0.95,(const int)(MAX_ITER-BURN_IN));
    }
    catch (...) {
      matlabError("Error computing summary statistics. Model output may still be valid");
    }

    rewind(fseed);
    for (int i=0;i<3;i++)
      fprintf(fseed,"%lu ",seed[i]);
    fclose(fseed);
    free(seed);

    fclose(flog);
    
    // printf("\n\nPlease see file Parameter_Estimates.log for parameter estimates and summaries\n");
    printToMatlab("\n\nPlease see file Parameter_Estimates.log for parameter estimates and summaries\n");
    
  };
};
