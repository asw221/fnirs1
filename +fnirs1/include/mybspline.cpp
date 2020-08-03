#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern FILE *flog;

double *create_bspline_basis(int N,double TIME,double freq,int int_knots,int *dim,double *knots,int type)
{
    int i,j;
    int datalength;
    int sdegree = 4,nknots;
    double *data,*Y;
    void bsplinebasis(int,int,int,double,double *,int,double *);
    void proj(double *out,double *v,double *u,int length);
   
    datalength = N-2;
    nknots = int_knots + 6;
    dim[0] = N;
    dim[1] = nknots - sdegree;
    
    fprintf(flog,"Time = %lf freq = %lf\n",TIME,freq);
    fprintf(flog,"int_knots = %d nknots = %d\n",int_knots,nknots);
    fprintf(flog,"dim = %d %d\n",dim[0],dim[1]);

    switch (type) {
        case 0: default: // Not a knot condition
            knots[0] = -0.00000111;
            knots[1] = -0.0000011;;
            knots[2] = -0.000001;
            j=3;
            for (i=0;i<=int_knots-1;i++) {
                knots[j] = (double)i*(N-1)/(int_knots-1);
                j++;
            }
            knots[j] = N-1 + 0.000001;
            knots[j+1] = N-1 + 0.0000011;
            knots[j+2] = N-1 + 0.00000111;
      
            break;
        case 1:  // end pts = 0 condition (not correct yet)
            knots[0] = -0.00000111;
            knots[1] = -0.0000011;;
            knots[2] = -0.000001;
            knots[3] = 0;
            knots[4] = 0.5*(N)/(int_knots-3);
            j = 5;
            for (i=1;i<int_knots-3;i++) {
                knots[j] = (double)i*(N)/(int_knots-3);
                j++;
            }
            knots[j] = knots[j-1] + 0.5*(N)/(int_knots-3);
            knots[j+1] = datalength;
            knots[j+2] = datalength + 0.000001;
            knots[j+3] = datalength + 0.0000011;
            knots[j+4] = datalength + 0.00000111;
            dim[1] -= 2;
            break;
        case 2: //  end pts + first derivative = 0 condition (not correct)
            knots[0] = -0.00000111;
            knots[1] = -0.0000011;;
            knots[2] = -0.000001;
            knots[3] = 0;
            knots[4] = 0.5*(N)/(int_knots-3);
            j = 5;
            for (i=1;i<int_knots-3;i++) {
                knots[j] = (double)i*(N)/(int_knots-3);
                j++;
            }
            knots[j] = knots[j-1] + 0.5*(N)/(int_knots-3);
            knots[j+1] = datalength;
            knots[j+2] = datalength + 0.000001;
            knots[j+3] = datalength + 0.0000011;
            knots[j+4] = datalength + 0.00000111;
            dim[1] -= 4;
            break;
        case 3: // end pts + first and second derivative = 0 condition
            knots[0] = -0.00000111;
            knots[1] = -0.0000011;;
            knots[2] = -0.000001;
            knots[3] = 0;
            knots[4] = 0.5*(N-1)/(int_knots-3);
            j = 5;
            for (i=1;i<int_knots-3;i++) {
                knots[j] = (double)i*(N-1)/(int_knots-3);
                j++;
            }
            knots[j] = knots[j-1] + 0.5*(N-1)/(int_knots-3);
            knots[j+1] = N-1;
            knots[j+2] = N-1 + 0.000001;
            knots[j+3] = N-1 + 0.0000011;
            knots[j+4] = N-1 + 0.00000111;
            dim[1] -= 6;
            break;
    }
    
    data = (double *)calloc(N,sizeof(double));
    j = 0;
    for (i=0;i<N;i++) {
        data[j] = i;
        j++;
    }

    Y = (double *)calloc(N*(nknots-sdegree),sizeof(double));
//    for (i=0;i<N;i++)
//        Y[i] = (double *)calloc(nknots-sdegree,sizeof(double));
    
    fprintf(flog,"%d %d %d\n",nknots-sdegree,dim[1],int_knots);
    for (i=0;i<N;i++) {
        bsplinebasis(0,0,sdegree,data[i],knots,nknots,&Y[i*(nknots-sdegree)]);
    }

    free(data);
 
    double *Y2;
    Y2 = (double *)calloc(N*dim[1],sizeof(double));
//    for (i=0;i<N;i++)
//        Y2[i] = (double *)calloc(dim[1],sizeof(double));
    switch (type) {
        case 0: default:
            for (i=0;i<N;i++)
                for (j=0;j<dim[1];j++)
                    Y2[i*dim[1]+j] = Y[i*(nknots-sdegree)+j];
 //           for (i=0;i<N;i++)
 //               free(Y[i]);
            free(Y);
            break;
        case 1:
            for (i=0;i<N;i++)
                for (j=0;j<dim[1];j++)
                    Y2[i*dim[1]+j] = Y[i*(nknots-sdegree)+j+1];
 //           for (i=0;i<N;i++)
 //               free(Y[i]);
            free(Y);
            break;
        case 2:
            for (i=0;i<N;i++)
                for (j=0;j<dim[1];j++)
                    Y2[i*dim[1]+j] = Y[i*(nknots-sdegree)+j+2];
//            for (i=0;i<N;i++)
//                free(Y[i]);
            free(Y);
            break;
        case 3:
            for (i=0;i<N;i++)
                for (j=0;j<dim[1];j++)
                    Y2[i*dim[1]+j] = Y[i*(nknots-sdegree)+j+3];
//           for (i=0;i<N;i++)
//                free(Y[i]);
            free(Y);
           break;
    }
/*    FILE *ff;
    ff = fopen("spline.log","w");
    for (i=0;i<N;i++) {
        for (j=0;j<dim[1];j++)
            fprintf(ff,"%g ",Y2[i*dim[1]+j]);
        fprintf(ff,"\n");
     }
     fclose(ff);*/
    return Y2;
}

void bsplinebasis(int j,int k,int SplineOrder,double x,double *knots,int Nknots,double *b)
{
	int i,start;
	if (k==0) {
		for (i=0;i<Nknots-SplineOrder;i++) 
			b[i] = 0;
		j = 0;
		while (x < knots[j] || knots[j+1] < x)
			j++;
		b[j] = 1;
	}
	else {
		start = (j-k > 0) ? j-k :0;
		for (i=start;i<j+1;i++) { 
			b[i] = (x-knots[i])*b[i]/(knots[i+k]-knots[i]);
			if (i < Nknots-SplineOrder-1)
				b[i] += (knots[i+k+1] - x)*b[i+1]/(knots[i+k+1]-knots[i+1]);
        }
    }
	
	if (k+1 < SplineOrder) 
		bsplinebasis(j,k+1,SplineOrder,x,knots,Nknots,b);
	
}

