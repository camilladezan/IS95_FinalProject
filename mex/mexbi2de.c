/*
  Replacement for bi2de, expects a binary matrix (one row for each
  decimal number, one column for each bit, the first column being
  the LSB). Returns a column vector containing the decimal representation of the numbers.
  
  $Log: mexbi2de.c,v $

*/

#ifndef __MEXHEADER__
#define __MEXHEADER__

#include <mex.h>
#include <math.h>
#include <string.h>

#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define NIL 0L
#define nil 0L

typedef mxArray * mxArrayPtr;
typedef const mxArray * mxArrayConstPtr;
typedef double * doublePtr;

#endif

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A;
    int rowsA,colsA;
    unsigned int a;
    int dec;
    double *R;
    int r, c;
    
    // Get the dimensions of A and the data pointer to A
    A = mxGetPr(prhs[0]);
    rowsA =  mxGetM(prhs[0]);
    colsA =  mxGetN(prhs[0]);

    // Create the result array and get a pointer to it
    plhs[0] = mxCreateDoubleMatrix(rowsA, 1, mxREAL);
    R = mxGetPr(plhs[0]);
    
    // Fill the result vector (bi2dec)
    for(r=0; r<rowsA; r++) {
       dec=0;
       for(c=colsA-1; c>=0; c--) {
          dec = dec << 1;
          dec += (int)(A[r+rowsA*c]);
       }
       R[r] = (double)dec;
    }
}
