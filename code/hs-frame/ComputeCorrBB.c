#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"        /* the algorithm is connect to matlab */
#include "math.h"
#include "matrix.h"
#define ABS(x) ((x)>0? (x):(-(x)))
#define MAX(x, y) ((x)>(y)? (x):(y))
#define MIN(x, y) ((x)<(y)? (x):(y))

/*
[B'B]=ComputeCorrBB(Corr,template,single(halfFilterSize));
*/

float ROUND(float d){
    
    return floor(d+0.5);
}

int px(int x, int y, int bx, int by)    
{            
   return (x + (y-1)*bx - 1); 
 }

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
  int N,i,j,bytes_to_copy,g;   
 
  const mxArray *f;
  const mwSize* dims;
  mxClassID datatype;
  float **Corr;
  const mxArray* Template;      /* template */
  const float *selectedRow, *selectedCol, *selectedFilter;
  int numSelectedFeature;
  float* halfFilterSize;
  float * result;
  
  mwSize dimsOutput[2];
  void* start_of_pr;
  int i_ind, i_row, i_col, j_ind, j_row, j_col, hi, hj;
   
  
          
  /* input variable 0: correlation matrix  */
  dims = mxGetDimensions(prhs[0]);
  N=dims[0];
  Corr = (float**) mxCalloc(N*N, sizeof(float*));    /* C: correlation/inhibition between filters */
  for (i=0; i<N; i++)
  {  
     for (j=0; j<N; j++)
     {
         f = mxGetCell(prhs[0], j*N+i); 
         datatype = mxGetClassID(f);
         if (datatype != mxSINGLE_CLASS)
            mexErrMsgTxt("warning !! single precision required for correlation matrix.");
         Corr[j*N+i] = (const float*)mxGetPr(f);         
       
     }   
   }
  
  /*
     * input variable 1: template
     */
    Template = prhs[1];
    f = mxGetField( Template, 0, "selectedRow" );
	if ( f == NULL )
	{
	  mexErrMsgTxt( "no such field 1 in struct" );
	}
    selectedRow = (const float*)mxGetPr(f);
       
    
    f = mxGetField( Template, 0, "selectedCol" );
    if ( f == NULL )
	{
	  mexErrMsgTxt( "no such field 2 in struct" );
	}
    selectedCol = (const float*)mxGetPr(f);
    
    f = ( mxGetField( Template, 0, "selectedFilter" ) );
    if ( f == NULL )
	{
	  mexErrMsgTxt( "no such field 3 in struct" );
    }
    selectedFilter = (const float*)mxGetPr(f);
    
    numSelectedFeature= mxGetM(f)*mxGetN(f);
    
    /*
     * input variable 2: halfFilterSize
     */
    datatype = mxGetClassID(prhs[2]);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    halfFilterSize = (const float*)mxGetPr(prhs[2]); 
        
      
        
    /* =============================================
     * Compute the B'B
     * ============================================= 
     */
     
    /*  initialization */
    result = (float*) mxCalloc(numSelectedFeature*numSelectedFeature, sizeof(float));
    
    for(i=0; i<numSelectedFeature; i++ ){
         for(j=0; j<numSelectedFeature; j++ ){
             result[i+j*numSelectedFeature]=0.0;
         }
    }   
    
   
    /* compute B'B*/
    for(i=0; i<numSelectedFeature; i++ ){    
    
        i_ind=selectedFilter[i]-1;   /*  start from 0 */
        i_row=selectedRow[i];    /*  start from 1 */
        i_col=selectedCol[i];
        hi=halfFilterSize[i_ind];
    
        result[i+i*numSelectedFeature]=1;
        for(j=i+1; j<numSelectedFeature; j++ ){        
        
            j_ind=selectedFilter[j]-1;
            j_row=selectedRow[j];
            j_col=selectedCol[j];
            hj=halfFilterSize[j_ind];       
     
           if( ABS(j_row-i_row)<=(hi+hj) && ABS(j_col-i_col)<=(hi+hj)  ){
               
                g = px(j_row-i_row+hi+hj+1, j_col-i_col+hi+hj+1, 2*(hi+hj)+1, 2*(hi+hj)+1);
                result[i+j*numSelectedFeature]= Corr[i_ind+j_ind*N][g]; 
                result[j+i*numSelectedFeature]= result[i+j*numSelectedFeature];   
           }   
   
        }      
    
    }
    
            
     /* =============================================
     * Handle output variables.
     * ============================================= 
     */
    
    dimsOutput[0] = numSelectedFeature; dimsOutput[1] = numSelectedFeature;
    plhs[0]  = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL);
    start_of_pr = (float*)mxGetData(plhs[0]);   
    bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(plhs[0]);
    memcpy( start_of_pr, result, bytes_to_copy );         
    
    
}





