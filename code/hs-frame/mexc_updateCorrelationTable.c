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
 *   [corrBB]=update_correlation(corrBB, Corr, single(halfFilterSize), filter_added, x_added, y_added, template);
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
  float *corrBB;
  const mxArray* Template;      /* template */
  const float *selectedRow, *selectedCol, *selectedFilter;
  int numSelectedFeature;
  float* halfFilterSize;
  float * result;
  
  mwSize dimsOutput[2];
  void* start_of_pr;
  int i_ind, i_row, i_col, add_ind, add_row, add_col, hi, h_add;
  int heightCorrBB, widthCorrBB; 
  
  
  /*
	 * input variable 0: corrBB
  */
  datatype = mxGetClassID(prhs[0]);
  if (datatype != mxSINGLE_CLASS)
     mexErrMsgTxt("warning !! single precision required.");
  corrBB = (const float*)mxGetPr(prhs[0]); 
  heightCorrBB = mxGetM(prhs[0]);    
  widthCorrBB = mxGetN(prhs[0]);
    
  
  
          
  /* input variable 1: correlation matrix  */
  dims = mxGetDimensions(prhs[1]);
  N=dims[0];
  Corr = (float**) mxCalloc(N*N, sizeof(float*));    /* C: correlation/inhibition between filters */
  for (i=0; i<N; i++)
  {  
     for (j=0; j<N; j++)
     {
         f = mxGetCell(prhs[1], j*N+i); 
         datatype = mxGetClassID(f);
         if (datatype != mxSINGLE_CLASS)
            mexErrMsgTxt("warning !! single precision required for correlation matrix.");
         Corr[j*N+i] = (const float*)mxGetPr(f);         
       
     }   
   }
  
 
    /*
     * input variable 2: halfFilterSize
     */
    datatype = mxGetClassID(prhs[2]);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    halfFilterSize = (const float*)mxGetPr(prhs[2]); 
        
    
     /*
     * input variable 3: template
    */
    Template = prhs[3];
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
     
           
    /* =============================================
     * Compute the B'B
     * ============================================= 
     */
      
    
   add_ind=selectedFilter[numSelectedFeature-1]-1;
   add_row=selectedRow[numSelectedFeature-1];
   add_col=selectedCol[numSelectedFeature-1];
   h_add=halfFilterSize[add_ind];   
        
        
    /* compute B'B only for last column and last row*/
    for(i=0; i<numSelectedFeature-1; i++ ){    
    
        i_ind=selectedFilter[i]-1;   /*  start from 0 */
        i_row=selectedRow[i];    /*  start from 1 */
        i_col=selectedCol[i];
        hi=halfFilterSize[i_ind];
    
       /* result[i+i*numSelectedFeature]=1; */
       
     
        if( ABS(add_row-i_row)<=(hi+h_add) && ABS(add_col-i_col)<=(hi+h_add)  ){
               
            g = px(add_row-i_row+hi+h_add+1, add_col-i_col+hi+h_add+1, 2*(hi+h_add)+1, 2*(hi+h_add)+1);
            corrBB[i+(numSelectedFeature-1)*heightCorrBB]= Corr[i_ind+add_ind*N][g]; 
            corrBB[(numSelectedFeature-1)+i*heightCorrBB]= corrBB[i+(numSelectedFeature-1)*heightCorrBB];   
        }   
   
              
    
    }
   
    corrBB[(numSelectedFeature-1)+(numSelectedFeature-1)*heightCorrBB]=1; 
   
    
            
     /* =============================================
     * Handle output variables.
     * ============================================= 
     
    
    dimsOutput[0] = heightCorrBB; dimsOutput[1] = widthCorrBB;
    plhs[0]  = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL);
    start_of_pr = (float*)mxGetData(plhs[0]);   
    bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(plhs[0]);
    memcpy( start_of_pr, corrBB, bytes_to_copy );         
    */
    
}





