/* mex-C: 
 * I = I + cB
 *
 * Usage:
 *    I = I_plus_cB(I, c, filters{template.selectedFilter(i)}, template.selectedRow(i), template.selectedCol(i));
 *	  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"        /* the algorithm is connect to matlab */
#include "math.h"
#include "matrix.h"
#define ABS(x) ((x)>0? (x):(-(x)))
#define MAX(x, y) ((x)>(y)? (x):(y))
#define MIN(x, y) ((x)<(y)? (x):(y))


float ROUND(float d){
    
    return floor(d+0.5);
}


void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
    mxClassID datatype;
    const float *image; 
    const float *basis;
    float c;
    int heightBasis, widthBasis;
    int heightImage, widthImage, nFilter;
    float* result;                  /* output result */
    int iRow, iCol;
    int imageStartingRow,imageStartingCol;
    int iLocation;
    int bytes_to_copy;    
    int row,col;  
    mwSize dimsOutput[2];      
    void* start_of_pr;       
     
    
    /*
	 * input variable 0: image
	 */
	datatype = mxGetClassID(prhs[0]);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    image = (const float*)mxGetPr(prhs[0]); 
    heightImage = mxGetM(prhs[0]);    
    widthImage = mxGetN(prhs[0]);
    
    /*
	 * input variable 1: coefficient (c)
	 */
    c = mxGetScalar(prhs[1]); 
     
    /*
	 * input variable 2: basis (B)
	 */
	datatype = mxGetClassID(prhs[2]);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    basis = (const float*)mxGetPr(prhs[2]); 
    heightBasis = mxGetM(prhs[2]);    
    widthBasis = mxGetN(prhs[2]);
    
    /*
	 * input variable 3 and 4: center of basis
	 */
    row = ROUND(mxGetScalar(prhs[3]))-1; 
    col = ROUND(mxGetScalar(prhs[4]))-1; 
  
       
    /* prepare the results by intializing it with image map I*/
    
    result = (float*)mxCalloc( heightImage * widthImage, sizeof(float) );
    
    for( iCol = 0; iCol < widthImage; ++iCol )
    {
        for( iRow = 0; iRow < heightImage; ++iRow )
        {
             result[iRow + iCol*heightImage] = image[iRow + iCol*heightImage];
        }
    }
    
    /*  compute: I+cB */
           
    imageStartingRow=row-ROUND(heightBasis/2.0-1);
    imageStartingCol=col-ROUND(widthBasis/2.0-1);
          
        
    iLocation=0;
    for (iCol=imageStartingCol; iCol<=imageStartingCol+widthBasis-1; iCol++)
    {         
           for(iRow=imageStartingRow; iRow<=imageStartingRow+heightBasis-1; iRow++)
           {
            
               if( (iRow>=0) && (iRow<heightImage) && (iCol>=0) && (iCol<widthImage)){
                   
                   result[iRow+iCol*heightImage] += c * basis[iLocation];  
               
               }           
            
               iLocation++;            
           }
   
    }
    
    
   
     dimsOutput[0] = heightImage; dimsOutput[1] = widthImage;
     plhs[0] = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
     start_of_pr = (float*)mxGetData(plhs[0]);
     bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(plhs[0]);
     memcpy( start_of_pr, result, bytes_to_copy );
    
    
}