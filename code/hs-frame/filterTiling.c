/* mex-C: 
 * tile the filter matrix into a big matrix with size imageSizeX and imageSizeY at center (centerX, centerY)
 *
 * Usage:
 *
 *    result = filterTiling(filter, imageSizeX, imageSizeY, centerX, centerY);
 *     
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
    float* result;                  /* output result */
    const float *filter; 
    int imageSizeX, imageSizeY,iCol,iRow;
    int centerX, centerY;
    int heightFilter, widthFilter;
    int iLocation;
    mwSize dimsOutput[2];
    mxClassID datatype;
    void* start_of_pr;
    int bytes_to_copy;
    int startingRow,startingCol;
    /*
	 * input variable 0: filter
	 */
	datatype = mxGetClassID(prhs[0]);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    filter = (const float*)mxGetPr(prhs[0]); 
    heightFilter = mxGetM(prhs[0]);    
    widthFilter = mxGetN(prhs[0]); 
    
    
    imageSizeX = (int)mxGetScalar(prhs[1]);
    imageSizeY = (int)mxGetScalar(prhs[2]);
    
    centerX = (int)mxGetScalar(prhs[3]);
    centerY = (int)mxGetScalar(prhs[4]);
    
    result = (float*)mxCalloc( imageSizeX * imageSizeY, sizeof(float) );
    
    for( iCol = 0; iCol < imageSizeY; ++iCol )
    {
        for( iRow = 0; iRow < imageSizeX; ++iRow )
        {
             result[iRow + iCol*imageSizeX] = 0.0;
        }
    }
    
    
    startingRow=(centerX-1)-ROUND(heightFilter/2.0-1);
    startingCol=(centerY-1)-ROUND(widthFilter/2.0-1);
    
    iLocation=0;
    for (iCol=startingCol; iCol<=startingCol+widthFilter-1; iCol++)
    {         
        for(iRow=startingRow; iRow<=startingRow+heightFilter-1; iRow++)
        {
            
            if( (iRow>=0) && (iRow<imageSizeX) && (iCol>=0) && (iCol<imageSizeY)){
                   
                   result[iRow+iCol*imageSizeX] = filter[iLocation];  
               }
            
             iLocation++;
    
        }
    }
    
       
     dimsOutput[0] = imageSizeX; dimsOutput[1] = imageSizeY;
     plhs[0] = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
     start_of_pr = (float*)mxGetData(plhs[0]);
     bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(plhs[0]);
     memcpy( start_of_pr, result, bytes_to_copy );
    
    
}