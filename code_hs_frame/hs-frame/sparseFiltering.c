/* mex-C: 
 * extract SUM1 feature from original image by sparse convolution
 *
 * Usage:
 *
 *    response = sparseFiltering(currImage, filters, index_Filter, X, Y);
 *     
 *
 *  template is a matlab struct with the following fields:
 *	selectedRow: single array
 *	selectedCol: single array
 *	selectedFilter: single array

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
      
    int i,bytes_to_copy;
    const mxArray *f;
    const mxArray *pFilter;
	const mwSize* dims;
    mxClassID datatype;
    mwSize dimsOutput[2];
    int iF;
    mxArray *pA;
    int iRow, iCol;
    void* start_of_pr;
      
    
    const float **FilterMap;            /* Filter maps */
    const float *image;                 /* image */
    int* heightFilterMap,* widthFilterMap;
    int heightImage, widthImage, nFilter;
    float selectedRow, selectedCol, selectedFilter;
      

    float result;
    int index_filter,h,w,row,col;
    int imageStartingRow,imageStartingCol;
    int iLocation;
    float *output_pointer;
  
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
	 * input variable 1: filter maps
	 */
    pFilter = prhs[1];
    dims = mxGetDimensions(pFilter);
    nFilter = dims[0] * dims[1];
     
 
    FilterMap = (const float**)mxCalloc( nFilter, sizeof(*FilterMap) );   /* Filter maps */
    heightFilterMap = (int*)mxCalloc( nFilter, sizeof(int) );
    widthFilterMap = (int*)mxCalloc( nFilter, sizeof(int) );
     
   
    
    for (i=0; i<nFilter; ++i)
    {
        f = mxGetCell(pFilter, i);
        datatype = mxGetClassID(f);
        if (datatype != mxSINGLE_CLASS)
            mexErrMsgTxt("warning !! single precision required.");
        FilterMap[i] = (const float*)mxGetPr(f);    /* get the pointer to cell content */
        heightFilterMap[i] = mxGetM(f);    
        widthFilterMap[i] = mxGetN(f);
    }
    
    
    /*
     * input variable 2: S2 template
     
    Template = prhs[2];
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
    */ 
    
      /*
     * input variable 2: selected filter
     */
    selectedFilter  = mxGetScalar(prhs[2]); 
     /*
     * input variable 3: selected row
     */
    selectedRow = mxGetScalar(prhs[3]); 
     /*
     * input variable 4: selected col
     */
    selectedCol = mxGetScalar(prhs[4]); 
    
  
     /* =============================================
     * filtering
     * ============================================= */
   
        
        index_filter=(int)selectedFilter-1;
        h=heightFilterMap[index_filter];
        w=widthFilterMap[index_filter];
       
        row=(int)selectedRow-1;
        col=(int)selectedCol-1;
        
        imageStartingRow=row-ROUND(h/2.0-1);
        imageStartingCol=col-ROUND(w/2.0-1);
        
               
        result=0.0;
        iLocation=0;
        
        for (iCol=imageStartingCol; iCol<=imageStartingCol+w-1; iCol++)
        {
                    
          for(iRow=imageStartingRow; iRow<=imageStartingRow+h-1; iRow++)
          {
              
               if( (iRow>=0) && (iRow<heightImage) && (iCol>=0) && (iCol<widthImage)){
                   
                   result = result + image[iRow+iCol*heightImage] * FilterMap[index_filter][iLocation];  
               }
               
               iLocation++;
          }
          
        }        
                       
   
    
    /* =============================================
     * Handle output variables.
     * ============================================= 
     */
     /*
     * output variable 0: response
     */
     
    dimsOutput[0] =  1; dimsOutput[1] = 1;
    plhs[0] = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL);
    output_pointer=mxGetPr(plhs[0]);
    *output_pointer=result;
   
   
    
    
}
