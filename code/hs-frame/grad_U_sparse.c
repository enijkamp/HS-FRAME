/* mex-C: 
 * dE(I)/d(I) = -sum_{lambdas * sign(<I,B>) B} + I
 *
 * Usage:
 *   
 *	  g_U = grad_U_sparse(currImage, filters, template)
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

int getSign(float x){
   
    if (x>0){
        
        return 1;
        
    }else if(x<0){
        
        return -1;
        
    }else{
        
        return 0;
    }
    
}


float ROUND(float d){
    
    return floor(d+0.5);
}


void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
    mxClassID datatype;
    const float *image;   
    int heightImage, widthImage, nFilter;
    const mxArray *pFilter;
	const mwSize* dims;
    const float **FilterMap;            /* Filter maps */
    int* heightFilterMap,* widthFilterMap;
    int i,bytes_to_copy;
    const mxArray *f;
    int numSelectedFeature;
    const mxArray* Template;      /* template */
    const float *selectedRow, *selectedCol, *selectedFilter, *selectedLambdas;
    int index_filter,h,w,row,col;
    int imageStartingRow,imageStartingCol;
    int iLocation;
    float sum;
    int iRow, iCol;
    float* response;
    float* result;                  /* output result */
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
     * input variable 2: template
     */
    Template = prhs[2];
    f = mxGetField( Template, 0, "selectedRow" );
	if ( f == NULL )
	{
	  mexErrMsgTxt( "no such field 1 in struct" );
	}
    datatype = mxGetClassID(f);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    
    selectedRow = (const float*)mxGetPr(f);
    
   
    
    f = mxGetField( Template, 0, "selectedCol" );
    if ( f == NULL )
	{
	  mexErrMsgTxt( "no such field 2 in struct" );
	}
    
    datatype = mxGetClassID(f);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    
    selectedCol = (const float*)mxGetPr(f);
    
    f = ( mxGetField( Template, 0, "selectedFilter" ) );
    if ( f == NULL )
	{
	  mexErrMsgTxt( "no such field 3 in struct" );
    }
    datatype = mxGetClassID(f);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    
    selectedFilter = (const float*)mxGetPr(f);
    
    f = ( mxGetField( Template, 0, "selectedLambdas" ) );
    if ( f == NULL )
	{
	  mexErrMsgTxt( "no such field 4 in struct" );
    }
    
    datatype = mxGetClassID(f);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    
    selectedLambdas = (const float*)mxGetPr(f);
    
    numSelectedFeature= mxGetM(f)*mxGetN(f);
    
    
    
    /* Step 1: compute the sparse filtering response, and save them into an array in the order by which features ared stored in template   */
    response = (float*)mxCalloc( numSelectedFeature, sizeof(float) );
    
    for (i=0; i<numSelectedFeature; i++)
    {
        
        index_filter=(int)selectedFilter[i]-1;
        h=heightFilterMap[index_filter];
        w=widthFilterMap[index_filter];
       
        row=(int)selectedRow[i]-1;
        col=(int)selectedCol[i]-1;
        
        imageStartingRow=row-ROUND(h/2.0-1);
        imageStartingCol=col-ROUND(w/2.0-1);
        
               
        sum=0.0;
        iLocation=0;
        
        for (iCol=imageStartingCol; iCol<=imageStartingCol+w-1; iCol++)
        {
                    
          for(iRow=imageStartingRow; iRow<=imageStartingRow+h-1; iRow++)
          {
              
               if( (iRow>=0) && (iRow<heightImage) && (iCol>=0) && (iCol<widthImage)){
                   
                   sum = sum + image[iRow+iCol*heightImage] * FilterMap[index_filter][iLocation];  
               }
               
               iLocation++;
          }
          
        }        
                       
        response[i]=sum;
     
    }
    
    
    /* Step 2: prepare the results by intializing it with zeros value*/
    
    result = (float*)mxCalloc( heightImage * widthImage, sizeof(float) );
    
    for( iCol = 0; iCol < widthImage; ++iCol )
    {
        for( iRow = 0; iRow < heightImage; ++iRow )
        {
             result[iRow + iCol*heightImage] = 0.0;
        }
    }
    
    /*  Step 3: compute: sum_{lambdas * sign(<I,B>) B} */
    for (i=0; i<numSelectedFeature; i++)
    { 
        
    
        index_filter=(int)selectedFilter[i]-1;
        h=heightFilterMap[index_filter];
        w=widthFilterMap[index_filter];
        
        row=(int)selectedRow[i]-1;
        col=(int)selectedCol[i]-1;
        
        imageStartingRow=row-ROUND(h/2.0-1);
        imageStartingCol=col-ROUND(w/2.0-1);
                
       
        
        iLocation=0;
        for (iCol=imageStartingCol; iCol<=imageStartingCol+w-1; iCol++)
        {         
           for(iRow=imageStartingRow; iRow<=imageStartingRow+h-1; iRow++)
           {
            
               if( (iRow>=0) && (iRow<heightImage) && (iCol>=0) && (iCol<widthImage)){
                   
                   result[iRow+iCol*heightImage] += selectedLambdas[i]* getSign(response[i])* FilterMap[index_filter][iLocation];  
               
               }           
            
               iLocation++;            
           }
   
        }
    
    
    }
    
    /*  Step 4: dE(I)/d(I) = -sum_{lambdas * sign(<I,B>) B} + I  */
    
    for (iCol=0; iCol<widthImage; iCol++)
    {         
        
        for(iRow=0; iRow<heightImage; iRow++)
        {
             
            result[iRow+iCol*heightImage] = -result[iRow+iCol*heightImage] + image[iRow+iCol*heightImage]; 
           
        }
          
    }
   
     dimsOutput[0] = heightImage; dimsOutput[1] = widthImage;
     plhs[0] = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
     start_of_pr = (float*)mxGetData(plhs[0]);
     bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(plhs[0]);
     memcpy( start_of_pr, result, bytes_to_copy );
    
    
}