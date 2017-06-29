/* mex-C: 
 * extract SUM1 feature from original image by sparse convolution
 *
 * Usage:
 *
 *    [rModel_selected]= sparseFeatureExtraction(currImage,filters, template,isAbsoluteValue);
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
    float** filteredMaps;                  /* filtered maps (output result)*/
    const float *image;                 /* image */
    int* heightFilterMap,* widthFilterMap;
    int heightImage, widthImage, nFilter;
    const mxArray* Template;      /* template */
    const float *selectedRow, *selectedCol, *selectedFilter;
    int numSelectedFeature;
   

    float result;
    int index_filter,h,w,row,col;
    int imageStartingRow,imageStartingCol;
    int iLocation;
    
    float* arrayFilteringResponse;
    const mxArray * isAbsoluteValue;

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
     */
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
        
    
    /*
     * input variable 3: absolute value or not
     */
    isAbsoluteValue = prhs[3];
    if ( ! mxIsLogicalScalar(isAbsoluteValue) )
    {
        mexErrMsgTxt( "not Logical Scalar!" );
    }
    
/*    
    filteredMaps = (float**)mxCalloc( nFilter, sizeof(*filteredMaps) );
    
    for( iF = 0; iF < nFilter; ++iF )
    {
        
        filteredMaps[iF] = (float*)mxCalloc( heightImage*widthImage, sizeof(**filteredMaps) );
        
        for( iCol = 0; iCol < widthImage; ++iCol )
        {
            for( iRow = 0; iRow < heightImage; ++iRow )
           {
                filteredMaps[iF][iRow + iCol*heightImage] = 0.0;
           }
        }
     }  

*/
    
   /* printf("%d", numSelectedFeature); */
     /* =============================================
     * sparse filtering
     * ============================================= */
   
   arrayFilteringResponse = (float*)mxCalloc( numSelectedFeature, sizeof(float) );
   
   
   if( mxIsLogicalScalarTrue(isAbsoluteValue)){

       /*printf("sparse filtering (absolute value).\n");*/
       
       for (i=0; i<numSelectedFeature; i++)
       {
        
        index_filter=(int)selectedFilter[i]-1;
        h=heightFilterMap[index_filter];
        w=widthFilterMap[index_filter];
       
        row=(int)selectedRow[i]-1;
        col=(int)selectedCol[i]-1;
        
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
                       
        
        /*   filteredMaps[index_filter][row+col*heightImage]=ABS(result);  */
        
           arrayFilteringResponse[i]=ABS(result);
        }

   }else{
       
       /*printf("sparse filtering (convolution value).\n");*/
       
       for (i=0; i<numSelectedFeature; i++)
       {
        
        index_filter=(int)selectedFilter[i]-1;
        h=heightFilterMap[index_filter];
        w=widthFilterMap[index_filter];
       
        row=(int)selectedRow[i]-1;
        col=(int)selectedCol[i]-1;
        
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
                       
        
         
     /*     filteredMaps[index_filter][row+col*heightImage]=result;   */    /*here is the difference*/  
        
          arrayFilteringResponse[i]=result;    /* here is the difference*/
       }

       
   }
    
    /* =============================================
     * Handle output variables.
     * ============================================= 
     */
     /*
     * output variable 0: array of filtering response
     */
     
     dimsOutput[0] = 1; dimsOutput[1] = numSelectedFeature;
     plhs[0] = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
     start_of_pr = (float*)mxGetData(plhs[0]);
     bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(plhs[0]);
     memcpy( start_of_pr, arrayFilteringResponse, bytes_to_copy );
     
     
     /*
     * output variable 1: filtered maps (sparse maps)
     */ 
/*
     dimsOutput[0] = nFilter; dimsOutput[1] = 1;
	 plhs[1] = mxCreateCellArray( 2, dimsOutput );
     dimsOutput[0] = heightImage; dimsOutput[1] = widthImage;
     for( iF = 0; iF < nFilter; ++iF )
     {
        pA = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
        
        start_of_pr = (float*)mxGetData(pA);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(pA);
        memcpy( start_of_pr, filteredMaps[iF], bytes_to_copy );
        mxSetCell( plhs[1], iF, pA );
     }
    
  */   
    
}
