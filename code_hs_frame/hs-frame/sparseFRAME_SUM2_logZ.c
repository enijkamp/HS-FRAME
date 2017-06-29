/* mex-C: 
 * compute SUM2 maps for a *multiple scales* image
 *
 * Usage:
 *    [SUM2] = sparseFRAME_SUM2_logZ(numResolution, allSizex, allSizey, numFilter, template, MAX1mapFind,halfTempSizex,halfTempSizey,MaxHalfFilterSize)
 * 
 *  template is a matlab struct with the following fields:
 *	selectedRow: single array
 *	selectedCol: single array
 *	selectedFilter: single array
 *	selectedLambdas: single array
 *  logZ: single
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


void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
       
    
    int iR,iF,i,bytes_to_copy;
    int numResolution;
    mxClassID datatype;
    const float *allSizex, *allSizey;
    const mxArray* Template;  
    const mxArray *f;
    const float *selectedRow, *selectedCol, *selectedFilter, *selectedLambdas;
    float logZ;
    int nElement,  numFilters;
    int heightM1Map, widthM1Map;
    int halfTempSizex,halfTempSizey,MaxHalfFilterSize;
    float** S2Map;                  /* SUM2 maps */
    int heightS2Map, widthS2Map;
    int iColS2, iRowS2,iColM1,iRowM1;
    int iFilterM1;
    mxArray *pA;
    mwSize dimsOutput[2];
    void* start_of_pr;
    const mxArray *pAM1Map;
    const mwSize* dims;
    const float **M1Map;            /* MAX 1 maps */
    
	/*
	 * input variable 0: numResolution
	 */
	numResolution = (int)mxGetScalar(prhs[0]);
 
    /*
	 * input variable 1: allSizex
	 */
    datatype = mxGetClassID(prhs[1]);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    allSizex = (const float*)mxGetPr(prhs[1]); 
    
    /*  
    *  input variable 2: allSizey
    */
    datatype = mxGetClassID(prhs[2]);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    allSizey = (const float*)mxGetPr(prhs[2]); 
    
    
     /*  
    *  input variable 3: numFilters
    */
    
    numFilters = (int)mxGetScalar(prhs[3]);
    /*  
    *  input variable 4: template
    */
    Template = prhs[4];
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
    
    nElement= mxGetM(f)*mxGetN(f);
    
    f = ( mxGetField( Template, 0, "logZ" ) );
    if ( f == NULL )
	{
	  mexErrMsgTxt( "no such field 5 in struct" );
    }
    datatype = mxGetClassID(f);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    
    logZ = (float)mxGetScalar(f);
    
    
 /*   for (i=0; i<nElement; i++)
    mexPrintf("%d: %f \n",i, selectedLambdas[i]);
   */ 
    
    
    /*  
    *  input variable 5: MAX1 map cell {iResulotion, iFilter}
    */
    pAM1Map = prhs[5];
    dims = mxGetDimensions(pAM1Map);
     
    M1Map = (const float**)mxCalloc( (dims[0] * dims[1]), sizeof(*M1Map) );   /* MAX1 maps */
   
        
    for (iF=0; iF<numFilters; ++iF){    
  
      for (iR=0; iR<numResolution; ++iR)            
      {
          i=iR+iF*numResolution;
          f = mxGetCell(pAM1Map, i);
          datatype = mxGetClassID(f);
          if (datatype != mxSINGLE_CLASS)
            mexErrMsgTxt("warning !! single precision required.");
          M1Map[i] = (const float*)mxGetPr(f);    /* get the pointer to cell content */
       
       }
       
    }
    
      
   
    /* initialize SUM2 maps for all resolutions*/
    S2Map = (float**)mxCalloc( numResolution, sizeof(*S2Map) );
        
    for (iR=0; iR<numResolution; ++iR){    
          
       heightS2Map = allSizex[iR];
       widthS2Map = allSizey[iR];
       
       S2Map[iR] = (float*)mxCalloc( heightS2Map*widthS2Map, sizeof(**S2Map) );
        
       for( iColS2 = 0; iColS2 < widthS2Map; ++iColS2 )
       {
          for( iRowS2 = 0; iRowS2 < heightS2Map; ++iRowS2 )
          {
             S2Map[iR][iRowS2 + iColS2*heightS2Map] = 0.0;
          }
        }
    
    }
    
    /* compute the SUM2 map for all resolutions*/
    
    for (iR = 0; iR < numResolution; ++iR)
    {
        
      heightS2Map = allSizex[iR];
      widthS2Map = allSizey[iR];
      
      heightM1Map = allSizex[iR];
      widthM1Map = allSizey[iR];      
     
            
            for( iColS2 = 0; iColS2 < widthS2Map; ++iColS2 )
            {
                
                for( iRowS2 = 0; iRowS2 < heightS2Map; ++iRowS2 )
                {
                    
                   
                  for( iF = 0; iF < nElement; ++iF )
                  {
                     iColM1 = iColS2 + ((int)selectedCol[iF]-1);
                     iRowM1 = iRowS2  + ((int)selectedRow[iF]-1);
                     iFilterM1 = (int)selectedFilter[iF]-1;
                    
                     if( iRowM1 >= 0 && iRowM1 < heightM1Map && iColM1 >= 0 && iColM1 < widthM1Map )
                     {
                        S2Map[iR][ iRowS2 + iColS2 * heightS2Map ] += 
                        selectedLambdas[iF] * M1Map[iR+iFilterM1*numResolution][iRowM1+iColM1*heightM1Map];                        
                       
                     }                        
                        
                   }
                    
                    S2Map[iR][ iRowS2 + iColS2 * heightS2Map ] -= logZ; 
                                   
                }
            }
       
      
    }
      
    
    /* =============================================
     * Handle output variables.
     * ============================================= 
     */
    /*
     * output variable 0: S2 maps
     */
    dimsOutput[0] = 1; dimsOutput[1] = numResolution;
	plhs[0] = mxCreateCellArray( 2, dimsOutput );
    
    for( iR = 0; iR < numResolution; ++iR )
    {
        dimsOutput[0] = allSizex[iR]; dimsOutput[1] = allSizey[iR];
        pA = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
        /* populate the real part of the created array */
        start_of_pr = (float*)mxGetData(pA);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(pA);
        memcpy( start_of_pr, S2Map[iR], bytes_to_copy );
        mxSetCell( plhs[0], iR, pA );
    }
}
