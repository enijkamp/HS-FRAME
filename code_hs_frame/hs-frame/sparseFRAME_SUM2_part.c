/* mex-C: 
 * compute SUM2 maps by all part template for one 'scale' image
 *
 * Usage:
 *   
 *    SUM2 = sparseFRAME_SUM2_part(single(allSizex), single(allSizey), single(numFilter), S2T{:}, MAX1mapFind(iRes,:));   
 *
 *  template is a matlab struct with the following fields:
 *	selectedRow: single array
 *	selectedCol: single array
 *	selectedFilter: single array
 *	selectedLambdas: single array
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
       
    
    int iT,iF,i,bytes_to_copy, numTemplate;    
    mxClassID datatype;    
    const mxArray* Template;  
    const mxArray *f, *f2;
    const float *selectedRow, *selectedCol, *selectedFilter, *selectedLambdas;
    int nElement,  numFilters;
    int heightM1Map, widthM1Map;    
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
	 * input variable 0: heightS2Map
	 */
    heightS2Map = (int)mxGetScalar(prhs[0]);
    heightM1Map = heightS2Map;
    /*  
    *  input variable 1: widthS2Map
    */
    widthS2Map = (int)mxGetScalar(prhs[1]);
    widthM1Map = widthS2Map;
    
     /*  
    *  input variable 2: numFilters
    */
    
    numFilters = (int)mxGetScalar(prhs[2]);
    /*  
    *  input variable 3: template
    */
    
    
    Template = prhs[3];
    numTemplate = mxGetM(Template)*mxGetN(Template);
    
    
    
    
 /*   for (i=0; i<nElement; i++)
    mexPrintf("%d: %f \n",i, selectedLambdas[i]);
   */ 
    
    
    /*  
    *  input variable 4: MAX1 map cell {1, iFilter} (one resolution)
    */
    pAM1Map = prhs[4];
    dims = mxGetDimensions(pAM1Map);
     
    M1Map = (const float**)mxCalloc( dims[0]*dims[1], sizeof(*M1Map) );   /* MAX1 maps */
   
        
    for (iF=0; iF<numFilters; ++iF){  
    
          f = mxGetCell(pAM1Map, iF);
          datatype = mxGetClassID(f);
          if (datatype != mxSINGLE_CLASS)
            mexErrMsgTxt("warning !! single precision required.");
          M1Map[iF] = (const float*)mxGetPr(f);    /* get the pointer to cell content */      
      
       
    }
    
      
   
    /* initialize SUM2 maps for all resolutions*/
    S2Map = (float**)mxCalloc( numTemplate, sizeof(*S2Map) );
        
    for (iT=0; iT<numTemplate; ++iT){    
          
       S2Map[iT] = (float*)mxCalloc( heightS2Map*widthS2Map, sizeof(**S2Map) );
        
       for( iColS2 = 0; iColS2 < widthS2Map; ++iColS2 )
       {
          for( iRowS2 = 0; iRowS2 < heightS2Map; ++iRowS2 )
          {
             S2Map[iT][iRowS2 + iColS2*heightS2Map] = 0.0;
          }
        }
    
    }
    
    
    
    /* compute the SUM2 map for all resolutions*/
    
    for (iT=0; iT<numTemplate; ++iT)
    {
      
      /* get each template*/  
      f = mxGetCell(Template, iT);  
      f2 = mxGetField( f, 0, "selectedRow" );
	  if ( f2 == NULL )
	  {
	    mexErrMsgTxt( "no such field 1 in struct" );
	  }
      datatype = mxGetClassID(f2);
      if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
      
      selectedRow = (const float*)mxGetPr(f2);
    
      f2 = mxGetField( f, 0, "selectedCol" );
      if ( f2 == NULL )
	  {
	    mexErrMsgTxt( "no such field 2 in struct" );
	  }
      datatype = mxGetClassID(f2);
      if (datatype != mxSINGLE_CLASS)
         mexErrMsgTxt("warning !! single precision required.");  
      
      selectedCol = (const float*)mxGetPr(f2);
    
      f2 = ( mxGetField( f, 0, "selectedFilter" ) );
      if ( f2 == NULL )
      {
	     mexErrMsgTxt( "no such field 3 in struct" );
      }
      datatype = mxGetClassID(f2);
      if (datatype != mxSINGLE_CLASS)
         mexErrMsgTxt("warning !! single precision required.");    
      
      selectedFilter = (const float*)mxGetPr(f2);       
    
      f2 = ( mxGetField( f, 0, "selectedLambdas" ) );
      if ( f2 == NULL )
	  {
	    mexErrMsgTxt( "no such field 4 in struct" );
      }
      datatype = mxGetClassID(f2);
      if (datatype != mxSINGLE_CLASS)
         mexErrMsgTxt("warning !! single precision required.");  
      
      selectedLambdas = (const float*)mxGetPr(f2);
    
      nElement= mxGetM(f2)*mxGetN(f2);  
        
      
      /* compute SUM2 map for each tempate */
      for( iF = 0; iF < nElement; ++iF )
      {
            
            iFilterM1 = (int)selectedFilter[iF]-1;
            
            for( iColS2 = 0; iColS2 < widthS2Map; ++iColS2 )
            {
                iColM1 = iColS2 + ((int)selectedCol[iF]-1);
                for( iRowS2 = 0; iRowS2 < heightS2Map; ++iRowS2 )
                {
                    iRowM1 = iRowS2  + ((int)selectedRow[iF]-1);
                   
                  /*  if( iRowM1 < 0 || iRowM1 >= heightM1Map[iR] || iColM1 < 0 || iColM1 >= widthM1Map[iR] )
                    {
                        S2Map[iR][ iRowS2 + iColS2 * heightS2Map ] += 
                            selectedLambdas[iF] * 0;
                        continue;
                    }
                    
                    
                    S2Map[iR][ iRowS2 + iColS2 * heightS2Map ] += 
                        selectedLambdas[iF] * M1Map[iR+iFilterM1*numResolution][iRowM1+iColM1*heightM1Map[iR]];
                    
                    */ 
                    
                    if( iRowM1 >= 0 && iRowM1 < heightM1Map && iColM1 >= 0 && iColM1 < widthM1Map )
                    {
                        
                        S2Map[iT][ iRowS2 + iColS2 * heightS2Map ] += 
                        selectedLambdas[iF] * M1Map[iFilterM1][iRowM1+iColM1*heightM1Map];
                        
                       
                    }                    
                           
                                   
                }
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
    dimsOutput[0] = numTemplate; dimsOutput[1] = 1;
	plhs[0] = mxCreateCellArray( 2, dimsOutput );
    dimsOutput[0] = heightS2Map; dimsOutput[1] = widthS2Map;
    
    for( iT = 0; iT < numTemplate; ++iT )
    {
        pA = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
        /* populate the real part of the created array */
        start_of_pr = (float*)mxGetData(pA);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(pA);
        memcpy( start_of_pr, S2Map[iT], bytes_to_copy );
        mxSetCell( plhs[0], iT, pA );
    }
}
