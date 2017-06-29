/* mex-C: 
 * compute SUM3 maps for a single image
 *
 * Usage:
 *    S3Map = mexc_ComputeSUM3_logZ( M2Map, S3Template, subsampleS3, nS2Trans, logZ , displacementRatio);
 *
 *eg. tmpS3 = mexc_ComputeSUM3_logZ( tmpM2(:), S3T(r), 1, numPartRotate, clusters(c).logZ);
 *
 *      For S3 template, we use a Bayes tree model on top of S2 templates.
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

/* variable declaration */
const float **M2Map;            /* MAX 1 maps */
float RatioDisplacementSUM3;
float** S3Map;                  /* SUM3 maps */
const mxArray* S3Template;      /* SUM3 templates */
int nS3Template;                /* number of S3 templates */
int* nElement;                  /* number of elements for each template */
int subsampleS3;                /* step length when scanning the SUM2 template over MAX1 map */
int heightS3Map, widthS3Map;    /* (down-sampled) size of SUM2 maps */
int heightM2Map, widthM2Map;    /* size of MAX1 maps */
int nTS2;
int nS2Trans;
float logZ;

void compute()
{
    int iT, iF;
    const mxArray *pA, *f2;
    /* Note: we assume that iRow and iCol is relative to (0,0). */
    int iRowS3, iColS3; /* iRowS3, iColS3 belongs to a subsampled lattice */
    int iRowM2, iColM2, iTemplateM2;
    const float *selectedRow, *selectedCol;
    const float *selectedS2Ind, *selectedS2Transform;
    const float *selectedLambda;
    int rowOffset, colOffset;
    float r;
    int jTS2, jTS2Trans;
    mxClassID datatype; 
    int displacementWidth, displacementHeight;
    
    heightS3Map = (int)floor((double)heightM2Map / subsampleS3);
    rowOffset = 0;
    widthS3Map = (int)floor((double)widthM2Map / subsampleS3);
    colOffset = 0;
    S3Map = (float**)mxCalloc( nS3Template, sizeof(*S3Map) );
    
    /*
     * compute the lowest value in each M2 map: to deal with out-of-boundary case
     */
    float *minValPerM2Map = (float*)mxCalloc( nTS2, sizeof(*minValPerM2Map) );
    for( iT = 0; iT < nTS2; ++iT )
    {
    	float minval = M2Map[iT][0];
    	for( int j = 1; j < heightM2Map * widthM2Map; ++j )
    	{
    		if( M2Map[iT][j] < minval )
    		{
    			minval = M2Map[iT][j];
    		}
    	}
    	minValPerM2Map[iT] = minval;
    }
    
    /* 
     * compute number of elements for the S3 templates
     */
    nElement = (int*)mxCalloc( nS3Template, sizeof(*nElement) );
    for( iT = 0; iT < nS3Template; ++iT )
    {
        pA = mxGetCell( S3Template, iT );
        pA = mxGetField( pA, 0, "selectedRow" );
        nElement[iT] = mxGetM( pA ) * mxGetN( pA );
    }
    
    
    
    displacementWidth = RatioDisplacementSUM3 * widthS3Map;
    displacementHeight = RatioDisplacementSUM3 * heightS3Map;
    
    
    
    
    /* About the visiting order in the FOR loop:
     *      The scan over M1 map positions should be inner-most.
     */
    for( iT = 0; iT < nS3Template; ++iT )
    {
        /* allocate memory for the iT-th S2 map */ 
        S3Map[iT] = (float*)mxCalloc( heightS3Map*widthS3Map, sizeof(**S3Map) );
        
        for( iColS3 = 0; iColS3 < widthS3Map; ++iColS3 )
        {
            for( iRowS3 = 0; iRowS3 < heightS3Map; ++iRowS3 )
            {
                S3Map[iT][iRowS3 + iColS3*heightS3Map] = 0.0;
            }
        }
        
        /* load the selected features from the template */
        
        pA = mxGetCell(S3Template, iT);  
        f2 = mxGetField( pA, 0, "selectedRow" );
	    if ( f2 == NULL )
	    {
	      mexErrMsgTxt( "no such field 1 in struct" );
	    }
        datatype = mxGetClassID(f2);
        if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
      
        selectedRow = (const float*)mxGetPr(f2);
        
        
        f2 = mxGetField( pA, 0, "selectedCol" );
	    if ( f2 == NULL )
	    {
	      mexErrMsgTxt( "no such field 2 in struct" );
	    }
        datatype = mxGetClassID(f2);
        if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
      
        selectedCol = (const float*)mxGetPr(f2);
        
        
        f2 = mxGetField( pA, 0, "selectedInd" );
	    if ( f2 == NULL )
	    {
	      mexErrMsgTxt( "no such field 3 in struct" );
	    }
        datatype = mxGetClassID(f2);
        if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
      
        selectedS2Ind = (const float*)mxGetPr(f2);
        
        
        f2 = mxGetField( pA, 0, "selectedTransform" );
	    if ( f2 == NULL )
	    {
	      mexErrMsgTxt( "no such field 4 in struct" );
	    }
        datatype = mxGetClassID(f2);
        if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
      
        selectedS2Transform = (const float*)mxGetPr(f2);
        
        
        f2 = mxGetField( pA, 0, "selectedLambdas" );
	    if ( f2 == NULL )
	    {
	      mexErrMsgTxt( "no such field 5 in struct" );
	    }
        datatype = mxGetClassID(f2);
        if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
      
        selectedLambda = (const float*)mxGetPr(f2);
    
          

        
        
        for( iColS3 = 0 + displacementWidth; iColS3 < widthS3Map - displacementWidth; ++iColS3 )
        {
                
                for( iRowS3 = 0 + displacementHeight; iRowS3 < heightS3Map - displacementHeight; ++iRowS3 )
                {
                    
                    
                    for( iF = 0; iF < nElement[iT]; ++iF )
                    {
                    
                       jTS2 = selectedS2Ind[iF];
                       jTS2Trans = selectedS2Transform[iF];
                       
                       
                       iColM2 = (int)(colOffset + iColS3 * subsampleS3 + selectedCol[iF]);
                       iRowM2 = (int)(rowOffset + iRowS3 * subsampleS3 + selectedRow[iF]);
                       
                       if( iRowM2 < 0 || iRowM2 >= heightM2Map || iColM2 < 0 || iColM2 >= widthM2Map )
                       {
                        S3Map[iT][ iRowS3 + iColS3 * heightS3Map ] += 
                            minValPerM2Map[jTS2Trans+nS2Trans*jTS2] * selectedLambda[iF];
                        continue;
                       }
                    
                       S3Map[iT][ iRowS3 + iColS3 * heightS3Map ] += 
                           selectedLambda[iF] * M2Map[jTS2Trans+nS2Trans*jTS2][iRowM2+iColM2*heightM2Map];
                    
                       /* if( M2Map[jTS2Trans+nS2Trans*jTS2][iRowM2+iColM2*heightM2Map] > 1e4 )
                       {
                        
                           mexErrMsgTxt("for debugging: stop now");
                        *
                        }*/
                    
                    }
                    
                  S3Map[iT][ iRowS3 + iColS3 * heightS3Map ] = S3Map[iT][ iRowS3 + iColS3 * heightS3Map ] - logZ; 
                  
                }
        }
 
                
    }
}

/* mex function is used to pass on the pointers and scalars from matlab, 
   so that heavy computation can be done by C, which puts the results into 
   some of the pointers. After that, matlab can then use these results. 
   
   So matlab is very much like a managing platform for organizing the 
   experiments, and mex C is like a work enginee for fast computation. */

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
    int ind, i, x, y, o, dataDim, j, bytes_to_copy, nGaborFilter;
    const mxArray *f;
    const mxArray *pAM2Map;
    mwSize ndim;
    const mwSize* dims;
    mwSize dimsOutput[2];
    void* start_of_pr;
    mxClassID datatype;
    int iT;
    
    /*
	 * input variable 0: M2 maps
	 */
    pAM2Map = prhs[0];
    dims = mxGetDimensions(pAM2Map);
    nTS2 = dims[0] * dims[1];
 
    M2Map = (const float**)mxCalloc( nTS2, sizeof(*M2Map) );   /* MAX1 maps */
    for (i=0; i<nTS2; ++i)
    {
        f = mxGetCell(pAM2Map, i);
        datatype = mxGetClassID(f);
        if (datatype != mxSINGLE_CLASS)
            mexErrMsgTxt("warning !! single precision required.");
        M2Map[i] = (const float*)mxGetPr(f);    /* get the pointer to cell content */
        heightM2Map = mxGetM(f);    /* overwriting is ok, since it is constant */
        widthM2Map = mxGetN(f);
    }

    /*
     * input variable 1: S3 templates
     */
    S3Template = prhs[1];
    nS3Template= mxGetM(S3Template) * mxGetN(S3Template);
    
    /*
     * input variable 2: subsampleS3
     */
    subsampleS3 = (int)mxGetScalar(prhs[2]);
    
    /*
     * input variable 3: nS2Trans (number of transformations for each S2 template)
     */
    nS2Trans = (int)mxGetScalar(prhs[3]);
    
    logZ = (float)mxGetScalar(prhs[4]);
    
   
	RatioDisplacementSUM3 = (float)mxGetScalar(prhs[5]);
    
    compute();
    
    /* =============================================
     * Handle output variables.
     * ============================================= 
     */
    
    /*
     * output variable 0: S3 maps
     */
    dimsOutput[0] = nS3Template; dimsOutput[1] = 1;
	plhs[0] = mxCreateCellArray( 2, dimsOutput );
    dimsOutput[0] = heightS3Map; dimsOutput[1] = widthS3Map;
    for( iT = 0; iT < nS3Template; ++iT )
    {
        mxArray* f = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
        /* populate the real part of the created array */
        start_of_pr = (float*)mxGetData(f);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(f);
        memcpy( start_of_pr, S3Map[iT], bytes_to_copy );
        mxSetCell( plhs[0], iT, f );
    }

}

