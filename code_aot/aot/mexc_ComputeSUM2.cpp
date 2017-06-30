/* mex-C: 
 * compute SUM2 maps for a *single* image
 *
 * Usage:
 *    S2Map = mexc_CcomputeSUM2( nGaborOri, M1Map, S2Template, subsampleS2 );
 *
 * S2Template is a matlab struct with the following fields:
 *	selectedRow: single array
 *	selectedCol: single array
 *	selectedOri: single array
 *	selectedScale: single array
 *	selectedLambda: single array
 *	selectedLogZ: single array
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

/* variable declaration */
const float **M1Map;            /* MAX 1 maps */
float** S2Map;                  /* SUM2 maps */
const mxArray* S2Template;      /* SUM2 templates */
int nT;                         /* number of S2 templates */
int* nElement;                  /* number of elements for each template */
int subsampleS2;                /* step length when scanning the SUM2 template over MAX1 map */
int heightS2Map, widthS2Map;    /* (down-sampled) size of SUM2 maps */
int heightM1Map, widthM1Map;    /* size of MAX1 maps */
int nGaborScale;
int nGaborOri;

void compute()
{
    int iT, iF;
    const mxArray* pA, *pA2;
    /* Note: we assume that iOri and iS starts from 0. */
    /* Note: we also assume that iRow and iCol is relative to (0,0). */
    int iRowS2, iColS2;
    int iOriM1, iRowM1, iColM1, iScaleM1;
    const float *selectedRow, *selectedCol;
    const float *selectedOri, *selectedScale;
    const float *selectedLambda, *selectedLogZ;
    mxClassID datatype;

    heightS2Map = (int)floor((double)heightM1Map / (double)subsampleS2);
    widthS2Map = (int)floor((double)widthM1Map / (double)subsampleS2);
    S2Map = (float**)mxCalloc( nT, sizeof(*S2Map) );

    /*
     * compute number of elements for the templates
     */
    nElement = (int*)mxCalloc( nT, sizeof(*nElement) );
    for( iT = 0; iT < nT; ++iT )
    {
        pA = mxGetCell( S2Template, iT );
        pA = mxGetField( pA, 0, "selectedRow" );
        nElement[iT] = mxGetM( pA ) * mxGetN( pA );
    }
    
    /* About the visiting order in the FOR loop:
     *      The scan over M1 map positions should be inner-most.
     */
    for( iT = 0; iT < nT; ++iT )
    {
        /* allocate memory for the iT-th S2 map */ 
        S2Map[iT] = (float*)mxCalloc( heightS2Map*widthS2Map, sizeof(**S2Map) );
        
        for( iColS2 = 0; iColS2 < widthS2Map; ++iColS2 )
        {
            for( iRowS2 = 0; iRowS2 < heightS2Map; ++iRowS2 )
            {
                S2Map[iT][iRowS2 + iColS2*heightS2Map] = 0.0;
            }
        }
        
        
        pA = mxGetCell( S2Template, iT );
		pA2 = mxGetField( pA, 0, "selectedRow" );
        selectedRow = (const float*)mxGetPr(pA2);
		pA2 = mxGetField( pA, 0, "selectedCol" );
        selectedCol = (const float*)mxGetPr(pA2);
        pA2 = ( mxGetField( pA, 0, "selectedOri" ) );
        selectedOri = (const float*)mxGetPr(pA2);
        pA2 = ( mxGetField( pA, 0, "selectedScale" ) );
        selectedScale = (const float*)mxGetPr(pA2);
        pA2 = ( mxGetField( pA, 0, "selectedLambda" ) );
        selectedLambda = (const float*)mxGetPr(pA2);
        pA2 = ( mxGetField( pA, 0, "selectedLogZ" ) );
        selectedLogZ = (const float*)mxGetPr(pA2);
        
        
        
        for( iF = 0; iF < nElement[iT]; ++iF )
        {
            iOriM1 = (int)selectedOri[iF];
            iScaleM1 = (int)selectedScale[iF];
            for( iColS2 = 0; iColS2 < widthS2Map; ++iColS2 )
            {
                iColM1 = (int)floor(.5+(iColS2+.5) * subsampleS2 + selectedCol[iF]);
                for( iRowS2 = 0; iRowS2 < heightS2Map; ++iRowS2 )
                {
                    iRowM1 = (int)floor(.5+(iRowS2+.5) * subsampleS2 + selectedRow[iF]);
                    if( iRowM1 < 0 || iRowM1 >= heightM1Map || iColM1 < 0 || iColM1 >= widthM1Map )
                    {
                        S2Map[iT][ iRowS2 + iColS2 * heightS2Map ] += 
                            -selectedLogZ[iF] + selectedLambda[iF] * 0;
                        continue;
                    }
                   
                    S2Map[iT][ iRowS2 + iColS2 * heightS2Map ] += 
                        -selectedLogZ[iF] + selectedLambda[iF] * M1Map[iOriM1+iScaleM1*nGaborOri][iRowM1+iColM1*heightM1Map];
                    
                    if( M1Map[iOriM1+iScaleM1*nGaborOri][iRowM1+iColM1*heightM1Map] > 1e4 )
                    {
                        mexPrintf("s2value=%.2f, iOri=%d,iScale=%d,iRow=%d,iCol=%d,lambda=%.2f, logZ=%.2f, m1value=%.2f\n", S2Map[iT][ iRowS2 + iColS2 * heightS2Map ], iOriM1,iScaleM1,iRowM1,iColM1,selectedLambda[iF],selectedLogZ[iF], M1Map[iOriM1+iScaleM1*nGaborOri][iRowM1+iColM1*heightM1Map]);
                        mexErrMsgTxt("for debugging: stop now");
                    }
                }
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
    const mxArray *pAM1Map;
	mxArray *pA;
    mwSize ndim;
    const mwSize* dims;
    mwSize dimsOutput[2];
    void* start_of_pr;
    mxClassID datatype;
    int iT;

	/*
	 * input variable 0: nGaborOri
	 */
	nGaborOri = (int)mxGetScalar(prhs[0]);
 
    /*
	 * input variable 1: M1 maps
	 */
    pAM1Map = prhs[1];
    dims = mxGetDimensions(pAM1Map);
    nGaborFilter = dims[0] * dims[1];
    nGaborScale = nGaborFilter / nGaborOri;
 
    M1Map = (const float**)mxCalloc( nGaborFilter, sizeof(*M1Map) );   /* MAX1 maps */
    for (i=0; i<nGaborFilter; ++i)
    {
        f = mxGetCell(pAM1Map, i);
        datatype = mxGetClassID(f);
        if (datatype != mxSINGLE_CLASS)
            mexErrMsgTxt("warning !! single precision required.");
        M1Map[i] = (const float*)mxGetPr(f);    /* get the pointer to cell content */
        heightM1Map = mxGetM(f);    /* overwriting is ok, since it is constant */
        widthM1Map = mxGetN(f);
    }

    /*
     * input variable 2: S2 templates
     */
    S2Template = prhs[2];
    nT = mxGetM(S2Template) * mxGetN(S2Template);
    
    /*
     * input variable 3: subsampleS2
     */
    subsampleS2 = (int)mxGetScalar(prhs[3]);
    

    compute();
    
    
    /* =============================================
     * Handle output variables.
     * ============================================= 
     */
    /*
     * output variable 0: S2 maps
     */
    dimsOutput[0] = nT; dimsOutput[1] = 1;
	plhs[0] = mxCreateCellArray( 2, dimsOutput );
    dimsOutput[0] = heightS2Map; dimsOutput[1] = widthS2Map;
    for( iT = 0; iT < nT; ++iT )
    {
        pA = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
        /* populate the real part of the created array */
        start_of_pr = (float*)mxGetData(pA);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(pA);
        memcpy( start_of_pr, S2Map[iT], bytes_to_copy );
        mxSetCell( plhs[0], iT, pA );
    }
}

