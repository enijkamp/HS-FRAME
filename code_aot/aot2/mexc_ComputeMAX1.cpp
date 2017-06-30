/* mex-C: 
 * compute MAX1 maps for a *single* image
 *
 * Usage:
 * 
 *    [M1Map M1Trace M1RowShift M1ColShift M1OriShifted] = 
 *      CcomputeMAX1( nGaborOri, S1Map, locationPerturb, 
 *      orientationPerturb, subsampleM1 );
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"        /* the algorithm is connect to matlab */
#include <cmath>
#include "matrix.h"
#define ABS(x) ((x)>0? (x):(-(x)))
#define MAX(x, y) ((x)>(y)? (x):(y))
#define MIN(x, y) ((x)<(y)? (x):(y))
# define PI 3.1415926
# define NEGMAX -1e10

/* variable declaration */
float **M1Map;                  /* MAX 1 maps */
int **M1Trace;                  /* we currently do not allow shift in scale, so no need to track original scale */
const float** S1Map;            /* SUM1 maps */
int height, width;              /* size of MAX1 maps */
int nGaborScale;                /* number of Gabor scales */
int nGaborOri;                  /* number of orientations for Gabor elements */
float *sizeGabor;               /* size of Gabor filters at different scales */ 
int numShift;                  /* number of local shifts, the same for all orientations and all scales */
int **rowShift, **colShift, **orientShift;  /* store the location/orientation shifts per scale and orientation */
int sizexSubsample, sizeySubsample;         /* (down sampled) size of M1 map */
int subsample;                              /* sub sampling step length */
int orientationPerturb;               /* allowed orientation perturbation */
int locationPerturb;                /* allowed location perturbation, relative to the Gabor filter size */

/* store all the shifts or perturbations in local maximum pooling */
void StoreShift()
{
    int orient, i, j, shift, ci, cj, si, sj, os, iS;
    float alpha;

    /* store all the possible shifts for all orientations */
    numShift = (locationPerturb * 2 + 1) * (orientationPerturb*2+1);
    rowShift = (int**)mxCalloc( nGaborScale * nGaborOri, sizeof(*rowShift) );
    colShift = (int**)mxCalloc( nGaborScale * nGaborOri, sizeof(*colShift) );
    orientShift = (int**)mxCalloc( nGaborScale * nGaborOri, sizeof(*orientShift) );

    for( iS = 0; iS < nGaborScale; ++iS ) /* for each scale */
    {
        for (orient=0; orient<nGaborOri; orient++)  /* for each orientation */
        {
            rowShift[orient+nGaborOri*iS] = (int*)mxCalloc( numShift, sizeof(**rowShift) );
            colShift[orient+nGaborOri*iS] = (int*)mxCalloc( numShift, sizeof(**colShift) );
            orientShift[orient+nGaborOri*iS] = (int*)mxCalloc( numShift, sizeof(**orientShift) );
            
            /* order the shifts from small to large */
            alpha = (float)PI*orient/nGaborOri;
            ci = 0;
            for( i=0; i<=locationPerturb; ++i ) // location shifts
            {
                for (si=0; si<=1; si++) /* sign of i */
                {
                    double location_shift = i;

                    cj = 0;
                    for( j = 0; j<=orientationPerturb; j++ ) // orientation shifts
                    {
                        for (sj=0; sj<=1; sj++) /* sign of j */
                        {
                            shift = ci*(2*orientationPerturb+1) + cj; /* hybrid counter for location and orientation */

                            //mexPrintf("locationshift=%.3f\n", location_shift);
                            //mexErrMsgTxt("for debug");
                            
                            rowShift[iS*nGaborOri+orient][shift] = (int)floor(.5+location_shift*(si*2-1)*cos(alpha));
                            colShift[iS*nGaborOri+orient][shift] = (int)floor(.5+location_shift*(si*2-1)*sin(alpha));
                            os = orient + j*(sj*2-1);
                            if (os<0)
                            {
                                os += nGaborOri;
                            }
                            else if (os>=nGaborOri)
                            {
                                os -= nGaborOri;
                            }
                            orientShift[iS*nGaborOri+orient][shift] = os;
                            if (j>0 || sj==1) /* triggers the orientation counter */
                            {
                                cj++;
                            }
                        }
                    }
                    if (i>0 || si==1)
                    {
                        ci++;
                    }
                }
            }
        }
        
    }
}

void Compute()
{
    int orient, x, y, iS;
    /* local maximum pooling at (x, y, orient) */
    float maxResponse, r;
    int shift, x1, y1, orient1, mshift;
    
    sizexSubsample = (int)floor((double)height / subsample);
    sizeySubsample = (int)floor((double)width / subsample);
    
    M1Map = (float**)mxCalloc( nGaborOri * nGaborScale, sizeof(*M1Map) );
    M1Trace = (int**)mxCalloc( nGaborOri * nGaborScale, sizeof(*M1Trace) );
    
    for( iS = 0; iS < nGaborScale; ++iS )
    {
        for (orient=0; orient<nGaborOri; orient++)
        {
            M1Map[iS*nGaborOri+orient] = (float*)mxCalloc( sizeySubsample * sizexSubsample, sizeof(**M1Map) );
            M1Trace[iS*nGaborOri+orient] = (int*)mxCalloc( sizeySubsample * sizexSubsample, sizeof(**M1Trace) );
            
            for (y=0; y<sizeySubsample; y++)
                for (x=0; x<sizexSubsample; x++)
                {
                    maxResponse = NEGMAX; mshift = 0;
                    for( shift=0; shift<numShift; shift++ )
                    {
                        x1 = x * subsample + rowShift[iS*nGaborOri+orient][shift];
                        y1 = y * subsample + colShift[iS*nGaborOri+orient][shift];
                        orient1 = orientShift[iS*nGaborOri+orient][shift];
                        if ((x1>=0)&&(x1<height)&&(y1>=0)&&(y1<width))
                        {
                            r = S1Map[iS*nGaborOri+orient1][x1+y1*height];
                            if (r>maxResponse)
                            {
                                maxResponse = r;
                                mshift = shift;
                            }
                        }
                    }
                    M1Map[iS*nGaborOri+orient][x+y*sizexSubsample] = maxResponse;
                    M1Trace[iS*nGaborOri+orient][x+y*sizexSubsample] = (float)mshift;
                }
        }
    }
    
}


/* entry point */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
    int ind, i, x, y, o, dataDim, j, bytes_to_copy, nGaborFilter;
    const mxArray *f;
    const mxArray *pAS1Map;
    mxArray *outPA;
    mwSize ndim;
    const mwSize* dims;
    mwSize dimsOutput[2];
    void* start_of_pr;
    mxClassID datatype;

	/*
	 * input variable 0: nGaborOri
	 */
	nGaborOri = (int)mxGetScalar(prhs[0]);
 
    /*
	 * input variable 1: S1 maps
	 */
    pAS1Map = prhs[1];
    dims = mxGetDimensions(pAS1Map);
    nGaborFilter = dims[0] * dims[1];
    nGaborScale = nGaborFilter / nGaborOri;
 
    S1Map = (const float**)mxCalloc( nGaborFilter, sizeof(*S1Map) );   /* SUM1 maps */
    for (i=0; i<nGaborFilter; ++i)
    {
        f = mxGetCell(pAS1Map, i);
        datatype = mxGetClassID(f);
        if (datatype != mxSINGLE_CLASS)
            mexErrMsgTxt("warning !! single precision required.");
        S1Map[i] = (const float*)mxGetPr(f);    /* get the pointer to cell content */
        height = mxGetM(f);    /* overwriting is ok, since it is constant */
        width = mxGetN(f);
    }
    
    /*
     * input variable 2: location shift radius
     */
    locationPerturb = (int)mxGetScalar(prhs[2]);
    
    /*
     * input variable 3: orientation shift radius
     */
    orientationPerturb = (int)mxGetScalar(prhs[3]);
    
    /*
     * input variable 4: sub sample step length
     */
    subsample = (int)mxGetScalar(prhs[4]);

    StoreShift();
    Compute();
    
    /* =============================================
     * Handle output variables.
     * ============================================= 
     */
    
    /*
     * output variable 0: M1 map
     */
    dimsOutput[0] = 1; dimsOutput[1] = nGaborScale * nGaborOri;
	plhs[0] = mxCreateCellArray( 2, dimsOutput );
    dimsOutput[0] = sizexSubsample; dimsOutput[1] = sizeySubsample;
    for( i = 0; i < nGaborOri * nGaborScale; ++i )
    {
        outPA = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
        /* populate the real part of the created array */
        start_of_pr = (float*)mxGetData(outPA);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(outPA);
        memcpy( start_of_pr, M1Map[i], bytes_to_copy );
        mxSetCell( plhs[0], i, outPA );
    }
    
    /*
     * output variable 1: M1 trace
     */
    dimsOutput[0] = 1; dimsOutput[1] = nGaborScale * nGaborOri;
	plhs[1] = mxCreateCellArray( 2, dimsOutput );
    dimsOutput[0] = sizexSubsample; dimsOutput[1] = sizeySubsample;
    for( i = 0; i < nGaborOri * nGaborScale; ++i )
    {
        outPA = mxCreateNumericArray( 2, dimsOutput, mxINT32_CLASS, mxREAL );
        /* populate the real part of the created array */
        start_of_pr = (int*)mxGetData(outPA);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(outPA);
        memcpy( start_of_pr, M1Trace[i], bytes_to_copy );
        mxSetCell( plhs[1], i, outPA );
    }
    
    /*
     * output variable 2: stored Gabor shifts : row shifts
     */
    dimsOutput[0] = nGaborScale * nGaborOri; dimsOutput[1] = 1;
    plhs[2] = mxCreateCellArray( 2, dimsOutput );
    for( i = 0; i < nGaborScale*nGaborOri; ++i )
    {
        dimsOutput[0] = numShift; dimsOutput[1] = 1;
        outPA = mxCreateNumericArray( 2, dimsOutput, mxINT32_CLASS, mxREAL );
        start_of_pr = (int*)mxGetData(outPA);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(outPA);
        memcpy( start_of_pr, rowShift[i], bytes_to_copy );
        mxSetCell( plhs[2], i, outPA );
    }
    
    /*
     * output variable 3: stored Gabor shifts : col shifts
     */
    dimsOutput[0] = nGaborScale * nGaborOri; dimsOutput[1] = 1;
    plhs[3] = mxCreateCellArray( 2, dimsOutput );
    for( i = 0; i < nGaborScale*nGaborOri; ++i )
    {
        dimsOutput[0] = numShift; dimsOutput[1] = 1;
        outPA = mxCreateNumericArray( 2, dimsOutput, mxINT32_CLASS, mxREAL );
        start_of_pr = (int*)mxGetData(outPA);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(outPA);
        memcpy( start_of_pr, colShift[i], bytes_to_copy );
        mxSetCell( plhs[3], i, outPA );
    }
    
    /*
     * output variable 4: stored Gabor shifts : orientation shifts
     */
    dimsOutput[0] = nGaborScale * nGaborOri; dimsOutput[1] = 1;
    plhs[4] = mxCreateCellArray( 2, dimsOutput );
    for( i = 0; i < nGaborScale*nGaborOri; ++i )
    {
        dimsOutput[0] = numShift; dimsOutput[1] = 1;
        outPA = mxCreateNumericArray( 2, dimsOutput, mxINT32_CLASS, mxREAL );
        start_of_pr = (int*)mxGetData(outPA);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(outPA);
        memcpy( start_of_pr, orientShift[i], bytes_to_copy );
        mxSetCell( plhs[4], i, outPA );
    }
}

