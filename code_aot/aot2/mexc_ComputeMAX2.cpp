/* mex-C: 
 * compute MAX2 maps for a single image
 *
 * Usage:
 * 
 *    [M2Map M2LocationTrace M2TemplateTrace M2RowShift M2ColShift] = 
 *      CcomputeMAX2( templateAffinityMatrix, S2Map, locationPerturbationFraction, 
 *          sizeTemplate, subsampleM2 );
 *
 *     <templateAffinityMatrix> is a cell array, which stores a sparse version
 *      of the template affinity matrix. The array should be in int32 format. 
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
# define PI 3.1415926
# define NEGMAX -1e10

/* variable declaration */
float **M2Map;                  /* MAX 1 maps */
int **M2LocationTrace;          /* ARGMAX map for row and column shifts */
int **M2TemplateTrace;          /* ARGMAX map for template index. This index may be hybrid index of S2 template and its transform */
const float** S2Map;            /* SUM1 maps */
int height, width;              /* size of MAX1 maps */
int *sizeTemplate;               /* size of templates at different scales */
float *rowShift, *colShift;           /* store the location/orientation shifts per scale */
int sizexsubsampleM2, sizeysubsampleM2; /* (down sampled) size of M1 map */
int *numLocationShift;              /* number of location shifts per scale */
int subsampleM2;                      /* sub sampling step length */
float locationPerturbRadius;        /* allowed location perturbation, relative to the template size */
const int **templateAffMat;         /* indices of similar templates (starting from 0), i.e. template affinity matrix in sparse format */
int *nNeighborTemplate;             /* number of neighbor templates per template */
int nLocationShift;                 /* number of location shifts for each template */
int nTemplate;                      /* number of S2 templates */

/* store all the shifts or perturbations in local maximum pooling */
void StoreShift()
{
    int orient, i, j, shift, si, sj, os, iS;
    float dc, dr;
    float alpha;
    float stepsize;
    
    /* store all the possible shifts for all templates */
    stepsize = 0.05;
    
    /* For each template, the number of location shifts is the same */
    /* We store the shifts as relative positions to the template size. */
    
    /* count the number of location shifts in a lazy way */
    j = 0;
    for( dc = -locationPerturbRadius; dc <= locationPerturbRadius; dc += stepsize )
    {
        for( dr = -locationPerturbRadius; dr <= locationPerturbRadius; dr += stepsize )
        {
            ++j;
        }
    }
    nLocationShift = j;
    
    rowShift = (float*)mxCalloc( nLocationShift, sizeof(*rowShift) );
    colShift = (float*)mxCalloc( nLocationShift, sizeof(*colShift) );
    
    /* store the location shifts in a lazy way */
    j = 0;
    for( dc = -locationPerturbRadius; dc <= locationPerturbRadius; dc += stepsize )
    {
        for( dr = -locationPerturbRadius; dr <= locationPerturbRadius; dr += stepsize )
        {
            rowShift[j] = dr;
            colShift[j] = dc;
            ++j;
        }
    }
}

void Compute()
{
    int orient, x, y, iS;
    /* local maximum pooling at (x, y, iT) */
    float maxResponse, r;
    int shift, x1, y1, orient1, mshift;
    int iT, jT, jjT;
    int bestLocationshift, bestTemplate;
    bool success;
    float *M2Map_toUpdate, *M2Map_neighbor;
    
    sizexsubsampleM2 = (int)floor((double)height / subsampleM2);
    sizeysubsampleM2 = (int)floor((double)width / subsampleM2);
    
    M2Map = (float**)mxCalloc( nTemplate, sizeof(*M2Map) );
    M2LocationTrace = (int**)mxCalloc( nTemplate, sizeof(*M2LocationTrace) );
    M2TemplateTrace = (int**)mxCalloc( nTemplate, sizeof(*M2TemplateTrace) );

    /* For each channel, perform local maximization on its own. */
    for( jT = 0; jT < nTemplate; ++jT )
    {
        M2Map[jT] = (float*)mxCalloc( sizeysubsampleM2 * sizexsubsampleM2, sizeof(**M2Map) );
        M2LocationTrace[jT] = (int*)mxCalloc( sizeysubsampleM2 * sizexsubsampleM2, sizeof(**M2LocationTrace) );
        M2TemplateTrace[jT] = (int*)mxCalloc( sizeysubsampleM2 * sizexsubsampleM2, sizeof(**M2TemplateTrace) );

        for (y=0; y<sizeysubsampleM2; y++)
            for (x=0; x<sizexsubsampleM2; x++)
            {
                success = false;
                maxResponse = NEGMAX;

                for( shift=0; shift<nLocationShift; ++shift )
                {
                    x1 = (int)floor( (x) * subsampleM2 + rowShift[shift] * sizeTemplate[jT]);
                    y1 = (int)floor( (y) * subsampleM2 + colShift[shift] * sizeTemplate[jT]);
                    if ((x1>=0)&&(x1<height)&&(y1>=0)&&(y1<width))
                    {
                        r = S2Map[jT][x1+y1*height];
                        if (r>maxResponse)
                        {
                            success = true;
                            maxResponse = r;
                            bestLocationshift = shift;
                            bestTemplate = jT;
                        }
                    }
                }
                
                M2Map[jT][x+y*sizexsubsampleM2] = maxResponse;
                M2LocationTrace[jT][x+y*sizexsubsampleM2] = bestLocationshift;
                M2TemplateTrace[jT][x+y*sizexsubsampleM2] = bestTemplate;
            }
    }
    
    /* Then, combine the channels. */
    for( iT = 0; iT < nTemplate; ++iT )
    {
        /* update the MAX2 maps for template iT */
        for( jjT = 0; jjT < nNeighborTemplate[iT]; ++jjT )
        {
            jT = templateAffMat[iT][jjT];
            for (y=0; y<sizeysubsampleM2; y++)
                for (x=0; x<sizexsubsampleM2; x++)
                {
                    if( M2Map[jT][x+y*sizexsubsampleM2] >
                        M2Map[iT][x+y*sizexsubsampleM2] )
                    {
                        M2Map[iT][x+y*sizexsubsampleM2] = M2Map[jT][x+y*sizexsubsampleM2];
                        M2TemplateTrace[iT][x+y*sizexsubsampleM2] = M2TemplateTrace[jT][x+y*sizexsubsampleM2];
                        M2LocationTrace[iT][x+y*sizexsubsampleM2] = M2LocationTrace[jT][x+y*sizexsubsampleM2];
                    }
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
    const mxArray *pATemplateAffMat;
    const mxArray *pAS2Map;
    mxArray *outPA;
    mwSize ndim;
    const mwSize* dims;
    mwSize dimsOutput[2];
    void* start_of_pr;
    mxClassID datatype;

	/*
	 * input variable 0: templateAffinityMatrix
	 */
    pATemplateAffMat = prhs[0];
    nTemplate = mxGetM(pATemplateAffMat) * mxGetN(pATemplateAffMat);
    templateAffMat = (const int**)mxCalloc( nTemplate, sizeof(*templateAffMat) );
    nNeighborTemplate = (int*)mxCalloc(nTemplate,sizeof(*nNeighborTemplate));
    for( i = 0; i < nTemplate; ++i )
    {
        f = mxGetCell(pATemplateAffMat, i);
        datatype = mxGetClassID(f);
        if (datatype != mxINT32_CLASS)
            mexErrMsgTxt("warning !! int32 data type required.");
        nNeighborTemplate[i] = mxGetM(f) * mxGetN(f);
        templateAffMat[i] = (const int*)mxGetPr(f);
    }
 
    /*
	 * input variable 1: S2 maps
	 */
    pAS2Map = prhs[1];
    S2Map = (const float**)mxCalloc( nTemplate, sizeof(*S2Map) );   /* SUM1 maps */
    for (i=0; i<nTemplate; ++i)
    {
        f = mxGetCell(pAS2Map, i);
        datatype = mxGetClassID(f);
        if (datatype != mxSINGLE_CLASS)
            mexErrMsgTxt("warning !! single precision required.");
        S2Map[i] = (const float*)mxGetPr(f);    /* get the pointer to cell content */
        height = mxGetM(f);    /* overwriting is ok, since it is constant, or we would only input one image */
        width = mxGetN(f);
    }
    
    /*
     * input variable 2: location shift radius
     */
    locationPerturbRadius = (float)mxGetScalar(prhs[2]);

    /*
     * input variable 3: size of the templates
     */
    datatype = mxGetClassID(prhs[3]);
    if( datatype != mxINT32_CLASS )
        mexErrMsgTxt("warning !! int32 data type required for sizeTemplates");
    sizeTemplate = (int*)mxGetPr(prhs[3]);
    
    /*
     * input variable 4: sub sample step length
     */
    subsampleM2 = (int)mxGetScalar(prhs[4]);

    StoreShift();
    Compute();
    
    /* =============================================
     * Handle output variables.
     * ============================================= 
     */
    
    /*
     * output variable 0: M2 map
     */
    dimsOutput[0] = 1; dimsOutput[1] = nTemplate;
	plhs[0] = mxCreateCellArray( 2, dimsOutput );
    dimsOutput[0] = sizexsubsampleM2; dimsOutput[1] = sizeysubsampleM2;
    for( i = 0; i < nTemplate; ++i )
    {
        outPA = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
        /* populate the real part of the created array */
        start_of_pr = (float*)mxGetData(outPA);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(outPA);
        memcpy( start_of_pr, M2Map[i], bytes_to_copy );
        mxSetCell( plhs[0], i, outPA );
    }
    
    /*
     * output variable 1: M2 location trace
     */
    dimsOutput[0] = 1; dimsOutput[1] = nTemplate;
	plhs[1] = mxCreateCellArray( 2, dimsOutput );
    dimsOutput[0] = sizexsubsampleM2; dimsOutput[1] = sizeysubsampleM2;
    for( i = 0; i < nTemplate; ++i )
    {
        outPA = mxCreateNumericArray( 2, dimsOutput, mxINT32_CLASS, mxREAL );
        /* populate the real part of the created array */
        start_of_pr = (int*)mxGetData(outPA);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(outPA);
        memcpy( start_of_pr, M2LocationTrace[i], bytes_to_copy );
        mxSetCell( plhs[1], i, outPA );
    }
    
    /*
     * output variable 1: M2 other transformation trace
     */
    dimsOutput[0] = 1; dimsOutput[1] = nTemplate;
	plhs[2] = mxCreateCellArray( 2, dimsOutput );
    dimsOutput[0] = sizexsubsampleM2; dimsOutput[1] = sizeysubsampleM2;
    for( i = 0; i < nTemplate; ++i )
    {
        outPA = mxCreateNumericArray( 2, dimsOutput, mxINT32_CLASS, mxREAL );
        /* populate the real part of the created array */
        start_of_pr = (int*)mxGetData(outPA);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(outPA);
        memcpy( start_of_pr, M2TemplateTrace[i], bytes_to_copy );
        mxSetCell( plhs[2], i, outPA );
    }
    
    /*
     * output variable 2: stored template location shifts : row, col
     */
    dimsOutput[0] = nLocationShift; dimsOutput[1] = 2;
    outPA = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
    start_of_pr = (float*)mxGetData(outPA);
    bytes_to_copy = dimsOutput[0] * 1 * mxGetElementSize(outPA);
    memcpy( start_of_pr, rowShift, bytes_to_copy );
    start_of_pr = (float*)mxGetData(outPA) + nLocationShift;
    memcpy( start_of_pr, colShift, bytes_to_copy );
    plhs[3] = outPA;

}

