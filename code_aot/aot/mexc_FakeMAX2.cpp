/* mex-C: 
 * compute MAX2 maps for a single image
 *
 * Usage:
 * 
 *    [M2Map] = 
 *      mexc_FakeMAX2( S2Map, LocTrace, TransformTrace, templateAffinityMatrix, sizeTemplate );
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
const int **LocTrace;          /* ARGMAX map for row and column shifts */
const int **TransformTrace;          /* ARGMAX map for template index. This index may be hybrid index of S2 template and its transform */
const float** S2Map;            /* SUM1 maps */
const int *sizeTemplate;               /* size of templates at different scales */
int *heightPerMap;					/* sizes of feature maps */
int *widthPerMap;					/* sizes of feature maps */
const int **templateAffMat;         /* indices of similar templates (starting from 0), i.e. template affinity matrix in sparse format */
int *nNeighborTemplate;             /* number of neighbor templates per template */
int nTemplate;                      /* number of S2 templates */
float *rowColShift;           /* store the location/orientation shifts */
int nShift;
int nMap,numPartRotate,numCandPart,numResolution;

void Compute()
{
    M2Map = (float**) mxCalloc( nMap, sizeof(*M2Map) );
	for( int i = 0; i < nMap; ++i )
	{
		M2Map[i] = (float*)mxCalloc( heightPerMap[i] * widthPerMap[i], sizeof(**M2Map) );
	}
	
    for( int iRes = 0; iRes < numResolution; ++iRes )
    {
    	for( int iPart = 0; iPart < numCandPart; ++iPart )
        {
    		for( int iRot = 0; iRot < numPartRotate; ++iRot )
    		{
				int ind1 = iRot+iPart*numPartRotate+iRes*numCandPart*numPartRotate;
                
    			for( int y = 0; y < widthPerMap[ind1]; ++y )
    				for( int x = 0; x < heightPerMap[ind1]; ++x )
    				{
    					// use TransformTrace iPart and iRot to find the rotation
						int ind2 = x + y * heightPerMap[ind1];
						int indTransform = TransformTrace[ind1][ind2]; // the index for the transformed template (hybrid of iRot and iPart)
						int ind3 = indTransform + iRes * numCandPart * numPartRotate;
                        
						// use x, y and iRes to find the translation
						int indTranslation = LocTrace[ind1][ind2]; // LocTrace[ind3][ind2];
                        int row = (int) floor(  x + rowColShift[indTranslation] * (float)sizeTemplate[ind3] );
						int col = (int) floor(  y + rowColShift[indTranslation+nShift] * (float)sizeTemplate[ind3] );
						
    					// fetch the point in S2Map
                        float val;
                        if( row < 0 || row >= heightPerMap[ind3] || col < 0 || col >= widthPerMap[ind3] )
                        {
                            mexPrintf("out of bound: row=%d,col=%d,height=%d,width=%d\n",row,col,heightPerMap[ind3],widthPerMap[ind3]);
                            // mexErrMsgTxt( "out of bound \n" );
                            val = -10000;
                        }
						val = S2Map[ind3][row+col*heightPerMap[ind3]];
						M2Map[ind1][ind2] = val;
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
    const mxArray *f, *pA;
    const mxArray *pATemplateAffMat;
    const mxArray *pAS2Map;
    mxArray *outPA;
    mwSize ndim;
    const mwSize* dims;
    mwSize dimsOutput[3];
    void* start_of_pr;
    mxClassID datatype;


	/*
	 * input variable 0: S2 maps
	 */
    pAS2Map = prhs[0];
    const int* dimsInput = mxGetDimensions(pAS2Map);
    numPartRotate = dimsInput[0];
    numCandPart = dimsInput[1];
    numResolution = dimsInput[2];

    nMap = numPartRotate*numCandPart*numResolution;
    heightPerMap = (int*)mxCalloc( nMap, sizeof(int) );
    widthPerMap = (int*)mxCalloc( nMap, sizeof(int) );
    
    S2Map = (const float**)mxCalloc( nMap, sizeof(*S2Map) );   /* SUM1 maps */
    for (i=0; i<nMap; ++i)
    {
        f = mxGetCell(pAS2Map, i);
        datatype = mxGetClassID(f);
        if (datatype != mxSINGLE_CLASS)
            mexErrMsgTxt("warning !! single precision required.");
        S2Map[i] = (const float*)mxGetPr(f);    /* get the pointer to cell content */
        heightPerMap[i] = mxGetM( f );
        widthPerMap[i] = mxGetN( f );
    }
    
    
    /*
	 * input variable 1: MAX2LocTrace
	 */
    pA = prhs[1];
    LocTrace = (const int**)mxCalloc( nMap, sizeof(*LocTrace) );   /* SUM1 maps */
    for (i=0; i<nMap; ++i)
    {
        f = mxGetCell(pA, i);
        datatype = mxGetClassID(f);
        if (datatype != mxINT32_CLASS)
            mexErrMsgTxt("warning !! int32 precision required.");
        LocTrace[i] = (const int*)mxGetPr(f);    /* get the pointer to cell content */
    }
    
    /*
	 * input variable 2: MAX2TransformTrace
	 */
    pA = prhs[2];
    TransformTrace = (const int**)mxCalloc( nMap, sizeof(*TransformTrace) );   /* SUM1 maps */
    for (i=0; i<nMap; ++i)
    {
        f = mxGetCell(pA, i);
        datatype = mxGetClassID(f);
        if (datatype != mxINT32_CLASS)
            mexErrMsgTxt("warning !! int32 precision required.");
        TransformTrace[i] = (const int*)mxGetPr(f);    /* get the pointer to cell content */
    }
    

	/*
	 * input variable 3: templateAffinityMatrix
	 */
    pATemplateAffMat = prhs[3];
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
     * input variable 4: size of the templates
     */
    datatype = mxGetClassID(prhs[4]);
    if( datatype != mxINT32_CLASS )
        mexErrMsgTxt("warning !! int32 data type required for sizeTemplates");
    sizeTemplate = (const int*)mxGetPr(prhs[4]);
	
	/*
	* input variable 5: row and column shifts
	*/
	rowColShift = (float*)mxGetPr( prhs[5] );
	nShift = mxGetM( prhs[5] ) * mxGetN( prhs[5] ) / 2;
    

    Compute();
    
    /* =============================================
     * Handle output variables.
     * ============================================= 
     */
    
    /*
     * output variable 0: M2 map
     */
    dimsOutput[0] = numPartRotate; dimsOutput[1] = numCandPart; dimsOutput[2] = numResolution;
	plhs[0] = mxCreateCellArray( 3, dimsOutput );
    for( int i = 0; i < nMap; ++i )
    {
		dimsOutput[0] = heightPerMap[i]; dimsOutput[1] = widthPerMap[i];
        outPA = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
        /* populate the real part of the created array */
        start_of_pr = (float*)mxGetData(outPA);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(outPA);
        memcpy( start_of_pr, M2Map[i], bytes_to_copy );
        mxSetCell( plhs[0], i, outPA );
    }
    
}

