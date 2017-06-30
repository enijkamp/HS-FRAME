/* 
 * Affine transform of a list of components. A component is a 4D point with its position, orientation and scale (length). 
 *
 * Usage:
 *    [outRow, outCol, outO, outS] =
            mexc_TemplateAffineTransform(tScale,rScale,cScale,rotation,inRow,inCol,inO,inS,nOri)
 *
 * Note: outO may be out of range. 
 * 
 */
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "mex.h"        
# include "math.h"
# define PI 3.1415926
# define ABS(x) ((x)>0? (x):(-(x)))
# define MAX(x, y) ((x)>(y)? (x):(y))
# define MIN(x, y) ((x)<(y)? (x):(y))
# define ROUND(x) (floor((x)+.5))

float **A1, **A2;       /* transformation matrices */
int nElement;               /* number of edgelet elements to be transformed */
float *outRow, *outCol, *outO, *outS;
float *inRow, *inCol, *inO, *inS;
float rotation, tScale, rScale, cScale;
int nOri;

void matrixMultiplication(int mA, int nA, float **inA,
        int nB, float **inB, float **outC)
{
    int c, r, k;
    float s;
    for( c = 0; c < nB; ++c )
        for( r = 0; r < mA; ++r )
        {
            s = 0;
            for( k = 0; k < nA; ++k )
            {
                s += inA[r][k] * inB[k][c];
            }
            outC[r][c] = s;
        }
}

void compute()
{
    int i, j;
    float angle;
    float **pt, **tmp;
    
    /* pt is 1 by 3 */
    pt = (float**)mxCalloc(1,sizeof(*pt));
    for( i = 0; i < 1; ++i )
    {
        pt[i] = (float*)mxCalloc(3,sizeof(**pt));
    }
    /* tmp is 1 by 3 */
    tmp = (float**)mxCalloc(1,sizeof(*tmp));
    for( i = 0; i < 1; ++i )
    {
        tmp[i] = (float*)mxCalloc(3,sizeof(**tmp));
    }
    
    /* set the transformation matrix */
    A1 = (float**)mxCalloc(3,sizeof(*A1));
    for( i = 0; i < 3; ++i )
    {
        A1[i] = (float*)mxCalloc(3,sizeof(**A1));
    }
    A2 = (float**)mxCalloc(3,sizeof(*A2));
    for( i = 0; i < 3; ++i )
    {
        A2[i] = (float*)mxCalloc(3,sizeof(**A2));
    }
    for( i = 0; i < 3; ++i )
        for( j = 0; j < 3; ++j )
        {
            A1[i][j] = 0;
            A2[i][j] = 0;
        }
    A1[0][0] = cScale * pow( (double)2, (double)(tScale/2) );
    A1[1][1] = rScale * pow( (double)2, (double)(tScale/2) );

    angle = rotation * PI/nOri;
    
    A2[0][0] = cos(angle); 
    A2[0][1] = sin(angle);
    A2[1][0] = -A2[0][1];
    A2[1][1] = A2[0][0];
        
    for( i = 0; i < nElement; ++i )
    {
        tmp[0][0] = inCol[i];
        tmp[0][1] = -inRow[i];
        tmp[0][2] = 1;
        matrixMultiplication(1,3,tmp,3,A1,pt);
        tmp[0][0] = pt[0][0];
        tmp[0][1] = pt[0][1];
        tmp[0][2] = pt[0][2];
        matrixMultiplication(1,3,tmp,3,A2,pt);
        
        outCol[i] = ROUND( pt[0][0] );
        outRow[i] = ROUND( -pt[0][1] );
        outO[i] = inO[i] + rotation;
        outS[i] = inS[i] + tScale;
        
        /*
        while( outO[i] < 0 )
            outO[i] = outO[i] + nOri;
        while( outO[i] >= nOri )
            outO[i] = outO[i] - nOri;
        */
    }
}


/* entry point: load input variables and run the algorithm */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
    mxArray *pA;
    mwSize dimsOutput[2];
    void* start_of_pr;
    mxClassID datatype;
    int bytes_to_copy;
    /* =============================================
     * Handle input variables.
     * ============================================= 
     */
    /*
	 * input variable 0: tScale
	 */
	tScale = (float)mxGetScalar(prhs[0]);
    /*
	 * input variable 1: rScale
	 */
	rScale = (float)mxGetScalar(prhs[1]);
    /*
	 * input variable 2: cScale
	 */
	cScale = (float)mxGetScalar(prhs[2]);
    /*
	 * input variable 3: rotation
	 */
	rotation = (float)mxGetScalar(prhs[3]);
    /*
	 * input variable 4: inRow
	 */
	inRow = (float*)mxGetPr(prhs[4]);
    nElement = mxGetM(prhs[4]) * mxGetN(prhs[4]);
    datatype = mxGetClassID(prhs[4]);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    /*
	 * input variable 5: inCol
	 */
	inCol = (float*)mxGetPr(prhs[5]);
    datatype = mxGetClassID(prhs[5]);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    /*
	 * input variable 6: inO
	 */
	inO = (float*)mxGetPr(prhs[6]);
    datatype = mxGetClassID(prhs[6]);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    /*
	 * input variable 7: inS
	 */
	inS = (float*)mxGetPr(prhs[7]);
    datatype = mxGetClassID(prhs[7]);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    /*
	 * input variable 8: nOri  (number of orientation levels within 0 ~ pi) 
	 */
	nOri = (float)mxGetScalar(prhs[8]);
    
    
    /* =============================================
     * Computation.
     * ============================================= 
     */
    outRow = (float*)mxCalloc( nElement, sizeof(*outRow) );
    outCol = (float*)mxCalloc( nElement, sizeof(*outCol) );
    outO = (float*)mxCalloc( nElement, sizeof(*outO) );
    outS = (float*)mxCalloc( nElement, sizeof(*outS) );
    compute();
    
    /* =============================================
     * Handle output variables.
     * ============================================= 
     */
    
    /*
     * output variable 0: outRow
     */
    dimsOutput[0] = nElement; dimsOutput[1] = 1;
	pA = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
    /* populate the real part of the created array */
    start_of_pr = (float*)mxGetData(pA);
    bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(pA);
    memcpy( start_of_pr, outRow, bytes_to_copy );
    plhs[0] = pA;
    
    /*
     * output variable 1: outCol
     */
    dimsOutput[0] = nElement; dimsOutput[1] = 1;
	pA = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
    /* populate the real part of the created array */
    start_of_pr = (float*)mxGetData(pA);
    bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(pA);
    memcpy( start_of_pr, outCol, bytes_to_copy );
    plhs[1] = pA;
    
    /*
     * output variable 2: outO
     */
    dimsOutput[0] = nElement; dimsOutput[1] = 1;
	pA = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
    /* populate the real part of the created array */
    start_of_pr = (float*)mxGetData(pA);
    bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(pA);
    memcpy( start_of_pr, outO, bytes_to_copy );
    plhs[2] = pA;
    
    /*
     * output variable 3: outS
     */
    dimsOutput[0] = nElement; dimsOutput[1] = 1;
	pA = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
    /* populate the real part of the created array */
    start_of_pr = (float*)mxGetData(pA);
    bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(pA);
    memcpy( start_of_pr, outS, bytes_to_copy );
    plhs[3] = pA;
    
}



