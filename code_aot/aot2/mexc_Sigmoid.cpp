/*
* Sigmoid transformation on S1 maps (pixels are updated in place). There is no output from this function. 
*
*
*/

# include <stdio.h>
# include <stdlib.h>
# include "mex.h"        
# include "math.h"
# define ROUND(x) (floor((x)+.5))
/* Compute pixel index in the vector that stores image */

 /* variables */
float **map;
int nMap, nPixel;
float saturation; /* saturation level for sigmoid transformation */

/* compute sigmoid transformation */
void SigmoidTransform(int i)
{
    for (int j=0; j<nPixel; ++j)
    {
        map[i][j] = saturation*(2./(1.+exp(-2.*map[i][j]/saturation))-1.);
    }
}

/* the entry function of this file */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
    int i;
    const mxArray *f;

    /* input 0: saturation value */
    saturation = mxGetScalar(prhs[0]);

    /* input 1: S1 maps */
    f = prhs[1];
    nMap = mxGetM( f ) * mxGetN( f );
    map = (float**)mxCalloc(nMap, sizeof(float*));
    for (i=0; i<nMap; i++)
    {
        mxArray* p;
        p = mxGetCell(f, i);
        map[i] = (float*)mxGetPr(p);
        nPixel = mxGetM(p) * mxGetN(p);
        SigmoidTransform(i);
    }


}





