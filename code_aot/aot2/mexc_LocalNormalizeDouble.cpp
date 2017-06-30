# include <stdio.h>
# include <stdlib.h>
# include "mex.h"        
# include "math.h"
# define ROUND(x) (floor((x)+.5))
# define NEGMAX -1e10
# define MAX(x, y) ((x)>(y)? (x):(y))
/* Generating double vector */
double *double_vector(int n)
{
    double *v; 
    v = (double*) mxCalloc (n, sizeof(double));    
    return v; 
}
/* Generating double matrix */
double **double_matrix(int m, int n)
{
    double **mat; 
    int i; 
    mat = (double**) mxCalloc(m, sizeof(double*)); 
    for (i=0; i<m; i++)
        mat[i] = double_vector(n); 
    return mat; 
}
/* Free matrix space */
void free_matrix(void **mat, int m, int n)
{
        int i;
        for (i=0; i<m; i++)
              mxFree(mat[i]);
        mxFree(mat);
}
/* Compute pixel index in the vector that stores image */
int px(int x, int y, int lengthx, int lengthy)  /* the image is lengthx*lengthy */
{            
   return (x + (y-1)*lengthx - 1); 
 }
 /* variables */
double **SUM1map, **SUM1mapAll, **integralMap, **averageMap, thresholdFactor;    
int numOrient, halfFilterSize, localHalfx, localHalfy, sizex, sizey;         
/* compute sigmoid transformation */
void LocalNormalize()
{
   int x, y, here, orient, leftx, rightx, upy, lowy, k, startx[8], endx[8], starty[8], endy[8], copyx[8], copyy[8], fx, fy;  
   double maxAverage, averageDivide; 
   /* compute the sum over all the orientations at each pixel */
   SUM1mapAll = double_matrix(sizex, sizey); 
   for (x=0; x<sizex; x++)
       for (y=0; y<sizey; y++)
       {
        SUM1mapAll[x][y] = 0.; 
        here = px(x+1, y+1, sizex, sizey); 
        for (orient=0; orient<numOrient; orient++)
           {    
                SUM1mapAll[x][y] += SUM1map[orient][here]; 
            }
       }
    /* compute the integral map */
    integralMap = double_matrix(sizex, sizey); 
    integralMap[0][0] = SUM1mapAll[0][0];
    for (x=1; x<sizex; x++)
        integralMap[x][0] = integralMap[x-1][0]+SUM1mapAll[x][0];
    for (y=1; y<sizey; y++)
        integralMap[0][y] = integralMap[0][y-1]+SUM1mapAll[0][y]; 
    for (x=1; x<sizex; x++)
       for (y=1; y<sizey; y++)
        integralMap[x][y] = integralMap[x][y-1]+integralMap[x-1][y]-integralMap[x-1][y-1]+SUM1mapAll[x][y]; 
    /* compute the local average around each pixel */
    averageMap = double_matrix(sizex, sizey); 
    leftx = halfFilterSize+localHalfx; rightx = sizex-halfFilterSize-localHalfx; 
    upy = halfFilterSize+localHalfy; lowy = sizey-halfFilterSize-localHalfy; 
    maxAverage = NEGMAX; 
    if ((leftx<rightx)&&(upy<lowy))
    {
    for (x=leftx; x<rightx; x++)
       for (y=upy; y<lowy; y++)
       {
            averageMap[x][y] = (integralMap[x+localHalfx][y+localHalfy] 
				- integralMap[x-localHalfx-1][y+localHalfy] - integralMap[x+localHalfx][y-localHalfy-1]
				+ integralMap[x-localHalfx-1][y-localHalfy-1])/(2.*localHalfx+1.)/(2.*localHalfy+1.)/numOrient; 
            if (maxAverage < averageMap[x][y])
                maxAverage = averageMap[x][y]; 
       }
    /* take care of the boundaries */
    k = 0;  
    /* four corners */
    startx[k] = 0; endx[k] = leftx; starty[k] = 0; endy[k] = upy; copyx[k] = leftx; copyy[k] = upy; k++; 
    startx[k] = 0; endx[k] = leftx; starty[k] = lowy; endy[k] = sizey; copyx[k] = leftx; copyy[k] = lowy-1; k++;  
    startx[k] = rightx; endx[k] = sizex; starty[k] = 0; endy[k] = upy; copyx[k] = rightx-1; copyy[k] = upy; k++; 
    startx[k] = rightx; endx[k] = sizex; starty[k] = lowy; endy[k] = sizey; copyx[k] = rightx-1; copyy[k] = lowy-1; k++; 
    /* four sides */
    startx[k] = 0; endx[k] = leftx; starty[k] = upy; endy[k] = lowy; copyx[k] = leftx; copyy[k] = -1; k++; 
    startx[k] = rightx; endx[k] = sizex; starty[k] = upy; endy[k] = lowy; copyx[k] = rightx-1; copyy[k] = -1; k++; 
    startx[k] = leftx; endx[k] = rightx; starty[k] = 0; endy[k] = upy; copyx[k] = -1; copyy[k] = upy; k++; 
    startx[k] = leftx; endx[k] = rightx; starty[k] = lowy; endy[k] = sizey; copyx[k] = -1; copyy[k] = lowy-1; k++; 
    /* propagate the average to the boundaries */
    for (k=0; k<8; k++) 
        for (x=startx[k]; x<endx[k]; x++)
            for (y=starty[k]; y<endy[k]; y++)
            {
                if (copyx[k]<0)
                    fx = x; 
                else 
                    fx = copyx[k]; 
                if (copyy[k]<0)
                    fy = y; 
                else 
                    fy = copyy[k]; 
                averageMap[x][y] = averageMap[fx][fy]; 
            }
     /* normalize the responses by local averages */
     for (x=0; x<sizex; x++)
       for (y=0; y<sizey; y++)
       {      
        here = px(x+1, y+1, sizex, sizey);
        averageDivide = MAX(averageMap[x][y], maxAverage*thresholdFactor); 
        for (orient=0; orient<numOrient; orient++)
                SUM1map[orient][here] /= averageDivide; 
       }
    }
    /* free intermediate matrices */
    free_matrix((void**)SUM1mapAll, sizex, sizey);  
    free_matrix((void**)integralMap, sizex, sizey); 
    free_matrix((void**)averageMap, sizex, sizey); 
}
/* load in variables */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
 int orient, c; 
 mxArray *f;   
 
 c = 0; 
 sizex = ROUND(mxGetScalar(prhs[c++])); 
 sizey = ROUND(mxGetScalar(prhs[c++])); 
 numOrient = ROUND(mxGetScalar(prhs[c++])); 
 halfFilterSize = ROUND(mxGetScalar(prhs[c++])); 
 localHalfx = ROUND(mxGetScalar(prhs[c++])); /* half size of the local area for averaging */
 localHalfy = ROUND(mxGetScalar(prhs[c++])); 
 SUM1map = (double**)mxCalloc(numOrient, sizeof(double*));  /* filter responses */
 for (orient=0; orient<numOrient; orient++)
      {  
       f = mxGetCell(prhs[c], orient); 
       SUM1map[orient] = (double*)mxGetPr(f);       
      }
 c++;
 thresholdFactor = mxGetScalar(prhs[c++]); 
 LocalNormalize(); 
}

     



 

                    
