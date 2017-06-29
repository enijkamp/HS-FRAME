/********************************************************************************************
  The Match Pursuit Algorithm V2 is also to generate the deformed template for every training images
*******************************************/

# include "mex.h"        /* the algorithm is connect to matlab */
# include "math.h"
# include "matrix.h"
# include "stdlib.h"
# define PI 3.1415926
# define ABS(x) ((x)>0? (x):(-(x)))
# define MAX(x, y) ((x)>(y)? (x):(y))
# define MIN(x, y) ((x)<(y)? (x):(y))
# define NEGMAX -1e10



int N;                      /* number of filters */
double *syma;
double **allsymbol;
double* filterSelected;
double* xSelected;
double* ySelected;
int numSketchSelected;
int sx, sy;                 /* sizes of image and range of max pooling */ 
double *half;                      /* halfsizes of filters */

int px(int x, int y, int bx, int by)    
{            
   return (x + (y-1)*bx - 1); 
}

void draw(double *symo, int x0, int y0, int ind)
{
  int x, y, h; 
   
  h = floor(half[ind]+.5);
  for (x=MAX(1, x0-h); x<=MIN(sx, x0+h); x++)
     for (y=MAX(1, y0-h); y<=MIN(sy, y0+h); y++)
       {
          symo[px(x, y, sx, sy)] = MAX(symo[px(x, y, sx, sy)], allsymbol[ind][px(x-x0+h+1, y-y0+h+1, 2*h+1, 2*h+1)]);
       }
}

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
 int ind, i, j, c; 
 mxArray *f;   
 
 c=0; 
 N = floor(mxGetScalar(prhs[c++])+.5);  /* number of total filters*/
 half = mxGetPr(prhs[c++]); 
 sx = floor(mxGetScalar(prhs[c++])+.5); 
 sy = floor(mxGetScalar(prhs[c++])+.5); 
 
 allsymbol = (double**) mxCalloc(N, sizeof(double*));    
 for (ind=0; ind<N; ind++)
 {       
     f = mxGetCell(prhs[c], ind); 
     allsymbol[ind] = mxGetPr(f);        
 }
 c += 1; 
 syma = mxGetPr(prhs[c++]); 
  
 filterSelected = mxGetPr(prhs[c++]); 
 xSelected = mxGetPr(prhs[c++]); 
 ySelected = mxGetPr(prhs[c++]); 
 numSketchSelected =(int)mxGetScalar(prhs[c++]);
 
 for (i=0; i< numSketchSelected; i++){

     draw(syma, xSelected[i], ySelected[i], filterSelected[i]-1);        /* one comprehensive template*/
     
 } 
 

}
                 