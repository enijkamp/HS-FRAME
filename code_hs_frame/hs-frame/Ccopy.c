# include <stdio.h>
# include <stdlib.h>
# include "mex.h"        
# include "math.h"
# define ROUND(x) (floor((x)+.5))
# define px(x, y, lengthx, lengthy) ((x) + ((y)-1)*(lengthx) - 1) 

 /* variables */
float *I, *J;    
double theta, sintheta, costheta; 
int sx, sy, sx1, sy1, Fx, Fy, Cx, Cy; 

/* copy images */
void CopyImage()
{
   int x, y, here, there, x1, y1; 
   double dx, dy, sintheta, costheta; 
   
   sintheta = sin(theta); costheta = cos(theta);  
   for (x=1; x<=sx; x++)
       for (y=1; y<=sy; y++)
       {
        here = px(x, y, sx, sy); 
        dx = x-Cx; dy = y-Cy; 
        x1 = ROUND(Fx + dx*costheta + dy*sintheta); 
        y1 = ROUND(Fy - dx*sintheta + dy*costheta); 
        if ((x1>=1)&&(x1<=sx1)&&(y1>=1)&&(y1<=sy1))
        {
           there = px(x1, y1, sx1, sy1); 
           I[here] = J[there];
        }
        else
            I[here] = 0.; 
        }
}

/* read in images */  
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
 int c; 
 
 c = 0;
 I = mxGetPr(prhs[c++]);
 J = mxGetPr(prhs[c++]);
 Fx = ROUND(mxGetScalar(prhs[c++]));
 Fy = ROUND(mxGetScalar(prhs[c++])); 
 Cx = ROUND(mxGetScalar(prhs[c++]));
 Cy = ROUND(mxGetScalar(prhs[c++]));
 sx = ROUND(mxGetScalar(prhs[c++]));
 sy = ROUND(mxGetScalar(prhs[c++]));   
 sx1 = ROUND(mxGetScalar(prhs[c++]));
 sy1 = ROUND(mxGetScalar(prhs[c++])); 
 theta = mxGetScalar(prhs[c++]); 

 CopyImage();
}

     



 

                    