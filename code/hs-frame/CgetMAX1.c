# include <stdio.h>
# include <stdlib.h>
# include "mex.h"        
# include "math.h"
# define PI 3.1415926
# define ROUND(x) (floor((x)+.5))
# define NEGMAX -1e10
/* Generating integer vector */
int *int_vector(int n)
{
    int *v; 
    v = (int*) mxCalloc (n, sizeof(int));
    return v; 
}
/* Generating integer matrix */
int **int_matrix(int m, int n)
{
    int **mat; 
    int i; 
    mat = (int**) mxCalloc(m, sizeof(int*)); 
    for (i=0; i<m; i++)
        mat[i] = int_vector(n); 
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
int numOrient, locationShiftLimit, orientShiftLimit; /* key parameters */
int numImage; /* number of resolutions */   
double *allSizex, *allSizey;  /* sizes of images at multiple resolutions */
float **SUM1map, **MAX1map;           
int sizex, sizey, subsample, sizexSubsample, sizeySubsample; /* MAX1 maps are smaller than SUM1 maps by subsample */ 
int numShift, **xShift, **yShift, **orientShifted; /* stored shifts in x and y, and the shifted orientations */      
/* store all the shifts or perturbations in local maximum pooling */
void StoreShift()  
{
    int orient, i, j, shift, ci, cj, si, sj, os;
    double alpha; 
    /* store all the possible shifts for all orientations */
    numShift = (locationShiftLimit*2+1)*(orientShiftLimit*2+1); 
    xShift = int_matrix(numOrient, numShift); 
    yShift = int_matrix(numOrient, numShift); 
    orientShifted = int_matrix(numOrient, numShift);
    /* order the shifts from small to large */
    for (orient=0; orient<numOrient; orient++)        
    {
        alpha = PI*orient/numOrient;
        ci = 0; 
        for (i=0; i<=locationShiftLimit; i++)
           for (si=0; si<=1; si++)
           {
             cj = 0;    
             for (j = 0; j<=orientShiftLimit; j++)
                 for (sj=0; sj<=1; sj++)
                  {
                    shift = ci*(2*orientShiftLimit+1) + cj; 
                    xShift[orient][shift] = ROUND(i*(si*2-1)*subsample*cos(alpha)); 
                    yShift[orient][shift] = ROUND(i*(si*2-1)*subsample*sin(alpha)); 
                    os = orient + j*(sj*2-1);
                    if (os<0)
                        os += numOrient; 
                    else if (os>=numOrient)
                        os -= numOrient; 
                    orientShifted[orient][shift] = os; 
                    
                    if (j>0) 
                        cj ++; 
                    else if (sj==1)
                        cj++; 
                 }
             if (i>0)
                 ci ++; 
             else if (si==1)
                 ci++;                     
             }
    }
}
/* local maximum pooling at (x, y, orient) */
float LocalMaximumPooling(int img, int orient, int x, int y, int *trace) 
{
   float maxResponse, r;
   int shift, x1, y1, orient1, mshift; 
 
   maxResponse = NEGMAX; mshift = 0;  
   for (shift=0; shift<numShift; shift++)
   {
       x1 = x + xShift[orient][shift]; 
       y1 = y + yShift[orient][shift]; 
       orient1 = orientShifted[orient][shift];
       if ((x1>=1)&&(x1<=sizex)&&(y1>=1)&&(y1<=sizey))
        {            
           r = SUM1map[orient1*numImage+img][px(x1, y1, sizex, sizey)]; 
           if (r>maxResponse)
                {
                   maxResponse = r;   
                   mshift = shift; 
                }
         }
     }
   trace[0] = mshift;  /* return the shift in local maximum pooling */
   return(maxResponse); 
}
   
/* computer MAX1 maps for image img */
void ComputeMAX1map(int img)
{
   int orient, x, y, here, i, trace[2];  
   
   for (x=1; x<=sizexSubsample; x++)
      for (y=1; y<=sizeySubsample; y++)
       {
        here = px(x, y, sizexSubsample, sizeySubsample); 
        for (orient=0; orient<numOrient; orient++)
              MAX1map[orient*numImage+img][here] = LocalMaximumPooling(img, orient, x*subsample, y*subsample, trace); 
       }
}

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
 int img, orient, c; 
 mxArray *f;  
 
 c = 0; 
 numImage = ROUND(mxGetScalar(prhs[c++]));
 allSizex = mxGetPr(prhs[c++]);
 allSizey = mxGetPr(prhs[c++]);
 numOrient = ROUND(mxGetScalar(prhs[c++]));   
 locationShiftLimit = ROUND(mxGetScalar(prhs[c++]));    
 orientShiftLimit = ROUND(mxGetScalar(prhs[c++]));  
 subsample = ROUND(mxGetScalar(prhs[c++]));
 SUM1map = mxCalloc(numImage*numOrient, sizeof(float*));   
 for (img=0; img<numImage; img++)
     for (orient=0; orient<numOrient; orient++)
      {  
       f = mxGetCell(prhs[c], orient*numImage+img); 
       SUM1map[orient*numImage+img] = mxGetPr(f);       
      }
 c++;
 MAX1map = mxCalloc(numImage*numOrient, sizeof(float*));   
 for (img=0; img<numImage; img++)
   {
     for (orient=0; orient<numOrient; orient++)
      {  
       f = mxGetCell(prhs[c], orient*numImage+img); 
       MAX1map[orient*numImage+img] = mxGetPr(f);         
      }
    }
 c++;
 StoreShift();  
 for (img=0; img<numImage; img++)
 {
 sizex = ROUND(allSizex[img]); 
 sizey = ROUND(allSizey[img]); 
 sizexSubsample = floor((double)sizex/subsample); 
 sizeySubsample = floor((double)sizey/subsample); 
 ComputeMAX1map(img);  
 }
 free_matrix(xShift, numOrient, numShift); 
 free_matrix(yShift, numOrient, numShift); 
 free_matrix(orientShifted, numOrient, numShift);
}
 


     



 

                    