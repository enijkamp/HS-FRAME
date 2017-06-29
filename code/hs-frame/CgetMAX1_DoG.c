# include <stdio.h>
# include <stdlib.h>
# include "mex.h"        
# include "math.h"
# define PI 3.1415926
# define ROUND(x) (floor((x)+.5))
# define NEGMAX -1e10

/* Compute pixel index in the vector that stores image */
int px(int x, int y, int lengthx, int lengthy)  /* the image is lengthx*lengthy */
{            
   return (x + (y-1)*lengthx - 1); 
 }
 /* variables */
int locationShiftLimit; /* key parameters */
int numImage; /* number of resolutions */   
double *allSizex, *allSizey;  /* sizes of images at multiple resolutions */
float **SUM1map, **MAX1map;           
int sizex, sizey, subsample, sizexSubsample, sizeySubsample; /* MAX1 maps are smaller than SUM1 maps by subsample */ 

/* local maximum pooling at (x, y) */
float LocalMaximumPoolingDoG(int img, int x, int y, int *trace) 
{
   float maxResponse, r;
   int x1, y1, lx, ly, mshift,shift; 
 
   maxResponse = NEGMAX;  
   shift=0; mshift=0;
   for (lx=-locationShiftLimit; lx<=locationShiftLimit; lx++)
     for (ly=-locationShiftLimit; ly<=locationShiftLimit; ly++)
     {
       x1 = x + lx; 
       y1 = y + ly; 
       shift=shift+1;
       if ((x1>=1)&&(x1<=sizex)&&(y1>=1)&&(y1<=sizey))
        {            
           r = SUM1map[img][px(x1, y1, sizex, sizey)]; 
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
          MAX1map[img][here] = LocalMaximumPoolingDoG(img, x*subsample, y*subsample, trace); 
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
 
 locationShiftLimit = ROUND(mxGetScalar(prhs[c++]));    
 
 subsample = ROUND(mxGetScalar(prhs[c++]));
 SUM1map = mxCalloc(numImage, sizeof(float*));   
 for (img=0; img<numImage; img++)
 {  
    f = mxGetCell(prhs[c], img); 
    SUM1map[img] = mxGetPr(f);       
 }
 c++;
 MAX1map = mxCalloc(numImage, sizeof(float*));   
 for (img=0; img<numImage; img++)
 {
    f = mxGetCell(prhs[c], img); 
    MAX1map[img] = mxGetPr(f);         
    
 }
 c++;
 
 for (img=0; img<numImage; img++)
 {
 sizex = ROUND(allSizex[img]); 
 sizey = ROUND(allSizey[img]); 
 sizexSubsample = floor((double)sizex/subsample); 
 sizeySubsample = floor((double)sizey/subsample); 
 ComputeMAX1map(img);  
 }
 
}
 


     


