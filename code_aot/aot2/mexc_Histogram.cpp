# include <stdio.h>
# include <stdlib.h>
# include "mex.h"       
# include "math.h"
# define ROUND(x) (floor((x)+.5))
# define MIN(x, y) ((x)<(y)? (x):(y))
/* Compute pixel index in the vector that stores image */
int px(int x, int y, int lengthx, int lengthy)  /* the image is lengthx*lengthy */
{            
   return (x + (y-1)*lengthx - 1); 
 }
/* variables */
int numImage, sizex, sizey, numOrient, halfFilterSize, numBin;                 
double **fI, binSize, *histog, saturation;                
/* compute sigmoid transformation */
void SigmoidTransform()
{
   int x, y, here, orient, img, i; 
   
   for (x=1; x<=sizex; x++)
       for (y=1; y<=sizey; y++)
       {
        here = px(x, y, sizex, sizey); 
        for (orient=0; orient<numOrient; orient++)
           {    
             for (img=0; img<numImage; img++)
             {
                i = orient*numImage+img; 
                fI[i][here] = saturation*(2./(1.+exp(-2.*fI[i][here]/saturation))-1.); 
             }
            }
       }
}
/* compute histogram of q() */
void Histogramq()
{
   int orient, img, x, y, here, b, tot; 
   
   tot = 0; 
   for (x=halfFilterSize+1; x<sizex-halfFilterSize; x++)
      for (y=halfFilterSize+1; y<sizey-halfFilterSize; y++)
       {
        here = px(x, y, sizex, sizey); 
        for (orient=0; orient<numOrient; orient++)
           {    
             for (img=0; img<numImage; img++)
             {
                 b = MIN(floor(fI[orient*numImage+img][here]/binSize), numBin-1); 
                 histog[b] += 1.; 
                 tot ++; 
             }
        }
      }
   for (b=0; b<numBin; b++)
       histog[b] /= (binSize*tot);
}

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
 int orient, img, c; 
 mxArray *f;  
 
 c = 0; 
 numImage = ROUND(mxGetScalar(prhs[c++]));  /* number of images */
 numOrient = ROUND(mxGetScalar(prhs[c++]));  /* number of orientations */
 fI = (double**)mxCalloc(numImage*numOrient, sizeof(double*));  /* fitered images */
 for (img=0; img<numImage; img++)
   {
     for (orient=0; orient<numOrient; orient++)
      {  
       f = mxGetCell(prhs[c], orient*numImage+img); 
       fI[orient*numImage+img] = mxGetPr(f);  /* get pointers to filtered images */       
      }
    }
 c++; 
 halfFilterSize = ROUND(mxGetScalar(prhs[c++]));  /* half size of filters */
 sizex = ROUND(mxGetScalar(prhs[c++]));
 sizey = ROUND(mxGetScalar(prhs[c++])); /* size of images */
 binSize = mxGetScalar(prhs[c++]); /* size of bin */
 numBin = ROUND(mxGetScalar(prhs[c++])); /* number of bins */
 histog = mxGetPr(prhs[c++]);  /* histogram of q() */
 saturation = mxGetScalar(prhs[c++]); /* saturation level of sigmoid transform */
 
 SigmoidTransform();
 Histogramq();
}

     



 

                    