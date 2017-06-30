/* SHARED SKETCH ALGORITHM FOR LEARNING ACTIVE BASIS */
# include <stdio.h>
# include <stdlib.h>
# include "mex.h"        
# include "math.h"
# define PI 3.1415926
# define ABS(x) ((x)>0? (x):(-(x)))
# define MAX(x, y) ((x)>(y)? (x):(y))
# define MIN(x, y) ((x)<(y)? (x):(y))
# define ROUND(x) (floor((x)+.5))
# define NEGMAX -1e10
/* Generating float vector */
float *float_vector(int n)
{
    float *v; 
    v = (float*) mxCalloc (n, sizeof(float));    
    return v; 
}
/* Generating integer vector */
int *int_vector(int n)
{
    int *v; 
    v = (int*) mxCalloc (n, sizeof(int));
    return v; 
}
/* Generating float matrix */
float **float_matrix(int m, int n)
{
    float **mat; 
    int i; 
    mat = (float**) mxCalloc(m, sizeof(float*)); 
    for (i=0; i<m; i++)
        mat[i] = float_vector(n); 
    return mat; 
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
/* key input and output variables */
int numOrient, locationShiftLimit, orientShiftLimit; /* key parameters */
int numImage, numElement; /* number of images and number of Gabors */   
float **SUM1map, **MAX1map, **pooledMax1map; /* pooledMax1map is summed over all images */ 
int **trackMap; /* track the shift or perturbation at each location and orientation */
double **Correlation, **allSymbol; /* correlation between Gabors, and symbols of Gabors */              
int halfFilterSize; /* filter size = 2*halfFilterSize + 1 */ 
int sizex, sizey, subsample, sizexSubsample, sizeySubsample; /* MAX1 maps are smaller than SUM1 maps by subsample */ 
int numShift, **xShift, **yShift, **orientShifted; /* stored shifts in x and y, and the shifted orientations */      
int numStoredPoint; /* number of stored points of lambda in exponential model */                     
double *storedlambda, *storedExpectation, *storedLogZ; /* stored lambda, mu, logZ in exponential model */    
int startx, endx, starty, endy; /* search within the interior of the image */
double *selectedOrient, *selectedx, *selectedy, *selectedlambda, *selectedLogZ; /* parameters of learned active basis */  
float *commonTemplate, **deformedTemplate; /* templates of input images and template matching scores */   

/* local maximum pooling at (x, y, orient) */
float LocalMaximumPooling(int img, int orient, int x, int y, int *trace) 
{
   float maxResponse, r;
   int shift, x1, y1, orient1, mshift; 
 
   maxResponse = NEGMAX;  mshift = 0;
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
/* the local maximal Gabor inhibits overlapping Gabors */
void NonMaximumSuppression(int img, int mo, int mx, int my) 
{
   int x, y, orient, x1, y1, orient1, i, here, shift, startx0, endx0, starty0, endy0, trace[2]; 
   float *f, maxResponse, maxResponse1; 
   double *fc; 
   /* inhibit on the SUM1 maps */
   for (orient=0; orient<numOrient; orient++)   
     {
       f = SUM1map[orient*numImage+img];   
       fc = Correlation[mo+orient*numOrient];   
       for (x=MAX(1, mx-2*halfFilterSize); x<=MIN(sizex, mx+2*halfFilterSize); x++)
         for (y=MAX(1, my-2*halfFilterSize); y<=MIN(sizey, my+2*halfFilterSize); y++)
         {
          f[px(x, y, sizex, sizey)] *= 
            fc[px(x-mx+2*halfFilterSize+1, y-my+2*halfFilterSize+1, 4*halfFilterSize+1, 4*halfFilterSize+1)];
         }
      }
   /* update the MAX1 maps */
   startx0 = floor((double)(mx-2*halfFilterSize)/subsample)-locationShiftLimit+1; 
   starty0 = floor((double)(my-2*halfFilterSize)/subsample)-locationShiftLimit+1; 
   endx0 =   floor((double)(mx+2*halfFilterSize)/subsample)+locationShiftLimit; 
   endy0 =   floor((double)(my+2*halfFilterSize)/subsample)+locationShiftLimit;

   for (orient=0; orient<numOrient; orient++)   
     {
      i = orient*numImage+img; 
      for (x=MAX(startx, startx0); x<=MIN(endx, endx0); x++)
         for (y=MAX(starty, starty0); y<=MIN(endy, endy0); y++)
         { /* go over the locations that may be affected */         
           here = px(x, y, sizexSubsample, sizeySubsample);        
           if( here < 0 || here >= sizexSubsample*sizeySubsample )
           {
               mexPrintf("error ! aaa\n");
           }
           if( i < 0 || i >= numImage*numOrient )
           {
               mexPrintf("error ! bbb\n");
           }
           maxResponse = MAX1map[i][here]; 
           shift = ROUND(trackMap[i][here]);
           if( shift < 0 || shift >= numShift )
           {
               mexPrintf("error ! sss\n");
               mexPrintf("shift=%d, numshift=%d\n",shift,numShift);
           }
           orient1 = orientShifted[orient][shift];    
           x1 = x*subsample + xShift[orient][shift]; 
           y1 = y*subsample + yShift[orient][shift];
           
           if ((x1-mx>=-2*halfFilterSize)&&(x1-mx<=2*halfFilterSize)&&
               (y1-my>=-2*halfFilterSize)&&(y1-my<=2*halfFilterSize)) 
             { /* if the previous local maximum is within the inhibition range */
                 
                 
              if(Correlation[mo+orient1*numOrient]
                [px(x1-mx+2*halfFilterSize+1,y1-my+2*halfFilterSize+1,4*halfFilterSize+1,4*halfFilterSize+1)]==0.) 
                 {   /* if it is indeed inhibited */
                     maxResponse1 = LocalMaximumPooling(img, orient, x*subsample, y*subsample, trace);
                     trackMap[i][here] = trace[0]; 
                     MAX1map[i][here] = maxResponse1;
                     pooledMax1map[orient][here] += (maxResponse1-maxResponse); 
                 }
             }
         }
      }
}
/* plot the bar for Gabor at (mx, my, mo) */
void DrawElement(float *Template, int mo, int mx, int my, double w) 
{
  int x, y, here; 
  float a;
          
  for (x=mx-halfFilterSize; x<=mx+halfFilterSize; x++)
     for (y=my-halfFilterSize; y<=my+halfFilterSize; y++)
      if ((x>=1)&&(x<=sizex)&&(y>=1)&&(y<=sizey)) 
        {
         a = allSymbol[mo][px(x-mx+halfFilterSize+1, y-my+halfFilterSize+1, 
                              2*halfFilterSize+1, 2*halfFilterSize+1)]*w; 
         here = px(x, y, sizex, sizey); 
         if (Template[here]<a)
             Template[here] = a;
        }
}
/* Initialize MAX1 maps */
void Initialize()
{
   int img, orient, x, y, here, i, trace[2];  
   
   pooledMax1map = float_matrix(numOrient, sizexSubsample*sizeySubsample);  
   startx = floor(((double)halfFilterSize+1)/subsample)+1+locationShiftLimit; 
   endx = floor((double)(sizex-halfFilterSize)/subsample)-1-locationShiftLimit; 
   starty = floor((double)(halfFilterSize+1)/subsample)+1+locationShiftLimit; 
   endy = floor((double)(sizey-halfFilterSize)/subsample)-1-locationShiftLimit; 
   for (x=startx; x<=endx; x++)
      for (y=starty; y<=endy; y++)
       {
        here = px(x, y, sizexSubsample, sizeySubsample); 
        for (orient=0; orient<numOrient; orient++)
           {    
             pooledMax1map[orient][here] = 0.; 
             for (img=0; img<numImage; img++)
                pooledMax1map[orient][here] += MAX1map[orient*numImage+img][here];     
            }
       }
}
/* Shared sketch algorithm for pursuing Gabors */
void PursueElement()
{
   int img, orient, x, y, besto, bestx, besty, mo, mx, my, t, here, shift, i, j, trace[2]; 
   float r, maxPooled, maxResponse, average, overShoot; 

   t = 0; /* t is for iterations */
   do
   { /* select the next Gabor */
     maxPooled = NEGMAX;
     for (x=startx; x<=endx; x++)
      for (y=starty; y<=endy; y++)
       {
        here = px(x, y, sizexSubsample, sizeySubsample); 
        for (orient=0; orient<numOrient; orient++)
           {
            r = pooledMax1map[orient][here]; 
            if (maxPooled<r)
             {
               maxPooled = r; 
               besto = orient; bestx = x; besty = y; 
              }
           }
       }
     /* estimate the parameter of exponential model */
     selectedOrient[t] = besto; selectedx[t] = bestx; selectedy[t] = besty; 
     average = maxPooled/numImage;    
     j = numStoredPoint-1; 
     while (storedExpectation[j]>average)
         j--; 
     if (j==numStoredPoint-1)
      {
        selectedlambda[t] = storedlambda[j]; 
        selectedLogZ[t] = storedLogZ[j]; 
       }
     else 
      { /* linear interpolation */
        overShoot = (average-storedExpectation[j])/(storedExpectation[j+1]-storedExpectation[j]); 
        selectedlambda[t] = storedlambda[j]+(storedlambda[j+1]-storedlambda[j])*overShoot; 
        selectedLogZ[t] = storedLogZ[j]+(storedLogZ[j+1]-storedLogZ[j])*overShoot; 
      }
     /* plot selected and perturbed Gabor and inhibit nearby Gabors */
     DrawElement(commonTemplate, besto, bestx*subsample, besty*subsample, sqrt(average));
     here = px(bestx, besty, sizexSubsample, sizeySubsample); 
     for (img=0; img<numImage; img++)
        {
          i = besto*numImage+img; 
          maxResponse = MAX1map[i][here]; 
          shift = ROUND(trackMap[i][here]); 
          mo = orientShifted[besto][shift];      
          mx = bestx*subsample + xShift[besto][shift]; 
          my = besty*subsample + yShift[besto][shift];      
         if (maxResponse>0.)
            {
               DrawElement(deformedTemplate[img], mo, mx, my, sqrt(maxResponse)); 
               NonMaximumSuppression(img, mo, mx, my);  
             }
         }
      t++; 
   }
  while (t<numElement);  
}

/* read in input variables and run the algorithm */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
 int orient, img, i, j, c, x, y; 
 mxArray *f;
 mxClassID datatype;
 
 c = 0; 
 /* about active basis */
 numOrient = ROUND(mxGetScalar(prhs[c++]));  
 locationShiftLimit = ROUND(mxGetScalar(prhs[c++]));    
 orientShiftLimit = ROUND(mxGetScalar(prhs[c++])); 
 subsample = ROUND(mxGetScalar(prhs[c++]));   
 numElement = ROUND(mxGetScalar(prhs[c++])); 
 /* about input images */
 numImage = ROUND(mxGetScalar(prhs[c++])); 
 sizex = ROUND(mxGetScalar(prhs[c++])); 
 sizey = ROUND(mxGetScalar(prhs[c++]));   
 SUM1map = (float**)mxCalloc(numImage*numOrient, sizeof(float*));  
 for (img=0; img<numImage; img++)
     for (orient=0; orient<numOrient; orient++)
     {  
       i = orient*numImage+img; 
       f = mxGetCell(prhs[c], i);
       datatype = mxGetClassID(f);
        if (datatype != mxSINGLE_CLASS)
            mexErrMsgTxt("warning !! single precision required.");
       SUM1map[i] = (float*)mxGetPr(f);    
    }
 c++; 
 MAX1map = (float**)mxCalloc(numImage*numOrient, sizeof(float*));  
 for (img=0; img<numImage; img++)
     for (orient=0; orient<numOrient; orient++)
      {
       i = orient*numImage+img; 
       f = mxGetCell(prhs[c], i);
       datatype = mxGetClassID(f);
        if (datatype != mxSINGLE_CLASS)
            mexErrMsgTxt("warning !! single precision required.");
       MAX1map[i] = (float*)mxGetPr(f);    
      }
 c++; 
 trackMap = (int**)mxCalloc(numImage*numOrient, sizeof(*trackMap));  
 for (img=0; img<numImage; img++)
     for (orient=0; orient<numOrient; orient++)
      {
       i = orient*numImage+img; 
       f = mxGetCell(prhs[c], i);
       datatype = mxGetClassID(f);
        if (datatype != mxINT32_CLASS)
            mexErrMsgTxt("warning !! int32 type required.");
       trackMap[i] = (int*)mxGetPr(f);
      }
 c++;
 /* about Gabor filters */
 halfFilterSize = ROUND(mxGetScalar(prhs[c++]));     
 Correlation = (double**)mxCalloc(numOrient*numOrient, sizeof(double*));   
 for (orient=0; orient<numOrient; orient++)
 {
       for (j=0; j<numOrient; j++)
        {
         f = mxGetCell(prhs[c], j*numOrient+orient);
         datatype = mxGetClassID(f);
        if (datatype != mxDOUBLE_CLASS)
            mexErrMsgTxt("warning !! double precision required.");
         Correlation[j*numOrient+orient] = (double*)mxGetPr(f); 
        }   
     }
 c++;  
 allSymbol = (double**)mxCalloc(numOrient, sizeof(double*));    
 for (orient=0; orient<numOrient; orient++)
     {  
       f = mxGetCell(prhs[c], orient);
       datatype = mxGetClassID(f);
        if (datatype != mxDOUBLE_CLASS)
            mexErrMsgTxt("warning !! double precision required.");
       allSymbol[orient] = (double*)mxGetPr(f);       
     }
 c++;
 /* about exponential model */ 
 numStoredPoint = ROUND(mxGetScalar(prhs[c++])); 
 storedlambda = (double*)mxGetPr(prhs[c++]);   
 storedExpectation = (double*)mxGetPr(prhs[c++]);   
 storedLogZ = (double*)mxGetPr(prhs[c++]);  
 /* learned parameters of active basis */
 selectedOrient = (double*)mxGetPr(prhs[c++]);                 
 selectedx = (double*)mxGetPr(prhs[c++]);         
 selectedy = (double*)mxGetPr(prhs[c++]);   
 selectedlambda = (double*)mxGetPr(prhs[c++]);       
 selectedLogZ = (double*)mxGetPr(prhs[c++]);  
 /* templates of images */
 commonTemplate = (float*)mxGetPr(prhs[c++]);     
 deformedTemplate = (float**)mxCalloc(numImage, sizeof(float*));   
 for (img=0; img<numImage; img++)
    {
        f = mxGetCell(prhs[c], img);
        datatype = mxGetClassID(f);
        if (datatype != mxSINGLE_CLASS)
            mexErrMsgTxt("warning !! single precision required.");
        deformedTemplate[img] = (float*)mxGetPr(f);  
     }
 c++;
 
 /* stored shifts */
 xShift = (int**)mxCalloc(numOrient, sizeof(*xShift));    
 for (orient=0; orient<numOrient; orient++)
     {
       f = mxGetCell(prhs[c], orient);
       numShift = mxGetNumberOfElements( f );
       datatype = mxGetClassID(f);
        if (datatype != mxINT32_CLASS)
            mexErrMsgTxt("warning !! int32 precision required.");
       xShift[orient] = (int*)mxGetPr(f);  
     }
 c++;
 
 yShift = (int**)mxCalloc(numOrient, sizeof(*yShift));    
 for (orient=0; orient<numOrient; orient++)
     {
       f = mxGetCell(prhs[c], orient);
       datatype = mxGetClassID(f);
        if (datatype != mxINT32_CLASS)
            mexErrMsgTxt("warning !! int32 precision required.");
       yShift[orient] = (int*)mxGetPr(f);       
     }
 c++;
 
 orientShifted = (int**)mxCalloc(numOrient, sizeof(*orientShifted));    
 for (orient=0; orient<numOrient; orient++)
     {  
       f = mxGetCell(prhs[c], orient);
       datatype = mxGetClassID(f);
        if (datatype != mxINT32_CLASS)
            mexErrMsgTxt("warning !! int32 precision required.");
       orientShifted[orient] = (int*)mxGetPr(f);       
     }
 c++;

 
 for (x=1; x<=sizex; x++)
      for (y=1; y<=sizey; y++)
          commonTemplate[px(x, y, sizex, sizey)] = 0.;   
 for (img=0; img<numImage; img++)
 {
     for (x=1; x<=sizex; x++)
         for (y=1; y<=sizey; y++)
         {
             deformedTemplate[img][px(x, y, sizex, sizey)] = 0.; 
         }
 }
 

 
 /* MAX1 maps are smaller than SUM1 maps */
 sizexSubsample = floor((double)sizex/subsample); 
 sizeySubsample = floor((double)sizey/subsample);
 /* run the shared sketch algorithm */

 Initialize();
 PursueElement();

 /* free matrices to avoid memory leak */  
// free_matrix((void**)pooledMax1map, numOrient, sizexSubsample*sizeySubsample);  // not necessary
}

