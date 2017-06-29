/* mex-C: 
 * compute MAX2 maps for a single image (combine Yi Hong and Zhangzhang version)
 *
 * Usage:
 * 
 *     [tmpMAX2, tmpMAX2LocTrace, tmpMAX2TransformTrace,M2RowColShift] = ...
 *          mexc_sparseFRAME_MAX2( allTemplateAffinityMatrix(:), SUM2map(:,:,iRes), relativePartLocationRange, argmaxMethod, 1);
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
float **tmpM2Map;               /* temporary max2map */
int **M2LocationTrace;          /* ARGMAX map for row and column shifts */
int **tmpM2LocationTrace;
int **M2TemplateTrace;          /* ARGMAX map for template index. This index may be hybrid index of S2 template and its transform */
int **tmpM2TemplateTrace;
const float** S2Map;            /* SUM1 maps */
int height, width;              /* size of MAX1 maps */
float *rowShift, *colShift;           /* store the location/orientation shifts per scale */
int sizexsubsampleM2, sizeysubsampleM2; /* (down sampled) size of M1 map */
int *numLocationShift;              /* number of location shifts per scale */
int subsampleM2;                      /* sub sampling step length */
//float locationPerturbRadius;        /* allowed location perturbation, relative to the template size */
float locationPerturb;
const int **templateAffMat;         /* indices of similar templates (starting from 0), i.e. template affinity matrix in sparse format */
int *nNeighborTemplate;             /* number of neighbor templates per template */
int nLocationShift;                 /* number of location shifts for each template */
int nTemplate;                      /* number of S2 templates */
float stepsize;                     

/* store all the shifts or perturbations in local maximum pooling */
void StoreShift()  /* squared region */
{
    int j;
    float dc, dr;
    
    /* store all the possible shifts for all templates */
    stepsize = 1.0;
    j = 0;
    for( dc = -locationPerturb; dc <= locationPerturb; dc += stepsize )
    {
        for( dr = -locationPerturb; dr <= locationPerturb; dr += stepsize )
        {
            ++j;
        }
    }
    
    nLocationShift = j;
    
    rowShift = (float*)mxCalloc( nLocationShift, sizeof(*rowShift) );
    colShift = (float*)mxCalloc( nLocationShift, sizeof(*colShift) );
    
    /* store the location shifts in a lazy way */
    j = 0;
    for( dc = -locationPerturb; dc <= locationPerturb; dc += stepsize )
    {
        for( dr = -locationPerturb; dr <= locationPerturb; dr += stepsize )
        {
          
            rowShift[j] = dr;
            colShift[j] = dc;
            ++j;
        }
    }
}

void StoreShift2()   /* cross-shaped region */
{
    int j;
    float dc, dr;
    
    /* store all the possible shifts for all templates */    
       
    nLocationShift = locationPerturb*4+1;
    
    rowShift = (float*)mxCalloc( nLocationShift, sizeof(*rowShift) );
    colShift = (float*)mxCalloc( nLocationShift, sizeof(*colShift) );
    
    /* store the location shifts in a lazy way */
    j = 0;
    for( dc = -locationPerturb; dc < 0; dc ++ )
    {
        rowShift[j] = dc;
        colShift[j] = 0;
        ++j;
        
    }
    
    for( dc = 1; dc<=locationPerturb; dc ++ )
    {
        rowShift[j] = dc;
        colShift[j] = 0;
        ++j;
        
    }
    
    rowShift[j] = 0;
    colShift[j] = 0;
    ++j;
    
    for( dc = -locationPerturb; dc < 0; dc ++ )
    {
         rowShift[j] = 0;
         colShift[j] = dc;
         ++j;
    }
    
    for( dc = 1; dc<=locationPerturb; dc ++ )
    {
        rowShift[j] = 0;
        colShift[j] = dc;
        ++j;
        
    }
    
}


/*Compute MAX2map*/
void ComMax2Map()
{
    int orient, x, y, iS;
    int shift, x1, y1, orient1, mshift;
    int iT, jT, jjT;
    int rd;
    float currentMax;
    
    sizexsubsampleM2 = (int) floor((double)height / subsampleM2);//size of MAX2map
    sizeysubsampleM2 = (int) floor((double)width / subsampleM2);//size of MAX2map
    M2Map = (float**) mxCalloc( nTemplate, sizeof(*M2Map) );//MAX2map
    tmpM2Map = (float**) mxCalloc( nTemplate, sizeof(*tmpM2Map) );//temporary MAX2map
    M2LocationTrace = (int**)mxCalloc( nTemplate, sizeof(*M2LocationTrace) );//location pertubation
    tmpM2LocationTrace = (int**)mxCalloc( nTemplate, sizeof(*tmpM2LocationTrace) );//location pertubation
    M2TemplateTrace = (int**)mxCalloc( nTemplate, sizeof(*M2TemplateTrace) );//rotation pertubation
    tmpM2TemplateTrace = (int**)mxCalloc( nTemplate, sizeof(*tmpM2TemplateTrace) );//rotation pertubation
    
    //set initial MAX2map
    for( jT = 0; jT < nTemplate; ++jT )
    {
        M2Map[jT] = (float*)mxCalloc( sizeysubsampleM2 * sizexsubsampleM2, sizeof(**M2Map) );
        tmpM2Map[jT] = (float*)mxCalloc( sizeysubsampleM2 * sizexsubsampleM2, sizeof(**tmpM2Map) );
        M2LocationTrace[jT] = (int*)mxCalloc( sizeysubsampleM2 * sizexsubsampleM2, sizeof(**M2LocationTrace) );
        tmpM2LocationTrace[jT] = (int*)mxCalloc( sizeysubsampleM2 * sizexsubsampleM2, sizeof(**tmpM2LocationTrace) );
        M2TemplateTrace[jT] = (int*)mxCalloc( sizeysubsampleM2 * sizexsubsampleM2, sizeof(**M2TemplateTrace) );
        tmpM2TemplateTrace[jT] = (int*)mxCalloc( sizeysubsampleM2 * sizexsubsampleM2, sizeof(**tmpM2TemplateTrace) );
        for (y=0; y<sizeysubsampleM2; y++)
            for (x=0; x<sizexsubsampleM2; x++)
            {
                x1 = x * subsampleM2;
                y1 = y * subsampleM2;
                if ((x1>=0)&&(x1<height)&&(y1>=0)&&(y1<width))
                {
                    M2Map[jT][x+y*sizexsubsampleM2] = S2Map[jT][x1+y1*height];
                    M2LocationTrace[jT][x+y*sizexsubsampleM2] = (int) floor ( (double)nLocationShift / 2.0 ); /* index starts from 0*/
                    M2TemplateTrace[jT][x+y*sizexsubsampleM2] = jT;  /* index starts from 0*/
                    
                    tmpM2Map[jT][x+y*sizexsubsampleM2] = M2Map[jT][x+y*sizexsubsampleM2];  /* copy of M2Map*/
                    tmpM2LocationTrace[jT][x+y*sizexsubsampleM2] = M2LocationTrace[jT][x+y*sizexsubsampleM2];
                    tmpM2TemplateTrace[jT][x+y*sizexsubsampleM2] = M2TemplateTrace[jT][x+y*sizexsubsampleM2];
                }
            }
    }
    
    //maximize with respect to rotations
    for( iT = 0; iT < nTemplate; ++iT )
    {
        for( jjT = 0; jjT < nNeighborTemplate[iT]; ++jjT )
        {
            jT = templateAffMat[iT][jjT];
            for (y=0; y<sizeysubsampleM2; y++)
                for (x=0; x<sizexsubsampleM2; x++)
                {
                    if( tmpM2Map[jT][x+y*sizexsubsampleM2] > tmpM2Map[iT][x+y*sizexsubsampleM2] )
                    {
                        M2Map[iT][x+y*sizexsubsampleM2] = tmpM2Map[jT][x+y*sizexsubsampleM2];
                        M2TemplateTrace[iT][x+y*sizexsubsampleM2] = tmpM2TemplateTrace[jT][x+y*sizexsubsampleM2];
                        M2LocationTrace[iT][x+y*sizexsubsampleM2] = tmpM2LocationTrace[jT][x+y*sizexsubsampleM2];
                    }
                }
        }
    }
    
    //maximize with respect to location by recusive algorithm
   //for (rd=0; rd<locationPerturb; rd++)
   // {
        for( iT = 0; iT < nTemplate; ++iT )
        {
            
            
            for (y=0; y<sizeysubsampleM2; y++)
            {
                for (x=0; x<sizexsubsampleM2; x++)
                {
                    
                        currentMax = M2Map[iT][x+y*sizexsubsampleM2];
                        
                          for( shift=0; shift<nLocationShift; ++shift )  /* go over all the shifts */
                          {
                                  x1 = (int)floor( (x) * subsampleM2 + rowShift[shift]);
                                  y1 = (int)floor( (y) * subsampleM2 + colShift[shift]);
                                  
                                  if ((x1>=0)&&(x1<sizexsubsampleM2)&&(y1>=0)&&(y1<sizeysubsampleM2)) {
                                           
                                   if (currentMax < M2Map[iT][ x1 + y1 * sizexsubsampleM2])                        
                                   {
                                           currentMax = M2Map[iT][x1 + y1 * sizexsubsampleM2];
                                           tmpM2Map[iT][x+y*sizexsubsampleM2] = M2Map[iT][x1 + y1 * sizexsubsampleM2];
                                           tmpM2TemplateTrace[iT][x+y*sizexsubsampleM2] = M2TemplateTrace[iT][x1 + y1 * sizexsubsampleM2];
                                           tmpM2LocationTrace[iT][x+y*sizexsubsampleM2] = shift;
                       
                                   }                      
                                   
                                  }
                                                  
                          }
                        
                        /* if ( x>0 && x<sizexsubsampleM2-1 && y>0 && y<sizeysubsampleM2-1 )
                    {
                       if (currentMax < M2Map[iT][x-1+y*sizexsubsampleM2])
                        {
                            currentMax = M2Map[iT][x-1+y*sizexsubsampleM2];
                            tmpM2Map[iT][x+y*sizexsubsampleM2] = M2Map[iT][x-1+y*sizexsubsampleM2];
                            tmpM2TemplateTrace[iT][x+y*sizexsubsampleM2] = M2TemplateTrace[iT][x-1+y*sizexsubsampleM2];
                            tmpM2LocationTrace[iT][x+y*sizexsubsampleM2] = 0;
                        }
                        
                        if (currentMax < M2Map[iT][x+1+y*sizexsubsampleM2])
                        {
                            currentMax = M2Map[iT][x+1+y*sizexsubsampleM2];
                            tmpM2Map[iT][x+y*sizexsubsampleM2] = M2Map[iT][x+1+y*sizexsubsampleM2];
                            tmpM2TemplateTrace[iT][x+y*sizexsubsampleM2] = M2TemplateTrace[iT][x+1+y*sizexsubsampleM2];
                            tmpM2LocationTrace[iT][x+y*sizexsubsampleM2] = 1;
                        }
                        
                        if (currentMax < M2Map[iT][x+(y-1)*sizexsubsampleM2])
                        {
                            currentMax = M2Map[iT][x+(y-1)*sizexsubsampleM2];
                            tmpM2Map[iT][x+y*sizexsubsampleM2] = M2Map[iT][x+(y-1)*sizexsubsampleM2];
                            tmpM2TemplateTrace[iT][x+y*sizexsubsampleM2] = M2TemplateTrace[iT][x+(y-1)*sizexsubsampleM2];
                            tmpM2LocationTrace[iT][x+y*sizexsubsampleM2] = 3; 
                        }
                        
                        if (currentMax < M2Map[iT][x+(y+1)*sizexsubsampleM2])
                        {
                            currentMax = M2Map[iT][x+(y+1)*sizexsubsampleM2];
                            tmpM2Map[iT][x+y*sizexsubsampleM2] = M2Map[iT][x+(y+1)*sizexsubsampleM2];
                            tmpM2TemplateTrace[iT][x+y*sizexsubsampleM2] = M2TemplateTrace[iT][x+(y+1)*sizexsubsampleM2];
                            tmpM2LocationTrace[iT][x+y*sizexsubsampleM2] = 4; 
                        }
                        
                    }*/
                }
            }
            //copy results from above iterations
            for (y=0; y<sizeysubsampleM2; y++)
            {
                for (x=0; x<sizexsubsampleM2; x++)
                {
                    M2Map[iT][x+y*sizexsubsampleM2] = tmpM2Map[iT][x+y*sizexsubsampleM2];
                    M2TemplateTrace[iT][x+y*sizexsubsampleM2] = tmpM2TemplateTrace[iT][x+y*sizexsubsampleM2];
                    M2LocationTrace[iT][x+y*sizexsubsampleM2] = tmpM2LocationTrace[iT][x+y*sizexsubsampleM2];
                }
            }
        }
   // }
    
    //release memory
    mxFree( tmpM2Map );
    mxFree( tmpM2TemplateTrace );
    mxFree( tmpM2LocationTrace );
}



/* entry point */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
    int ind, i, x, y, o, dataDim, j, bytes_to_copy, nGaborFilter, toUseRegion;
    const mxArray *f;
    const mxArray *pATemplateAffMat;
    const mxArray *pAS2Map;
    mxArray *outPA;
    mwSize ndim;
    const mwSize* dims;
    mwSize dimsOutput[2];
    void* start_of_pr;
    mxClassID datatype;

	/*
	 * input variable 0: templateAffinityMatrix
	 */
    pATemplateAffMat = prhs[0];
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
	 * input variable 1: S2 maps
	 */
    pAS2Map = prhs[1];
    S2Map = (const float**)mxCalloc( nTemplate, sizeof(*S2Map) );   /* SUM1 maps */
    for (i=0; i<nTemplate; ++i)
    {
        f = mxGetCell(pAS2Map, i);
        datatype = mxGetClassID(f);
        if (datatype != mxSINGLE_CLASS)
            mexErrMsgTxt("warning !! single precision required.");
        S2Map[i] = (const float*)mxGetPr(f);    /* get the pointer to cell content */
        height = mxGetM(f);    /* overwriting is ok, since it is constant, or we would only input one image */
        width = mxGetN(f);
    }
    
    /*
     * input variable 2: location shift radius
     */
    locationPerturb = (float)mxGetScalar(prhs[2]);
    
    /*
     * input variable 3: method of local max for parts (1: squared region; 0: cross-shaped region)
     */
    toUseRegion = (int)mxGetScalar(prhs[3]);
 
    
    /*
     * input variable 4: sub sample step length
     */
    subsampleM2 = (int)mxGetScalar(prhs[4]);

    if(toUseRegion==1){
      StoreShift();
    }
    else{
      StoreShift2();
    }
    ComMax2Map();
    
    /* =============================================
     * Handle output variables.
     * ============================================= 
     */
    
    /*
     * output variable 0: M2 map
     */
    dimsOutput[0] = 1; dimsOutput[1] = nTemplate;
	plhs[0] = mxCreateCellArray( 2, dimsOutput );
    dimsOutput[0] = sizexsubsampleM2; dimsOutput[1] = sizeysubsampleM2;
    for( i = 0; i < nTemplate; ++i )
    {
        outPA = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
        /* populate the real part of the created array */
        start_of_pr = (float*)mxGetData(outPA);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(outPA);
        memcpy( start_of_pr, M2Map[i], bytes_to_copy );
        mxSetCell( plhs[0], i, outPA );
    }
    
    /*
     * output variable 1: M2 location trace
     */
    dimsOutput[0] = 1; dimsOutput[1] = nTemplate;
	plhs[1] = mxCreateCellArray( 2, dimsOutput );
    dimsOutput[0] = sizexsubsampleM2; dimsOutput[1] = sizeysubsampleM2;
    for( i = 0; i < nTemplate; ++i )
    {
        outPA = mxCreateNumericArray( 2, dimsOutput, mxINT32_CLASS, mxREAL );
        /* populate the real part of the created array */
        start_of_pr = (int*)mxGetData(outPA);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(outPA);
        memcpy( start_of_pr, M2LocationTrace[i], bytes_to_copy );
        mxSetCell( plhs[1], i, outPA );
    }
    
    /*
     * output variable 2: M2 other transformation trace
     */
    dimsOutput[0] = 1; dimsOutput[1] = nTemplate;
	plhs[2] = mxCreateCellArray( 2, dimsOutput );
    dimsOutput[0] = sizexsubsampleM2; dimsOutput[1] = sizeysubsampleM2;
    for( i = 0; i < nTemplate; ++i )
    {
        outPA = mxCreateNumericArray( 2, dimsOutput, mxINT32_CLASS, mxREAL );
        /* populate the real part of the created array */
        start_of_pr = (int*)mxGetData(outPA);
        bytes_to_copy = dimsOutput[0] * dimsOutput[1] * mxGetElementSize(outPA);
        memcpy( start_of_pr, M2TemplateTrace[i], bytes_to_copy );
        mxSetCell( plhs[2], i, outPA );
    }
    
    /*
     * output variable 3: stored template location shifts : row, col
     */
    dimsOutput[0] = nLocationShift; dimsOutput[1] = 2;
    outPA = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL );
    start_of_pr = (float*)mxGetData(outPA);
    bytes_to_copy = dimsOutput[0] * 1 * mxGetElementSize(outPA);
    memcpy( start_of_pr, rowShift, bytes_to_copy );
    start_of_pr = (float*)mxGetData(outPA) + nLocationShift;
    memcpy( start_of_pr, colShift, bytes_to_copy );
    plhs[3] = outPA;

}

