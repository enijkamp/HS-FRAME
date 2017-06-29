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

double *double_vector(int n)
{
    double *v; 
    v = (double*) mxCalloc (n, sizeof(double));
    return v; 
}

int *int_vector(int n)
{
    int *v; 
    v = (int*) mxCalloc (n, sizeof(int));
    return v; 
}

double **double_matrix(int m, int n)
{
    double **mat; 
    int i; 
    mat = (double**) mxCalloc(m, sizeof(double*)); 
    for (i=0; i<m; i++)
        mat[i] = double_vector(n); 
    return mat; 
}

int **int_matrix(int m, int n)
{
    int **mat; 
    int i; 
    mat = (int**) mxCalloc(m, sizeof(int*)); 
    for (i=0; i<m; i++)
        mat[i] = int_vector(n); 
    return mat; 
}



void free_matrix(void **a, int nrow, int ncol)
{
	int i;
	for(i = 0; i < nrow; i++)
		mxFree(a[i]);
	mxFree(a);
}


/* getting the index of matlab image, (x,y) location, (sx,sy) sizes */
int px(int x, int y, int bx, int by)    
{            
   return (x + (y-1)*bx - 1); 
 }


int n;                      /* number of images */
int N;                      /* number of filters */
int M;                      /* number of orientations */
int LAYER;                  /* number of scales */
int K;                      /* number of DoGs */
double **fIr, **fIi, **fI;                /* filtered images and max image */
double **Crr, **Cri, **Cir, **Cii;    /* inhibition coefficients */
double *half;                      /* halfsizes of filters */
int sx, sy;                 /* sizes of image and range of max pooling */ 
double **allsymbol, **allfilterr, **allfilteri;
double *syma, **sym, **Asym, **Asym0, **Asyma;
double *syma_top, **Asyma_top;
int L, ore; 
double TOTAL; 
double epsilon, Upperbound, SHUTUP; 
int **sinsh, **cossh;       /* store the shift  */
int **pI;                /* indicator */
int** MAX1mapML; /*relative perturbation in location */
int** MAX1mapMX; /*x position */
int** MAX1mapMY; /*y position */
int** MAX1mapMO; /*template index */
int** MAX1mapMD; /*relative perturbation in orientation */
double** MAX1mapMC; /* MAX1 map*/

double* filterSelected;
double* xSelected;
double* ySelected;
//double* numSketchSelected;


double* deformed_filterSelected;
double* deformed_xSelected;
double* deformed_ySelected;

int numTopFeatures;

int DoGLayer;  /* if DoG is not used, DoGLayer will be 0, which mean no template for illustrating DoGs */

/* store the shift values so that we do not need to repeat the sin and cos computation */ 
void storeshift()
{
    int ind, l;
    double theta; 
    sinsh = int_matrix(M, L+L+1); 
    cossh = int_matrix(M, L+L+1); 
    for (ind=0; ind<M; ind++)        
    {
        theta = PI*ind/M; 
        for (l=-L; l<=L; l++)
         {
            sinsh[ind][l+L] = floor(l*sin(theta)+.5); 
            cossh[ind][l+L] = floor(l*cos(theta)+.5); 
        }   
    }
}

/* local max in X or Y direction */
double shiftmax(int i, int ind, int x, int y, int *rm)  
{
   double m, en;
   int layer, ind1, x1, y1, l, ml, mx, my, de, d, d1, g, r, mo, md; 
   ml = 0;
   mx = 0;
   my = 0;
   mo = 0;
   md = 0;
   layer = floor(ind/M);  /* ind: index of filter; M: number of orientation */
   ind1 = ind - layer*M;  /* orientation */
   m = NEGMAX;  
   for (l=-L; l<=L; l++)
     {
      x1 = x + cossh[ind1][l+L]; 
      y1 = y + sinsh[ind1][l+L]; 
      if ((x1>=1)&&(x1<=sx)&&(y1>=1)&&(y1<=sy))
      {
        g = px(x1, y1, sx, sy);
        for (de=-ore; de<=ore; de++)
            {
              d = de+ind1;
              d1 = d;
              if (d<0)
                 d1 = d+M;
              if (d>=M)
                  d1 = d-M;
              r = (d1+layer*M)*n+i;
              en = fI[r][g];   
              if (en>m)
                {
                   m = en; 
                   ml = l; mx = x1; my = y1; mo = d1+layer*M; md = de; 
                }
             }
        }
     }
   rm[0] = ml; rm[1] = mx; rm[2] = my; rm[3] = mo; rm[4] = md; 
   return(m); 
   /*return(Upperbound*(2./(1.+exp(-2.*m/Upperbound))-1.)); */
}



/* local max in orientation */
double shiftmaxo(int i, int ind, int x, int y, int *rm)  
{
    double m, en;
    int x1, y1, lx, ly, ml, mx, my, mo, md, g; 
    ml = 0;
    mx = 0;
    my = 0;
    mo = 0;
    md = 0;
    m = NEGMAX;  
    for (lx=-L; lx<=L; lx++)
        for (ly=-L; ly<=L; ly++)
        {
            x1 = x + lx; 
            y1 = y + ly; 
            if ((x1>=1)&&(x1<=sx)&&(y1>=1)&&(y1<=sy))
            {
                g = px(x1, y1, sx, sy);
                en = fI[ind*n+i][g];   
                if (en>m)
                {
                    m = en; 
                    ml = lx; mx = x1; my = y1; mo = ind; md = ly; 
                }
            }
        }
   rm[0] = ml; rm[1] = mx; rm[2] = my; rm[3] = mo; rm[4] = md; 
   return(m); 
    /* return(Upperbound*(2./(1.+exp(-2.*m/Upperbound))-1.)); */
}


/* inhibitation */
void inhibit(int i, int mi, int mx, int my)   
{
   int x, y, ind, mh, h, g, r; 
   double *fr, *fi, *f, *fcrr, *fcri, *fcir, *fcii, ar, ai, mlayer, layer;
   /* mlayer = floor(mi/M); */
   mh = floor(half[mi]+.5); 
   ar = fIr[mi*n+i][px(mx, my, sx, sy)]; 
   ai = fIi[mi*n+i][px(mx, my, sx, sy)]; 
   for (ind=0; ind<N; ind++)   /* number of filters */
     {
     /*  layer = floor(ind/M); */
       h = floor(half[ind]+.5);  
       fr = fIr[ind*n+i]; 
       fi = fIi[ind*n+i]; 
       f = fI[ind*n+i];
       fcrr = Crr[mi+ind*N]; fcri = Cri[mi+ind*N];
       fcir = Cir[mi+ind*N]; fcii = Cii[mi+ind*N];
       for (x=MAX(1, mx-mh-h); x<=MIN(sx, mx+mh+h); x++)
         for (y=MAX(1, my-mh-h); y<=MIN(sy, my+mh+h); y++)
           {
              r = px(x, y, sx, sy);
              g = px(x-mx+mh+h+1, y-my+mh+h+1, 2*(mh+h)+1, 2*(mh+h)+1);
              fr[r] -= (fcrr[g]*ar+fcir[g]*ai);
              fi[r] -= (fcri[g]*ar+fcii[g]*ai);
              {
                   f[r] = fr[r]*fr[r] +fi[r]*fi[r];
                    
                   /*if (fcrr[g]*fcrr[g]+fcri[g]*fcri[g]+fcir[g]*fcir[g]+fcii[g]*fcii[g]>epsilon) f[r] = 0.; */
             
              }
           }
      }  
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

void draw_withoutShade(double *symo, int x0, int y0, int ind) /* because the allsymbol has shade property, we use this function to draw template without using different shades for different scales Gabor*/
{
  int x, y, h; 
   
  h = floor(half[ind]+.5);
  for (x=MAX(1, x0-h); x<=MIN(sx, x0+h); x++)
     for (y=MAX(1, y0-h); y<=MIN(sy, y0+h); y++)
       {
          if (allsymbol[ind][px(x-x0+h+1, y-y0+h+1, 2*h+1, 2*h+1)]!=0) {                 
                    
              symo[px(x, y, sx, sy)] = MAX(symo[px(x, y, sx, sy)], allsymbol[ind][px(x-x0+h+1, y-y0+h+1, 2*h+1, 2*h+1)] / allsymbol[ind][px(x-x0+h+1, y-y0+h+1, 2*h+1, 2*h+1)]);
          }
          else{
                   
              symo[px(x, y, sx, sy)] = MAX(symo[px(x, y, sx, sy)], allsymbol[ind][px(x-x0+h+1, y-y0+h+1, 2*h+1, 2*h+1)]);
          }
       }
}

void draw1(double *symo, int x0, int y0, int ind, int i)
{
  int x, y, h, r;
  double ar, ai; 

  ar = fIr[ind*n+i][px(x0, y0, sx, sy)]; 
  ai = fIi[ind*n+i][px(x0, y0, sx, sy)]; 
  h = floor(half[ind]+.5);
  for (x=MAX(1, x0-h); x<=MIN(sx, x0+h); x++)
     for (y=MAX(1, y0-h); y<=MIN(sy, y0+h); y++)
       {
          r = px(x-x0+h+1, y-y0+h+1, 2*h+1, 2*h+1); 
          symo[px(x, y, sx, sy)] += ar*allfilterr[ind][r]+ai*allfilteri[ind][r];
       }
}


/* new code here, calculate MAX1map */
void calculatemax()
{
    int k,i, j, x, y, h, here, rm[5], ind, sizeX, sizeY;
    double mm;
    MAX1mapML = (int**) mxCalloc(n*N, sizeof(int*));/*relative perturbation in location */
    MAX1mapMX = (int**) mxCalloc(n*N, sizeof(int*));/*x position */
    MAX1mapMY = (int**) mxCalloc(n*N, sizeof(int*));/*y position */
    MAX1mapMO = (int**) mxCalloc(n*N, sizeof(int*));/*template index */
    MAX1mapMD = (int**) mxCalloc(n*N, sizeof(int*));/*relative perturbation in orientation */
    MAX1mapMC = (double**) mxCalloc(n*N, sizeof(double*)); /* maximal MAX1 score */
    for (i=0; i<n; i++)  /* number of images */
    {
        for (ind=0; ind<N; ind++)  /* number of filters */
        {
            h = floor(half[ind]+.5);
            sizeX = sx-h-L-(h+1+L)+1;
            sizeY = sy-h-L-(h+1+L)+1;
            
             if ((sizeX<=0)||(sizeY<=0)){    
              mexPrintf("The valid region for matching pursuit is non-positive: sizeX=%d, sizeY=%d\n",sizeX, sizeY);
              mexErrMsgTxt("Error! Please use larger template size or smaller filter size!\n"); 
            }
            
            MAX1mapML[i*N+ind] = (int*) mxCalloc(sizeX*sizeY, sizeof(int));
            MAX1mapMX[i*N+ind] = (int*) mxCalloc(sizeX*sizeY, sizeof(int));
            MAX1mapMY[i*N+ind] = (int*) mxCalloc(sizeX*sizeY, sizeof(int));
            MAX1mapMO[i*N+ind] = (int*) mxCalloc(sizeX*sizeY, sizeof(int));
            MAX1mapMD[i*N+ind] = (int*) mxCalloc(sizeX*sizeY, sizeof(int));
            MAX1mapMC[i*N+ind] = (double*) mxCalloc(sizeX*sizeY, sizeof(double));
            for (x=h+1+L; x<=sx-h-L; x++)
            {
                for (y=h+1+L; y<=sy-h-L; y++)
                {
                    if (ind<N-K)
                        mm = shiftmax(i, ind, x, y, rm);   /* for Gabor filter */ 
                    else 
                        mm = shiftmaxo(i, ind, x, y, rm);  /* for DoG*/
                    j = x-h-1-L;
                    k = y-h-1-L; 
                    MAX1mapML[i*N+ind][j*sizeY+k] = rm[0];
                    MAX1mapMX[i*N+ind][j*sizeY+k] = rm[1];
                    MAX1mapMY[i*N+ind][j*sizeY+k] = rm[2];
                    MAX1mapMO[i*N+ind][j*sizeY+k] = rm[3];
                    MAX1mapMD[i*N+ind][j*sizeY+k] = rm[4];
                    MAX1mapMC[i*N+ind][j*sizeY+k] = mm;
                }
            }
        }
    }
}


/* new code here, update MAX1map */
/* you can update the parameter alpha to adjust inhabit neighbors */
void updateMAX1map(int mi, int orgX, int orgY)
{
    int halfX, halfY;
    float alpha = 0.7;
    int k,i, j, x, y, h, here, rm[5], ind, sizeX, sizeY, mh;
    double mm;
    mh = floor(half[mi]+.5); 
    for (i=0; i<n; i++)  /* all images */
    {
        for (ind=0; ind<N; ind++)  /* all filters */
        {
            h = floor(half[ind]+.5);
            sizeX = sx-h-L-(h+1+L)+1;
            sizeY = sy-h-L-(h+1+L)+1;
            
            /* neighbor sizes. */
            /* option 1: using a parameter alpha to control the size of the update region */
            /*
            halfX = floor(h*alpha); 
            halfY = floor(h*alpha); 
            */
            /* option 2: go over all the locations that may be affected */ 
            
            halfX = floor(h+mh+L); 
            halfY = floor(h+mh+L); 
            
            for (x = MAX(orgX-halfX,h+1+L); x<=MIN(sx-h-L,orgX+halfX);x++)
            {
                for (y = MAX(orgY-halfY,h+1+L); y<=MIN(sy-h-L,orgY+halfY);y++)
                {
                    if (ind<N-K)
                        mm = shiftmax(i, ind, x, y, rm);    
                    else 
                        mm = shiftmaxo(i, ind, x, y, rm);
                    j = x-h-1-L;
                    k = y-h-1-L; 
                    MAX1mapML[i*N+ind][j*sizeY+k] = rm[0];
                    MAX1mapMX[i*N+ind][j*sizeY+k] = rm[1];
                    MAX1mapMY[i*N+ind][j*sizeY+k] = rm[2];
                    MAX1mapMO[i*N+ind][j*sizeY+k] = rm[3];
                    MAX1mapMD[i*N+ind][j*sizeY+k] = rm[4];
                    MAX1mapMC[i*N+ind][j*sizeY+k] = mm;
                }
            }
        }
    }  
}



/* the match pursuit algorithm */
void Cmp()     
{
   int i, j, k, ind, x, y, mi, mx, my, mj, mk, t, rm[5], h, here, orie, layer, mlayer,sizeX, sizeY; 
   double m, s, mm; 
   pI = int_matrix(N, sx*sy);  
   for (ind=0; ind<N; ind++)
     for (here=0; here<sx*sy; here++)
         pI[ind][here] = 1; 
   storeshift();
   calculatemax();
   
   t = 0; 
   do
     {
      m = NEGMAX;  
      for (ind=0; ind<N; ind++)  /* for all filters*/
      { 
        h = floor(half[ind]+.5);
        sizeX = sx-h-L-(h+1+L)+1;
        sizeY = sy-h-L-(h+1+L)+1;
        for (x=h+1+L; x<=sx-h-L; x++)
           for (y=h+1+L; y<=sy-h-L; y++)
            {
             here = px(x, y, sx, sy); 
             if (pI[ind][here])
             {
               s = 0.; 
               for (i=0; i<n; i++)  /* for all images */
                 {  
                   /* original code here*/
                   /*
                     if (ind<N-K)
                       mm = shiftmax(i, ind, x, y, rm);    
                     else 
                       mm = shiftmaxo(i, ind, x, y, rm);
                   */
                   j = x-h-1-L;
                   k = y-h-1-L; 
                   mm = MAX1mapMC[i*N+ind][j*sizeY+k];
                   
                   s += mm; 
                  }   
                if (m<s)
                  {
                    m = s; 
                    mi = ind; 
                    mx = x; 
                    my = y; 
                    mj = j;   /* j is the actual index in MAX1mapMC, which can be computed by j=x-h-1-L   */
                    mk = k;   /* k is the actual index in MAX1mapMC, which can be computed by k=y-h-1-L;  */
                   }
                if (s/n<SHUTUP)
                   pI[ind][here] = 0; 
              }
             }
      }
    
      if (m/n>SHUTUP)     
      {
          
          t++; 
          mlayer = floor(mi/M);      /* M: number of orientation*/
          draw(syma, mx, my, mi);        /* one comprehensive template*/
          if (t<=numTopFeatures){
              draw_withoutShade(syma_top, mx, my, mi);  
          }           
          draw(sym[mlayer], mx, my, mi);   /* all scales of templates and one DoG template*/
     
            for (i=0; i<n; i++)
            { 
              /*original code here*/
              /*if (mi<N-K)
                   mm = shiftmax(i, mi, mx, my, rm);   
                else 
                   mm = shiftmaxo(i, mi, mx, my, rm); 
              */
         
              /* new code here*/
              h = floor(half[mi]+.5);
              sizeX = sx-h-L-(h+1+L)+1;  /* effective size of x */
              sizeY = sy-h-L-(h+1+L)+1;  /* effective size of y */
              rm[0] = MAX1mapML[i*N+mi][mj*sizeY+mk];
              rm[1] = MAX1mapMX[i*N+mi][mj*sizeY+mk];
              rm[2] = MAX1mapMY[i*N+mi][mj*sizeY+mk];
              rm[3] = MAX1mapMO[i*N+mi][mj*sizeY+mk];
              rm[4] = MAX1mapMD[i*N+mi][mj*sizeY+mk];
              
              deformed_filterSelected[i+n*(t-1)]=(double)(rm[3]+1);
              deformed_xSelected[i+n*(t-1)]=(double)rm[1];
              deformed_ySelected[i+n*(t-1)]=(double)rm[2];
                      
              draw1(Asym[i], rm[1], rm[2], rm[3], i);  /* draw constructed image */
              draw(Asym0[i*(LAYER+DoGLayer)+mlayer], rm[1], rm[2], rm[3]); 
              draw(Asyma[i], rm[1], rm[2], rm[3]);  /* draw sketches */
              if (t<=numTopFeatures){
                draw_withoutShade(Asyma_top[i], rm[1], rm[2], rm[3]);  
              } 
              inhibit(i, rm[3], rm[1], rm[2]);
     
            } 
     
          
            /*new code here */
           /*updateMAX1map(mi,mx,my); */
           updateMAX1map(mi,mx,my);     
                     
      }
        
      printf("No. %d:  (filter = %d   x = %d   y = %d) --> average response = %f\n", t, mi+1, mx, my, m/n); 
      filterSelected[t-1]=(double)(mi+1);
      xSelected[t-1]=(double)mx;
      ySelected[t-1]=(double)my;
   }
   while ((m/n>SHUTUP)&&(t<TOTAL));  /* stopping criterion */  
   
   //numSketchSelected[0]=(double)t;
}




void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
 int ind, i, j, c; 
 mxArray *f;  
 
 
 c=0;
 n = floor(mxGetScalar(prhs[c++])+.5);  /* number of images */
 N = floor(mxGetScalar(prhs[c++])+.5);  /* number of total filters*/
 M = floor(mxGetScalar(prhs[c++])+.5);  /* number of orientation */
 LAYER = floor((N+1.0)/M); 
 K = floor(mxGetScalar(prhs[c++])+.5);  /* number of DoG*/
 
 fIr = (double**) mxCalloc(n*N, sizeof(double*));    
 fIi = (double**) mxCalloc(n*N, sizeof(double*)); 
 fI =  (double**) mxCalloc(n*N, sizeof(double*));  
 for (i=0; i<n; i++)
   {
     for (ind=0; ind<N; ind++)
      {  
       f = mxGetCell(prhs[c], ind*n+i); 
       fIr[ind*n+i] = mxGetPr(f); 
       f = mxGetCell(prhs[c+1], ind*n+i); 
       fIi[ind*n+i] = mxGetPr(f);   
       f = mxGetCell(prhs[c+2], ind*n+i); 
       fI[ind*n+i] = mxGetPr(f);       
      }
    }
 c += 3;
 Crr = (double**) mxCalloc(N*N, sizeof(double*));    /* C: correlation/inhibition between filters */
 Cri = (double**) mxCalloc(N*N, sizeof(double*)); 
 Cir = (double**) mxCalloc(N*N, sizeof(double*)); 
 Cii = (double**) mxCalloc(N*N, sizeof(double*)); 
 for (ind=0; ind<N; ind++)
     {  
       for (j=0; j<N; j++)
        {
         f = mxGetCell(prhs[c], j*N+ind); 
         Crr[j*N+ind] = mxGetPr(f);         /* get correlation/inhibition */
         f = mxGetCell(prhs[c+1], j*N+ind); 
         Cri[j*N+ind] = mxGetPr(f);   
         f = mxGetCell(prhs[c+2], j*N+ind); 
         Cir[j*N+ind] = mxGetPr(f);   
         f = mxGetCell(prhs[c+3], j*N+ind); 
         Cii[j*N+ind] = mxGetPr(f);   
        }   
     }
 c += 4;     
 half = mxGetPr(prhs[c++]); 
 sx = floor(mxGetScalar(prhs[c++])+.5); 
 sy = floor(mxGetScalar(prhs[c++])+.5); 

 allfilterr = (double**) mxCalloc(N, sizeof(double*));  
 allfilteri = (double**) mxCalloc(N, sizeof(double*));  
 allsymbol = (double**) mxCalloc(N, sizeof(double*));    
 for (ind=0; ind<N; ind++)
     {  
       f = mxGetCell(prhs[c], ind); 
       allfilterr[ind] = mxGetPr(f);  
       f = mxGetCell(prhs[c+1], ind); 
       allfilteri[ind] = mxGetPr(f);        
       f = mxGetCell(prhs[c+2], ind); 
       allsymbol[ind] = mxGetPr(f);        
     }
 c += 3; 
 
 syma = mxGetPr(prhs[c++]); 
 
 syma_top = mxGetPr(prhs[c++]); 
 
 if (K>0){
     DoGLayer=1;
 }
 else{
     DoGLayer=0;
 }
     
 sym = (double**) mxCalloc(LAYER+DoGLayer, sizeof(double*)); 
 for (i=0; i<LAYER + DoGLayer; i++)
 {
        f = mxGetCell(prhs[c], i);
        sym[i] = mxGetPr(f);  
 }
 c++; 
 
 Asym = (double**) mxCalloc(n, sizeof(double*)); 
 for (i=0; i<n; i++)
 {
        f = mxGetCell(prhs[c], i);
        Asym[i] = mxGetPr(f);  
 }
 c++; 
 
 Asym0 = (double**) mxCalloc(n*(LAYER+DoGLayer), sizeof(double*)); 
 for (i=0; i<n*(LAYER+DoGLayer); i++)
    {
        f = mxGetCell(prhs[c], i);
        Asym0[i] = mxGetPr(f);  
     }
 c++;
 
 Asyma = (double**) mxCalloc(n, sizeof(double*)); 
 for (i=0; i<n; i++)
 {
        f = mxGetCell(prhs[c], i);
        Asyma[i] = mxGetPr(f);  
 }
 c++; 
 
 Asyma_top = (double**) mxCalloc(n, sizeof(double*)); 
 for (i=0; i<n; i++)
 {
        f = mxGetCell(prhs[c], i);
        Asyma_top[i] = mxGetPr(f);  
 }
 c++; 
 
 L = floor(mxGetScalar(prhs[c++])+.5);    /* location limit */
 ore = floor(mxGetScalar(prhs[c++])+.5);  /* orientation limit  */
 TOTAL = mxGetScalar(prhs[c++]);          /* number of sketches */
 epsilon = mxGetScalar(prhs[c++]); 
 Upperbound = mxGetScalar(prhs[c++]);  
 SHUTUP = mxGetScalar(prhs[c++]);  
 
 filterSelected = mxGetPr(prhs[c++]); 
 xSelected = mxGetPr(prhs[c++]); 
 ySelected = mxGetPr(prhs[c++]); 
 
 
 deformed_filterSelected = mxGetPr(prhs[c++]); 
 deformed_xSelected = mxGetPr(prhs[c++]); 
 deformed_ySelected = mxGetPr(prhs[c++]); 
         
 //numSketchSelected = mxGetPr(prhs[c++]); 
 numTopFeatures = (int)mxGetScalar(prhs[c++]);
  
 Cmp();
 
   
 
 /* release memory */
 mxFree(fIr);
 mxFree(fIi);
 mxFree(fI);
 mxFree(Crr);
 mxFree(Cri);
 mxFree(Cir);
 mxFree(Cii);
 mxFree(allfilterr);
 mxFree(allfilteri);
 mxFree(allsymbol);
 mxFree(sym);
 mxFree(Asym);
 mxFree(Asym0);
 mxFree(Asyma);
 mxFree(Asyma_top);
 for (i=0; i<n*N; i++)
 {
     mxFree(MAX1mapML[i]);
     mxFree(MAX1mapMX[i]);
     mxFree(MAX1mapMY[i]);
     mxFree(MAX1mapMO[i]);
     mxFree(MAX1mapMD[i]);
     mxFree(MAX1mapMC[i]);  
 }
mxFree(MAX1mapML);
mxFree(MAX1mapMX);
mxFree(MAX1mapMY);
mxFree(MAX1mapMO);
mxFree(MAX1mapMD);
mxFree(MAX1mapMC);
for (i=0; i<M; i++)
{
	mxFree(sinsh[i]);
	mxFree(cossh[i]);
}
mxFree(sinsh);
mxFree(cossh);
for (i=0; i<N; i++)
{
	mxFree(pI[i]);
}
mxFree(pI);
}
                 