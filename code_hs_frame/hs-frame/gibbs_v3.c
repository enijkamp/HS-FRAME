#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"        /* the algorithm is connect to matlab */
#include "math.h"
#include "matrix.h"
#define ABS(x) ((x)>0? (x):(-(x)))
#define MAX(x, y) ((x)>(y)? (x):(y))
#define MIN(x, y) ((x)<(y)? (x):(y))

/* 
 c=gibbs_v3(template.selectedLambdas, response_IB, single(c_val_list), corrBB, j, single(I_squared_sum), response_IB(j), single(rand_num), nonzero_corrBB_List{j});    
*/

float ROUND(float d){
    
    return floor(d+0.5);
}

int px(int x, int y, int bx, int by)    
{            
   return (x + (y-1)*bx - 1); 
 }

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])                
{
  
 
  mxClassID datatype;  
  const float *selectedLambdas;
  int numSelectedFeature;
  float *response_IB;
  const float *c_val_list;
  int num_c_val_list;
  const float *corrBB;
  int heightcorrBB, widthcorrBB;
  float I_squared_sum, I_times_B_sum; 
  float *p, *cdf_p;
  float max_p, c_val, sum, log_sum_exp, c, ra;
  int i,j, iter, N, index;    
  mwSize dimsOutput[2];
  float* output_pointer;  
  float *nonzero_corrBB_List;
  int num_nonzero_corrBB_List;
  float sum2;
 
  /*
     * input variable 0: selectedLambdas
     */
    datatype = mxGetClassID(prhs[0]);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    selectedLambdas = (const float*)mxGetPr(prhs[0]); 
    numSelectedFeature = mxGetM(prhs[0])*mxGetN(prhs[0]);       
   
     /*
	 * input variable 1: response_IB
	 */
    datatype = mxGetClassID(prhs[1]);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    response_IB = (const float*)mxGetPr(prhs[1]); 
       
    
     /*
	 * input variable 2: c_val_list
	 */
    datatype = mxGetClassID(prhs[2]);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    c_val_list = (const float*)mxGetPr(prhs[2]); 
    num_c_val_list = mxGetM(prhs[2])*mxGetN(prhs[2]);        
        
    /*
	 * input variable 3: corrBB
	 */
    datatype = mxGetClassID(prhs[3]);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    corrBB = (const float*)mxGetPr(prhs[3]); 
    heightcorrBB = mxGetM(prhs[3]);    
    widthcorrBB = mxGetN(prhs[3]);
    
     /*
     * input variable 4: j
     */
    j = mxGetScalar(prhs[4])-1; 
     /*
	 * input variable 5, 6 and 7: 
	 */
    I_squared_sum = (float)mxGetScalar(prhs[5]);  /* sum_x{I^2} */
    /* B_squared_sum = mxGetScalar(prhs[6]);*/  /* sum_x{B^2}  */
    I_times_B_sum = (float)mxGetScalar(prhs[6]);  /* sum_x{IB} */
    
    ra = (float)mxGetScalar(prhs[7]);  /* random number */   
    
    
    datatype = mxGetClassID(prhs[8]);
    if (datatype != mxSINGLE_CLASS)
        mexErrMsgTxt("warning !! single precision required.");
    nonzero_corrBB_List = (const float*)mxGetPr(prhs[8]); 
    num_nonzero_corrBB_List = mxGetM(prhs[8])*mxGetN(prhs[8]);  
     
     
    /*
     *
     *start to compute
     *
     *
     */
    
    /* p=zeros(size(c_val_list));  */
    p = (float*)mxCalloc( num_c_val_list, sizeof(float) );
    for(iter=0; iter<num_c_val_list; iter++ ){
        p[iter]=0;        
    }   
    
    
    max_p=-10000000000;     
    for(iter=0; iter<num_c_val_list; iter++ ){
        c_val=c_val_list[iter];
        sum=0;        
        
        for(i=0; i<num_nonzero_corrBB_List; i++ ){     /*  go over all the elements in neighborhood system of basis j */
            index=(int)nonzero_corrBB_List[i]-1;
            sum += selectedLambdas[index] * ABS( response_IB[index] + c_val * corrBB[j+index*heightcorrBB]  );   /* sum_i{ lambda_i * |<I, B_i> + c <B_j, B_i>|}  i in the neighborhood system  */ 
           
        }
 
        /* add reference probability q(I)*/ 
        sum = sum -  0.5* ( I_squared_sum + c_val*c_val + 2 * c_val * I_times_B_sum);     /* sum_i{ lambda_i * |<I, B_i> + c <B_j, B_i>|} - 0.5 sum_x {I+cB_j} */     
        p[iter]=sum;    
        
        if (sum>max_p){ max_p=sum; }  /* store the maximum of p for safe softmax*/
        
    }
       
    /* p = exp (p - ( mx+log(sum(exp(p-mx))) ));  */
    /* compute log_sum_exp(p) = max + log_sum_exp(p-max)  */
    log_sum_exp=0;
    for(iter=0; iter<num_c_val_list; iter++ ){
        log_sum_exp += exp(p[iter] - max_p);        
    } 
    log_sum_exp = log(log_sum_exp) + max_p;
   
    
    /*  softmax and cdf */
    for(iter=0; iter<num_c_val_list; iter++ ){
        p[iter] = exp(p[iter] - log_sum_exp);        
    } 
  
    /* cdf  */
    /* cdf_p=cumsum(p); */  
    cdf_p = (float*)mxCalloc( num_c_val_list, sizeof(float) );
    cdf_p[0]=p[0];
    for(iter=1; iter<num_c_val_list; iter++ ){
        cdf_p[iter]=cdf_p[iter-1] + p[iter];        
    }   
    
    /* inverted cdf */
    /*rnum=rand();
	ra = rnum/RAND_MAX;
    */
    /* srand( time(0));*/    /*set seed */
    /*ra = ((float)rand()+1.0)/((float)RAND_MAX+2.0);*/    
	iter=0;
	while( (cdf_p[iter]<ra) && (iter<num_c_val_list-1) ){  /* we assume cdf_p is sorted in ascending order */
		iter++;
	}
	c=c_val_list[iter];
        
    
    /*printf("ra=%f, p[0]=%f, iter=%d, c=%f\n", ra, p[0], iter, c);*/ 
     /* =============================================
     * Handle output variables.
     * ============================================= 
     */
    
    dimsOutput[0] =  1; dimsOutput[1] = 1;
    plhs[0] = mxCreateNumericArray( 2, dimsOutput, mxSINGLE_CLASS, mxREAL);
    output_pointer=mxGetPr(plhs[0]);
    *output_pointer=c;
        
  
}