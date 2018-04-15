/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * VITERBI MEX-C ALGORITHM
 ***************************
 * Usage of the module:
 * [u_hat u_win_hat]=viterbi(r,[p1 p2],w);
 * r is the decoder input array
 * u_hat is the output array
 * u_win_hat is the output array obtained with a window of w (optional,
 *   if w=0 no such output will be produced, see later)
 * [p1 p2] defines the puncturing scheme:
 * - if any of the two values are zero it means "no puncturing"
 * - otherwise, the puncturing scheme is considered to be such that,
 *   defining q1/q2 the simplified fraction with respect to (n*p1)/p2,
 *   every q1 symbols we transmit only q2 symbols. In this way, the
 *   overall code rate is p1/p2
 * w defines the length of the window used for the decoding (w=0 means
 *   no window, output only the word decoded with the traditional method)
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include "codeStructure.h"

/* prevWeight: weight of previous node
 * symbol: output symbol in the path from the previous node and the current one
 * r: received symbol
 * returns the weight of the path
 */
double getWeight(double prevWeight, int symbol[], double r[])
{
    double weight=0;
    int i;
    for(i=0;i<n;i++)
    {
        weight+=fabs((2*symbol[i]-1)-r[i]);
    }
    return prevWeight+weight;
}

/* a: array of double elements
 * aLength: length of the array
 * returns the index of the minimum value of the array
 */
int getMin(double a[], int aLength)
{
    int index=0;
    double min=a[0];
    int j;
    for(j=1;j<aLength;j++)
    {
        if(min>a[j])
        {
            min=a[j];
            index=j;
        }
    }
    return index;
}

/* serial to parallel conversion
 * a: input array
 * aLen: length of the array
 * b: output matrix of aLen/n rows and n columns
 */
void reshape(double a[], int aLen, double b[][n])
{
    int i;
    for(i=0;i<aLen;i++)
        b[i/n][i%n]=a[i];
}

/* return the gcd of two numbers */
int gcd(int a, int b)
{
    int remainder;
    while(a!=0)
    {
        remainder=b%a;
        b=a;
        a=remainder;
    }
    return b;
}

/* MAIN ALGORITHM */
void viterbi(double r0[], int r0_len, double u_hat[], int u_len,
        double punct[], int window, double u_hat_win[])
{
    /* initialize the structure of the code */
    initializeT();
    
    /* compensate the puncturing */
    if(punct[0]!=0 && punct[1]!=0)
    {
        int gcdVal=gcd(n*(int)punct[0],(int)punct[1]);
        int q1=((int)(n*punct[0]))/gcdVal;
        int q2=(int)punct[1]/gcdVal;
        
        int i=0;
        int j=0;
        int r1_len=r0_len*q1/q2;
        double r1[r1_len];
        while(i<r1_len)
        {
            if(i%q1>=q2)
            {
                r1[i]=0;
            }
            else
            {
                r1[i]=r0[j];
                j++;
            }
            i++;
        }
        r0=r1;
        r0_len=r1_len;
    }
    
    /* serial to parallel conversion */
    int r_len=r0_len/n;
    double r[r_len][n];
    reshape(r0, r0_len, r);
    
    /* initialize the trellis results container */
    struct StatesCollection
    {
        double weight[r_len+1][ns];
        int u_hat[r_len+1][ns];
        int prec_state[r_len+1][ns];
    };
    struct StatesCollection states;
    states.weight[0][0]=0;
    int i;
    for(i=1;i<ns;i++)
        states.weight[0][i]=INFINITY;
    
    /* -- main cycle -- */
    int a,b,c;
    for(a=1;a<=r_len;a++) /* cycle through trellis columns */
    {
        for(b=0;b<ns;b++) /* cycle through states */
        {
            /* get the candidate nodes */
            struct CodeStructure candidates=N(b);
            double candidatesWeight[candidates.aLength];
            
            /* cycle through candidates and calculate the respective weight */
            for(c=0;c<candidates.aLength;c++)
            {
                int m=candidates.prec_state[c];
                candidatesWeight[c]=getWeight(states.weight[a-1][m],
                    candidates.out_symbol[c], r[a-1]);
            }
            
            /* choose the best candidate */
            int minIndex=getMin(candidatesWeight, candidates.aLength);
            states.weight[a][b]=candidatesWeight[minIndex];
            states.u_hat[a][b]=candidates.in_symbol[minIndex];
            states.prec_state[a][b]=candidates.prec_state[minIndex];
        }
        
        /* normalize weight */
        double wMin=states.weight[a][getMin(states.weight[a],ns)];
        for(b=0;b<ns;b++)
        {
            states.weight[a][b]=states.weight[a][b]-wMin;
        }
        
        /* windowed decision */
        if(window>0 && a>window)
        {
            int state=getMin(states.weight[a],ns);
            for(b=a;b>a-window;b--)
                state=states.prec_state[b][state];
            
            u_hat_win[a-window-1]=states.u_hat[a-window][state];
        }
    }
    
    /* -- backtrack phase -- */
    int state=0;
    for(a=r_len;a>0;a--)
    {
        if(a<u_len)
        {
            u_hat[a-1]=states.u_hat[a][state];
            if(a>=r_len-window)
                u_hat_win[a-1]=states.u_hat[a][state];
        }
        state=states.prec_state[a][state];
    }
}


/* ENTRY POINT */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *r;
    double *u_hat;
    double *punct;
    double *u_hat_win=NULL;
    int window;
    int r_l,u_l;
    
    /* get input */
    r=mxGetPr(prhs[0]);
    r_l=mxGetN(prhs[0]);
    punct=mxGetPr(prhs[1]);
    if(mxGetN(prhs[1])!=2)
        mexErrMsgTxt("The length of the second paramenter must be 2!");
    window=(int)mxGetScalar(prhs[2]);
    
    /* prepare output */
    if(punct[0]!=0 && punct[1]!=0)
        u_l=r_l*punct[0]/punct[1]-deg;
    else
        u_l=r_l/n-deg;
    plhs[0]=mxCreateDoubleMatrix(1, u_l, mxREAL);
    u_hat=mxGetPr(plhs[0]);
    if(nlhs==2)
    {
        plhs[1]=mxCreateDoubleMatrix(1, u_l, mxREAL);
        u_hat_win=mxGetPr(plhs[1]);
    }
    
    /* call the main algorithm */
    viterbi(r, r_l, u_hat, u_l, punct, window, u_hat_win);
    
    return;
}
