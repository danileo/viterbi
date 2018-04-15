#define N(i) T[i]
/* n (assuming k=1 => n=1/R): */
#define n 2
/* nu: */
#define deg 2
/* number of states (2^deg): */
#define ns 4

struct CodeStructure
{
    int prec_state[2]; /* preceding states */
    int in_symbol[2]; /* respective input symbols */
    int out_symbol[2][n]; /* output symbols (one symbol in one row) */
    int aLength; /* number of incoming edges */
};

struct CodeStructure T[4];

void initializeT()
{
    /* [1 (1+D^2)/(1+D+D^2)] */
    T[0].prec_state[0]=0;
    T[0].prec_state[1]=1;
    T[0].in_symbol[0]=0;
    T[0].in_symbol[1]=1;
    T[0].out_symbol[0][0]=0;
    T[0].out_symbol[0][1]=0;
    T[0].out_symbol[1][0]=1;
    T[0].out_symbol[1][1]=1;
    T[0].aLength=2;
    
    T[1].prec_state[0]=2;
    T[1].prec_state[1]=3;
    T[1].in_symbol[0]=1;
    T[1].in_symbol[1]=0;
    T[1].out_symbol[0][0]=1;
    T[1].out_symbol[0][1]=0;
    T[1].out_symbol[1][0]=0;
    T[1].out_symbol[1][1]=1;
    T[1].aLength=2;
    
    T[2].prec_state[0]=0;
    T[2].prec_state[1]=1;
    T[2].in_symbol[0]=1;
    T[2].in_symbol[1]=0;
    T[2].out_symbol[0][0]=1;
    T[2].out_symbol[0][1]=1;
    T[2].out_symbol[1][0]=0;
    T[2].out_symbol[1][1]=0;
    T[2].aLength=2;
    
    T[3].prec_state[0]=2;
    T[3].prec_state[1]=3;
    T[3].in_symbol[0]=0;
    T[3].in_symbol[1]=1;
    T[3].out_symbol[0][0]=0;
    T[3].out_symbol[0][1]=1;
    T[3].out_symbol[1][0]=1;
    T[3].out_symbol[1][1]=0;
    T[3].aLength=2;
}
