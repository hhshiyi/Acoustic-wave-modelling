/*
  Compute the 2D acoustic wave fields with perfect match layer.
  Author : Hu Hao
  Date   : 2013-10-04
  Revise : 2015-12-16
  Adress : HangZhou, China

  Notation : Using this code must have the agreement of author.

  Revise information : Refine the grid, memory and displaying.
*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "GL/glut.h"

#define  PI  3.141592654
#define  Xn  101
#define  Zn  101
#define  Tn  900
#define  Sx  50
#define  Sz  50
#define  ML  4
#define  PL  100

double** array_2(int x,int z);
void free_2(double** m,int x,int z);

double*** array_3(int x,int y,int z);
void free_3(double*** m,int x,int y,int z);

void Forward();
void DisplayWave(double **signal,int M,int N,double gain);

void vectorFlip(double *x,int N);

double **u0,**u1,**u2,**v,**coef_r,*ricker,*alpha,*cm;
double **receiver;

int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_DOUBLE| GLUT_RGBA|GLUT_DEPTH);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(500,500);
    glutCreateWindow("DrawWave");

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glClearColor(1,1,1,1);
    glClearDepth(1);

    glutDisplayFunc(&Forward);
    glutMainLoop();

    return 0;
}

void Forward()
{
    int i,j,k,m,n;
    double Delta_hx,Delta_hz,Delta_t,Is_Center,r,f,alpha1,alpha2,beta1,beta2,gain,delay;

    receiver=array_2(Tn,Xn+4*ML+2*PL);

    coef_r= array_2(Zn+4*ML+2*PL,Xn+4*ML+2*PL);  /* Malloc the memery for Zn rows and Xn columns. */
    u0 = array_2(Zn+4*ML+2*PL,Xn+4*ML+2*PL);     /* As above. */
    u1 = array_2(Zn+4*ML+2*PL,Xn+4*ML+2*PL);     /* As above. */
    u2 = array_2(Zn+4*ML+2*PL,Xn+4*ML+2*PL);     /* As above. */
    v  = array_2(Zn+4*ML+2*PL,Xn+4*ML+2*PL);     /* As above. */
    ricker = (double*)calloc(Tn,sizeof(double)); /* Discrete source function use the ricker wavelet. */
    alpha=(double*)calloc(PL,sizeof(double));    /* Store the reduce factors temporarily */
    cm=(double*)calloc(ML,sizeof(double));       /* Store the difference coefficient temporarily */

    /* Initialize the parameters */
    Delta_hx = 5.0;
    Delta_hz = 5.0;
    Delta_t = 0.001;
    Is_Center = 0;
    r = 3.7;
    f = 20;
    gain=0.05;
    delay=200.0;

    /* Initialize the velocity */
    for(i=0;i<Zn+4*ML+2*PL;i++)
        for(j=0;j<Xn+4*ML+2*PL;j++)
        {
            v[i][j]=1000.0;
        }

    /* Initializing the difference coefficient */
    cm[0]=1.196289,cm[1]=-0.079753,cm[2]=0.009570,cm[3]=-0.000698;

    /* Calculate the reduce factors */
    for(i=0;i<PL;i++)
    {
        alpha[i]=1.0-cos(PI*(PL-i)/(2.0*PL));
        alpha[i]*=100.0;
    }

    /* UP */
    for(i=ML;i<ML+PL;i++)
        for(j=ML;j<ML+PL+ML+Xn+ML+PL;j++)
            coef_r[i][j]=alpha[i-ML];

    /* LEFT */
    for(i=ML;i<ML+PL+ML+Zn+ML+PL;i++)
        for(j=ML;j<ML+PL;j++)
            coef_r[i][j]=alpha[j-ML];

    vectorFlip(alpha,PL);

    /* RIGHT */
    for(i=ML;i<ML+PL+ML+Zn+ML+PL;i++)
        for(j=ML+PL+ML+Xn+ML;j<ML+PL+ML+Xn+ML+PL;j++)
            coef_r[i][j]=alpha[j-ML-PL-ML-Xn-ML];

    /* DOWN */
    for(i=ML+PL+ML+Zn+ML;i<ML+PL+ML+Zn+ML+PL;i++)
        for(j=ML;j<ML+PL+ML+Xn+ML+PL;j++)
            coef_r[i][j]=alpha[i-ML-PL-ML-Zn-ML];

    /* LEFT-UP */
    for(i=ML;i<ML+PL;i++)
        for(j=ML;j<ML+PL;j++)
            coef_r[i][j]=coef_r[ML+PL][j]+coef_r[i][ML+PL];

    /* RIGHT-UP */
    for(i=ML;i<ML+PL;i++)
        for(j=ML+PL+ML+Xn+ML;j<ML+PL+ML+Xn+ML+PL;j++)
            coef_r[i][j]=coef_r[ML+PL][j]+coef_r[i][ML+PL+ML+Xn+ML-1];

    /* RIGHT-DOWN */
    for(i=ML+PL+ML+Zn+ML;i<ML+PL+ML+Zn+ML+PL;i++)
        for(j=ML+PL+ML+Xn+ML;j<ML+PL+ML+Xn+ML+PL;j++)
            coef_r[i][j]=coef_r[ML+PL+ML+Zn+ML-1][j]+coef_r[i][ML+PL+ML+Xn+ML-1];

    /* LEFT-DOWN */
    for(i=ML+PL+ML+Zn+ML;i<ML+PL+ML+Zn+ML+PL;i++)
        for(j=ML;j<ML+PL;j++)
            coef_r[i][j]=coef_r[ML+PL+ML+Zn+ML-1][j]+coef_r[i][ML+PL];


    /* Compute the 2D acoustic wave field. */
    for(k=1;k<Tn;k++)
    {
        /* Discrete ricker wavelet. */
        ricker[k] = exp((-4*PI*PI*f*f*(k-delay)*(k-delay)*Delta_t*Delta_t)/(r*r))*cos(2*PI*f*Delta_t*(k-delay));

        /* Display the field. */
        DisplayWave(u1,Zn+4*ML+2*PL,Xn+4*ML+2*PL,gain);

        for(i=ML;i<ML+PL+ML+Zn+ML+PL;i++)
            for(j=ML;j<ML+PL+ML+Xn+ML+PL;j++)
            {
                Is_Center=0.0;
                if(i==Sz+ML+PL+ML&&j==Sx+ML+PL+ML) Is_Center=1.0;

                beta1=0.0;
                beta2=0.0;
                for(m=1;m<=ML;m++)
                {
                    beta1+=cm[m-1]*(u1[i+m][j]+u1[i-m][j]-2.0*u1[i][j]); /* X axis diff */
                    beta2+=cm[m-1]*(u1[i][j+m]+u1[i][j-m]-2.0*u1[i][j]); /* Z axis diff */
                }
                /* With PML */
                alpha1=(2.0-coef_r[i][j]*coef_r[i][j]*Delta_t*Delta_t)*u1[i][j];
                alpha2=(coef_r[i][j]*Delta_t-1.0)*u0[i][j];
                u2[i][j]=1.0/(1.0+coef_r[i][j]*Delta_t)*(alpha1+alpha2+v[i][j]*v[i][j]*Delta_t*Delta_t*(beta1/(Delta_hz*Delta_hz)+beta2/(Delta_hx*Delta_hx)))+ricker[k]*Is_Center;

                /* store the wave field on free surface (i==0)
                   NOTATION : if you want to use the data to draw a image,
                   please output the receiver variation or u1 variation.*/
                receiver[k][j]=u2[ML+PL+ML][j];
            }

            for(m=ML;m<ML+PL+ML+Zn+ML+PL;m++)
                for(n=ML;n<ML+PL+ML+Xn+ML+PL;n++)
                {
                    u0[m][n]=u1[m][n];
                    u1[m][n]=u2[m][n];
                }
    }
    for(m=ML;m<ML+PL+ML+Zn+ML+PL;m++)
        for(n=ML;n<ML+PL+ML+Xn+ML+PL;n++)
        {
            u2[m][n]=u1[m][n];
            u1[m][n]=u0[m][n];
        }
}

/*
  Reverse order of vector
  <1,2,3...>  =>  <...3,2,1>
*/
void vectorFlip(double *x,int N)
{
    int i;
    double *t;
    t=(double*)calloc(N,sizeof(double));

    for(i=0;i<N;i++)
        t[i]=x[i];

    for(i=0;i<N;i++)
        x[i]=t[N-1-i];
}
/*
  Malloc memery for 2D array (z rows,x columns).
*/
double** array_2(int z,int x)
{
    double **m=NULL;
    int i;
    m=(double **)calloc(z,sizeof(double*));
    if(m==NULL) fprintf(stderr,"Can not malloc memory!\n");
    for(i=0;i<z;i++)
    {
        m[i]=(double*)calloc(x,sizeof(double));
        if(m[i]==NULL) fprintf(stderr,"Can not malloc memory!\n");
    }
    return m;
}
/*
  Malloc memery for 3D array.
*/
double ***array_3(int z,int y,int x)
{
    int i,j;
    double  ***m;
    m=(double***)malloc(z*sizeof(double**));
    if (!m)
    {
        fprintf(stderr,"Can not malloc memory for 3D array !\n");
        exit(1);
    }
    for (i=0;i<z;i++)
    {
        m[i]=(double **)malloc(y*sizeof(double*));
        if (!m[i])
        {
            fprintf(stderr,"Can not malloc memory for 3D array !\n");
            exit(1);
        };
        for (j=0;j<y;j++)
        {
            m[i][j]=(double *)malloc(x*sizeof(double));
            if (!m[i][j])
            {
                fprintf(stderr,"Can not malloc memory for 3D array !\n");
                exit(1);
            }
        }
    }
    return m;
}
/*
  Free the 2D array memery.
*/
void free_2(double** m,int z,int x)
{
    int i;
    for(i=0;i<z;i++)
        free(m[i]);
    free(m);
    m=NULL;
}
/*
  Free the 3D array memery .
*/
void free_3(double ***m,int z,int y,int x)
{
    int i,j;
    for (i=0;i<z;i++)
    {
        for (j=0;j<y;j++)
        {
            free(m[i][j]);
            m[i][j]=NULL;
        }
        free(m[i]);
        m[i]=NULL;
    }
    free(m);
    m=NULL;
}
/*
  Display the wave field. OpenGL is used in follow codes.
*/
void DisplayWave(double **signal,int N,int M,double gain)
{
    int i,j;
    double x,z,t;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //X-O-Z
    glColor3f(0.0f, 0.0f, 1.0f);
    for(i=ML+PL+ML;i<ML+PL+ML+Xn;i++)
    {
        x=(i-ML-PL-ML)*(2.0/Xn)-1.0;
        glBegin(GL_LINE_STRIP);
        for(j=ML+PL+ML;j<ML+PL+ML+Zn;j++)
        {
            t=signal[j][i]/1.0;
            z=1-(2.0f*(j-ML-PL-ML))/(1.0f*Zn);
            glVertex2f(x+t*gain,z);
        }
        glEnd();
    }

    glFlush();
    glutSwapBuffers();
}
