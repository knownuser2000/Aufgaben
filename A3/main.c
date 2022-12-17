#include "MyNumerikBib.h"
#include "gnuplot.h"

double a[3];
    //Zentraldiff
void koef_ZD(double D, double delta_x, double v){
	a[0]=(D/delta_x -v/2)/delta_x;  //u_j_l+1
	a[1]=-2*D/pow(delta_x,2);         //u_j_l
	a[2]=(D/delta_x +v/2)/delta_x;  //u_j_l-1
}; 
    //Forward differentiaton

void koef_FD(double D, double delta_x, double v){
	a[0]=(D/delta_x -v)/delta_x;        //u_j_l+1
	a[1]= (-2*D/delta_x +v)/delta_x;    //u_j_l
	a[2]=(D/delta_x)/delta_x;           //u_j_l-1
}; 

//Backward differentiation
void koef_BD(double D, double delta_x, double v){
	a[0]=(D/delta_x)/delta_x;        //u_j_l+1
	a[1]=(-2*D/delta_x -v)/delta_x;     //u_j_l
	a[2]=(D/delta_x +v)/delta_x;        //u_j_l-1
}; 

void function(int n, double *x, double *f){
    f[0]=x[1]*a[0]+x[0]*a[1];
    for (int i = 1; i < n-1; i++)
    {
       f[i]= x[i+1]*a[0]+ x[i]*a[1]+x[i-1]*a[2];
    }
    f[n-1]=x[n-1]*a[1]+x[n-2]*a[2];
}

void function2(int n, double *x, double *f){
    f[0]=0;
    f[1]=x[1]*a[0]+x[0]*a[1];
    for (int i = 2; i < n-2; i++)
    {
       f[i]= x[i+1]*a[0]+ x[i]*a[1]+x[i-1]*a[2];
    }
    f[n-2]=x[n-2]*a[1]+x[n-3]*a[2];
    f[n-1]=0;
}

//n=13 By_dot=B*y'
void koeff_Matrix(int n, double *By_dot){
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; i < n; j++)
        {
            if (i==j-1)
            {
                By_dot[i*n+j]=a[2];     //u_j_(l-1) Koeff
            }
            else if (i==j)
            {
                By_dot[i*n+j]=a[1];     //u_j_(l) Koeff
            }
            else if (i==j+1)
            {
                By_dot[i*n+j]=a[2];     //u_j_(l-1) Koeff
            }
            else{
                By_dot[i*n+j]=0;
            }
        }
    }
};


//n= anzahl gl, l=zeit-zeile
void test(int n, int l, double *u, double *D, double delta_x, double n_x, double v){
    u[l*n]=u[l*n+n]*(D[0]/pow(delta_x,2)-v/(delta_x*2)) -2*u[l*n]/pow(delta_x,2) 
             + u[l*n-n]*(D[0]/pow(delta_x,2)+v/(delta_x*2));
    
    u[l*n+1]=u[l*n+n+1]*(D[1]/pow(delta_x,2)-v/(delta_x*2)) -2*u[l*n+1]/pow(delta_x,2) 
             + u[l*n-n+1]*(D[1]/pow(delta_x,2)+v/(delta_x*2));

    u[l*n+2]=u[l*n+n+2]*(D[2]/pow(delta_x,2)-v/(delta_x*2)) -2*u[l*n+2]/pow(delta_x,2) 
             + u[l*n-n+2]*(D[2]/pow(delta_x,2)+v/(delta_x*2));
};

int main(){
    int n_t,n_x;
    double *t, *x;
    double tmax,xmax;
    xmax=30,tmax=30;
    double delta_t,delta_x;
    delta_t=1, delta_x=1;    
    n_t=tmax/delta_t;
    n_x=xmax/delta_x;   
    int n=n_x-2; 
    double D=0.5;
    double v=1;
    koef_ZD( D, delta_x, v);

    t=calloc(n_t,sizeof(double));
    //Anfangsvektor
    t[0]=delta_t;
    for (int i = 1; i < n; i++)
    {
        t[i]=t[i-1]+delta_t;
    }
    
    x=calloc(n_x,sizeof(double));
    x[0]=20; //Anfangsvektor

    double *B;
    B=calloc((n_x-2)*(n_x-2),sizeof(double));
    
    for (int i = 0; i <(n_x-2); i++)
    {
        B[i*(n_x-2)+i]=1;
    }

    double *z; 
    z=calloc((n_x-2)*n_t,sizeof(double));
    int z_counter=0;
    int j_t;

    //void euler_SI(int n,int *j_t, int *z_counter,double delta_T, double *B, double *t,double tmax, double *y_i, double *z, void(fcn)(int, double*,double*)){
    euler_SI(n,&j_t,&z_counter,delta_t,B,t,tmax,x,z,function);
    //		euler_semi_implizit_Schritt(n,delta_T,B,t,y_t,y_t_2,z,fcn);

    //Spalte= konzentration
    //Zeile:Zeitpunkt
    double *a_a;
    a_a=calloc(n,sizeof(double));

    for (int i = 0; i < n; i++)
    {
        a_a[i]=z[i*n];
        /* code */
    }
    
    

    matrix_trans(n,n,z);
    struct gnuplot_arg plot={.grid = 1,.plotTitle="Aufgabe 3.1 D_1",.length=n,.inputFile="A2_3_1.txt"};
    gnuplot_multi_2D_in_3d(plot,t,t,n,z);
    char *test;
    gnuplot(plot,test,x,1,a_a);
}