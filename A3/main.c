#include "MyNumerikBib.h"
#include "gnuplot.h"

double a[3];    //Koeff für D1
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

//Bestimme Peclet Zahl, und entscheide das Diskretisierungsverfahren anhand der Zahl
void Peclet(double D, double v, double delta_x ){
     double Pe;
     Pe=v*delta_x/D;  //Peclet Zahl
    printf("Pe:%lf\n", Pe);
        if (fabs(Pe)>1.9)
    {
        if (v>0)
        {
            koef_BD(D,delta_x,v);
        }else
        {
            koef_FD(D,delta_x,v);
        }
    }else
    {
        koef_ZD(D, delta_x, v);
    }
}

//RB: Direclet für x=0 und x[xmax]=0 werden ignoriert
void function(int n, double *x, double *f){
    f[0]=x[1]*a[0]+x[0]*a[1];
    for (int i = 1; i < n-1; i++)
    {
       f[i]= x[i+1]*a[0]+ x[i]*a[1]+x[i-1]*a[2];
    }
    f[n-1]=x[n-1]*a[1]+x[n-2]*a[2];
}

//RB: Direclet für x=0 und x[xmax]=0 beachtet [funktioniert noch nicht]
void function3(int n, double *x, double *f){
    f[0]=0;
    for (int i = 1; i < n-1; i++)
    {
       f[i]= x[i+1]*a[0]+ x[i]*a[1]+x[i-1]*a[2];
    }
    f[n-1]=0;
}

//Aufgabe 3. Ignorier
void function4(int n_x, int n_t, double v, double D, double *x, double *t, double *f){
double x0=0.5;
double x_n[n_x];
double *z;
z=calloc(n_t*n_x,sizeof(double));

for (int i = 0; i < n_x; i++)
{
    x_n[i]=-x0+x[i];
}

x_n[0]=x[0];            //für u[0] x0=0;
x_n[n_x-1]=x[n_x-1];    //für u[n_x-1] x0=0;

for (int i = 0; i < n_t; i++)
{
    for (int j = 0; j < n_x; j++)
    {
        z[i*n_x+j]=x_n[j]-v*t[i];
        z[i*n_x+j]= -pow(z[i*n_x+j],2)/(4*D*t[i]);
        z[i*n_x+j]= exp(z[i*n_x+j])/sqrt(4*PI*D*t[i]);
    }
}
}

//Erstellt die Konzentrationsvektoren für die t=10/50/75
void Vektor_aufgabe_1(int n, double *z, double *l1, double *l2, double *l3){
    for (int i = 0; i < n; i++)
    {
        l1[i]=z[n*(10)+i];
        l2[i]=z[n*(50)+i];
        l3[i]=z[n*(75)+i];
    }
}

void Aufgabe1(int xmax, int tmax, int delta_x, int delta_t, double *t, double *x, double D, double v, double *z){
    int n_t,n_x;        //Schrittzahl
    n_t=tmax/delta_t+1; //wegen startwert
    n_x=xmax/delta_x+1; 
    int n=n_x;
    double *B;  //Einheitsmatrix
    B=calloc((n)*(n),sizeof(double));   
    for (int i = 0; i <(n); i++)
    {
        B[i*(n)+i]=1;
    }
    //B[0]=0,B[n*n-1]=0;
    int z_counter=0;
    int j_t=0;
    Peclet(D,v,delta_x);
    euler_SI(n,&j_t,&z_counter,delta_t,B,t,tmax,x,z,function3);
}

void Aufgabe3(){
    double D[3];
    D[0]=0.01,D[1]=0.0008,D[2]=0.00005;
    double v=0.2;
    double tmax=1,xmax=1;
    double delta_x=0.04,delta_t=0.02;
    int n_t,n_x;
    n_t=tmax/delta_t+1;
    n_x=xmax/delta_x+1;
    int n=n_x;

    Peclet(D[0],v,delta_x);
    
}

int main(){
    double v=1,D[3];    //v=geschw. D=Diffussionskonstante
    D[0]=0.5, D[1]=0.05,D[2]=0.0005;

    double tmax=76,xmax=30;
    double delta_t=1,delta_x=1;
    int n_t=tmax/delta_t+1;
    int n_x=xmax/delta_x+1;        //Schrittzahl
    //n_x sollte xmax/delta_x +1 sein, weil x=0 und x=xmax
    int n=n_x;

    double *t, *x, *c;      //Zeit-; Orts-; Konzentrationsvektor
    t=calloc(n_t,sizeof(double));   //Zeitvektor
    x=calloc(n,sizeof(double));   
    c=calloc(n,sizeof(double));   //Konzentrationsvektor des ortes

    //x[0]=delta_x;
    c[1]=20; //Anfangsvektor
    //t[0]=delta_t;

    for (int i = 1; i < n; i++)
    {
        x[i]=x[i-1]+delta_x;
    }

    for (int i = 1; i < n_t; i++)
    {
        t[i]=t[i-1]+delta_t;
    }
    
    double *z1,*z2,*z3; //Konzentrationsmatrizen für die unterschiedlichen diffussionen
    z1=calloc((n)*n_t,sizeof(double));
    z2=calloc((n)*n_t,sizeof(double));
    z3=calloc((n)*n_t,sizeof(double));
    double *l1,*l2,*l3,*l4,*l5,*l6,*l7,*l8,*l9;
    
    //euler_SI(n,&j_t,&z_counter,delta_t,B,t,tmax,x,z1,function);
    //void Aufgabe1(int xmax, int tmax, int delta_x, int delta_t, double *t, double *x, double D, double v, double *z){
    Aufgabe1(xmax,tmax,delta_x,delta_t,t,c,D[0],v,z1);
    l1=calloc(n,sizeof(double));
    l2=calloc(n,sizeof(double));
    l3=calloc(n,sizeof(double));
    Vektor_aufgabe_1(n,z1,l1,l2,l3);

    Aufgabe1(xmax,tmax,delta_x,delta_t,t,c,D[1],v,z2);
    l4=calloc(n,sizeof(double));
    l5=calloc(n,sizeof(double));
    l6=calloc(n,sizeof(double));
    Vektor_aufgabe_1(n,z2,l4,l5,l6);

    Aufgabe1(xmax,tmax,delta_x,delta_t,t,c,D[2],v,z3);
    l7=calloc(n,sizeof(double));
    l8=calloc(n,sizeof(double));
    l9=calloc(n,sizeof(double));
    Vektor_aufgabe_1(n,z3,l7,l8,l9);
    
    char *legende[]={"D_1", "D_2", "D_3"};
    struct gnuplot_arg plot1={.grid = 1,.plotTitle="Aufgabe 3.2 t=10",.length=n,.inputFile="A2_3_t1.txt"};
    struct gnuplot_arg plot2={.grid = 1,.plotTitle="Aufgabe 3.2 t=50",.length=n,.inputFile="A2_3_t2.txt"};
    struct gnuplot_arg plot3={.grid = 1,.plotTitle="Aufgabe 3.2 t=75",.length=n,.inputFile="A2_3_t3.txt"};

    gnuplot(plot1,legende,x,3,l1,l4,l7);
    gnuplot(plot2,legende,x,3,l2,l5,l8);
    gnuplot(plot3,legende,x,3,l3,l6,l9);


    //euler_SI(n,&j_t,&z_counter,delta_t,B,t,tmax,x,z,function);
    //euler_SI(n_x,&j_t,&z_counter,delta_t,B,t,tmax,x,z,function);

    //		euler_semi_implizit_Schritt(n,delta_T,B,t,y_t,y_t_2,z,fcn);

    //matrix_trans(n,n,z);
    /* 
    struct gnuplot_arg plot={.grid = 1,.plotTitle="Aufgabe 3.1 D_1",.length=n,.inputFile="A2_3_1.txt"};
    gnuplot_multi_2D_in_3d(plot,t,t,n,z);
    */

    /* 
    char *test;
    gnuplot(plot,test,x,1,a_a);
    */
}