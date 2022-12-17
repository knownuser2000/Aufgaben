#include "MyNumerikBib.h"


void fcn1(int n, double *x,double *f){
 f[0]= (-1.)*x[0];
 f[1]= x[0] - 10*x[1];
 f[2]= 10.*x[1];
}

void fcn2(int n, double *x,double *f){
 f[0]= (-1.)*x[0];
 f[1]= x[0] - pow(10.,6)*x[1];
 f[2]= pow(10.,6)*x[1];
}

//Erstelle Matrix mit Spalte 1=Zeit, Spalte 2,3,4:Konzentrationen
void Matrix_t_c(int n, int n_i, double *z, double *Z_,double *t_1){
    for(int i=0;i<n_i+1;i++){
        Z_[i*(n+1)]=t_1[i];         //Spalte 1=Zeit
        Z_[i*(n+1)+1]=z[i*n];       //Spalte 2: c_A
        Z_[i*(n+1)+2]=z[i*n+1];     //Spalte 3: c_B
        Z_[i*(n+1)+3]=z[i*n+2];     //Spalte 4: c_C
        }
}

//Erstelle Matrix mit Spalte 1=Zeit, Spalte 2,3, mit Fehlernormen
void Matrix_t_frob(int n, int n_i, double *error_Matrix, double *error_frob,double *error_y_j_frob,double *t_1){
    for(int i=0;i<n_i+1;i++){
        error_Matrix[i*(n)]=t_1[i];                 //Spalte 1=Zeit
        error_Matrix[i*(n)+1]=error_frob[i];        //Spalte 2: Frob von |y_t_2-y_t|/|y_t|
        error_Matrix[i*(n)+2]=error_y_j_frob[i];    //Spalte 3: Frob von |y_t|
        }
}


void writeFile(int n, int q,double *z_, double n_i, char *dateiname1){
    FILE *fd1;
    fd1 = fopen(dateiname1,"w");
    switch (q)
    {
    case 1:
        fprintf(fd1, "#t\t\t\t\t#c_A\t\t\t\t#c_B\t\t\t\t#c_C\n");
        for(int i=0;i<(n+1)*(n_i);i=i+n+1)
        {
        fprintf(fd1, "%e\t%e\t%e\t%e\n", z_[i],z_[i+1],z_[i+2],z_[i+3]); //z_[i]=c_A,z_[i+1]=c_B,z_[i+2]=c_C
        }
        break;
    
    default:
        fprintf(fd1, "#t\t\t\t\t#error_frob\t\t\t\t#cerror_y_j_frob\n");
        for(int i=0;i<n*(n_i+1);i=i+n+1)
        {
        fprintf(fd1, "%e\t%e\t%e\n", z_[i],z_[i+1],z_[i+2]);
        }
        break;
    }
    fclose(fd1);
}

int main(){
    int q=1;
    int x;
    int y;
    double tmax;
    double h;

    printf("Euler-Verfahren\n");
    printf("Geben sie für…\n\t\tEuler explizit:1\n\t\tEuler imlizit:2\n\t\tEuler semi implizite irgendeine Zahl ein\n");
    scanf("%d", &x);
    printf("Schrittgröße?\t\tFür a) h=0.05\t\t");
    scanf("%lf", &h);
    printf("Zeitbereich?\t\tFür a) tmax=20\t\t Für b) tmax=1e-04\n");
    scanf("%lf", &tmax);
    printf("Fall?\t\tFür Fall 1 (k_1=1, k_2=10)\t\t Gebe 1 ein\n");
    printf("\t\tFür Fall 2 (k_1=1, k_2=1e06)\t\t Gebe 2 ein\n");
    scanf("%d", &y);

    int i; 
    int n=3;
    double tmin=0;
    double t=tmin;
    double n_i=(tmax-tmin)/h;   //Schrittzahl
    double *t_1;

    t_1=calloc(n_i+1,sizeof(double));   //Zeitvektor, mit teilschritten
    t_1[0]=0;
    for(i=0;i<n_i;i++){
     t_1[i+1]=t_1[i]+h; 
    }

    double *y_i;
    y_i=calloc(n,sizeof(double));
    y_i[0]=1;
    y_i[1]=0;
    y_i[2]=0;

    double *z;
    z=calloc(n*n_i+n,sizeof(double));

    double *Z_;
    Z_=calloc((n+1)*n_i+n+1,sizeof(double));

    //Für semi implizit
    double *error_frob;
    error_frob=calloc(n_i+1,sizeof(double));

    double *error_y_j_frob;
    error_y_j_frob=calloc(n_i+1,sizeof(double));

    double *error_Matrix;
    error_Matrix=calloc(n*(n_i+1),sizeof(double));

    char *dateiname;        //Für die Konzentrationen
    char *dateiname2;       //Für die Fehlerrate

    switch (x)
    {
    case 1:         //Explizit
        switch (y)
        {
        case 1:         //k_1=1 k_2=10
            euler_explizit_2(n,h,t,tmax,y_i,z,fcn1);
            break;
        
        default:        //k_1=1 k_2=10^6
            euler_explizit_2(n,h,t,tmax,y_i,z,fcn2);
            break;
        }
        dateiname="Euler_explizit.txt";
        Matrix_t_c(n,n_i,z,Z_,t_1);
        writeFile(n,q,Z_,n_i,dateiname);
        system("gnuplot gnuplot_explizit.gpl -p");
        break;

    case 2:         //Implizit
        switch (y)
        {
        case 1:     //k_1=1 k_2=10
            euler_implizit(n,h,t,tmax,y_i,z,fcn1);
            break;
        
        default:    //k_1=1 k_2=10^6
            euler_implizit(n,h,t,tmax,y_i,z,fcn2);
            break;
        }
        dateiname="Euler_implizit.txt";
        Matrix_t_c(n,n_i,z,Z_,t_1);
        writeFile(n,q,Z_,n_i,dateiname);
        system("gnuplot gnuplot_implizit.gpl -p");
        break;
    
    default:        //Semi-implizit
        switch (y)
        {
        case 1:     //k_1=1 k_2=10
            euler_semi_implizit(n,h,t_1,tmax,y_i, z,error_frob,error_y_j_frob,fcn1);
            break;
        
        default:    //k_1=1 k_2=10^6
            euler_semi_implizit(n,h,t_1,tmax,y_i, z,error_frob,error_y_j_frob,fcn2);

            break;
        }
        dateiname="Euler_semi_implizit.txt";
        Matrix_t_c(n,n_i,z,Z_,t_1);
        writeFile(n,q,Z_,n_i,dateiname);
        system("gnuplot gnuplot_semi_implizit.gpl -p");
        dateiname2="Euler_semi_implizit_Fehler.txt";
        q=2;
        Matrix_t_frob(n,n_i,error_Matrix,error_frob,error_y_j_frob,t_1);
        writeFile(n,q,Z_,n_i,dateiname2);
        system("gnuplot gnuplot__Fehler.gpl -p");
        break;
    }

    free(t_1);
    free(z);
    free(Z_);
    free(error_frob);
    free(error_y_j_frob);
    free(error_Matrix);
}