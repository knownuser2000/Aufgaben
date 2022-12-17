#include "MyNumerikBib.h"
#include "gnuplot.h"


int main(){
    double R=8.314;
    double T=293;
    double z=1;
    double F= 96485;
    double beta;
    double D=10E-9;
    double rho=50*10E-6;
    beta=D/rho;
    double i_0=10E-5;
    double i_2;
    double c_1;//mol/l=1000 mol/m^3;
    double *i, *c;
    i=calloc(10E3,sizeof(double));
    c=calloc(10E3,sizeof(double));



}