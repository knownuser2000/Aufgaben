#include<stdlib.h>     
#include<stdio.h>     
#include<math.h>       
#include<unistd.h>
#include<string.h>
//#include "C:\Users\admin\Desktop\NumerikII_test\Library\MyNumerikBib.h" //oder 
#include "MyNumerikBib.h"

//externe Prototypen (aus MyNumerikBib oder anderen C-Datei)
extern void read_array(int ,double *); //andere C-Datei
 
int main(){
    
   //Deklarieren/Initialisieren von Variablen
    int dim  = 3; //Dimension von Array b
    double *b;

    b   =   calloc(   dim,sizeof(double));

    //Einlesen von Vektor aus Datei
    read_array(dim,b);

    //Ausgabe von Vektor in Konsole
   // print_vec(b, dim,"Vektor b");

    return 0;
}