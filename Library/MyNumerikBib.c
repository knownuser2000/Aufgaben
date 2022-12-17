//Musterlösung MyNumerikBib "Numerische Methoden I - ICVT Uni Stuttgart"
//Autoren: Nikos Alexiadis, Jan Brunner, Giuliano Lombardo

#include "MyNumerikBib.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////Euler-Verfahren/////////////////////////////////////

void euler_explizit_2(int n, double h, double t, double tmax, double *y_i, double *z,void(fcn)(int, double*,double*)){
    int i;
    int k=n;
    double *f;
    f=calloc(n,sizeof(double));
    copy_array(z,y_i,n);
    do{
        fcn(n,y_i,f);
        for(i=0;i<n;i++){
        y_i[i]=y_i[i]+h*f[i];
        z[k]=y_i[i];
        k=k+1;
        }
		t=t+h;      //Erhöhung von t
    }while(t<tmax);
    free(f);
}

void euler_modified(int n, double h, double t, double tmax, double *y_i, double *z,void(fcn)(int, double*,double*)){
    //basic equation. y^(n+1)=y^n+0.5*h*[f(x_n,y_n)+f(x_n+h,y_n+h*f(x_n,y_n))]
    //y^(n+1)=y^n+0.5*h*[f(x_n,y_n)+f(x_n+h,y_n_dach)]
    //y_n_dach=y_i[i]+h*f[i]        //explizites verfahren. 
    int i;
    int k=n;
    double *f;
    f=calloc(n,sizeof(double));
    double *f_plus;
    f_plus=calloc(n,sizeof(double));
    double y_dach[3];
    double y_i_plus[3];
    copy_array(z,y_i,n);
    do{
        fcn(n,y_i,f);
        for(i=0;i<n;i++){
        y_dach[i]=y_i[i]+h*f[i];    //euler explizit
        }
        fcn(n,y_dach,f_plus);
        for(i=0;i<n;i++){
        y_i_plus[i]=y_i[i]+(h/2)*(f[i]+f_plus[i]);
        z[k]=y_i_plus[i];
        }
        t=t+h;      //Erhöhung von t
        k=k+1;
    }while(t<tmax);
    free(f);
    free(f_plus);
}

void euler_explizit(int n,int k, double h, double t, double tmax, double *y_i, double *z,void(fcn)(int, double*,double*)){
    int i;
    double *f;
    f=calloc(n,sizeof(double));
    //copy_array(double *arr1, double *arr2,int size); //arr1[i] = arr2[i] 
    //copy_array(z,y_i,n);
    for(i=0;i<n;i++){
        z[i]=y_i[i];
    }
    
    do{
        t=t+h;      //Erhöhung von t
        fcn(n,y_i,f);
        for(i=0;i<n;i++){
        y_i[i]=y_i[i]+h*f[i];
        z[k]=y_i[i];
        k=k+1;
        }
    }while(t<tmax);
//     for(i=0;i<k;i++){
//         printf("k[%i]=%lf\n",i,z[i]);
//     }
    free(f);
    }

void euler_implizit(int n, double h,double t, double tmax, double *y_i, double *z, void(fcn)(int, double*,double*)){
    //n=Dimension des Vektors, t=Startbeginn,aktueller Ort beim Integral, tmax=Integrationsende, 
    //y_i ist der Vektor mit y(0) anfangs, aber ist eig. der Vektor für y^(i+1).
    //y_0 ist eig. y^(i).
    int i;   
    int k=0;
    double h_jac=pow(10.,-8.);          //Schrittweite jacobi Matrix
    double tol= pow(10.,-6);            //Toleranz
    double *DF;                         //Jacobi Matrix
    double *L;                          //Left triangular Matrix LR-Zerlegung
    double *R;                          //Right triangular matrix
    DF = calloc(n*n,sizeof(double));
    L =  calloc(n*n,sizeof(double));
    R =  calloc(n*n,sizeof(double));
    double y_0[n];
    double F[n];
    double y[n];
    double delta_y[n];
    double *f;
    f=calloc(n,sizeof(double));
    double frob_delta_y;
    double frob_y_i;
    int iteration;
    int max_iteration=1000;
    double I_3[n*n];                    //soll einheitsmatrix sein
    
    for(i=0;i<n*n;i++){         //Nullmatrix I_3
        I_3[i]=0;
    }
    for(i=0;i<n;i++){
        I_3[i*n+i]=1.;     
    }         //Einheitsmatrix I_3

    do{
        t=t+h;                    //Änderung der Integralstelle
        fcn(n,y_i,f);                         //f(y_(i+1))
        for(i=0;i<n;i++){
                y_0[i]=y_i[i];                  //y_von i
                y_i[i]=y_i[i]+h*f[i];           //y_von_i+1
                z[k]=y_0[i]; 
                k=k+1;
                iteration=0;
        }                         //In der oberen do loop wird y^(i+1) über das explizite Verfahren gewonnen.
        do
        {
            fcn(n,y_i,f);
            fjac_cd(n,h_jac,f,DF,fcn);          //DF(y[i+1])
            for(i=0;i<n;i++){
                F[i]=y_0[i]+h*f[i]-y_i[i];      //
            }
            for(i=0;i<n*n;i++){
                DF[i]=h*DF[i]-I_3[i];
            }
            //LR_seperation(DF,F,L,R,n);
            //void LR_seperation(double *A, double *b, double *L, double *R,int n){
            //void LR_separation(int n, double *A,int lda, double *b, double *L,double *R){
            LR_separation(n,DF,n,F,L,R);
            array_minus(F,n);
            /*for(i=0;i<n;i++){
                F[i]=-F[i];    //xk als erster 
            }*/
            forward_substitution(L,F,y,n);
            back_substitution(R,y,delta_y,n);
            for(i=0;i<n;i++){
                y_i[i]=y_i[i]+delta_y[i];
            }
            //Frobeniusnorm(n,1,delta_y,&frob_delta_y);
            frob_delta_y=nmrc_mnorm(n,1,delta_y,n,NmrcNormFrob);
            //Frobeniusnorm(n,1,y_i,&frob_y_i);
            frob_y_i=nmrc_mnorm(n,1,y_i,n,NmrcNormFrob);
            iteration++;
            //innerer Loop 
        } while (frob_delta_y>tol*frob_y_i && iteration<max_iteration);
        
    }while(t<tmax);
    k=k-1;
    free(DF);
    free(L);
    free(R);
}

//norm = nmrc_mnorm(n, m, a[0], lda, NmrcNormRS);
//norm = nmrc_mnorm(n, m, a[0], lda, NmrcNormCS);
//norm = nmrc_mnorm(n, m, a[0], lda, NmrcNormFrob);
//double nmrc_mnorm(int n, int m, const double *a, int lda, enum NmrcNorm type)
//void Frobeniusnorm(int n, int m, double *A, double *f){
//nmrc_mnorm(n,1,delta_y,n,NmrcNormFrob)

void euler_Schrittweitensteuerung(int n,double h, int *s, double *t, double tmax, double *y_i, double *z, double P,double rho, double tol, double h_min, double n_w, void(fcn)(int, double*,double*)){    
    double a, b, Maxnorm;
    double EST;
    int i;
    double *f;
    double *y_s, *y_h, *y_d;
    y_s=calloc(n,sizeof(double));
    y_h=calloc(n,sizeof(double));
    y_d=calloc(n,sizeof(double));
    f=calloc(n,sizeof(double));
    double *delta_y;
    delta_y=calloc(n,sizeof(double));
    printf("Jonni ist cool\n");
    int j=0;
    for(i=0;i<n;i++){
        z[i]=y_i[i];    //Startwert draufschreiben;
    }
    
    do{
        h=h+h;
        do{
            h=h/2.;
            fcn(n,y_i,f);
            for(i=0;i<n;i++){
                y_s[i]=y_i[i]+h*f[i];       //y_s y_ schräg
                y_h[i]=y_i[i]+(h/2.)*f[i];       //y_h y halbe
            }
            fcn(n,y_h,f);
            for(i=0;i<n;i++){
                y_d[i]=y_h[i]+(h/2.)*f[i];         //y_dach
            }
            for(i=0;i<n;i++){
                delta_y[i]=y_s[i]-y_d[i]; 
            }
//             printf("Jonni war cool\n");
            //Maximumnorm(n,1,delta_y,&Maxnorm); nmrc_mnorm(n,1,y_i,n,NmrcNormFrob);
            //en = nmrc_condA(n,  A, lda, NmrcNormRS);    // für Maximumnorm

            Maxnorm=nmrc_condA(n,y_i,n,NmrcNormRS);
            EST=Maxnorm/(1-pow(2.,-P));
        }while((EST/h)>tol && h>h_min);
        t[j+1]=t[j]+h;
        j=j+1;
        a=n_w*h;                 //Wachstumsschranke*Schrittweite
        b=rho*pow((tol/EST)*pow(h,P+1),1/P);
        if(b<a){
            a=b;
        }
        if(h_min<a){
            h=a;
        }else{
            h=h_min;
        }
        if(t[j]+h>=tmax){
            h=tmax-t[j];
        }
        for(i=0;i<n;i++){
            y_i[i]=y_d[i];          //y_dach sei neuer y-Wert
            z[j*n+i]=y_i[i];            //Vektor mit allen y Werten
        }
    }while(t[j]<tmax);
     *s=j;
     printf("s=%i\n",j);
free(y_d);
free(y_s);
free(y_h);
free(f);
free(delta_y);
}

void euler_semi_implizit(int n,double delta_T, double *t, double tmax,double *y_i, double *z, double *error_frob,double *error_y_j_frob, void(fcn)(int, double*,double*)){

int i;
int k_steps=2;
int j_t=0;		//t[j_t+1]=t[j_t]+delta_T
int z_counter=0;
int error_counter=0;

double delta_frob;			//Frob von y_t_2-y_t
//double *error_y_j_frob,
double delta_frob_V[n];
double y_j_frob_V[n];	//Abspeicherung von y_t bevor frob

double *y_h;
y_h=calloc(n,sizeof(double));
double *y_t;
y_t=calloc(n,sizeof(double));
double *y_t_2;
y_t_2=calloc(n,sizeof(double));
double *delta_y;
delta_y= calloc(n,sizeof(double));

double h;
//delta_T;	//Schrittweite eines Basisschrittes
double *h_Werte;
h_Werte=calloc(k_steps,sizeof(double));
double h_jac=pow(10.,-8.);          //Schrittweite jacobi Matrix
double tol= pow(10.,-6);            //Toleranz
double B[n*n];                    //soll einheitsmatrix sein
double J_h[n*n]; 

double *DF;                         //Jacobi Matrix
DF = calloc(n*n,sizeof(double));

double *P;
P = calloc(n*k_steps,sizeof(double));

double *f;					//funktionswerte
double f_h[n];				//f_h[i]=h*f[i];
f=calloc(n,sizeof(double));

for(i=0;i<n*n;i++){         //Nullmatrix B
        B[i]=0;
		J_h[i]=0;
    }
for(i=0;i<n;i++){
    B[i*n+i]=1.;     //Einheitsmatrix B
	z[i]=y_i[i];	//Gesamtvektor mit allen y_werten
}         
copy_array(y_t,y_i,n);	//alles auf y_t

do
{
    fjac_cd(n,h_jac,y_t,DF,fcn);          //f'(y^t)=DF
	for (int k = 1; k <= k_steps; k++)
	{
		i=0;
		copy_array(y_h,y_t,n);			//y_h=y_t
		do
		{
			h=delta_T/k;	//Update innere Schrittweite
			for (int j = 0; j <n*n; j++)
			{
				J_h[j]=B[j]-h*DF[j];	//Update linke Seite
			}
			fcn(n,y_h,f);	
			for (int j = 0; j <n; j++)
			{
				f_h[j]=h*f[j];	//Update rechte Seite
			}
			lu_complete(J_h,delta_y,f_h,n);
			for (int j = 0; j <n; j++)
			{
				y_h[j]=y_h[j]+delta_y[j];		//y an aktueller Stelle h berechnen
			}
			i++;
		} while (i<=k);
		
		h_Werte[k-1]=h;
		for ( int j = 0; j < n; j++)
		{
			P[(k-1)*n+j]=y_h[j];
		}
	}

	z_counter=z_counter+n;
	for (int j = 0; j < n; j++)
	{
		y_t_2[j]= P[(k_steps-1)*n+j]+ ((P[(k_steps-2)*n+j] - P[(k_steps-1)*n+j]) * -h_Werte[k_steps-1])
											/(h_Werte[k_steps-2]-h_Werte[k_steps-1]);
		delta_frob_V[j]=y_t_2[j]-y_t[j];		//Differenz vektor von y_t_2-y_t
		y_j_frob_V[i]= y_t[i];					
		y_t[j]= y_t_2[j];		//überschreiben
		z[z_counter+j]=y_t[j];	//Vektor mit allen y_t werten
	}

	delta_frob= nmrc_mnorm(n,1,delta_frob_V,n,NmrcNormFrob);	//Frob von y_t_2-y_t
	error_y_j_frob[error_counter]=	nmrc_mnorm(n,1,y_j_frob_V,n,NmrcNormFrob);	//Frob von y_t
	error_frob[error_counter]= delta_frob/error_y_j_frob[error_counter];		//error_frob |y_t_2-y_t|/|y_t|
	error_counter=error_counter + 1;

	t[j_t+1]=t[j_t]+delta_T;
	j_t=j_t+1;

} while (t[j_t]<tmax);

free(y_h);
free(y_t);
free(y_t_2);
free(delta_y);
free(h_Werte);
free(DF);
free(P);
free(f);

}


//Abgabe 2 beginnt hier
//Test
//Für genau implementierung siehe euler_semi_implizit_Schritt
////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////DAE-SOLVER/////////////////////////////////////
void Richardson_extrapolation(int n, int k_steps, double *h_Werte, double *P, double *y_t_2){
for (int j = 0; j < n; j++)
	{
		y_t_2[j]= P[(k_steps-1)*n+j]+ ((P[(k_steps-2)*n+j] - P[(k_steps-1)*n+j]) * -h_Werte[k_steps-1])
											/(h_Werte[k_steps-2]-h_Werte[k_steps-1]);
	}
}


void euler_semi_implizit_Schritt(int n,double delta_T, double *B, double *t,double *y_t, double *y_t_2, double *z, void(fcn)(int, double*,double*)){

int i, k_steps=2;
double h; //delta_T;	//Schrittweite eines Basisschrittes
double h_jac=pow(10.,-8.);          //Schrittweite jacobi Matrix
double tol= pow(10.,-6);            //Toleranz

double *y_h, *delta_y, *J_h, *DF, *P, *f, *h_Werte;
//double *y_t, *y_t_2;
y_h=calloc(n,sizeof(double));
//y_t=calloc(n,sizeof(double));
//y_t_2=calloc(n,sizeof(double));
delta_y= calloc(n,sizeof(double));
J_h=calloc(n*n,sizeof(double));
h_Werte=calloc(k_steps,sizeof(double));            
DF = calloc(n*n,sizeof(double));		//Jacobi Matrix
P = calloc(n*k_steps,sizeof(double));
f=calloc(n,sizeof(double));//funktionswerte
double *a=calloc(n,sizeof(double));

double f_h[n];				//f_h[i]=h*f[i];
fjac_cd(n,h_jac,y_t,DF,fcn);          //f'(y^t)=DF
	for (int k = 1; k <= k_steps; k++)
	{
		i=0;
		copy_array(y_h,y_t,n);			//y_h=y_t
		do
		{
			h=delta_T/k;	//Update innere Schrittweite
			for (int j = 0; j <n*n; j++)
			{
				J_h[j]=B[j]-h*DF[j];	//Update linke Seite
			}
			fcn(n,y_h,f);	
			for (int j = 0; j <n; j++)
			{
				f_h[j]=h*f[j];	//Update rechte Seite
			}
			lu_complete(J_h,a,f_h,n);
			for (int j = 0; j <n; j++)
			{
				y_h[j]=y_h[j]+a[j];		//y an aktueller Stelle h berechnen
			}
			i++;
		} while (i<=k);
		
		h_Werte[k-1]=h;
		for ( int j = 0; j < n; j++)
		{
			P[(k-1)*n+j]=y_h[j];
		}
	}
	//Richardson_extrapolation(n,k_steps,h_Werte,P,y_t_2);
	for (int j = 0; j < n; j++)
	{
		y_t_2[j]= P[(k_steps-1)*n+j]+ ((P[(k_steps-2)*n+j] - P[(k_steps-1)*n+j]) * -h_Werte[k_steps-1])
											/(h_Werte[k_steps-2]-h_Werte[k_steps-1]);
	}


free(J_h);
free(y_h);
free(delta_y);
free(h_Werte);
free(DF);
free(P);
free(f);
}


//noch bearbeiten
void euler_SI(int n,int *j_t, int *z_counter,double delta_T, double *B, double *t,double tmax, double *y_i, double *z, void(fcn)(int, double*,double*)){
	
	int j_t1,z_counter1;
	j_t1=*j_t;
	z_counter1=*z_counter;
	double *y_t,*y_t_2;
	y_t=calloc(n,sizeof(double));
	y_t_2=calloc(n,sizeof(double));
	copy_array(y_t,y_i,n);	//alles auf y_t
	
	for(int i=0;i<n;i++){
	z[z_counter1+i]=y_i[i];	//Gesamtvektor mit allen y_werten
	}
	z_counter1+=n;  
	do
	{
		euler_semi_implizit_Schritt(n,delta_T,B,t,y_t,y_t_2,z,fcn);

		for(int i=0;i<n;i++)
		{
			y_t[i]=y_t_2[i];
			z[z_counter1+i]=y_t[i];	//Gesamtvektor mit allen y_werten
		}	
		z_counter1 +=n;
		t[j_t1+1]=t[j_t1]+delta_T;
		j_t1++;
	} while (t[j_t1]<tmax);
	*j_t=j_t1;
	*z_counter=z_counter1;
}
//Funktionert
void euler_SI_dense_output(int n,int anzahl,int *j_t, int *z_counter, double delta_T, double *B, double *t,double tmax, double *y_i, double *z, void(fcn)(int, double*,double*)){
	
	int j_t1,z_counter1;
	j_t1=*j_t;
	z_counter1=*z_counter;
	double *y_t,*y_t_2, *f_t,*f_t_2;
	y_t=calloc(n,sizeof(double));
	y_t_2=calloc(n,sizeof(double));
	//f_t=calloc(n,sizeof(double));
	//f_t_2=calloc(n,sizeof(double));
	copy_array(y_t,y_i,n);	//alles auf y_t
	
	for(int i=0;i<n;i++){
	z[j_t1+i]=y_i[i];	//Gesamtvektor mit allen y_werten
	}

	t[j_t1+1]=delta_T/anzahl;
	z_counter1+=n; 

	do
	{
		euler_semi_implizit_Schritt(n,delta_T,B,t,y_t,y_t_2,z,fcn);
		//Ändere Zeile drei,der weitergegeben Funktion;
		//Funktion von Modell eins soll weitergegeben werden.

		hermit_interpolation_Schritt(n,j_t1,z_counter1,anzahl,delta_T,t,y_t,y_t_2,z,fcn);

		z_counter1+=n*anzahl;
		j_t1+=anzahl;

		for(int i=0;i<n;i++)
		{
			y_t[i]=y_t_2[i];
			z[z_counter1+i]=y_t[i];	//Gesamtvektor mit allen y_werten
		}	
		z_counter1 +=n;
		t[j_t1+1]=t[j_t1-anzahl]+delta_T;
		j_t1++;
	} while (t[j_t1]<tmax);
	*j_t=j_t1;
	*z_counter=z_counter1;

} 

//bearbeiten
void euler_SI_DO_abbruch(int n, int *j_t, int *z_counter, double delta_T, double *B, double *t,double tmax, double ymax, double *y_i, double *z, void(fcn)(int, double*,double*)){
//x=Angabe Modell; 
	
	int j_t1,z_counter1;
	j_t1=*j_t;
	z_counter1=*z_counter;
	double *y_t,*y_t_2, *f_t,*f_t_2;
	y_t=calloc(n,sizeof(double));
	y_t_2=calloc(n,sizeof(double));
	//f_t=calloc(n,sizeof(double));
	//f_t_2=calloc(n,sizeof(double));
	copy_array(y_t,y_i,n);	//alles auf y_t
	
	for(int i=0;i<n;i++){
	z[j_t1+i]=y_i[i];	//Gesamtvektor mit allen y_werten
	}

	t[j_t1+1]+=delta_T;
	z_counter1+=n;

	int a_switch=0;
	//&&Regler==0
	do
	{
		euler_semi_implizit_Schritt(n,delta_T,B,t,y_t,y_t_2,z,fcn);

		if (y_t_2[n-1]>ymax)
		{
			a_switch=1;
	//void hermit_interpolation_t(int n,int j_t,double ymax,double delta_T, double *t, double *y_t, double *y_t_2, double *y,void(fcn)(int, double*,double*)){
			hermit_interpolation_t(n,j_t1,ymax,delta_T,t,y_t,y_t_2,y_i,fcn);
			for(int i=0;i<n;i++)
			{
			z[z_counter1+i]=y_i[i];	//Gesamtvektor mit allen y_werten
			}
			//copy_array(y_t,y_i,n);
		}else
		{
			for(int i=0;i<n;i++)
			{
				y_t[i]=y_t_2[i];
				z[z_counter1+i]=y_t[i];	//Gesamtvektor mit allen y_werten
				t[j_t1+1]=t[j_t1]+delta_T;
			}	
		}
			
		z_counter1 +=n;
		j_t1++;
	} while (t[j_t1]<tmax&&a_switch==0);
	*j_t=j_t1;
	*z_counter=z_counter1;
};

void euler_SI_DAE_In_2(int n, int *j_t, int *z_counter, double delta_T, double *B, double *t,double tmax, double ymax, double *y_i, double *z, void(fcn)(int, double*,double*)){
	int j_t1,z_counter1;
	j_t1=*j_t;
	z_counter1=*z_counter;
	double *y_t,*y_t_2,*f_t_2;
	y_t=calloc(n,sizeof(double));
	y_t_2=calloc(n,sizeof(double));
	//f_t=calloc(n,sizeof(double));
	f_t_2=calloc(n,sizeof(double));
	copy_array(y_t,y_i,n);	//alles auf y_t
	double tol= pow(10.,-8);            //Toleranz

	for(int i=0;i<n;i++){
	z[z_counter1+i]=y_i[i];	//Gesamtvektor mit allen y_werten
	}
	z_counter1+=n;  
	do
	{
		euler_semi_implizit_Schritt(n,delta_T,B,t,y_t,y_t_2,z,fcn);

		//fcn1(n,y_t_2,f_t_2);
		f_t_2[2] = 100.*(y_t_2[0]/800. + y_t_2[1]/1100.) + 10.;
		//f_t_2[2]=f_t_2[2]-20.;
	printf("f: %lf\n",fabs(f_t_2[2]));

		for(int i=0;i<n;i++)
		{
			y_t[i]=y_t_2[i];
			z[z_counter1+i]=y_t[i];	//Gesamtvektor mit allen y_werten
		}	
		z_counter1 +=n;
		t[j_t1+1]=t[j_t1]+delta_T;
		j_t1++;
	} while (t[j_t1]<tmax);
	*j_t=j_t1;
	*z_counter=z_counter1;
}
	



void dae_semi_implizit_Index_1(int n,int *j_t,int *z_counter,double delta_T, double *B, double *t, double tmax,double *y_i, double *z, void(fcn)(int, double*,double*)){

int i;
int k_steps=2;
int j_t1;
j_t1=*j_t;
int z_counter1;
z_counter1=*z_counter;

double *y_h;
y_h=calloc(n,sizeof(double));
double *y_t;
y_t=calloc(n,sizeof(double));
double *y_t_2;
y_t_2=calloc(n,sizeof(double));
double *delta_y;
delta_y= calloc(n,sizeof(double));
double *J_h;
J_h=calloc(n*n,sizeof(double));
double h;
//delta_T;	//Schrittweite eines Basisschrittes
double *h_Werte;
h_Werte=calloc(k_steps,sizeof(double));
double h_jac=pow(10.,-8.);          //Schrittweite jacobi Matrix
double tol= pow(10.,-6);            //Toleranz

double *DF;                         //Jacobi Matrix
DF = calloc(n*n,sizeof(double));

double *P;
P = calloc(n*k_steps,sizeof(double));

double *f;					//funktionswerte
double f_h[n];				//f_h[i]=h*f[i];
f=calloc(n,sizeof(double));


for(i=0;i<n;i++){
	z[j_t1+i]=y_i[i];	//Gesamtvektor mit allen y_werten
}
z_counter1=z_counter1+n;         

copy_array(y_t,y_i,n);	//alles auf y_t

do
{
    fjac_cd(n,h_jac,y_t,DF,fcn);          //f'(y^t)=DF
	for (int k = 1; k <= k_steps; k++)
	{
		i=0;
		copy_array(y_h,y_t,n);			//y_h=y_t
		do
		{
			h=delta_T/k;	//Update innere Schrittweite
			for (int j = 0; j <n*n; j++)
			{
				J_h[j]=B[j]-h*DF[j];	//Update linke Seite
			}
			fcn(n,y_h,f);	
			for (int j = 0; j <n; j++)
			{
				f_h[j]=h*f[j];	//Update rechte Seite
			}
			lu_complete(J_h,delta_y,f_h,n);
			for (int j = 0; j <n; j++)
			{
				y_h[j]=y_h[j]+delta_y[j];		//y an aktueller Stelle h berechnen
			}
			i++;
		} while (i<=k);
		
		h_Werte[k-1]=h;
		for ( int j = 0; j < n; j++)
		{
			P[(k-1)*n+j]=y_h[j];
		}
	}
	z_counter1= z_counter1+n;

	//*z_counter=*z_counter+n;
	for (int j = 0; j < n; j++)
	{
		y_t_2[j]= P[(k_steps-1)*n+j]+ ((P[(k_steps-2)*n+j] - P[(k_steps-1)*n+j]) * -h_Werte[k_steps-1])
											/(h_Werte[k_steps-2]-h_Werte[k_steps-1]);
		y_t[j]= y_t_2[j];		//überschreiben
		z[z_counter1+j]=y_t[j];	//Vektor mit allen y_t werten
	}

	t[j_t1+1]=t[j_t1]+delta_T;
	j_t1++;

} while (t[j_t1]<tmax);

*j_t=j_t1;
*z_counter=z_counter1;
free(J_h);
free(y_h);
free(y_t);
free(y_t_2);
free(delta_y);
free(h_Werte);
free(DF);
free(P);
free(f);
}

void dae_si_dense_output(int n,int *j_t,int *z_counter,double delta_T, double *B, double *t, double tmax,double *y_i, double ymax, double *z, void(fcn)(int, double*,double*)){
int i;
int k_steps=2;
int j_t1;
j_t1=*j_t;
int z_counter1;
z_counter1=*z_counter;
//int j_t=0;		//t[j_t+1]=t[j_t]+delta_T
//int z_counter=0;
int Regler=0;

double *y_h;
y_h=calloc(n,sizeof(double));
double *y_t;
y_t=calloc(n,sizeof(double));
double *y_t_2;
y_t_2=calloc(n,sizeof(double));
double *delta_y;
delta_y= calloc(n,sizeof(double));
double *J_h;
J_h=calloc(n*n,sizeof(double));
double h;
//delta_T;	//Schrittweite eines Basisschrittes
double *h_Werte;
h_Werte=calloc(k_steps,sizeof(double));
double h_jac=pow(10.,-8.);          //Schrittweite jacobi Matrix
double tol= pow(10.,-6);            //Toleranz

double *DF;                         //Jacobi Matrix
DF = calloc(n*n,sizeof(double));

double *P;
P = calloc(n*k_steps,sizeof(double));

double *f;					//funktionswerte
double f_h[n];				//f_h[i]=h*f[i];
f=calloc(n,sizeof(double));

double *y;
y=calloc(n,sizeof(double));

for(i=0;i<n;i++){
	z[z_counter1+i]=y_i[i];	//Gesamtvektor mit allen y_werten
}         
z_counter1=z_counter1+n;
copy_array(y_t,y_i,n);	//alles auf y_t

do
{
    fjac_cd(n,h_jac,y_t,DF,fcn);          //f'(y^t)=DF
	for (int k = 1; k <= k_steps; k++)
	{
		i=0;
		copy_array(y_h,y_t,n);			//y_h=y_t
		do
		{
			h=delta_T/k;	//Update innere Schrittweite
			for (int j = 0; j <n*n; j++)
			{
				J_h[j]=B[j]-h*DF[j];	//Update linke Seite
			}
			fcn(n,y_h,f);	
			for (int j = 0; j <n; j++)
			{
				f_h[j]=h*f[j];	//Update rechte Seite
			}
			lu_complete(J_h,delta_y,f_h,n);
			for (int j = 0; j <n; j++)
			{
				y_h[j]=y_h[j]+delta_y[j];		//y an aktueller Stelle h berechnen
			}
			i++;
		} while (i<=k);
		
		h_Werte[k-1]=h;
		for ( int j = 0; j < n; j++)
		{
			P[(k-1)*n+j]=y_h[j];
		}
	}
	//z_counter1=z_counter1+n;
	for (int j = 0; j < n; j++)
	{
		y_t_2[j]= P[(k_steps-1)*n+j]+ ((P[(k_steps-2)*n+j] - P[(k_steps-1)*n+j]) * -h_Werte[k_steps-1])
											/(h_Werte[k_steps-2]-h_Werte[k_steps-1]);
	}
	if (y_t_2[n-1]>ymax)
	{
		//hermit_interpolation(n,j_t1,ymax,delta_T,t,y_t,y_t_2,y,fcn);
		for (int j = 0; j < n; j++)
		{
			z[z_counter1+j]=y[j];	//Vektor mit allen y_t werten
			//printf("z[%i]=%lf\n",j,z[z_counter1+j]);
		}
		Regler=1;
		}else if (y_t_2[n-1]<ymax)
		{
			for (int j = 0; j < n; j++)
			{
				y_t[j]= y_t_2[j];		//überschreiben
				z[z_counter1+j]=y_t[j];	//Vektor mit allen y_t werten
			}
				t[j_t1+1]=t[j_t1]+delta_T;
		}
		j_t1=j_t1+1;
		z_counter1=z_counter1+n;

} while (t[j_t1]<tmax&&Regler==0);
*j_t=j_t1;
*z_counter=z_counter1;
free(J_h);
free(y_h);
free(y_t);
free(y_t_2);
free(delta_y);
free(h_Werte);
free(DF);
free(P);
free(f);

}

void dae_si(int n,int *j_t,int *z_counter,double delta_T, double *B, double *t, double tmax,double *y_i, double *z, void(fcn)(int, double*,double*)){

int i;
int k_steps=2;
int j_t1;
j_t1=*j_t;
int z_counter1;
z_counter1=*z_counter;
double *y_h;
y_h=calloc(n,sizeof(double));
double *y_t;
y_t=calloc(n,sizeof(double));
double *y_t_2;
y_t_2=calloc(n,sizeof(double));
double *delta_y;
delta_y= calloc(n,sizeof(double));
double *J_h;
J_h=calloc(n*n,sizeof(double));
double h;
double *h_Werte;
h_Werte=calloc(k_steps,sizeof(double));
double h_jac=pow(10.,-8.);          //Schrittweite jacobi Matrix
double tol= pow(10.,-6);            //Toleranz

double *DF;                         //Jacobi Matrix
DF = calloc(n*n,sizeof(double));

double *P;
P = calloc(n*k_steps,sizeof(double));

double *f;					//funktionswerte
double f_h[n];				//f_h[i]=h*f[i];
f=calloc(n,sizeof(double));

double *y;
y=calloc(n,sizeof(double));

copy_array(y_t,y_i,n);	//alles auf y_t

do
{
    fjac_cd(n,h_jac,y_t,DF,fcn);          //f'(y^t)=DF
	for (int k = 1; k <= k_steps; k++)
	{
		i=0;
		copy_array(y_h,y_t,n);			//y_h=y_t
		do
		{
			h=delta_T/k;	//Update innere Schrittweite
			for (int j = 0; j <n*n; j++)
			{
				J_h[j]=B[j]-h*DF[j];	//Update linke Seite
			}
			fcn(n,y_h,f);	
			for (int j = 0; j <n; j++)
			{
				f_h[j]=h*f[j];	//Update rechte Seite
			}
			lu_complete(J_h,delta_y,f_h,n);
			for (int j = 0; j <n; j++)
			{
				y_h[j]=y_h[j]+delta_y[j];		//y an aktueller Stelle h berechnen
			}
			i++;
		} while (i<=k);
		
		h_Werte[k-1]=h;
		for ( int j = 0; j < n; j++)
		{
			P[(k-1)*n+j]=y_h[j];
		}
	}
	//z_counter1=z_counter1+n;
	for (int j = 0; j < n; j++)
	{
		y_t_2[j]= P[(k_steps-1)*n+j]+ ((P[(k_steps-2)*n+j] - P[(k_steps-1)*n+j]) * -h_Werte[k_steps-1])
											/(h_Werte[k_steps-2]-h_Werte[k_steps-1]);
	}

	for (int j = 0; j < n; j++)
	{
		y_t[j]= y_t_2[j];		//überschreiben
		z[z_counter1+j]=y_t[j];	//Vektor mit allen y_t werten
	}
	t[j_t1+1]=t[j_t1]+delta_T;
		
	j_t1=j_t1+1;
	z_counter1=z_counter1+n;

} while (t[j_t1]<tmax);
*j_t=j_t1;
*z_counter=z_counter1;
free(J_h);
free(y_h);
free(y_t);
free(y_t_2);
free(delta_y);
free(h_Werte);
free(DF);
free(P);
free(f);
}


void dae_si_Schritt(int n,int *j_t,int *z_counter,double delta_T, double *B, double *t, double tmax,double *y_i, double *z, void(fcn)(int, double*,double*),void(fcn1)(int, double*,double*)){

int i;
int k_steps=2;
int Regler=0;
int j_t1;
j_t1=*j_t;
int z_counter1;
z_counter1=*z_counter;
double *y_h;
y_h=calloc(n,sizeof(double));
double *y_t;
y_t=calloc(n,sizeof(double));
double *y_t_2;
y_t_2=calloc(n,sizeof(double));
double *delta_y;
delta_y= calloc(n,sizeof(double));
double *J_h;
J_h=calloc(n*n,sizeof(double));
double h;
double *h_Werte;
h_Werte=calloc(k_steps,sizeof(double));
double h_jac=pow(10.,-8.);          //Schrittweite jacobi Matrix
double tol= pow(10.,-8);            //Toleranz

double *DF;                         //Jacobi Matrix
DF = calloc(n*n,sizeof(double));

double *P;
P = calloc(n*k_steps,sizeof(double));

double *f;					//funktionswerte
double f_h[n];				//f_h[i]=h*f[i];
f=calloc(n,sizeof(double));

double *y;
y=calloc(n,sizeof(double));

double *f_2;
f_2=calloc(n,sizeof(double));
copy_array(y_t,y_i,n);	//alles auf y_t

do
{
    fjac_cd(n,h_jac,y_t,DF,fcn);          //f'(y^t)=DF
	for (int k = 1; k <= k_steps; k++)
	{
		i=0;
		copy_array(y_h,y_t,n);			//y_h=y_t
		do
		{
			h=delta_T/k;	//Update innere Schrittweite
			for (int j = 0; j <n*n; j++)
			{
				J_h[j]=B[j]-h*DF[j];	//Update linke Seite
			}
			fcn(n,y_h,f);	
			for (int j = 0; j <n; j++)
			{
				f_h[j]=h*f[j];	//Update rechte Seite
			}
			lu_complete(J_h,delta_y,f_h,n);
			for (int j = 0; j <n; j++)
			{
				y_h[j]=y_h[j]+delta_y[j];		//y an aktueller Stelle h berechnen
			}
			i++;
		} while (i<=k);
		
		h_Werte[k-1]=h;
		for ( int j = 0; j < n; j++)
		{
			P[(k-1)*n+j]=y_h[j];
		}
	}

	for (int j = 0; j < n; j++)
	{
		y_t_2[j]= P[(k_steps-1)*n+j]+ ((P[(k_steps-2)*n+j] - P[(k_steps-1)*n+j]) * -h_Werte[k_steps-1])
											/(h_Werte[k_steps-2]-h_Werte[k_steps-1]);
	}
	fcn1(n,y_t_2,f_2);
	f_2[2] = 100.*(y_t_2[0]/800. + y_t_2[1]/1100.) + 10;
	f_2[2]=f_2[2]-20.;
	printf("f: %lf\n",fabs(f_2[2]));
/*

	if(20-f_2[2] > tol){
		printf("f: %lf\n",f_2[2]);
	}
	*/

	/*
	

	if (fabs(f_2[2])>tol)
	{
		printf("Fehler: %lf\n",f_2[2]);
	
	}
	*/

	//wichtig
	//fcn1(n,y_t,f_2);

	for (int j = 0; j < n; j++)
	{
		y_t[j]= y_t_2[j];		//überschreiben
		z[z_counter1+j]=y_t[j];	//Vektor mit allen y_t werten
	}
	t[j_t1+1]=t[j_t1]+delta_T;
		
	j_t1=j_t1+1;
	z_counter1=z_counter1+n;

} while (t[j_t1]<tmax&&Regler==0);
*j_t=j_t1;
*z_counter=z_counter1;
free(J_h);
free(y_h);
free(y_t);
free(y_t_2);
free(delta_y);
free(h_Werte);
free(DF);
free(P);
free(f);
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////Hermit Interpolation/////////////////////////////////////

void hermit_interpolation_Schritt(int n,int j_t1, int z_counter1, int anzahl,double delta_T, double *t, double *y_t, double *y_t_2, double *z,void(fcn)(int, double*,double*)){
	//theta=(t_t-t[j_t])/delta_T;
	//(t[j_t+1]-t[j_t])/Anzahl/delta_T=1/Anzahl
	int i_1=0;	//Schrittzähler für Anzahl	
	double *y,*theta,*f_t,*f_t_2;
	y=calloc(n,sizeof(double));
	theta=calloc(anzahl,sizeof(double));
	f_t=calloc(n,sizeof(double));
	f_t_2=calloc(n,sizeof(double));	
	fcn(n,y_t,f_t);
	f_t[n-1]= 1. + 0.1 * 100 * y_t[n-3]*(1/1100 - 1/800);
	fcn(n,y_t_2,f_t_2);
	f_t_2[n-1]= 1. + 0.1 * 100 * y_t_2[n-3]*(1/1100 - 1/800);


		do
		{
			theta[i_1]=(i_1+1.)/(anzahl+1.);
			t[j_t1+1]=t[j_t1]+theta[0]*delta_T;
			for (int j = 0; j < n; j++)
			{
				y[j] = (1-theta[i_1])*y_t[j] + theta[i_1]*y_t_2[j] + theta[i_1]*(theta[i_1]+1)*((1-2*theta[i_1])*(y_t_2[j]-y_t[j]) 
               		 		+ (theta[i_1]-1)*(delta_T)*f_t[j] + theta[i_1]*(delta_T)*f_t_2[j]);
				z[z_counter1+j]=y[j];
			}
			i_1++;
			j_t1++;
			z_counter1+=n;
		} while (i_1<anzahl);

	free(y);
	free(theta);
	free(f_t);
	free(f_t_2);
}

void hermit_interpolation_t(int n,int j_t,double ymax,double delta_T, double *t, double *y_t, double *y_t_2, double *y,void(fcn)(int, double*,double*)){

	int Regler=0;
	int k=0;
	double ak=0.1;
	double tol=10e-9;
	int l_b=0;		//lower boundery
	int u_b=10;		//upper boundery
	
	double theta=0;
	double *f_t,*f_t_2,*y_alt;
	f_t=calloc(n,sizeof(double));
	f_t_2=calloc(n,sizeof(double));
	y_alt=calloc(n,sizeof(double));	//temporär

	fcn(n,y_t,f_t);
	fcn(n,y_t_2,f_t_2);
	
	do
	{
		for (int j = l_b; j < u_b; j++)
		{
			theta=1.0*j*ak;	//anfangs theta=0
			copy_array(y_alt,y,n);
			for (int i = 0; i < n; i++)
			{
				y[i] = (1-2*theta)*(y_t_2[i]-y_t[i]) + delta_T*(f_t[i]*(theta-1)+f_t_2[i]*theta);
				y[i] *=	theta*(theta+1);
				y[i] += (1-theta)*y_t[i] + theta*y_t_2[i];
			}

			//wenn y_alt[2]<ymax<y[2]
			if (ymax<y[n-1]&&ymax>y_alt[n-1])
			{
				k++;
				if (ymax>y_alt[n-1]&&ymax-tol<y_alt[n-1])
				{
					Regler=1;
				}else
				{
					ak *=0.1;
					l_b=10*(j-1);
					u_b=10*j;
				}
			}
		}
	} while (k<50&&Regler==0);

	for (int i = 0; i < n; i++)
			{
				y[i]=y_alt[i];
			}
	t[j_t+2]=t[j_t]+ delta_T;
	t[j_t+1]=t[j_t]+ theta*delta_T;
	printf("t[j_t+1]=%lf\n",t[j_t+1]);
	printf("k=%i\n",k);

	free(f_t);
	free(f_t_2);
	free(y_alt);
}

void hermit_interpolation(int n,int j_t1, int z_counter1, int anzahl,double delta_T, double *t, double *y_t, double *y_t_2, double *z,void(fcn)(int, double*,double*)){
	//theta=(t_t-t[j_t])/delta_T;
	//(t[j_t+1]-t[j_t])/Anzahl/delta_T=1/Anzahl
	int i_1=0;	//Schrittzähler für Anzahl	
	double *y,*theta,*f_t,*f_t_2;
	y=calloc(n,sizeof(double));
	theta=calloc(anzahl,sizeof(double));
	f_t=calloc(n,sizeof(double));
	f_t_2=calloc(n,sizeof(double));	
	fcn(n,y_t,f_t);
	fcn(n,y_t_2,f_t_2);

		do
		{
			theta[i_1]=(i_1+1.)/(anzahl+1.);
			t[j_t1+1]=t[j_t1]+theta[0]*delta_T;
			for (int j = 0; j < n; j++)
			{
				y[j] = (1-theta[i_1])*y_t[j] + theta[i_1]*y_t_2[j] + theta[i_1]*(theta[i_1]+1)*((1-2*theta[i_1])*(y_t_2[j]-y_t[j]) 
               		 		+ (theta[i_1]-1)*(delta_T)*f_t[j] + theta[i_1]*(delta_T)*f_t_2[j]);
				z[z_counter1+j]=y[j];
			}
			i_1++;
			j_t1++;
			z_counter1+=n;
		} while (i_1<anzahl);

	free(y);
	free(theta);
	free(f_t);
	free(f_t_2);
}



////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////Abgabe 3/////////////////////////////////////
void matrix_trans(int n, int m, double *matrix){
	double *matrix_2;
	matrix_2=calloc(n*m,sizeof(double));
	copy_array(matrix_2,matrix,n*m);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            matrix[j*n+i] = matrix_2[i*m+j];
        }
    }
	free(matrix_2);
}






////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////Jakobi/////////////////////////////////////

// Funktionen zur Numerischen Berechnung der Jacobi-Matrix (Zentraldifferenz)
void fjac_cd(int n, const double h, double *x, double *jac, void (*fcn)(int, double *, double *)){
    
    int i,j;
    //double help;
    
    double *fhelp1, *fhelp2, *xhelp1 ,*xhelp2;
    
    fhelp1=malloc(2*n*sizeof(fhelp1));
    fhelp2=fhelp1+n;
    xhelp1=malloc(2*n*sizeof(xhelp1));
    xhelp2=xhelp1+n;
    

        //Jakobi-Matrix wird spaltenweise durchlaufen
        for ( j = 0; j < n ; j++ ) {
                
                //Update des x-Vektors
                for ( i = 0; i < n ; i++){
                   xhelp1[i] = x[i];
                   xhelp2[i] = x[i];
                }
                 xhelp1[j] =  xhelp1[j]-h;
                 xhelp2[j] =  xhelp2[j]+h;
                
                
                //Auswertung der Funktionen an der Stelle x+h
                fcn(n, xhelp2, fhelp2);
                
                //Auswertung der Funktionen an der Stelle x+h
                fcn(n, xhelp1, fhelp1);
                
                
                //Zeile der Jacobi-Matrix df_i/dx_k
                for ( i = 0; i < n ; i++){
                    
                    jac[j+i*n] = (fhelp2[i]-fhelp1[i])/(2*h);
                    
                }
        }
        
    free(fhelp1);
    free(xhelp1);
}




/////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Newton-Verfahren //////////////////////////////////

//Newton-Verfahren mit Betrag des Residuums als Abbruchkriterium
void newton_iter(int n, int *ni, int kmax, double tol, double *x_start, double *residuum, double *lsg, void(fcn)(int, double*, double*)){

    int k = 0;

    double *A;      //Matrix A des LGS (Jacobi-Matrix-Auswertung)
    double *b;      //rechte Seite des LGS (Funktionsauswertung)
    double *x;      //Lösungsvektor des LGS Ax=b
    double res;     //Residuum für Abbruch der Newton-Iteration
    
    
    A  = calloc(n*n,sizeof(double));
    b  = calloc(n,sizeof(double));
    x  = calloc(n,sizeof(double));

    //Startwert: Newton-Iteration festlegen
    copy_array(lsg, x_start, n);

    //Beginn: Newton-Iteration
    do{

        //Aufruf: Berechnung f(x1,x2)            
        fcn(n, lsg, b);

        //Residuum: Betrag von f(x1,x2)
        res = vec_fabs(b, n);
        residuum[k] = res;
        
        //Aufruf: Berechnung Jacobi-Matrix von f(x1,x2)
        fjac_cd(n, 10e-5, lsg, A, fcn);
        
        //Lu-Zerlegung für Newton-Verfahren
        lu_complete(A, x, b, n);


        //Vorschrift: Newton-Verfahren
        vec_diff(lsg, x , lsg, n);
        
        //Counter: Anzahl an Iterationen
        k++;   
        //printf("Starten: %lf datei\n", res);
    }while((res > tol) && (k<kmax));
    
    //Adressübergabe
    *ni     = k;
}



////////////////////////////////////////////////////////////////////////////////
/////////////////////////////Matrix Vektor Operationen//////////////////////////




//Ändert Vorzeichen aller Arrayeinträge
void array_minus(double *arr, int size){
    int i;
    for(i=0; i<size; i++){
        arr[i] = -arr[i];
    }  
}



//Kopiert Array 2 auf Array 1
void copy_array(double *arr1, double *arr2,int size){
    int i;
    for (i = 0; i < size; i++){
		arr1[i] = arr2[i];
	}
    
}


//Multpilizert zwei Matrizen
void mat_vec_multiplication(int n, int m, double *A, double *x,double *b)
{
    //Übergabe in main.c
    //print_vec_scientific(lsg_temp[0],k, "neuer vektor" );
    int i,j;
    double sum;
    for(i=0;i<n;i++)
    {
        sum=0;
        for(j=0;j< m;j++)
            sum = sum + A[j +i*m] * x[j];
        b[i] = sum;
    }
}


void matrix(double *A, double *B,int n){  
    int i,j;

    for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
                    A[j+i*n] = B[j+i*n];
        }
}


//Allokiert eine Matrix
void matrix_alloc(double ***A, int size1, int size2){
    int i;
    double **a;

    a=malloc(size1*sizeof(**a));
    a[0] = malloc(size1*size2*sizeof(*a));
    for(i=0;i<size1;i++)
        a[i]=(*a+size2*i);

    *A = a;
}


//Funktion für Matrixmultiplikation    
void matrix_multiplication(double *A,double *B,int n, double *C){
    int i,d,j;
    double sum;
    for(i=0;i< n;i++)
    {
        for(d=0;d< n;d++)
        {
            sum=0;
            for(j=0;j< n;j++)
                sum = sum + A[j +i*n] * B[d + j*n ];
            C[d + i*n] = sum;
        }
    }
    
}


//Subtrahiert vektor arr2 von arr1 -> erg = arr1 - arr2
void vec_diff(double *arr1, double *arr2, double *erg, int size){
    int i;
    for(i=0; i<size; i++){
        erg[i] = arr1[i]- arr2[i];
    }
}




//Betrag eines Vektors
double vec_fabs(double *vec, int size){
    double sum=0;
    int i;
    for(i=0; i<size; i++){
        sum += vec[i] * vec[i];
    }
    return sqrt(sum);
}


///////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////Print Befehle////////////////////////////////////

// //Voraussetzung für Gaußeliminierung: reguläre, also quadratische Matrix mit vollen Rang: m=n
//void LR_separation(int n, double *A,int lda, double *b, double *L,double *R, double *x){
void LR_separation(int n, double *A,int lda, double *b, double *L,double *R){
    
    int k,i,j,p,s;
    double  Maximum,
            sum, *A_initial,
            *A_help, *b_help,
            *L_help, *L_help_multiplication, *L_help_pivot, *L_help_invers, *L_invers,
            *P, *P_help, *P_help_multiplication; //Matrix für Pivotisierung
    
     //Allokieren
     
    A_initial = calloc( n * n, sizeof(double));
    b_help = calloc( n, sizeof(double));
    A_help = calloc( n * n, sizeof(double));
    L_help = calloc( n * n, sizeof(double));
    L_help_invers = calloc( n * n, sizeof(double));
    L_help_multiplication = calloc( n * n, sizeof(double));
    L_help_pivot = calloc( n * n, sizeof(double));
    P_help = calloc( n * n, sizeof(double));
    P_help_multiplication = calloc( n * n, sizeof(double));
    L_invers = calloc( n * n, sizeof(double));
    
    
    //Anfangszustand von Matrix A
    matrix(A_initial,A,n);
    
    
    //Ausgangspunkt: L ist zunächst Einheitsmatrix
    for(i=0;i<n;i++)
            {
                for(j=0;j<n;j++)
                {
                    if(i==j){

                        L[j+i*n] = 1;
                        P_help_multiplication[j+i*n] = 1;
                        
                        
                    }
                    else{
           
                        L[j+i*n] = 0;
                        P_help_multiplication[j+i*n] = 0;
                    }
                        
                }
            }
    

        /*Suche des Pivotelementes: betragsmäßig größte Zahl in einer Spalte muss auf der Hauptdiagonalen liegen*/
        for (k=0;k< n-1;k++) {
            
            
            //Hilfsmatrizen zur Bestimmung von R und L
    for(i=0;i<n;i++)
            {
                for(j=0;j<n;j++)
                {
                    if(i==j){
                        L_help[j+i*n] = 1;
                        L_help_invers[j+i*n] = 1;
                        L_help_pivot[j+i*n] = 1;
                        
                    }
                    else{
                        L_help[j+i*n] = 0;
                        L_help_invers[j+i*n] = 0;
                        L_help_pivot[j+i*n] = 0;
                    }
                        
                }
            }
            
            
            p=k;
            Maximum = fabs(A[k + n *(k)]);
            for (s=k+1;s< n;s++)
            {   
                if (fabs(A[k + n *s]) > Maximum )
                {
                    Maximum = fabs(A[k + n * s]);
                    p=s;
                }
            
            }
            
        //Definition der Umtauschmatrix P -> Pivotisierung über Matrixmultiplikation mit
             P= calloc(n * n, sizeof(double)); //erschafft Nullmatrix
        
        /*Beispiel: wird erste Zeile mit dritten Zeile getauscht wegen Pivotisierung:
         P = 0 0 1
             0 1 0
             1 0 0 
             vertauschen der Spalten eins und drei in Umtauschmatrix*/
        for (i=0;i< n;i++)
        {
            if(i==p)
            {
                        
                        P[k + n *i] = 1;
            }
        }
//             printf("In Spalte %d befindet sich das Pivolelement in Zeile %d \n", k+1,p+1);
//             printf(" \n");
                    
                for(j=0;j< n;j++) {
            for(i=0;i< n;i++)
                {
                    if(j==p)
                    { 
                        P[j + n * k] = 1;
                    }
                    else if ( j!=k && j!=p)
                    {
                        P[j + n * j] =1;
                    }
                }
                }
            
            
    //Zeilen vertauschen über Matrixmultiplikation nach dem Prinzip A=L*A -> wenn man direkt A = L*A resultiert Fehler (deshalb der Weg über B)
    matrix_multiplication(P,A,n,A_help);
    
              
    //Überschreiben von A durch Matrix A_help
    matrix(A,A_help,n);
    
    
    //Zeilen des Vektors b vertauschen nach dem Prinzip u=L b 
    for(i=0;i< n;i++)
    {
            sum=0;
            
            for(j=0;j< n;j++)
            {
                sum = sum + P[j + i* n] * b[j];
            }
            
            b_help[i] = sum;
        }
    
    //Überschreiben von b durch b_help
    for(i=0;i< n;i++)
    {
        b[i] = b_help[i];
    }
  
    
    //Ermittlung der Pivotisierungsmatrix P
    matrix_multiplication(P,P_help_multiplication,n,P_help);
    
              
    //Überschreiben von P_help_multiplication durch Matrix P
    
    matrix(P_help_multiplication,P_help,n);
            
            if(k!=0){
                //Pivotisierung auch auf L-Matrix anwenden
                matrix_multiplication(P,L_help_multiplication,n,L_help_pivot);
            
                //Änderung der l_ij Einträge in L-Matrix
                        
                        for(j=0;j<k;j++)
                {
                        for(i=j+1;i<n ;i++){
            
                            L[j+i*n] = L_help_pivot[j+i*n];
                            
                                
                        }
        
                }
            
            }
            
  
    
    //Ende der Pivotisierung:
    //Gleichungssystem der Form P*A = P*b
    
    //LR-Zerlegung der pivotisierten Matrix A
    
    for (i=k+1;i<n;i++){
        L_help[k + i*n] = A[k + i*n]/A[k + k*n]; 
        L_help_invers[k + i*n] = - A[k + i*n]/A[k + k*n]; 
    }
    
    //Multiplikn,ation der L-Matrix mit A, da gilt:
    // R = L_(n-1) * ... * L_1 * A 
    
    
    matrix_multiplication(L_help_invers,A,n,A_help);
              
    //Überschreiben von A durch Matrix A_help
    matrix(A,A_help,n);
    
    
    //Bestimmung der L-Matrix
    matrix_multiplication(L,L_help,n,L_help_multiplication);
    
              
    //Überschreiben von L durch Matrix L_help_multiplication
    matrix(L,L_help_multiplication,n);
            
} 
        
        
        //Bestimmung der inversen L-Matrix
          for(j=0;j< n;j++)
        {
            L_invers[j+j*n] = L[j+j*n];
            
    
                for(i=j+1;i<n ;i++){
    
                    L_invers[j+i*n] = -L[j+i*n];
                    
                        
                }
        
        }
        
        //Ende der LR-Zerlegung
        
        //Bestimmung von R = L^(-1)*P*A --> entspricht A
        matrix(R,A,n);
        
        //A wird wieder auf Initialbedingung resetet
        matrix(A,A_initial,n);
        

}    
  

//Vorwärtssubstitution des LGS L (R x) = L y = b mit y = R x
void forward_substitution(double *L, double *b, double *y,int n){
    
     int i,j;
     double a;
           
     
        y[0] = (b[0])/(L[0]);
        for(i= 1;i<n;i++){
            a=0;
            for(j= 0;j<i;j++){
                a = a + L[j + n *i] * y[j];
            }
             y[i] = (b[i] - a)/(L[i + n *i]);
        }
        
}    

 
//Rücksubstitution der LGS R x = y z
void back_substitution(double *R, double *y, double *x, int n){
    
    int i,j;
    double a;

        x[n -1] = (y[n -1])/(R[(n -1) + n *(n -1)]);
        for(i= n-2;i>=0;i--){
            a=0;
            for(j= n-1;j>i;j--){
                a = a + R[j + n *i] * x[j];
            }
             x[i] = (y[i] - a)/(R[i + n *i]);
        }
        


}    

//LU-Zerlegung
void lu_complete(double *A, double *x, double *b, int n){
    
    double *L  = calloc(n*n,sizeof(double)); 
    double *R  = calloc(n*n,sizeof(double));
    double *y  = calloc(n,sizeof(double));
    
    //Aufruf: LR-Zerlegung aus Assignment1
    LR_separation(n,A,n,b,L,R);
        
    //Aufruf: Vorwärtssubstitution
    forward_substitution(L, b, y, n);
        
    //Aufruf: Rückwärtssubstitution
    back_substitution(R, y, x, n);
    
    free(L);
    free(R);
    free(y);
}



/////////////////////////////////////// NMRC-LIBRAY //////////////////////////////////////////////////

/* -------------------------------------------------------------------------direct linear systems solver -------------------------------------------------------------------------*/

//	n	Dimension des Gleichungssystems
 // lda	führende Dimension von A
 // a	Koeffizientenmatrix A
 // b	Rechte Seite b
 //	x	Lösungsvektor

int nmrc_gauss(int n, double *a, int lda, double *b, double *x)
{

	int i, k, j, p, nm1, print_flg;
	double t, big;

	if (n < 2 || lda < n) return -2;
	if (!a || !b || !x) return -3;
	
	// print_flg = print_iter_flag & PRINT_ITER_LS;
	   
	/* Gauss-Elimination mit Spaltenpivotierung */

	nm1 = n - 1;
        for (k = 0; k<nm1 ; k++) {

		/* Pivotsuche */
		big = 0.0;
		p = k;
		for (i=k; i<n; i++) {
			t = fabs(a[i*lda+k]);
			if (t > big) {
				p = i;
				big = t;
			}
		}
		
		/* Test fuer die Singularitaet */
		if ( big < TOL ) {
			if (big == ZERO)
				return -1; /* Matrix ist singulaer */
			
			printf("Warning %s: Matrix is close to singular\n", __func__);
		}
		/* Zeile p und k tauschen */
		if (p != k) {
			for (j=k; j<n; j++) {
				t = a[k*lda +j];
				a[k*lda +j] = a[p*lda +j];
				a[p*lda +j] = t;
			}
			t = b[k];
			b[k] = b[p];
			b[p] = t;
		}
		
		/* Aufdatierung */
		for (i = k + 1; i<n; i++) {
	    		t = a[i*lda + k] / a[k*lda+k];
	    		for (j = k + 1; j<n; j++) 
				a[i*lda + j] -= t * a[k*lda + j];
	    		b[i] -= t * b[k];
		}
		if (print_flg) {
                        printf("\n%d. Umformungsschritt:\n", k+1);
                       // mat_print(" Matrix A", "%13g", stdout, n, n, lda, a);
                        //vec_print(" Vektor b", "%13g", stdout, n, n, b);			
		}
	}
	/* Test fuer die Singularitaet */
	if ((big = fabs(a[lda*nm1 + nm1])) < TOL) {
		if (big == ZERO)
			return -1;
		printf("Warning %s: Matrix is close to singular\n", __func__);
	}

	/* Rueckwaertssubstitution */
	for (i=nm1; i>=0; i--) {
		t = b[i];
		for (j=i+1; j<n; j++) 
	    		t -= a[i*lda + j] * x[j];
		x[i] = t / a[i*lda +i];
	}
	return 0;
	
}	/* gauss */


/*-------------------------------------------------------------------Hilfsfunktionen-----------------------------------------------------------------------*/

 //   ls = Parameter fuer lineares System
 //     .kmax   = Maximale Anzahl von Iterationen
 //     .rtol   = gewuenschte Genauigkeit Norm2(r_k)/Norm2(r_0)
 //     .relax  = Relaxationsfaktor (0, 2]
 //   n  = Zeilenanzahl der Koeffizientenmatrix
  //  b  = Rechte Seite b
  //  x  = Startvektor fuer die gesuchte Loesung
    
//  Parameter nur fuer vollbesetzte Matrix A:
//    a   = Matrix A (n x lda)
//    lda = fuehrende Dimension der Matrix A

//  Parameter nur fuer duennbesetzte Matrix A:
//    a   = Vektor Nichtnullelemente von A, in CRS-Format gespeichert.
//    ir  = "row-pointer" Zeiger auf den Beginn der i-ten Zeile. Vektor der dimension n+1.
//    ic  = "col-index" Splatenindizes der entsprechenden Wert in A-Vektor.



int testIvArgs(int n, const double *a, int km, double tol, const double *b, const double *x)
{
	if (n < 2 || km < 0) return -3;
	if (a == NULL || b == NULL || x == NULL) return -4;
	if (tol < ZERO) return -5;
	return 0;
}

/* Kehrwerte der Diagonalelemente fuer CRS-Matrix */
int crsInvDiag(int n, const double *a, const int *ir, const int *ic, double *ia)
{
	int	i, l;
	double	diag;
	
	for (i = 0; i < n; i++) {
		for (l = ir[i]; l < ir[i+1]; l++) {
			if (ic[l] != i) continue; /* kein Diagonalelement */
			diag = a[l];
			if (fabs(diag) == ZERO) return 2;
			ia[i] = 1.0 / diag;
			break;
		}
		/* kein Diagonalelement in aktueller Zeile */
		if ( l == ir[i+1] ) return 1;
	}
	return 0;	
}
/* Kehrwerte der Diagonalelemente fuer vollbesetzte Matrix */
static int denseInvDiag(int n, const double *a, int lda, double *ia)
{
	register int i, np1;
	register double t;
	
	np1 = lda + 1;
	for (i = 0; i < n; i++) {
		t = a[i*np1];
		/* Test auf Singularitaet */
		if (fabs(t) == ZERO)  return 1;	
		ia[i] = 1. / t; /* omega / a(i,i) */
	}
	return 0;
}


  /*--------------------------------------------------------------- dense/sparse splitting linear system solvers--------------------------------------------------------------- */

 //   ls = Parameter fuer lineares System
 //     .kmax   = Maximale Anzahl von Iterationen
 //     .rtol   = gewuenschte Genauigkeit Norm2(r_k)/Norm2(r_0)
 //     .relax  = Relaxationsfaktor (0, 2]
 //   n  = Zeilenanzahl der Koeffizientenmatrix
  //  b  = Rechte Seite b
  //  x  = Startvektor fuer die gesuchte Loesung
    
//  Parameter nur fuer vollbesetzte Matrix A:
//    a   = Matrix A (n x lda)
//    lda = fuehrende Dimension der Matrix A

//  Parameter nur fuer duennbesetzte Matrix A:
//    a   = Vektor Nichtnullelemente von A, in CRS-Format gespeichert.
//    ir  = "row-pointer" Zeiger auf den Beginn der i-ten Zeile. Vektor der dimension n+1.
//    ic  = "col-index" Splatenindizes der entsprechenden Wert in A-Vektor.


// SOR-Verfahrens

int nmrc_sor(const struct leqSys *ls, int n, const double *a, int lda, const double *b, double *x)
{
	int i, k, err = 0;
	double *r, *v, res, tol;
	
	/* Validierung der Argumente */
	if ((err = testIvArgs(n, a, ls->kmax, ls->rtol, b, x)))
		return err;
	
	if (lda < n) return -3;
	
	r = nmrc_alloc(2*n, sizeof(*r));
	v = r + n;
	
	/* Kehrwerte der Diagonalelemente */
	if ((err = denseInvDiag(n, a, lda, v))) {
		nmrc_free(r); return -1;	/* A ist singulaer */
	}
	/* Skalierung des Schrittweitenvektors */
	if (ls->relax != 1.0) nmrc_vscale(n, ls->relax, v, v);
	
	/* Residuum am Startwert */
	nmrc_res(n, n, a, lda, x, b, r);
	tol = ls->rtol * (res = nmrc_vnorm(n, r, 1, NmrcNormL2));
	
	if (ls->pout) ls->pout(n, x, res, 0);
	
	/* Iterative Berechnung des Loesungsvektors */ 
	k = 0;
	while (res > tol) {
		/* keine Konvergenz */
		if ( k >= ls->kmax ) { k = -2; break; }
		
		/* Neues x berechnen */
		for (i=0; i<n; i++) {
			res = b[i] - nmrc_dot(n, a+i*lda, 1, x, 1);
			x[i] += res * v[i];
		}
		/* Residuum mit endgueltigem x */
		nmrc_res(n, n, a, lda, x, b, r);
		res = nmrc_vnorm(n, r, 1, NmrcNormL2);
		++k;
		
		if (ls->pout) ls->pout(n, x, res, k);
	}
	nmrc_free(r);
	return k;
}

//Gauss-Seidel-Verfahrens
int nmrc_gs(const struct leqSys *ls, int n, const double *a, int lda, const double *b, double *x)
{
	if (ls->relax != 1.0) return -6;
	return nmrc_sor(ls, n, a, lda, b, x);
}

//Jacobi-Verfahrens

int nmrc_jac(const struct leqSys *ls, int n, const double *a, int lda, const double *b, double *x)
{
	int	i, k, err = 0;
	double	*v, *r, res, tol;
	
	/* Validierung der Argumente */
	if ((err = testIvArgs(n, a, ls->kmax, ls->rtol, b, x)))
		return err;
	
	if (lda < n) return -3;

	r = nmrc_alloc(2*n, sizeof(*r));
	v = r + n;
	
	/* Kehrwerte der Diagonalelemente */
	if ((err = denseInvDiag(n, a, lda, v))) {
		nmrc_free(r); return -1;	/* A ist singulaer */
	}
	/* Skalierung des Schrittweitenvektors */
	nmrc_vscale(n, ls->relax, v, v);
	
	/* Residuum am Startwert */
	nmrc_res(n, n, a, lda, x, b, r);	/* r = b - A*x */
	tol = ls->rtol * (res = nmrc_vnorm(n, r, 1, NmrcNormL2));
	
	if (ls->pout) ls->pout(n, x, res, 0);
	
	/* Iterative Berechnung des Loesungsvektors */ 
	k = 0;
	while (res > tol) {
		/* keine Konvergenz */
		if ( k >= ls->kmax ) { k = -2; break; }
		
		/* Neues x berechnen */
		for (i=0; i<n; i++)
			x[i] = x[i] + r[i] * v[i];
		
		/* Residuum mit endgueltigem x */
		nmrc_res(n, n, a, lda, x, b, r);
		res = nmrc_vnorm(n, r, 1, NmrcNormL2);
		++k;
		
		if (ls->pout) ls->pout(n, x, res, k);
	}
	nmrc_free(r);
	return k;
}

//acobi-Verfahrens
int nmrc_jacs(const struct leqSys *ls, int n, const double *a, const int *ir, const int *ic, const double *b, double *x)
{
	int	i, k, err = 0;
	double	*xn, *r, *work, *v, res, tol;
	
	/* Validierung der Argumente */
	if ((err = testIvArgs(n, a, ls->kmax, ls->rtol, b, x)))
		return err;
	
	if (ir == NULL || ic == NULL) return -4;
	
	work = nmrc_alloc(4*n, sizeof(*work));
	r = work + n; xn = r + n; v = xn + n;
	
	/* Kehrwerte der Diagonalelemente */
	if ((err = crsInvDiag(n, a, ir, ic, v))) {
		nmrc_free(work); return -1;
	}
	/* Skalierung des Schrittweitenvektors */
	nmrc_vscale(n, ls->relax, v, v);
	
	/* Residuum am Startwert */
	nmrc_crsmvp(n, a, ir, ic, x, r);
	nmrc_vsub(n, b, r, r);
	res = nmrc_vnorm(n, r, 1, NmrcNormL2);
	
	if (ls->pout) ls->pout(n, x, res, 0);
	
	/* Toleranz als Vielfaches des Anfangsresiduums */
	tol = ls->rtol * res;
	
	/* Iterative Berechnung des Loesungsvektors */
	k = 0;
	while (res > tol) {
		/* keine Konvergenz */
		if (k >= ls->kmax) { k = -2; break; }
		
		/* neues x berechene */
		for (i=0; i<n; i++)
			x[i] = x[i] + r[i] * v[i];
		++k;
		
		/* Residuum r = b - A*x */
		nmrc_crsmvp(n, a, ir, ic, x, r);
		nmrc_vsub(n, b, r, r);
		res = nmrc_vnorm(n, r, 1, NmrcNormL2);
		
		/* Ausgabe der Werte am Ende des Iterationsschrittes */
		if (ls->pout) ls->pout(n, x, res, k);
	}
	nmrc_free(work);
	return k; /* Konvergenz erreicht? */
}

//SOR-Verfahrens
int nmrc_sors(const struct leqSys *ls, int n, const double *a, const int *ir, const int *ic, const double *b, double *x)
{
	int i, k, err = 0;
	double *r, *v, res, tol;
	
	/* Validierung der Argumente */
	if ((err = testIvArgs(n, a, ls->kmax, ls->rtol, b, x)))
		return err;
	
	if (ir == NULL || ic == NULL) return -4;
	
	r = nmrc_alloc(2*n, sizeof(*r));
	v = r + n;
	
	/* Kehrwerte der Diagonalelemente */
	if ((err = crsInvDiag(n, a, ir, ic, v))) {
		nmrc_free(r); return err;
	}
	/* Skalierung des Schrittweitenvektors */
	if (ls->relax != 1.0) nmrc_vscale(n, ls->relax, v, v);
	
	
	/* Residuum am Startwert */
	nmrc_crsmvp(n, a, ir, ic, x, r);
	nmrc_vsub(n, b, r, r);
	res = nmrc_vnorm(n, r, 1, NmrcNormL2);
	
	if (ls->pout) ls->pout(n, x, res, 0);
	
	tol = ls->rtol * res; /* tol <= res(i)/res(0) */
	
	/* Iterative Berechnung des Loesungsvektors */ 
	k = 0;
	while (res > tol) {
		int	l;
		
		/* keine Konvergenz */
		if (k >= ls->kmax) { k = -2; break; }
		
		/* Berechnen der neuen X-Werte */
		for (i=0, l=ir[0]; i<n; i++) {
			int	ln = ir[i+1];
			res = b[i] - nmrc_crsdot(ln-l, a+l, ic+l, x);
			x[i] += res * v[i];
			l = ln;
		}
		/* Residuum am Ende des Iterationsschrittes */
		nmrc_crsmvp(n, a, ir, ic, x, r);
		nmrc_vsub(n, b, r, r);
		res = nmrc_vnorm(n, r, 1, NmrcNormL2);
		k++;
		
		if (ls->pout) ls->pout (n, x, res, k);
	}
	nmrc_free(r);
	return k;
}

//Gauss-Seidel-Verfahren
int nmrc_gss(const struct leqSys *ls, int n, const double *a, const int *ir, const int *ic, const double *b, double *x)
{
	if (ls->relax != 1.0) return -6;
	return nmrc_sors(ls, n, a, ir, ic, b, x);
}

/* ----------------------------------------------------dense/sparse gradient linear system solvers--------------------------------------------------------------------------- */

 //   ls = Parameter fuer lineares System
 //     .kmax   = Maximale Anzahl von Iterationen
 //     .rtol   = gewuenschte Genauigkeit Norm2(r_k)/Norm2(r_0)
 //     .relax  = Relaxationsfaktor (0, 2]
 //   n  = Zeilenanzahl der Koeffizientenmatrix
  //  b  = Rechte Seite b
  //  x  = Startvektor fuer die gesuchte Loesung
    
//  Parameter nur fuer vollbesetzte Matrix A:
//    a   = Matrix A (n x lda)
//    lda = fuehrende Dimension der Matrix A

//  Parameter nur fuer duennbesetzte Matrix A:
//    a   = Vektor Nichtnullelemente von A, in CRS-Format gespeichert.
//    ir  = "row-pointer" Zeiger auf den Beginn der i-ten Zeile. Vektor der dimension n+1.
//    ic  = "col-index" Splatenindizes der entsprechenden Wert in A-Vektor.

//Gradientenmethode 
int nmrc_grad(const struct leqSys *ls, int n, const double *a, int lda, const double *b, double *x)
{
	int	np1, i, k, err = 0;
	double	alpha, tol, res, rr;
	double	*q, *r;
	
	/* Validierung der Argumente */
	if ((err = testIvArgs(n, a, ls->kmax, ls->rtol, b, x)))
		return err;
	
	if (lda < n) return -3;
	
	/* Test auf Singularitaet */
	for (i = 0, np1 = lda + 1; i < n; i++)
		if (fabs(a[i*np1]) == ZERO) return -1;
		
	q = nmrc_alloc(2*n, sizeof(*q));
	r = q + n;
	
	nmrc_res(n, n, a, lda, x, b, r);	/* r = b - A * x */
	rr = nmrc_dot(n, r, 1, r, 1);
	res = sqrt(rr);
	tol = ls->rtol * res;
	
	if (ls->pout) ls->pout(n, x, res, 0);
	
	/* Iterative Berechnung des Loesungsvektors */ 
	k = 0;
	while (res > tol) {
		/* keine Konvergenz */
		if (k >= ls->kmax) { k = -2; break; }
		
		nmrc_mvp(n, n, a, lda, r, q);	/* q = A * r */
		alpha = rr / nmrc_dot(n, r, 1, q, 1);	/* alpha = <r',r> / <r',q> */
		
		k++;
		nmrc_axpy(n, alpha, r, x, x);	/* x = x + alpha * p */
		nmrc_res(n, n, a, lda, x, b, r); /* r = b - A * x */
		rr = nmrc_dot(n, r, 1, r, 1);
		res = sqrt(rr);
		
		if (ls->pout) ls->pout(n, x, res, k);
	}
	nmrc_free(q);
	return k;
}

//Konjugierte Gradientenmethode 

int nmrc_cg(const struct leqSys *ls, int n, const double *a, int lda, const double *b, double *x)
{
	int	i, k, np1, err = 0;
	double	rr, rro, res, t, alpha, tol;
	double	*p, *r, *q;
	
	/* Validierung der Argumente */
	if ((err = testIvArgs(n, a, ls->kmax, ls->rtol, b, x)))
		return err;
	
	if (lda < n) return -3;
	
	/* Test auf Singularitaet */
	for (i = 0, np1 = lda + 1; i < n; i++)
		if (fabs(a[i*np1]) == ZERO) return -1;	

	q = nmrc_alloc(3*n, sizeof(*q));
	r = q + n; p = r + n;
	
	/* Initialisierung */
	nmrc_res(n, n, a, lda, x, b, r);	/* r0 = b - A*x0 */
	nmrc_vcopy(n, r, 1, p, 1);	/* s0 = r0 */

	rr  = rro = nmrc_dot(n, r, 1, r, 1);		/*  <r0' , r0> */
	tol = ls->rtol * (res = sqrt(rr));	/*  tol = rtol * Norm2(r0) */
	
	if (ls->pout) ls->pout(n, x, res, 0);
	
	/* Iterative Berechnung des Loesungsvektors */ 
	k = 0;
	while (res > tol) {
		/* keine Konvergenz */
		if (k >= ls->kmax) { k = -2; break; }
		
		nmrc_mvp(n, n, a, lda, p, q);	/* q = A * p */
		alpha = rr / nmrc_dot(n, p, 1, q, 1);	/* alpha = <r',r> / <p',q> */
		
		nmrc_axpy(n, alpha, p, x, x);		/* x = x + alpha * p */
		
		nmrc_axpy(n, -alpha, q, r, r);		/* r = r - alpha * q */
		
		rro = rr;
		rr = nmrc_dot(n, r, 1, r, 1) ;	/* rr = (r ' * r ) */
		t = rr / rro;			/* rr / rr_old */
		res = sqrt(rr);
		
		nmrc_axpy(n, t, p, r, p);	/* p = r + t * p */
		
		++k;
		
		if (ls->pout) ls->pout(n, x, res, k);
	}
	nmrc_free(q);
	return k;
}

// Konjugierte Gradientenmethode mit Vorkonditionierung
int nmrc_pcg(const struct leqSys *ls, int n, const double *a, int lda, const double *b, double *x)
{
	int	np1, i, k, err = 0;
	double	rz, res, tol, beta, alpha;
	double	*q, *r, *p, *z;
	void	(*pcond)(), *ppar;
		
	/* Validierung der Argumente */
	if ((err = testIvArgs(n, a, ls->kmax, ls->rtol, b, x)))
		return err;
	
	ppar = ls->ppar;
	if (!(pcond = ls->prec)) { pcond = nmrc_pc_ident; ppar = &n; }
	
	if (lda < n) return -3;
		
	/* Test auf Singularitaet */
	for (i = 0, np1 = lda + 1; i < n; i++)
		if (fabs(a[i*np1]) == ZERO) return -1;
		
	q = nmrc_alloc(4*n, sizeof(*q));
	p = q + n;	/* Suchrichtung */
	r = p + n;	/* Residuum */
	z = r + n;	/* Arbeitsvektor */
	
	nmrc_res(n, n, a, lda, x, b, r); /*  r0 = b - A*x0 */
	tol  = ls->rtol * (res = nmrc_vnorm(n, r, 1, NmrcNormL2));
	beta = 0.0;			/* startwert fuer 1/rz_old */
	
	if (ls->pout) ls->pout(n, x, res, 0);
	
	/* Iterative Berechnung des Loesungsvektors */ 
	k = 0;
	while (res > tol) {
		/* keine Konvergenz */
		if (k >= ls->kmax) { k = -2; break; }
		
		pcond(ppar, r, z, ls->relax);
		rz = nmrc_dot(n, r, 1, z, 1);		/* rz = r' * z */
		
		beta *= rz;			/* beta = <r',z> / <r',z>old */
		nmrc_axpy(n, beta, p, z, p);	/* p = z + beta * p (beta = 0 fuer k = 0) */
		
		nmrc_mvp(n, n, a, lda, p, q);	/* q = A * p */
		alpha = rz / nmrc_dot(n, p, 1, q, 1);	/* <r',z> / <p',q> */ 
		nmrc_axpy(n, alpha, p, x, x);	/* x = x + alpha * p */ 
		nmrc_axpy(n, -alpha, q, r, r);	/* r = r - alpha * q */ 
		
		res = nmrc_vnorm(n, r, 1, NmrcNormL2);
		beta = 1.0/rz;			/* beta = 1 / <r',z>old */
		++k;
		
		if (ls->pout) ls->pout(n, x, res, k);
	}
	nmrc_free(q);
	return k;
}

//Konjugierte Gradientenmethode

int nmrc_cgs(const struct leqSys *ls, int n, const double *a, const int *ir, const int *ic, const double *b, double *x)
{
	int	k, err = 0;
	double	rr, rro, t, alpha, res, tol;
	double	*p, *r, *q;
	
	/* Test the input arguments */
	if ((err = testIvArgs(n, a, ls->kmax, ls->rtol, b, x)))
		return err;
	
	if (ir == NULL || ic == NULL) return -3;
	
	q = nmrc_alloc(3*n, sizeof(*q));
	r = q + n; p = r + n;
	
	if ((err = crsInvDiag(n, a, ir, ic, q))) {
		nmrc_free(q); return -1;	/* A ist singulaer */
	}
	
	/* Initialisierung */
	nmrc_crsmvp(n, a, ir, ic, x, r);
	nmrc_vsub(n, b, r, r);			/* r = b - A * x */
	rr = nmrc_dot(n, r, 1, r, 1);		/* <r', r> */
	nmrc_vcopy(n, r, 1, p, 1);		/* s0 = r0 */
	tol = ls->rtol * (res = sqrt(rr));	/* tol = rtol * Norm2(r) */
	
	if (ls->pout) ls->pout(n, x, res, 0);
	
	/* Iterative Berechnung des Loesungsvektors */ 
	k = 0;
	while (res > tol) {
		/* keine Konvergenz */
		if (k >= ls->kmax) { k = -2; break; }
		
		nmrc_crsmvp(n, a, ir, ic, p, q);	/* q = A * p */
		alpha = rr / nmrc_dot(n, p, 1, q, 1);	/* alpha = <r',r> / <p',q> */
		
		nmrc_axpy(n,  alpha, p, x, x);	/* x = x + alpha * p */
		nmrc_axpy(n, -alpha, q, r, r);	/* r = r - alpha * q */
		
		rro = rr;
		rr = nmrc_dot(n, r, 1, r, 1);	/* rr = (r ' * r ) */
		t = rr / rro;			/* rr / rr_old */
		res = sqrt(rr);
		
		nmrc_axpy(n, t, p, r, p);	/* p = r + t * p */
		
		++k;
		
		if (ls->pout) ls->pout(n, x, res, k);
	}
	nmrc_free(q);
	return k;
}

//Konjugierte Gradientenmethode mit Vorkonditionierung

int nmrc_pcgs(const struct leqSys *ls, int n, const double *a, const int *ir, const int *ic, const double *b, double *x)
{
	int	k, err = 0;
	double	rz, res, tol, beta, alpha;
	double	*q, *r, *p, *z;
	void	(*pcond)(), *ppar;
	
	/* Validierung der Argumente */
	if ((err = testIvArgs(n, a, ls->kmax, ls->rtol, b, x)))
		return err;
	
	if (ir == NULL || ic == NULL) return -3;
	
	ppar = ls->ppar;
	if (!(pcond = ls->prec)) { pcond = nmrc_pc_ident; ppar = &n; }
	
	q = nmrc_alloc(4*n, sizeof(*q));
	p = q + n;	/* Search direction */
	r = p + n;	/* Residual */
	z = r + n;	/* Temporary vector */
	
	/* Suchrichtung: Kehrwert der Diagonalelemente */
	if ((err = crsInvDiag(n, a, ir, ic, q))) {
		nmrc_free(q); return -1;	/* A ist singulaer */
	}
	
	nmrc_crsmvp(n, a, ir, ic, x, r);
	nmrc_vsub(n, b, r, r);		/* r = b - A * x */
	tol  = ls->rtol * (res = nmrc_vnorm(n, r, 1, NmrcNormL2));
	beta = 0.0;			/* startwert fuer 1 / rz_old */
	
	if (ls->pout) ls->pout(n, x, res, 0);
	
	/* Iterative Berechnung des Loesungsvektors */ 
	k = 0;
	while (res > tol) {
		/* keine Konvergenz */
		if (k >= ls->kmax) { k = -2; break; }
		
		pcond(ppar, r, z, ls->relax);		/* setup for z */
		rz = nmrc_dot(n, r, 1, z, 1);		/* rz = r' * z */
		
		beta *= rz;				/* beta = <r',z>neu / <r',z>alt */
		nmrc_axpy(n, beta, p, z, p);		/* p = z + beta * p */
		
		nmrc_crsmvp(n, a, ir, ic, p, q);	/* q = A * p */
		alpha = rz / nmrc_dot(n, p, 1, q, 1);	/* <r',z> / <p',q> */ 
		nmrc_axpy(n,  alpha, p, x, x);		/* x = x + alpha * p */ 
		nmrc_axpy(n, -alpha, q, r, r);		/* r = r - alpha * q */ 
		
		res = nmrc_vnorm(n, r, 1, NmrcNormL2);
		beta = 1.0 / rz;
		++k;
		
		if (ls->pout) ls->pout(n, x, res, k);
	}
	nmrc_free(q);
	return k;
}


/* create sparse from dense matrix */

//    nnz = non-null(A)
//    n   = Dimension von A
//    ir  = "row-pointer"- oder "row_index" - Vektor.(CRS/COO)
//   ic  = "col-index" - Vektor
//    val = "value" - Vektor

 /* 
   Fuer die 3 Vektoren wird Speicher allokiert
   Rueckgabewert ist die Elementanzahl von (a, ic) 
 */

int nmrc_crs(int m, int n, const double *a, int lda, int **ir, int **ic, double **val)
{
	int k, i, j, nnz, *ptr_r, *idx_c;
	double *va, t;
	
	if (!ir ||  !ic) {
		errno = EFAULT; return -1;
	}
	
	/* count the non zero elements*/
	for(i = 0, nnz = 0; i < m ; i++)
		for(j = 0; j < n ; j++)
			if(a[i*lda + j] != 0.0) nnz++;
	
	/* storage in compressed sparse row format */
	if (!(va = nmrc_alloc(nnz, sizeof(*va)))) {
		return -1;
	}
	if (!(ptr_r = nmrc_alloc(m+1, sizeof(*ptr_r)))) {
		nmrc_free(va); return -1;
	}
	if (!(idx_c = nmrc_alloc(nnz, sizeof(*idx_c)))) {
		nmrc_free(va); nmrc_free(ptr_r); return -1;
	}
	
#ifdef _OPENMP
# pragma omp parallel for
#endif
	for(i = 0, k = 0; i < m ; i++) {
		ptr_r[i] = k;
		for(j = 0; j < n ; j++) {
		 	t = a[i*lda + j];
			if(t != 0.0) {
				va[k] = t;
				idx_c[k++] = j;
			}
		}
	}
	ptr_r[n] = nnz;
	*val = va;
	*ir = ptr_r;
	*ic = idx_c;
	
	return nnz;
}

/*---------------------------------------------------------------------- initialize solver parameter--------------------------------------------------------------------- */

void nmrc_leq_init(struct leqSys *ls)
{
	ls->kmax	= 500;	/* max. iterations */
	ls->relax	= 1.0;	/* relaxation factor */
	ls->pout	= 0;	/* iteration data output */
	ls->prec	= 0;	/* preconditioner */
	ls->ppar	= 0;	/* preconditioner parameter */
	ls->rtol	= 1e-7;	/* relative residual resired */
}

/*------------------------------------------------------------------------- vektor operations---------------------------------------------------------------------------- */



 /* Vektor-Initialisierung
  
  Belegen der Elemente eines Vektors mit einem Skalar
  
  n	Anzahl der Elemente
  sk	Skalar
  v	Vektor (output)
 */

void nmrc_vset(int n, double sk, double *v)
{
	int	i;
	
	if (!v) nmrc_die(EFAULT, __func__, "Nullzeiger im Argument");
	
	for (i = 0 ; i < n; i++)
		v[i] = sk;
}


/* Vektor-Skalierung
 * 
 * Skalierung eines Vektors durch Multiplikation mit einem Skalar. 
 
 * 
 	n	Dimension der Vektoren
 	s	Skalierungsfaktor
 	x	Eingabevektor (input)
 	y	Ausgebavektor (input)
 */

void nmrc_vscale(int n, double s, const double *x, double *y)
{
	int	i;
	
	if ( !x || !y )
		nmrc_die(EFAULT, __func__, "Nullzeiger im Argument");
	
	if ( n < 0 )
		nmrc_die(EDOM, __func__, "Negative Dimension");
	
	for (i = 0; i < n; i++)
		y[i] = s * x[i];
}


 /* brief Kopieren eines Vektors
 * 
 * Kopieren der Elemente des Vektors x auf den Vektor y 
 * 
	Anzahl der Elemente
 x	ursprünglicher Vektor (input)
 ldx	Abstand der Elemete in x
 y	Zielvektor (output)
 ldy	Abstand der Elemete in y
 */

void nmrc_vcopy(int n, const double *x, int ldx, double *y, int ldy)
{
	int	i;
	
	if (!x || !y)
		nmrc_die(EFAULT, __func__, "Nullzeiger im Argument");
	
	for (i = 0; i < n; i += 1, x += ldx, y += ldy)
		*y = *x;
}


 /* Addieren eines skalierten mit einem unskalierten Vektor

 n	Länge der Vektoren
 a	Faktor zur Skalierung von x
 x	zu skalierender Vektor x
 y	zu addierender Vektor
 r	Ergebnisvektor
*/

void nmrc_axpy(int n, double a, const double *x, const double *y, double *r)
{
	int i, m = n % 4;
	
	if (!x || !y)
		nmrc_die(EFAULT, __func__, "Nullzeiger im Argument");
	
	if (n < 0)
		nmrc_die(EDOM, __func__, "negative Dimensionn");
	
	if (m != 0) {
		for (i = 0; i < m; i++)
			r[i] = y[i] + a * x[i];
		
		x += m;
		y += m;
		r += m;
	}
	for (i = m; i < n; i += 4, x += 4, y += 4, r += 4) {
		r[0] = y[0] + a * x[0];
		r[1] = y[1] + a * x[1];
		r[2] = y[2] + a * x[2];
		r[3] = y[3] + a * x[3];
	}
}



/* Vektor-Subtraktion 
  n	Dimension der Vektoren
  x	Erster vektor (input)
  y	Zweiter Vektor (input)
  z	Ergebnis der Subtraktion (output)
 */

void nmrc_vsub(int n, const double *x, const double *y, double *z)
{
	int	i;
	
	if ( !x || !y )
		nmrc_die(EFAULT, __func__, "Nullzeiger im Argument");
	
	if ( n < 0 )
		nmrc_die(EDOM, __func__, "Negative Dimension");
	
	for (i = 0; i < n; i++)
		z[i] = x[i] - y[i];
}


/* -----------------------------------------------------------dense preconditioner-------------------------------------------------- */


 /*Residuen als Korrekturwert
 * \param par	Parameter für Vorkonditionierer
 * \param r	aktuelles Residuum
 * \param z	neue Korrektur
 * \param omega	Relaxationsfaktor \c ppar.relax
 */


void nmrc_pc_ident(const int *par, const double *r, double *z, double omega)
{
	nmrc_vcopy(*par, r, 1, z, 1);	/* z = r */
}

/*---------------------------------------------------------- heap memory management --------------------------------------------------*/


 /* Die Funktion allokiert Speicherplatz für einen Vektor
 
  Aufruf
   \code vec_alloc(&v, sizeof(*v), dim); \endcode
 
 rows		Dimension des Vektors
 elSize	Größe des Datentyps in Byte
 pA		Zeiger auf den allokierten Speicherbereich (O)
 */




void vec_alloc(void *pA, size_t elSize, int rows)
{
	size_t	rowBytes;
	char	*a;
	
	if ( rows < 1 ) {
		fprintf(stderr, "Fehlerhafte Vektordimension (rows = %d)\n", rows);
		exit (EXIT_FAILURE);
	}
	
	rowBytes = rows * elSize;
	
	if ( (a = malloc(rowBytes)) == NULL ) {
		fputs("Fehler bei Anforderung vom Speicherplatz\n", stderr);
		exit (EXIT_FAILURE);
	}
	
	*(char **)pA = a;
	
}	/* vec_alloc */


 /* Die Funktion allokiert Speicherplatz für eine
 \f$ n \times m \f$ Matrix und deren Zeilenzeiger.
  

double **a = (double **) nmrc_matalloc(n, m, sizeof(**a));

  
  rows	Anzahl der Zeilen
  cols	Anzahl der Spalten
  esze	Größe des Datentyps in Byte
 */


void mat_alloc(void *pA, size_t elSize, int rows, int cols)
{
	int	i;
	size_t	rowBytes;
	char	*p, **a;
	
	if ( rows < 1 || cols < 1 ) {
		fprintf(stderr, "Fehlerhafte Matrixdimensionen (rows = %d, cols = %d)\n", rows, cols);
		exit (EXIT_FAILURE);
	}
	
	rowBytes = cols * elSize;
	
	/* Zeiger auf Feld von Zeigern ( Zeilen ) und Zeilendaten */
	if ( (a = malloc(rows * (sizeof(*a) + rowBytes))) == NULL ) {
		fputs("Fehler bei Anforderung vom Speicherplatz\n", stderr);
		exit (EXIT_FAILURE);
	}
	
	/* Verweis auf Indexfeld wird gespeichert */
	*(char ***)pA = a;
	
	/* Feld beginnt nach letzter Referenz */
	for (i = 0, p = (char *)(a+rows); i < rows; i++, p += rowBytes)
		a[i] = p;
	memset(a[0], 0, rows*rowBytes);
}	/*mat_alloc*/


void vec_free(void *pV)
{
	if ( pV != NULL ) free(pV);
}	/* vec_free */


  
/*Die Funktion allokiert Speicherplatz für Datenelemente
double *a = nmrc_alloc(5, sizeof(*a));
len	Anzahl der Elemente
 esze	Größe des Datentyps in Byte
 */

extern void *nmrc_alloc(int len, size_t esze)
{
	if (len < 0 || !esze) {
		errno = EINVAL; return 0;
	}
	return calloc(len, esze);
}



 /* Freigabe von Speicherplatz der über
 ::nmrc_alloc oder ::nmrc_matalloc bereitgestellt wurde.
  

 double *a = nmrc_alloc(5, sizeof(*a));
 // arbeiten mit a
 nmrc_free(a);
*/

extern void nmrc_free(void *ptr)
{
	free(ptr);
}

/*--------------------------------------------------------- BLAS-like vector operations -----------------------------------------------------*/


 
 /* Berechnung des Skalarprodukts \f$ x \cdot y \f$
n	Länge der Vektoren
x	erster Vektor
ldx	führende Dimension des ersten Vektors
y	zweiter Vektor
ldy	führende Dimension des zweiten Vektors
 * \return Ergebnis des Skalarproduktes*/

double nmrc_dot(int n, const double *x, int ldx, const double *y, int ldy)
{
	int	i;
	double	d;
	
	if (!x || !y) {
		nmrc_die(EFAULT, __func__, "Nullzeiger im Argument");
	}
	for (i = 0, d = 0.0; i < n; i++)
		d += x[i*ldx] * y[i*ldy];
	
	return d;
}




/*------------------------------------------------------------------------------ norms --------------------------------------------------*/




 // n	Anzahl der Zeilen
 // m	Anzahl der Spalten
 // a	Matrix
 // lda	Führende Dimension der Matrix A
 // type	Art der Norm
 
 // \return 	die gewählte Matrixnorm

//norm = nmrc_mnorm(n, m, a[0], lda, NmrcNormRS);
//norm = nmrc_mnorm(n, m, a[0], lda, NmrcNormCS);
//norm = nmrc_mnorm(n, m, a[0], lda, NmrcNormFrob);




double nmrc_mnorm(int n, int m, const double *a, int lda, enum NmrcNorm type)
{
	int	i;
	double	en;
	
	if (!a) {
		errno = EFAULT; return -1;
	}
	if (n < 1 || !InRange(m,0,lda)) {
		errno = EDOM; return -2;
	}
	
	switch (type) {
	  case NmrcNormFrob:
		for (i = 0, en = 0.0; i < n; i++, a += lda)
			en += nmrc_dot(m, a, 1, a, 1);
		return sqrt(en);
	  case NmrcNormRS:
		for (i = 0, en = 0.0; i < n; i++, a += lda) {
			double s = nmrc_vnorm(m, a, 1, NmrcNormL1);
			if (s > en) en = s;
		}
		return en;
	  case NmrcNormCS:
		for (i = 0, en = 0.0; i < n; i++, ++a) {
			double s = nmrc_vnorm(m, a, lda, NmrcNormL1);
			if (s > en) en = s;
		}
		return en;
	  default:
	//	errno = EBADRQC;
		return -2;
	}
}

// n	Anzahl der Elemente
// v	Vektor
// ld	Abstand der Vektorelemente 
// type	Art Norm
  
// return 	die gewählte Vektornorm


//en = nmrc_vnorm(n, v, 1, NmrcNormInf); // für Maximumnorm
//en = nmrc_vnorm(n, v, 1, NmrcNormL2);  // für Euklidische Norm



double nmrc_vnorm(int n, const double *v, int ld, enum NmrcNorm type)
{
	int i;
	double en = 0.0;
	
	if (!v) {
		errno = EFAULT; return -3;
	}
	if (n <= 0) {
		errno = ERANGE; return -2;
	}
	
	switch (type) {
	  case NmrcNormL2:
		en = nmrc_dot(n, v, ld, v, ld);
		return sqrt(en);
	  case NmrcNormL1:
		for (i = 0, en = 0.0; i < n; i++, v += ld)
			en += fabs(*v);
		return en;
	  case NmrcNormInf:
		for (i = 0, en = 0; i < n; i++, v += ld) {
			double val = fabs(*v);
			if (en < val) en = val;
		}
		return en;
	  default:
//		errno = EBADRQC;
		return -1;
	}
}





// n	Anzahl der Elemente
// A	Matrix
// lda	Führende Dimension der Matrix
// type	Art der Norm
  
//return die Konditionszahl zur gewählten Matrixnorm
 
//en = nmrc_condA(n,  A, lda, NmrcNormRS);    // für Maximumnorm
//en = nmrc_condA(n, *A, lda, NmrcNormFrob);  // für Euklidische Norm

double nmrc_condA(int n, const double *A, int lda, enum NmrcNorm norm)
{
	double	*iA, ka, kai; /* Inverse Matrix + Konditionszahl*/
	
	ka = nmrc_mnorm(n, n, A, lda, norm);
	
	if (ka < 0.0) return NAN;
	
	vec_alloc(&iA, sizeof(*iA), n*n);
	if (nmrc_cinverse(n, A, lda, iA, n) < 0) ka = NAN;
	else if ((kai = nmrc_mnorm(n, n, iA, n, norm)) < 0) ka = NAN;
	else ka *= kai;
	
	vec_free(iA);
	
	return ka;
}


/* -------------------------------------------------------------------- inverse matrix -------------------------------------------*/
/*
 * Das Unterprogramm berechnet die inverse Matrix \f$ A^{-1}\f$ von A.
 * Die Koeffizienten von A werden überschrieben.
  
 	n	Dimension der Matrix A
 	lda	führende Dimension der Matrix A lda >= n
 	a	Koeffizientenmatrix A
 	a	Inverse von A
 * */


int nmrc_inverse(int n, double *a, int lda)
{
	int	*pv, j, ldai = n, err;	
	double	*ai, *b;
	
	if (n < 2 || lda < n) return -2;
	
	vec_alloc(&pv, sizeof(*pv), n);
	
	if ((err = nmrc_ludec(n, a, lda, pv)) < 0 ) {
		vec_free(pv);
		return err;
	}
	
	vec_alloc(&ai, sizeof(*ai), n*n + n);
	b = ai + n*n;
	
	for (j = 0; j < n; j++) {
	
		/* Belegung der rechten Seite */
		nmrc_vset(n, 0.0, b);
		b[j] = 1.0;
		
		/* Berechnung der Vektoren der inversen Matrix */
		nmrc_lusol(n, a, lda, pv, b, ai + j*n);
	}
	
	nmrc_trans(n, n, ai, ldai, a, lda);
	
	vec_free(ai);
	vec_free(pv);
	
	return 0;
}



/* Inverse
 * 
 Das Unterprogramm berechnet die inverse Matrix \f$ A^{-1}\f$ von A.
 Die Matrix A bleibt unverändert.
  
 n	Dimension der Matrix A
 lda	führende Dimension der Matrix A lda >= n
 a	Koeffizientenmatrix A
 ai	Koeffizientenmatrix der Inverse von A
 ldai	führende Dimension der Inversen
  */

int nmrc_cinverse(int n, const double *a, int lda, double *ai, int ldai)
{
	double *ca;
	int ldca = n, err = 0;
	
	vec_alloc(&ca, sizeof(*ca), n*n);
	nmrc_mcopy(n, n, a, lda, ca, ldca);
	
	if ((err = nmrc_inverse(n, ca, lda)) == 0)
		nmrc_mcopy(n, n, ca, lda, ai, ldai);
	
	vec_free(ca);
	
	return err;
}

/* ---------------------------------------------------------- matrix/vektor operations -------------------------------------------------------- */


 /* Matrix-Vektor-Multiplikation
 *

 * 
  n	Anzahl der Zeilen von A
  m	Anzahl der Spalten von A
  a	Matrix A (input)
  lda	führende Dimension der Matrix A
  x	Vektor x (input)
  y	Ergebnisvektor (output)
 */

void nmrc_mvp(int n, int m, const double *a, int lda, const double *x, double *y)
{
	int	i;
	
	if ( !a || !x || !y )
		nmrc_die(EFAULT, __func__, "Nullzeiger im Argument");
	
	if ( !InRange(m,0,lda) )
		nmrc_die(EDOM, __func__, "Parameter nicht in Grenzen");
#ifdef _OPENMP
# pragma omp parallel for
#endif
	for (i = 0; i < n; i++)
		y[i]= nmrc_dot(m, a + i*lda, 1, x, 1);
}


 /* Residuenvektor eines Linearen Gleichungssystems \n
 \f$ r = b - A \cdot x \f$
  
  n	Anzahl der Zeilen von A bzw. Länge von b
  m	Anzahl der Spalten von A bzw. Länge von x
  a	Matrix A der Dimension (n x m)
  lda	führende Dimension der Matrix A
  x	Vektor x der Dimension (m)
  b	Vektor b der Dimension (n)
  r	Ergebnisvektor der Dimension (n)
*/


void nmrc_res(int n, int m, const double *a, int lda, const double *x, const double *b, double *r)
{
	register int i;
	
	if (!a || !x || !b || !r)
		nmrc_die(EFAULT, __func__, "Nullzeiger im Argument");
	
	if (!InRange(m,0,lda))
		nmrc_die(EDOM, __func__, "Negative Dimension");
#ifdef _OPENMP
# pragma omp parallel for
#endif
	for (i = 0; i < n; i++)
		r[i] = b[i] - nmrc_dot(m, a + i*lda, 1, x, 1);
}


/*--------------------------------------------------------- multiply sparse matrix and vektor (y = Sx) -------------------------------------------*/
 /*  Matrix/Vektor-Produkt: y = A * x
  
  Berechnung des inneren Produkts: d = x' * y
  A ist in val, ir und ic gepackt. x und y sind voller Vektoren der Laenge n. */

void nmrc_crsmvp(int n, const double *val, const int *ir, const int *ic, const double *x, double *y)
{
	int i;
#ifdef _OPENMP
# pragma omp parallel for
#endif
	for (i = 0; i < n; i++) {
		int	j = ir[i];
		y[i] = nmrc_crsdot(ir[i+1]-j, val+j, ic+j, x);
	}
}

/*-------------------------------------------------------- dot product of indexed and dense vector--------------------------------------- */


/* inneres Produkt d = x' * y
 *  Berechnung des inneren Produkts: d = x' * y
 *  x ist in val und idcol gepackt, y ist voll besetzter Vektor. 
*/


double nmrc_crsdot(int lenr, const double *val, const int *idcol, const double *y)
{
	int	i;
	double	d;
	
	for (i = 0, d = 0.; i < lenr; i++)
		d += val[i] * y[idcol[i]];
	return d;
}



/*---------------------------------------------------------------- matrix operations-------------------------------------------------------  */

 /* Transponieren einer rechteckigen Matrix
 * 
 * Ziel und Quelle duerfen keine ueberlappenden Bereiche
 * enthalten oder muessen identisch sein.
 * 
 n	Anzahl der Zeilen von A bzw. Anzahl der Spalten von A^t
 m	Anzahl der Spalten von A bzw. Anzahl der Zeilen von A^t
 a	Matrix A (input)
 lda	führende Dimension der Matrix A
 aT	Transponierte Matrix A^t (output)
 ldt	führende Dimension der transponierten Matrix*/


void nmrc_trans(int n, int m, const double *a, int lda, double *aT, int ldt)
{
	int	i, j, k;
	
	if (!a || !aT)
		nmrc_die(EFAULT, __func__, "Nullzeiger im Argument");
	
	if (!InRange(m,0,lda) || !InRange(n,0,ldt))
		nmrc_die(EDOM, __func__, "Parameter nicht in Grenzen");
	
	k = (n < m) ? n : m;
	
	for (i = 0; i < k; i++, a += lda + 1, aT += ldt + 1) {
		const double *ar = a, *ac = a;
		double	*tr = aT, *tc = aT;
		
		/* transponiere quadratischen Teilbereich */
		for (j = i; j < k; j++, ar += lda, ac++, tr += ldt, tc++) {
			double	tmp = *ar;
			*tr = *ac;
			*tc = tmp;
		}
		
		/* Schleife wird bei m < n aktiv */
		for ( ; j < n; j++, ar += lda, tc++)
			*tc = *ar;
		
		/* Schleife wird bei n < m akiv */
		for ( ; j < m; j++, ac++, tr += ldt)
			*tr = *ac;
	}
}


/* brief Kopieren einer Matrix
 * 
 * Kopieren der Koeffizienten der Matrix A auf die Matrix B \n
 * \f$ B_{ij} = A_{ij} \f$
 * 
  n	Anzahl der Zeilen von A bzw. von B
  m	Anzahl der Spalten von A bzw. von B
  a	Matrix A (input)
  lda	führende Dimension der Matrix A
  b	Die Kopie von A, also Matrix B (output)
  ldb	führende Dimension der Matrix B
 */

void nmrc_mcopy(int n, int m, const double *a, int lda, double *b, int ldb)
{
	int	i;
	
	if (!a || !b)
		nmrc_die(EFAULT, __func__, "Nullzeiger im Argument");
	
	if (!InRange(m,0,lda) || !InRange(m,0,ldb))
		nmrc_die(EDOM, __func__, "Dimension nicht in Grenzen");
#ifdef _OPENMP
# pragma omp parallel for
#endif
	for (i = 0; i < n; i++)
		(void) memcpy(b + i*ldb, a + i*lda, m*sizeof(*b));
}


/*------------------------------------------------------------------- LAPACK wrapper ---------------------------------------------------*/

 /* Zerlegung einer Matrix \a A in eine untere (linke) Dreiecksmatrix \a L 
 und eine obere (rechte) Dreiecksmatrix \a U, wobei die Diagonalelemente
  von \a L gleich 1 sind.
  
n	Dimension der Matrix A
a	Koeffizientenmatrix A
a	Koeffizientenmatrix LU
lda	führende Dimension der Matrix A
pv	Permutationsvektor

 */
int nmrc_ludec(int n, double *a, int lda, int *pv)
{
	int	k, j, i, p, print_flg;
	double	*v, sum, t, big;
	
	if (n < 2 ||  lda < n) return -2;
	if (!a || !pv) return -3;
	
	// print_flg = print_iter_flag & PRINT_ITER_LS;
	vec_alloc(&v, sizeof(*v), n);
	
	/* Initialisierung des Permutationsvektors */
	if (pv) {
		for (i = 0; i < n; i++)
			pv[i] = i;
	}
	
	/* LU-Zerlegung */
	for (i = 0; i < n; i++) {
		for (k = i; k < n; k++) {
			for (j = 0, sum = 0.0; j < i; j++)
				sum += a[k*lda + j] * a[j*lda + i];
			v[k] = a[k*lda + i] - sum;
		}
		big = fabs(v[i]);
		p = i;
		/* Pivotsuche */
		if (pv) {
			for (k = i+1; k < n; k++) {
				t = fabs(v[k]);
				if (t > big) {
					p = k;
					big = t;
				}
			}
		}
		if (big < TOL) {
			if (big == ZERO) {
				vec_free(v);
				return -1;
			}
			v[p] = sign(v[p]) * TOL;
			fprintf(stderr, "Warning %s: Matrix is close to singular\n", __func__);
		}
		/* Zeilentausch */
		if (p != i) {
			for (k = 0; k < n; k++) {
				t = a[i*lda + k];
				a[i*lda + k] = a[p*lda + k];
				a[p*lda + k] = t;
			}
			k = pv[i];
			pv[i] = pv[p];
			pv[p] = k;
			t = v[i];
			v[i] = v[p];
			v[p] = t;
		}
		/* Aktualisierung der Diagonalelemente */
		t = v[i];
		a[i*lda + i] = t;
		for (k = i+1; k < n; k++) {
			for (j = 0, sum = 0.0; j < i; j++)
				sum += a[i*lda + j] * a[j*lda + k];
			a[i*lda + k] -= sum;
			a[k*lda + i] = v[k] / t;
		}
		if (print_flg) {
			printf("%d. Umformungsschritt:\n", i+1);
			//mat_print(" Matrix A", "% 13g", stdout, n, n, lda, a);
		}
	}
	vec_free(v);
	return 0;
}


void nmrc_lusol(int n, const double *a, int lda, const int *pv, const double *b, double *x)
{
	int	i, j;
	double	t;
	
	/* Reorganisieren der rechten Seite */
	if (!pv) {
		nmrc_vcopy(n, b, 1, x, 1);
	}
	else for (i = 0; i < n; i++) {
		x[i] = b[pv[i]];
	}
	/* Vorwaertssubstitution */
	for (i = 0; i < n; i++) {
		for (j = 0, t = x[i]; j<i; j++)
			t -= a[i*lda + j] * x[j];
		x[i] = t;
	}
	/* Rueckwaertssubstitution */
	for (i = n-1; i >= 0; i--) {
		for (j = i+1, t = x[i]; j < n; j++) {
			t -= a[i*lda + j] * x[j];
		}
		x[i] = t / a[i*lda + i];
	}
}

/*-------------------------------------------------------------------- miscellaneous-------------------------------------------------------- */


extern void nmrc_die(int type, const char *fcn, const char *msg)
{
	errno = type;
	
	fputs(fcn, stderr);
	fputs("(): ", stderr);
	fputs(msg, stderr);
	fputc('\n', stderr);
	
	abort();
}
























