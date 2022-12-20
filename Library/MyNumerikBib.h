//Autoren: Nikos Alexiadis, Jan Brunner, Giuliano Lombardo

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>



# define MACHEPS 1.2e-16
# define MACH_EPSILON 2.2204460492503131E-16
# define EPSFCN 1.0e-15
# define TOL  1.0E-12
# define MIN_STEPSIZE  1.0E-8
# define ZERO  0.0
#define PI 3.142857



# define sign(X)	((X)>= 0 ? 1 :-1)
# define max(a, b)	(((a) > (b)) ? (a) : (b))
# define min(a, b)	(((a) < (b)) ? (a) : (b))
# define InRange(val, min, max)	((min < max) ? ((val >= min) ? (val <= max) : 0) : 0)

//Euler-Verfahren
//Abgabe 1
void euler_explizit_2(int n, double h, double t, double tmax, double *y_i, double *z,void(fcn)(int, double*,double*));
void euler_explizit(int n,int k, double h, double t, double tmax, double *y_i, double *z,void(fcn)(int, double*,double*));
void euler_implizit(int n, double h,double t, double tmax, double *y_i, double *z, void(fcn)(int, double*,double*));
void euler_Schrittweitensteuerung(int n,double h, int *s, double *t, double tmax, double *y_i, double *z, double P,double rho, double tol, double h_min, double n_w, void(fcn)(int, double*,double*));
void euler_semi_implizit(int n,double delta_T, double *t, double tmax,double *y_i, double *z, double *error_frob,double *error_y_j_frob, void(fcn)(int, double*,double*));

//Abgabe 2
void dae_semi_implizit_Index_1(int n, int *j_t,int *z_counter, double delta_T,double *B, double *t, double tmax,double *y_i, double *z, void(fcn)(int, double*,double*));
void dae_si_dense_output(int n,int *j_t,int *z_counter,double delta_T, double *B, double *t, double tmax,double *y_i, double ymax, double *z, void(fcn)(int, double*,double*));
void dae_si(int n,int *j_t,int *z_counter,double delta_T, double *B, double *t, double tmax,double *y_i, double *z, void(fcn)(int, double*,double*));
void hermit_interpolation(int n,int j_t1, int z_counter1, int anzahl,double delta_T, double *t, double *y_t, double *y_t_2, double *z,void(fcn)(int, double*,double*));
void dae_si_Schritt(int n,int *j_t,int *z_counter,double delta_T, double *B, double *t, double tmax,double *y_i, double *z, void(fcn)(int, double*,double*),void(fcn1)(int, double*,double*));



//neue Funktionen
void Richardson_extrapolation(int n, int k_steps, double *h_Werte, double *P, double *y_t_2);
void euler_semi_implizit_Schritt(int n,double delta_T, double *B, double *t,double *y_t, double *y_t_2, double *z, void(fcn)(int, double*,double*));
void euler_SI(int n,int *j_t, int *z_counter,double delta_T, double *B, double *t,double tmax, double *y_i, double *z, void(fcn)(int, double*,double*));
void euler_SI_dense_output(int n,int anzahl,int *j_t, int *z_counter, double delta_T, double *B, double *t,double tmax, double *y_i, double *z, void(fcn)(int, double*,double*));
void euler_SI_DO_abbruch(int n, int *j_t, int *z_counter, double delta_T, double *B, double *t,double tmax, double ymax, double *y_i, double *z, void(fcn)(int, double*,double*));
void euler_SI_DAE_In_2(int n, int *j_t, int *z_counter, double delta_T, double *B, double *t,double tmax, double ymax, double *y_i, double *z, void(fcn)(int, double*,double*));
void hermit_interpolation_Schritt(int n,int j_t1, int z_counter1, int anzahl,double delta_T, double *t, double *y_t, double *y_t_2, double *z,void(fcn)(int, double*,double*));
void hermit_interpolation_t(int n,int j_t,double ymax,double delta_T, double *t, double *y_t, double *y_t_2, double *y,void(fcn)(int, double*,double*));

//Matrix-Vektor-Operationen
void mat_vec_multiplication(int n, int m, double *A, double *x,double *b);
void matrix_multiplication(double *,double *,int , double *);

// Funktionen zur Numerischen Berechnung der Jacobi-Matrix (fd:Vorwaerts, cd: Central)
void fjac_cd(int n, double h, double *x, double *jac, void (*fcn)(int, double *, double *));

//Newton -Verfahren//
//Abbruchkriterium: Betrag des Residuums
void newton_iter(int n, int *ni, int kmax, double tol, double *x_start, double *residuum,  double *lsg, void(fcn)(int, double*, double*));


//Array kopieren -> kopiert arr 2 auf arr 1 -> arr1 = arr2
void copy_array(double *arr1, double *arr2,int size);
//Array mal minus 1 nehmen
void array_minus(double *arr, int size);

//[size1][size2] Matrix allokieren
void matrix_alloc(double ***A, int size1, int size2);
void matrix(double *A, double *B,int n);


//Vektor betrag
double vec_fabs(double *vec, int size);
//vektor differenz: erg = arr1 - arr2
void vec_diff(double *arr1, double *arr2, double *erg, int size);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//Matrix transponieren
void matrix_trans(int n, int m, double *matrix);

//LR-Zerlegung
void LR_separation(int, double*, int, double*, double*, double*);
void forward_substitution(double*, double*, double*, int);
void back_substitution(double* , double*, double*, int);
//Lu_complete --> Komplette LU-Zerlegung mit Spaltenpivotisierung
void lu_complete(double *A, double *x, double *b, int n);







////
// Funktionen aus der Numerik - Bibliothek des ICVT 
// NMRC-LIBRARY

/* linear systems solver parameter */
struct leqSys {
	int	kmax;
	float	relax;
	void	(*pout)(int n, double *x, double res, int iter);
	void	(*prec)();
	void	*ppar;
	double	rtol;
};

/* linear systems preconditioner parameter */
struct precPar
{
	int	n, lda;
	const double *a, *av;
	const int *ir, *ic;
};

/* Norms */
enum NmrcNorm {
	NmrcNormL2,	/*!< L2-Norm */
	NmrcNormL1,	/*!< L1-Norm */
	NmrcNormInf,	/*!< Infinimum Norm */
	
	NmrcNormFrob,	/*!< Frobenius Norm */
	NmrcNormRS,	/*!< Zeilensummennorm */
	NmrcNormCS	/*!< Spaltensummennorm */
};


/* direct linear systems solver */
extern int  nmrc_gauss(int , double *, int , double *, double *);

/* dense/sparse splitting linear system solvers */
extern int nmrc_gs  (const struct leqSys *, int , const double *, int , const double *, double *);//
extern int nmrc_sor (const struct leqSys *, int , const double *, int , const double *, double *);//
extern int nmrc_jac (const struct leqSys *, int , const double *, int , const double *, double *);//
extern int nmrc_gss (const struct leqSys *, int , const double *, const int *, const int *, const double *, double *);//
extern int nmrc_sors(const struct leqSys *, int , const double *, const int *, const int *, const double *, double *);//
extern int nmrc_jacs(const struct leqSys *, int , const double *, const int *, const int *, const double *, double *);//

/* dense/sparse gradient linear system solvers */
extern int nmrc_grad(const struct leqSys *, int , const double *, int , const double *, double *);  //
extern int nmrc_cg  (const struct leqSys *, int , const double *, int , const double *, double *); //
extern int nmrc_pcg (const struct leqSys *, int , const double *, int , const double *, double *);//
extern int nmrc_cgs (const struct leqSys *, int , const double *, const int *, const int *, const double *, double *);//
extern int nmrc_pcgs(const struct leqSys *, int , const double *, const int *, const int *, const double *, double *);//


/* create sparse from dense matrix */
extern int nmrc_crs(int , int , const double *, int , int **, int **, double **);


/* initialize solver parameter */
extern void nmrc_leq_init(struct leqSys *);

/* vektor operations */
extern void   nmrc_vset  (int , double , double *);
extern void   nmrc_vscale(int , double , const double *, double *);
extern void   nmrc_vcopy (int , const double *, int , double *, int );
extern void   nmrc_axpy  (int , double , const double *, const double *, double *);
extern void   nmrc_vsub  (int , const double *, const double *, double *);

/* dense preconditioner */
extern void   nmrc_pc_ident(const int *,            const double *r, double *z, double relax);

/* heap memory management */
extern void vec_alloc(void *, size_t , int );
extern void mat_alloc(void *, size_t , int , int );
extern void vec_free(void *);
extern void   *nmrc_alloc(int len, size_t size);
extern void   nmrc_free(void *);

/* BLAS-like vector operations */
extern double nmrc_dot  (int , const double *, int , const double *, int);

/* norms */
extern double nmrc_mnorm  (int , int , const double *, int , enum NmrcNorm );
extern double nmrc_vnorm  (int , const double *, int , enum NmrcNorm );
extern double nmrc_condA  (int , const double *, int , enum NmrcNorm );


extern int    nmrc_inverse (int , double *, int );
extern int    nmrc_cinverse(int , const double * , int , double *, int );

/* matrix/vektor operations  */
extern void   nmrc_mvp (int , int , const double *, int , const double *, double *);
extern void   nmrc_res (int , int , const double *, int , const double *, const double *, double *);
/* multiply sparse matrix and vektor (y = Sx) */
extern void   nmrc_crsmvp(int , const double *, const int *, const int *, const double *, double *);

/* dot product of indexed and dense vector */
extern double nmrc_crsdot(int , const double *, const int *, const double *);
/* matrix operations  */
extern void   nmrc_trans(int , int , const double *, int , double *, int );
extern void   nmrc_mcopy(int , int , const double *, int , double *, int );

/* LAPACK wrapper */
extern int    nmrc_ludec(int , double *, int , int *);
extern void   nmrc_lusol(int , const double *, int , const int *, const double *, double *);

/* miscellaneous */
extern void   nmrc_die(int , const char *, const char *);