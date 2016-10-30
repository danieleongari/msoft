/*******************************************************************************
File MSOFT.h is a header file for program MSOFT.c.
*******************************************************************************/
#define NX 1024   /* Number of mesh points */
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

/* Function prototypes ********************************************************/
void init_param();
void init_prop();
void init_wavefn();
void single_step(int step);
void pot_prop();
void kin_prop();
void periodic_bc();
void calc_energy();
void four1(double data[], unsigned long nn, int isign);
void create_C1f();
void create_C2f();
void update_C1();
void update_C2();
void calc_ekin();
void calc_epot();
void calc_norm();
void calc_energy(); 
void calc_eigenvalues(int i);
void calc_De_and_Det(int i);
void pop_states();
void print_pop(int step, FILE *f5);
void print_energy(int step, FILE *f1);
void print_pot_ad(FILE *f4);
void print_pot_di(FILE *f6);
void print_wavefn(int step, FILE *f2, FILE *f3);

/* Input parameters ***********************************************************/
double LX;       /* Simulation box length */
double DT;       /* Time discretization unit */ 
double M;	 /* Mass of the system */
int NSTEP;       /* Number of simulation steps */
int NECAL;       /* Interval to calculate energies */
int NNCAL;       /* Interval to calculate norm     */
int pot_type;    /* Potential type */
double X0,S0,P0; /* Center-of-mass, spread & momentum of initial wave packet */
double A,B,C,D;  /* Parameters of potential */

/* Arrays **********************************************************************
C1[NX+2][2]:     C1[i][0|1] is the real|imaginary part of the first component of the 
		 wave function on mesh point i
C2[NX+2][2]:     C2[i][0|1] is the real|imaginary part of the second component of the 
		 wave function on mesh point i
Cf[(NX+1)*2]:    C[2i|2i+1] is the real|imaginary part of the wave function
		 on mesh point i 
T[NX+2]:	 T[i] is the kinetic energy at mesh point i
t[NX+2][2]:	 t[i][] is the kinetic propagator on i (real|imaginary part)
u[NX+2][2][2]:   u[i][j][] is the jth component of the diagonal potential propagator on i (real|imaginary part)
h[NX+2][2][2]:   h11[i][j][k] is the element j,k of the matrix h at mesh point i
E[NX+2][2]:      E[i][j] is the j-th eigenvalue of the matrix at mesh point i
De[NX+2][2][2]:  De[i][j][k] is the element j,k of the matrix De, which has the 
		 eigenvectors (of the matrix h) as columns, at mesh point i
Det[NX+2][2][2]: Det[i][j][k] is the element j,k of the matrix Det, which is the 
		 transposed of De, at mesh point i

*******************************************************************************/
double C1[NX+2][2];
double C2[NX+2][2];
double Cf[(NX+2)*2];
double T[NX+2];
double t[NX+2][2];
double u[NX+2][2][2];
double h[NX+2][2][2];
double E[NX+2][2];
double De[NX+2][2][2];
double Det[NX+2][2][2];

/* Variables *******************************************************************
dx   = Mesh spacing
ekin = Kinetic energy
epot = Potential energy
etot = Total energy
P1 = Population of state 1
P2 = Population of state 2
*******************************************************************************/
double dx;
double norm;
double ekin,epot,etot; 
double P1,P2,d1,d2;
