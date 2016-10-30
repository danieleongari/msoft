#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MSOFT.h"



int main(int argc, char **argv) {

  int step; /* Simulation loop iteration index */

  FILE *f1=fopen("energy.dat","w");
  FILE *f2=fopen("psisq.dat","w");
  FILE *f3=fopen("norm.dat","w");
  FILE *f4=fopen("potential_ad.dat","w");
  FILE *f5=fopen("population.dat","w");
  FILE *f6=fopen("potential_di.dat","w");
  printf("  completed      norm           etot     pop1    pop2\n") ;

  init_param();  /* Initialize input parameters */
  init_prop();   /* Initialize the kinetic & potential propagators */
  init_wavefn(); /* Initialize the electron wave function */
  print_pot_ad(f4);
  print_pot_di(f6);
  print_wavefn(0,f2,f3); /* print initial conditions */

  for (step=1; step<=NSTEP; step++) {
    single_step(step); /* Time propagation for one step, DT */
    if (step==1 || step%NECAL==0) {
      calc_energy();
      print_energy(step,f1);
      pop_states();
      print_pop(step,f5);
    }
    if (step==1 || step%NNCAL==0) {
      calc_norm();
      print_wavefn(step,f2,f3);
    }
    if (step==1|| step%(NSTEP/10)==0){
      int progress=step*100/NSTEP;
      printf("  %3d %%         %6.5f      %8.4f   %5.4f  %5.4f \n", progress, norm, etot,P1,P2) ; 
    }
  }

  return 0;
}
/*------------------------------------------------------------------------------*/
void init_param() {
/*------------------------------------------------------------------------------
  Initializes parameters by reading them from standard input.
------------------------------------------------------------------------------*/
  /* Initialize control parameters */
  LX=50.0;
  DT=0.1;
  NSTEP=100000;
  NECAL=1;
  NNCAL=1000;

  /* Calculate the mesh size */
  dx = LX/NX;
}

/*----------------------------------------------------------------------------*/
void init_prop() {
/*------------------------------------------------------------------------------
  Initializes the kinetic & potential propagators.
------------------------------------------------------------------------------*/
  int i;
  double x, k;
  

  M=2000;
  pot_type=1; 

  /* Set up kinetic propagators */
  for (i=0; i<=NX; i++) {
    if (i < NX/2)
      k = 2*M_PI*i/LX;
    else
      k = 2*M_PI*(i-NX)/LX;
    
    /* kinetic operator */ 
    T[i] = k*k*0.5/M;
    /* kinetic propagator */
    t[i][0] = cos(-DT*T[i]);
    t[i][1] = sin(-DT*T[i]);
  }
  
  /* Set up potential propagator */
  if (pot_type==1) {         // exp up 6 exp down

     A=0.01;
     B=1.6;
     C=0.005;  //SB: 0.005
     D=1.0;  

     for (i=1; i<=NX; i++) {
      x = -0.5*LX + dx*i;
      h[i][0][0] = A*(1-exp(-B*fabs(x)))*x/fabs(x);
      h[i][0][1] = C*exp(-D*x*x);
      h[i][1][0] = h[i][0][1];
      h[i][1][1] = -h[i][0][0];
      if(i==NX/2) {
       h[i][0][0] = 0;
       h[i][0][1] = C;
       h[i][1][0] = C;
       h[i][1][1] = 0; 
      }                 
     }
  }
  else if (pot_type==2) {     //double harmonic  
     A=0.001;
     B=20;
     C=30;
     D=1.0; 

     for (i=1; i<=NX; i++) {
      x = i*dx;
      h[i][0][0] = A*(x-B)*(x-B);
      h[i][0][1] = C*exp(-D*(x-0.5*LX)*(x-0.5*LX));
      h[i][1][0] = h[i][0][1];
      h[i][1][1] = A*(x-C)*(x-C);
     }
   }   //end definition of potentials 

    /* calc eigenvalues of h */
    calc_eigenvalues(i);
    /* calc De and Det */
    calc_De_and_Det(i); 
    /* Half-step diagonal propagator */
    /* 1st component */
    u[i][0][0] = cos(-0.5*DT*E[i][0]);
    u[i][0][1] = sin(-0.5*DT*E[i][0]);
    /* 2nd component */
    u[i][1][0] = cos(-0.5*DT*E[i][1]);
    u[i][1][1] = sin(-0.5*DT*E[i][1]);
}

/*----------------------------------------------------------------------------*/
void init_wavefn() {
/*------------------------------------------------------------------------------
  Initializes the components of the wave function.
------------------------------------------------------------------------------*/
  int sx,s;
  double x,gauss,Csq,norm_fac;

  X0= 17;
  S0= 1.0;    //sb:0.5
  P0= 20;   //sb: 55, do:10

  /* Calculate the the wave function value mesh point-by-point */
  for (sx=1; sx<=NX; sx++) {
    x = dx*sx;
    gauss = exp(-(x-X0)*(x-X0)/4.0/(S0*S0));
    C1[sx][0] = gauss*cos(P0*(x-X0));  	        /* wf on surface 1 */
    C1[sx][1] = gauss*sin(P0*(x-X0));  
    C2[sx][0] = 0;			/* wf on surface 2 */
    C2[sx][1] = 0;
  }
  
  // Normalize C1 if not null 
  Csq=0.0;
  for (sx=1; sx<=NX; sx++)
    for (s=0; s<2; s++)
      Csq += C1[sx][s]*C1[sx][s];
  Csq *= dx;
  if (Csq >0) {
   norm_fac = 1.0/sqrt(Csq);
   for (sx=1; sx<=NX; sx++)
     for (s=0; s<2; s++)
       C1[sx][s] *= norm_fac;
  }

  // Normalize C2 if not null
    Csq=0.0;
    for (sx=1; sx<=NX; sx++)
    for (s=0; s<2; s++)
    Csq += C2[sx][s]*C2[sx][s];
    Csq *= dx;
    if (Csq >0) {
     norm_fac = 1.0/sqrt(Csq);
     for (sx=1; sx<=NX; sx++)
       for (s=0; s<2; s++)
         C2[sx][s] *= norm_fac;  
    }
  periodic_bc();
}

/*----------------------------------------------------------------------------*/
void periodic_bc() {
/*------------------------------------------------------------------------------
  Applies the periodic boundary condition to wave function, by copying
  the boundary values to the auxiliary array positions at the other ends.
------------------------------------------------------------------------------*/
  int s;

  /* Copy boundary wave function values */
  for (s=0; s<=1; s++) {
    C1[0][s] = C1[NX][s];
    C1[NX+1][s] = C1[1][s];
    C2[0][s] = C2[NX][s];
    C2[NX+1][s] = C2[1][s];
  }
}

/*----------------------------------------------------------------------------*/
void four1(double data[], unsigned long nn, int isign)
/*******************************************************************************
Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as
1; or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform,
if isign is input as -1.  data is a complex array of length nn or, equivalently,
a real array of length 2*nn.  nn MUST be an integer power of 2 (this is not
checked for!).
*******************************************************************************/
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;

  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) { /* This is the bit-reversal section of the routine. */
    if (j > i) {
      SWAP(data[j],data[i]); /* Exchange the two complex numbers. */
      SWAP(data[j+1],data[i+1]);
    }
    m=nn;
    while (m >= 2 && j > m) { 
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  
  mmax=2;
  while (n > mmax) { /* Outer loop executed log2 nn times. */
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax); /* Initialize the trigonometric recurrence. */
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) { /* Here are the two nested inner loops. */
      for (i=m;i<=n;i+=istep) {
	j=i+mmax; /* This is the Danielson-Lanczos formula. */
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr; /* Trigonometric recurrence. */
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
/*----------------------------------------------------------------------------*/
void create_C1f() {
  int sx;
  for (sx=1; sx <= NX+1; sx++) {
    Cf[2*sx-1] = C1[sx-1][0];
    Cf[2*sx] = C1[sx-1][1];
  }
}
/*----------------------------------------------------------------------------*/
void create_C2f() {
  int sx;
  for (sx=1; sx <= NX+1; sx++) {
    Cf[2*sx-1] = C2[sx-1][0];
    Cf[2*sx] = C2[sx-1][1];
  }
}
/*----------------------------------------------------------------------------*/
void update_C1() {
  int sx;
  for (sx=1; sx <= NX+1; sx++) {
    C1[sx-1][0] = Cf[2*sx-1];
    C1[sx-1][1] = Cf[2*sx];
  }
}
/*----------------------------------------------------------------------------*/
void update_C2() {
  int sx;
  for (sx=1; sx <= NX+1; sx++) {
    C2[sx-1][0] = Cf[2*sx-1];
    C2[sx-1][1] = Cf[2*sx];
  }
}
/*----------------------------------------------------------------------------*/
void calc_eigenvalues(int i) {
  //E1 = (h11+h22+sqrt((h11-h22)*(h11-h22)+4*h12*h21))/2;  //nb1
  //E2 = (h11+h22-sqrt((h11-h22)*(h11-h22)+4*h12*h21))/2;
  E[i][0] = (h[i][0][0]+h[i][1][1]+sqrt((h[i][0][0]-h[i][1][1])*(h[i][0][0]-h[i][1][1])+4*h[i][0][1]*h[i][1][0]))/2;
  E[i][1] = (h[i][0][0]+h[i][1][1]-sqrt((h[i][0][0]-h[i][1][1])*(h[i][0][0]-h[i][1][1])+4*h[i][0][1]*h[i][1][0]))/2;

}
/*----------------------------------------------------------------------------*/
void calc_De_and_Det(int i) {
  double norm_fac;
  /* De */
  if(i<=NX/2) {
    /* 1st column = 1st eigenvector */
    De[i][0][0] = h[i][0][1]/(E[i][0]-h[i][0][0]);
    De[i][1][0] = 1;
    norm_fac = sqrt(De[i][0][0]*De[i][0][0]+1);
    De[i][0][0] /= norm_fac;
    De[i][1][0] /= norm_fac;
    /* 2nd column = 2nd eigenvector */
    De[i][0][1] = -De[i][1][0];	
    De[i][1][1] = De[i][0][0];		
  }
  else {
    /* 2nd column = 2nd eigenvector */
    De[i][0][1] = h[i][0][1]/(E[i][1]-h[i][0][0]);
    De[i][1][1] = 1;
    norm_fac = sqrt(De[i][0][1]*De[i][0][1]+1);
    De[i][0][1] /= norm_fac;
    De[i][1][1] /= norm_fac;
    /* 1st column = 1st eigenvector */
    De[i][1][0] = -De[i][0][1];	
    De[i][0][0] = De[i][1][1];		
  }
  Det[i][0][0] = De[i][0][0];
  Det[i][1][1] = De[i][1][1];
  Det[i][1][0] = De[i][0][1];
  Det[i][0][1] = De[i][1][0];
}

/*----------------------------------------------------------------------------*/
void pot_prop() {
/*------------------------------------------------------------------------------
  Potential propagator for a half time step, DT/2.
------------------------------------------------------------------------------*/
  int sx;
  double wr1,wi1,wr2,wi2;

  for (sx=1; sx<=NX; sx++) {
    /* De<r|psi> */
    wr1=De[sx][0][0]*C1[sx][0]+De[sx][0][1]*C2[sx][0];
    wi1=De[sx][0][0]*C1[sx][1]+De[sx][0][1]*C2[sx][1];
    wr2=De[sx][1][0]*C1[sx][0]+De[sx][1][1]*C2[sx][0];
    wi2=De[sx][1][0]*C1[sx][1]+De[sx][1][1]*C2[sx][1];
    C1[sx][0]=wr1;
    C1[sx][1]=wi1;
    C2[sx][0]=wr2;
    C2[sx][1]=wi2;
    /* exp(-iE1dt)*C1, exp(-iE2dt)*C2 */
    /* 1st component */
    wr1=u[sx][0][0]*C1[sx][0]-u[sx][0][1]*C1[sx][1];
    wi1=u[sx][0][0]*C1[sx][1]+u[sx][0][1]*C1[sx][0];
    /* 2nd component */
    wr2=u[sx][1][0]*C2[sx][0]-u[sx][1][1]*C2[sx][1];
    wi2=u[sx][1][0]*C2[sx][1]+u[sx][1][1]*C2[sx][0];
    
    C1[sx][0]=wr1;
    C1[sx][1]=wi1;
    C2[sx][0]=wr2;
    C2[sx][1]=wi2;
    /* Det<r|fi> */
    wr1=Det[sx][0][0]*C1[sx][0]+Det[sx][0][1]*C2[sx][0];
    wi1=Det[sx][0][0]*C1[sx][1]+Det[sx][0][1]*C2[sx][1];
    wr2=Det[sx][1][0]*C1[sx][0]+Det[sx][1][1]*C2[sx][0];
    wi2=Det[sx][1][0]*C1[sx][1]+Det[sx][1][1]*C2[sx][1];
    C1[sx][0]=wr1;
    C1[sx][1]=wi1;
    C2[sx][0]=wr2;
    C2[sx][1]=wi2;
  }
  periodic_bc();
}
/*----------------------------------------------------------------------------*/
void kin_prop() {
/*------------------------------------------------------------------------------
  Kinetic propagation for t step.
-------------------------------------------------------------------------------*/
  int sx,s;
  double wr,wi;
  
  for (sx=1; sx<=NX; sx++) {
    /* 1st component C1 */
    wr=t[sx][0]*C1[sx][0]-t[sx][1]*C1[sx][1];
    wi=t[sx][0]*C1[sx][1]+t[sx][1]*C1[sx][0];
    C1[sx][0]=wr;
    C1[sx][1]=wi;
    /* 2nd component C2 */
    wr=t[sx][0]*C2[sx][0]-t[sx][1]*C2[sx][1];
    wi=t[sx][0]*C2[sx][1]+t[sx][1]*C2[sx][0];
    C2[sx][0]=wr;
    C2[sx][1]=wi;
  }	
}
/*----------------------------------------------------------------------------*/
void single_step(int step) {
/*------------------------------------------------------------------------------
  Propagates the wave function for a unit time step, DT.
------------------------------------------------------------------------------*/
  int j;
	
  /* half step potential propagation */
  pot_prop();  
  /* fft of the 2 components of <r|psi> */	
  /* 1st component */	
  create_C1f();
  four1(Cf, (unsigned long) NX, -1);
  for (j=0; j <= 2*(NX+1); j++) 
    Cf[j] /= NX;
  update_C1();
  /* 2nd component */	
  create_C2f();
  four1(Cf, (unsigned long) NX, -1);
  for (j=0; j <= 2*(NX+1); j++) 
    Cf[j] /= NX;
  update_C2();
  /* step kinetic propagation   */
  kin_prop(); 
  if (step%NECAL==0)
    calc_ekin();
  /* fft^(-1) */
  /* 1st component */
  create_C1f();
  four1(Cf, (unsigned long) NX, 1);
  update_C1();
  /* 2nd component */
  create_C2f();
  four1(Cf, (unsigned long) NX, 1);
  update_C2();
  /* half step potential propagation */
  pot_prop();  
  if (step%NECAL==0)		
    calc_epot();	
}

/*----------------------------------------------------------------------------*/
void pop_states() {
  int i;
  P1=0.0;
  P2=0.0;
  for (i=1; i<=NX;i++) {
    P1 += C1[i][0]*C1[i][0]+C1[i][1]*C1[i][1];
    P2 += C2[i][0]*C2[i][0]+C2[i][1]*C2[i][1];
  }
  P1 *= dx;
  P2 *= dx;
}
/*----------------------------------------------------------------------------*/

void print_pop(int step, FILE *f5) {
  int i;

    fprintf(f5,"%8i %15.10f %15.10f\n",step,P1,P2);
} 
      
/*----------------------------------------------------------------------------*/
void calc_ekin() {
  int sx;
  double k;
  /* kinetic energy */
  ekin = 0.0;
  for (sx=1; sx<=NX; sx++) {
    if (sx < NX/2)
      k = 2*M_PI*sx/LX;
    else
      k = 2*M_PI*(sx-NX)/LX;
      //printf("%f\n",k);
  ekin += k*k*(C1[sx][0]*C1[sx][0]+C1[sx][1]*C1[sx][1])*0.5/M +
          k*k*(C2[sx][0]*C2[sx][0]+C2[sx][1]*C2[sx][1])*0.5/M;     //domod "/M" 
  }
  ekin *= dx;
  ekin *= NX;
}

/*----------------------------------------------------------------------------*/
void calc_epot() {
  int sx;
  /* Potential energy */
  epot = 0.0;
  for (sx=1; sx<=NX; sx++) {
    epot += h[sx][0][0]*(C1[sx][0]*C1[sx][0]+C1[sx][1]*C1[sx][1])-
            2.0*h[sx][0][1]*(C1[sx][0]*C2[sx][0]+C1[sx][1]*C2[sx][1])+
            h[sx][1][1]*(C2[sx][0]*C2[sx][0]+C2[sx][1]*C2[sx][1]);   //nb1 

  }
  epot *= dx;
}

/*----------------------------------------------------------------------------*/
void calc_energy() {
  /* Total energy */
  etot = ekin+epot;
}

/*----------------------------------------------------------------------------*/
void print_energy(int step, FILE *f1) {
/*------------------------------------------------------------------------------
  Print energy values in file energy.dat
-------------------------------------------------------------------------------*/
  
 fprintf(f1, "%8i %15.10f %15.10f %15.10f\n", step,ekin,epot,etot); 
  
}
/*-------------------------------------------------------------------------------
  Print the potential
-------------------------------------------------------------------------------*/

void print_pot_ad(FILE *f4) {
  
  int i;
  double x;

  for (i=1; i<=NX; i++) {
    x = dx*i;
    fprintf(f4,"%8i %15.10f %15.10f %15.10f\n",i,x,E[i][1],E[i][0]);  //nb1
  }
}
/*-------------------------------------------------------------------------------*/
void print_pot_di(FILE *f6) {
  
  int i;
  double x;

  for (i=1; i<=NX; i++) {
    x = dx*i;
    fprintf(f6,"%8i %15.10f %15.10f %15.10f %15.10f %15.10f\n"
                ,i,x,h[i][0][0],h[i][1][1],h[i][0][1],h[i][1][0]);
  }
}
/*-------------------------------------------------------------------------------*/
void calc_norm() {
/*------------------------------------------------------------------------------
  Calculate the norm                        
-------------------------------------------------------------------------------*/
  
 int sx;
 double psisq,psisq2;

 norm=0.0;

 for (sx=1; sx<=NX; sx++)                              // domod from 1 to NX (not NX+1) 
 {
  psisq=C1[sx][0]*C1[sx][0]+C1[sx][1]*C1[sx][1]+
        C2[sx][0]*C2[sx][0]+C2[sx][1]*C2[sx][1];
  psisq2=C1[sx+1][0]*C1[sx+1][0]+C1[sx+1][1]*C1[sx+1][1]+
         C2[sx+1][0]*C2[sx+1][0]+C2[sx+1][1]*C2[sx+1][1];
  norm=norm+dx*((psisq2+psisq)/2.0);
 }
  norm=sqrt(norm);
}

/*-------------------------------------------------------------------------------*/
void print_wavefn(int step, FILE *f2, FILE *f3) {
/*------------------------------------------------------------------------------
  Print energy values in file energy.dat
-------------------------------------------------------------------------------*/
  
 int sx;
 double x;

 fprintf(f2,"\n"); 
 fprintf(f2,"\n"); 
 for (sx=1; sx<=NX; sx++)
 {
  x=dx*sx;
  fprintf(f2,"%8i %15.10f %15.10f %15.10f\n",sx,x,
    ((C1[sx][0]*C1[sx][0]+C1[sx][1]*C1[sx][1]))+h[sx][0][0],  
    ((C2[sx][0]*C2[sx][0]+C2[sx][1]*C2[sx][1]))+h[sx][1][1]);  
 }

  if (step>0)
  {
  fprintf(f3,"%8i %15.10f\n",step,norm);
  }
}
/*----------------------------------------------------------------------------*/
