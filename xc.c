/* Calculate XC pot */

/* Copyright (c) 2007-2019 MJ Rutter 
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3
 * of the Licence, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see http://www.gnu.org/licenses/
 */ 


#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "c2xsf.h"

/* Calculate XC pot according to Perdew and Zunger
 * 
 * PRB 23 5048 (1981) 
 */

double lda(double dens){
  double r_ws,Ec,uc,ux;
  double B1,B2,gamma,A,B,C,D;
  double alpha,rscfcv;

  if (dens<1e-30) return 0;
  
  A= 0.0311;
  B=-0.048;
  C= 0.0020;
  D=-0.0116;
  B1=1.0529;
  B2=0.3334;
  gamma=-0.1423;  /* Castep has -3.782211 = 27.211 * gamma */
  alpha=pow(9*M_PI/4,1./3.);
  rscfcv=3*alpha/(4*M_PI);  /* = 0.4581653 */
  /* Castep has 6.59747 = 27.211 * 0.5291 * 0.4581653 */

  /* Find Wigner-Seitz radius */

  r_ws=pow(3/(4*M_PI*dens),1./3.);

  /* and convert to Bohr */

  r_ws/=BOHR;

  if (debug>2) fprintf(stderr,"Wigner-Seitz radius: %lf Bohr\n",r_ws);

  ux=-4*rscfcv/(3*r_ws);
  
  if (r_ws>=1.0){
    Ec=gamma/(1+B1*sqrt(r_ws)+B2*r_ws);
    uc=Ec*(1+(7*B1/6)*sqrt(r_ws)+(4*B2/3)*r_ws)/(1+B1*sqrt(r_ws)+B2*r_ws);
  }
  else
    uc=A*log(r_ws)+(B-A/3)+(2*C*r_ws/3)*log(r_ws)+(2*D-C)*r_ws/3;

  return (ux+uc)*H_eV;

}

#ifdef TEST
int debug;

int main(int argc, char **argv)
{
  double dens,pot;

  debug=1;
  
  sscanf(argv[1],"%lf",&dens);

  printf("Density: %lf eA^-3\n",dens);

  pot=lda(dens);

  printf("Potential: %lf V\n",pot);

  return 0;
}
#endif
