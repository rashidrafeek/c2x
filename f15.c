/* Copyright (c) 2013, 2020 MJ Rutter 
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
#include<string.h>
#include<stdlib.h>
#include<ctype.h>

#include "c2xsf.h"

void f15_write(struct unit_cell *c, struct contents *m, struct kpts *k){
  int i,j,nspec;
  int *natomsp,*spchg;
  FILE *fort15,*fort14,*fort4;

  fort15=fopen("fort.15","w");

  /* Write basis, transposed, twice */
  for(j=0;j<2;j++)
    for(i=0;i<3;i++)
      fprintf(fort15,"%.7f %.7f %.7f\n",c->basis[0][i],c->basis[1][i],
	      c->basis[2][i]);

  if (!m->spec) addspec(m);
  nspec=m->nspec;

  natomsp=malloc(nspec*sizeof(int));
  spchg=malloc(nspec*sizeof(int));
  if ((!natomsp)||(!spchg)) error_exit("Malloc error in f15_write");
  for(i=0;i<nspec;i++) natomsp[i]=0;
  for(i=0;i<nspec;i++)
    for(j=0;j<m->n;j++)
      if (m->atoms[j].atno==m->spec[i].atno){
	natomsp[i]++;
	spchg[i]=(int)m->atoms[j].chg;
      }
  
  /* Write out atoms, sorted by species */
  for(i=0;i<nspec;i++){
    for(j=0;j<m->n;j++)
      if (m->atoms[j].atno==m->spec[i].atno)
        fprintf(fort15," %f %f %f 0.0\n",m->atoms[j].frac[0],
                m->atoms[j].frac[1],m->atoms[j].frac[2]);
    for(j=0;j<m->n;j++)
      if (m->atoms[j].atno==m->spec[i].atno)
        fprintf(fort15," %f %f %f\n",m->atoms[j].frac[0],
                m->atoms[j].frac[1],m->atoms[j].frac[2]);
  }

  for(i=0;i<k->n;i++)
    fprintf(fort15,"%.8f %.8f %.8f\n",k->kpts[i].frac[0],
	    k->kpts[i].frac[1],k->kpts[i].frac[2]);
  for(i=0;i<k->n;i++)
    fprintf(fort15,"%.8f\n",k->kpts[i].wt);

  fclose(fort15);

  fort14=fopen("fort.14","w");
  fprintf(fort14,"1  NITER\n"
	  "1  NPRINT\n"
	  "1  ENMAX\n"
	  "1  NLPLOT\n"
	  "1  NITMAX\n"
	  "1  NDELAY\n"
	  "1  ISTART\n"
	  "1  ISBROT\n"
	  "1  IOCCUP\n"
	  "1  DELMIN\n"
	  "1  DELMAX\n"
	  "1  NDEL\n"
	  "1  IION\n"
	  "1  IBOX\n"
	  "1  IPRINT\n"
	  "1  NIONCG\n"
	  "1  NITFIX\n"
	  "1  INRAND\n"
	  "1  ICLOCK\n"
	  "1  SITIM\n"
	  "1  SIDAMP\n"
	  "1  SIMASS\n"
	  "1  SIDISP\n"
	  "1  POTIM\n"
	  "1  PODISP\n");
  for(i=0;i<nspec;i++)
    fprintf(fort14,"%d ICHARGE\n",spchg[i]);
  for(i=0;i<nspec;i++)
    fprintf(fort14,"1 POMASS\n");
  for(i=0;i<nspec;i++)
    fprintf(fort14,"%d NIONSP\n",natomsp[i]);
  fclose(fort14);

  fort4=fopen("fort.4","w");
  fprintf(fort4,"%d       ! NSPEC\n",nspec);
  fprintf(fort4,"800 15 15 15 0            ! A2MAX NA1 NA2 NA3 EPSILO\n"
	  "0 0                       ! No more detail of direct lattice\n"
	  "-1                        ! -1 No repeat generation\n"
	  "0 %d 4 0 0 0               ! IQ1 IQ2 IQ3 WVK0\n"
	  "-1                        ! Symmetrize mesh\n"
	  "1                         ! Save special points on file\n"
	  "0 0 0 0 0 0               ! Exit K290\n",k->n);
  fprintf(fort4,"\n"
	  "This fort.4 is for MJR's version of k290.\n"
	  "The first line must be NSPEC - the value in fort.2 is \n"
	  "  treated as a dummy\n"
	  "If IQ1=0, then IQ2 k-points will be read from fort.15\n");
  fclose(fort4);
  free(spchg);
  free(natomsp);
}

