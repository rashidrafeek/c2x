/* If passed a non-null pointer, assume is pointer to three integers
 * Shift everything by these numbers of FFT grid cells
 *
 * If pointer is null, attempt to work out what shift would be required
 * in x,y and z independently.
 *
 * Assumes all 3D data are on same FFT grid
 *
 * Revised 2019 so that automatic mode tries to centre molecule in cell,
 *  and does so precisely if there is no grid data to worry about
 */


/* Copyright (c) 2007,2019 MJ Rutter 
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

void molecule_fix(int m_abc[3], double m_rel[3], struct unit_cell *c,
		  struct contents *m,
                  struct grid *gptr){
  int i,j,k,ii,jj,kk,off1,off2;
  int fft[3],shift[3];
  double amin,amax,ashift[3],*grid3;

  if (gptr->data){
    fft[0]=gptr->size[0];
    fft[1]=gptr->size[1];
    fft[2]=gptr->size[2];
  }
  else fft[0]=fft[1]=fft[2]=1000;


  if(!m_abc){ /* We are expected to determine the shift automatically */

    if (m_rel){
      for(i=0;i<3;i++) ashift[i]=m_rel[i];
    }
    else{
      /* Have we one data grid with an offset origin? */
      if ((gptr)&&(gptr->origin_abs)&&((!gptr->next)||(!gptr->next->data))){
	if (debug) fprintf(stderr,"Shifting atoms by (%f,%f,%f) A\n",
			   -gptr->origin_abs[0],-gptr->origin_abs[1],
			   -gptr->origin_abs[2]);
	for(i=0;i<m->n;i++)
	  for(k=0;k<3;k++)
	    m->atoms[i].abs[k]-=gptr->origin_abs[k];
	addfrac(m->atoms,m->n,c->recip);
	free(gptr->origin_abs);
	gptr->origin_abs=NULL;
	return;
      } /* General case */
      else{
	for(k=0;k<3;k++){
	  amin=1000;
	  amax=-1000;
	  for(i=0;i<m->n;i++){
	    if (m->atoms[i].frac[k]>amax)
	      amax=m->atoms[i].frac[k];
	    if (m->atoms[i].frac[k]<amin)
	      amin=m->atoms[i].frac[k];
	  }
	  ashift[k]=0.5-0.5*(amax+amin);
	}
      }
    }
    if (!gptr->data){
      if ((!m_rel)||(debug))
	fprintf(stderr,"Shifting by (%.6f,%.6f,%.6f)\n",ashift[0],
		ashift[1],ashift[2]);
      for(i=0;i<m->n;i++)
        for(k=0;k<3;k++)
          m->atoms[i].frac[k]+=ashift[k];

      addabs(m->atoms,m->n,c->basis);
      return;
    }
    else{
      m_abc=malloc(3*sizeof(int));
      if (!m_abc) error_exit("Malloc error in molecule_fix");
      for(k=0;k<3;k++) m_abc[k]=floor(ashift[k]*fft[k]+0.5);
    }
  }

  if (debug)
      fprintf(stderr,"Will translate grid by (%d,%d,%d) gridpoints\n"
                     " which is (%f,%f,%f)\n",
              m_abc[0],m_abc[1],m_abc[2],(double)m_abc[0]/fft[0],
              (double)m_abc[1]/fft[1],(double)m_abc[2]/fft[2]);

  /* Ions are easy */

  for(i=0;i<m->n;i++)
    for(k=0;k<3;k++)
      m->atoms[i].frac[k]+=(double)m_abc[k]/fft[k];

  addabs(m->atoms,m->n,c->basis);

    /* Grids are harder */

  while((gptr)&&(gptr->data)){

    grid3=malloc(sizeof(double)*gptr->size[0]*gptr->size[1]*gptr->size[2]);
    if (!grid3) error_exit("Memory allocation error in grid shift");
    for(k=0;k<3;k++){
      shift[k]=(gptr->size[k]-(m_abc[k]*gptr->size[k])/fft[k])%gptr->size[k];
      while(shift[k]<0) shift[k]+=gptr->size[k];
    }

    for(k=0;k<gptr->size[0];k++){
      kk=(k+shift[0])%gptr->size[0];
      for(j=0;j<gptr->size[1];j++){
        jj=(j+shift[1])%gptr->size[1];
        off1=(k*gptr->size[1]+j)*gptr->size[2];
        off2=(kk*gptr->size[1]+jj)*gptr->size[2];
        for(i=0;i<gptr->size[2];i++){
          ii=(i+shift[2])%gptr->size[2];
          *(grid3+off1+i)=*(gptr->data+off2+ii);
        }
      }
    }
    free(gptr->data);
    gptr->data=grid3;
    gptr=gptr->next;
  }
}
