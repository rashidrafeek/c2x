/* Write a Gaussian cube file, one density only
 *       Units must be Bohr, not Angstoms
 */


/* Copyright (c) 2007-2017 MJ Rutter 
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

#include "c2xsf.h"

void cube_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                struct grid *g){
  int i,j,k;
  double x,y,z,*dptr1,*dptr2;
  char *fmt1,*fmt2,*fmt3;
  
  if (!g->data){
    fprintf(stderr,"Cannot write .cube file when no 3D data requested.\n");
    exit(1);
  }

  fprintf(outfile,"DENSITY:\n\n");

  fprintf(outfile,"%d 0.0 0.0 0.0\n",m->n);

  if (flags&HIPREC){
    fmt1="%d %.15f %.15f %.15f\n";
    fmt2="%d 0.0 %.15f %.15f %.15f\n";
    fmt3="%.15f\n";
  }
  else{
    fmt1="%d %f %f %f\n";
    fmt2="%d 0.0 %f %f %f\n";
    fmt3="%f\n";
  }
      
  
  for(i=0;i<3;i++) fprintf(outfile,fmt1,g->size[i],
                                c->basis[i][0]/g->size[i]/BOHR,
                                c->basis[i][1]/g->size[i]/BOHR,
                                c->basis[i][2]/g->size[i]/BOHR);

/* Need to write coords in Cartesian basis */
  for(i=0;i<m->n;i++){
    x=m->atoms[i].abs[0]/BOHR;
    y=m->atoms[i].abs[1]/BOHR;
    z=m->atoms[i].abs[2]/BOHR;
    fprintf(outfile,fmt2,m->atoms[i].atno,x,y,z);
  }

  dptr2=g->data;
  for(k=0;k<g->size[0];k++){
    for(j=0;j<g->size[1];j++){
      dptr1=dptr2+((k*g->size[1])+j)*g->size[2];
      for(i=0;i<g->size[2];i++)
        fprintf(outfile,fmt3,*(dptr1+i));
    }
  }
}
