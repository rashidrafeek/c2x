/* Write a CIF file */

/* Copyright (c) 2014 MJ Rutter 
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

/* Write a file in format similar to that which Castep would use 
 * for CIF output. Though this format may be compatible with CIF / mmCIF,
 * it is very far from making use of all the features of that format. */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

#include "c2xsf.h"

#define MAX_ELS 103

void cif_write(FILE* outfile, struct unit_cell *c, struct contents *m, 
	       struct symmetry *s, int mm){
  int i,j,prec;
  unsigned char sep;
  double abc[6];

  sep='_';
  if (mm) sep='.';

  prec=7;
  if (flags&HIPREC) prec=15;

  fprintf(outfile,"# written by c2x\n");

  fprintf(outfile,"data_c2x\n");

  if (m->title) fprintf(outfile,"_struct%ctitle %s\n",sep,m->title);

  cart2abc(c,m,abc,NULL,1);
  fprintf(outfile,"\n");
  fprintf(outfile,"_cell%clength_a      %.*f\n",sep,prec,abc[0]);
  fprintf(outfile,"_cell%clength_b      %.*f\n",sep,prec,abc[1]);
  fprintf(outfile,"_cell%clength_c      %.*f\n",sep,prec,abc[2]);
  fprintf(outfile,"_cell%cangle_alpha   %.*f\n",sep,prec,abc[3]);
  fprintf(outfile,"_cell%cangle_beta    %.*f\n",sep,prec,abc[4]);
  fprintf(outfile,"_cell%cangle_gamma   %.*f\n",sep,prec,abc[5]);
  fprintf(outfile,"\n");

  fprintf(outfile,"\n");
  fprintf(outfile,"loop_\n"
	  "_atom_site%ctype_symbol\n"
	  "_atom_site%cfract_x\n"
	  "_atom_site%cfract_y\n"
	  "_atom_site%cfract_z\n"
	  "_atom_site%cU_iso_or_equiv\n"
	  "_atom_site%coccupancy\n",sep,sep,sep,sep,sep,sep);

  if (prec<10) prec=10;

  for(i=0;i<m->n;i++){
    fprintf(outfile,"%s ",atno2sym(m->atoms[i].atno));
    for(j=0;j<3;j++)
      fprintf(outfile,"%.*f ",prec,m->atoms[i].frac[j]);
    fprintf(outfile,"0.01 1.00\n");
  }

  if (s->n){
    fprintf(outfile,"\nloop_\n");
    if (mm==0){
      fprintf(outfile,"_symmetry_equiv_pos_as_xyz\n");
      for(i=0;i<s->n;i++)
	equiv_sym(s->ops+i,c,outfile);
    }
    else{
      fprintf(outfile,"_symmetry_equiv.id\n");
      fprintf(outfile,"_symmetry_equiv.pos_as_xyz\n");
      for(i=0;i<s->n;i++){
	fprintf(outfile,"%d ",i+1);
	equiv_sym(s->ops+i,c,outfile);
      }
    }
  }
}
