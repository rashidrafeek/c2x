/* Write an xsf file (XCrySDen). */

/* Copyright (c) 2007 MJ Rutter 
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

/* The format is:

CRYSTAL
PRIMVEC
[lattice_vector_1_x] [y] [z]
[lattice_vector_2_x] [y] [z]
[lattice_vector_3_x] [y] [z]
CONVVEC
[lattice_vector_1_x] [y] [z]
[lattice_vector_2_x] [y] [z]
[lattice_vector_3_x] [y] [z]
PRIMCOORD
[number_of_atoms] 1
[atom_symbol] [pos_x] [pos_y] [pos_z]
....
BEGIN_BLOCK_DATAGRID_3D
Densities
BEGIN_DATAGRID_3D_chden
[grid_points_vec1+1] [grid_points_vec2+1] [grid_points_vec3+1]
[offset_x] [offset_y] [offset_z] 
[grid_vector_1_x] [y] [z]
[grid_vector_2_x] [y] [z]
[grid_vector_3_x] [y] [z]
[data x=0 y=0 z=0]
[data x=1 y=0 z=0]
...
[data x=grid_points_vec1 y=0 z=0]
[data x=0 y=0 z=0]
[data x=0 y=1 z=0]
...
END_DATAGRID_3D
BEGIN_DATAGRID_3D_spin (optional, as
END_DATAGRID_3D         for chden above)
END_BLOCK_DATAGRID_3D

Comments start with a #

For added confusion, and to convert from Fortran to C ordering, this
code used to output the grid vectors with x and z exchanged, and thus exchanged
x and z when writing the data...

No longer, for VESTA reads the XSF format, and silently misinterprets it
unless the vectors in the datagrid section are the same as in the cell
description. It also gets upset if the first line is blank. A line
containing a single '#' is fine.

*/

#include<stdio.h>
#include<stdlib.h>

#include "c2xsf.h"

void xsf_write(FILE* outfile, struct unit_cell *c, struct contents *m,
               int molecule, struct grid *gptr){
  int i,j,k;
  double x,y,z,*dptr1,*dptr2;
  char *fmtf,*fmtd6f,*fmtd3f,*fmt3f;
  struct cmt *comment;
  
  if (flags&HIPREC){
    fmtf="% 19.15f";
    fmt3f="% 19.15f % 19.15f % 19.15f\n";
    fmtd3f="%d % 19.15f % 19.15f % 19.15f\n";
    fmtd6f="%d % 19.15f % 19.15f % 19.15f % 19.15f % 19.15f % 19.15f\n";
  }
  else{
    fmtf="% 11.7f\n";
    fmt3f="% 11.7f % 11.7f % 11.7f\n";
    fmtd3f="%d % 11.7f % 11.7f % 11.7f\n";
    fmtd6f="%d % 11.7f % 11.7f % 11.7f % 11.7f % 11.7f % 11.7f\n";
  }
  
/* Atomic positions are required to lie within unit cell, so: */
  reduce_cell(m->atoms,m->n,c->basis);

  if (m->title) fprintf(outfile,"# %s\n\n",m->title);
  if (m->comment->txt){
    fprintf(outfile,"#\n");
    comment=m->comment;
    while((comment)&&(comment->txt)){
      fprintf(outfile,"# %s\n",comment->txt);
      comment=comment->next;
    }
    fprintf(outfile,"#\n");
  }
  
  if (molecule==0){
    fprintf(outfile,"CRYSTAL\nPRIMVEC\n");
    for(i=0;i<=2;i++) fprintf(outfile,fmt3f,c->basis[i][0],
                              c->basis[i][1],c->basis[i][2]);
    fprintf(outfile,"CONVVEC\n");
    for(i=0;i<=2;i++) fprintf(outfile,fmt3f,c->basis[i][0],
                              c->basis[i][1],c->basis[i][2]);
    fprintf(outfile,"PRIMCOORD\n");
    fprintf(outfile,"%d 1\n",m->n);
  }else{
    fprintf(outfile,"MOLECULE\nATOMS\n");
  }

    
/* Need to write coords in Cartesian basis */
  for(i=0;i<m->n;i++){
    x=m->atoms[i].abs[0];
    y=m->atoms[i].abs[1];
    z=m->atoms[i].abs[2];
    if (m->forces)
      fprintf(outfile,fmtd6f,m->atoms[i].atno,x,y,z,
              m->atoms[i].force[0],m->atoms[i].force[1],m->atoms[i].force[2]);
    else
      fprintf(outfile,fmtd3f,m->atoms[i].atno,x,y,z);
  }

  if((gptr)&&(gptr->data)){
    fprintf(outfile,"BEGIN_BLOCK_DATAGRID_3D\n"
                    "Densities\n");
    while((gptr)&&(gptr->data)){
      if (debug>1) fprintf(stderr,"Writing %s\n",gptr->name);
      fprintf(outfile,"BEGIN_DATAGRID_3D_%s\n",gptr->name);
      fprintf(outfile,"%d %d %d\n",gptr->size[0]+1,
                      gptr->size[1]+1,gptr->size[2]+1);
      fprintf(outfile,"0.0 0.0 0.0\n");
      for(i=0;i<3;i++) fprintf(outfile,fmt3f,c->basis[i][0],
                                  c->basis[i][1],c->basis[i][2]);

      dptr2=gptr->data;
#if 0
      for(k=0;k<=gptr->size[0];k++){
        for(j=0;j<=gptr->size[1];j++){
          dptr1=dptr2+
           (((k%gptr->size[0])*gptr->size[1])+(j%gptr->size[1]))*gptr->size[2];
          for(i=0;i<=gptr->size[2];i++){
            fprintf(outfile,"%f\n",*(dptr1+i%gptr->size[2]));
          }
        }
      }
#endif
      for(i=0;i<=gptr->size[2];i++){
        for(j=0;j<=gptr->size[1];j++){
          dptr1=dptr2+(j%gptr->size[1])*gptr->size[2]+i%gptr->size[2];
          for(k=0;k<=gptr->size[0];k++){
            fprintf(outfile,fmtf,
                    *(dptr1+(k%gptr->size[0])*gptr->size[1]*gptr->size[2]));
          }
        }
      }
      fprintf(outfile,"END_DATAGRID_3D\n");
      gptr=gptr->next;
    }

    fprintf(outfile,"END_BLOCK_DATAGRID_3D\n");
  }
}
