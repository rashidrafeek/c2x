/* A reader for Crystal's fort.34 files */

/* Copyright (c) 2018 MJ Rutter 
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
#include<stdlib.h> /* malloc */

#include "c2xsf.h"

#define LINE_SIZE 2049

void fort34_read(FILE* infile, struct unit_cell *c, struct contents *m,
              struct symmetry *s){
  int i,j,k,nsym,natoms;
  char buffer[LINE_SIZE+1];
  double *ptr,x[3];

  buffer[0]=0;
  i=0;
  
  fgets(buffer,LINE_SIZE,infile);

  sscanf(buffer,"%d %d %d",&i,&j,&k);

  if (i!=3) error_exit("Not a valid fort.34 describing a crystalline system");

  if (j!=1) fprintf(stderr,"Warning: ignoring centring information\n");

  if (!(c->basis=malloc(9*sizeof(double))))
    error_exit("Malloc error in fort34_read for c->basis");

  ptr=(double*)c->basis;
  for(i=0;i<3;i++){
    fgets(buffer,LINE_SIZE,infile);
    if (sscanf(buffer,"%lf %lf %lf",ptr,ptr+1,ptr+2)!=3)
      error_exit("Error reading basis in fort.34");
    ptr+=3;
  }

  real2rec(c);
  
  nsym=0;
  fgets(buffer,LINE_SIZE,infile);
  if (sscanf(buffer,"%d",&nsym)!=1)
    error_exit("Error reading number of symmetry ops in fort.34");
    
  s->ops=malloc(nsym*sizeof(struct sym_op));
  if (!s->ops) error_exit("Malloc error for symops");

  for(i=0;i<nsym;i++){
    ptr=(double*)s->ops[i].mat;
    for(j=0;j<3;j++){
      fgets(buffer,LINE_SIZE,infile);
      if (sscanf(buffer,"%lf %lf %lf",ptr,ptr+1,ptr+2)!=3)
        error_exit("Error reading sym_op in fort.34");
      ptr+=3;
    }
    fgets(buffer,LINE_SIZE,infile);
    if (sscanf(buffer,"%lf %lf %lf",x,x+1,x+2)!=3)
        error_exit("Error reading sym_op in fort.34");
    if((x[0]==0)&&(x[1]==0)&&(x[2]==0))
      s->ops[i].tr=NULL;
    else{
      s->ops[i].tr=malloc(3*sizeof(double));
      if (!s->ops[i].tr) error_exit("Malloc error for symop.tr");
      for(j=0;j<3;j++)
        s->ops[i].tr[j]=x[j];
    }
  }
  s->n=nsym;

  natoms=0;
  fgets(buffer,LINE_SIZE,infile);
  if (sscanf(buffer,"%d",&natoms)!=1)
    error_exit("Error reading number of atoms in fort.34");

  m->atoms=malloc(natoms*sizeof(struct atom));
  if (!m->atoms) error_exit("Malloc error for atoms");
  m->n=natoms;
  init_atoms(m->atoms,m->n);
  
  for(i=0;i<natoms;i++){
    fgets(buffer,LINE_SIZE,infile);
    ptr=m->atoms[i].abs;
    if (sscanf(buffer,"%d %lf %lf %lf",&m->atoms[i].atno,ptr,ptr+1,ptr+2)!=4)
      error_exit("Error parsing atom in fort.34");
    m->atoms[i].atno=m->atoms[i].atno%100;
  }
  addfrac(m->atoms,natoms,c->recip);

  sym_expand(c,m,s);
  
}
