/* Reader for a Siesta .XV and .KP files */

/* Copyright (c) 2019-2021 MJ Rutter 
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

/* XV files are in Bohr
   3 x (cell vector, cell vector velocity)
   natoms
   natoms x (spec no, atomic no, abs coords, velocity)
 */

#include<stdio.h>
#include<stdlib.h> /* malloc */

#include "c2xsf.h"

#define LINE_SIZE 200

double siesta_pot_charge(char *prefix, char *label); /* From fdf_read.c */

void xv_read(FILE* infile, struct unit_cell *c, struct contents *m){
  int i,j,species;
  double dummy,*species_to_charge,x;
  char buffer[LINE_SIZE+1],*prefix;

  species_to_charge=(double*)dict_get(m->dict,"Siesta_species_to_charge");
  
  if (!(c->basis=malloc(9*sizeof(double))))
    error_exit("Malloc error in xv_read for c->basis");

  for(i=0;i<3;i++){
    if (!fgets(buffer,LINE_SIZE,infile))
      error_exit("Unexpected EOF reading basis");
    j=sscanf(buffer,"%lf %lf %lf %lf %lf %lf",c->basis[i],c->basis[i]+1,
	     c->basis[i]+2,&dummy,&dummy,&dummy);
    /* Be tolerant of some unexpected line breaks */
    if ((j>=3)&&(j<6)){
      if (!fgets(buffer,LINE_SIZE,infile))
	error_exit("Unexpected EOF reading basis");
      j+=sscanf(buffer,"%lf %lf %lf",&dummy,&dummy,&dummy);
    }
    if (j!=6) error_exit("Error reading basis");
  }

  /* Basis was in Bohr */

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[i][j]*=BOHR;

  if (!fgets(buffer,LINE_SIZE,infile))
    error_exit("Unexpected EOF reading natoms");
  i=sscanf(buffer,"%d",&(m->n));
  if (i!=1) error_exit("Failed to read natoms");

  m->atoms=malloc(m->n*sizeof(struct atom));
  if (!m->atoms) error_exit("Malloc error for m->atoms in xv_read");
  init_atoms(m->atoms,m->n);
  
  for(i=0;i<m->n;i++){
    if (!fgets(buffer,LINE_SIZE,infile))
      error_exit("Unexpected EOF reading atoms");
    j=sscanf(buffer,"%d %d %lf %lf %lf %lf %lf %lf",&species,
	     &(m->atoms[i].atno),m->atoms[i].abs,
	     m->atoms[i].abs+1,m->atoms[i].abs+2,
             m->atoms[i].v,m->atoms[i].v+1,m->atoms[i].v+2);
    /* Be tolerant of some unexpected line breaks */
    if (j==5){
      if (!fgets(buffer,LINE_SIZE,infile))
	error_exit("Unexpected EOF reading atoms");
      j+=sscanf(buffer,"%lf %lf %lf",
                m->atoms[i].v,m->atoms[i].v+1,m->atoms[i].v+2);
    }
    if (j!=8) error_exit("Error parsing atom entry for 8 fields");
    for(j=0;j<3;j++)
      m->atoms[i].abs[j]*=BOHR;
    for(j=0;j<3;j++) /* velocities are in Bohr/fs, we want A/ps */
      m->atoms[i].v[j]*=1000*BOHR;
    if (species_to_charge) m->atoms[i].chg=species_to_charge[species];
  }

  real2rec(c);
  addfrac(m->atoms,m->n,c->recip);
  addspec(m);

 /* Don't set m->velocities if all v's are zero */
  x=0;
  for(i=0;i<m->n;i++)
    x+=m->atoms[i].v[0]*m->atoms[i].v[0]+
      m->atoms[i].v[1]*m->atoms[i].v[1]+m->atoms[i].v[2]*m->atoms[i].v[2];
  if (x>0) m->velocities=1;

  if ((!species_to_charge)&&((flags&CHDEN)||(flags&OCCUPANCIES))){
    prefix=dict_get(m->dict,"in_dir");
    for(i=0;i<m->nspec;i++){
      x=siesta_pot_charge(prefix,atno2sym(m->spec[i].atno));
      for(j=0;j<m->n;j++)
	if (m->atoms[j].atno==m->spec[i].atno) m->atoms[j].chg=x;
    }
  }

}

void siesta_kp_read(FILE* infile, struct kpts *k){
  int i,j;
  char buffer[LINE_SIZE+1];
  
  if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
  if (sscanf(buffer,"%d",&(k->n))!=1)
    error_exit("parse error on first line of KP file");

  if (debug) fprintf(stderr,"Reading %d kpoints from Siesta KP file\n",
		     k->n);
  
  k->kpts=malloc(k->n*sizeof(struct atom));
  if (!k->kpts) error_exit("malloc error for kpoints");
  init_atoms(k->kpts,k->n);
  
  for(i=0;i<k->n;i++){
    if (!fgets(buffer,LINE_SIZE,infile)){
      fprintf(stderr,"Unexpected EOF after %d lines\n",i+1);
      exit(1);
    }
    if (sscanf(buffer,"%d %lf %lf %lf %lf",&j,k->kpts[i].frac,
	       k->kpts[i].frac+1,k->kpts[i].frac+2,&(k->kpts[i].wt))!=5){
      fprintf(stderr,"Error parsing line %d\n",i+1);
      exit(1);
    }
    if (j!=i+1){
      fprintf(stderr,"Confused by kpt ordering. Expected %d, found %d\n",
	      i+1,j);
      exit(1);
    }
  }
  
}
