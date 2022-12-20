/* Reader for a Siesta .XV file */

/* Copyright (c) 2019 MJ Rutter 
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

void xv_read(FILE* infile, struct unit_cell *c, struct contents *m){
  int i,j,species;
  double dummy,*species_to_charge;
  char buffer[LINE_SIZE+1];

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
	     m->atoms[i].abs+1,m->atoms[i].abs+2,&dummy,&dummy,&dummy);
    /* Be tolerant of some unexpected line breaks */
    if ((j>=5)&&(j<8)){
      if (!fgets(buffer,LINE_SIZE,infile))
	error_exit("Unexpected EOF reading atoms");
      j+=sscanf(buffer,"%lf %lf %lf",&dummy,&dummy,&dummy);
    }
    for(j=0;j<3;j++)
      m->atoms[i].abs[j]*=BOHR;
    if (species_to_charge) m->atoms[i].chg=species_to_charge[species];
  }

  real2rec(c);
  addfrac(m->atoms,m->n,c->recip);

}
