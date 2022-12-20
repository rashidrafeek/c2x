/* Read an extended xyz format file */

/* Copyright (c) 2020 MJ Rutter 
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
#include<string.h>
#include<ctype.h>

#include "c2xsf.h"

#define LINE_SIZE 300

void xyz_read(FILE* infile, struct unit_cell *c, struct contents *m){
  int i,coda;
  double *dptr;
  char buffer[LINE_SIZE+1],*ptr,*ptr2;

  if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
  if (sscanf(buffer,"%d",&m->n)!=1)
    error_exit("Unable to read number of atoms");
  if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");

  /* We insist on a lattice being present, either in a Lattice string
   * on this line, or at the end of the file with "%PBC" here
   */
  coda=0;
  if (!strncmp(buffer,"%PBC",4))
    coda=1;
  else{
    ptr=strstr(buffer,"Lattice=");
    if (!ptr)
      error_exit("XYZ file does not include lattice specification");
    ptr+=strlen("Lattice=")+1;
    c->basis=malloc(9*sizeof(double));
    if (!c->basis) error_exit("Malloc error for basis");
    dptr=(double*)c->basis;
    if (sscanf(ptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",dptr,dptr+1,
	       dptr+2,dptr+3,dptr+4,dptr+5,dptr+6,dptr+7,dptr+8)!=9)
      error_exit("Failed to parse lattice");
  }

  m->atoms=malloc(m->n*sizeof(struct atom));
  if (!m->atoms) error_exit("Malloc error for atoms");
  
  for(i=0;i<m->n;i++){
    if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
    ptr=buffer;
    while(isspace(*ptr)) ptr++;
    if (isdigit(*ptr)){
      if (sscanf(ptr,"%d %lf %lf %lf",&(m->atoms[i].atno),m->atoms[i].abs,
		 m->atoms[i].abs+1,m->atoms[i].abs+2)!=4)
	error_exit("Error parsing atom line");
    }
    else{
      ptr2=ptr;
      while((*ptr2)&&(!isspace(*ptr2))) ptr2++;
      *ptr2=0;
      m->atoms[i].atno=atsym2no(ptr);
      if (sscanf(ptr2+1,"%lf %lf %lf",m->atoms[i].abs,
		 m->atoms[i].abs+1,m->atoms[i].abs+2)!=3)
	error_exit("Error parsing atom line");
    }
  }
  
  if (coda){
    c->basis=malloc(9*sizeof(double));
    if (!c->basis) error_exit("Malloc error for basis");
    i=0;
    while(1){
      if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
      ptr=buffer;
      while(isspace(*ptr)) ptr++;
      if (!(strncasecmp(ptr,"vec",3))){
	ptr2=ptr;
	while((*ptr2)&&(!isspace(*ptr2))) ptr2++;
	if (sscanf(ptr2+1,"%lf %lf %lf",c->basis[i],c->basis[i]+1,
		   c->basis[i]+2)!=3)
	  error_exit("Error parsing lattice vector");
	i++;
	if (i==3) break;
      }
    }
  }

  real2rec(c);
  addfrac(m->atoms,m->n,c->recip);
  
}
  
  
