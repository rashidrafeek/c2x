/* Write a shelx file. See http://shelx.uni-ac.gwdg.de/SHELX/shelx97.pdf */


/* Copyright (c) 2013 MJ Rutter 
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

extern int periodic_max_el;

void shelx_write(FILE* outfile, struct unit_cell *c, struct contents *m){
  int i,j,nspec;
  int *atom_present;
  int *atom_list;
  int *atom_spec;
  double abc[6];
  char cspec[4];

  char *fmt1,*fmt2;

  atom_present=calloc(periodic_max_el+1,sizeof(int));
  atom_list=calloc(periodic_max_el+1,sizeof(int));
  atom_spec=calloc(periodic_max_el+1,sizeof(int));

  if ((!atom_present)||(!atom_list)||(!atom_spec))
    error_exit("malloc error in shelx_write");
  
  if (flags&HIPREC){
    fmt1="CELL 1.0 %16.14f %16.14f %16.14f %16.14f %16.14f %16.14f\n";
    fmt2="%-4s %d %16.14f %16.14f %16.14f";

  }
  else{
    fmt1="CELL 1.0 %12.8f %12.8f %12.8f %8.4f %8.4f %8.4f\n";
    fmt2="%-4s %d %12.8f %12.8f %12.8f";
  }

  /* Find which atoms we have in our cell */

  for(i=0;i<m->n;i++)
    atom_present[m->atoms[i].atno]++;

  nspec=1;
  for(i=0;i<periodic_max_el;i++)
    if (atom_present[i]){
      atom_list[nspec]=i;
      atom_spec[i]=nspec;
      nspec++;
    }
      
  /* Write header */

  if (m->title)
    fprintf(outfile,"TITL %s\n",m->title);
  else
    fprintf(outfile,"TITL shelx file written by c2x\n");

  cart2abc(c,m,abc,NULL);
  fprintf(outfile,fmt1,
                  abc[0],abc[1],abc[2],abc[3],abc[4],abc[5]);

  if (!(flags&ALT_OUT)) fprintf(outfile,"ZERR 1 .0 .0 .0 0 0 0\n");
  fprintf(outfile,"LATT -1\n");

  fprintf(outfile,"SFAC");
  for(i=1;i<nspec;i++){
    strcpy(cspec,atno2sym(atom_list[i]));
    /* Make cspec upper case */
    for(j=0;cspec[j];j++) cspec[j]=toupper(cspec[j]);
    fprintf(outfile," %s",cspec);
  }

  if (!(flags&ALT_OUT)){
    fprintf(outfile,"\nUNIT");
    for(i=1;i<nspec;i++){
      fprintf(outfile," %d",atom_present[atom_list[i]]);
    }
  }
  fprintf(outfile,"\n");

  for(i=0;i<m->n;i++){
    fprintf(outfile,fmt2,atno2sym(m->atoms[i].atno),
	    atom_spec[m->atoms[i].atno],m->atoms[i].frac[0],m->atoms[i].frac[1],
	    m->atoms[i].frac[2]);
    if (flags&ALT_OUT) fprintf(outfile," 11\n");
    else fprintf(outfile,"\n");
  }

  fprintf(outfile,"END\n");

  free(atom_spec);
  free(atom_list);
  free(atom_present);
	   
}
