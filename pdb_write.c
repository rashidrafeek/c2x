/* Write a PDB file */

/* Copyright (c) 2007, 2019 MJ Rutter 
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

/* Changed 4/07 to start atom numbering at one, not zero
 * and to put atomic symbols in the correct place */

/* Excerpted from PDB file format (www.pdb.org):
 *
 * Columns    Thing
 *
 *  1-6       "ATOM  "
 *  7-11      Atom serial no (int)
 * 13-16      Atom name (lable)
 *  ...
 * 31-38      x co-ord, as Fortran's real(8.3), Angstroms
 * 39-46      y co-ord, ditto
 * 47-54      z co-ord, ditto
 *  ...
 * 77-78      Element symbol, right justified
 *
 *  1-6       "CRYST1"
 *  7-15      a, as Fortran's real(9.3), Angstroms
 * 16-24      b, ditto
 * 25-33      c, ditto
 * 34-40      alpha, as Fortran's real(7.2), degrees
 * 41-47      beta, ditto
 * 48-54      gamma, ditto
 * 56-66      space group, left justified string
 * 67-70      number of polmeric chains per cell (z)
 *
 * Set z to one, and space group to "P 1" if irrelevant / undetermined
 */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

#include "c2xsf.h"

extern int periodic_max_el;

void pdb_write(FILE* outfile, struct unit_cell *c, struct contents *m){
  int i,*n_in_el=NULL;
  double abc[6];

  if (m->n>99999)
    error_exit("Cannot write pdb file with more than 99,999 atoms");
  
  if (flags&ALT_OUT){
    n_in_el=calloc(periodic_max_el+1,sizeof(int));
    if (!n_in_el) error_exit("Calloc error in pdbn_write");
  }

  fprintf(outfile,"REMARK written by c2x\n");

  cart2abc(c,m,abc,NULL);
  fprintf(outfile,"CRYST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2f P 1"
                  "           1\n",
                  abc[0],abc[1],abc[2],abc[3],abc[4],abc[5]);

/* Ditch chain identifier (X)? */

  for(i=0;i<m->n;i++)
    if (flags&ALT_OUT){
      fprintf(outfile,"ATOM  %5d ",i+1);
      if(++n_in_el[min(m->atoms[i].atno,periodic_max_el)]<=99)
          fprintf(outfile,"%2s%02d",atno2sym(m->atoms[i].atno),
                  n_in_el[min(m->atoms[i].atno,periodic_max_el)]);
      else if((n_in_el[min(m->atoms[i].atno,periodic_max_el)]<=999)&&
              (strlen(atno2sym(m->atoms[i].atno))==1))
          fprintf(outfile,"%1s%d",atno2sym(m->atoms[i].atno),
                  n_in_el[min(m->atoms[i].atno,periodic_max_el)]);
      else fprintf(outfile,"%2s00",atno2sym(m->atoms[i].atno));
      fprintf(outfile,"   X     1     %7.3f %7.3f %7.3f"
                   "                      %2s\n",
                     m->atoms[i].abs[0],
                     m->atoms[i].abs[1],m->atoms[i].abs[2],atno2sym(m->atoms[i].atno));
    }
    else fprintf(outfile,"ATOM  %5d %2s     X     1     %7.3f %7.3f %7.3f"
                   "                      %2s\n",
                     i+1,atno2sym(m->atoms[i].atno),m->atoms[i].abs[0],
                     m->atoms[i].abs[1],m->atoms[i].abs[2],atno2sym(m->atoms[i].atno));

  fprintf(outfile,"END\n");
  if (flags&ALT_OUT) free(n_in_el);
}
