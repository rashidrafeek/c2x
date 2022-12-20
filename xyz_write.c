/* Write an xyz format file */

/* Copyright (c) 2007, 2018 MJ Rutter 
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

/* Each frame (one here) has the format of
 * NATOMS
 * Comment line
 * Atomic symbol x y z (repeated natoms times)
 */

#include<stdio.h>
#include<stdlib.h>

#include "c2xsf.h"

void xyz_write(FILE* outfile, struct unit_cell *c, struct contents *m){
  int i;
  char *fmt;

  if (flags&HIPREC)
    fmt="%3s %19.14f %19.14f %19.14f\n";
  else
    fmt="%3s %12.8f %12.8f %12.8f\n";

  fprintf(outfile,"%d\n",m->n);

  if (m->title)
    fprintf(outfile,"%s\n",m->title);
  else
    fprintf(outfile,"\n");

  
  for(i=0;i<m->n;i++)
    fprintf(outfile,fmt,
                     atno2sym(m->atoms[i].atno),m->atoms[i].abs[0],
                     m->atoms[i].abs[1],m->atoms[i].abs[2]);
}
