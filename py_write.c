/* Produce output in the form of python's ASE Atoms class */

/* Copyright (c) 2017 MJ Rutter 
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

#include "c2xsf.h"

/* Type 'a' -- use ASE module explicitly
 *      'd' -- write dictionary compatible with
 *
 *               from ase import Atoms
 *               import foo
 *               st=Atoms(**foo.structure)
 */

void py_write(FILE* outfile, struct unit_cell *c, struct contents *m,
              char type){
  int i;
  char *fmt1,*fmt2,*fmt3;
  
  if(flags&HIPREC){
    fmt1="[(%.15f,%.15f,%.15f),\n"
      "                     (%.15f,%.15f,%.15f),\n"
      "                     (%.15f,%.15f,%.15f)],\n";
    fmt2="    (%.15f,%.15f,%.15f),\n";
    fmt3="    (%.15f,%.15f,%.15f)],\n";
  }
  else{
    fmt1="[(%f,%f,%f),\n"
      "                     (%f,%f,%f),\n"
      "                     (%f,%f,%f)],\n";
    fmt2="    (%f,%f,%f),\n";
    fmt3="    (%f,%f,%f)],\n";
  }
  
  if (m->title)
    fprintf(outfile,"# %s\n",m->title);
  
  if (type=='a') fprintf(outfile,"from ase import Atoms\n");

  if (type=='a') fprintf(outfile,"structure = Atoms(cell=");
  else if (type=='d') fprintf(outfile,"structure = {'cell':");
  
  fprintf(outfile,fmt1,
	  c->basis[0][0],c->basis[0][1],c->basis[0][2],
  	  c->basis[1][0],c->basis[1][1],c->basis[1][2],
	  c->basis[2][0],c->basis[2][1],c->basis[2][2]);


  if (type=='a') fprintf(outfile,"  symbols=[");
  else if (type=='d') fprintf(outfile,"  'symbols':[");
  
  for(i=0;i<m->n-1;i++)
    fprintf(outfile,"'%s',",atno2sym(m->atoms[i].atno));

  fprintf(outfile,"'%s'],\n",atno2sym(m->atoms[m->n-1].atno));

  if (type=='a') fprintf(outfile,"  scaled_positions=[\n");
  else if (type=='d') fprintf(outfile,"  'scaled_positions':[\n");
  for(i=0;i<m->n-1;i++)
    fprintf(outfile,fmt2,m->atoms[i].frac[0],
	    m->atoms[i].frac[1],m->atoms[i].frac[2]);
  i=m->n-1;
  fprintf(outfile,fmt3,m->atoms[i].frac[0],
	  m->atoms[i].frac[1],m->atoms[i].frac[2]);
  
  if (type=='a') fprintf(outfile,"  pbc=True)\n");
  else if (type=='d') fprintf(outfile,"  'pbc':True}\n");
}
