/* Write an xyz or extended xyz format file */

/* Copyright (c) 2007, 2018-2020 MJ Rutter 
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

/* Each frame has the format of
 * NATOMS
 * Comment line
 * Atomic symbol x y z (repeated natoms times)
 *
 * if ALT_OUT set, write an extended xyz format
 */

#include<stdio.h>
#include<stdlib.h>

#include "c2xsf.h"

void xyz_write(FILE* outfile, struct unit_cell *c, struct contents *m,
	       struct time_series *ts){
  int i;
  char *fmt;

  if ((ts)&&(ts->nsteps>1)){
    if ((ts->nm<ts->nsteps)||(ts->nc<ts->nsteps)){
      error_exit("Too few steps for atoms or cell. Aborting.\n");
    }
    if (!ts->m[0].title) ts->m[0].title=m->title;
    if (!ts->m[0].dict) ts->m[0].dict=m->dict;
    for(i=0;i<ts->nsteps;i++)
      xyz_write(outfile,ts->cells+i,ts->m+i,NULL);
    return;
  }

  
  if (flags&HIPREC)
    fmt="%3s % 19.14f % 19.14f % 19.14f\n";
  else
    fmt="%3s % 14.8f % 14.8f % 14.8f\n";

  fprintf(outfile,"%d\n",m->n);

  if (flags&ALT_OUT){
#if 0
    fprintf(outfile,"%%PBC\n");
#else
    fprintf(outfile,"Lattice=\"");
    for (i=0;i<3;i++)
      fprintf(outfile,"%.8f %.8f %.8f%s",c->basis[i][0],c->basis[i][1],
              c->basis[i][2],(i==2)?"":" ");
    fprintf(outfile,"\" Properties=species:S:1:pos:R:3\n");
#endif
  }
  else if (m->title)
    fprintf(outfile,"%s\n",m->title);
  else{
    if (dict_get(m->dict,"in_file"))
      fprintf(outfile,"Converted from %s by c2x\n",
	      (char*)dict_get(m->dict,"in_file"));
    else
      fprintf(outfile,"\n");
  }
  
  for(i=0;i<m->n;i++)
    fprintf(outfile,fmt,
                     atno2sym(m->atoms[i].atno),m->atoms[i].abs[0],
                     m->atoms[i].abs[1],m->atoms[i].abs[2]);

#if 0  
  if (flags&ALT_OUT){
    fprintf(outfile,"\n");
    for(i=0;i<3;i++){
      fprintf(outfile,"Vector%d",i+1);
      fprintf(outfile,fmt,"",c->basis[i][0],c->basis[i][1],c->basis[i][2]);
    }
    fprintf(outfile,"Offset ");
    fprintf(outfile,fmt,"",0.0,0.0,0.0);
  }
#endif
}
