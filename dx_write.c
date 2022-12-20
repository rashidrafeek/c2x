/* Write an dx file, single density only. Updated 2020 to make VMD-compatible
 */


/* Copyright (c) 2007, 2020 MJ Rutter 
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

#include "c2xsf.h"

void dx_write(FILE* outfile,struct unit_cell *c, struct grid *g){
  int i,j,k;
  double *dptr1,*dptr2;

  fprintf(outfile,"# %s\n\n",g->name);

  fprintf(outfile,"object 1 class gridpositions counts %d %d %d\n",
                   g->size[0],g->size[1],g->size[2]);
  if (!g->origin_abs)
    fprintf(outfile," origin 0.0 0.0 0.0\n");
  else
    fprintf(outfile," origin %f %f %f\n",
	    g->origin_abs[0],g->origin_abs[1],g->origin_abs[2]);
  for(i=0;i<3;i++)
   fprintf(outfile," delta %f %f %f\n",c->basis[i][0]/g->size[i],
                                       c->basis[i][1]/g->size[i],
                                       c->basis[i][2]/g->size[i]);

  fprintf(outfile,"object 2 class gridconnections counts %d %d %d\n",
                  g->size[0],g->size[1],g->size[2]);

  fprintf(outfile,"object 3 class array items %d data follows\n",
                  g->size[0]*g->size[1]*g->size[2]);

  dptr2=g->data;
  for(k=0;k<g->size[0];k++){
    for(j=0;j<g->size[1];j++){
      dptr1=dptr2+((k*g->size[1])+j)*g->size[2];
      for(i=0;i<g->size[2];i++)
        fprintf(outfile,"%g\n",*(dptr1+i));
    }
  }

  fprintf(outfile,"object \"%s\" class field\n",g->name);

}
