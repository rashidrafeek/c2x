/* Write a Gaussian cube file, one density only
 *       Units must be Bohr, not Angstoms
 * If no density given, as density is a required part of a Guassian
 * file, writes a 1x1x1 grid containing the value zero (in common with
 * some other programs).
 *
 * If ALT_OUT flag set, use molecular orbital form (but just one orbital)
 */


/* Copyright (c) 2007-2019 MJ Rutter 
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

void cube_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                struct grid *g){
  int i,j,k;
  double x,y,z,zero,*dptr1,*dptr2;
  char *fmt1,*fmt2,*fmt3;

  zero=0;
  
  if (!g->data){ /* Make dummy 1x1x1 grid with value zero */
    g->data=&zero;
    g->size[0]=g->size[1]=g->size[2]=1;
  }

  /* Cube files start with two lines of comments */
  
  if (m->title)
    fprintf(outfile,"%s\n",m->title);
  else
    fprintf(outfile,"\n");

  if (g->name)
    fprintf(outfile,"%s\n",g->name);
  else
    fprintf(outfile,"\n");

  if (flags&ALT_OUT)
    fprintf(outfile,"%d 0.0 0.0 0.0\n",-m->n);
  else
    fprintf(outfile,"%d 0.0 0.0 0.0\n",m->n);

  if (flags&HIPREC){
    fmt1="%d %.15f %.15f %.15f\n";
    fmt2="%d %.4f %.15f %.15f %.15f\n";
    fmt3="%22.15E";
  }
  else{
    fmt1="%d %f %f %f\n";
    fmt2="%d %.3f %f %f %f\n";
    fmt3="%12.5E";
  }
      
  
  for(i=0;i<3;i++) fprintf(outfile,fmt1,g->size[i],
                                c->basis[i][0]/g->size[i]/BOHR,
                                c->basis[i][1]/g->size[i]/BOHR,
                                c->basis[i][2]/g->size[i]/BOHR);

/* Need to write coords in Cartesian basis */
  for(i=0;i<m->n;i++){
    x=m->atoms[i].abs[0]/BOHR;
    y=m->atoms[i].abs[1]/BOHR;
    z=m->atoms[i].abs[2]/BOHR;
    fprintf(outfile,fmt2,m->atoms[i].atno,m->atoms[i].chg,x,y,z);
  }

  if (flags&ALT_OUT) fprintf(outfile," 1 1\n");
  
  dptr2=g->data;
  for(k=0;k<g->size[0];k++){
    for(j=0;j<g->size[1];j++){
      dptr1=dptr2+((k*g->size[1])+j)*g->size[2];
      for(i=0;i<g->size[2];i++){
        if (flags&AU){
          fprintf(outfile,fmt3,(*(dptr1+i))*BOHR*BOHR*BOHR);
        }
        else{
          fprintf(outfile,fmt3,*(dptr1+i));
        }
        if ((i%6)==5) fwrite("\n",1,1,outfile);
        else fwrite(" ",1,1,outfile);
      }
      if ((g->size[2]%6)!=0) fwrite("\n",1,1,outfile);
    }
  }
}
