/* Write a CASTEP den_fmt file
 * will write grids until the first change in grid size
 */


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
#include<stdlib.h>
#include<string.h>

#include "c2xsf.h"

void denfmt_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                struct grid *g){
  double abc[6];
  unsigned int i,j,k,ii,offset,ngrids,nfft[3],same;
  char *fmt;
  struct grid *gptr;

  if (flags&HIPREC){
    fmt="  %.15f";
  }
  else{
    fmt="  % 12.6f";
  }
  
  cart2abc(c,NULL,abc,NULL,0);
  
  fprintf(outfile," BEGIN header\n \n");

  fprintf(outfile,"           Real Lattice(A)               "
          "Lattice parameters(A)    Cell Angles\n");

  fprintf(outfile," %11.7f %11.7f %11.7f     a = %11.6f  alpha = %11.6f\n",
          c->basis[0][0],c->basis[0][1],c->basis[0][2],abc[0],abc[3]);
  fprintf(outfile," %11.7f %11.7f %11.7f     b = %11.6f  beta  = %11.6f\n",
          c->basis[1][0],c->basis[1][1],c->basis[1][2],abc[1],abc[4]);
  fprintf(outfile," %11.7f %11.7f %11.7f     c = %11.6f  gamma = %11.6f\n",
          c->basis[2][0],c->basis[2][1],c->basis[2][2],abc[2],abc[5]);

  fprintf(outfile," \n");

  ngrids=0;
  gptr=g;
  for(i=0;i<3;i++) nfft[i]=g->size[i];
  while(gptr){
    if (gptr->data){
      same=1;
      for(i=0;i<3;i++) if (gptr->size[i]!=nfft[i]) same=0;
      if (same) ngrids++;
      else break;
    }
    gptr=gptr->next;
  }
  if (!same) fprintf(stderr,"Some grids missed as different sizes present\n");
/* Jmol detects den_fmt files by the existence of the string "! nspins",
 * so we cannot use the more accurate comment "! ngrids"
 */
/*  fprintf(outfile,"   %d             ! ngrids\n",ngrids); */
  fprintf(outfile,"   %d             ! nspins (really ngrids)\n",ngrids);
  fprintf(outfile,"  %d   %d   %d   ! grid size along <a,b,c>\n",nfft[0],
          nfft[1],nfft[2]);
  fprintf(outfile," END header: data is \"<a b c>");

  if (debug) fprintf(stderr,"Writing %d grids:",ngrids);
  gptr=g;
  for (i=0;i<ngrids;i++){
    if (gptr->name){
      if (!strncasecmp(gptr->name,"Density",7)){
        fprintf(outfile," charge");
        if (debug) fprintf(stderr," charge");
      }
      else if (!strncasecmp(gptr->name,"Spin",4)){
        fprintf(outfile," spin");
        if (debug) fprintf(stderr," spin");
      }
      else{
        fprintf(outfile," (%s)",gptr->name);
        if (debug) fprintf(stderr," (%s)",gptr->name);
      }
    }
    else{
        fprintf(outfile," (unknown)");
        if (debug) fprintf(stderr," (unknown)");
    }
    gptr=gptr->next;
  }
  fprintf(outfile,"\"\n");
  if (debug) fprintf(stderr,"\n");
  fprintf(outfile," \n");

  /* Now write the data */

  for(i=1;i<=nfft[2];i++)
    for(j=1;j<=nfft[1];j++)
      for(k=1;k<=nfft[0];k++){
        fprintf(outfile,"  %4d  %4d  %4d  ",k,j,i);
        gptr=g;
        offset=(i-1)+(j-1)*nfft[2]+(k-1)*nfft[2]*nfft[1];
        for(ii=0;ii<ngrids;ii++){
          fprintf(outfile,fmt,gptr->data[offset]*c->vol);
          gptr=gptr->next;
        }
        fprintf(outfile,"\n");
      }
}
