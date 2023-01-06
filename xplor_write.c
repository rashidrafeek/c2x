/* Write an xplor file, charge density only */

/* Copyright (c) 2007 MJ Rutter 
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


/* Some programs (pymol) are very fussy about the precise format statements
 *
 */

/* The xplor file format is:
 *
 * Blank line
 * I8      ntitles (must be >0)
 * ntitles lines of text as titles
 * 9I8     na,amin,amax,nb,bmin,bmax,nc,cmin,cmax
 * 6E12.5  a,b,c,alpha,beta,gamma  (Angstroms and degrees)
 * ZYX     precisely those three characters
 * do  c=cmin,cmax
 * I8      section number (1..nc ?)
 * 6E12.5  ((map(a,b,c),a=amin,amax),b=bmin,bmax)
 * enddo
 *
 * The stupid thing has a grid offset of half a grid cell compared to
 * anyone else's idea of sanity. So we have an option to shift all
 * atoms by this amount, rather than using inexact grid interpolation
 * schemes
 */

#include<stdio.h>
#include<stdlib.h>

#include "c2xsf.h"

void xplor_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                 struct grid *g){
  int i,j,k;
  double *dptr1,*dptr2,abc[6];

  if(!g->data) error_exit("Xplor output requested, but no grid data");

  fprintf(outfile,"\n      1 !NTITLE\n%s\n",g->name);

  fprintf(outfile," %7d %7d %7d %7d %7d %7d %7d %7d %7d\n",
                                                 g->size[0],0,g->size[0]-1,
                                                 g->size[1],0,g->size[1]-1,
                                                 g->size[2],0,g->size[2]-1);

  cart2abc(c,m,abc,NULL);
  fprintf(outfile," %11f %11f %11f %11f %11f %11f\n",abc[0],abc[1],abc[2],
                                                     abc[3],abc[4],abc[5]);

  fprintf(outfile,"ZYX\n");

  dptr2=g->data;
  for(k=0;k<g->size[2];k++){
    fprintf(outfile,"%d\n",k+1);
    for(j=0;j<g->size[1];j++){
      dptr1=dptr2+k+j*g->size[2];              ;
      for(i=0;i<g->size[0];i++)
        fprintf(outfile,"%f\n",*(dptr1+i*g->size[1]*g->size[2]));
    }
  }
}

void xplor_fudge(struct unit_cell *c, struct contents *m, struct grid *g){
  double sx,sy,sz;
  int i;

  if((g->size[0]==0)||(g->size[1]==0)||(g->size[2]==0))
    error_exit("Shift of half grid cell requested, but no grid found.");

  if (debug>0) fprintf(stderr,"Shifting by half grid cell as requested\n");

  sx=0.5/g->size[0];
  sy=0.5/g->size[1];
  sz=0.5/g->size[2];

  for(i=0;i<m->n;i++){
    m->atoms[i].frac[0]+=sx;
    m->atoms[i].frac[1]+=sy;
    m->atoms[i].frac[2]+=sz;
  }

  addabs(m->atoms,m->n,c->basis);
}
