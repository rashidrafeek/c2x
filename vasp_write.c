/* Write a VASP-style output file, one density only
 */

/* Copyright (c) 2007, 2017 MJ Rutter 
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
#include<math.h>
#include<string.h>

#include "c2xsf.h"

void vasp_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                struct grid *g){
  int i,j,k,count,line;
  double *dptr1,scale,x;
  int nspec,*natomsp,*spatno;
  char *fmt;
  
  if (m->title) fprintf(outfile,"%s\n",m->title);
  else fprintf(outfile,"\n");
  
  fprintf(outfile," %20.16f\n",1.0);  /* Scale factor */

  
  if (flags&HIPREC)
    fmt="% 19.15f % 19.15f % 19.15f\n";
  else
    /* Some VASP readers eschew free format, so are fussy about details */
    fmt="% 12.6f % 12.6f % 12.6f\n";

  for(i=0;i<3;i++) fprintf(outfile,fmt,
                                c->basis[i][0],
                                c->basis[i][1],
                                c->basis[i][2]);

  if(!(flags&HIPREC)) fmt="% 8.6f  % 8.6f  % 8.6f\n";

  /* Now we need to know the number of species.
     It must be fewer than the number of atoms...
  */

  nspec=0;
  natomsp=malloc(m->n*sizeof(int));
  if (!natomsp) error_exit("Malloc error in vasp_write");
  spatno=malloc(m->n*sizeof(int));
  if (!spatno) error_exit("Malloc error in vasp_write");


  for(i=0;i<m->n;i++){
    for(j=0;j<nspec;j++) if (m->atoms[i].atno==spatno[j]) break;
    if (j==nspec){  /* new species */
      spatno[j]=m->atoms[i].atno;
      natomsp[j]=1;
      nspec++;
    }else{          /* existing species */
      natomsp[j]++;
    }
  }

  /* VASP now puts atomic symbols here */
  for(i=0;i<nspec;i++) fprintf(outfile,"   %s",atno2sym(spatno[i]));
  fprintf(outfile,"\n");

  
  /* Must not have trailing space on this line... */
  for(i=0;i<nspec;i++) fprintf(outfile,"  %d",natomsp[i]);
  fprintf(outfile,"\n");

  fprintf(outfile,"Direct\n");  /* i.e. fractional co-ords */

  /* Write out atoms, sorted by species */
  for(i=0;i<nspec;i++)
    for(j=0;j<m->n;j++)
      if (m->atoms[j].atno==spatno[i])
        fprintf(outfile,fmt,m->atoms[j].frac[0],
                m->atoms[j].frac[1],m->atoms[j].frac[2]);

  /* And now write density */

  if (!g->data) return;

  fprintf(outfile,"\n %d %d %d\n",g->size[0],g->size[1],g->size[2]);

  if ((flags&HIPREC)||(flags&CHGCAR)){
    /* Some VASP readers eschew free format, so are fussy about details */
    fmt=" %17.11E";
    line=3;
  }
  else{
    fmt=" %g";
    line=7;
  }
  scale=fabs(c->vol);

  count=-1;
  for(k=0;k<g->size[2];k++){
    for(j=0;j<g->size[1];j++){
      dptr1=g->data+k+j*g->size[2];
      for(i=0;i<g->size[0];i++){
        x=*(dptr1+i*g->size[1]*g->size[2])*scale;
        if (fabs(x)>1e-99) fprintf(outfile,fmt,x);
        else fprintf(outfile,fmt,0.0);
        count++;
        if ((count&line)==line) fprintf(outfile,"\n");
      }
    }
  }
  if ((count&line)!=line) fprintf(outfile,"\n");

  /* Should we also write a spin density? */

  g=g->next;
  if ((!g)||(!g->data)) return;
  if ((!g->name)||(strncasecmp(g->name,"Spin",4))) return;

  /* Write out mag moments, same order as atoms */

  k=0;
  for(i=0;i<nspec;i++)
    for(j=0;j<m->n;j++)
      if (m->atoms[j].atno==spatno[i]){
        fprintf(outfile,"  %f",m->atoms[j].spin);
        k++;
        if (k%4==0) fprintf(outfile,"\n");
      }
  if (k%4) fprintf(outfile,"\n");
  
  /* And now write spin */

  fprintf(outfile," %d %d %d\n",g->size[0],g->size[1],g->size[2]);
  count=-1;
  for(k=0;k<g->size[2];k++){
    for(j=0;j<g->size[1];j++){
      dptr1=g->data+k+j*g->size[2];
      for(i=0;i<g->size[0];i++){
        x=*(dptr1+i*g->size[1]*g->size[2])*scale;
        if (fabs(x)>1e-99) fprintf(outfile,fmt,x);
        else fprintf(outfile,fmt,0.0);
        count++;
        if ((count&line)==line) fprintf(outfile,"\n");
      }
    }
  }
  if ((count&line)!=line) fprintf(outfile,"\n");

}
