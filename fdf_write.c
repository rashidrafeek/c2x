/* Write an fdf Siesta file in a fashion that gdis understands 
 *
 * write a density file if there is a density to write
 * and a filename ending .fdf was specified
 * Density is written to a file ending .RHO
 *
 * Messy: work in progress...
 *
 * MJR 5/05
 */


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


#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "c2xsf.h"

void fdf_write(FILE* fdf, char* filename, struct unit_cell *c,
               struct contents *m, struct grid *g){
  int i,j,k,nspec,*natomsp,*spatno,len;
  double *dptr1;
  float *fptr;
  FILE *bin;
  char *cptr;

  fprintf(fdf,"# fdf file created by c2x\n");

  fprintf(fdf,"LatticeConstant  1.0 Ang\n");

  fprintf(fdf,"%%block LatticeVectors\n");
  for(i=0;i<3;i++)fprintf(fdf,"%f %f %f\n",c->basis[i][0],
        c->basis[i][1],c->basis[i][2]);
  fprintf(fdf,"%%endblock LatticeVectors\n\n");

  /* Now we need to know the number of species.
     It must be fewer than the number of atoms...
     This is horribly inefficient, but I intend to
     use it for ethene only... */

  nspec=0;
  natomsp=malloc(m->n*sizeof(int));
  if (!natomsp) error_exit("Malloc error in fdf_write");
  spatno=malloc(m->n*sizeof(int));
  if (!spatno) error_exit("Malloc error in fdf_write");


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

  fprintf(fdf,"NumberOfAtoms  %d\n",m->n);
  fprintf(fdf,"NumberOfSpecies  %d\n",nspec);

  fprintf(fdf,"%%block ChemicalSpeciesLabel\n");
  for(i=0;i<nspec;i++)
    fprintf(fdf,"%d  %d   %s\n",i+1,spatno[i],atno2sym(spatno[i]));
  fprintf(fdf,"%%endblock ChemicalSpeciesLabel\n\n");

  fprintf(fdf,"AtomicCoordinatesFormat Fractional\n");
    /* Write out atoms, sorted by species */
  fprintf(fdf,"%%block AtomicCoordinatesAndAtomicSpecies\n");
  for(i=0;i<nspec;i++)
    for(j=0;j<m->n;j++)
      if (m->atoms[j].atno==spatno[i])
        fprintf(fdf," %f %f %f %d\n",m->atoms[j].frac[0],
                m->atoms[j].frac[1],m->atoms[j].frac[2],i+1);
  fprintf(fdf,"%%endblock AtomicCoordinatesAndAtomicSpecies\n");


  /* And now write density */

  if (g->data==0) return;

  if (filename==NULL) return;

  if (strlen(filename)<4) return;

  cptr=filename+strlen(filename)-4;
  
  if (strcmp(cptr,".fdf")) return;

  cptr[1]='R';
  cptr[2]='H';
  cptr[3]='O';

  bin=fopen(filename,"wb");

  if(!bin){
    fprintf(stderr,"Error: unable to open %s for writing\n",filename);
    exit(1);
  }

  len=9*sizeof(double);
  fwrite(&len,sizeof(int),1,bin);
  fwrite(c->basis,sizeof(double),9,bin);
  fwrite(&len,sizeof(int),1,bin);

  len=4*sizeof(int);
  fwrite(&len,sizeof(int),1,bin);
  fwrite(g->size,sizeof(int),3,bin);
  i=1;
  fwrite(&i,sizeof(int),1,bin);
  fwrite(&len,sizeof(int),1,bin);

  if(!(fptr=malloc(g->size[0]*sizeof(float))))
      error_exit("Malloc error in fdf_write\n");

  len=g->size[0]*sizeof(float);
  for (k=0;k<g->size[2];k++)
    for (j=0;j<g->size[1];j++){
      dptr1=g->data+k+j*g->size[2];
      for (i=0;i<g->size[0];i++)
        fptr[i]=*(dptr1+i*g->size[1]*g->size[2]);
      fwrite(&len,sizeof(int),1,bin);
      fwrite(fptr,sizeof(float),g->size[0],bin);
      fwrite(&len,sizeof(int),1,bin);
    }

}



