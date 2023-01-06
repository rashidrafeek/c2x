/* Reader for a Siesta .RHO file */

/* Copyright (c) 2019 MJ Rutter 
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
#include<stdlib.h> /* malloc */
#include<string.h>
#include<math.h>

#include "c2xsf.h"


void rho_read(FILE* infile, struct unit_cell *c, struct contents *m,
              struct grid *gptr, struct es *elect){
  int tmp,nspins,i,j,ix,iy,iz,grid[3],okay,csize,single;
  int is_pot,junk;
  void *column;
  double scale,basis[3][3],*dptr1,*dptr2;
  char *filename,*name,*cptr;
  FILE *coords;

  dptr1=dptr2=NULL;
  scale=1;
  is_pot=0;
  filename=dict_get(m->dict,"in_file");
  
  tmp=0;
  fread(&tmp,4,1,infile);
  if (tmp!=72) error_exit("Unexpected start to RHO file");

  fread(basis,72,1,infile);
  fread(&junk,4,1,infile);

  fread(&tmp,4,1,infile);
  if (tmp!=16) error_exit("Unexpected second record length in RHO file");
  
  if ((!strcasecmp(filename+strlen(filename)-4,".rho"))||
      (!strcasecmp(filename+strlen(filename)-6,".rhoxc"))){
    if (!(flags&RAW)){
      if (flags&AU){
        add_cmt(m->comment,"Densities in e Bohr^-3");
      }
      else{
        scale=1/(BOHR*BOHR*BOHR);
        add_cmt(m->comment,"Densities in e A^-3");
      }
    }
  }
  
  if (!strcasecmp(filename+strlen(filename)-5,".toch")){
    if (!(flags&RAW)){
      if (flags&AU){
        scale=-1;
        add_cmt(m->comment,"Densities in e Bohr^-3");
      }
      else{
        scale=-1/(BOHR*BOHR*BOHR);
        add_cmt(m->comment,"Densities in e A^-3");
      }
    }
  }

  if ((!strcasecmp(filename+strlen(filename)-3,".vh"))||
      (!strcasecmp(filename+strlen(filename)-3,".vt"))){
    is_pot=1;
    if (!(flags&RAW)){
      scale=0.5*H_eV;
      add_cmt(m->comment,"Potential in V");
    }
  }
  
  fread(grid,12,1,infile);
  nspins=0;
  fread(&nspins,4,1,infile);
  fread(&junk,4,1,infile);

  if (debug) {
    if (is_pot)
      fprintf(stderr,"Grid size %dx%dx%d\n",
                     grid[0],grid[1],grid[2]);
    else
      fprintf(stderr,"Grid size %dx%dx%d spins=%d\n",
                     grid[0],grid[1],grid[2],nspins);
  }

  if ((flags&SPINDEN)&&(nspins!=2))
    fprintf(stderr,"Spin density requested, but nspins!=2. Ignored.\n");

  if (nspins>2)
    fprintf(stderr,"Confused: nspins>2. Rubbish may result.\n");
  
  csize=0;
  fread(&csize,4,1,infile);
  fseek(infile,-4,SEEK_CUR);
  if ((csize!=4*grid[0])&&(csize!=8*grid[0]))
        error_exit("Unexpected column size");
  single=1;
  if (csize==8*grid[0]) single=0;

  if (single)
    fprintf(stderr,"Recorded data are single precision\n");
  else
    fprintf(stderr,"Recorded data are double precision\n");
  
  column=malloc(csize);
  if (!column) error_exit("Malloc error for column");

  dptr1=malloc(grid[0]*grid[1]*grid[2]*sizeof(double));
  if (!dptr1) error_exit("Malloc error for grid data");
  
  for(iz=0;iz<grid[2];iz++){
    for(iy=0;iy<grid[1];iy++){
      fread(&tmp,4,1,infile);
      if (tmp!=csize)
        error_exit("Unexpected column size");
      fread(column,csize,1,infile);
      fread(&junk,4,1,infile);
      if (single)
        for(ix=0;ix<grid[0];ix++)
          dptr1[iz+grid[2]*(iy+ix*grid[1])]=scale*((float*)column)[ix];
      else
        for(ix=0;ix<grid[0];ix++)
          dptr1[iz+grid[2]*(iy+ix*grid[1])]=scale*((double*)column)[ix];
    }
  }

  if (nspins==2){
    dptr2=malloc(grid[0]*grid[1]*grid[2]*sizeof(double));
    if (!dptr2) error_exit("Malloc error for grid data");

    for(iz=0;iz<grid[2];iz++){
      for(iy=0;iy<grid[1];iy++){
        fread(&tmp,4,1,infile);
        if (tmp!=csize)
          error_exit("Unexpected column size");
        fread(column,csize,1,infile);
	fread(&junk,4,1,infile);
      if (single)
        for(ix=0;ix<grid[0];ix++)
          dptr2[iz+grid[2]*(iy+ix*grid[1])]=scale*((float*)column)[ix];
      else
        for(ix=0;ix<grid[0];ix++)
          dptr2[iz+grid[2]*(iy+ix*grid[1])]=scale*((double*)column)[ix];
      }
    }
  }

  free(column);
  
  if (is_pot){
    gptr=grid_new(gptr);
    gptr->name="Potential";
    for(i=0;i<3;i++) gptr->size[i]=grid[i];
    gptr->data=dptr1;
  }
  else{
    if (flags&CHDEN){
      gptr=grid_new(gptr);
      gptr->name=NULL;
      for(i=0;i<3;i++) gptr->size[i]=grid[i];
      gptr->data=dptr1;
      if (!strcasecmp(filename+strlen(filename)-4,".rho")){
        gptr->name="Density";
        if (nspins==2)
          for(i=0;i<grid[0]*grid[1]*grid[2];i++)
            dptr1[i]+=dptr2[i];
      }
    }

    if ((flags&SPINDEN)&&(nspins==2)){
      if (gptr->next) gptr=gptr->next;
      gptr->next=malloc(sizeof(struct grid));
      if (!gptr->next) error_exit("Malloc error for struct grid");
      gptr->next->next=NULL;
      gptr->next->data=NULL;
      gptr->name=NULL;
      for(i=0;i<3;i++) gptr->size[i]=grid[i];
      gptr->name="Spin";
      gptr->data=dptr2;
      if (flags&CHDEN)
        for(i=0;i<grid[0]*grid[1]*grid[2];i++)
          dptr2[i]=dptr1[i]-2*dptr2[i];
      else{
        for(i=0;i<grid[0]*grid[1]*grid[2];i++)
          dptr2[i]=dptr1[i]-dptr2[i];
        free(dptr1);
        dptr1=NULL;
      }
    }
    else if (dptr2){
      free(dptr2);
      dptr2=NULL;
    }
  }

  /* Must now read .fdf first, for atomic charges, and then .XV second,
   * for updated positions
   */
  
  name=malloc(strlen(filename)+5);
  if (!name) error_exit("Malloc error for name");
  strcpy(name,filename);
  cptr=name+strlen(name);
  while ((cptr>=name)&&(*cptr!='.')) cptr--;
  cptr++;

  *cptr='f';
  *(cptr+1)='d';
  *(cptr+2)='f';
  *(cptr+3)=0;
  coords=fopen(name,"rb");
  if (debug) fprintf(stderr,"Attempting to read %s %s\n",name,
                     (coords)?"":"...failed");
  if (coords){
    fprintf(stderr,"Also reading %s\n",name);
    fdf_read(coords,c,m,NULL,elect);
  }

  *cptr='X';
  *(cptr+1)='V';
  *(cptr+2)=0;
  coords=fopen(name,"rb");
  if (debug) fprintf(stderr,"Attempting to read %s %s\n",name,
                     (coords)?"":"...failed");
  if (coords){
    fprintf(stderr,"Also reading %s\n",name);
    xv_read(coords,c,m);
  }
  else {
    *cptr='x';
    *(cptr+1)='v';
    coords=fopen(name,"rb");
    if (debug) fprintf(stderr,"Attempting to read %s %s\n",name,
                       (coords)?"":"...failed");
    if (coords){
      fprintf(stderr,"Also reading %s\n",name);
      xv_read(coords,c,m);
    }
  }
  free(name);

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      basis[i][j]*=BOHR;

  if (c->basis){
    okay=1;
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
        if (!aeq(c->basis[i][j],basis[i][j])) okay=0;
    if (!okay)
      fprintf(stderr,"Warning: basis vectors in the two files differ!\n");
  }
  if (!c->basis) c->basis=malloc(9*sizeof(double));
  if (!c->basis) error_exit("Malloc error for basis");
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[i][j]=basis[i][j];
  real2rec(c);

  elect->nspins=nspins;
  
}
