/* Read data set from a .cube file */

/* Copyright (c) 2014, 2018 MJ Rutter 
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

#include "c2xsf.h"

#define LINE_SIZE 100

void cube_read(FILE* infile, struct unit_cell *c, struct contents *m,
               struct grid *gptr){
  static char buffer[LINE_SIZE+1];
  char *ptr;
  int i,j,n;
  double x,y,z,t[3],scale;

  if (debug>2) fprintf(stderr,"cube_read called\n");

  fgets(buffer,LINE_SIZE,infile);
  /* Was the title more than spaces and an end-of-line? */
  ptr=buffer;
  while((*ptr==' ')||(*ptr=='\n')) ptr++;
  if (*ptr){
    m->title=malloc(strlen(buffer)+1);
    if (!m->title) error_exit("Malloc error for title in cube_read()");
    strcpy(m->title,buffer);
    /* Remove newline */
    m->title[strlen(m->title)-1]=0;
  }

  fgets(buffer,LINE_SIZE,infile);

  fgets(buffer,LINE_SIZE,infile);
  if (sscanf(buffer,"%d %lf %lf %lf",&n,&x,&y,&z)!=4)
    error_exit("Scan error 3rd line");

  if ((n<=0)||(x!=0)||(y!=0)||(z!=0))
    error_exit("Unexpected data 3rd line");

  if (!(c->basis=malloc(72)))
    error_exit("Malloc error in cube_read for c->basis");

  gptr=grid_new(gptr);
  
  for(i=0;i<3;i++){
    fgets(buffer,LINE_SIZE,infile);
    if (sscanf(buffer,"%d %lf %lf %lf",gptr->size+i,t,t+1,t+2)!=4)
      error_exit("Scan error in grid description");
    for(j=0;j<3;j++) c->basis[i][j]=t[j]*gptr->size[i]*BOHR;
  }

  if(!(m->atoms=malloc(n*sizeof(struct atom))))
    error_exit("Malloc error in cube_read for atoms");
  m->n=n;
  init_atoms(m->atoms,n);
  
  for(i=0;i<n;i++){
    fgets(buffer,LINE_SIZE,infile);
    if (sscanf(buffer,"%d %lf %lf %lf %lf",&(m->atoms[i].atno),
               &(m->atoms[i].chg),t,t+1,t+2)!=5)
      error_exit("Scan error in atom description");
    for(j=0;j<3;j++) m->atoms[i].abs[j]=t[j]*BOHR;
  }

  real2rec(c);
  addfrac(m->atoms,m->n,c->recip);

  if(!(gptr->data=malloc(gptr->size[0]*gptr->size[1]*
			 gptr->size[2]*sizeof(double))))
    error_exit("Malloc error in cube_read for grid data");

  n=0;
  scale=1/(BOHR*BOHR*BOHR);
  for(i=0;i<gptr->size[0]*gptr->size[1]*gptr->size[2];i++){
    if (fscanf(infile,"%lf",gptr->data+i)!=1){n=1;break;}
    if (flags&DE_AU) gptr->data[i]*=scale;
  }

  if (n) error_exit("Read error for grid data");

  /* Discard volumetic data if it is a dummy single point of zero set */
  if ((gptr->size[0]==1)&&(gptr->size[1]==1)&&(gptr->size[2]==1)&&
      (*gptr->data==0)){
    free(gptr->data);
  }
  else
    gptr->name="Density";

}
