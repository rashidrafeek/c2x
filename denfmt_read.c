/* Read a single data set from a CASTEP .den_fmt file 
 *
 * The format seems to be undocumented,
 * so it is hard to work out how it should be parsed...
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
#include<stdlib.h> /* malloc */

#include "c2xsf.h"

#define LINE_SIZE 100

void denfmt_read(FILE* infile, struct unit_cell *c,
               struct grid *gptr){
  static char buffer[LINE_SIZE+1];
  int i,nx,ny,nz,ngx,ngy,ngz;
  double x,scale;

  if (debug>2) fprintf(stderr,"denfmt_read called\n");
  
  fgets(buffer,LINE_SIZE,infile);
  fgets(buffer,LINE_SIZE,infile);
  fgets(buffer,LINE_SIZE,infile);

  if (!(c->basis=malloc(72)))
    error_exit("Malloc error in denfmt_read for c->basis");
 
  for(i=0;i<3;i++){
    fgets(buffer,LINE_SIZE,infile);
    if (sscanf(buffer,"%lf %lf %lf",&(c->basis[i][0]),
               &(c->basis[i][1]),&(c->basis[i][2]))!=3)
      error_exit("Scan error in lattice description");
  }
  real2rec(c);
   
  fgets(buffer,LINE_SIZE,infile);
  fgets(buffer,LINE_SIZE,infile);

  if (gptr->next) gptr=gptr->next;
  gptr->next=malloc(sizeof(struct grid));
  if (!gptr->next) error_exit("Malloc error for struct grid");
  gptr->next->next=NULL;
  gptr->next->data=NULL;

  fgets(buffer,LINE_SIZE,infile);
  if (sscanf(buffer,"%d %d %d",&ngx,&ngy,&ngz)!=3)
    error_exit("Scan error in grid dimensions");

  gptr->size[0]=ngx;
  gptr->size[1]=ngy;
  gptr->size[2]=ngz;

  if (debug>1) fprintf(stderr,"denfmt read grid size to be %dx%dx%d\n",
                       ngx,ngy,ngz);
  
  if(!(gptr->data=malloc(ngx*ngy*ngz*sizeof(double))))
    error_exit("Malloc error in denfmt_read for grid data");

  fgets(buffer,LINE_SIZE,infile);
  fgets(buffer,LINE_SIZE,infile);

  while(fgets(buffer,LINE_SIZE,infile)){
    if(sscanf(buffer,"%d %d %d",&nx,&ny,&nz)==3) break;
  }
  
  scale=1/c->vol;
  for(i=0;i<ngx*ngy*ngz;i++){
    if (i) if (!fgets(buffer,LINE_SIZE,infile)) break;
    if (sscanf(buffer,"%d %d %d %lf",&nx,&ny,&nz,&x)!=4){
      fprintf(stderr,"Parse error for data item %d\n",i);
      exit(1);
    }
    nx--;ny--;nz--;
    gptr->data[nx*ngy*ngz+ny*ngz+nz]=x*scale;
  }

  if (i!=ngx*ngy*ngz) error_exit("Too few data points read in denfmt_read");
  
  gptr->name="Unknown";

}
