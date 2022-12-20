/* Read a single data set from a CASTEP .den_fmt (or .pot_fmt) file 
 *
 * The format seems to be undocumented,
 * so it is hard to work out how it should be parsed...
 */

/* Copyright (c) 2017-2019 MJ Rutter 
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

#define LINE_SIZE 120

void denfmt_read(FILE* infile, struct unit_cell *c, struct contents *m,
                 struct grid *gptr, struct es *elect, int rescale){
  static char buffer[LINE_SIZE+1];
  char *spin_range,*filename;
  int i,nx,ny,nz,ngx,ngy,ngz,nspins;
  double x,x2,scale,*data2;
  struct grid *gptr2;

  if (debug>2) fprintf(stderr,"denfmt_read called\n");

  data2=NULL;
  spin_range=elect->spin_range;

  filename=dict_get(m->dict,"in_file");
  i=strlen(filename);
  if ((i>7)&&(!strcmp(filename+i-7,"den_fmt"))){
    if ((flags&CHDEN)&&(flags&SPINDEN)) spin_range="-";
   else if (flags&SPINDEN) spin_range="1";
    else spin_range="0";
  }
  
  while(fgets(buffer,LINE_SIZE,infile))
    if (strstr(buffer,"BEGIN header")) break;

  if (!strstr(buffer,"BEGIN header")) error_exit("BEGIN header not found");
  
  while(fgets(buffer,LINE_SIZE,infile))
    if (strstr(buffer,"Real Lattice")) break;

  if (!strstr(buffer,"Real Lattice")) error_exit("Real Lattice not found");
  
  if (!(c->basis=malloc(9*sizeof(double))))
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
  if (sscanf(buffer,"%d",&nspins)!=1)
    error_exit("Scan error for nspins");

  if ((nspins!=1)&&(nspins!=2)){
    fprintf(stderr,"Unexpected value for nspins, assuming one\n");
    nspins=1;
  }
  
  fgets(buffer,LINE_SIZE,infile);
  if (sscanf(buffer,"%d %d %d",&ngx,&ngy,&ngz)!=3)
    error_exit("Scan error in grid dimensions");

  if (debug>1) fprintf(stderr,"denfmt read grid size to be %dx%dx%d\n",
                       ngx,ngy,ngz);

  if (gptr->next) gptr=gptr->next;
  gptr->next=malloc(sizeof(struct grid));
  if (!gptr->next) error_exit("Malloc error for struct grid");
  gptr->next->next=NULL;
  gptr->next->data=NULL;
  gptr->size[0]=ngx;
  gptr->size[1]=ngy;
  gptr->size[2]=ngz;

  if(!(gptr->data=malloc(ngx*ngy*ngz*sizeof(double))))
    error_exit("Malloc error in denfmt_read for grid data");

  if (nspins==2){
    if(!(data2=malloc(ngx*ngy*ngz*sizeof(double))))
      error_exit("Malloc error in denfmt_read for grid data");
  }
  
  while(fgets(buffer,LINE_SIZE,infile))
    if (strstr(buffer,"END header")) break;

  if (!strstr(buffer,"END header")) error_exit("END header not found");

  fgets(buffer,LINE_SIZE,infile);

  while(fgets(buffer,LINE_SIZE,infile)){
    if(sscanf(buffer,"%d %d %d",&nx,&ny,&nz)==3) break;
  }

  scale=1;
  if (!(flags&RAW)){
    if (rescale==1)
      scale=1/c->vol;
    else if (rescale==2)
      scale=H_eV;
  }
  
  for(i=0;i<ngx*ngy*ngz;i++){
    if (i) if (!fgets(buffer,LINE_SIZE,infile)) break;
    if (nspins==2){
      if (sscanf(buffer,"%d %d %d %lf %lf",&nx,&ny,&nz,&x,&x2)!=5){
	fprintf(stderr,"Parse error for data item %d\n",i);
	exit(1);
      }
    }
    else{
      if (sscanf(buffer,"%d %d %d %lf",&nx,&ny,&nz,&x)!=4){
	fprintf(stderr,"Parse error for data item %d\n",i);
	exit(1);
      }
    }
    if ((nx>ngx)||(ny>ngy)||(nz>ngz)){
      fprintf(stderr,"Found point (%d,%d,%d) but grid is %dx%dx%d\n",
              nx,ny,nz,ngx,ngy,ngz);
      exit(1);
    }
    nx--;ny--;nz--;
    gptr->data[nx*ngy*ngz+ny*ngz+nz]=x*scale;
    if (nspins==2)
      data2[nx*ngy*ngz+ny*ngz+nz]=x2*scale;
  }

  if (i!=ngx*ngy*ngz) error_exit("Too few data points read in denfmt_read");

  if (nspins==2){
    if (inrange(0,spin_range)){
      if (rescale==1){
	if (flags&RAW)
	  gptr->name="Density_raw";
	else
	  gptr->name="Density";
      }
      else if (rescale==2){
	if (flags&RAW)
	  gptr->name="Potential_raw_s0";
	else
	  gptr->name="Potential_s0";
      }
      else if (rescale==1)
	gptr->name="ELF_s0";
      else
	gptr->name="Unknown_s0";
    
      if (inrange(1,spin_range)){
	gptr2=gptr->next;
	gptr2->next=malloc(sizeof(struct grid));
	if (!gptr2->next) error_exit("Malloc error for struct grid");
	gptr2->next->next=NULL;
	gptr2->next->data=NULL;
	gptr2->size[0]=ngx;
	gptr2->size[1]=ngy;
	gptr2->size[2]=ngz;
	gptr2->data=data2;
	if (rescale==1){
	  if (flags&RAW)
	    gptr2->name="Spin_raw";
	  else
	    gptr2->name="Spin";
	}
	else if (rescale==2){
	  if (flags&RAW)
	    gptr2->name="Potential_raw_s1";
	  else
	    gptr2->name="Potential_s1";
	}
	else if (rescale==1)
	  gptr2->name="ELF_s1";
	else
	  gptr2->name="Unknown_s1";
      }
      else free(data2);
    }
    else{ /* if (inrange(0,spin_range) */
      if (inrange(1,spin_range)){
	free(gptr->data);
	gptr->data=data2;
	if (rescale==1){
	  if (flags&RAW)
	    gptr->name="Spin_raw";
	  else
	    gptr->name="Spin";
	}
	else if (rescale==2){
	  if (flags&RAW)
	    gptr->name="Potential_raw_s1";
	  else
	    gptr->name="Potential_s1";
	}
	else if (rescale==1)
	  gptr->name="ELF_s1";
	else
	  gptr->name="Unknown_s1";
      }
      else error_exit("Spin range excludes 0 and 1");
    }
  } /* if spins==2 */
  else{
    if (rescale==1){
      if (flags&RAW)
	gptr->name="Density_raw";
      else
	gptr->name="Density";
    }
    else if (rescale==2){
      if (flags&RAW)
	gptr->name="Potential_raw";
      else
	gptr->name="Potential";
    }
    else if (rescale==1)
      gptr->name="ELF";
    else
      gptr->name="Unknown";
  }

}
