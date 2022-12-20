/* Write a binary Fortran file
 *
 * Read with
 *
 *  open(unit=u,file='foo',form='unformatted')
 *  read(u) ngx,ngy,ngz
 *  allocate (g(ngx,ngy,ngz))
 *  read(u) g
 *
 *  If multiple grids present, the grid size is given just once, so
 *  the next read should be "read(u) g2".
 *
 *  The file does not indicate whether the grid is stored as real or
 *  complex -- you are expected to know!
 *
 */

/* Copyright (c) 2013 MJ Rutter 
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

void reverse4(void *data);
void reverse4n(int *data,int n);
void reverse8n(double *data,int n);

void fbin_write(FILE* outfile, struct grid *g){
  int i,endian,self_little,out_little;
  char *p;
  int buff[5];

  if (sizeof(int)!=4){
    fprintf(stderr,"Sorry: fbin_write() cannot run if sizeof(int)!=4\n");
    exit(1);
  }

  if (!g->data){
    fprintf(stderr,"Cannot write fbin file when no 3D data requested.\n");
    exit(1);
  }

  i=1;
  p=(char*)&i;
  self_little=*p;

  out_little=0;
  if (flags&LITTLE_OUT) out_little=1;

  endian=0; /* Output requested is same as our architecture */
  if (out_little!=self_little) endian=1;

  buff[0]=12;
  buff[1]=g->size[0];
  buff[2]=g->size[1];
  buff[3]=g->size[2];
  buff[4]=12;

  if (endian) reverse4n(buff,5);
  fwrite(buff,4,5,outfile);

  while((g)&&(g->data)){
    i=g->size[0]*g->size[1]*g->size[2];
    buff[0]=8*i;
    if (endian) reverse4(buff);
    fwrite(buff,4,1,outfile);
    if (endian) reverse8n(g->data,i);
    fwrite(g->data,8,i,outfile);
    fwrite(buff,4,1,outfile);
    g=g->next;
  }
}
