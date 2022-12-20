/* Read a CASTEP .chdiff file
 *
 * Cope with either endianness
 *
 * Much code in common with esp_read.c -- these should really be combined
 */

/* Copyright (c) 2008 MJ Rutter 
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

static void reverse4(void *data){
/* reverse endian a single 4 byte int */
   int out;
   char *p1,*p2;

   p1=(char*)data;
   p2=((char*)&out)+3;

   *(p2--)=*(p1++);
   *(p2--)=*(p1++);
   *(p2--)=*(p1++);
   *p2=*p1;

   *((int*)data)=out;
}

static void reverse8n(double *data, int n){
/* reverse endian n words of 8 byte data */
   int i;
   double out;
   char *p1,*p2;

   for(i=0;i<n;i++){
     p1=(char*)(data+i);
     p2=(char*)&out;
  
     p2=p2+7;
  
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
   
     *(data+i)=out;
   }
}

void chdiff_read(FILE* infile, struct grid *gptr){
  int endian,tmp,i,j,k;
  int fft[3];
  double *column,*dptr1,*dptr2;

  if (debug>2) fprintf(stderr,"chdiff_read called\n");

  endian=0; /* Keep compiler quiet */
  fread(&tmp,4,1,infile);

  /* The first record is of length 12. Being Fortran, the first
   * item will therefore be an integer, 4 bytes, of value 12. If we
   * have an endian problem, we will see this as an integer of value
   * 12*(1<<24)
   */

  if (tmp==12) endian=0;
  else if (tmp==12*(1<<24)) endian=1;
  else error_exit("Not a .chdiff file");

  fread(fft,12,1,infile);
  if (endian) {reverse4(fft);reverse4(fft+1);reverse4(fft+2);}
  if (debug>1) fprintf(stderr,"Grid of %dx%dx%d in .chdiff\n",fft[0],
		       fft[1],fft[2]);
  fseek(infile,4,SEEK_CUR);

  /* Set up grid to receive data */

  gptr->next=malloc(sizeof(struct grid));
  if (!gptr->next) error_exit("Malloc error for struct grid");
  gptr->next->next=NULL;
  gptr->next->data=NULL;
  gptr->name="Density_diff";
  for(i=0;i<3;i++) gptr->size[i]=fft[i];
  gptr->data=malloc(8*fft[0]*fft[1]*fft[2]);
  if (!gptr->data) error_exit("Malloc error for chdiff data grid");

  /* The data will be presented in columns of complexes */

  column=malloc(16*fft[2]);
  if (!column) error_exit("Malloc error for chdiff column");

  fread(&tmp,4,1,infile);
  if (endian) reverse4(&tmp);
  if (tmp!=16*fft[2]+8) error_exit("Error reading .chdiff file");
  while(tmp==16*fft[2]+8){
    fread(&i,4,1,infile);
    if (endian) reverse4(&i);
    fread(&j,4,1,infile);
    if (endian) reverse4(&j);
    fread(column,16*fft[2],1,infile);
    if (endian) reverse8n(column,2*fft[2]);
    /* Copy the complex column into our real grid, at specified location */
    dptr1=gptr->data+((i-1)*fft[1]+(j-1))*fft[2];
    dptr2=column;
    for(k=0;k<fft[2];k++){*dptr1++=*dptr2; dptr2+=2;}
    fseek(infile,4,SEEK_CUR);
    if (fread(&tmp,4,1,infile)==0) break;
    if (endian) reverse4(&tmp) ;
    if (tmp!=16*fft[2]+8) error_exit("Error reading .chdiff file");
  }
  free(column);
}

