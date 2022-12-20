/* Read a CASTEP .esp file
 *
 * Cope with either endianness
 *
 * Note that Castep always uses units of -Hartrees for this file
 * (at least up to Castep 18.1)
 */

/* Copyright (c) 2008-2018 MJ Rutter 
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
#include<math.h>

#include "c2xsf.h"

int inrange(int x, char *range); /* from check_read.c */

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

void esp_read(FILE* infile, struct contents *m, struct grid *gptr,
	      struct es *elect){
  int endian,tmp,i,j,k,nspins,ns,c;
  int fft[3];
  double *column,*dptr1,*dptr2,csum,scale;

  endian=0; /* Keep compiler quiet */
  csum=0;

  if (debug>2) fprintf(stderr,"esp_read called\n");

  fread(&tmp,4,1,infile);

  /* The first record is of length 4. Being Fortran, the first
   * item will therefore be an integer, 4 bytes, of value 12. If we
   * have an endian problem, we will see this as an integer of value
   * 4*(1<<24)
   */

  if (tmp==4) endian=0;
  else if (tmp==4*(1<<24)) endian=1;
  else error_exit("Not a .cst_esp file");

  fread(&nspins,4,1,infile);

  if(endian) reverse4(&nspins);

  if ((nspins!=1)&&(nspins!=2)){
    fprintf(stderr,"Unexpected number of spins %d. Aborting\n",nspins);
    exit(1);
  }

  fseek(infile,4,SEEK_CUR);

  fread(&tmp,4,1,infile);
  if (endian) reverse4(&tmp);
  if (tmp!=12) error_exit("Error reading FFT grid size in cst_esp file");

  fread(fft,12,1,infile);
  if (endian) {reverse4(fft);reverse4(fft+1);reverse4(fft+2);}
  if (debug>1) fprintf(stderr,"Grid of %dx%dx%d in .cst_esp\n",fft[0],
                       fft[1],fft[2]);
  fseek(infile,4,SEEK_CUR);


  /* The data will be presented in columns of complexes */

  column=malloc(16*fft[2]);
  if (!column) error_exit("Malloc error for cst_esp column");

  for(ns=0;ns<nspins;ns++){
    if (inrange(ns,elect->spin_range)){
  /* Set up grid to receive data */

      if (gptr->next) gptr=gptr->next;
      gptr->next=malloc(sizeof(struct grid));
      if (!gptr->next) error_exit("Malloc error for struct grid");
      gptr->next->next=NULL;
      gptr->next->data=NULL;
      if (nspins==1) {
	gptr->name=dict_get(m->dict,"grid_name");
	if (!gptr->name) gptr->name="ESP";
      }
      else{
        gptr->name=malloc(40);
        if (!gptr->name) error_exit("Malloc error for grid name");
	if (dict_get(m->dict,"grid_name"))
	  sprintf(gptr->name,"%s_s%d",(char*)dict_get(m->dict,"grid_name"),ns);
	else
	  sprintf(gptr->name,"ESP_s%d",ns);
      }
      for(i=0;i<3;i++) gptr->size[i]=fft[i];
      gptr->data=malloc(8*fft[0]*fft[1]*fft[2]);
      if (!gptr->data) error_exit("Malloc error for pot data grid");

 
      for(c=0;c<fft[0]*fft[1];c++){
        fread(&tmp,4,1,infile);
        if (endian) reverse4(&tmp);
        if (tmp!=16*fft[2]+8) error_exit("Error reading .cst_esp file");
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
        if (debug){
          dptr2=column+1;
          for(k=0;k<fft[2];k++){csum+=fabs(*dptr2); dptr2+=2;}
        }
        fseek(infile,4,SEEK_CUR);
      }
      if (debug){
        if(nspins==1)
          fprintf(stderr,
                  "Sum of absolute values of imaginary parts of data"
                  " is %g\n",
                  csum);
        else
          fprintf(stderr,"For spin %d, "
                  "sum of absolute values of imaginary parts of data"
                  " is %g\n",
                  ns,csum);
      }
      /* Castep stores the potential in -Hartrees */
      /* From 2.32 we show all potentials in -V */
      if (!(flags&RAW)){
        scale=H_eV;
        for(i=0;i<fft[0]*fft[1]*fft[2];i++) gptr->data[i]*=scale;
        if (debug>1) fprintf(stderr,"Potential rescaled by %f\n",scale);
      }
    }
    else fseek(infile,fft[0]*fft[1]*(16*fft[2]+16),SEEK_CUR);
  }
  free(column);
}
