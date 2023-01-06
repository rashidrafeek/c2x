/* Write a CCP4 / .map file. See
 *  http://www.ccpem.ac.uk/mrc_format/mrc2014.php and
 *  http://www.sciencedirect.com/science/article/pii/S104784771500074X
 */

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
#include<stdlib.h>
#include<math.h>

#include "c2xsf.h"

void ccp4_write(FILE* outfile,struct unit_cell *c, struct grid *g){
  int i,j,k,npts;
  double *dptr1,abc[6],mx,mn,sum,sumsq;
  float fabc[6],ftmp,*fptr,*fptr1;
  struct contents m;
  
  if (sizeof(int)!=4) error_exit("sizeof(int)!=4 in ccp4_write");
  if (sizeof(float)!=4) error_exit("sizeof(float)!=4 in ccp4_write");

  if ((!g)||(!g->data)) error_exit("No data to write");

  npts=g->size[0]*g->size[1]*g->size[2];
  m.n=0;
  
  dptr1=g->data;
  sum=mx=mn=*dptr1;
  sumsq=*dptr1*(*dptr1);
  for(i=1;i<npts;i++){
    mx=max(mx,dptr1[i]);
    mn=min(mn,dptr1[i]);
    sum+=dptr1[i];
    sumsq+=dptr1[i]*dptr1[i];
  }

  /* This may change g->size if a lhs to rhs conversion is req. */

  cart2abc(c,&m,abc,g);
  
  /* NX, NY, NZ */
  
  fwrite(g->size,4,3,outfile);

  /* MODE */
  
  i=2;
  fwrite(&i,4,1,outfile);

  /* NXSTART, NYSTART, NZSTART */
  
  i=0;
  for(j=0;j<3;j++)
    fwrite(&i,4,1,outfile);

  /* MX, MY, MZ */

  fwrite(g->size,4,3,outfile);

  /* CELL */

  for(i=0;i<6;i++)
    fabc[i]=abc[i];
  fwrite(fabc,4,6,outfile);

  /* MAPC, MAPR, MAPS */
  
  i=1;
  fwrite(&i,4,1,outfile);
  i=2;
  fwrite(&i,4,1,outfile);
  i=3;
  fwrite(&i,4,1,outfile);

  /* DMIN, DMAX, DMEAN */
  
  ftmp=mn;
  fwrite(&ftmp,4,1,outfile);
  ftmp=mx;
  fwrite(&ftmp,4,1,outfile);
  ftmp=sum/npts;
  fwrite(&ftmp,4,1,outfile);

  if (debug) fprintf(stderr,"Min: %lf  max:  %lf   mean:  %lf\n",mn,
		     mx,sum/npts);

  
  /* ISPG */
  
  i=1;
  fwrite(&i,4,1,outfile);

  /* NSYBT */
  
  i=0;
  fwrite(&i,4,1,outfile);

  /* EXTRA */
  
  for(j=25;j<=49;j++)
    fwrite(&i,4,1,outfile);

  /* ORIGIN */
  
  ftmp=0;
  for(j=0;j<3;j++)
    fwrite(&ftmp,4,1,outfile);

  /* MAP */
  
  fwrite("MAP",1,4,outfile);

  /* MACHST */
  
  fputc(0x44,outfile);
  fputc(0x44,outfile);
  fputc(0,outfile);
  fputc(0,outfile);

  /* RMS */
  
  ftmp=sqrt((sumsq-sum*sum)/npts);
  fwrite(&ftmp,4,1,outfile);

  if (debug) fprintf(stderr,"RMS deviation from mean: %f\n",ftmp);
  
  /* NLABL */
  
  i=0;
  fwrite(&i,4,1,outfile);

  /* 800 bytes of spaces */
  
  for(i=225;i<=1024;i++)
    fputc(' ',outfile);

  fptr=malloc(npts*4);
  if (!fptr) error_exit("Malloc error for data in ccp4_write");

  dptr1=g->data;
  fptr1=fptr;

  /* Transpose and change precision */
  for(i=0;i<g->size[2];i++)
    for(j=0;j<g->size[1];j++)
      for(k=0;k<g->size[0];k++)
	*(fptr1++)=dptr1[i+j*g->size[2]+k*g->size[2]*g->size[1]];

  fwrite(fptr,4,npts,outfile);
      

  free(fptr);
  
}
