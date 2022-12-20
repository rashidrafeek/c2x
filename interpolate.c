/* FFT interpolation of 3D real grid data */

/* Copyright (c) 2014 MJ Rutter 
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

void fft3d(double *c, int *ngptar, int dir);
void pad_recip(double *o, int fft[3], double **ptr, int nfft[3]);

/* Do FFT interpolation of real space real data */
void interpolate3d(struct grid *old_grid, struct grid *new_grid){
  int old_size,new_size,ffft[3],fft[3],nfft[3];
  int i;
  double *o,*n,scale;

  /* A new grid dimension of zero means leave as was */
  for(i=0;i<3;i++)
    if (new_grid->size[i]==0) new_grid->size[i]=old_grid->size[i];

  for(i=0;i<3;i++) fft[i]=old_grid->size[i];
  for(i=0;i<3;i++) nfft[i]=new_grid->size[i];

  if (debug>1)
    fprintf(stderr,"Interpolating real data from %dx%dx%d to %dx%dx%d\n",
            fft[0],fft[1],fft[2],nfft[0],nfft[1],nfft[2]);

  old_size=fft[0]*fft[1]*fft[2];
  new_size=nfft[0]*nfft[1]*nfft[2];
  new_grid->data=malloc(new_size*sizeof(double));
  if (!new_grid->data) error_exit("Malloc error for final grid in interpolate");

  if((fft[0]==nfft[0])&&(fft[1]==nfft[1])&&(fft[2]==nfft[2])){
    if (debug>1) fprintf(stderr,"Null interpolation reduced to copy.\n");
    for(i=0;i<old_size;i++) new_grid->data[i]=old_grid->data[i];
    return;
  }

  /* Pad real data to complex */

  o=malloc(2*old_size*sizeof(double));
  if (!o) error_exit("Malloc error for first grid in interpolate");

  for(i=0;i<old_size;i++){
    o[2*i]=old_grid->data[i];
    o[2*i+1]=0.0;
  }

  /* FFT to reciprocal space */

  /* A FORTRAN data order ... */

  ffft[0]=fft[2];
  ffft[1]=fft[1];
  ffft[2]=fft[0];

  if (debug>1) fprintf(stderr,"first FFT in interpolate\n");

  fft3d(o,ffft,-1);

  /* Pad onto interpolated reciprocal space grid */

  /* Assume all bits zero is a double zero */

  if (debug>1) fprintf(stderr,"padding in interpolate\n");

  pad_recip(o,fft,&n,nfft);

  free(o);

  /* FFT back to real space */

  if (debug>1) fprintf(stderr,"second FFT in interpolate\n");

  ffft[0]=nfft[2];
  ffft[1]=nfft[1];
  ffft[2]=nfft[0];

  fft3d(n,ffft,1);

  if (debug>1) fprintf(stderr,"end of second FFT in interpolate\n");

  /* Convert back to real and rescale */

  scale=1.0/old_size;

  for(i=0;i<new_size;i++)
    new_grid->data[i]=scale*n[2*i];

  free(n);

}


double interpolate0d(struct grid *gptr,double x_in[3]){
  int i,j,n;
  int ii,jj,kk;
  int ny,nz;
  int v[2][2][2][3];
  double c1[2][2],c2[2];
  double t,z,x[3];

  ny=gptr->size[1];
  nz=gptr->size[2];

  for(i=0;i<3;i++){
    x_in[i]=fmod(x_in[i],1.0);
    if (x_in[i]<0) x_in[i]++;
  }
  
  for(i=0;i<3;i++) x[i]=x_in[i]*gptr->size[i];

  for(i=0;i<3;i++) {
    v[0][0][0][i]=(int)floor(x[i])%gptr->size[i];
    if (v[0][0][0][i]<0) v[0][0][0][i]+=gptr->size[i];
    x[i]-=floor(x[i]);
  }

  /* Fill in x co-ords of all cube vertices */

  for(i=0;i<2;i++)
    for(j=0;j<2;j++){
      v[0][i][j][0]=v[0][0][0][0];
      n=(v[0][0][0][0]+1)%gptr->size[0];
      v[1][i][j][0]=n;
    }

  /* y co-ords */

  for(i=0;i<2;i++)
    for(j=0;j<2;j++){
      v[i][0][j][1]=v[0][0][0][1];
      n=(v[0][0][0][1]+1)%gptr->size[1];
      v[i][1][j][1]=n;
    }

  /* z co-ords */

  for(i=0;i<2;i++)
    for(j=0;j<2;j++){
      v[i][j][0][2]=v[0][0][0][2];
      n=(v[0][0][0][2]+1)%gptr->size[2];
      v[i][j][1][2]=n;
    }

  /* First step of interpolation: push out z */

  for(i=0;i<2;i++)
    for(j=0;j<2;j++){
      ii=v[i][j][0][0];
      jj=v[i][j][0][1];
      kk=v[i][j][0][2];
      t=(1-x[2])*gptr->data[kk+nz*(jj+ny*ii)];
      ii=v[i][j][1][0];
      jj=v[i][j][1][1];
      kk=v[i][j][1][2];
      t+=x[2]*gptr->data[kk+nz*(jj+ny*ii)];
      c1[i][j]=t;
    }

  /* Second step, push out y */

  for(i=0;i<2;i++){
    t=(1-x[1])*c1[i][0];
    t+=x[1]*c1[i][1];
    c2[i]=t;
  }

  z=(1-x[0])*c2[0]+x[0]*c2[1];

  return(z);
      
}

void interpolate1d(struct grid *gptr, double st[3], double end[3],
                   int npts, double *out){
  int i,j;
  double x[3];

  for(i=0;i<npts;i++){
    for(j=0;j<3;j++)
      x[j]=st[j]+(i/(double)(npts-1))*(end[j]-st[j]);
    out[i]=interpolate0d(gptr,x);
  }
}

/* Pad into array assumed to be zeroed */
/* Deal with case of target being both larger and smaller than source */
void pad_recip(double *o, int fft[3], double **nptr, int nfft[3]){
  int i,j,k,ii,jj,kk,nii,njj,nkk;
  int imin,imax,jmin,jmax,kmin,kmax;
  int ind,nind;
  double *n;

  for(i=0;i<3;i++)
    if (nfft[i]==0) nfft[i]=fft[i];

  for(i=0;i<3;i++)
    if (nfft[i]<1) {
      fprintf(stderr,"Invalid grid size %dx%dx%d\n",nfft[0],nfft[1],nfft[2]);
      exit(1);
    }

  if (debug) fprintf(stderr,"Moving from %dx%dx%d to %dx%dx%d recip grid\n",
                     fft[0],fft[1],fft[2],nfft[0],nfft[1],nfft[2]);

  *nptr=calloc(nfft[0]*nfft[1]*nfft[2],2*sizeof(double));
  if (!*nptr) error_exit("Malloc error for interpolated grid\n");
  n=*nptr;

  imax=min(fft[0]/2,nfft[0]/2);
  imin=max((1-fft[0])/2,(1-nfft[0])/2);

  jmax=min( fft[1]/2,nfft[1]/2);
  jmin=max((1-fft[1])/2,(1-nfft[1])/2);

  kmax=min(fft[2]/2,nfft[2]/2);
  kmin=max((1-fft[2])/2,(1-nfft[2])/2);
  
  for(i=imin;i<=imax;i++){
    if (i>=0) {
      ii=i;
      nii=i;
    }
    else{
      ii=fft[0]+i;
      nii=nfft[0]+i;
    }
    for(j=jmin;j<=jmax;j++){
      if (j>=0) {
        jj=j;
        njj=j;
      }
      else{
        jj=fft[1]+j;
        njj=nfft[1]+j;
      }
      for(k=kmin;k<=kmax;k++){
        if (k>=0) {
          kk=k;
          nkk=k;
        }
        else{
          kk=fft[2]+k;
          nkk=nfft[2]+k;
        }
        ind=2*(kk+fft[2]*(jj+ii*fft[1]));
        nind=2*(nkk+nfft[2]*(njj+nii*nfft[1]));
        n[nind]=o[ind];
        n[nind+1]=o[ind+1];
      }
    }
  }
}

