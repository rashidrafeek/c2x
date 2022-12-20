/* A slow, unoptimised GPFA FFT in ANSI C89 for smallish 3D FFTs
 *
 * As written here it has very little benefit, save that it does the job
 *
 * This code is copyright MJ Rutter
 *
 * The author is well aware of faster ways of doing FFTs than this version.
 *
 * Sign convention is 
 *   dir=-1, exponential is -ikx, generally regarded as forward
 *   dir=+1, exponential is +ikx, generally regarded as backwards
 *
 * No scaling in either case.
 */


/* Copyright (c) 2007, 2015 MJ Rutter 
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

#include <math.h>
#include <stdio.h> /* fprintf */
#include <stdlib.h> /* malloc */

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

/* Big should be approximately half the size of the L1 cache, in units
 * of 2*sizeof(double). This code will fail (with an error) if it
 * encounters a prime factor >BIG
 */

#define BIG 2000

#define min(a,b) (a)<(b)?(a):(b)

void m1fft(double* a, int n, int str, int dir, int rep, int rstr);
static void m1spfar(double* a, int p, int n, int dir, int off,int str,
             int max, int r, int rep, int rstr);
static void m1ftrw(double* a, int n, int nn, int dir, int off, int max, int r,
            int rep, int rstr, double* w);
static void m1tr(double* a,int n,int nn1,int nn2,int off,int max,int rep,
            int rstr);
static void factor(int n, int* factors, int* powers);
static int nfactors(int n);

void fft3d(double* a, int* nn, int dir){
  int i,offset,lots,remainder,off;

  off=0;
  remainder=nn[1]*nn[2];
  lots=min(BIG/nn[0],remainder);
  if (lots==0) lots=1;
  while(remainder>0){
    m1fft(a+off,nn[0],1,dir,lots,nn[0]);
    remainder=remainder-lots;
    off=off+2*lots*nn[0];
    lots=min(lots,remainder);
  }

  for(i=0;i<nn[2];i++){
    offset=nn[0]*nn[1]*i;
    m1fft(a+2*offset,nn[1],nn[0],dir,nn[0],1);
  }

  off=0;
  remainder=nn[0]*nn[1];
  lots=min(BIG/nn[2],remainder);
  if (lots==0) lots=1;
  while(remainder>0){
    m1fft(a+off,nn[2],nn[0]*nn[1],dir,lots,1);
    remainder=remainder-lots;
    off=off+2*lots;
    lots=min(lots,remainder);
  }
}

void m1fft(double* a, int n, int str, int dir, int rep, int rstr){
/* Multiple FFTs, GPFA,
 *
 * a -- matrix
 * n -- FFT length
 * str -- stride between elements
 * dir -- direction
 * rep -- number of FFTs to perform
 * rstr -- stride between interleaved FFTs
 */

  int nf,f,f1,f2,i,j;
  int *factors=0,*powers=0;

  nf=nfactors(n);
  factors=malloc(nf*sizeof(int));
  powers=malloc(nf*sizeof(int));
  if ((!factors)||(!powers)){
    fprintf(stderr,"Malloc error in FFT");
    exit(1);
  }
  factor(n,factors,powers);

  for(f=0;f<nf;f++){
    f1=pow(factors[f],powers[f]);
    f2=n/f1;
    for(i=0;i<f2;i++){
      j=(i*f1)%n;
      m1spfar(a,factors[f],powers[f],dir,j*str,f2*str,n*str,f2,rep,rstr);
    }
  }

  free(powers);
  free(factors);

}

static void m1spfar(double* a, int p, int n, int dir, int off,int str,
             int max, int r, int rep, int rstr){
/* Multiple FFTs, single prime factor, using PFA,
 *
 * a -- matrix
 * p -- prime factor
 * n -- its power
 * dir -- direction
 * off -- offset of first element in a
 * max -- calculate indices modulo this
 * r -- rotate by this
 * rep -- number of FFTs to perform
 * rstr -- stride between interleaved FFTs
 */

  int j,k,len,stride,blocks,ind1,ind2;
  int pass,stride_in,stride_butt,butt,nbutt,nsubblock;
  int sb,butt1,butt2;

  double wr[2],wrk[2],twopibyn,tmp;

  len=pow(p,n);
  twopibyn=dir*2*M_PI/len;

/* Simple passes, first half */

  stride=len/p;
  blocks=1;

  for(pass=0;pass<=(n-1)/2;pass++){
    for(j=0;j<blocks;j++){
      ind1=(j*p*stride*str+off)%max;
      wr[0]=cos(((blocks*r)%len)*twopibyn);
      wr[1]=sin(((blocks*r)%len)*twopibyn);
      wrk[0]=1.0;
      wrk[1]=0.0;
      for(k=0;k<stride;k++){
        m1ftrw(a,p,stride*str,dir,ind1,max,r,rep,rstr,wrk);
        tmp=wrk[0];
        wrk[0]=wrk[0]*wr[0]-wrk[1]*wr[1];
        wrk[1]=tmp*wr[1]+wrk[1]*wr[0];
        ind1=ind1+str;
        if (ind1>=max) ind1-=max;
      }
    }
    stride=stride/p;
    blocks=blocks*p;
  }

/* Exchanging passes, second half */

  for(pass=(n-1)/2+1;pass<n;pass++){
    nsubblock=len/pow(p,n-pass-1);
    nbutt=nsubblock/(p*p);
    stride_in=pow(p,n-pass-1);
    stride_butt=pow(p,pass);
    for(sb=0;sb<stride_in;sb++){
      for(butt1=0;butt1<stride_in;butt1++){
        wr[0]=cos(((stride_butt*butt1*r)%len)*twopibyn);
        wr[1]=sin(((stride_butt*butt1*r)%len)*twopibyn);
        for(butt2=0;butt2<nbutt/stride_in;butt2++){
          butt=butt1+p*stride_in*butt2;
          ind1=sb*nsubblock+butt;  
          ind1=(ind1*str+off)%max;
          ind2=ind1;
          for(j=0;j<p;j++){
            m1ftrw(a,p,stride_in*str,dir,ind2,max,r,rep,rstr,wr);
            ind2=ind2+stride_butt*str;
            if (ind2>=max) ind2-=max;
          }
          /* Transpose */
          m1tr(a,p,stride_in*str,stride_butt*str,ind1,max,rep,rstr);
        }
      }
    }
  }
}

static void m1ftrw(double* a, int n, int nn, int dir, int off, int max, int r,
           int rep, int rstr, double* w){
/* Multiple DFTs
 *
 * a -- matrix
 * n -- size
 * nn -- increment
 * dir -- direction
 * off -- offset of first element in a
 * max -- calculate indices modulo this
 * r -- rotate by this
 * rep -- number of DFTs to perform
 * rstr -- stride between interleaved DFTs
 * w -- twiddle factor
 */

  int i,j,k,ind,remainder,lots,off2;
  double b[2*BIG];
  double cr,ci,c0r,c0i,c1r,c1i,tmp;

  /*
  if (n*rep>BIG){fprintf(stderr,"Increase BIG in m1ftr\n"); exit(1);}
  */

  /* That wasn't friendly, so... */

  if (n*rep>BIG){
    remainder=rep;
    lots=BIG/n;
    if (lots==0){
      fprintf(stderr,"Increase BIG in m1ftr\n");
      exit(1);
    }
    off2=off;
    while(remainder>0){
      lots=min(lots,remainder);
      m1ftrw(a,n,nn,dir,off2,max,r,lots,rstr,w);
      remainder-=lots;
      off2+=lots*rstr;
    }
    return;
  }

/* Correct for indexing complexes as doubles */

  off*=2;
  max*=2;
  nn*=2;
  rstr*=2;


  c0r=cos(2*(r%n)*M_PI/n);
  c0i=dir*sin(2*(r%n)*M_PI/n);

  c1r=1.0;
  c1i=0.0;

  for(i=0;i<2*n*rep;i+=2*rep){
    for(j=0;j<2*rep;j++) b[i+j]=0.0;
    cr=1.0;
    ci=0.0;
    ind=off;
    for(k=0;k<rep;k++){
      b[i+2*k]=0.0;
      b[i+2*k+1]=0.0;
    }
    for(j=0;j<n;j++){
      for(k=0;k<rep;k++){
        b[i+2*k]+=a[ind+k*rstr]*cr-a[ind+k*rstr+1]*ci;
        b[i+2*k+1]+=a[ind+k*rstr]*ci+a[ind+k*rstr+1]*cr;
      }
      tmp=cr;
      cr=cr*c1r-ci*c1i;
      ci=tmp*c1i+ci*c1r;
      ind+=nn;
      if (ind>max) ind=ind-max;
    }
    tmp=c1r;
    c1r=c1r*c0r-c1i*c0i;
    c1i=tmp*c0i+c1i*c0r;
  }

  cr=1.0;
  ci=0.0;
  ind=off;
  for(i=0;i<2*n*rep;i+=2*rep){
    for(k=0;k<rep;k++){
      a[ind+k*rstr]=b[i+2*k]*cr-b[i+2*k+1]*ci;
      a[ind+k*rstr+1]=b[i+2*k]*ci+b[i+2*k+1]*cr;
    }
    tmp=cr;
    cr=cr*w[0]-ci*w[1];
    ci=tmp*w[1]+ci*w[0];
    ind+=nn;
    if (ind>max) ind=ind-max;
  }
}



static void m1tr(double* a,int n,int nn1,int nn2,int off,int max,int rep, int rstr){
/* Transpose complex (sub)matrix, stored as real vector
 *
 *  a -- matrix
 *  n -- size
 *  nn1 -- column increment
 *  nn2 -- row increment
 *  off -- offset of first element
 *  max -- calculate indices modulo this
 *  rep -- number of interleaved matrices
 *  rstr -- stride between interleaved matrices
 */

  int i,j,k,ind1,ind2;
  double tmp1,tmp2;

/* Correct for indexing a double array rather than a complex one */
  nn1*=2;
  nn2*=2;
  off*=2;
  max*=2;
  rstr*=2;

  for(i=0;i<=n-2;i++){
    ind1=(i*(nn1+nn2)+off)%max;
    ind2=ind1;
    for(j=i+1;j<=n-1;j++){
      ind1=ind1+nn2;
      if (ind1>=max) ind1-=max;
      ind2=ind2+nn1;
      if (ind2>=max) ind2-=max;
      for(k=0;k<rep*rstr;k+=rstr){
        tmp1=a[ind1+k];
        tmp2=a[ind1+k+1];
        a[ind1+k]=a[ind2+k];
        a[ind1+k+1]=a[ind2+k+1];
        a[ind2+k]=tmp1;
        a[ind2+k+1]=tmp2;
      }
    }
  }
}

static void factor(int n, int* factors, int* powers){
/* Find prime factors, and their powers, in n */
  int i,j;

  j=-1;
  i=2;

  while(n>=i){
    if ((n%i)==0){
      j++;
      factors[j]=i;
      powers[j]=0;
      while ((n%i)==0){
        n/=i;
        powers[j]++;
      }
    }
    i+=2;
    if (i==4) i=3;
  }
}

static int nfactors(int n){
/* Find number of prime factors of n */
  int i,nfact;

  nfact=0;
  i=2;

  while(n>=i){
    if ((n%i)==0){
      nfact++;
      while ((n%i)==0) n/=i;
    }
    i+=2;
    if (i==4) i=3;
  }
  return nfact;
}
