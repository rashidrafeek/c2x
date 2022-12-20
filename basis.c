/* Various utility functions for dealing with basis sets */


/* Copyright (c) 2007, 2014 MJ Rutter 
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

/* cart2abc requires this funct from basis.c */
void grid_interp(struct grid *grid,int old_fft[3],
   double old_basis[3][3],double old_recip[3][3]);

/* Check whether set is right- or left-handed */
int is_rhs(double b[3][3]){
  if ((b[0][0]*(b[1][1]*b[2][2]-b[1][2]*b[2][1])
      -b[0][1]*(b[1][0]*b[2][2]-b[1][2]*b[2][0])
      +b[0][2]*(b[1][0]*b[2][1]-b[1][1]*b[2][0]))<0) return(0);
  else return(1);
}

/* Update reciprocal basis set to reflect real basis set */
/* return a negative volume if fed a lhs of axes */
void real2rec(struct unit_cell *c){
  int i,j,k;
  double rec[3][3],v;
  double (*b)[3];

  b=c->basis;

  v=0.0;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      v+=b[i][j]*b[i][j];

  if (v==0.0) error_exit("all cell axes are length zero: no basis found?");

  for(i=0;i<3;i++){
    j=(i+1)%3;
    k=(i+2)%3;
    rec[i][0]=b[j][1]*b[k][2]-b[j][2]*b[k][1];
    rec[i][1]=-b[j][0]*b[k][2]+b[j][2]*b[k][0];
    rec[i][2]=b[j][0]*b[k][1]-b[j][1]*b[k][0];
  }

  v=0.0;
  for(i=0;i<3;i++) v+=b[0][i]*rec[0][i];

  if (v==0.0) error_exit("unit cell volume is zero in real2rec.");

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->recip[i][j]=rec[i][j]/v;

  c->vol=v;
}

/* Update global fractional atomic co-ordinates to reflect global
 * absolute atomic co-ordinates. Assumes global reciprocal basis
 * is correct. Does NOT reduce co-ords to be less than 1 --
 * sym_vec relies on this.
 */
void addfrac(struct atom *a,int na, double rec[3][3]){
  int i,j,k;

  for(i=0;i<na;i++){
    for(j=0;j<3;j++){
      a[i].frac[j]=0;
      for(k=0;k<3;k++)
         a[i].frac[j]+=a[i].abs[k]*rec[j][k];
      if ((a[i].frac[j]<0.0)&&(a[i].frac[j]>-1e-15)) a[i].frac[j]=0.0;
    }
  }
}

/* Update global absolute atomic co-ordinates to reflect global
 * fractional atomic co-ordinates
 */
void addabs(struct atom *a,int na, double b[3][3]){
  int i,j,k;

  for(i=0;i<na;i++){
    for(j=0;j<3;j++){
      a[i].abs[j]=0;
      for(k=0;k<3;k++)
         a[i].abs[j]+=a[i].frac[k]*b[k][j];
    }
  }
}

/* Fix co-ordinates to force everything into 1st unit cell */

void reduce_cell(struct atom *a, int na, double b[3][3]){
  int i,j,fixed=0;

  for(i=0;i<na;i++){
    for(j=0;j<3;j++){
      if ((a[i].frac[j]>=0)&&(a[i].frac[j]<1)) continue;
      fixed=1;
      a[i].frac[j]=fmod(a[i].frac[j],1.0);
      if (a[i].frac[j]<0) a[i].frac[j]+=1.0;
    }
  }

  if (fixed) addabs(a,na,b);
}

/* Fix co-ordinates to force everything into 1st unit cell (with tolerance) */

void reduce_cell_tol(struct atom *a, int na, double b[3][3]){
  int i,j,fixed=0;

  for(i=0;i<na;i++){
    for(j=0;j<3;j++){
      if ((a[i].frac[j]>=0)&&(a[i].frac[j]<1-tol)) continue;
      fixed=1;
      a[i].frac[j]=fmod(a[i].frac[j],1.0);
      if (a[i].frac[j]<0) a[i].frac[j]+=1.0;
      if (a[i].frac[j]>1-tol) a[i].frac[j]=0;
    }
  }

  if (fixed) addabs(a,na,b);
}

void abc2cart(double *abc, struct unit_cell *c){
/* convert from a,b,c,alpha,beta,gamma to Cartesian basis */
  double alpha,beta,gamma,x;
  double (*b)[3];
  int i;

  b=c->basis;

  if (debug>2)fprintf(stderr,"abc2cart: original basis:\n%f %f %f\n%f %f %f\n",
     abc[0],abc[1],abc[2],abc[3],abc[4],abc[5]);

  alpha=abc[3]*M_PI/180;
  beta=abc[4]*M_PI/180;
  gamma=abc[5]*M_PI/180;

/* a lies along x axis */

  b[0][0]=abc[0];
  b[0][1]=0.0;
  b[0][2]=0.0;

/* b is in xy plane and angle gamma to x */

  b[1][0]=abc[1]*cos(gamma);
  b[1][1]=abc[1]*sin(gamma);
  b[1][2]=0.0;

/* a,b,c is a right-hand set */

  b[2][0]=abc[2]*cos(beta);

/* basis[2].basis[1] = abc[2] abc[1] cos(alpha) */

  x=abc[2]*abc[1]*cos(alpha)-b[1][0]*b[2][0];
  b[2][1]=0;
  if (fabs(x)>1e-20){
    if (fabs(b[1][1])>1e-30){
      b[2][1]=x/b[1][1];
    }else{
      error_exit("impossible problem in abc2cart, perhaps gamma is zero.\n");
    }
  }

/* And mod(basis[2][])=abc[2] */

  b[2][2]=sqrt(abc[2]*abc[2]-b[2][0]*b[2][0]
                                -b[2][1]*b[2][1]);

/* Worry about handedness */

  real2rec(c);
  if (c->vol<0){
    for(i=0;i<3;i++) b[2][i]*=-1;
    c->vol*=-1;
  }

  if(debug>2){
    int i;
    fprintf(stderr,"abc2cart: final basis:\n");
      for(i=0;i<=2;i++)
        fprintf(stderr,"%f %f %f\n",b[i][0],b[i][1],b[i][2]);
  }

}

/* Convert from cartesian basis set to a,b,c,alpha,beta,gamma */
/* If m and gptr are NULL, and fix==0, just calculate a,b,c,alpha,beta,gamma */
/* If m not NULL, rotate axes to a along x, b in x-y plane, etc */

void cart2abc(struct unit_cell *c, struct contents *m, double *abc, 
              struct grid *gptr, int fix){
  int i,j,k,itmp,old_fft[3],preserve_c;
  double dtmp,b[3][3];
  double old_basis[3][3],old_recip[3][3];

  preserve_c=(flags&PRESERVE_C)>>PC_SHIFT;

  if (debug>2) fprintf(stderr,"cart2abc called with motif %s\n",
                       m?"present":"absent");

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      b[i][j]=c->basis[i][j];

  for(i=0;i<3;i++)
    abc[i]=sqrt(b[i][0]*b[i][0]+b[i][1]*b[i][1]+
                b[i][2]*b[i][2]);

  for(i=3;i<6;i++){
    j=(i+1)%3;
    k=(i+2)%3;
    abc[i]=acos((b[j][0]*b[k][0]+b[j][1]*b[k][1]+
                 b[j][2]*b[k][2])/(abc[j]*abc[k]))*180/M_PI;
  }

  /* We may now have a different orientation to the original, so,
     if we were passed a motif: */

  if (m){
    if (is_rhs(b)==0){ /* We are trying to convert a lhs of vectors
                            * to abc format.
               * We should warn people,
               * Exchange second and third axes
               * Exchange second and third fractional co-ords
               */
      if (flags&LHS_FUDGE) {
         fprintf(stderr,"Need rh set of vectors for abc output.\n"
                        "Have lh set, and exchange prohibitted...\n");
      }
      else {
        if (debug) fprintf(stderr,"Need rh set of vectors for abc output. "
                           "Exchanging basis vectors %d and %d.\n",
                           (1+preserve_c)%3+1,(2+preserve_c)%3+1);

        real2rec(c);
        for(i=0;i<3;i++)
          for(j=0;j<3;j++)
            old_basis[i][j]=b[i][j];
        for(i=0;i<3;i++)
          for(j=0;j<3;j++)
            old_recip[i][j]=c->recip[i][j];

        for(i=0;i<3;i++){
          dtmp=b[(1+preserve_c)%3][i];
          b[(1+preserve_c)%3][i]=b[(2+preserve_c)%3][i];
          b[(2+preserve_c)%3][i]=dtmp;
        }
        real2rec(c);
        for(i=0;i<m->n;i++){
          dtmp=m->atoms[i].frac[(1+preserve_c)%3];
          m->atoms[i].frac[(1+preserve_c)%3]=m->atoms[i].frac[(2+preserve_c)%3];
          m->atoms[i].frac[(2+preserve_c)%3]=dtmp;
        }

        dtmp=abc[(1+preserve_c)%3];
        abc[(1+preserve_c)%3]=abc[(2+preserve_c)%3];
        abc[(2+preserve_c)%3]=dtmp;

        dtmp=abc[(1+preserve_c)%3+3];
        abc[(1+preserve_c)%3+3]=abc[(2+preserve_c)%3+3];
        abc[(2+preserve_c)%3+3]=dtmp;

        while((gptr)&&(gptr->data)){
          if (debug>1) fprintf(stderr,"Exchanging axes for 3D grid %dx%dx%d\n",
			     gptr->size[0],gptr->size[1],gptr->size[2]);
          for(i=0;i<3;i++) old_fft[i]=gptr->size[i];
          itmp=gptr->size[(1+preserve_c)%3];
          gptr->size[(1+preserve_c)%3]=gptr->size[(2+preserve_c)%3];
          gptr->size[(2+preserve_c)%3]=itmp;
          grid_interp(gptr,old_fft,old_basis,old_recip);
          gptr=gptr->next;
        }
      }
    } /* if (m) */

    abc2cart(abc,c);
    real2rec(c);
    addabs(m->atoms,m->n,b);
  }

}


/* Convert from cartesian basis set to a,b,c,alpha,beta,gamma */
void basis2abc(double b[3][3], double abc[6]){
  int i,j,k;

  for(i=0;i<3;i++)
    abc[i]=sqrt(b[i][0]*b[i][0]+b[i][1]*b[i][1]+
                b[i][2]*b[i][2]);

  for(i=3;i<6;i++){
    j=(i+1)%3;
    k=(i+2)%3;
    abc[i]=acos((b[j][0]*b[k][0]+b[j][1]*b[k][1]+
                 b[j][2]*b[k][2])/(abc[j]*abc[k]))*180/M_PI;
  }
}

int atom_in_list(struct atom *b, struct atom *a, int n){
  int hit,i;

  hit=-1;
  for(i=0;i<n;i++){
    if ((a[i].atno==b->atno)&&
        (fabs(a[i].frac[0]-b->frac[0])<tol)&&
        (fabs(a[i].frac[1]-b->frac[1])<tol)&&
        (fabs(a[i].frac[2]-b->frac[2])<tol)){
      hit=i;
      break;
    }
  }

  return hit;
}

