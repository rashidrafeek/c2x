/* Function is passed two vectors. It applies the rotation
 * which moves first to second to the basis set, returning
 * that in new_cell. For obfuscation, the input vectors are
 * also passed in new_cell.
 */


/* Copyright (c) 2007 MJ Rutter 
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

#include <stdio.h>
#include <math.h>

#include "c2xsf.h"

int is_rhs(double b[3][3]); /* Found in basis.c */

void rotation(struct unit_cell *c, struct contents *m, double new_cell[3][3]){
  double v1[3],v2[3],e[3],angle,tmp;
  double rot_mat[3][3];
  int i,j,k;
  
  /* If new_cell[2] is not the zero vector, use it as Euler vector
   * else determine our own Euler vector
   */

  if ((new_cell[2][0]*new_cell[2][0]+new_cell[2][1]*new_cell[2][1]+
       new_cell[2][2]*new_cell[2][2])==0.0){


    /* Extract two vectors in absolute co-ords*/
    for(i=0;i<3;i++){
      v1[i]=new_cell[0][i];
      v2[i]=new_cell[1][i];
    }

    /* euler is v1xv2 */

    e[0]= v1[1]*v2[2]-v1[2]*v2[1];
    e[1]=-v1[0]*v2[2]+v1[2]*v2[0];
    e[2]= v1[0]*v2[1]-v1[1]*v2[0];

    /* Which needs normalising */

    tmp=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);

    if (tmp<1e-20) error_exit("Impossibly small cross product in rotation");

    e[0]/=tmp;
    e[1]/=tmp;
    e[2]/=tmp;

  }else{
    e[0]=new_cell[2][0];
    e[1]=new_cell[2][1];
    e[2]=new_cell[2][2];

    tmp=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);

    if (tmp<1e-20) error_exit("Impossibly small euler axis in rotation");

    e[0]/=tmp;
    e[1]/=tmp;
    e[2]/=tmp;

    /* project the vectors given onto plane perpendicular to euler axis */

    tmp=e[0]*new_cell[0][0]+e[1]*new_cell[0][1]+e[2]*new_cell[0][2];
    for(i=0;i<3;i++) v1[i]=new_cell[0][i]-tmp*e[i];

    tmp=e[0]*new_cell[1][0]+e[1]*new_cell[1][1]+e[2]*new_cell[1][2];
    for(i=0;i<3;i++) v2[i]=new_cell[1][i]-tmp*e[i];
  }

    /* angle is arccos(v1.v2 / mod(v1).mod(v2)) */
  angle=acos((v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])/
              sqrt((v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]) *
                   (v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2])));

  /* We have a sign ambiguity here, which we should try to resolve */

  for(i=0;i<3;i++) rot_mat[0][i]=v1[i];
  for(i=0;i<3;i++) rot_mat[1][i]=v2[i];
  for(i=0;i<3;i++) rot_mat[2][i]=e[i];

  if (debug>2) fprintf(stderr,"Sign is %s\n",is_rhs(rot_mat)?"+":"-");

  if (!is_rhs(rot_mat)) angle*=-1;

  if (debug>1) fprintf(stderr,"Rotating by %g degrees about (%g,%g,%g)\n",
                        angle*180/M_PI,e[0],e[1],e[2]);

  /* Wikipedia says that the rotation matrix is:
   *
   *  I cos(theta) + (1-cos(theta))ee^t - E sin(theta)
   *
   *                 ( 0  -e3   e2 )
   * where       E = ( e3  0   -e1 )
   *                 ( -e2 e1   0  )
   *
   * with e=(e1,e2,e3) being the Euler vector (axis of rotation)
   */

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      rot_mat[i][j]=0;

  for(i=0;i<3;i++) rot_mat[i][i]=cos(angle);

  tmp=1-cos(angle);

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      rot_mat[i][j]+=tmp*e[i]*e[j];

  tmp=sin(angle);

  rot_mat[0][1]-=tmp*e[2];
  rot_mat[0][2]+=tmp*e[1];
  rot_mat[1][0]+=tmp*e[2];
  rot_mat[1][2]-=tmp*e[0];
  rot_mat[2][0]-=tmp*e[1];
  rot_mat[2][1]+=tmp*e[0];

  if (debug>2){
    fprintf(stderr,"Rotation matrix\n");
    for(i=0;i<=2;i++)
      fprintf(stderr,"%f %f %f\n",
              rot_mat[i][0],rot_mat[i][1],rot_mat[i][2]);
    if (is_rhs(rot_mat)) fprintf(stderr,"(Determinant is >=0)\n");
    else fprintf(stderr,"(Warning: determinant is <0\n");
    fprintf(stderr,"(Old basis is a %shs)\n",is_rhs(c->basis)?"r":"l");
  }


  /* Now apply rotation to basis */

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      new_cell[i][j]=0;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      for(k=0;k<3;k++)
        new_cell[i][j]+=rot_mat[j][k]*c->basis[i][k];

  if (debug>1){
    fprintf(stderr,"New basis set\n");
    for(i=0;i<=2;i++)
      fprintf(stderr,"%f %f %f\n",
              new_cell[i][0],new_cell[i][1],new_cell[i][2]);
  }


  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[i][j]=new_cell[i][j];

  if (debug>2)
    fprintf(stderr,"(New basis is a %shs)\n",is_rhs(c->basis)?"r":"l");

  /* Everything remains as was in relative co-ordinates, but the
   * absolute co-ords of the atoms need recalculating, as does the
   * reciprocal basis set
   */

  addabs(m->atoms,m->n,c->basis);
  real2rec(c);

}
