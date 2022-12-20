/* Copyright (c) 2013 -- 2017 MJ Rutter 
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
#include<string.h>

#include "c2xsf.h"


/* routine to identify a symmetry operation */

void vcross(double a[3],double b[3],double c[3]);
int minvert(double m[3][3]);
double vmod2(double v[3]);

void frac_print(FILE *out, double x){
  int i;

  if (fabs(x)<tol){
    fprintf(out,"0");
    return;
  }

  if (x<0){
    x*=-1;
    fprintf(out,"-");
  }
  else fprintf(out,"+");

  if (aeq(x,1)){
    fprintf(out,"1");
    return;
  }
  
  for(i=2;i<=13;i++){
    if (aeq(floor(x*i+tol),x*i)) break;
  }
  if (i==13) fprintf(out,"%f",x);
  else fprintf(out,"%d/%d",(int)(x*i+tol),i);
}

/* Print symmetry matrix in the "equivalent" form of
 * -x,-y,z+1/2 etc. Must not ever print a space.
 */

void equiv_sym(struct sym_op *s, struct unit_cell *c, FILE *out){
  int i,j;
  double m[3][3],v[3];


  /* Do the trivial bit of re-expressing in frac co-ords */

  mat_a2f(s->mat,m,c->basis,c->recip);
  for(i=0;i<3;i++){
    v[i]=0;
    if (s->tr)
      for(j=0;j<3;j++)
	v[i]+=s->tr[j]*c->recip[i][j];
  }
  

  for(i=0;i<3;i++){
    for(j=0;j<3;j++)
      if(!aeq(m[i][j],0)){
	if(aeq(m[i][j],1))fprintf(out,"+%c",'x'+j);
	else if(aeq(m[i][j],-1))fprintf(out,"-%c",'x'+j);
	else if(aeq(m[i][j],rint(m[i][j])))fprintf(out,"%+d%c",
						   (int)(rint(m[i][j])),'x'+j);
	else fprintf(out,"%+f%c",m[i][j],'x'+j);
      }
    if ((!aeq(v[i],0))&&(!aeq(v[i],1))) frac_print(out,v[i]);
    if(i<2) fprintf(out,",");
  }
  fprintf(out,"\n");
}

/* Read one of the three components of a symmetry matrix
 * in the "equivalent" form of -x,-y,z+1/2 etc.
 */

static void equiv_read_comp(char *str_in, double x[4]){
  double num,denom,sign,t;
  int i;
  char *str;

  for(i=0;i<4;i++) x[i]=0;
  denom=sign=1;
  num=0;
  str=str_in;

  while(*str){
    if (*str==' ') {str++; continue;}
    if (*str=='-'){sign=-1; str++; continue;}
    else if (*str=='+') {str++; continue;}
    /* Start parsing a number, as xxx.yyy, xxx/yyy or xxx */
    else if ((*str>='0')&&(*str<='9')){
      num=*str-'0';
      str++;
      while ((*str>='0')&&(*str<='9')){
	num=10*num+*str-'0';
	str++;
      }
      if (*str=='.'){
	t=0.1;
	str++;
	while ((*str>='0')&&(*str<='9')){
	  num+=t*(*str-'0');
	  t*=0.1;
	  str++;
	}
      }
      while (*str==' ') str++;
      if (*str=='/'){
	while (*str==' ') str++;
	str++;
	denom=0;
	while ((*str>='0')&&(*str<='9')){
	  denom=10*denom+*str-'0';
	  str++;
	}
      }
      while (*str==' ') str++;
    }
    /* End parsing number */
    i=3;
    if ((*str=='x')||(*str=='X')) {i=0; str++;}
    else if ((*str=='y')||(*str=='Y')) {i=1; str++;}
    else if ((*str=='z')||(*str=='Z')) {i=2; str++;}
    /* Next thing must be a number, or x-z, or end of string */
    while (*str==' ') str++;
    if ((*str)&(!strchr("0123456789+-xyzXYZ",*str))){
      fprintf(stderr,"Parsing of %s failed at %c\n",str_in,*str);
      exit(1);
    }
    if (num==0) x[i]=sign;
    else x[i]=sign*num/denom;
    denom=sign=1;
    num=0;
    
  }

  if (debug>2) fprintf(stderr,"Parsed %s as %lf %lf %lf %lf\n",str_in,
		       x[0],x[1],x[2],x[3]);

}

/* Read symmetry matrix in the "equivalent" form of
 * -x,-y,z+1/2 etc.
 */

void equiv_read(struct sym_op *s, struct unit_cell *c, char *line){
  double x[4],m[3][3],v[3];
  int i,j;
  char *p1,*p2;

  if (debug>2) fprintf(stderr,"equiv_read called with %s\n",line);

  p1=p2=line;
  if ((*p1=='"')||(*p1=='\'')){p1++;p2++;}
  for(i=0;i<3;i++){
    while((*p2!=',')&&(*p2!='\n')&&(*p2!='"')&&(*p2!='\'')&&(*p2)) p2++;
    if (((*p2==0)||(*p2=='"')||(*p2=='\''))&&(i!=2))
      error_exit("Unexpected end of string in equiv_read");
    *p2=0;
    equiv_read_comp(p1,x);
    for(j=0;j<3;j++) m[i][j]=x[j];
    v[i]=x[3];
    p2++;
    p1=p2;
  }

  if((v[0]!=0)||(v[1]!=0)||(v[2]!=0)){
    if (!s->tr) s->tr=malloc(3*sizeof(double));
    if (!s->tr) error_exit("malloc error in equiv_read");
    for(j=0;j<3;j++) s->tr[j]=v[0]*c->basis[0][j]+
		       v[1]*c->basis[1][j]+v[2]*c->basis[2][j];
  }

  mat_f2a(m,s->mat,c->basis,c->recip);

}

void ident_sym(struct sym_op *s, struct unit_cell *c, FILE *out){
  int i,j,inv,mult,screw;
  double m[3][3];
  double a[3],af[3],off[3];
  double v[3],v2[3];
  double det,angle,x;


  inv=0;
  screw=0;
  for(i=0;i<3;i++) off[i]=0;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      m[i][j]=s->mat[i][j];

  /* Calculate determinant */

  vcross(m[1],m[2],v);
  det=m[0][0]*v[0]+m[0][1]*v[1]+m[0][2]*v[2];

  if(!aeq(fabs(det),1)){
    fprintf(stderr,"Surprise: determinant is %f\n",det);
    fprintf(out,"(Error)\n");
    return;
  }

  if (aeq(det,1)){
    if (debug>2) fprintf(stderr,"Looks like a rotation matrix\n");
  }

  if (aeq(det,-1)){
    if (debug>2) fprintf(stderr,"Looks like an inverse rotation matrix\n");
    inv=1;
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	m[i][j]*=-1;
  }

  x=0.5*(m[0][0]+m[1][1]+m[2][2]-1);
  if (aeq(x,1.0)) x=1.0; /* Else rounding errors lead to acos(1+epsilon) */
  angle=acos(x);
  angle*=180/M_PI;

  if (debug>2) fprintf(stderr,"Angle looks like %f\n",angle);

  mult=0;
  if(aeq(angle,0)) mult=1;
  if(aeq(angle,180)) mult=2;
  if(aeq(angle,120)) mult=3;
  if(aeq(angle,90)) mult=4;
  if(aeq(angle,60)) mult=6;
  if (mult==0){
    fprintf(stderr,"Impossible angle %f in ident_sym\n",angle);
    fprintf(out,"(Error)\n");
    return;
  }

  /* Now need axis. It will be e-vector with e-value +1 */

  if (mult==1){
    v[0]=v[1]=0;
    v[2]=1;
  }
  else{
    for(i=0;i<3;i++)
      m[i][i]-=1;
    vcross(m[0],m[1],v);
    if (aeq(vmod2(v),0)) vcross(m[1],m[2],v);
    if (aeq(vmod2(v),0)) vcross(m[2],m[0],v);
    if (aeq(vmod2(v),0)){
      fprintf(stderr,"Surprise! 1\n");
      fprintf(out,"(Error)\n");
      return;
    }
  }

  for(i=0;i<3;i++) a[i]=v[i];
  if (debug>2) fprintf(stderr,"Axis looks like (%f,%f,%f)\n",a[0],a[1],a[2]);

  /* Try that in a fractional basis */

  for(i=0;i<3;i++)
    v2[i]=v[0]*c->recip[i][0]+v[1]*c->recip[i][1]+v[2]*c->recip[i][2];
  x=1e10;
  if ((fabs(v2[0])>tol)&&(fabs(v2[0])<x)) x=fabs(v2[0]);
  if ((fabs(v2[1])>tol)&&(fabs(v2[1])<x)) x=fabs(v2[1]);
  if ((fabs(v2[2])>tol)&&(fabs(v2[2])<x)) x=fabs(v2[2]);
  for(i=0;i<3;i++) v2[i]/=x;
  for(i=0;i<3;i++) a[i]/=x;
    
  for(i=0;i<3;i++) af[i]=v2[i];
  if (debug>2)
    fprintf(stderr,"Axis looks like (%f,%f,%f)\n",af[0],af[1],af[2]);

  /* Now deal with any offset */

  if (s->tr){
    for(i=0;i<3;i++) v[i]=s->tr[i];
    if(!(aeq(v[0],0)&&aeq(v[1],0)&&aeq(v[2],0))){
      x=v[0]*a[0]+v[1]*a[1]+v[2]*a[2];
      if((inv==0)&&(!aeq(x,0))){   /* we have a screw */
	x/=vmod2(a);  /* dot with unit vector parallel to axis, then divide by
		       * repeat length in that direction. Hence no sqrt. */
	for(i=0;i<3;i++) v[i]-=x*a[i];
	x*=mult;
	x=fmod(x,1);
	if (x<0) x+=1;
	if (x>1-tol) x=0;
	if (!aeq(x,floor(x+tol))){
	  fprintf(stderr,"Screw problem, x=%f mult=%d\n",x,mult);
	  fprintf(stderr,"Axis (%f,%f,%f) Disp (%f,%f,%f)\n",
		  a[0],a[1],a[2],s->tr[0],s->tr[1],s->tr[2]);
          fprintf(stderr,"Is (%f,%f,%f) a lattice vector?\n",
                  x*a[0],x*a[1],x*a[2]);
          screw=-1;
	}
        else screw=x+tol;
      }
    }
    if(!(aeq(v[0],0)&&aeq(v[1],0)&&aeq(v[2],0))){ /* we have a translation */
      //      fprintf(stderr,"Have tr=(%f,%f,%f)\n",v[0],v[1],v[2]);
      if((inv==1)&&(aeq(angle,0))){ /* An inversion point */
	for(i=0;i<3;i++) off[i]=0.5*v[i];
      }
      else if((inv==1)&&(aeq(angle,180))){ /* A mirror -- two evals are one */
	for(i=0;i<3;i++) off[i]=0.5*v[i];  /* Hmmm. Should take only bit
					    * parallel to a */
      }
      else if (aeq(angle,0)){ /* A pure translation */
	for(i=0;i<3;i++) off[i]=v[i];
      }
      else {
	for(i=0;i<3;i++)
	  for(j=0;j<3;j++)
	    m[i][j]=-s->mat[i][j];
	for(i=0;i<3;i++) m[i][i]+=1.0;
	/* We now have 1-R, and the translation is we have is the result
         * of (1-R)t, where t is the translation we want. Unfortunately,
         * if R is a proper rotation, it has an evalue of 1, so 1-R
         * has a zero evalue and is uninvertable. Fix this by adding a matrix
	 * with two zero evals, and one non-zero eval whose evec is along
	 * the rotation axis. This squashes the zero evalue, and leaves
	 * the others and their evectors alone.
	 */
	if (inv==0){
	  for(i=0;i<3;i++)
	    for(j=0;j<3;j++)
	      m[i][j]+=a[i]*a[j];
	}

	if(minvert(m)){
	  fprintf(stderr,"Surprise! Can't invert matrix\n");
	  fprintf(stderr,"Seem to have a %d axis along (%f,%f,%f)"
		  " inv=%d\n",mult,
		  af[0],af[1],af[2],inv);
	  fprintf(stderr,"Original matrix:\n");
	  for(i=0;i<3;i++)
	    fprintf(stderr,"( %6f %6f %6f )\n",s->mat[i][0],
		    s->mat[i][1],s->mat[i][2]);
	  fprintf(stderr,"Original translation:\n");
	  fprintf(stderr,"( %6f %6f %6f )\n",s->tr[0],s->tr[1],s->tr[2]);

	  return;
	}
	for(i=0;i<3;i++)
	  off[i]=m[i][0]*v[0]+m[i][1]*v[1]+m[i][2]*v[2];
      }
    }
  }

  if (inv) mult*=-1;
  for(i=0;i<3;i++){
    v[i]=off[0]*c->recip[i][0]+off[1]*c->recip[i][1]+off[2]*c->recip[i][2];
    v[i]=fmod(v[i],1);
    if (v[i]<0) v[i]+=1;
    if (v[i]>1-tol) v[i]=0;
  }

  if (mult==1){
    if ((s->tr)&&(!aeq(v[0]*v[0]+v[1]*v[1]+v[2]*v[2],0)))
	fprintf(out,"Translation of (%6.3f,%6.3f,%6.3f)\n",v[0],v[1],v[2]);
    else fprintf(out,"Identity\n");
    return;
  }

  if (screw)
    fprintf(out,"%d_%d axis along (%6.3f,%6.3f,%6.3f)",mult,screw,
	    af[0],af[1],af[2]);
  else
    fprintf(out,"%2d  axis along (%6.3f,%6.3f,%6.3f)",mult,
	    af[0],af[1],af[2]);

  if (aeq(v[0],0)&&aeq(v[1],0)&&aeq(v[2],0))
    fprintf(out," through (0,0,0)\n");
  else
    fprintf(out," through (%f,%f,%f)\n",v[0],v[1],v[2]);
}

void mpr(double m[3][3]){
  int i;
  for(i=0;i<3;i++)
    fprintf(stderr,"[%7f %7f %7f]\n",m[i][0],m[i][1],m[i][2]);
}

void vcross(double a[3],double b[3],double c[3]){
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-b[2]*a[0];
  c[2]=a[0]*b[1]-b[0]*a[1];
}

double vmod2(double v[3]){
  return (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

int minvert(double m[3][3]){
  int i,j;
  double det;
  double v[3],c[3][3];
#if 0
  double t[3][3];
#endif

  vcross(m[1],m[2],v);
  det=m[0][0]*v[0]+m[0][1]*v[1]+m[0][2]*v[2];
  if (aeq(det,0)) {
#if 0
    if (debug){
      fprintf(stderr,"Singular matrix:\n");
      mpr(m);
    }
#endif
    return 1;
  }
  /* Transpose of cofactor matrix */
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c[j][i]=m[(1+i)%3][(1+j)%3]*m[(2+i)%3][(2+j)%3]-
	m[(2+i)%3][(1+j)%3]*m[(1+i)%3][(2+j)%3];

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c[i][j]=c[i][j]/det;

#if 0
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      t[i][j]=c[i][0]*m[0][j]+c[i][1]*m[1][j]+c[i][2]*m[2][j];

  fprintf(stderr,"Identify in minv:\n");
  mpr(t);
#endif

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      m[i][j]=c[i][j];

  return 0;
}
