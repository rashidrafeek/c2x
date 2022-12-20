/* Various utility functions for dealing with symmetry and k points */


/* Copyright (c) 2007, 2014, 2015 MJ Rutter 
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

void shorten(double vfinal[3][3]); /* from primitive.c */

/* Compare symmetry operations, rotation only */

int sym_in_list(struct sym_op *e, struct sym_op *s, int n){
  int i,j,k,hit;

  hit=0;
  for(i=0;i<n;i++){
    hit=1;
    for(j=0;j<3;j++)
      for(k=0;k<3;k++)
        if (fabs(s[i].mat[j][k]-e->mat[j][k])>tol) hit=0;

    if (hit==1) {break;}
  }
  return(hit); 
}


/* Real space syms to k-space syms -- strip translations and
 * add inversion */

void sym2ksym(struct symmetry *rs, struct symmetry *ks){
  int i,j,k,n;
  int inv=0;
  struct sym_op invert;

  for(j=0;j<3;j++)
    for(k=0;k<3;k++)
      if (j==k) invert.mat[j][k]=-1;
      else invert.mat[j][k]=0;
  
  /* If no sym ops given, initialise to identity and inversion and return */
  if (rs->ops==NULL){
    ks->ops=malloc(2*sizeof(struct sym_op));
    if (!ks->ops) error_exit("Malloc error in sym2ksym");
    ks->n=2;
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
        ks->ops[0].mat[i][j]=0;
    for(i=0;i<3;i++) ks->ops[0].mat[i][i]=1;
    ks->ops[0].tr=NULL;
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
        ks->ops[1].mat[i][j]=0;
    for(i=0;i<3;i++) ks->ops[1].mat[i][i]=-1;
    ks->ops[1].tr=NULL;
    if (debug>1) fprintf(stderr,"No symops passed to sym2ksym\n");
    return;
  }

  /* Does real space symmetry contain inversion? */

  inv=sym_in_list(&invert,rs->ops,rs->n);

  if (debug>1)
    fprintf(stderr,"Real space sym of %d ops does%s contain inversion\n",
	    rs->n,inv?"":" not");

  ks->ops=malloc((inv?2:1)*rs->n*sizeof(struct sym_op));
  if (!ks->ops) error_exit("Malloc error in sym2ksym");
  
  n=0;

  for(i=0;i<rs->n;i++){
    if (sym_in_list(rs->ops+i,ks->ops,n-1)==0){
      memcpy(ks->ops+n,rs->ops+i,sizeof(struct sym_op));
      ks->ops[n].tr=NULL;
      n++;
    }
    if (inv==0){
      for(j=0;j<3;j++)
	for(k=0;k<3;k++)
	  invert.mat[j][k]=-rs->ops[i].mat[j][k];
      if (sym_in_list(&invert,ks->ops,n-1)==0){
        memcpy(&invert,rs->ops+i,sizeof(struct sym_op));
	ks->ops[n].tr=NULL;
        n++;
      }
    }
  }
  
  ks->n=n;
  ks->tol=NULL;
  ks->gen=NULL;

  if (debug>1)
    fprintf(stderr,"Reciprocal space sym of %d ops\n",ks->n);

}

void sym_vec(struct atom *a, struct atom *b, struct sym_op *s,
	     double recip[3][3], int tr){

  int i,j;

  for(i=0;i<3;i++){
    b->abs[i]=0;
    for(j=0;j<3;j++)
      b->abs[i]+=s->mat[i][j]*a->abs[j];
  }

  if (tr && s->tr)
    for(i=0;i<3;i++)
      b->abs[i]+=s->tr[i];

  addfrac(b,1,recip);
}

void sym_atom(struct atom *a, struct atom *b, struct sym_op *s,
              double recip[3][3]){
  int i;

  b->atno=a->atno;
/* Unclear what the best thing is to do with spin. Leaving uninitialised
 * is not best */
  b->spin=a->spin;
  b->label=a->label;
  b->chg=a->chg;
  sym_vec(a,b,s,recip,1);

  for(i=0;i<3;i++){
    b->frac[i]=fmod(b->frac[i],1.0);
    if (b->frac[i]<0.0) b->frac[i]+=1.0;
    if (b->frac[i]>1.0-tol) b->frac[i]=0.0;
  }
}

int k_in_list(struct atom *k,struct kpts *kl){
  int hit,i;

  hit=-1;
  for(i=0;i<kl->n;i++){
    if ((dist(k->frac[0],kl->kpts[i].frac[0])<tol)&&
        (dist(k->frac[1],kl->kpts[i].frac[1])<tol)&&
        (dist(k->frac[2],kl->kpts[i].frac[2])<tol)){
      hit=i;
      break;
    }
  }

  return hit;
}

/* Remove symmetry-related kpoints from list */

void sym_kpts(struct kpts *k_in, struct kpts *k_out, struct symmetry *s,
              double basis[3][3]){
  struct atom k;
  int i,j,ii,n,hit;

  hit=0; /* Stop warning re uninitialised variable */

  if ((!s)||(s->n==0)) return;
  
  k_out->kpts=malloc(k_in->n*sizeof(struct atom));
  if (!k_out) error_exit("Malloc error in sym_kpts");

  for(i=0;i<3;i++) k_out->kpts[0].abs[i]=k_in->kpts[0].abs[i];
  for(i=0;i<3;i++) k_out->kpts[0].frac[i]=k_in->kpts[0].frac[i];
  for(i=0;i<3;i++){
    k_out->kpts[0].frac[i]=fmod(k_out->kpts[0].frac[i],1.0);
    if (k_out->kpts[0].frac[i]<0.0) k_out->kpts[0].frac[i]+=1.0;
    if (k_out->kpts[0].frac[i]>1.0-tol) k_out->kpts[0].frac[i]=0.0;
  }
  k_out->kpts[0].wt=k_in->kpts[0].wt;
  n=1;
  k_out->n=1;  

  for (i=1;i<k_in->n;i++){
    for(j=0;j<s->n;j++){
      sym_atom(k_in->kpts+i,&k,s->ops+j,basis);
      hit=k_in_list(&k,k_out);
      if (hit>=0) break;
    }
    if (hit>=0)
      k_out->kpts[hit].wt+=k_in->kpts[i].wt;
    else{
      for(ii=0;ii<3;ii++) k_out->kpts[n].abs[ii]=k_in->kpts[i].abs[ii];
      for(ii=0;ii<3;ii++) k_out->kpts[n].frac[ii]=k_in->kpts[i].frac[ii];
      k_out->kpts[n].wt=k_in->kpts[ii].wt;
      n++;
      k_out->n=n;
    }
  }

  k_out->n=n;
  k_out->mp=k_in->mp;

}

/* Generate "Monkhurst Pack" set according to Castep's rules */
void mp_gen(struct kpts *ks, struct unit_cell *c){
  int i,j,k,sym,sym2,max_x,hit;
  int ngx,ngy,ngz;
  double kx,ky,kz;
  double d[3];
  struct atom kinv;
  
  ngx=ks->mp->grid[0];
  ngy=ks->mp->grid[1];
  ngz=ks->mp->grid[2];

  for(i=0;i<3;i++){
    d[i]=ks->mp->disp[i];
    if ((ks->mp->grid[i]&1)==0) /* Grid even */
      d[i]+=1.0/(2*ks->mp->grid[i]); /* For thus says Castep */
    d[i]=fmod(d[i],1.0/ks->mp->grid[i]);
    if (d[i]<0) d[i]+=1.0/ks->mp->grid[i];
  }

  /* inversion symmetry will map points onto each other only
   * if all displacements are of zero or half a grid cell */
  sym=1;
  for(i=0;i<3;i++){
    if ((!aeq(d[i],0.0))&&(!aeq(ks->mp->grid[i]*d[i],0.5))&&
        (!aeq(ks->mp->grid[i]*d[i],1.0))) sym=0;
  }

  max_x=ngx;
  if (sym) max_x=floor((0.5-d[0])*ngx)+1;
  
  ks->kpts=malloc(max_x*ngy*ngz*sizeof(struct atom));
  if (!ks->kpts) error_exit("Malloc error in mp_gen");

  ks->n=0;
  for(i=0;i<max_x;i++){
    kx=(i)/(double)ngx+d[0];
    sym2=0;
    if ((sym)&&((aeq(kx,0.0))||(aeq(kx,0.5)))) sym2=1;
    for(j=0;j<ngy;j++){
      ky=(j)/(double)ngy+d[1];
      for(k=0;k<ngz;k++){
	kz=(k)/(double)ngz+d[2];
        if (sym2){
          kinv.frac[0]=-kx;
          kinv.frac[1]=-ky;
          kinv.frac[2]=-kz;
          hit=k_in_list(&kinv,ks);
          if (hit!=-1) {
            ks->kpts[hit].wt+=1.0/(ngx*ngy*ngz);
            continue;
          }
        }
        ks->kpts[ks->n].frac[0]=kx;
        ks->kpts[ks->n].frac[1]=ky;
        ks->kpts[ks->n].frac[2]=kz;
        ks->kpts[ks->n].wt=1.0/(ngx*ngy*ngz);
        if ((sym)&&(sym2==0))
          ks->kpts[ks->n].wt=2.0/(ngx*ngy*ngz);
	ks->n++;
      }
    }
  }

  addabs(ks->kpts,ks->n,c->recip);

  if (debug>1) fprintf(stderr,"Generated %d kpts, after inversion %d\n",
		       ngx*ngy*ngz,ks->n);
  
}

/* Integrate a given k-space wave using unsymmetrised weighted kpoints */
double kint(double x[3], struct kpts *ks){
  int i;
  double f,sum;

  if(debug>3) fprintf(stderr,"kint called with x=(%f,%f,%f)\n",x[0],x[1],x[2]);

  sum=0;
  for(i=0;i<ks->n;i++){
    f=cos(x[0]*ks->kpts[i].frac[0]+
          x[1]*ks->kpts[i].frac[1]+
          x[2]*ks->kpts[i].frac[2]);
    sum+=ks->kpts[i].wt*f;
  }
  return(sum);
}

/* Integrate a given k-space wave and all of its symmetry-related waves
 * using unsymmetrised weighted kpoints */
double kstar(double x[3],struct kpts *k, struct unit_cell *c,
	     struct symmetry *ks){
  int i,j,ii;
  double xa[3],t;
  double xsa[3],xs[3];

  if(debug>3) fprintf(stderr,"kstar called with x=(%f,%f,%f)\n",x[0],x[1],x[2]);

  for(i=0;i<3;i++){
    xa[i]=0;
    for(j=0;j<3;j++) xa[i]+=x[j]*c->recip[j][i];
  }

  t=0;
  for(ii=0;ii<ks->n;ii++){

    for(i=0;i<3;i++){
      xsa[i]=0;
      for(j=0;j<3;j++) xsa[i]+=ks->ops[ii].mat[i][j]*xa[j];
    }

    for(i=0;i<3;i++){
      xs[i]=0;
      for(j=0;j<3;j++) xs[i]+=xsa[j]*c->basis[i][j];
    }

    t+=kint(xs,k);
  }

  return(t/ks->n);
}

/* Calculate first failure star */
void fstar(struct kpts *kp, struct unit_cell *c_in, struct symmetry *rs){
  int i,j,k,ii;
  int nx,ny,nz,f[3];
  double x[3],t;
  double f2,x2,tf;
  struct symmetry ks;
  struct unit_cell *c;

  f[0]=f[1]=f[2]=0; /* Stop warning re uninitialised variable */
  tf=1e300;

  if (kp->n==0){
    fprintf(stderr,"No k-points given: no failure star\n");
    return;
  }

  /* Work in a sanitised basis
   * N.B. sym ops are independent of basis choice
   */

  c=malloc(sizeof(struct unit_cell));
  if (!c) error_exit ("Malloc error in fstar");
  c->basis=malloc(9*sizeof(double));
  if (!c->basis) error_exit ("Malloc error in fstar");
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[i][j]=c_in->basis[i][j];
  shorten(c->basis);
  real2rec(c);

  sym2ksym(rs,&ks);


  x[0]=x[1]=x[2]=0;
  for(nx=1;;nx++) {
    x[0]=2*nx*M_PI;
    t=kstar(x,kp,c,&ks);
    if (fabs(t)>tol) break;
  }
  if (debug>2) fprintf(stderr,"Fail at nx=%d, int=%f\n",nx,t);

  x[0]=x[1]=x[2]=0;
  for(ny=1;;ny++) {
    x[1]=2*ny*M_PI;
    t=kstar(x,kp,c,&ks);
    if (fabs(t)>tol) break;
  }
  if (debug>2) fprintf(stderr,"Fail at ny=%d, int=%f\n",ny,t);

  x[0]=x[1]=x[2]=0;
  for(nz=1;;nz++) {
    x[2]=2*nz*M_PI;
    t=kstar(x,kp,c,&ks);
    if (fabs(t)>tol) break;
  }
  if (debug>2) fprintf(stderr,"Fail at nz=%d, int=%f\n",nz,t);

  /* I believe that, if the cell is reduced to a compact form,
   * then iterating over all points with the first index between
   * -a*nx and a*nx, second index between -a*ny and a*ny, etc, with
   * a=1/cos(30) will include all vectors in the elipsoid formed by
   * nx, ny and nz. */ 

  nx*=1.155;
  ny*=1.155;
  nz*=1.155;

  f2=1e300;
  for(i=-nx;i<=nx;i++){
    x[0]=2*i*M_PI;
    for(j=-ny;j<=ny;j++){
      x[1]=2*j*M_PI;
      for(k=-nz;k<=nz;k++){
        if(i+j+k<0) continue;
        if((i==0)&&(j==0)&&(k==0)) continue;
        x[2]=2*k*M_PI;
	t=kstar(x,kp,c,&ks);
	if (fabs(t)>tol){
	  /* Need mod2(x) in absolute co-ords */
	  x2=0;
	  for(ii=0;ii<3;ii++)
	    x2+=pow(c->basis[0][ii]*i+c->basis[1][ii]*j+c->basis[2][ii]*k,2);
          if (x2<f2){
	    f2=x2;
	    f[0]=i;
	    f[1]=j;
	    f[2]=k;
            tf=t;
	  }
	}
      }
    }
  }

  if (f2>1e250) fprintf(stderr,"Impossible: search for failure star failed!\n");

  for(i=0;i<3;i++)
    x[i]=c->basis[0][i]*f[0]+c->basis[1][i]*f[1]+c->basis[2][i]*f[2];

  fprintf(stderr,"First failure star, integral %f\n",tf);
  fprintf(stderr,"Typical vector: (%.6f,%.6f,%.6f), length**2=%f\n",x[0],
	  x[1],x[2],f2);

  free(c->basis);
  free(c);
}

/* Check that all permutations of (1,0,0), (1,1,0) and (1,1,1) map
 * under symmetry (excluding translations) to vectors whose relative
 * components are zero or +/-1
 */

void sym_basis(struct symmetry *s, struct unit_cell *c){
  int i,j,k,del;
  struct atom v1,v2;

  for(i=0;i<s->n;i++){
    del=0;
    for(j=0;j<3;j++){
      for(k=0;k<3;k++) v1.frac[k]=0;
      v1.frac[j]=1;
      addabs(&v1,1,c->basis);
      sym_vec(&v1,&v2,s->ops+i,c->recip,0);
      for(k=0;k<3;k++){
	/* Components of v2 must be 0 or +/-1 */
	v2.frac[k]=fabs(v2.frac[k]);
	if ((v2.frac[k]<tol)||(fabs(1-v2.frac[k])<tol)){;}
        else del=1;
      }
      if (del) break;
      for(k=0;k<3;k++) v1.frac[k]=1;
      v1.frac[j]=0;
      addabs(&v1,1,c->basis);
      sym_vec(&v1,&v2,s->ops+i,c->recip,0);
      for(k=0;k<3;k++){
	/* Components of v2 must be 0 or +/-1 */
	v2.frac[k]=fabs(v2.frac[k]);
	if ((v2.frac[k]<tol)||(fabs(1-v2.frac[k])<tol)){;}
        else {del=1; break;}
      }
    }
    if (del==0){
      for(k=0;k<3;k++) v1.frac[k]=1;
      addabs(&v1,1,c->basis);
      sym_vec(&v1,&v2,s->ops+i,c->recip,0);
      for(k=0;k<3;k++){
	/* Components of v2 must be 0 or +/-1 */
	v2.frac[k]=fabs(v2.frac[k]);
	if ((v2.frac[k]<tol)||(fabs(1-v2.frac[k])<tol)){;}
        else {del=1; break;}
      }
    }
    if (del){
      if (debug) fprintf(stderr,
			 "Removing sym op as incompatible with unit cell\n");
      if (debug>1){
	ident_sym(s->ops+i,c,stderr);
        fprintf(stderr,"(%f,%f,%f) -> (%f,%f,%f)\n",
		v1.frac[0],v1.frac[1],v1.frac[2],
		v2.frac[0],v2.frac[1],v2.frac[2]);
      }

      for(j=i+1;j<s->n;j++){
	memcpy(s->ops+j-1,s->ops+j,sizeof(struct sym_op));
      }
      s->n--;
      i--;
    }
  }
}

/* Convert m1 in fractional co-ords to m2 in absolute */

void mat_f2a(double m1[3][3], double m2[3][3], double basis[3][3],
	    double recip[3][3]){
  int i,j,k;
  double m[3][3];
  
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      m[i][j]=m2[i][j]=0.0;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      for(k=0;k<3;k++)
	m[i][j]+=m1[i][k]*recip[k][j];

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      for(k=0;k<3;k++)
	m2[i][j]+=basis[k][i]*m[k][j];

}

/* Convert m1 in absolute co-ords to m2 in fractional */

void mat_a2f(double m1[3][3], double m2[3][3], double basis[3][3],
	    double recip[3][3]){
  int i,j,k;
  double m[3][3];
  
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      m[i][j]=m2[i][j]=0.0;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      for(k=0;k<3;k++)
	m[i][j]+=m1[i][k]*basis[j][k];

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      for(k=0;k<3;k++)
	m2[i][j]+=recip[i][k]*m[k][j];

}
