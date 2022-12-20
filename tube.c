/* Make nanotube
 * c must be perpendicular to a and b
 * the circumference vector must be in the ab plane, and is
 * rpt (in fractional coordinates)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "c2xsf.h"

void tube(struct unit_cell *c, struct contents *m, int rpt[3], double spacing){
  double new_basis[3][3],mn,mx,dtmp,amax;
  double r, theta, z, r0, abc[6],v1[2],v2[2],v2f[2],vtmp[3],vtmp2[3],ratio;
  double bsum,bwt,dot,cvec[3],stretch;
  int i,j,k,ii,jj,rep,orth[3];
  struct atom *atoms,*a;

  if (rpt[2]!=0) error_exit("Circumference vector not in ab plane\n");
  
  basis2abc(c->basis,abc);
  if ((!(aeq(abc[3],90)))||(!(aeq(abc[4],90))))
    error_exit("c not perpendicular to ab plane\n");

  for(i=0;i<3;i++)
    vtmp[i]=rpt[0]*c->basis[0][i]+rpt[1]*c->basis[1][i];
  
  if (debug)
    fprintf(stderr,"Circumference length = %lf A\n",sqrt(vmod2(vtmp)));

  amax=max(abc[0],abc[1]);
   
  /* See if system really is in ab plane, and, if so, where */

  mn=2;
  mx=-1;
  for(i=0;i<m->n;i++){
    dtmp=m->atoms[i].frac[2];
    dtmp=fmod(dtmp,1.0);
    if (dtmp<0.0) dtmp=dtmp+1.0;
    m->atoms[i].frac[2]=dtmp;
    if (dtmp>mx) mx=dtmp;
    if (dtmp<mn) mn=dtmp;
  }
  if ((1.0-(mx-mn))*abc[2]>4){  /* Look for 4A vacuum gap */
    if ((mx>=0.5)&&(mn<=0.5))
      ;
    else{
      dtmp=0.5-(mx+mn)/2;
      if (debug) fprintf(stderr,"Shifing atoms by %lf in c\n",dtmp);
      for(i=0;i<m->n;i++)
        m->atoms[i].frac[2]+=dtmp;
      addabs(m->atoms,m->n,c->basis);
    }
  }
  else{
    mn=2;
    mx=-1;
    for(i=0;i<m->n;i++){
      dtmp=m->atoms[i].frac[2]+0.5;
      dtmp=fmod(dtmp,1.0);
      if (dtmp<0.0) dtmp=dtmp+1.0;
      if (dtmp>mx) mx=dtmp;
      if (dtmp<mn) mn=dtmp;
    }
    if ((1.0-(mx-mn))*abc[2]>4){  /* Look for 4A vacuum gap */
      if ((mx>0.5)&&(mn<0.5)){
	if (debug) fprintf(stderr,"Shifing atoms by %lf in c\n",0.5);
        for(i=0;i<m->n;i++){
          dtmp=m->atoms[i].frac[2]+0.5;
          dtmp=fmod(dtmp,1.0);
          m->atoms[i].frac[2]=dtmp;
        }
        addabs(m->atoms,m->n,c->basis);
      }
      else{
        if (debug) fprintf(stderr,"Shifing atoms by %lf in c\n",dtmp);
        dtmp=1.0-(mx+mn)/2;
        for(i=0;i<m->n;i++)
          m->atoms[i].frac[2]=fmod(m->atoms[i].frac[2]+dtmp,1.0);
      }
      addabs(m->atoms,m->n,c->basis);
    }
    else{
      error_exit("Confused: cannot find vacuum gap in c\n");
    }
  }

  if (debug>2) print_cell(c,m);
  
  if (debug)
    fprintf(stderr,"Requested circumference vector is (%d,%d,%d)\n",
            rpt[0],rpt[1],rpt[2]);

  rep=1;
  j=min(rpt[0],rpt[1]);
  if (j==0) j=max(rpt[0],rpt[1]);
  for(i=2;i<=j;i++)
    if (((rpt[0]%i)==0)&&((rpt[1]%i)==0)){
      rpt[0]=rpt[0]/i;
      rpt[1]=rpt[1]/i;
      rep=rep*i;
      i=i-1;
    }

  if (debug)
    fprintf(stderr,"Factorised as %dx(%d,%d,%d)\n",rep,rpt[0],rpt[1],rpt[2]);


  /* Try to find orthogonal vec in ab plane */

  v1[0]=rpt[0]*c->basis[0][0]+rpt[1]*c->basis[1][0];
  v1[1]=rpt[0]*c->basis[0][1]+rpt[1]*c->basis[1][1];
  v2[0]=-v1[1];
  v2[1]=v1[0];
  v2f[0]=v2[0]*c->recip[0][0]+v2[1]*c->recip[0][1];
  v2f[1]=v2[0]*c->recip[1][0]+v2[1]*c->recip[1][1];

  if (debug)
    fprintf(stderr,"Orthogonal vector (fractional): (%lf,%lf,0)\n",
            v2f[0],v2f[1]);
  
  if (fabs(v2f[0])<tol) v2f[0]=0;
  if (fabs(v2f[1])<tol) v2f[1]=0;

  if ((v2f[0]!=0)&&(v2f[1]!=0)){ /* Need to reduce to integers */
    ratio=fabs(v2f[1]/v2f[0]);
    for(i=1;i<12;i++)
      if (aeq(i*ratio,(int)(i*ratio+.4))) break;
    if (i==12) error_exit("Search for orthogonal vector failed\n");
    if (v2f[0]<0) i=-i;
    orth[0]=i;
    orth[1]=floor(i*v2f[1]/v2f[0]+.4);
    if (debug)
      fprintf(stderr,"Cell vector parallel to tube (fractional): (%d,%d,0)\n",
              orth[0],orth[1]);
  }
  else{
    orth[0]=orth[1]=0;
    if (v2f[0]==0)
      orth[1]=(v2f[1]>0)?1:-1;
    else orth[0]=(v2f[0]>0)?1:-1;
  }
  if (((rpt[0]!=0)&&(rpt[1]!=0))||((orth[0]!=0)&&(orth[1]!=0))){
    /* Now make supercell using new vectors */
    for(j=0;j<3;j++)
      new_basis[0][j]=rpt[0]*c->basis[0][j]+rpt[1]*c->basis[1][j];
    for(j=0;j<3;j++)
      new_basis[1][j]=orth[0]*c->basis[0][j]+orth[1]*c->basis[1][j];
    for(j=0;j<3;j++)
      new_basis[2][j]=c->basis[2][j];
    if (debug>1){
      fprintf(stderr,"New basis\n");
      print_basis(new_basis);
    }
    
    super(c,m,new_basis,NULL,NULL,NULL,0);
  }
  
  basis2abc(c->basis,abc);
    
  if (spacing==0) spacing=ceil(rep*abc[0]/1.5);
  
  if (rep*abc[0]>2*M_PI*spacing)
    error_exit("Spacing too low for tube transform\n");

  if (debug>1)
    fprintf(stderr,"Circumference length = %lf A\n",rep*abc[0]);
  
  r0=rep*abc[0]/(2*M_PI*spacing);

  if (!(flags&RAW)){
  
    /* Now we'd like to find the average projection of a bond onto
     * the circumference so as to work out how to stretch the radius
     * to keep bond lengths mostly the same.
     */

    /* Assume an atom is bonded to all atoms with 1.3x the distance
     * of its closest neighbour
     */

    bsum=0;
    bwt=0;
    /* make cvec unit vector in dir of circumference */
    for(i=0;i<3;i++) cvec[i]=c->basis[0][i]/abc[0];
  
    for(i=0;i<m->n;i++){
      for(j=0;j<3;j++)
	vtmp[j]=m->atoms[i].abs[j]-0.5*c->basis[2][j];
      if (vmod2(vtmp)>2*amax*amax) continue;
      mn=1e30;
      for(j=0;j<m->n;j++){
	if (i!=j){
	  for(k=0;k<3;k++) vtmp[k]=m->atoms[i].abs[k]-m->atoms[j].abs[k];
	  if (vmod2(vtmp)<mn) mn=vmod2(vtmp);
	}
	else{
	  for(k=0;k<3;k++) vtmp[k]=0;
	}
	for(ii=-1;ii<=1;ii++)
	  for(jj=-1;jj<=1;jj++){
	    if ((ii==0)&&(jj==0)) continue;
	    for(k=0;k<3;k++)
	      vtmp2[k]=vtmp[k]+ii*c->basis[0][k]+jj*c->basis[1][k];
	    if (vmod2(vtmp2)<mn) mn=vmod2(vtmp2);
	  }
      }
      if (debug>2)
	fprintf(stderr,"Atom %d shortest neighbour distance %lf\n",i,sqrt(mn));
      for(j=0;j<m->n;j++){
	for(k=0;k<3;k++) vtmp[k]=m->atoms[i].abs[k]-m->atoms[j].abs[k];
	for(ii=-1;ii<=1;ii++){
	  for(jj=-1;jj<=1;jj++){
	    if ((i==j)&&(ii==0)&&(jj==0)) continue;
	    for(k=0;k<3;k++)
	      vtmp2[k]=vtmp[k]+ii*c->basis[0][k]+jj*c->basis[1][k];
	    if (vmod2(vtmp2)<1.3*1.3*mn){ /* Have bond */
	      dot=0;
	      for(k=0;k<3;k++)
		dot+=vtmp2[k]*cvec[k];
	      dot=fabs(dot);
	      /* dot is now projected length onto circumference */
	      bsum+=dot*dot/sqrt(vmod2(vtmp2));
	      bwt+=dot/sqrt(vmod2(vtmp2));
	    }
	  }
	}
      }
    }
    if (debug>1)
      fprintf(stderr,"Av weighted bond proj onto circ is %lf\n",bsum/bwt);
    theta=2*M_PI*(bsum/bwt)/(rep*abc[0]);
    stretch=theta/(2*sin(0.5*theta));
    if (debug) fprintf(stderr,"Stretch = %f\n",stretch);
    r0*=stretch;
  }
  
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      new_basis[i][j]=0;

  new_basis[0][0]=spacing;
  new_basis[1][1]=spacing;
  new_basis[2][2]=abc[1];

  atoms=malloc(m->n*rep*sizeof(struct atom));
  init_atoms(atoms,m->n*rep);
  
  for(i=0;i<rep;i++){
    for(j=0;j<m->n;j++){
      a=atoms+m->n*i+j;
      theta=2*M_PI*(m->atoms[j].frac[0]+i)/rep;
      z=m->atoms[j].frac[1];
      r=r0+(m->atoms[j].frac[2]-0.5)*(abc[2]/spacing);
      
      a->frac[0]=0.5+r*cos(theta);
      a->frac[1]=0.5+r*sin(theta);
      a->frac[2]=z;

      a->atno=m->atoms[j].atno;
      a->chg=m->atoms[j].chg;
      a->spin=m->atoms[j].spin;
      a->label=m->atoms[j].label;
    }
  }

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[i][j]=new_basis[i][j];

  free(m->atoms);
  m->n=m->n*rep;
  m->atoms=atoms;
  addabs(m->atoms,m->n,c->basis);
  fprintf(stderr,"Number of atoms in tube: %d\n",m->n);
  fprintf(stderr,"Unit cell length parallel to tube: %lf A \n",abc[1]);
  cell_check(c,m);
}
