/* Find primitive unit cell */


/* Copyright (c) 2010 MJ Rutter 
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

double volume(double b[3][3]);
void basis_signs(double basis[3][3]);
void basis_equal(double basis[3][3]);
int is_rhs(double b[3][3]); /* function found in cell.basis.c */
int rarest(struct unit_cell *c, struct contents *m, int *spec);

int lv_sort(const void *a, const void *b){ /* This ghastly declaration */
  struct vector *v1 = (struct vector *) a; /* stops compilers moaning */
  struct vector *v2 = (struct vector *) b; /* about qsort's prototype */

  if (v1->mod2<v2->mod2) return -1;
  if (v1->mod2>v2->mod2) return 1;
  return 0;
}


void primitive(struct unit_cell *c, struct contents *m, double vfinal[3][3]){
  double vec[3],tvec[3];
  double ab,a2,b2,prim_vol;
  struct vector *lv;
  int i,j,k,found,ok;
  int rare_spec,nrare,nlattice_vecs;
  int *lattice_vecs;

  /* Fill in default answer */

  for (i=0;i<3;i++)
    for(j=0;j<3;j++)
      vfinal[i][j]=c->basis[i][j];

  prim_vol=c->vol;

  found=rarest(c,m,&rare_spec);

  nrare=0;
  for(i=0;i<m->n;i++) if (m->atoms[i].atno==rare_spec) nrare++;

  for(i=0;i<3;i++) tvec[i]=fmod(m->atoms[found].frac[i]+1.5,1.0)-0.5;

  if (debug>1)
    fprintf(stderr,"Shift of origin (fractional co-ords) "
	    "(%6f,%6f,%6f)\n",tvec[0],tvec[1],tvec[2]);


  for(i=0;i<m->n;i++)
    for(j=0;j<3;j++)
      m->atoms[i].frac[j]=fmod(1.0+m->atoms[i].frac[j]-tvec[j],1.0);

  /* Remove all fract co-ords which are 1.0 */

  for(i=0;i<m->n;i++)
    for(j=0;j<3;j++)
      if (m->atoms[i].frac[j]>1-tol) m->atoms[i].frac[j]=0.0;

  addabs(m->atoms,m->n,c->basis);

  if (nrare==1) return;

  lattice_vecs=malloc(nrare*sizeof(int));
  nlattice_vecs=0;

  for(i=0;i<m->n;i++){
    if (m->atoms[i].atno!=rare_spec) continue;
    if (i==found) continue;
    for(j=0;j<3;j++) vec[j]=m->atoms[i].frac[j];
    if (debug>2) fprintf(stderr,"Trying (%f,%f,%f) as lattice vector\n",
			 vec[0],vec[1],vec[2]);
    ok=1;
    for(k=0;k<m->n;k++){
      for(j=0;j<3;j++) {
	tvec[j]=fmod(1+vec[j]+m->atoms[k].frac[j],1.0);
        if (tvec[j]>1.0-tol) tvec[j]=0.0;
      }
      /* If this is a lattice vector, tvec[] should have hit another atom */
      for(j=0;j<m->n;j++)
        if(((m->atoms[j].frac[0]-tvec[0])*(m->atoms[j].frac[0]-tvec[0])+
	    (m->atoms[j].frac[1]-tvec[1])*(m->atoms[j].frac[1]-tvec[1])+
	    (m->atoms[j].frac[2]-tvec[2])*(m->atoms[j].frac[2]-tvec[2]))
	   <tol*tol) break;
      if (j==m->n) {
        ok=0;
        if (debug>2) fprintf(stderr,"Atom %d maps to nothing\n",k);
        if (debug>3) fprintf(stderr,"nothing found at (%f,%f,%f)\n",
			     tvec[0],tvec[1],tvec[2]);
        break;
      }
      if (m->atoms[j].atno!=m->atoms[k].atno){
        ok=0;
        if (debug>2) fprintf(stderr,"Atom %d species %d maps to species %d\n",
			     k,m->atoms[k].atno,m->atoms[j].atno);
        if (debug>3) fprintf(stderr,"at (%f,%f,%f)\n",
			     tvec[0],tvec[1],tvec[2]);
	break;
      }
      if (fabs(m->atoms[j].spin-m->atoms[k].spin)>tol){
	ok=0;
	if (debug>1)
	  fprintf(stderr,"Atom %d species %d maps to atom of same species"
		  " but different spin\n",k,m->atoms[k].atno);
	break;
      }
      /* Atoms the same if their label pointers are equal (possibly null),
         or their label pointers are not null but point to strings which
         are equal to strcmp */
      if (!((m->atoms[j].label==m->atoms[i].label)||
            ((m->atoms[j].label)&&(m->atoms[i].label)&&
             (!strcmp(m->atoms[j].label,m->atoms[i].label))))){
	ok=0;
	if (debug>1)
	  fprintf(stderr,"Atom %d species %d maps to atom of same species"
		  " but different label\n",k,m->atoms[k].atno);
	break;
      }
    }
    if (ok){
      lattice_vecs[nlattice_vecs]=i;
      nlattice_vecs++;
    }
  }

  if (nlattice_vecs>0) {

    prim_vol=c->vol/(nlattice_vecs+1);
    if (debug>1){
      fprintf(stderr,"Found extra lattice vectors\n");
      for(i=0;i<nlattice_vecs;i++) fprintf(stderr,"(%f,%f,%f)\n",
  				 m->atoms[lattice_vecs[i]].frac[0],
  				 m->atoms[lattice_vecs[i]].frac[1],
  				 m->atoms[lattice_vecs[i]].frac[2]);
      fprintf(stderr,"Initial cell volume %f\n",c->vol);
      fprintf(stderr,"Primitive cell volume should be %f\n",prim_vol);
    }
    /* Form array of lattice vecs in absolute co-ords */
  
    lv=malloc((3+nlattice_vecs)*sizeof(struct vector));
    for(i=0;i<nlattice_vecs;i++){
      for(j=0;j<3;j++){
        lv[i].v[j]=0;
  	for(k=0;k<3;k++)
  	  lv[i].v[j]+=m->atoms[lattice_vecs[i]].frac[k]*c->basis[k][j];
      }
      lv[i].mod2=0;
      for(j=0;j<3;j++) lv[i].mod2+=lv[i].v[j]*lv[i].v[j];
    }
  
    /* And add original c->basis vectors */
  
    for(i=0;i<3;i++){
      for(j=0;j<3;j++) lv[nlattice_vecs+i].v[j]=c->basis[i][j];
      lv[nlattice_vecs+i].mod2=0;
      for(j=0;j<3;j++) lv[nlattice_vecs+i].mod2+=lv[nlattice_vecs+i].v[j]*
  		       lv[nlattice_vecs+i].v[j];
    }
  
    nlattice_vecs+=3;
  
    /* Sort them */
    if (debug>2){
      fprintf(stderr,"Unsorted lattice vectors\n");
      for(i=0;i<nlattice_vecs;i++) fprintf(stderr,"(%f,%f,%f)  %f\n",
  					 lv[i].v[0],lv[i].v[1],lv[i].v[2],
  					 lv[i].mod2);
    }
  
    qsort(lv,(size_t)nlattice_vecs ,sizeof(struct vector),lv_sort);
  
    if (debug>2){
      fprintf(stderr,"Sorted lattice vectors\n");
      for(i=0;i<nlattice_vecs;i++) fprintf(stderr,"(%f,%f,%f)  %f\n",
  					 lv[i].v[0],lv[i].v[1],lv[i].v[2],
  					 lv[i].mod2);
    }
    /* Take first */
    for(j=0;j<3;j++) vfinal[0][j]=lv[0].v[j];
  
    /* Take next which is not parallel to first */
  
    for(i=1;i<nlattice_vecs;i++){
      ab=0;
      for(j=0;j<3;j++)ab+=vfinal[0][j]*lv[i].v[j];
      a2=lv[0].mod2;
      b2=lv[i].mod2;
      if ((a2*b2-ab*ab)>tol*a2*b2) break;
    }
  
    for(j=0;j<3;j++) vfinal[1][j]=lv[i].v[j];
  
    /* Take next which gives correct volume */
  
    for(i++;i<nlattice_vecs;i++){
      for(j=0;j<3;j++) vfinal[2][j]=lv[i].v[j];
      if (aeq(volume(vfinal),prim_vol)) break;
    }
  
    if (!aeq(volume(vfinal),prim_vol)){
      fprintf(stderr,
	      "Horror: failed to find lattice vectors for primitive cell\n");
      exit(1);
    }

    /* Success */
  
    if (debug>2){
      fprintf(stderr,"Final lattice vectors are:\n");
      for(i=0;i<3;i++) fprintf(stderr,"(%f,%f,%f)\n",
  			     vfinal[i][0],vfinal[i][1],vfinal[i][2]);
      fprintf(stderr,"Final volume is %f\n",volume(vfinal));
    }

  }  /* end if (nlattice_vecs>0) */  


}

void shorten(double vfinal[3][3]){
  double fix,ab,b2,len[3],tvec[3];
  int iter,i,j,k,ii;

  /* Worry about shortening vectors */
  fix=1;
  iter=0;
  while((fix)&&(iter<20)){
    fix=0;
    iter++;
    for(i=0;i<3;i++){
    /* Look at other two vectors */
      for(j=0;j<3;j++){
	if (i==j) continue;
	ab=0;
	for(k=0;k<3;k++) ab+=vfinal[i][k]*vfinal[j][k];
	b2=0;
	for(k=0;k<3;k++) b2+=vfinal[j][k]*vfinal[j][k];
	if (ab>0.5*b2) {
	  fix=(int)(ab/b2+0.499);
	  if (debug>2) fprintf(stderr,"Correcting %d-=%f * %d\n",i,fix,j);
	  for(k=0;k<3;k++) vfinal[i][k]-=fix*vfinal[j][k];
	}
	else if (-ab>0.5*b2) {
	  fix=(int)(-ab/b2+0.499);
	  if (debug>2) fprintf(stderr,"Correcting %d+=%f * %d\n",i,fix,j);
	  for(k=0;k<3;k++) vfinal[i][k]+=fix*vfinal[j][k];
	}
      }

      /* Then their sum/difference as req. for non-orthogonal vectors */
      j=(i+1)%3;
      k=(i+2)%3;
      ab=0;
      for(ii=0;ii<3;ii++) ab+=vfinal[j][ii]*vfinal[k][ii];
      if (ab>0)
	for(ii=0;ii<3;ii++) tvec[ii]=vfinal[j][ii]-vfinal[k][ii];
      else
	for(ii=0;ii<3;ii++) tvec[ii]=vfinal[j][ii]+vfinal[k][ii];
      
      ab=0;
      for(ii=0;ii<3;ii++) ab+=vfinal[i][ii]*tvec[ii];
      b2=0;
      for(ii=0;ii<3;ii++) b2+=tvec[ii]*tvec[ii];
      if (ab>0.5*b2) {
	fix=(int)(ab/b2+0.499);
	if (debug>2) fprintf(stderr,"Correcting %d-=%f * (%d %d)\n",i,fix,j,k);
	for(k=0;k<3;k++) vfinal[i][k]-=fix*tvec[k];
      }
      else if (-ab>0.5*b2) {
	fix=(int)(-ab/b2+0.499);
	if (debug>2) fprintf(stderr,"Correcting %d+=%f * (%d %d)\n",i,fix,j,j);
	for(k=0;k<3;k++) vfinal[i][k]+=fix*tvec[k];
      }

    }
  }
  if (iter>=20){
    fprintf(stderr,"Problem: infinite loop in vector shortening code\n");
    fprintf(stderr,"Troublesome vectors are:\n");
    for(i=0;i<3;i++) fprintf(stderr,"(%f,%f,%f)\n",
			     vfinal[i][0],vfinal[i][1],vfinal[i][2]);
    fprintf(stderr,"Carrying on regardless.\n");
  }


  if (debug>2){
    fprintf(stderr,"Final, final lattice vectors are:\n");
    for(i=0;i<3;i++) fprintf(stderr,"(%f,%f,%f)\n",
			     vfinal[i][0],vfinal[i][1],vfinal[i][2]);
    fprintf(stderr,"Final,final volume is %f\n",volume(vfinal));
  }

  /* Canonicalise slightly */

  basis_signs(vfinal);
  for(i=0;i<3;i++){
    len[i]=0;
    for(j=0;j<3;j++) len[i]+=vfinal[i][j]*vfinal[i][j];
  }

  if((aeq(len[0],len[1]))&&(aeq(len[1],len[2]))) basis_equal(vfinal);



  if (debug>2){
    fprintf(stderr,"Final, final,final lattice vectors are:\n");
    for(i=0;i<3;i++) fprintf(stderr,"(%f,%f,%f)\n",
			     vfinal[i][0],vfinal[i][1],vfinal[i][2]);
    fprintf(stderr,"Final, final, final volume is %f\n",volume(vfinal));
  }


}

int rarest(struct unit_cell *c, struct contents *m,int *spec){
  double tvec[3],dtmp,dist2;
  int i,j,k,found;
  int nspec,rare_spec,nrare;
  int *spec_at,*spec_n,*atom_spec;

  /* Find least frequently occuring species */

  spec_n=malloc(m->n*sizeof(int));
  spec_at=malloc(m->n*sizeof(int));
  atom_spec=malloc(m->n*sizeof(int));

  for(i=0;i<m->n;i++) spec_n[i]=spec_at[i]=0;
  nspec=0;

  for(i=0;i<m->n;i++){
    found=0;
    for(j=0;j<nspec;j++){
      if (m->atoms[i].atno==spec_at[j]){
        spec_n[j]++;
        found=1;
        break;
      }
    }
    if (found==0){
      spec_at[nspec]=m->atoms[i].atno;
      spec_n[nspec]=1;
      atom_spec[i]=nspec;
      nspec++;
    }
    else atom_spec[i]=j;
  }

  if (debug>2){
    fprintf(stderr,"Cell contains %d species\n",nspec);
    for(i=0;i<nspec;i++)
      fprintf(stderr,"Species %3d  atomic number %3d  m->atoms %d\n",
	      i,spec_at[i],spec_n[i]);
  }

  nrare=spec_n[0];
  rare_spec=spec_at[0];
  for(i=1;i<nspec;i++){
    if(spec_n[i]<nrare){
      nrare=spec_n[i];
      rare_spec=spec_at[i];
    }
  }

  if (debug>1) fprintf(stderr,"(Joint?) rarest species is atomic number %d"
		       " occuring %d times in the cell\n",rare_spec,nrare);

  /* Place rarest atom at origin */

  /* Find one closest to origin */

  found=-1;
  dist2=1e100;
  for (i=0;i<m->n;i++){
    if (spec_n[atom_spec[i]]==nrare){
      dtmp=0;
      for(j=0;j<3;j++){
	tvec[j]=0;
	for(k=0;k<3;k++)
	  tvec[j]+=(fmod(m->atoms[i].frac[k]+11.5,1.0)-0.5)*c->basis[k][j];
      }
      for(j=0;j<3;j++) dtmp+=tvec[j]*tvec[j];
      if (dtmp<dist2){
	dist2=dtmp;
	found=i;
      }
    }
  }

  if (spec_at[atom_spec[found]]!=rare_spec){
    if (debug>1) fprintf(stderr,"Changing rarest species to atomic no %d\n",
			 spec_at[atom_spec[found]]);
    rare_spec=spec_at[atom_spec[found]];
  }

  if (found==-1) error_exit("Impossible in primitive, found==-1");

  if (debug>1) fprintf(stderr,"Closest rare to origin atom no %d "
		       "atno %d dist2 %f\n",found,atom_spec[found],dist2);

  *spec=rare_spec;

  free(spec_n);
  free(spec_at);
  free(atom_spec);

  return found;

}


void basis_signs(double b[3][3]){
  int i;
  double a2,b2,c2,ab,tvec[3];

  if (b[0][0]*b[1][0]+b[0][1]*b[1][1]+
      b[0][2]*b[1][2]<0) for(i=0;i<3;i++) b[1][i]*=-1;

  a2=b2=c2=0;
  for (i=0;i<3;i++) {
    a2+=b[0][i]*b[0][i];
    b2+=b[1][i]*b[1][i];
  }
  a2=sqrt(a2);
  b2=sqrt(b2);

  for(i=0;i<3;i++) tvec[i]=b[0][i]/a2+b[1][i]/b2;

  if (tvec[0]*b[2][0]+tvec[1]*b[2][1]+
      tvec[2]*b[2][2]<0) for(i=0;i<3;i++) b[2][i]*=-1;

  if (!is_rhs(b)){
    for(i=0;i<3;i++) {
      ab=b[1][i];
      b[1][i]=b[2][i];
      b[2][i]=ab;
    }
  }

}


void basis_equal(double b[3][3]){
  /* See if we can find another set of equilengthed vectors which
   * lead to alpha=beta=gamma */
  int i,j,k,ii;
  double tvec[3],t2;
  double abc[6],nabc[6],nb[3][3];

  if (debug>1) fprintf(stderr,"basis_equal called\n");

  /* First convert to a,b,c alpha, beta, gamma */
  basis2abc(b,abc);

  if (!(aeq(abc[0],abc[1])&&aeq(abc[1],abc[2]))) return;

  if ((aeq(abc[3],abc[4]))&&(aeq(abc[4],abc[5]))) return;

  if (debug>2){
    for(i=0;i<3;i++) fprintf(stderr,"(%f,%f,%f)\n",
			     b[i][0],b[i][1],b[i][2]);
    fprintf(stderr,"Basis equal length: %f\n",abc[0]);
    fprintf(stderr,"Basis equal alpha, beta, gamma: %f %f %f\n",
	    abc[3],abc[4],abc[5]);
  }

  for(i=0;i<3;i++) tvec[i]=b[0][i]+b[1][i]+b[2][i];
  t2=0;
  for(i=0;i<3;i++) t2+=tvec[i]*tvec[i];
  t2=sqrt(t2);

  if(debug>2) fprintf(stderr,"t2=%f\n",t2);

  if(aeq(t2,abc[0])) return;

  for(i=0;i<3;i++){
    j=(i+1)%3;
    k=(i+2)%3;
    if (!(aeq(abc[j+3],abc[k+3]))) continue;

    for (j=0;j<3;j++)
      for(k=0;k<3;k++)
	nb[j][k]=b[j][k];
    for(j=0;j<3;j++) nb[i][j]=-nb[i][j];

    if (debug>2)
      for(ii=0;ii<3;ii++) fprintf(stderr,"(%f,%f,%f)\n",
				  nb[ii][0],nb[ii][1],nb[ii][2]);

    basis2abc(nb,nabc);
    if (debug>2) fprintf(stderr,"Basis equal t alpha, beta, gamma: %f %f %f\n",
	    nabc[3],nabc[4],nabc[5]);

    for(ii=0;ii<3;ii++) tvec[ii]=nb[0][ii]+nb[1][ii]+nb[2][ii];
    t2=0;
    for(ii=0;ii<3;ii++) t2+=tvec[ii]*tvec[ii];
    t2=sqrt(t2);

    if(debug>2) fprintf(stderr,"t t2=%f\n",t2);

    if (aeq(nabc[0],t2)){
    for (j=0;j<3;j++)
      for(k=0;k<3;k++)
	b[j][k]=nb[j][k];
    if (debug>1) fprintf(stderr,"Change in bravais_equal\n");
    return;
    }
  }
}


double volume(double b[3][3]){
  /* Return cell volume, i.e. a.(bxc) */
  double vol;

  vol=b[0][0]*(b[1][1]*b[2][2]-b[1][2]*b[2][1])+
      b[0][1]*(-b[1][0]*b[2][2]+b[1][2]*b[2][0])+
      b[0][2]*(b[1][0]*b[2][1]-b[1][1]*b[2][0]);

  return fabs(vol);
}

 
