/* Interface between c2x and spglib */

/* Copyright (c) 2014 -- 2022 MJ Rutter 
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

#ifdef SPGLIB
#include <spglib.h>
#endif

#include "c2xsf.h"

int is_rhs(double b[3][3]); /* From basis.c */

void rotate(double v[3], double axis[3], double theta);
struct sym_op *sym_frac2abs(int spg_rot[][3][3],double spg_tr[][3],
			    struct unit_cell *c,int nsym);

int igr2hall[]={ 0,    1,  2,  3,  6,  9, 18, 21, 30, 39,
                 57,  60, 63, 72, 81, 90,108,109,112,115,
                 116,119,122,123,124,125,128,134,137,143,
                 149,155,161,164,170,173,176,182,185,191,
                 197,203,209,212,215,218,221,227,228,230,
                 233,239,245,251,257,263,266,269,275,278,
                 284,290,292,298,304,310,313,316,322,334,
                 335,337,338,341,343,349,350,351,352,353,
                 354,355,356,357,358,359,361,363,364,366,
                 367,368,369,370,371,372,373,374,375,376,
                 377,378,379,380,381,382,383,384,385,386,
                 387,388,389,390,391,392,393,394,395,396,
                 397,398,399,400,401,402,404,406,407,408,
                 410,412,413,414,416,418,419,420,422,424,
                 425,426,428,430,431,432,433,435,436,438,
                 439,440,441,442,443,444,446,447,448,449,
                 450,452,454,455,456,457,458,460,462,463,
                 464,465,466,467,468,469,470,471,472,473,
                 474,475,476,477,478,479,480,481,482,483,
                 484,485,486,487,488,489,490,491,492,493,
                 494,495,497,498,500,501,502,503,504,505,
                 506,507,508,509,510,511,512,513,514,515,
                 516,517,518,520,521,523,524,525,527,529,
                 530};

#ifdef SPGLIB

static void spg_sym_print(SpglibDataset *spg,double spg_symprec, int op);
static void spg_scan(double spg_latt[3][3], double spg_pos[][3], int *spg_type,
		     int natoms, int op, double tolmin, double tolmax,
		     int nsmin, int nsmax, int first);

int cspg_op(struct unit_cell *c, struct contents *m, struct symmetry *s,
            struct kpts *kp, int op, double tolmin){
  int i,j,k,natoms2;
  int rotation;
  double dot,mod1,mod2,theta,axis[3],vtmp[3],vtmp1[3],vtmp2[3];
  double old_basis[3][3],old_frac[3][3],new_basis[3][3];
  double rhs[3][3],tr[3];
  struct contents *m2;
  struct unit_cell *c2;
  struct symmetry stmp;
  
  /* SPG variables */
  double spg_latt[3][3];
  double (*spg_pos)[3],spg_symprec;
  int *spg_type;
  SpglibDataset *spg;

  struct uniq {int atno; double spin; char *label;} *auid;
  int hit,n_uniq;

  if (debug>1) fprintf(stderr,"Calling spglib with op=0x%x\n",op);

  if ((op&CSPG_SNAP)||(op&CSPG_PRIM_NR)){
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	old_basis[i][j]=c->basis[i][j];
  }
  
  if (op&CSPG_SNAP){
    m2=malloc(sizeof(struct contents));
    if (!m2) error_exit("Malloc error for m2");
    memcpy(m2,m,sizeof(struct contents));
    m2->atoms=malloc(m2->n*sizeof(struct atom));
    if (!m2->atoms) error_exit("Malloc error for m2->atoms");
    memcpy(m2->atoms,m->atoms,m2->n*sizeof(struct atom));
    c2=malloc(sizeof(struct unit_cell));
    if (!c2) error_exit("Malloc error for c2");
    memcpy(c2,c,sizeof(struct unit_cell));
    c2->basis=malloc(9*sizeof(double));
    if (!c2->basis) error_exit("Malloc error for c2->basis");
    memcpy(c2->basis,c->basis,9*sizeof(double));

    /* This version does not correct cell vector lengths and angles */
    cspg_op(c2,m2,s,NULL,CSPG_PRIM_NR+CSPG_NO_SORT,tolmin);
    if (debug>2){
      fprintf(stderr,"CSPG_PRIM_NR returns:\n");
      print_cell(c2,m2);
    }

    /* This version does, but may add an unwanted rotation and translation */
    cspg_op(c,m,s,NULL,CSPG_PRIM+CSPG_NO_SORT,tolmin);
    if (debug>2){
      fprintf(stderr,"CSPG_PRIM returns:\n");
      print_cell(c,m);
    }
    if (debug){
      fprintf(stderr,"In primitive cell, ");
      cspg_op(c,m,&stmp,NULL,CSPG_INT,tolmin);
      fprintf(stderr,"(%d symmetry operations)\n",stmp.n);
    }      

    /* Wonder about translation */
    for(i=0;i<3;i++)
      tr[i]=m->atoms[0].frac[i]-m2->atoms[0].frac[i];

    if (debug>2)
      fprintf(stderr,"Translation is %12.8f %12.8f %12.8f\n",
              tr[0],tr[1],tr[2]);

    if (!(op&CSPG_SNAP_TR)){
      if (debug>1) fprintf(stderr,"Discarding translation\n");
      if (vmod2(tr)>0.0){
	for(i=0;i<m->n;i++)
	  for(j=0;j<3;j++)
	    m->atoms[i].frac[j]-=tr[j];

	reduce_cell(m->atoms,m->n,c->basis);
	real2rec(c);
	addabs(m->atoms,m->n,c->basis);
      }
    }
    

    /* Wonder if a axis has been rotated */
    /* Could try evaluating theta as arccos(a.b), but this has accuracy
     *  issues if a.b is close to 1. So try 2 arcsin |0.5(a-b)|,
     *  for normalised a and b 
     */
    rotation=0;

    mod1=sqrt(vmod2(c->basis[0]));
    mod2=sqrt(vmod2(c2->basis[0]));
    if (!aeq(mod1,mod2)) error_exit("Unexpected length differences in a");
    if (debug>1) fprintf(stderr,"length difference in a: %e\n",mod1-mod2);
    for(i=0;i<3;i++)
      vtmp[i]=0.5*(c->basis[0][i]/mod1-c2->basis[0][i]/mod2);
    dot=sqrt(vmod2(vtmp));
    
    if (!aeq(dot,0.0)){ /* Need to rotate a */
      rotation=1;
      theta=2*asin(dot);
      vcross(c->basis[0],c2->basis[0],axis);

      /* It is possible for axis to be the zero vector, if theta=pi.
       * In this case, need to find a different axis perpendicular to a.
       * The cross product of a with the other basis axis which needs
       * the greatest rotation is a good answer.
       */

      if (aeq(theta,M_PI)&&(vmod2(axis)<1e-16)){
        mod1=0;
        mod2=0;
        for(i=0;i<3;i++){
          mod1+=c->basis[1][i]*c2->basis[1][i];
          mod2+=c->basis[2][i]*c2->basis[2][i];
        }
        k=1;
        if (mod2<mod1) k=2;
        vcross(c->basis[0],c->basis[k],axis);
      }
      
#if 0
      /* Will be rhs by construction */
      for(i=0;i<3;i++){
	rhs[0][i]=c->basis[0][i];
	rhs[1][i]=c2->basis[0][i];
	rhs[2][i]=axis[i];
      }
      if (is_rhs(rhs)) fprintf(stderr,"rhs\n");
      else fprintf(stderr,"lhs\n");
#endif
      
      if (debug>1) fprintf(stderr,"Rotating a by %lf\n",theta*180/M_PI);
      for(i=0;i<3;i++)
	rotate(c->basis[i],axis,theta);
      /* Did this work? */
      dot=0;
      for(i=0;i<3;i++) dot+=c2->basis[0][i]*c->basis[0][i];
      dot/=(mod1*mod2);
      if ((dot-1)<1e-12) dot=1;
      theta=acos(dot);
      if (debug>1) fprintf(stderr,"Required rotation now %lf\n",theta*180/M_PI);
      if (debug>2){
        fprintf(stderr,"New cell vectors of:\n");
	print_cell(c,NULL);
      }
    }

    /* a axis now okay, what about b? */

    mod1=sqrt(vmod2(c->basis[1]));
    mod2=sqrt(vmod2(c2->basis[1]));
    if (!aeq(mod1,mod2)) error_exit("Unexpected length differences in b");
    if (debug>1) fprintf(stderr,"length difference in b: %e\n",mod1-mod2);
    for(i=0;i<3;i++)
      vtmp[i]=0.5*(c->basis[1][i]/mod1-c2->basis[1][i]/mod2);
    dot=sqrt(vmod2(vtmp));

    if (!aeq(dot,0.0)){/* Need to rotate b */
      rotation=1;
      /* a should be preserved: make it the axis */
      for(i=0;i<3;i++) axis[i]=c->basis[0][i];
      mod1=sqrt(vmod2(axis));
      for(i=0;i<3;i++) axis[i]/=mod1;

      dot=0;
      for(i=0;i<3;i++) dot+=c->basis[1][i]*axis[i];
      for(i=0;i<3;i++) vtmp1[i]=c->basis[1][i]-dot*axis[i];
      mod1=sqrt(vmod2(vtmp1));

      dot=0;
      for(i=0;i<3;i++) dot+=c2->basis[1][i]*axis[i];
      for(i=0;i<3;i++) vtmp2[i]=c2->basis[1][i]-dot*axis[i];
      mod2=sqrt(vmod2(vtmp2));

      for(i=0;i<3;i++)
        vtmp[i]=0.5*(vtmp1[i]/mod1-vtmp2[i]/mod2);
      dot=sqrt(vmod2(vtmp));

      for(i=0;i<3;i++){
	rhs[0][i]=c->basis[1][i];
	rhs[1][i]=c2->basis[1][i];
	rhs[2][i]=axis[i];
      }
      if (is_rhs(rhs)) theta=2*asin(dot);
      else theta=-2*asin(dot);
      if (debug>1) fprintf(stderr,"Rotating b by %lf\n",theta*180/M_PI);

      for(i=0;i<3;i++)
	rotate(c->basis[i],axis,theta);

      /* Did this work? */
      dot=0;
      for(i=0;i<3;i++) dot+=c2->basis[1][i]*c->basis[1][i];
      mod1=sqrt(vmod2(c->basis[1]));
      mod2=sqrt(vmod2(c2->basis[1]));
      dot/=(mod1*mod2);
      if ((dot-1)<1e-12) dot=1;
      theta=acos(dot);
      if (debug>1) fprintf(stderr,"Required rotation now %lf\n",theta*180/M_PI);
      
    }    

    if (rotation){
      if (debug>1){
	fprintf(stderr,"*****New basis:\n");
	print_cell(c,NULL);
      }
      real2rec(c);
      addabs(m->atoms,m->n,c->basis);
    }
    
    /* Work out old basis in terms of new */

    for(i=0;i<3;i++){
      for(j=0;j<3;j++){
	old_frac[i][j]=0;
	for(k=0;k<3;k++)
	  old_frac[i][j]+=old_basis[i][k]*c->recip[j][k];
      }
    }

    /* These should all be integer... */

    for(i=0;i<3;i++){
      for(j=0;j<3;j++){
	if (!aeq(old_frac[i][j],floor(old_frac[i][j]+0.5))){
	  fprintf(stderr,"Unexpected transformation from primitive cell\n");
	  for(k=0;k<3;k++)
	    fprintf(stderr,"%6lf %6lf %6lf\n",old_frac[k][0],
		    old_frac[k][1],old_frac[k][2]);
	  exit(1);
	}
      }
    }

    /* Work out snapped old basis */

    for(i=0;i<3;i++){
      for(j=0;j<3;j++){
	old_basis[i][j]=0;
	for(k=0;k<3;k++)
	  old_basis[i][j]+=floor(old_frac[i][k]+0.5)*c->basis[k][j];
      }
    }

    if (debug>2){
      fprintf(stderr,"Snapped old basis:\n");
      for(i=0;i<3;i++)
	fprintf(stderr,"% 15.11f % 15.11f % 15.11f\n",old_basis[i][0],
		old_basis[i][1],old_basis[i][2]);
    }

    super(c,m,old_basis,NULL,s,NULL,0);

    free(c2->basis);
    free(c2);
    free(m2->atoms);
    free(m2);
    
    return 0;
    
  } /* end if (op&CSPG_SNAP) */

  spg=NULL;

  /* spg_refine_cell can return four times as many atoms as is passed */
  if ((op&CSPG_REF)||(op&CSPG_STD)||(op&CSPG_STD_IDEAL)){
    spg_pos=malloc(4*3*m->n*sizeof(double));
    spg_type=malloc(4*m->n*sizeof(int));
  }
  else{
    spg_pos=malloc(3*m->n*sizeof(double));
    spg_type=malloc(m->n*sizeof(int));
  }
  auid=malloc(m->n*sizeof(struct uniq));
  n_uniq=0;

  if ((!spg_pos)||(!spg_type)||(!auid)) {
    fprintf(stderr,"Malloc error in cspg_op\n");
    exit(1);
  }

  natoms2=m->n;
  spg_symprec=tol;

  n_uniq=0;
  for(i=0;i<m->n;i++){
    for(j=0;j<3;j++)
      spg_pos[i][j]=m->atoms[i].frac[j];
    hit=0;
    for(j=0;j<n_uniq;j++){
      /* Atoms the same if their atomic numbers are equal,
	 their spins are equal,
	 their label pointers are equal (possibly null), or
	 their label pointers are not null but point to strings which
	 are equal to strcmp
      */
      if((auid[j].atno==m->atoms[i].atno)&&
	 aeq(auid[j].spin,m->atoms[i].spin)&&
	 ((auid[j].label==m->atoms[i].label)||
	  ((auid[j].label)&&(m->atoms[i].label)&&
	   (!strcmp(auid[j].label,m->atoms[i].label))))){
	spg_type[i]=j;
	hit=1;
	break;
      }
    }
    if (hit==0){
      auid[n_uniq].atno=m->atoms[i].atno;
      auid[n_uniq].spin=m->atoms[i].spin;
      auid[n_uniq].label=m->atoms[i].label;
      spg_type[i]=n_uniq;
      n_uniq++;
    }
  }
  if(debug>1) fprintf(stderr,"%d unique atom types in c2x2spg\n",n_uniq);


  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      spg_latt[i][j]=c->basis[j][i];

  if ((tolmin!=tol)&&(op&(CSPG_INT+CSPG_SYM+CSPG_PNT+CSPG_SCH))){
    spg_scan(spg_latt,spg_pos,spg_type,m->n,op,tolmin,tol,0,0,1);
  }
    
  if (op&CSPG_REF){
    natoms2=spg_refine_cell(spg_latt,spg_pos,spg_type,m->n,spg_symprec);
    if (natoms2==0){
      fprintf(stderr,"spglib finds no refinement\n");
      natoms2=m->n;
      if (op==CSPG_REF) return 0;
    }
    else{
      m->n=natoms2;
      s->n=0;
      if (c->primitive){
	free_cell(c->primitive);
	c->primitive=NULL;
      }
    }
  }

  if ((op&CSPG_STD)||(op&CSPG_STD_IDEAL)){
    natoms2=spg_standardize_cell(spg_latt,spg_pos,spg_type,m->n,
				 0,(op&CSPG_STD_IDEAL)?0:1,spg_symprec);
    if (natoms2==0){
      fprintf(stderr,"spglib finds no standardisation\n");
      natoms2=m->n;
      if ((op==CSPG_STD)||(op==CSPG_STD_IDEAL)) return 0;
    }
    else{
      m->n=natoms2;
      if (op&CSPG_STD_IDEAL) {
	s->n=0;
	if (c->primitive){
	  free_cell(c->primitive);
	  c->primitive=NULL;
	}
      }
    }
  }
    
  if (op&CSPG_PRIM){
    natoms2=spg_find_primitive(spg_latt,spg_pos,spg_type,m->n,spg_symprec);
    if (natoms2==0){
      fprintf(stderr,"spglib finds no primitive cell\n");
      natoms2=m->n;
      if (op==CSPG_PRIM) return 0;
    }
    else{
      m->n=natoms2;
      s->n=0;
      if (c->primitive){
	free_cell(c->primitive);
	c->primitive=NULL;
      }
    }
  }

  /* This one will not rotate basis vectors in Cartesian space */
  if (op&CSPG_PRIM_NR){
    natoms2=spg_standardize_cell(spg_latt,spg_pos,spg_type,m->n,
				 1,1,spg_symprec);
    if (natoms2==0)
      error_exit("spglib finds no primitive cell\n");
    else{
      m->n=natoms2;
    }
  }
      
  if (op&(CSPG_SYM|CSPG_PNT|CSPG_INT|CSPG_SCH)){
    spg=spg_get_dataset(spg_latt,spg_pos,spg_type,m->n,spg_symprec);
    if ((tolmin==tol)||(debug))
      spg_sym_print(spg,spg_symprec,op);
    s->n=spg->n_operations;
  }
  

  
  /* Convert back to c2x-style variables */

  /* If we were asked for info, not a conversion, simply return */
  if (!(op&(~(CSPG_PNT|CSPG_INT|CSPG_SCH)))){ 
    free(auid);
    free(spg_type);
    free(spg_pos);
    if (spg) spg_free_dataset(spg);
    return 0;
  }

  if (debug) fprintf(stderr,"SPGlib returns %d atoms\n",m->n);

  
  /* First k-points */

  if (kp){
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	new_basis[j][i]=spg_latt[i][j];
    super(c,m,new_basis,kp,NULL,NULL,4);
  }

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[j][i]=spg_latt[i][j];
  real2rec(c);

  if (op==CSPG_PRIM_LATT){
    free(auid);
    free(spg_type);
    free(spg_pos);
    return 0;
  }
  
  if (debug&&(op&CSPG_PRIM_NR))
    print_old_in_new(old_basis,c->basis);

  /* Should we retain old sym-ops? */
  
  if ((op&(CSPG_PRIM_NR|CSPG_STD))&&(s->n)){
    sym_basis(s,c);
  }
  
  free(m->atoms);
  m->atoms=malloc(m->n*sizeof(struct atom));
  for(i=0;i<m->n;i++){
    for(j=0;j<3;j++)
      m->atoms[i].frac[j]=spg_pos[i][j];
    m->atoms[i].atno=auid[spg_type[i]].atno;
    m->atoms[i].spin=auid[spg_type[i]].spin;
    m->atoms[i].chg=0;
    m->atoms[i].label=auid[spg_type[i]].label;
  }
  addabs(m->atoms,m->n,c->basis);

  /* SPG is capable of returning fractional co-ordinates of 1,
   * whereas we like zeros */

  reduce_cell_tol(m->atoms,m->n,c->basis,1e-9);
  if (!(op&CSPG_NO_SORT)) sort_atoms(m,1);

  /* SPG has sym-ops in fractional co-ords, Castep in cartesians */
  if ((op&(CSPG_SYM|CSPG_PNT))&&(s->n)) {
    s->ops=sym_frac2abs(spg->rotations,spg->translations,c,s->n);
    if (op&CSPG_LST)
      for(i=0;i<s->n;i++) ident_sym(&s->ops[i],c,m,stderr);
  }

  free(auid);
  free(spg_type);
  free(spg_pos);
  if (spg) spg_free_dataset(spg);
  return 0;

}

  struct sym_op *sym_frac2abs(int spg_rot[][3][3],double spg_tr[][3],
			      struct unit_cell *c,int nsym){
    int i,j,k,kk;
    struct sym_op *s;
    double mat[3][3];

    s=malloc(nsym*sizeof(struct sym_op));
    if(!s) error_exit("Malloc error in sym_frac2abs");


    for(i=0;i<nsym;i++){
      if (debug>2){
	fprintf(stderr,"Sym op on entry:\n");
	for(j=0;j<3;j++)
	  fprintf(stderr,"%2d %2d %2d\n",spg_rot[i][j][0],
		  spg_rot[i][j][1],spg_rot[i][j][2]);
	fprintf(stderr,"%13.9f %13.9f %13.9f\n",spg_tr[i][0],spg_tr[i][1],
		spg_tr[i][2]);
      }
      /* Translations are "easy" */
      /* Sometimes SPGlib returns a translation vector component of 1 */
      for(j=0;j<3;j++)
	if (aeq(fabs(spg_tr[i][j]),1.0)) spg_tr[i][j]=0;
      if ((spg_tr[i][0]!=0)||(spg_tr[i][1]!=0)||(spg_tr[i][2]!=0)){
	s[i].tr=malloc(3*sizeof(double));
	if(!s[i].tr) error_exit("Malloc error for tr in sym_frac2abs");
	for(j=0;j<3;j++) s[i].tr[j]=spg_tr[i][0]*c->basis[0][j]+
			   spg_tr[i][1]*c->basis[1][j]+
			   spg_tr[i][2]*c->basis[2][j];
      }
      else s[i].tr=NULL;
      /* Matrices are harder. Permute indices until answer looks right? */
      for(j=0;j<3;j++)
	for(k=0;k<3;k++)
	  mat[j][k]=s[i].mat[j][k]=0.0;
      for(j=0;j<3;j++)
	for(k=0;k<3;k++)
	  for(kk=0;kk<3;kk++)
	    mat[j][k]+=spg_rot[i][j][kk]*c->recip[kk][k];
      for(j=0;j<3;j++)
	for(k=0;k<3;k++)
	  for(kk=0;kk<3;kk++)
	    s[i].mat[j][k]+=c->basis[kk][j]*mat[kk][k];
      /* We don't like matrix entries of "-0" */
      for(j=0;j<3;j++)
	for(k=0;k<3;k++)
	  if (fabs(s[i].mat[j][k])<1e-10) s[i].mat[j][k]=0;
   
      if (debug>2){
	fprintf(stderr,"Sym op on exit:\n");
	for(j=0;j<3;j++)
	  fprintf(stderr,"%13.9f %13.9f %13.9f\n",s[i].mat[j][0],
		  s[i].mat[j][1],s[i].mat[j][2]);
	if (s[i].tr)
	  fprintf(stderr,"%13.9f %13.9f %13.9f\n",s[i].tr[0],
		  s[i].tr[1],s[i].tr[2]);
	else
	  fprintf(stderr,"0.0 0.0 0.0\n");
      }
    }

    return s;
  }

  void cspg_hall2sym(int hall, struct unit_cell *c, struct symmetry *s){
    int rotations[192][3][3];
    double translations[192][3];

    s->n=spg_get_symmetry_from_database(rotations,translations,hall);
    s->ops=sym_frac2abs(rotations,translations,c,s->n);
  }

#else
int cspg_op(struct unit_cell *c, struct contents *m, struct symmetry *s,
            struct kpts *kp, int op, double tolmin){
  fprintf(stderr,"Error: spglib functionality not available"
	  " in this version of c2x\n");
  exit(1);
}

void cspg_hall2sym(int hall, struct unit_cell *c, struct symmetry *s){
  if (hall==1) return;
  fprintf(stderr,"Error: spglib functionality not available"
	  " so cannot expand spacegroup\n");
  exit(1);
}

#endif

/* Does space group have two possible origins? */
int spgr_is_double(int spgr){
  int doubles[]={48,50,59,70,85,86,88,125,126,129,130,133,
    134,137,138,141,142,201,203,222,224,227,228,0};
  int i,hit;
  
  hit=0;
  i=0;
  while(doubles[i]){
    if (doubles[i]==spgr){
      hit=1;
      break;
    }
    i++;
  }
  return hit;
}


  void rotate(double v[3], double axis[3], double theta){
    int i,j;
    double basis[3][3],frac[3];
    double tmp;

    if (debug>1)
      fprintf(stderr,"Rotating (%lf,%lf,%lf) by %lf about (%lf,%lf,%lf)\n",
	      v[0],v[1],v[2],theta,axis[0],axis[1],axis[2]);
  

    /* First basis vector is rotation axis */
  
    tmp=sqrt(vmod2(axis));
    for(i=0;i<3;i++) basis[0][i]=axis[i]/tmp;

    /* Second basis vector is residual from vector to be rotated */

    tmp=0;
    for(i=0;i<3;i++) tmp+=v[i]*basis[0][i];

    for(i=0;i<3;i++) basis[1][i]=v[i]-tmp*basis[0][i];

    tmp=sqrt(vmod2(basis[1]));

    if (tmp<1e-10) return; /* v and axis are parallel -- nothing to do */
    for(i=0;i<3;i++) basis[1][i]/=tmp;

    /* Third basis vector is cross product of other two */
  
    vcross(basis[0],basis[1],basis[2]);
    tmp=sqrt(vmod2(basis[2]));  /* Unnecessary -- should be 1 anyway */
    for(i=0;i<3;i++) basis[2][i]/=tmp;

    /* New basis is orthogonal, so easy to convert v into it */

    for(i=0;i<3;i++){
      frac[i]=0;
      for(j=0;j<3;j++)
	frac[i]+=v[j]*basis[i][j];
    }

    /* Now rotate -- frac[2] will be zero at this point */
    /* and frac[0] will be unchanged */
  
    frac[2]=frac[1]*sin(theta);
    frac[1]=frac[1]*cos(theta);

    /* And convert back to v */

    for(i=0;i<3;i++){
      v[i]=0;
      for(j=0;j<3;j++)
	v[i]+=frac[j]*basis[j][i];
    }
  
  }

#ifdef SPGLIB
  static void spg_scan(double spg_latt[3][3], double spg_pos[][3],
		       int *spg_type, int natoms, int op, double tolmin,
		       double tolmax, int hallmin, int hallmax, int first){
    int i,j,nsym,nhall;
    double tmp_latt[3][3];
    double (*tmp_pos)[3],spg_symprec;
    int *tmp_type;
    SpglibDataset *spg;

    if (debug>2) fprintf(stderr,"spg_scan called, tolmin=%g tolmax=%g "
			 "hallmin=%d hallmax=%d\n",tolmin,tolmax,
                         hallmin,hallmax);
    
    tmp_pos=malloc(3*natoms*sizeof(double));
    tmp_type=malloc(natoms*sizeof(int));
    if ((!tmp_pos)||(!tmp_type)) error_exit("Malloc error in spg_scan");
    spg=NULL;

    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	tmp_latt[i][j]=spg_latt[i][j];
    memcpy(tmp_pos,spg_pos,3*natoms*sizeof(double));
    memcpy(tmp_type,spg_type,natoms*sizeof(int));

    if (first){
      spg_symprec=tolmin;
      spg=spg_get_dataset(tmp_latt,tmp_pos,tmp_type,natoms,spg_symprec);
      spg_sym_print(spg,spg_symprec,op);
      if (spg) {
	hallmin=spg->hall_number;
	spg_free_dataset(spg);
      }
      else hallmin=0;

      for(i=0;i<3;i++)
	for(j=0;j<3;j++)
	  tmp_latt[i][j]=spg_latt[i][j];
      memcpy(tmp_pos,spg_pos,3*natoms*sizeof(double));
      memcpy(tmp_type,spg_type,natoms*sizeof(int));
      spg_symprec=tolmax;
      spg=spg_get_dataset(tmp_latt,tmp_pos,tmp_type,natoms,spg_symprec);
      if (spg)
	hallmax=spg->hall_number;
      else
	hallmax=0;

      if (hallmin==hallmax){
	free(tmp_pos);
	free(tmp_type);
	if (spg) spg_free_dataset(spg);
	return;
      }
      else{
	if (spg) spg_free_dataset(spg);
      }
    }
    
    spg_symprec=sqrt(tolmin*tolmax);
  
    spg=spg_get_dataset(spg_latt,spg_pos,spg_type,natoms,spg_symprec);
    if (spg){
      nsym=spg->n_operations;
      nhall=spg->hall_number;
    }
    else{
      nsym=0;
      nhall=0;
    }

    if (debug>2) fprintf(stderr,"spg_get_dataset returns nsym=%d,"
                         " int=%s and tol=%g\n",
			 nsym,spg->international_symbol,spg_symprec);
    
    if ((tolmax/spg_symprec>1.0003)&&(spg_symprec/tolmin>1.0003)){ /* Refine */
      if (nhall!=hallmin)
	spg_scan(spg_latt, spg_pos, spg_type, natoms,
		 op, tolmin, spg_symprec, hallmin, nhall, 0);
      if (nhall!=hallmax)
	spg_scan(spg_latt, spg_pos, spg_type, natoms,
		 op, spg_symprec, tolmax, nhall, hallmax, 0);
    }
    else{
      if (nhall!=hallmin)
	spg_sym_print(spg,spg_symprec,op);
      else if (nhall!=hallmax){
	if (spg) spg_free_dataset(spg);
	for(i=0;i<3;i++)
	  for(j=0;j<3;j++)
	    tmp_latt[i][j]=spg_latt[i][j];
	memcpy(tmp_pos,spg_pos,3*natoms*sizeof(double));
	memcpy(tmp_type,spg_type,natoms*sizeof(int));
	spg_symprec=tolmax;
	spg=spg_get_dataset(tmp_latt,tmp_pos,tmp_type,natoms,spg_symprec);
	spg_sym_print(spg,spg_symprec,op);
      }
    }

    free(tmp_type);
    free(tmp_pos);
    if (spg) spg_free_dataset(spg);
  
  }

  static void spg_sym_print(SpglibDataset *spg,double spg_symprec, int op){
  
    if (!spg){
      fprintf(stderr,"Tol=%-7.3g attempt to find symmetry failed\n",
	      spg_symprec);
      return;
    }
    
    if (op&CSPG_INT)
      fprintf(stderr,"Tol=%-7.3g International symmetry is %s\n",spg_symprec,
	      spg->international_symbol);
    if (op&CSPG_SYM)
      fprintf(stderr,"Tol=%-7.3g %d symmetry operations found.\n",spg_symprec,
	      spg->n_operations);
    if (op&CSPG_PNT)
      fprintf(stderr,"Tol=%-7.3g Point group is %s\n",spg_symprec,
	      spg->pointgroup_symbol);
    if (op&CSPG_SCH)
      fprintf(stderr,"Tol=%-7.3g Schoenflies symmetry is %s\n",spg_symprec,
	      spg_get_spacegroup_type(spg->hall_number).schoenflies);

  }
#endif


