/* Interface between c2x and spglib */

/* Copyright (c) 2014, 2016 MJ Rutter 
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
int cspq_op(struct unit_cell *c, struct contents *m, struct symmetry *s,
            int op, double tolmin){
  int i,j,natoms2;
  int dummy[3][3];
  int result,results,resulti,resultp;
  double stol;

  /* SPG variables */
  double spg_latt[3][3];
  double (*spg_pos)[3],spg_symprec;
  int *spg_type;
  char spg_sym[22];
  SpglibDataset *spg;

  struct uniq {int atno; double spin; char *label;} *auid;
  int hit,n_uniq;

  if (debug>1) fprintf(stderr,"Calling spglib with op=%x\n",op);
  else if(debug) fprintf(stderr,"Calling spglib\n");


  spg=NULL;
  result=0;

  /* spg_refine_cell can return four times as many atoms as is passed */
  if (op&CSPG_REF){
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

  results=resulti=resultp=-1;
  natoms2=m->n;
  stol=tolmin;

  while(stol<=tol){

    spg_symprec=stol;
    stol*=2;
    if ((stol>tol)&&(stol<1.99*tol)) stol=tol;

    n_uniq=0;
    for(i=0;i<m->n;i++){
      for(j=0;j<3;j++)
        spg_pos[i][j]=m->atoms[i].frac[j];
      hit=0;
      for(j=0;j<n_uniq;j++){
        /* Atoms the same if there atomic numbers are equal,
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
      }
    }

    if (op&CSPG_SCH){
      result=spg_get_schoenflies(spg_sym,spg_latt,spg_pos,spg_type,m->n,
                                 spg_symprec);

      if (!result)
        fprintf(stderr,"Attempt to find Schoenflies symmetry failed\n");
      else if (result!=results){
	if ((debug)||(2*tolmin<tol))
	  fprintf(stderr,"Tol=%g ",spg_symprec);
        fprintf(stderr,"Schoenflies symmetry is %s\n",spg_sym);
        results=result;
      }
    }
  
    if (op&(CSPG_SYM|CSPG_PNT|CSPG_INT)){
      spg=spg_get_dataset(spg_latt,spg_pos,spg_type,m->n,spg_symprec);
      if ((!spg->rotations)||(!spg->translations)){
	spg->n_operations=0;
      }
      if (!spg){
        fprintf(stderr,"Attempt to find symmetry failed\n");
        exit(1);
      }
      if (op&CSPG_INT){
        result=spg->spacegroup_number;
        if (!result)
          fprintf(stderr,"Attempt to find international symmetry failed\n");
        else if (result!=resulti){
          if ((debug)||(2*tolmin<tol))
	    fprintf(stderr,"Tol=%g ",spg_symprec);
          fprintf(stderr,"International symmetry is %s\n",
		  spg->international_symbol);
          resulti=result;
        }
      }
      if (op&(CSPG_SYM|CSPG_PNT)){
        result=spg->n_operations;
//        fprintf(stderr,"%d symmetry operations found.\n",spg->n_operations);
//        fprintf(stderr,"spg->rotations is %s null\n",(spg->rotations)?"not":"");
        if (result!=s->n){
          s->n=result;
          if ((debug)||(2*tolmin<tol))
	    fprintf(stderr,"Tol=%g ",spg_symprec);
	  fprintf(stderr,"%d symmetry operations found.\n",s->n);
	}
        if (op&CSPG_PNT){
          result=spg_get_pointgroup(spg_sym,dummy,spg->rotations,s->n);
          if (!result)
            fprintf(stderr,"Attempt to find point group failed\n");
          else if (result!=resultp){
            fprintf(stderr,"Tol=%g Point group is %s\n",spg_symprec,spg_sym);
            resultp=result;
          }
        }
      }
    }
  }
    
  m->n=natoms2;

  /* Convert back to c2x-style variables */
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[j][i]=spg_latt[i][j];
  real2rec(c);

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

  /* SPG has sym-ops in fractional co-ords, Castep in cartesians */
  if ((op&(CSPG_SYM|CSPG_PNT))&&(s->n)) {
    s->ops=sym_frac2abs(spg->rotations,spg->translations,c,s->n);
    if (op&CSPG_LST)
      for(i=0;i<s->n;i++) ident_sym(&s->ops[i],c,stderr);
  }

  free(auid);
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
int cspq_op(struct unit_cell *c, struct contents *m, struct symmetry *s,
            int op){
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
