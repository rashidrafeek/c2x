
/* Copyright (c) 2015 MJ Rutter 
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

/* This ordering is Human-friendly */
int atom_sort(const void *a, const void *b){ /* This ghastly declaration */
  struct atom *a1 = (struct atom *) a;       /* stops compilers moaning */
  struct atom *a2 = (struct atom *) b;       /* about qsort's prototype */
  double stol;

  stol=tol;
  if (stol>0.01) stol=0.01;

  if (a1->atno < a2->atno) return(1);
  if (a1->atno > a2->atno) return(-1);
  if (a1->frac[2] < a2->frac[2]-stol) return(-1);
  if (a1->frac[2] > a2->frac[2]+stol) return(1);
  if (a1->frac[1] < a2->frac[1]-stol) return(-1);
  if (a1->frac[1] > a2->frac[1]+stol) return(1);
  if (a1->frac[0] < a2->frac[0]-stol) return(-1);
  if (a1->frac[0] > a2->frac[0]+stol) return(1);
  return(0); /* All co-ords equal! */
}

/* This ordering is more Castep-friendly (but only "more") */
int atom_sort2(const void *a, const void *b){ /* This ghastly declaration */
  struct atom *a1 = (struct atom *) a;       /* stops compilers moaning */
  struct atom *a2 = (struct atom *) b;       /* about qsort's prototype */
  double stol;

  stol=tol;
  if (stol>0.01) stol=0.01;

  if (a1->atno > a2->atno) return(1);
  if (a1->atno < a2->atno) return(-1);
  if (a1->frac[2] < a2->frac[2]-stol) return(-1);
  if (a1->frac[2] > a2->frac[2]+stol) return(1);
  if (a1->frac[1] < a2->frac[1]-stol) return(-1);
  if (a1->frac[1] > a2->frac[1]+stol) return(1);
  if (a1->frac[0] < a2->frac[0]-stol) return(-1);
  if (a1->frac[0] > a2->frac[0]+stol) return(1);
  return(0); /* All co-ords equal! */
}

void sort_atoms(struct contents *mtf, int sort_style){
  int na;

  na=mtf->n;

  if (sort_style==1)
    qsort(mtf->atoms,(size_t)na,sizeof(struct atom),atom_sort);
  else if (sort_style==2)
    qsort(mtf->atoms,(size_t)na,sizeof(struct atom),atom_sort2);
  else
    fprintf(stderr,"Ignoring undefined sort type %d\n",sort_style);

}

int symop_cmp(const void *a1, const void *b1){
  struct sym_op *a = (struct sym_op *) a1;
  struct sym_op *b = (struct sym_op *) b1;

  double mult_a,mult_b;

  mult_a=ident_sym(a,NULL,NULL,NULL);
  if (mult_a<0) mult_a=fabs(mult_a)+0.5;

  mult_b=ident_sym(b,NULL,NULL,NULL);
  if (mult_b<0) mult_b=fabs(mult_b)+0.5;
  
  if (mult_a>mult_b) return(1);
  if (mult_a<mult_b) return(-1);
  if ((!a->tr)&&(!b->tr)) return(0);
  if ((a->tr)&&(!b->tr)) return(1);
  if ((b->tr)&&(!a->tr)) return(-1);
  if (vmod2(a->tr)>vmod2(b->tr)) return(1);
  if (vmod2(b->tr)>vmod2(a->tr)) return(-1);
  if (a->tr[0]>b->tr[0]) return(1);
  if (b->tr[0]>a->tr[0]) return(-1);
  if (a->tr[1]>b->tr[1]) return(1);
  if (b->tr[1]>a->tr[1]) return(-1);
  if (a->tr[2]>b->tr[2]) return(1);
  if (b->tr[2]>a->tr[2]) return(-1);
  return 0;
}
  
void sort_symops(struct symmetry *s, int sort_style){

  qsort(s->ops,s->n,sizeof(struct sym_op),symop_cmp);
  
}
