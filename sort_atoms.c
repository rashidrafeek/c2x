/* Conversion to arbitrary, specified supercells, with grid interpolation */

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

