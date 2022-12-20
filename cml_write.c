/* Someone will want CML output, even if I cannot imagine why.
 * This might amuse him.
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

#include<stdio.h>
#include<stdlib.h> /* calloc */

#include "c2xsf.h"

extern int periodic_max_el;

void cml_write(FILE* outfile, struct unit_cell *c, struct contents *m){
  int i,*n_in_el;
  double abc[6];

  n_in_el=calloc(periodic_max_el+1,sizeof(int));
  if (!n_in_el) error_exit("Calloc error in cml_write");

  fprintf(outfile,"<?xml version=\"1.0\"?>\n"
         "<molecule xmlns=\"http://www.xml-cml.org/schema\">\n");

  cart2abc(c,m,abc,NULL,1);

  fprintf(outfile,
      "<crystal xmlns:cmldict=\"http://www.xml-cml.org/dict/cmlDict\">\n");
  fprintf(outfile,"  <scalar title=\"a\">%f</scalar>\n",abc[0]);
  fprintf(outfile,"  <scalar title=\"b\">%f</scalar>\n",abc[1]);
  fprintf(outfile,"  <scalar title=\"c\">%f</scalar>\n",abc[2]);
  fprintf(outfile,"  <scalar title=\"alpha\">%f</scalar>\n",abc[3]);
  fprintf(outfile,"  <scalar title=\"beta\">%f</scalar>\n",abc[4]);
  fprintf(outfile,"  <scalar title=\"gamma\">%f</scalar>\n",abc[5]);
  fprintf(outfile,"</crystal>\n");

  fprintf(outfile,"<atomArray>\n");
  for(i=0;i<m->n;i++)
      fprintf(outfile,"  <atom id=\"%s%d\" xFract=\"%f\" yFract=\"%f\""
              " zFract=\"%f\" elementType=\"%s\"/>\n",atno2sym(m->atoms[i].atno),
              ++n_in_el[min(m->atoms[i].atno,periodic_max_el)],
              m->atoms[i].frac[0],m->atoms[i].frac[1],m->atoms[i].frac[2],
              atno2sym(m->atoms[i].atno));
  fprintf(outfile,"</atomArray>\n");
  fprintf(outfile,"</molecule>\n");

  free(n_in_el);
}
