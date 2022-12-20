/* Copyright (c) 2020 MJ Rutter 
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
#include<string.h>
#include<math.h>

#include "c2xsf.h"

void elk_write(FILE* outfile, struct unit_cell *c, struct contents *m,
              struct kpts *k, struct es *e){
  int i,j,*natomsp;
  double akpt_disp[3];
  char *fmt;
  
  /* Tasks */
  fprintf(outfile,"tasks\n");
  if (dict_get(m->dict,"Elk_tasks"))
    fprintf(outfile,"%s\n",(char*)dict_get(m->dict,"Elk_tasks"));
  else
    fprintf(outfile,"  0\n");
  
  /* Unit cell, in Bohr */

  if (flags&HIPREC)
    fmt="% 19.15f % 19.15f % 19.15f\n";
  else
    fmt="% 11.7f % 11.7f % 11.7f\n";
  
  fprintf(outfile,"\navec\n");
  for(i=0;i<3;i++)
    fprintf(outfile,fmt,
	    c->basis[i][0]/BOHR,c->basis[i][1]/BOHR,c->basis[i][2]/BOHR);

  /* Now atoms */
  
  if (dict_get(m->dict,"Elk_sppath"))
    fprintf(outfile,"\nsppath\n  %s\n",(char*)dict_get(m->dict,"Elk_sppath"));

  if (!m->spec) addspec(m);
  fprintf(outfile,"\natoms\n");
  fprintf(outfile,"  %d\n",m->nspec);
  natomsp=malloc(m->nspec*sizeof(int));
  if (!natomsp) error_exit("malloc error in elk_write");
  for(i=0;i<m->nspec;i++) natomsp[i]=0;
  for(i=0;i<m->n;i++){
    for(j=0;j<m->nspec;j++){
      if (m->atoms[i].atno==m->spec[j].atno){
        natomsp[j]++;
        break;
      }
    }
  }

  if (flags&HIPREC)
    fmt="% 19.15f % 19.15f % 19.15f  0.0 0.0 0.0\n";
  else
    fmt="% 11.7f % 11.7f % 11.7f  0.0 0.0 0.0\n";


  for(i=0;i<m->nspec;i++){
    if (dict_get(m->dict,"Elk_spfnames"))
      fprintf(outfile,"  %s\n",((char**)dict_get(m->dict,"Elk_spfnames"))[i]);
    else
      fprintf(outfile,"  '%s.in'\n",atno2sym(m->spec[i].atno));
    fprintf(outfile,"  %d\n",natomsp[i]);
    for(j=0;j<m->n;j++){
      if (m->atoms[j].atno==m->spec[i].atno){
	fprintf(outfile,fmt,m->atoms[j].frac[0],
                m->atoms[j].frac[1],m->atoms[j].frac[2]);
      }
    }
  }

  /* kpoints */

  if ((k->mp)&&(k->mp->grid[0]>0)){
    fprintf(outfile,"\nngridk\n");
    fprintf(outfile,"  %d %d %d\n",k->mp->grid[0],
            k->mp->grid[1],k->mp->grid[2]);

        /* Convert from Castep's convention to Elk's */
    for(i=0;i<3;i++){
      akpt_disp[i]=k->mp->disp[i]*k->mp->grid[i];
      if ((k->mp->grid[i]&1)==0) akpt_disp[i]+=0.5; /* is grid even? */
      akpt_disp[i]=fmod(akpt_disp[i],1.0);
    }
    if (flags&HIPREC)
      fmt="\nvlkoff\n  %.15g %.15g %.15g\n";
    else
      fmt="\nvkloff\n  %.11g %.11g %.11g\n";
    fprintf(outfile,fmt,akpt_disp[0],akpt_disp[1],akpt_disp[2]);
  }

  free(natomsp);
  
}
