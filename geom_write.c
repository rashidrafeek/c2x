/* Write a CASTEP .geom file */

/* Units are always atomic units (Ha, Bohr) */

/* Copyright (c) 2019-2021 MJ Rutter 
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
#include<string.h>
#include<stdlib.h>

#include "c2xsf.h"

void geom_write(FILE* outfile, struct unit_cell *c, struct contents *m_in,
		struct es *e, struct time_series *ts){
  int ii,i,j,free_ts;
  char *sp20="                    ";
  int nspec,*natomsp,*spatno,*ninspec;
  struct contents *m;

  free_ts=0;
  if (!ts){
    ts=malloc(sizeof(struct time_series));
    if (!ts) error_exit("malloc error in geom_write");
    free_ts=1;
    init_tseries(ts);
  }

  if (ts->nsteps==0){
    ts->nsteps=ts->nc=ts->nm=ts->nen=1;
    ts->m=m_in;
    ts->cells=c;
    ts->energies=e->energy;
    if (!ts->energies) ts->nen=0;
  }

  if ((ts->nm<ts->nsteps)||(ts->nc<ts->nsteps)){
    error_exit("Too few steps for atoms or cell. Aborting.\n");
  }
  
  fprintf(outfile," BEGIN header\n");
  fprintf(outfile," produced by c2x");
  if (dict_get(m_in->dict,"in_file"))
    fprintf(outfile," from %s",(char*)dict_get(m_in->dict,"in_file"));
  fprintf(outfile,"\n END header\n");
  fprintf(outfile," \n");

  for(ii=0;ii<ts->nsteps;ii++){
    if (ii==(ts->nsteps-1))
      fprintf(outfile,"%s         %10d      "
              "%s           T   T   T   T            <-- c\n",
              sp20,ii,sp20);
    else
      fprintf(outfile,"%s         %10d      "
              "%s           F   F   F   F            <-- c\n",
              sp20,ii,sp20);

    if (ts->nen==ts->nsteps)
      fprintf(outfile,"%s %23.16e    %23.16e"
              "%s          <-- E\n",sp20,
              ts->energies[ii]/H_eV,
              (ts->enthalpies)?ts->enthalpies[ii]/H_eV:ts->energies[ii]/H_eV,
              sp20);
    /* Axes */
    for(j=0;j<3;j++){
      fprintf(outfile,"%s %23.16e    %23.16e    %23.16e   <-- h\n",sp20,
              ts->cells[ii].basis[j][0]/BOHR,
              ts->cells[ii].basis[j][1]/BOHR,
              ts->cells[ii].basis[j][2]/BOHR);
    }
    /* Stress should be here */
    /* Atoms */
    m=ts->m+ii;

    nspec=0;
    natomsp=malloc(m->n*sizeof(int));
    if (!natomsp) error_exit("Malloc error in geom_write");
    spatno=malloc(m->n*sizeof(int));
    if (!spatno) error_exit("Malloc error in geom_write");
    ninspec=malloc(m->n*sizeof(int));
    if (!ninspec) error_exit("Malloc error in geom_write");
    
    for(i=0;i<m->n;i++){
      for(j=0;j<nspec;j++) if (m->atoms[i].atno==spatno[j]) break;
      if (j==nspec){  /* new species */
        spatno[j]=m->atoms[i].atno;
        natomsp[j]=1;
        ninspec[i]=1;
        nspec++;
      }else{          /* existing species */
        natomsp[j]++;
        ninspec[i]=natomsp[j];
      }
    }
    
    /* Write out atoms, sorted by species */
    for(i=0;i<nspec;i++){
      for(j=0;j<m->n;j++){
        if (m->atoms[j].atno==spatno[i]){
          fprintf(outfile,"%3s       %8d   %23.16e    %23.16e    %23.16e"
                  "   <-- R\n",
                  atno2sym(m->atoms[j].atno),ninspec[j],
                  m->atoms[j].abs[0]/BOHR,
                  m->atoms[j].abs[1]/BOHR,
                  m->atoms[j].abs[2]/BOHR);
        }
      }
    }

    /* Write out forces, sorted by species */
    if (m->forces){
      for(i=0;i<nspec;i++){
        for(j=0;j<m->n;j++){
          if (m->atoms[j].atno==spatno[i]){
            fprintf(outfile,"%3s       %8d   %23.16e    %23.16e    %23.16e"
                    "   <-- F\n",
                    atno2sym(m->atoms[j].atno),ninspec[j],
                    m->atoms[j].force[0]*=(BOHR/H_eV),
                    m->atoms[j].force[1]*=(BOHR/H_eV),
                    m->atoms[j].force[2]*=(BOHR/H_eV));
          }
        }
      }
    }
    fprintf(outfile,"\n");
    free(ninspec);
    free(spatno);
    free(natomsp);
  }

  if (free_ts) free(ts);
}
