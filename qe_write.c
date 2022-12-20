/* Write a Quantum Espresso format file */

/* Fails to deal with spin: QE expects spins to align with species,
 * not atoms...
 *
 * There is confusion over whether brackets are necessary for options
 * in cards, with
 *
 * K_POINTS gamma
 * K_POINTS {gamma}
 * K_POINTS (gamma)
 *
 * all seen. This code takes the no punctuation approach, which
 * pw.x 6.2 accepts.
 */

/* Copyright (c) 2018 MJ Rutter 
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

void qe_write(FILE* outfile, struct unit_cell *c, struct contents *m,
	      struct kpts *k, struct es *e){
  int nspec,i,j,*spatno,okay;
  char *fmt,*ptr;

  if (m->title) fprintf(outfile,"# %s\n\n",m->title);
  
  fprintf(outfile,"&CONTROL\n");
  if (m->title){
    /* remove single quotes */
    ptr=m->title;
    while (*ptr){
      if (*ptr=='\'') *ptr=' ';
      ptr++;
    }
    fprintf(outfile,"  title = '%s',\n",m->title);
  }
  fprintf(outfile,"  calculation = 'scf'\n");
  fprintf(outfile,"/\n");

  fprintf(outfile,"&SYSTEM\n");
  fprintf(outfile,"  ibrav = 0,\n");
  fprintf(outfile,"  nat   = %d,\n",m->n);

 /* Now we need to know the number of species.
     It must be fewer than the number of atoms...
  */

  spatno=malloc(m->n*sizeof(int));
  if (!spatno) error_exit("Malloc error in qe_write");
  nspec=0;
  
  for(i=0;i<m->n;i++){
    for(j=0;j<nspec;j++) if (m->atoms[i].atno==spatno[j]) break;
    if (j==nspec){  /* new species */
      spatno[j]=m->atoms[i].atno;
      nspec++;
    }
  }

  fprintf(outfile,"  ntyp = %d,\n",nspec);

  if (e->cut_off)
    fprintf(outfile,"  ecutwfc = %f\n",2*e->cut_off/H_eV);
  else
    fprintf(outfile,"  ecutwfc = 10\n");

  fprintf(outfile,"/\n");

  fprintf(outfile,"&ELECTRONS\n");
  fprintf(outfile,"/\n");

  if (flags&HIPREC)
    fmt="% 19.15f % 19.15f % 19.15f\n";
  else
    fmt="% 11.7f % 11.7f % 11.7f\n";

  fprintf(outfile,"\nCELL_PARAMETERS ");

  if (flags&AU){
    fprintf(outfile,"bohr\n");
    for(i=0;i<3;i++)
      fprintf(outfile,fmt,
                  c->basis[i][0]/BOHR,c->basis[i][1]/BOHR,c->basis[i][2]/BOHR);
  }
  else{
    fprintf(outfile,"angstrom\n");
    for(i=0;i<3;i++)
      fprintf(outfile,fmt,
                    c->basis[i][0],c->basis[i][1],c->basis[i][2]);
  }


  fprintf(outfile,"\nATOMIC_SPECIES\n");
  for(i=0;i<nspec;i++)
    fprintf(outfile,"  %s  1.0\n",atno2sym(spatno[i]));
  
  if (flags&HIPREC)
    fmt="%3s % .15f % .15f % .15f\n";
  else
    fmt="%3s % .9f % .9f % .9f\n"; 
  
  fprintf(outfile,"\nATOMIC_POSITIONS ");
  if (flags&AU){
    fprintf(outfile,"bohr\n");
    for(i=0;i<m->n;i++)
      fprintf(outfile,fmt,
	      atno2sym(m->atoms[i].atno),m->atoms[i].abs[0]/BOHR,
	      m->atoms[i].abs[1]/BOHR,m->atoms[i].abs[2]/BOHR);
  }
  else{
    fprintf(outfile,"angstrom\n");
    for(i=0;i<m->n;i++)
      fprintf(outfile,fmt,
	      atno2sym(m->atoms[i].atno),m->atoms[i].abs[0],
	      m->atoms[i].abs[1],m->atoms[i].abs[2]);
  }

  if ((k->n)||(k->mp)){
    fprintf(outfile,"\nK_POINTS ");
    if (k->mp){
      okay=1;
      for(i=0;i<3;i++)
        if ((k->mp->disp[i]!=0)&&(!aeq(k->mp->disp[i],0.5/k->mp->grid[i])))
          okay=0;
      if (okay){
        fprintf(outfile,"automatic\n");
        fprintf(outfile,"  %d %d %d",
	      k->mp->grid[0],k->mp->grid[1],k->mp->grid[2]);
        /* See comments in qe_read.c */
        for (i=0;i<3;i++)
          if (((k->mp->disp[i]==0)&&((k->mp->grid[i]&1)==1))||
              ((k->mp->disp[i]!=0)&&((k->mp->grid[i]&1)==0)))
            fprintf(outfile," 0");
          else
            fprintf(outfile," 1");
        fprintf(outfile,"\n");
      }
      else{ /* Need to list k-points explicitly */
        mp_gen(k,c);
      }
    }
    if ((k->n==1)&&(k->kpts[0].frac[0]==0)
	&&(k->kpts[0].frac[1]==0)&&(k->kpts[0].frac[2]==0))
      fprintf(outfile,"gamma\n");
    else if (k->n>0){
      fprintf(outfile,"crystal\n");
      fprintf(outfile,"  %d",k->n);
      if (flags&HIPREC){
	fmt="% 19.15f % 19.15f % 19.15f     %19.15f\n";
      }
      else{
	fmt="% 13.9f % 13.9f % 13.9f     %12.9f\n";
      }
      for(i=0;i<k->n;i++)
      fprintf(outfile,fmt,k->kpts[i].frac[0],
              k->kpts[i].frac[1],k->kpts[i].frac[2],k->kpts[i].wt);
    }
  }

}
