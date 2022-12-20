/* Write a Quantum Espresso format file */


/* There is confusion over whether brackets are necessary for options
 * in cards, with
 *
 * K_POINTS gamma
 * K_POINTS {gamma}
 * K_POINTS (gamma)
 *
 * all seen. This code takes the no punctuation approach, which
 * pw.x 6.2 accepts.
 *
 * https://www.quantum-espresso.org/Doc/INPUT_PW.html
 */

/* Copyright (c) 2018, 2019 MJ Rutter 
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
#include<ctype.h>

#include "c2xsf.h"

void qe_write(FILE* outfile, struct unit_cell *c, struct contents *m,
	      struct kpts *k, struct es *e){
  int nspec,i,j,*spatno,okay,xtra;
  char *fmt,*ptr;
  double *spspin;
  int *spxtra,*atxtra;

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
  if (dict_get(m->dict,"QE_prefix"))
    fprintf(outfile,"  prefix = '%s',\n",(char*)dict_get(m->dict,"QE_prefix"));
  fprintf(outfile,"  calculation = 'scf',\n");
  fprintf(outfile,"  outdir = '.',\n");
  fprintf(outfile,"  pseudo_dir = '%s'\n",dict_get(m->dict,"QE_pseudo_dir")?
	  (char*)dict_get(m->dict,"QE_pseudo_dir"):".");
  fprintf(outfile,"/\n");

  fprintf(outfile,"&SYSTEM\n");
  fprintf(outfile,"  ibrav = 0,\n");
  fprintf(outfile,"  nat   = %d,\n",m->n);

 /* Now we need to know the number of species.
     It must be fewer than the number of atoms...
  */

  spatno=malloc(m->n*sizeof(int));
  spxtra=malloc(m->n*sizeof(int));
  atxtra=malloc(m->n*sizeof(int));
  spspin=malloc(m->n*sizeof(double));
  if ((!spatno)||(!spxtra)||(!atxtra)||(!spspin))
    error_exit("Malloc error in qe_write");
  nspec=0;
  
  for(i=0;i<m->n;i++){
    xtra=0;
    for(j=0;j<nspec;j++) {
      if (m->atoms[i].atno==spatno[j]) {
	if (m->atoms[i].spin==spspin[j]) break;
	else xtra++;
      }
    }
    if (j==nspec){  /* new species */
      spatno[j]=m->atoms[i].atno;
      spspin[j]=m->atoms[i].spin;
      spxtra[j]=xtra;
      nspec++;
    }
    atxtra[i]=spxtra[j];
  }

  fprintf(outfile,"  ntyp = %d,\n",nspec);

  okay=1;  /* Do we have spin? */
  for(i=0;i<nspec;i++)
    if (spspin[i]!=0) okay=0;
  if (!okay){
    fprintf(outfile,"  nspin = 2,\n");
    for(i=0;i<nspec;i++)
      fprintf(outfile,"  starting_magnetization(%d) = %lf,\n",i+1,spspin[i]);
  }

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
  if (dict_get(m->dict,"QE_atomic_species"))
    fprintf(outfile,"%s",(char*)dict_get(m->dict,"QE_atomic_species"));
  else{
    for(i=0;i<nspec;i++){
      if (spxtra[i]==0)
        fprintf(outfile,"  %-3s",atno2sym(spatno[i]));
      else
        fprintf(outfile,"  %2s%1d",atno2sym(spatno[i]),spxtra[i]);
      fprintf(outfile,"  1.0   %s.UPF\n",atno2sym(spatno[i]));
    }
  }
  
  if (flags&HIPREC)
    fmt=" % .15f % .15f % .15f\n";
  else
    fmt=" % .9f % .9f % .9f\n"; 
  
  fprintf(outfile,"\nATOMIC_POSITIONS ");
  if (flags&FRAC)
    fprintf(outfile,"crystal\n");
  else if (flags&AU)
    fprintf(outfile,"bohr\n");
  else
    fprintf(outfile,"angstrom\n");
  for(i=0;i<m->n;i++){
    okay=0;
    if (m->atoms[i].label&&(strlen(m->atoms[i].label)<=3)){
      j=strlen(atno2sym(m->atoms[i].atno));
      if ((!strncasecmp(atno2sym(m->atoms[i].atno),m->atoms[i].label,j))&&
	  (!isalpha(m->atoms[i].label[j]))) okay=1;
    }
    if (okay)
      fprintf(outfile," %3s",m->atoms[i].label);
    else{
      fprintf(outfile,"%3s",atno2sym(m->atoms[i].atno));
      if (atxtra[i]) fprintf(outfile,"%d",atxtra[i]);
      else fprintf(outfile," ");
    }
    if (flags&FRAC)
      fprintf(outfile,fmt,m->atoms[i].frac[0],
	      m->atoms[i].frac[1],m->atoms[i].frac[2]);
    else if (flags&AU)
      fprintf(outfile,fmt,m->atoms[i].abs[0]/BOHR,
	      m->atoms[i].abs[1]/BOHR,m->atoms[i].abs[2]/BOHR);
    else
      fprintf(outfile,fmt,m->atoms[i].abs[0],
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
      fprintf(outfile,"  %d\n",k->n);
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
