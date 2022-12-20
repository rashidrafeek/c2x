/* Units are always atomic units (Ha, Bohr) */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ctype.h>

#include "c2xsf.h"

#define LINE_SIZE 200

void geom_read(FILE* infile, struct unit_cell *cf, struct contents *mf,
		struct time_series *ts){
  char buffer[LINE_SIZE+1],species[5],*cptr;
  int i,j;
  int natoms,nforces,nsteps,rl,dummy,n;
  struct contents *m;
  struct unit_cell *c;
  
  nsteps=0;
  buffer[0]=0;
  rl=0;
  
  while (rl||fgets(buffer,LINE_SIZE,infile)){
    rl=0;
    if (strstr(buffer,"<-- c")){
      nsteps++;
      ts->cells=realloc(ts->cells,nsteps*sizeof(struct unit_cell));
      if (!ts->cells) error_exit("realloc error for cells");
      ts->m=realloc(ts->m,nsteps*sizeof(struct contents));
      if (!ts->m) error_exit("realloc error for motifs");
      c=ts->cells+(nsteps-1);
      m=ts->m+(nsteps-1);
      
      c->basis=malloc(9*sizeof(double));
      if (!c->basis) error_exit("realloc error for basis");
      i=0;
      m->atoms=NULL;
      m->n=natoms=nforces=0;
      m->title=NULL;
      m->dict=NULL;
      m->comment=NULL;
      n=0;
      
      while(fgets(buffer,LINE_SIZE,infile)){
	cptr=buffer;
	while(isspace(*cptr)) cptr++;
	if (*cptr==0) break;
	if (strstr(buffer,"<-- c")) {rl=1; break;}
	if (buffer[0]=='\n') break;
	if (strstr(buffer,"<-- E")){
	  ts->energies=realloc(ts->energies,nsteps*sizeof(double));
	  ts->enthalpies=realloc(ts->enthalpies,nsteps*sizeof(double));
	  if (sscanf(buffer,"%lf %lf",ts->energies+(nsteps-1),
		      ts->enthalpies+(nsteps-1))!=2)
	    error_exit("Error parsing energies");
	  ts->energies[nsteps-1]*=H_eV;
	  ts->enthalpies[nsteps-1]*=H_eV;
	  ts->nen++;
	}
	else if (strstr(buffer,"<-- h")){
	  if (sscanf(buffer,"%lf %lf %lf",c->basis[i],c->basis[i]+1,
		     c->basis[i]+2)!=3) error_exit("Error parsing cell axes");
	  for(j=0;j<3;j++)
	    c->basis[i][j]*=BOHR;
	  i++;
	  if (i==3){
	    real2rec(c);
	    ts->nc++;
	  }
	}
	else if (strstr(buffer,"<-- R")){
	  if (natoms==0) ts->nm++;
	  natoms++;
	  m->atoms=realloc(m->atoms,natoms*sizeof(struct atom));
	  init_atoms(m->atoms+(natoms-1),1);
	  
	  if (sscanf(buffer,"%4s %d %lf %lf %lf",species,&dummy,
		     m->atoms[natoms-1].abs,
		     m->atoms[natoms-1].abs+1,m->atoms[natoms-1].abs+2)!=5)
	    error_exit("Error parsing atomic co-ords");
	  m->atoms[natoms-1].atno=atsym2no(species);
	  if ((natoms==1)||(m->atoms[natoms-1].atno!=m->atoms[natoms-2].atno))
	    n=1;
	  else
	    n++;
	  if (dummy!=n) error_exit("Unexpected atom ordering in .geom file");
	  for(j=0;j<3;j++)
	    m->atoms[natoms-1].abs[j]*=BOHR;
	  addfrac(m->atoms+(natoms-1),1,c->recip);
	  m->n=natoms;
	}
	else if (strstr(buffer,"<-- F")){
	  if (sscanf(buffer,"%4s %d %lf %lf %lf",species,&dummy,
		     m->atoms[nforces].force,
		     m->atoms[nforces].force+1,m->atoms[nforces].force+2)!=5)
	    error_exit("Error parsing atomic forces");
	  if (atsym2no(species)!=m->atoms[nforces].atno)
	    error_exit("Inconsistent ordering between forces and positions");
	  if ((nforces==0)||(m->atoms[nforces].atno!=m->atoms[nforces-1].atno))
	    n=1;
	  else
	    n++;
	  if (dummy!=n) error_exit("Unexpected force ordering in .geom file");
	  for(j=0;j<3;j++)
	    m->atoms[nforces].force[j]*=H_eV/BOHR;
	  m->forces=1;
	  nforces++;
	}
      }
    }
  }

  ts->nsteps=nsteps;

  if (debug) fprintf(stderr,"Read %d steps\n",nsteps);
  if (debug>1) fprintf(stderr,"Read %d cell steps, %d atom steps\n",
		       ts->nc,ts->nm);
  
  if (ts->cells[nsteps-1].basis){
    memcpy(cf,ts->cells+(nsteps-1),sizeof(struct unit_cell));
    cf->basis=malloc(9*sizeof(double));
    if (!cf->basis) error_exit("malloc error for basis");
    memcpy(cf->basis,ts->cells[nsteps-1].basis,9*sizeof(double));
  }

  if (ts->m[nsteps-1].atoms){
    mf->n=ts->m[nsteps-1].n;
    mf->atoms=malloc(mf->n*sizeof(struct atom));
    if (!mf->atoms) error_exit("malloc erro for atoms");
    memcpy(mf->atoms,ts->m[nsteps-1].atoms,mf->n*sizeof(struct atom));
    mf->n=ts->m[nsteps-1].n;
    mf->forces=ts->m[nsteps-1].forces;
  }

}
