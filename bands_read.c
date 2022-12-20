/* Read a CASTEP .bands file. Also a .cell file if found */

#include<stdio.h>
#include<stdlib.h> /* malloc */
#include<string.h>
#include<ctype.h>

#include "c2xsf.h"

#define LINE_SIZE 100
void cell_read(FILE* in, struct unit_cell *c, struct contents *m,
               struct kpts *kp, struct symmetry *s);

static void siesta_bands_read(FILE* infile, struct unit_cell *c,
                              struct contents *m, struct kpts *k,
                              struct symmetry *s,struct es *e);

void bands_read(FILE* infile, struct unit_cell *c, struct contents *m,
                 struct kpts *k, struct symmetry *s,struct es *e){
  int i,j,itmp,itmp2,first_line;
  int ik,nb,ns;
  double dtmp,dtmp2;
  char buffer[LINE_SIZE+1],*cptr,*newfile;
  FILE *f;
  struct kpts kdummy;
  
  /* Read header */
  first_line=1;
  while(1){
    if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
    if (first_line){
      if (sscanf(buffer,"%lf",&dtmp)==1){
        e->e_fermi=malloc(sizeof(double));
        if (!e->e_fermi) error_exit("malloc error");
        *(e->e_fermi)=dtmp;
        siesta_bands_read(infile, c, m, k, s, e);
        return;
      }
    }
    cptr=buffer;
    while(isspace(*cptr)) cptr++;
    if (*cptr==0) continue;
    if (*cptr=='#') continue;
    if (!strncasecmp(cptr,"Number of k-points",18)){
      if (sscanf(cptr+18,"%d",&itmp)!=1)
	error_exit("Error parsing number of kpoints");
      k->n=itmp;
      k->kpts=malloc(itmp*sizeof(struct atom));
      if (!k->kpts) error_exit("Malloc error for kpoints");
    }
    else if (!strncasecmp(cptr,"Number of spin components",25)){
      if (sscanf(cptr+25,"%d",&itmp)!=1)
	error_exit("Error parsing number of spins");
      if ((itmp<1)||(itmp>2))
	error_exit("Error, number of spins not one or two");
      e->nspins=itmp;
      e->nbspins=itmp;
    }
    else if (!strncasecmp(cptr,"Number of electrons",19)){
      if (e->nspins==1){
	if (sscanf(cptr+19,"%lf",&(e->nel))!=1)
	  error_exit("Error parsing number of electrons");
      }
      else if (e->nspins==2){
	if (sscanf(cptr+19,"%lf %lf",&dtmp,&dtmp2)!=2)
	  error_exit("Error parsing number of electrons");
	e->nel=dtmp+dtmp2;
	e->nup_minus_down=dtmp-dtmp2;
      }
      else
	error_exit("Found number of electrons, but nspins not set to 1 or 2");
    }
    else if (!strncasecmp(cptr,"Number of eigenvalues",21)){
      if (e->nspins==1){
	if (sscanf(cptr+21,"%d",&(e->nbands))!=1)
	  error_exit("Error parsing number of bands");
      }
      else if (e->nspins==2){
	if (sscanf(cptr+21,"%d %d",&itmp,&itmp2)!=2)
	  error_exit("Error parsing number of bands");
	if (itmp!=itmp2)
	  error_exit("Confused: number of bands differ for different spins");
	e->nbands=itmp;
      }
    }
    else if (!strncasecmp(cptr,"Fermi energy (in atomic units)",30)){
      e->e_fermi=malloc(sizeof(double));
      if (!e->e_fermi) error_exit("malloc error for a double!");
      if (sscanf(cptr+30,"%lf",e->e_fermi)!=1)
	error_exit("Error parsing Fermi energy");
      *(e->e_fermi)*=H_eV;
    }
    else if (!strncasecmp(cptr,"Fermi energies (in atomic units)",32)){
      e->e_fermi=malloc(sizeof(double));
      if (!e->e_fermi) error_exit("malloc error for a double!");
      if (e->nspins==2){
	if (sscanf(cptr+32,"%lf %lf",&dtmp,&dtmp2)!=2)
	  error_exit("Error parsing Fermi energies");
	*(e->e_fermi)=0.5*(dtmp+dtmp2)*H_eV;
      }
      else
	error_exit("Found Fermi energies, but nspins not set to 1 or 2");
    }
    else if (!strncasecmp(cptr,"Unit cell vectors",17)){
      c->basis=malloc(9*sizeof(double));
      if (!c->basis) error_exit("Malloc error for basis");
      for(i=0;i<3;i++){
	if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
	if (sscanf(buffer,"%lf %lf %lf",c->basis[i],
		   c->basis[i]+1,c->basis[i]+2)!=3)
	  error_exit("Error parsing basis");
	for(j=0;j<3;j++)
	  c->basis[i][j]*=BOHR;
      }
    }
    else if (!strncasecmp(cptr,"K-point",7)) break;
    else fprintf(stderr,"Ignoring: %s",cptr);
    first_line=0;
  }

  if (!c->basis)
    error_exit("Basis not found");
  real2rec(c);
  
  itmp=e->nbands*e->nspins*k->n;
  if (itmp==0) error_exit("One of nkpt, nbands and kspins is zero!");
  e->eval=malloc(itmp*sizeof(double));
  if (!e->eval) error_exit("Malloc error for eigenvalues");
  
  /* Read k-points */

  for(i=0;i<k->n;i++){
    if (i!=0)
      if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
    if (sscanf(buffer,"K-point %d%n",&ik,&itmp)!=1)
      error_exit("Error parsing kpoint line");
    if ((ik<1)||(ik>k->n))
      error_exit("Unexpected kpoint number");
    ik--;
    if (sscanf(buffer+itmp,"%lf %lf %lf %lf",
	       k->kpts[ik].frac,k->kpts[ik].frac+1,k->kpts[ik].frac+2,
	       &(k->kpts[ik].wt))!=4)
      error_exit("Error parsing kpoint line");
    for(ns=0;ns<e->nspins;ns++){
      if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
      if (sscanf(buffer,"Spin component %d",&itmp)!=1)
	error_exit("Error parsing spin component line");
      if (itmp!=ns+1)
	error_exit("Unexpected spin component number");
      for(nb=0;nb<e->nbands;nb++){
	if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
	if (sscanf(buffer,"%lf",&dtmp)!=1)
	  error_exit("Error parsing eigenvalue");
	e->eval[nb+e->nbands*ns+e->nbands*e->nspins*ik]=dtmp*H_eV;
      }
    }
  }

  /* See if we can find a .cell file */

  if (dict_get(m->dict,"in_file")){
    newfile=strrsubs(dict_get(m->dict,"in_file"),"bands","cell");
    if (newfile){
      f=fopen(newfile,"r");
      if (f){
	fprintf(stderr,"Additionally reading %s\n",newfile);
	kdummy.n=0;
	kdummy.kpts=NULL;
	kdummy.mp=NULL;
	kdummy.spacing=NULL;
	cell_read(f,c,m,&kdummy,s);
	fclose(f);
	if (kdummy.kpts) free(kdummy.kpts);
	if (kdummy.mp) free(kdummy.mp);
	if (kdummy.spacing) free(kdummy.spacing);
      }
    }
    free(newfile);
  }
    
}

/* Assume first line of file already read */
static void siesta_bands_read(FILE* infile, struct unit_cell *c,
                              struct contents *m, struct kpts *k,
                              struct symmetry *s,struct es *e){
  int nbands,nspins,nkpts;
  int ik,nb,is,off,i,j;
  char buffer[LINE_SIZE+1],*cptr,*newfile;
  double dtmp;
  FILE *f;
  struct kpts kdummy;
  
  if (debug) fprintf(stderr,"Reading Siesta bands file\n");
  
  if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");

  if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
  if (sscanf(buffer,"%d %d %d",&nbands,&nspins,&nkpts)!=3)
    error_exit("Error parsing third line");

  e->eval=malloc(nbands*nspins*nkpts*sizeof(double));
  e->nbands=nbands;
  e->nspins=nspins;
  if (!e->eval) error_exit("malloc error for evals");
  k->n=nkpts;
  k->kpts=malloc(nkpts*sizeof(struct atom));
  if (!e->eval) error_exit("malloc error for kpoints");
  
  for(ik=0;ik<nkpts;ik++){
    if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
    if (sscanf(buffer,"%lf %lf %lf %n",k->kpts[ik].abs,k->kpts[ik].abs+1,
               k->kpts[ik].abs+2,&off)!=3)
      error_exit("parse error for kpt");
    cptr=buffer+off;
    for(nb=0;nb<nbands;nb++){
      for(is=0;is<nspins;is++){
        if (sscanf(cptr,"%lf%n",&dtmp,&off)!=1){
          if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
          cptr=buffer;
          if (sscanf(cptr,"%lf%n",&dtmp,&off)!=1) error_exit("Parse error");
        }
        e->eval[nb+nbands*is+nbands*nspins*ik]=dtmp;
        cptr+=off;
      }
    }
  }

  if (debug) fprintf(stderr,"Read %d bands, %d spins, %d kpts\n",
                     nbands,nspins,nkpts);

  if (debug) print_bandwidths(e,k);
  
  /* See if we can find a .XV file */
  if(!c->basis){
    if (dict_get(m->dict,"in_file")){
      newfile=strrsubs(dict_get(m->dict,"in_file"),"bands","XV");
      if (newfile){
        f=fopen(newfile,"r");
        if (f){
          fprintf(stderr,"Additionally reading %s\n",newfile);
          xv_read(f,c,m);
          fclose(f);
        }
        else{
          free(newfile);
          newfile=strrsubs(dict_get(m->dict,"in_file"),"bands","fdf");
          if (newfile){
            f=fopen(newfile,"r");
            if (f){
              fprintf(stderr,"Additionally reading %s\n",newfile);
              kdummy.n=0;
              kdummy.kpts=NULL;
              kdummy.mp=NULL;
              kdummy.spacing=NULL;
              fdf_read(f,c,m,&kdummy,e);
              fclose(f);
              if (kdummy.kpts) free(kdummy.kpts);
              if (kdummy.mp) free(kdummy.mp);
              if (kdummy.spacing) free(kdummy.spacing);
            }
          }
        }
      }
      free(newfile);
    }
  } /* if (!c->basis) */

  if (!c->basis)
    error_exit("Unable to find basis");

  /* kpoint co-ordinates are absolute, and in reciprocal Bohr */

  for(i=0;i<k->n;i++)
    for(j=0;j<3;j++)
      k->kpts[i].abs[j]/=(2*M_PI*BOHR);
  
  addfrac(k->kpts,k->n,c->basis);

  if (tol==1e-4){
    fprintf(stderr,"Resetting tol to 5e-4 due to low precision in input\n");
    tol=5e-4;
  }
  
}

