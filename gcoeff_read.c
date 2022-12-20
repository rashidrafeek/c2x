/* Read a gcoeff.txt wavefunction file
 *
 * (c) 2020 MJ Rutter
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "c2xsf.h"

#define LINE_SIZE 100


void gcoeff_read(FILE *infile, struct unit_cell *c,
                 struct contents *m, struct kpts *k, struct grid *g,
                 struct es *e, int *i_grid){
  int i,j,nspins,nkpt,nbands,nplwv,junk,itmp,off;
  int ns,ikpt,nb,fft[3];
  long seek_save;
  char buffer[LINE_SIZE+1];
  double *psi,djunk;
  int *pwgrid;

  dict_add(m->dict,"band_read_order",NULL); /* Delete any old entry */
  dict_strcat(m->dict,"band_read_order","skb"); /* Malloc for new */
  
  psi=NULL;
  pwgrid=NULL;
  fft[0]=fft[1]=fft[2]=0;
  
  if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
  if (sscanf(buffer,"%d",&nspins)!=1)
    error_exit("Unable to read number of spins");
  e->nspins=e->nbspins=nspins;

  if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
  if (sscanf(buffer,"%d",&nkpt)!=1)
    error_exit("Unable to read number of kpoints");
  k->n=nkpt;
  k->kpts=malloc(nkpt*sizeof(struct atom));
  if (!k->kpts) error_exit("Malloc error for kpoints");

  if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
  i=sscanf(buffer,"%d %lf %lf",&nbands,&e->cut_off,&djunk);
  if (i==0) error_exit("Unable to read number of bands");
  if (i==3){
    e->e_fermi=malloc(sizeof(double));
    if (!e->e_fermi) error_exit("malloc error");
    *e->e_fermi=djunk;
  }
  e->nbands=nbands;
  e->occ=malloc(nspins*nkpt*nbands*sizeof(double));
  if (!e->occ) error_exit("Malloc error for occupancies");
  e->eval=malloc(nspins*nkpt*nbands*sizeof(double));
  if (!e->eval) error_exit("Malloc error for evals");

  c->basis=malloc(9*sizeof(double));
  if (!c->basis) error_exit("Malloc error for basis");
  for(i=0;i<3;i++){
    if (!fgets(buffer,LINE_SIZE,infile))
      error_exit("Unexpected EOF reading basis");
    if (sscanf(buffer,"%lf %lf %lf",c->basis[i],c->basis[i]+1,
	       c->basis[i]+2)!=3) error_exit("Parse error for basis");
  }
  real2rec(c);
  /* Skip offered recip basis vectors */
  for(i=0;i<3;i++){
    if (!fgets(buffer,LINE_SIZE,infile))
      error_exit("Unexpected EOF reading recip basis");
  }

  for(ns=0;ns<nspins;ns++){
    for(ikpt=0;ikpt<nkpt;ikpt++){
      if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
      if (sscanf(buffer,"%lf %lf %lf",k->kpts[ikpt].frac,
		 k->kpts[ikpt].frac+1,k->kpts[ikpt].frac+2)!=3)
	error_exit("Error parsing kpoint");
      k->kpts[ikpt].wt=0;
      off=(ikpt*e->nbspins+ns)*e->nbands;
      for(nb=0;nb<nbands;nb++){
	if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
	if (sscanf(buffer,"%d %d",&junk,&itmp)!=2)
	  error_exit("Parse error");
	if (nb==0) nplwv=itmp;
	else if (nplwv!=itmp) error_exit("nplwv changes within a kpt");
	if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
	if (sscanf(buffer,"( %lf , %lf ) %lf",e->eval+off+nb,&djunk,
		   e->occ+off+nb)!=3)
	  error_exit("Parse error");
	if (nb==0){
	  if (psi) free(psi);
	  if (pwgrid) free(pwgrid);
	  psi=malloc(2*nplwv*sizeof(double));
	  if(!psi) error_exit("Malloc error for psi");
	  pwgrid=malloc(3*nplwv*sizeof(int));
	  if(!pwgrid) error_exit("Malloc error for pwgrid");
          if (dict_get(m->dict,"wavecar_output")){ /* Need to collect all
                                                    * evals and occs for kpt
                                                    */
            seek_save=ftell(infile);
            for(i=1;i<nbands;i++){
              for(j=0;j<nplwv+2;j++)
                if (!fgets(buffer,LINE_SIZE,infile))
                  error_exit("Unexpected EOF");
              if (sscanf(buffer,"( %lf , %lf ) %lf",e->eval+off+i,&djunk,
                         e->occ+off+i)!=3)
                error_exit("Parse error");
            }
            fseek(infile,seek_save,SEEK_SET);
          }
	}
	for(i=0;i<nplwv;i++){
	  if (!fgets(buffer,LINE_SIZE,infile)) error_exit("Unexpected EOF");
	  if (sscanf(buffer,"%d %d %d ( %lf , %lf )",
		     pwgrid+3*i,pwgrid+3*i+1,pwgrid+3*i+2,
		     psi+2*i,psi+2*i+1)!=5)
	    error_exit("Parse error for psi component");
	}

	if (fft[0]==0){
	  for(i=0;i<nplwv;i++){
	    for(j=0;j<3;j++){
	      fft[j]=max(fft[j],pwgrid[3*i+j]);
	    }
	  }
	  for(j=0;j<3;j++)
	    fft[j]=4*fft[j]+2;
          if (debug) fprintf(stderr,"FFT grid of %dx%dx%d chosen\n",
                             fft[0],fft[1],fft[2]);
	}

        for(i=0;i<nplwv;i++){
          for(j=0;j<3;j++){
            if (pwgrid[3*i+j]<0) pwgrid[3*i+j]+=fft[j];
          }
        }
	
	if ((inrange(ns,e->spin_range))&&
	    (inrange(nb+1,e->band_range))&&
	    (inrange(ikpt+1,e->kpt_range)))
	  band_process(psi,fft,pwgrid,nplwv,0,
                         c,&g,e,k,m,ikpt,0,ns,nb,i_grid);

	
      }
      free(psi);
      psi=NULL;
      free(pwgrid);
      pwgrid=NULL;
    }
  }
}
