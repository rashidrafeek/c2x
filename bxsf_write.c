/* Write an XCrysDen .bxsf Fermi surface file */

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
#include<math.h>

#include "c2xsf.h"

int k_in_list(struct atom *k,struct kpts *kl);
void sym_vec(struct atom *a, struct atom *b, struct sym_op *s,
             double recip[3][3], int tr);

void bxsf_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                struct es *elect, struct kpts *kp,
                struct symmetry *rs){
  int i,j,k,ii,jj,kk,ic,off,ind,nb,ns,hit,okay,warn1,warn2;
  double mins[3],maxes[3],scale,lscale,frac1[3],frac2[3],chg,efermi;
  int grid[3],gpts,nbands,nspins;
  int *mapping;
  struct symmetry *ks;
  struct atom ak1;

  scale=1;
  if (flags&AU) scale=1/H_eV;
  lscale=1;
  if (flags&AU) lscale=BOHR;
  
  
  if (!elect->eval)
    error_exit("No evalues found, so cannot write bxsf");
  
  /* Force kpoint co-ords to 0<=x<1 */

  for(i=0;i<kp->n;i++)
    for(j=0;j<3;j++)
      kp->kpts[i].frac[j]=fmod(kp->kpts[i].frac[j]+1.0,1.0);

  addabs(kp->kpts,kp->n,c->basis);

  /* First check that the kpoint grid includes the gamma point */

  hit=0;
  for(i=0;i<kp->n;i++){
    if ((aeq(kp->kpts[i].frac[0],0))&&
        (aeq(kp->kpts[i].frac[1],0))&&
        (aeq(kp->kpts[i].frac[2],0))){
      hit=1;
      break;
    }
  }

  if (hit==0) error_exit("Gamma point kpt must be included for .bxsf output");

  if (!elect->e_fermi){
    efermi=NAN;
    chg=elect->nel;
    if (chg==0){
      for(i=0;i<m->n;i++)
	chg+=m->atoms[i].chg;
      if ((chg)&&(elect->charge)) chg-=*elect->charge;
      if (chg)
	fprintf(stderr,"Assuming nel=%lf from total ionic charge\n",chg);
    }
    if (chg) efermi=calc_efermi(elect,kp,chg);
    if (isnan(efermi)){
      efermi=0;
      fprintf(stderr,"Warning: no Fermi energy found\n");
    }
    else
      fprintf(stderr,"Warning: using estimate of Fermi energy\n");
  }
  else efermi=*elect->e_fermi;
  

  if (kp->mp){
    for(i=0;i<3;i++) grid[i]=kp->mp->grid[i];
  }
  else{
  
    mins[0]=mins[1]=mins[2]=1;
    maxes[0]=maxes[1]=maxes[2]=0;
    for(i=0;i<kp->n;i++){
      for(j=0;j<3;j++){
	if (!aeq(kp->kpts[i].frac[j],0)){
	  mins[j]=min(mins[j],fabs(kp->kpts[i].frac[j]));
	}
	if (!aeq(fabs(kp->kpts[i].frac[j]),1)){
	  maxes[j]=max(maxes[j],fabs(kp->kpts[i].frac[j]));
	}
      }
    }
    
    for(i=0;i<3;i++) mins[i]=min(mins[i],fabs(1-maxes[i]));

    if (debug>1)
      fprintf(stderr,"kpt minima: %f %f %f\n",mins[0],mins[1],mins[2]);
  
    for(i=0;i<3;i++){
      if (aeq(1.0/mins[i],floor(1.0/mins[i]+0.5)))
	grid[i]=floor(1.0/mins[i]+0.5);
      else
	error_exit("Unable to deduce k-point mesh parameters");
    }
  }

  fprintf(stderr,"Deduced kpoint grid: %dx%dx%d\n",grid[0],grid[1],grid[2]);
  
  gpts=grid[0]*grid[1]*grid[2];

  if ((rs->n==0)&&(gpts>2*kp->n)){
    fprintf(stderr,"No symmetry operations given, "
            "and kpoints appear to be symmetrised\n");
#ifdef SPGLIB
    if (m->n>0){
      fprintf(stderr,"Attempting to find symmetry operations\n");
      cspg_op(c,m,rs,kp,CSPG_SYM,tol);
    }
    else
      error_exit("Unable to find symmetry operations as no atoms found");
#else
    error_exit("Unable to find symmetry operations as not linked with spglib");
#endif
  }

  mapping=malloc(gpts*sizeof(int));
  if (!mapping) error_exit("Malloc error for mapping");

  ks=malloc(sizeof(struct symmetry));
  if (!ks) error_exit("malloc error");
  sym2ksym(rs,ks);

  warn1=warn2=0;
  for(i=0;i<gpts;i++) mapping[i]=-1;
  for(i=0;i<kp->n;i++){
    addabs(kp->kpts+i,1,c->recip);
    for(ii=0;ii<ks->n;ii++){
      sym_vec(kp->kpts+i,&ak1,ks->ops+ii,c->basis,0);
      okay=1;
      for(j=0;j<3;j++){
	frac1[j]=fmod(ak1.frac[j],1.0);
	if (frac1[j]<0) frac1[j]+=1.0;
	if (frac1[j]>1-tol) frac1[j]=0;
	frac2[j]=floor(frac1[j]*grid[j]+0.5);
	if (!aeq(frac1[j],(frac2[j]/grid[j]))) okay=0;
      }
      if (okay==0)
	warn1=1;
      else{
	off=frac2[2]+grid[2]*frac2[1]+grid[1]*grid[2]*frac2[0];
	if ((off<0)||(off>=gpts))
	  fprintf(stderr,"Error: %lf %lf %lf\n",frac2[0],frac2[1],frac2[2]);
	if (mapping[off]==-1) mapping[off]=i;
	else if (mapping[off]!=i) {
	  warn2=1;
	  fprintf(stderr,"Kpts %d and %d appear identical\n",i,mapping[off]);
	}
      }
    }
  }

  hit=0;
  for(i=0;i<gpts;i++)
    if (mapping[i]==-1) hit++;

  if (hit) {
    fprintf(stderr,"Error: no mapping found for %d grid points\n",hit);
    exit(1);
  }
  
  if (warn1)
    fprintf(stderr,"Warning: symmetry generates extra positions outside grid\n");
  if (warn2)
    fprintf(stderr,"Warning: superfluous k-points\n");
  
  fprintf(stderr,"nbspins=%d\n",elect->nbspins);
  /* Count bands */
  nbands=0;
  for(i=0;i<elect->nbands;i++)
    if (inrange(i+1,elect->band_range)) nbands++;
  /* Count spins */
  nspins=0;
  for(i=0;i<elect->nbspins;i++)
    if (inrange(i,elect->spin_range)) nspins++;
  
  fprintf(outfile,"BEGIN_INFO\n");
  fprintf(outfile,"# Generated by c2x\n");
  if (m->title) fprintf(outfile,"# %s\n",m->title);
  if (flags&AU)
    fprintf(outfile,"# Units: Hartrees and Bohr^-1\n");
  else
    fprintf(outfile,"# Units: eV and Angstrom^-1\n");
  if ((elect->nbspins>1)||(nbands!=elect->nbands)){
    fprintf(outfile,"#  Band number mapping\n");
    fprintf(outfile,"# number  label   band spin\n");
    i=1;
    for(ns=0;ns<elect->nbspins;ns++){
      if (!inrange(ns,elect->spin_range)) continue;
      for(nb=0;nb<elect->nbands;nb++){
        if (!inrange(nb+1,elect->band_range)) continue;
        fprintf(outfile,"# %3d     %4d   %3d    %1d\n",
                i,nb+1+elect->nbands*ns,nb+1,ns);
        i++;
      }
    }
  }
  fprintf(outfile,"  Fermi Energy: %.6f\n",efermi*scale);
  fprintf(outfile,"END_INFO\n\n");

  fprintf(outfile,"BEGIN_BLOCK_BANDGRID_3D\n");
  fprintf(outfile,"Comment\n");
  if (flags&AU)
    fprintf(outfile,"BEGIN_BANDGRID_3D_Ha\n");
  else
    fprintf(outfile,"BEGIN_BANDGRID_3D_eV\n");
  fprintf(outfile,"  %d\n",nbands*nspins);
  fprintf(outfile,"  %4d %4d %4d\n",grid[0]+1,grid[1]+1,grid[2]+1);
  fprintf(outfile,"  0.0 0.0 0.0\n");
  for(i=0;i<3;i++)
    fprintf(outfile,"  % 8f  % 8f  %8f\n",c->recip[i][0]*lscale,
            c->recip[i][1]*lscale,c->recip[i][2]*lscale);
  for(ns=0;ns<elect->nbspins;ns++){
    if (!inrange(ns,elect->spin_range)) continue;
    for(nb=0;nb<elect->nbands;nb++){
      if (!inrange(nb+1,elect->band_range)) continue;
      fprintf(outfile,"  BAND: %3d",nb+1+elect->nbands*ns);

 /* XCrysDen likes bizarre "general" grids */
      ic=0;
      for(i=0;i<=grid[0];i++){
        ii=i%grid[0];
        for(j=0;j<=grid[1];j++){
          jj=j%grid[1];
          off=jj*grid[2]+ii*grid[2]*grid[1];
          for(k=0;k<=grid[2];k++){
            kk=k%grid[2];
            ind=kk+off;
            if (ic%5==0) fprintf(outfile,"\n  ");
            fprintf(outfile,"  % 8f",elect->eval[nb+mapping[ind]*elect->nbands*elect->nbspins+elect->nbands*ns]*scale);
            ic++;
          }
        }
      }
      fprintf(outfile,"\n");
    }
  }
  fprintf(outfile,"END_BANDGRID_3D\n");
  fprintf(outfile,"END_BLOCK_BANDGRID_3D\n");

  free(mapping);
  print_bandwidths(elect,kp);
}
