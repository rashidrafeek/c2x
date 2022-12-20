/* Copyright (c) 2019-2020 MJ Rutter 
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
#include "c2xsf.h"

void fft3d(double *c, int *ngptar, int dir);

void ion_recrho(struct unit_cell *c, struct contents *m,
                int fft[3], double *cgrid, double musq);

void es_pot(struct unit_cell *c, struct contents *m,
            struct grid *g, struct es *elect, double musq){
  int i,j,k,ii,jj,kk,ngx,ngy,ngz;
  int grid_size,fft[3],ffft[3],ind[3],dipole_slab_dir;
  double *cgrid,*cigrid,gvec[3],scale,gsq,net_charge,ion_charge;
  double dpole[3],field,off,*dipole_ctr,q_off,qa,qc,ql;
  double vec[3],mag,abc[6];

  dipole_ctr=elect->dip_ctr;
  cgrid=NULL;
  
  if (debug) fprintf(stderr,"Calculating ES potential\n");

  if ((g->name)&&(!strncasecmp(g->name,"pot",3)))
    fprintf(stderr,"Warning: need electron density to calculate potential,\n"
            "         but appear to have potential in grid.\n");

  /* First calculate electronic contribution to the ES potential */

  for(i=0;i<3;i++) fft[i]=g->size[i];
  grid_size=fft[0]*fft[1]*fft[2];
  ffft[0]=fft[2];
  ffft[1]=fft[1];
  ffft[2]=fft[0];

  if (g->data){
    cgrid=malloc(2*sizeof(double)*grid_size);

    if (!cgrid) error_exit("Malloc error in es_pot() for cgrid");

    for(i=0;i<grid_size;i++){
      cgrid[2*i]=-g->data[i];
      cgrid[2*i+1]=0;
    }

    /* FFT to reciprocal space */


    if (debug>1) fprintf(stderr,"first FFT in esp_pot\n");

    fft3d(cgrid,ffft,-1);

    ion_charge=0;
    for(i=0;i<m->n;i++) ion_charge+=m->atoms[i].chg;
  }
  else{
    if (!grid_size){
      basis2abc(c->basis,abc);
      for(i=0;i<3;i++) fft[i]=10*abc[i]+1;
      if (debug)
	fprintf(stderr,"Using grid size of %dx%dx%d\n",fft[0],fft[1],fft[2]);
      for(i=0;i<3;i++) g->size[i]=fft[i];
      grid_size=fft[0]*fft[1]*fft[2];
      ffft[0]=fft[2];
      ffft[1]=fft[1];
      ffft[2]=fft[0];
    }
      
    g->data=malloc(grid_size*sizeof(double));
    if(!g->data) error_exit("malloc error for grid");

    net_charge=0;
    mag=0;
    ion_charge=0;
    for(i=0;i<m->n;i++){
      net_charge+=m->atoms[i].site_chg;
      mag+=m->atoms[i].site_chg*m->atoms[i].site_chg;
    }
    if (mag==0){
      fprintf(stderr,"Cannot calculate ES potential with no electron density\n"
	      "and no site charges\n");
      return;
    }
    if (fabs(net_charge)>tol)
      fprintf(stderr,"Warning: net charge of %lf e from site charges!\n",
	      net_charge);
    dict_strcat(m->dict,"pot_from_site_chg","");
  }
  
  /* Now add ionic potentials */

  cigrid=malloc(2*sizeof(double)*grid_size);
  if (!cigrid) error_exit("Malloc error in es_pot() for cigrid");
  ion_recrho(c,m,fft,cigrid,musq);
  if (cgrid){
    for(i=0;i<2*grid_size;i++)
      cgrid[i]+=cigrid[i];
    free(cigrid);
  }
  else cgrid=cigrid;

  net_charge=cgrid[0]*c->vol/(fft[0]*fft[1]*fft[2]);
  
  if ((debug>1)||(fabs(net_charge)>1e-7*ion_charge)) {
    fprintf(stderr,
            "g=0 component of total charge in esp: %lf+%lfi\n",
            cgrid[0],cgrid[1]);
    fprintf(stderr,"          corresponding net charge is %g\n",net_charge);
  }    

  /*  Now need to scale by 1/g^2, or, more precisely,
   *  1/(EPS0*g^2)
   */

  scale=1/(EPS0*grid_size);

  ngx=g->size[0];
  ngy=g->size[1];
  ngz=g->size[2];

  for(i=0;i<ngx;i++){
    ii=i;
    if (ii>ngx/2) ii=ii-ngx;
    for(j=0;j<ngy;j++){
      jj=j;
      if (jj>ngy/2) jj=jj-ngy;
      for(k=0;k<ngz;k++){
        kk=k;
        if (kk>ngz/2) kk=kk-ngz;
        if ((ii==0)&&(jj==0)&&(kk==0)){
          cgrid[2*(kk+fft[2]*(jj+ii*fft[1]))]=0;
          cgrid[2*(kk+fft[2]*(jj+ii*fft[1]))+1]=0;
          continue;
        }
        gvec[0]=ii*c->recip[0][0]+jj*c->recip[1][0]+kk*c->recip[2][0];
        gvec[1]=ii*c->recip[0][1]+jj*c->recip[1][1]+kk*c->recip[2][1];
        gvec[2]=ii*c->recip[0][2]+jj*c->recip[1][2]+kk*c->recip[2][2];
        gsq=4*M_PI*M_PI*(gvec[0]*gvec[0]+gvec[1]*gvec[1]+gvec[2]*gvec[2]);
        cgrid[2*(k+fft[2]*(j+i*fft[1]))]*=scale/gsq;
        cgrid[2*(k+fft[2]*(j+i*fft[1]))+1]*=scale/gsq;
      }
    }
  }

  /* FFT back to real space */

  if (debug>1) fprintf(stderr,"second FFT in es_pot\n");

  fft3d(cgrid,ffft,1);

  
  if (dipole_ctr){ /* Need to correct for dipole */
    dipole_calc(c,m,g,dipole_ctr,dpole);
    dipole_slab_dir=-1;
    if (elect->dip_corr_dir){
      dipole_slab_dir=elect->dip_corr_dir[0]-'a';
    }
    if ((dipole_slab_dir>=0)&&(dipole_slab_dir<3)){
      /* Find unit vector in dipole_slab_dir */
      for(i=0;i<3;i++)
        vec[i]=c->basis[dipole_slab_dir][i];
      mag=0;
      for(i=0;i<3;i++)
        mag+=vec[i]*vec[i];
      mag=sqrt(mag);
      /* Dot dipole with this to find field */
      field=0;
      for(i=0;i<3;i++)
        field+=dpole[i]*vec[i]/mag;
      field=field/(EPS0*c->vol);
      fprintf(stderr,"Adding dipole correction of %.3f V/A to potential\n",
              field);
      field=field*mag; /* in volt/unit cell */
      field=field/g->size[dipole_slab_dir];  /* in volt/grid cell */
      fprintf(stderr,"g->size=[%d,%d,%d]\n",g->size[0],g->size[1],
              g->size[2]);
      for(i=0;i<ngx;i++){
        ind[0]=i;
        for(j=0;j<ngy;j++){
          ind[1]=j;
          for(k=0;k<ngz;k++){
            ind[2]=k;
            off=(ind[dipole_slab_dir]-0.5*g->size[dipole_slab_dir])*field;
            cgrid[2*(k+fft[2]*(j+i*fft[1]))]+=off;
          }
        }
      }
    }
  }

  if ((dict_get(m->dict,"charge_correction"))&&(*elect->dip_corr!='m')){
    if (!dipole_ctr){
      dipole_ctr=malloc(3*sizeof(double));
      if (!dipole_ctr) error_exit("Malloc error for three doubles");
      for(i=0;i<3;i++) dipole_ctr[i]=0.5;
    }
    dipole_calc(c,m,g,dipole_ctr,dpole);
    dipole_slab_dir=elect->dip_corr_dir[0]-'a';
    if ((dipole_slab_dir<0)||(dipole_slab_dir>2))
      error_exit("Impossible dipole slab dir");
    ql=0;
    for(i=0;i<3;i++)
      ql+=c->basis[dipole_slab_dir][i]*c->basis[dipole_slab_dir][i];
    ql=sqrt(ql);
    q_off=dipole_ctr[dipole_slab_dir]+
      dpole[dipole_slab_dir]/(net_charge*ql);
    if (debug) fprintf(stderr,"Centre of charge correction: %lf\n",q_off);
    qa=net_charge/(2*c->vol*EPS0);
    qc=-net_charge*ql*ql/(24*c->vol*EPS0);
    for(i=0;i<ngx;i++){
      ind[0]=i;
      for(j=0;j<ngy;j++){
	ind[1]=j;
	for(k=0;k<ngz;k++){
	  ind[2]=k;
	  /* in grid cells */
	  off=(ind[dipole_slab_dir]-0.5*g->size[dipole_slab_dir]);
	  /* in A */
	  off*=ql/g->size[dipole_slab_dir];
	  /* in eV */
	  off=qa*off*off+qc;
	  cgrid[2*(k+fft[2]*(j+i*fft[1]))]-=off;
	}
      }
    }
  }
  
  /* copy back to real grid */
  
  for(i=0;i<grid_size;i++){
    g->data[i]=cgrid[2*i];
  }

  free(cgrid);
  
  g->name=malloc(30);
  if (!g->name) error_exit("Malloc error for grid name");
  sprintf(g->name,"Potential_Volts");

  if (debug){
    double min,max;
    min=1e20;
    max=-1e20;
    for(i=0;i<grid_size;i++){
      if (g->data[i]>max) max=g->data[i];
      if (g->data[i]<min) min=g->data[i];
    }
    fprintf(stderr,"Potential in volts on %dx%dx%d grid: min=%g  max=%g  \n",
            fft[0],fft[1],fft[2],min,max);
  }
  
}


/* Treat ions as Gaussian blobs of charge */
void ion_recrho(struct unit_cell *c, struct contents *m,
                int fft[3], double *cgrid, double musq){
  int i,j,k,ii,jj,kk,ngx,ngy,ngz,ion,bare;
  int grid_size,partial,offset;
  double scale,scale_r,scale_i,phi,cg,cg1;
  double gvec[3],gsq;
  double *cgrid_ps,*cgrid_bare;
  double dtmp,phase_r,phase_i,scale_r_org,scale_i_org,*pot;
  
  bare=0;
  if (musq<0){
    musq=-musq;
    bare=1;
  }

  partial=0;
  if (dict_get(m->dict,"pot_from_site_chg")) partial=1;
  
  if(debug>1){
    fprintf(stderr,"Using smearing of 1/%.3f A",sqrt(musq));
    if (bare) fprintf(stderr," and bare Coulomb potential\n");
    else fprintf(stderr," and both bare and pseudo Coulomb potentials\n");
  }
  
  /* Consider unit charge at the origin, but with Gaussian smearing
   * over exp(-mu^2 r^2)
   * Its FFT will be exp(-g^2/4mu^2) (and is real)
   * For consistency with the units in the charge density,
   * the normalisation is that the g=0 component is
   * grid size / cell volume
   */

  grid_size=fft[0]*fft[1]*fft[2];
  ngx=fft[0];
  ngy=fft[1];
  ngz=fft[2];

  cgrid_bare=malloc(grid_size*sizeof(double));
  if (!cgrid_bare) error_exit("malloc error in ion_recrho");
  cgrid_ps=malloc(grid_size*sizeof(double));
  if (!cgrid_ps) error_exit("malloc error in ion_recrho");
  
  scale=grid_size/c->vol;

  if (debug>1){
    fprintf(stderr,"Aliasing test in ion_recrho():");
    for(i=0;i<3;i++){
      gsq=c->recip[i][0]*c->recip[i][0]+
        c->recip[i][1]*c->recip[i][1]+
        c->recip[i][2]*c->recip[i][2];
      gsq*=M_PI*M_PI*fft[i]*fft[i];
      fprintf(stderr," %lg",exp(-gsq/(4*musq)));
    }
    fprintf(stderr,"\n");
    if (!bare){
      fprintf(stderr,"                              ");
      for(i=0;i<3;i++){
        gsq=c->recip[i][0]*c->recip[i][0]+
          c->recip[i][1]*c->recip[i][1]+
          c->recip[i][2]*c->recip[i][2];
        gsq*=M_PI*M_PI*fft[i]*fft[i];
        fprintf(stderr," %lg",(1-gsq/(6*musq))*exp(-gsq/(4*musq)));
      }
      fprintf(stderr,"\n");
    }
  }
  
  for(i=0;i<ngx;i++){
    ii=i;
    if (ii>ngx/2) ii=ii-ngx;
    for(j=0;j<ngy;j++){
      jj=j;
      if (jj>ngy/2) jj=jj-ngy;
      for(k=0;k<ngz;k++){
        kk=k;
        if (k>ngz/2) kk=kk-ngz;
        gvec[0]=ii*c->recip[0][0]+jj*c->recip[1][0]+kk*c->recip[2][0];
        gvec[1]=ii*c->recip[0][1]+jj*c->recip[1][1]+kk*c->recip[2][1];
        gvec[2]=ii*c->recip[0][2]+jj*c->recip[1][2]+kk*c->recip[2][2];
        gsq=4*M_PI*M_PI*(gvec[0]*gvec[0]+gvec[1]*gvec[1]+gvec[2]*gvec[2]);
        cgrid_bare[k+fft[2]*(j+i*fft[1])]=scale*exp(-gsq/(4*musq));
        cgrid_ps[k+fft[2]*(j+i*fft[1])]=scale*(1-gsq/(6*musq))*exp(-gsq/(4*musq));
      }
    }
  }

  /* For the cell, we need to add contributions from each ion.
   * Each ion contributes a scale factor of its charge, and a phase
   * factor arising from its position of exp(igr)
   */

  if ((m->n<1000)||(flags&HIPREC)){
    for(i=0;i<ngx;i++){
      ii=i;
      if (ii>ngx/2) ii=ii-ngx;
      for(j=0;j<ngy;j++){
	jj=j;
	if (jj>ngy/2) jj=jj-ngy;
	for(k=0;k<ngz;k++){
	  kk=k;
	  if (k>ngz/2) kk=kk-ngz;
	  cg=0;
	  cg1=0;
	  offset=k+fft[2]*(j+i*fft[1]);
	  for(ion=0;ion<m->n;ion++){
	    phi=m->atoms[ion].frac[0]*ii+
	      m->atoms[ion].frac[1]*jj+
	      m->atoms[ion].frac[2]*kk;
	    phi*=-2*M_PI;
	    /* Hope compiler uses a combined sincos call here... */
	    scale_r=cos(phi);
	    scale_i=sin(phi);
	    if (partial){
	      scale_r*=m->atoms[ion].site_chg;
	      scale_i*=m->atoms[ion].site_chg;
	    }
	    else{
	      scale_r*=m->atoms[ion].chg;
	      scale_i*=m->atoms[ion].chg;
	    }
	    if ((bare)||(m->atoms[ion].chg==m->atoms[ion].atno)){
	      cg+=scale_r*cgrid_bare[offset];
	      cg1+=scale_i*cgrid_bare[offset];
	    }
	    else{
	      cg+=scale_r*cgrid_ps[offset];
	      cg1+=scale_i*cgrid_ps[offset];
	    }
	  }
	  cgrid[2*offset]=cg;
	  cgrid[2*offset+1]=cg1;
	}
      }
    }

  }
  else{
    fprintf(stderr,"Using fast ionic potential calc as >=1000 ions\n");
    for(i=0;i<2*ngx*ngy*ngz;i++) cgrid[i]=0;

    if (bare) pot=cgrid_bare;
    else pot=cgrid_ps;
  
    for(i=0;i<ngx;i++){
      ii=i;
      if (ii>ngx/2) ii=ii-ngx;
      for(j=0;j<ngy;j++){
	jj=j;
	if (jj>ngy/2) jj=jj-ngy;
	for(ion=0;ion<m->n;ion++){
	  offset=fft[2]*(j+i*fft[1]);
	  phi=m->atoms[ion].frac[0]*ii+
	    m->atoms[ion].frac[1]*jj;
	  phi*=-2*M_PI;
	  scale_r=cos(phi);
	  scale_i=sin(phi);
	  if (partial){
	    scale_r*=m->atoms[ion].site_chg;
	    scale_i*=m->atoms[ion].site_chg;
	  }
	  else{
	    scale_r*=m->atoms[ion].chg;
	    scale_i*=m->atoms[ion].chg;
	  }
	  cgrid[2*offset]+=scale_r*pot[offset];
	  cgrid[2*offset+1]+=scale_i*pot[offset];

	  phi=-2*M_PI*m->atoms[ion].frac[2];
	  phase_r=cos(phi);
	  phase_i=sin(phi);
	  scale_r_org=scale_r;
	  scale_i_org=scale_i;
	  for(k=1;k<=ngz/2;k++){
	    dtmp=scale_r;
	    scale_r=scale_r*phase_r-scale_i*phase_i;
	    scale_i=scale_i*phase_r+dtmp*phase_i;
	    offset+=1;
	    cgrid[2*offset]+=scale_r*pot[offset];
	    cgrid[2*offset+1]+=scale_i*pot[offset];
	  }
	
	  scale_r=scale_r_org;
	  scale_i=scale_i_org;
	  phase_i=-phase_i;
	  offset=fft[2]*(j+1+i*fft[1]);
	  for(k=-1;k>-ngz/2;k--){
	    dtmp=scale_r;
	    scale_r=scale_r*phase_r-scale_i*phase_i;
	    scale_i=scale_i*phase_r+dtmp*phase_i;
	    offset-=1;
	    cgrid[2*offset]+=scale_r*pot[offset];
	    cgrid[2*offset+1]+=scale_i*pot[offset];
	  }
	}
      }
    }
  }
  
  if(debug>2) fprintf(stderr,"g=0 term in ion_recrho is %lf+%lfi\n",
		      cgrid[0],cgrid[1]);

  free(cgrid_ps);
  free(cgrid_bare);
  
}
