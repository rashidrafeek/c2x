#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "c2xsf.h"

void fft3d(double *c, int *ngptar, int dir);

void ion_recrho(struct unit_cell *c, struct contents *m,
                int fft[3], double *cgrid, double musq);

void es_pot(struct unit_cell *c, struct contents *m,
            struct grid *g, double musq){
  int i,j,k,ii,jj,kk,ngx,ngy,ngz;
  int grid_size,fft[3],ffft[3];
  double *cgrid,*cigrid,gvec[3],scale,gsq;

  if (debug) fprintf(stderr,"Calculating ES potential\n");
  
  /* First calculate electronic contribution to the ES potential */

  for(i=0;i<3;i++) fft[i]=g->size[i];
  grid_size=fft[0]*fft[1]*fft[2];
  cgrid=malloc(2*sizeof(double)*grid_size);

  if (!cgrid) error_exit("Malloc error in es_pot() for cgrid");

  for(i=0;i<grid_size;i++){
    cgrid[2*i]=-g->data[i];
    cgrid[2*i+1]=0;
  }

  /* FFT to reciprocal space */

  ffft[0]=fft[2];
  ffft[1]=fft[1];
  ffft[2]=fft[0];

  if (debug>1) fprintf(stderr,"first FFT in esp_pot\n");

  fft3d(cgrid,ffft,-1);

  /* Now add ionic potentials */

  cigrid=malloc(2*sizeof(double)*grid_size);
  if (!cigrid) error_exit("Malloc error in es_pot() for cigrid");
  ion_recrho(c,m,fft,cigrid,musq);
  for(i=0;i<2*grid_size;i++)
    cgrid[i]+=cigrid[i];
  free(cigrid);
  
  if (debug>1) fprintf(stderr,
                       "g=0 component of total charge in esp: %lf+%lfi\n",
                       cgrid[0],cgrid[1]);

  if (fabs(cgrid[0])>1e-5){
    fprintf(stderr,"Warning: cell not neutral generating ES pot."
            "g=0 component is %lf\n",cgrid[0]);
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
  int grid_size;
  double scale,scale_r,scale_i,phi;
  double gvec[3],gsq;
  double *cgrid_ps,*cgrid_bare;

  bare=0;
  if (musq<0){
    musq=-musq;
    bare=1;
  }

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

  for(i=0;i<2*grid_size;i++) cgrid[i]=0;
  
  for(i=0;i<ngx;i++){
    ii=i;
    if (ii>ngx/2) ii=ii-ngx;
    for(j=0;j<ngy;j++){
      jj=j;
      if (jj>ngy/2) jj=jj-ngy;
      for(k=0;k<ngz;k++){
        kk=k;
        if (k>ngz/2) kk=kk-ngz;
        for(ion=0;ion<m->n;ion++){
          phi=m->atoms[ion].frac[0]*ii+
            m->atoms[ion].frac[1]*jj+
            m->atoms[ion].frac[2]*kk;
          phi*=-2*M_PI;
          scale_r=m->atoms[ion].chg*cos(phi);
          scale_i=m->atoms[ion].chg*sin(phi);
          if ((bare)||(m->atoms[ion].chg==m->atoms[ion].atno)){
            cgrid[2*(k+fft[2]*(j+i*fft[1]))]+=
              scale_r*cgrid_bare[k+fft[2]*(j+i*fft[1])];
            cgrid[2*(k+fft[2]*(j+i*fft[1]))+1]+=
              scale_i*cgrid_bare[k+fft[2]*(j+i*fft[1])];
          }
          else{
            cgrid[2*(k+fft[2]*(j+i*fft[1]))]+=
              scale_r*cgrid_ps[k+fft[2]*(j+i*fft[1])];
            cgrid[2*(k+fft[2]*(j+i*fft[1]))+1]+=
              scale_i*cgrid_ps[k+fft[2]*(j+i*fft[1])];
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
