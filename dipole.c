#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "c2xsf.h"

int dipole_slab_dir=0;

void fft3d(double *c, int *ngptar, int dir);


void dipole(struct unit_cell *c, struct contents *m,
            struct grid *g, double *dipole_ctr){
  int i,j,ngx,ngy,ngz,ix,iy,iz;
  int nn[3],ii[3];
  double i_charge,i_dipole[3],e_charge,e_dipole[3],e_dipole_g[3],d[3];
  double ed[3],scale,disp,*ptr,*ptr2,mag,E;
  double *rgrid;
  double theta;


  i_charge=0;
  d[0]=d[1]=d[2]=0;
  for(i=0;i<m->n;i++){
    i_charge+=m->atoms[i].chg;
    for(j=0;j<3;j++){
      /* force disp to range 0.5<=disp<=0.5 */
      disp=fmod(m->atoms[i].frac[j]-dipole_ctr[j]+0.5,1.0);
      if (disp<0) disp+=1;
      disp-=0.5;
      /* ignore ions split between +/- 0.5 */
      if ((aeq(disp,0.5))||(aeq(disp,-0.5))) disp=0;
      d[j]+=m->atoms[i].chg*disp;
    }
  }

  /* Convert ionic dipole from fractional to absolute units */

  for(i=0;i<3;i++){
    i_dipole[i]=0;
    for(j=0;j<3;j++)
      i_dipole[i]+=d[j]*c->basis[j][i];
  }

  if (debug){
    fprintf(stderr,"Ionic charge is %f\n",i_charge);
    fprintf(stderr,"Ionic dipole is (%f,%f,%f)\n",
            i_dipole[0],i_dipole[1],i_dipole[2]);
  }

  ptr=g->data;
  if (!ptr) return;
  
  e_charge=0;
  ed[0]=ed[1]=ed[2]=0;
  scale=c->vol/(g->size[0]*g->size[1]*g->size[2]);

  for(ii[0]=0;ii[0]<g->size[0];ii[0]++){
    for(ii[1]=0;ii[1]<g->size[1];ii[1]++){
      for(ii[2]=0;ii[2]<g->size[2];ii[2]++){
        e_charge-=scale*(*ptr);
        for(i=0;i<3;i++){
          disp=(double)ii[i]/g->size[i]-dipole_ctr[i]; /* BUG!!! */
      /* force disp to range 0.5<=disp<=0.5 */
          disp=fmod(disp+0.5,1.0);
          if (disp<0) disp+=1;
          disp-=0.5;
          ed[i]-=disp*scale*(*ptr);
        }
        ptr++;
      }
    }
  }

  /* Convert electronic dipole from fractional to absolute units */

  for(i=0;i<3;i++){
    e_dipole[i]=0;
    for(j=0;j<3;j++)
      e_dipole[i]+=ed[j]*c->basis[j][i];
  }

  if (debug>1){
    fprintf(stderr,"Electronic charge is %f\n",e_charge);
    fprintf(stderr,"Electronic dipole is (%f,%f,%f)\n",
            e_dipole[0],e_dipole[1],e_dipole[2]);
  }

  /* Try calculating the electronic dipole again in reciprocal space */

  rgrid=malloc(2*g->size[0]*g->size[1]*g->size[2]*sizeof(double));
  if(!rgrid) error_exit("Malloc error in dipole.c");

  /* Copy real grid data to complex rgrid prior to FFT */

  ptr=g->data;
  ptr2=rgrid;
  for(i=0;i<g->size[0]*g->size[1]*g->size[2];i++){
    *(ptr2++)=*(ptr++);
    *(ptr2++)=0;
  }

  /* FFT */

  nn[0]=g->size[2];
  nn[1]=g->size[1];
  nn[2]=g->size[0];
  
  fft3d(rgrid,nn,-1);

  ngx=g->size[0];
  ngy=g->size[1];
  ngz=g->size[2];

  ed[0]=ed[1]=ed[2]=0;

  for(i=1;i<ngx;i++){
    ix=i;
    if (ix>ngx/2) ix=ix-ngx;
    theta=2*M_PI*(dipole_ctr[0]-0.5)*ix;
    ed[0]-=(rgrid[2*i*ngy*ngz]*sin(theta)+rgrid[2*i*ngy*ngz+1]*cos(theta))/ix;
    //    ed[0]+=(rgrid[2*i]*sin(theta)+rgrid[2*i+1]*cos(theta))/ix;
  }

  for(i=1;i<ngy;i++){
    iy=i;
    if (iy>ngy/2) iy=iy-ngy;
    theta=2*M_PI*(dipole_ctr[1]-0.5)*iy;
    ed[1]-=(rgrid[2*i*ngz]*sin(theta)+rgrid[2*i*ngz+1]*cos(theta))/iy;
    //    ed[1]+=(rgrid[2*i*ngx]*sin(theta)+rgrid[2*i*ngx+1]*cos(theta))/iy;
  }

  for(i=1;i<ngz;i++){
    iz=i;
    if (iz>ngz/2) iz=iz-ngz;
    theta=2*M_PI*(dipole_ctr[2]-0.5)*iz;
    ed[2]-=(rgrid[2*i]*sin(theta)+rgrid[2*i+1]*cos(theta))/iz;
  }

  for(i=0;i<3;i++) ed[i]*=scale/(2*M_PI);

  fprintf(stderr,"Ed: %f %f %f\n",ed[0],ed[1],ed[2]);
  
  /* Convert electronic dipole from fractional to absolute units */

  for(i=0;i<3;i++){
    e_dipole_g[i]=0;
    for(j=0;j<3;j++)
      e_dipole_g[i]+=ed[j]*c->basis[j][i];
  }

  if (debug){
    fprintf(stderr,"Elect (g)  dipole is (%f,%f,%f)\n",
            e_dipole_g[0],e_dipole_g[1],e_dipole_g[2]);
    if (debug>1) fprintf(stderr,"g=0 change density is %f+%fi\n",
                         scale*rgrid[0],scale*rgrid[1]);
  }


  fprintf(stderr,"Total charge: %f\n",i_charge-scale*rgrid[0]);
  if (fabs(i_charge-scale*rgrid[0])>tol)
    fprintf(stderr,"Warning: total charge not zero\n");
  
  fprintf(stderr,"Total dipole (eA): (");
  for(i=0;i<3;i++) fprintf(stderr,"%f%s",i_dipole[i]+e_dipole_g[i],
                           (i==2)?"":",");
  mag=0;  
  for(i=0;i<3;i++)
    mag+=(i_dipole[i]+e_dipole_g[i])*(i_dipole[i]+e_dipole_g[i]);
  mag=sqrt(mag);
  fprintf(stderr,")  magnitude %f\n",mag);

  if (debug){
    fprintf(stderr,"Extra dipole field in V/A: (");
    for(i=0;i<3;i++) fprintf(stderr,"%f%s",
                             (i_dipole[i]+e_dipole_g[i])/(EPS0*c->vol),
                             (i==2)?"":",");
    fprintf(stderr,")\n");
  }

  if ((dipole_slab_dir>=1)&&(dipole_slab_dir<=3)){
    mag=i_dipole[dipole_slab_dir-1]+e_dipole_g[dipole_slab_dir-1];
    E=0.5*mag*mag/(EPS0*c->vol);
    fprintf(stderr,"Calculated dipole energy correction (%c axis): %f eV\n",
            'w'+dipole_slab_dir,E);
  }
  
}
