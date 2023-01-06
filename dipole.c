/* Copyright (c) 2019, 2021 MJ Rutter 
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

/* References:
 *
 * Slab correction: https://doi.org/10.1103/PhysRevB.46.16067
 * Molecule in cube: https://doi.org/10.1103/PhysRevB.51.4014
 *               and https://doi.org/10.1103/PhysRevB.60.15476
 * Molecule in tetragon: https://doi.org/10.1088/1361-648X/ab20e1
 */


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "c2xsf.h"


void fft3d(double *c, int *ngptar, int dir);

/* From charge.c */
double quadrupole(struct unit_cell *c, struct contents *m,
                  struct grid *g, double *ctr);
double quadrupole_ii(struct unit_cell *c, struct contents *m,
                     struct grid *g, double *ctr, int dir);

void dipole_calc(struct unit_cell *c, struct contents *m,
                 struct grid *g, double *dipole_ctr, double *dpole){
  int i,j,ngx,ngy,ngz,ix,iy,iz;
  int nn[3],ii[3];
  double i_charge,i_dipole[3],e_charge,e_dipole[3],e_dipole_g[3],d[3];
  double ed[3],scale,disp,*ptr,*ptr2,mag;
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

  if (debug>2) fprintf(stderr,"Ed: %f %f %f\n",ed[0],ed[1],ed[2]);

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
  if ((fabs(i_charge-scale*rgrid[0])>tol)&&
      (!dict_get(m->dict,"dipole_no_charge_warn")))
    fprintf(stderr,"Warning: total charge not zero\n");
  
  fprintf(stderr,"Total dipole (eA): (");
  for(i=0;i<3;i++) fprintf(stderr,"%f%s",i_dipole[i]+e_dipole_g[i],
                           (i==2)?"":",");
  mag=0;  
  for(i=0;i<3;i++)
    mag+=(i_dipole[i]+e_dipole_g[i])*(i_dipole[i]+e_dipole_g[i]);
  mag=sqrt(mag);
  fprintf(stderr,")  magnitude %f\n",mag);

  for(i=0;i<3;i++) dpole[i]=i_dipole[i]+e_dipole_g[i];

  free(rgrid);
  
}


void dipole(struct unit_cell *c, struct contents *m,
            struct grid *g, struct es *elect){

  double dpole[3],mag,E,*dipole_ctr,v[3];
  int i,dipole_slab_dir,c_dir;
  double va,vc,dE,mag2,quad;

  dipole_ctr=elect->dip_ctr;

  fprintf(stderr,"Centre used for moments: (%f,%f,%f)\n",
	  dipole_ctr[0],dipole_ctr[1],dipole_ctr[2]);
  
  dipole_calc(c,m,g,dipole_ctr,dpole);
  
  if (debug){
    fprintf(stderr,"Extra dipole field in V/A: (");
    for(i=0;i<3;i++) fprintf(stderr,"%f%s",
                             (dpole[i])/(EPS0*c->vol),
                             (i==2)?"":",");
    fprintf(stderr,")\n");
  }

  quad=quadrupole(c,m,g,dipole_ctr);
  fprintf(stderr,"Trace of quadrupole tensor (eA^2): %f\n",quad);
  for(i=0;i<3;i++){
    quad=quadrupole_ii(c,m,g,dipole_ctr,i);
    fprintf(stderr,"         quadrupole_%c%c (eA^2):     %f\n",
	    i+'a',i+'a',quad);
  }
  
  if (!elect->dip_corr_dir) return; /* Report, don't correct */
  
  dipole_slab_dir=-1;
  if (elect->dip_corr_dir){
    dipole_slab_dir=elect->dip_corr_dir[0]-'a';
  }
  if ((dipole_slab_dir>=0)&&(dipole_slab_dir<=2)){
    /* Find unit vector corresponding to dipole slab axis */
    for(i=0;i<3;i++)
      v[i]=c->basis[dipole_slab_dir][i];
    mag=0;
    for(i=0;i<3;i++)
      mag+=v[i]*v[i];
    mag=sqrt(mag);
    for(i=0;i<3;i++)
      v[i]/=mag;
    /* Now have unit vector, so dot with dipole */
    mag=0;
    for(i=0;i<3;i++)
      mag+=v[i]*dpole[i];
    E=0.5*mag*mag/(EPS0*c->vol);
    fprintf(stderr,"Calculated dipole energy correction (%c axis): %f eV\n",
            'a'+dipole_slab_dir,E);
    if (elect->energy){
      if ((elect->dip_corr==NULL)||(elect->dip_corr[0]=='N'))
        fprintf(stderr,"Corrected energy %.6f + %.6f = %.6f eV\n",
                *elect->energy,E,*elect->energy+E);
      else
        fprintf(stderr,"Reported energy of %.6f appears to include "
                "a correction already\n",*elect->energy);
    }
  }

  /* Molecules are harder */

  if ((elect->dip_corr_dir)&&(*elect->dip_corr_dir=='m')){
    double abc[6];
    cart2abc(c,NULL,abc,NULL);
    if ((aeq(abc[3],90))&&(aeq(abc[4],90))&&(aeq(abc[3],90))){
      if ((aeq(abc[0],abc[1]))&&(aeq(abc[1],abc[2]))){
        mag=0;
        for(i=0;i<3;i++)
          mag+=dpole[i]*dpole[i];
        E=mag/(6*EPS0*c->vol);
        fprintf(stderr,"Calculated dipole energy correction (molecule in cube):"
                " %f eV\n",E);
        if (elect->energy){
          if ((elect->dip_corr==NULL)||(elect->dip_corr[0]=='N'))
            fprintf(stderr,"Corrected energy %.6f + %.6f = %.6f eV\n",
                    *elect->energy,E,*elect->energy+E);
          else
            fprintf(stderr,"Reported energy of %.6f appears to include "
                "a correction already\n",*elect->energy);      
        }
      }
      else{
	c_dir=-1;
	if (aeq(abc[0],abc[1])) c_dir=2;
	if (aeq(abc[0],abc[2])) c_dir=1;
	if (aeq(abc[1],abc[2])) c_dir=0;
	if (c_dir!=-1){
	  /* Find unit vector corresponding to tetragonal c axis */
	  for(i=0;i<3;i++)
	    v[i]=c->basis[c_dir][i];
	  mag=0;
	  for(i=0;i<3;i++)
	    mag+=v[i]*v[i];
	  mag=sqrt(mag);
	  for(i=0;i<3;i++)
	    v[i]/=mag;
	  /* Now have unit vector, so dot with dipole */
	  mag=0;
	  for(i=0;i<3;i++)
	    mag+=v[i]*dpole[i];
	  /* Calculate also dipole magnitude */
	  mag2=0;
	  for(i=0;i<3;i++)
	    mag2+=dpole[i]*dpole[i];
          mag2=sqrt(mag2);
	  if (mag2>1.05*fabs(mag))
	    fprintf(stderr,"Warning: only dipole moment parallel "
		    "to %c axis corrected\n",'a'+c_dir);
	  vc=abc[c_dir];
	  va=abc[0];
	  if (c_dir==0) va=abc[1];
	  E=((mag*mag)/(EPS0*va*va))*(1/(2*vc)-1/(3*va));
	  fprintf(stderr,
		  "Calculated dipole energy correction (%c axis):"
                  "      %10.6f eV\n",
		  'a'+c_dir,E);
	  if (elect->energy){
	    if ((elect->dip_corr==NULL)||(elect->dip_corr[0]=='N'))
	      fprintf(stderr,"Corrected energy %14.6f + %10.6f = %14.6f eV\n",
		      *elect->energy,E,*elect->energy+E);
	    else
	      fprintf(stderr,"Reported energy of %.6f appears to include "
		      "a correction already\n",*elect->energy);
          }
	  dE=-4*M_PI*mag*mag*(exp(-2*M_PI)-exp(-2*M_PI*vc/va))/(EPS0*va*va*va);
	  fprintf(stderr,"Additional correction if point dipole:"
                  "         %14.6f eV\n",dE);
	  if (elect->energy){
	    fprintf(stderr,"Leading to corrected energy of:"
		    "                %14.6f eV\n",
		    *elect->energy+E+dE);
	  }
	}
        else fprintf(stderr,"This version cannot "
                "calculate dipole corrections for non-tetragonal cells.\n");
      }
    }
    else
        fprintf(stderr,"This version cannot "
                "calculate dipole corrections for non-tetragonal cells.\n");
  }
  
}
