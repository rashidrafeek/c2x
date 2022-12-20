/* Calculate Madelung constant and leading term in correction for
 * charged cells. Also calculate quadrupole moment as needed for
 * 0D charged cell correction (Makov and Payne, PRB 51 4014 (1995))
 * and 2D charged cell correction (Rutter, Electron Struct (2021)) */

/* Copyright (c) 2019 -- 2021 MJ Rutter 
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

double madelung(struct unit_cell *c){
  double M, M_real,M_recip,M_self,M_charged;
  double E_real,E_recip,E_self,E_charged;
  double dist,disp[3];
  double sigma;
  int i,j,k,ii,jj,kk,m;

  sigma=M_PI*pow(0.7/(c->vol*c->vol),1.0/3.0);
  if (debug>1) fprintf(stderr,"sigma=%f\n",sigma);

  /* Real space sum */

  M_real=0.0;

  /* erfc(x) is approx exp(-x*x),
   * we don't care about terms smaller than 1e-16,
   * which is exp(-6*6) */
  ii=1+6.0/sqrt(sigma)/sqrt(vmod2(c->basis[0]));
  jj=1+6.0/sqrt(sigma)/sqrt(vmod2(c->basis[1]));
  kk=1+6.0/sqrt(sigma)/sqrt(vmod2(c->basis[2]));
  if (debug>1)
    fprintf(stderr,"Real space Ewald grid=%dx%dx%d\n",2*ii+1,2*jj+1,2*kk+1);
  for(i=-ii;i<=ii;i++){
    for(j=-jj;j<=jj;j++){
      for(k=-kk;k<=kk;k++){
	if ((i==0)&&(j==0)&&(k==0)) continue;
	for(m=0;m<3;m++) disp[m]=c->basis[0][m]*i+
			   c->basis[1][m]*j+c->basis[2][m]*k;
	dist=sqrt(vmod2(disp));
	M_real+=(1.0/dist)*erfc(dist*sqrt(sigma));
      }
    }
  }

  M_self=-2.0*sqrt(sigma/M_PI);

  E_real=M_real/(4*M_PI*EPS0);
  E_self=M_self/(4*M_PI*EPS0);

  ii=1+sqrt(4*sigma)/sqrt(vmod2(c->recip[0]));
  jj=1+sqrt(4*sigma)/sqrt(vmod2(c->recip[1]));
  kk=1+sqrt(4*sigma)/sqrt(vmod2(c->recip[2]));
  if (debug>1)
    fprintf(stderr,"Rec space Ewald grid=%dx%dx%d\n",2*ii+1,2*jj+1,2*kk+1);
  M_recip=0.0;
  for(i=-ii;i<=ii;i++){
    for(j=-jj;j<=jj;j++){
      for(k=-kk;k<=kk;k++){
	if ((i==0)&&(j==0)&&(k==0)) continue;
	for(m=0;m<3;m++) disp[m]=c->recip[0][m]*i+
			   c->recip[1][m]*j+c->recip[2][m]*k;
	dist=4*M_PI*M_PI*vmod2(disp);
	M_recip+=exp(-0.25*dist/sigma)/dist;
      }
    }
  }

  M_recip*=4*M_PI/c->vol;
  E_recip=M_recip/(4*M_PI*EPS0);

  M_charged=-M_PI/(c->vol*sigma);
  E_charged=M_charged/(4*M_PI*EPS0);
  
  if(debug>2){
	  /* For comparison with pymatgen */
    fprintf(stderr,"Volume:      %f\n",c->vol);
    fprintf(stderr,"Real space:  %f\n",E_real);
    fprintf(stderr,"Recip space: %f\n",E_recip);
    fprintf(stderr,"Self:        %f\n",E_self);
    fprintf(stderr,"Charged:     %f\n",E_charged);
    fprintf(stderr,"Total:       %f\n",E_real+E_recip+E_self+E_charged);
  }

  M=M_real+M_recip+M_self+M_charged;
  M=-M*sqrt(vmod2(c->basis[0]));

  if (debug) fprintf(stderr,"Madelung constant of lattice: %f\n",M);

  return M;
  
}

double quadrupole(struct unit_cell *c, struct contents *m,
                  struct grid *g, double *ctr){
  double q,Q_e,Q_i,vec[3],*ptr;
  int i,ii[3],j;
  double rvec[3],disp2,scale;

  scale=c->vol/(g->size[0]*g->size[1]*g->size[2]);
  Q_e=Q_i=0;
  q=0;

  ptr=g->data;
  if (!ptr) return 0;

  for(i=0;i<m->n;i++){
    for(j=0;j<3;j++){
      rvec[j]=m->atoms[i].frac[j]-ctr[j];
      rvec[j]=fmod(rvec[j]+0.5,1.0);
      if (rvec[j]<0) rvec[j]+=1;
      rvec[j]-=0.5;
    }
    for(j=0;j<3;j++)
      vec[j]=rvec[0]*c->basis[0][j]+
        rvec[1]*c->basis[1][j]+rvec[2]*c->basis[2][j];
    disp2=vmod2(vec);
    Q_i+=m->atoms[i].chg*disp2;
  }

  if (debug) fprintf(stderr,"Ionic quadrupole: Q %f eA^2\n",Q_i);

  
  for(ii[0]=0;ii[0]<g->size[0];ii[0]++){
    for(ii[1]=0;ii[1]<g->size[1];ii[1]++){
      for(ii[2]=0;ii[2]<g->size[2];ii[2]++){
        for(i=0;i<3;i++){
          rvec[i]=(double)ii[i]/g->size[i]-ctr[i]; /* BUG!!! */
      /* force disp to range 0.5<=disp<=0.5 */
          rvec[i]=fmod(rvec[i]+0.5,1.0);
          if (rvec[i]<0) rvec[i]+=1;
          rvec[i]-=0.5;
        }
        for(i=0;i<3;i++)
          vec[i]=rvec[0]*c->basis[0][i]+
            rvec[1]*c->basis[1][i]+rvec[2]*c->basis[2][i];
        disp2=vmod2(vec);
        q-=scale*(*ptr);
        Q_e-=scale*(*ptr)*disp2;
        ptr++;
      }
    }
  }

  if (debug) fprintf(stderr,"Electric quadrupole: q=%fe, Q %f eA^2\n",q,Q_e);

  if (debug) fprintf(stderr,"Total quadrupole: %f eA^2\n",Q_e+Q_i);
  return Q_e+Q_i;
  
}

double quadrupole_ii(struct unit_cell *c, struct contents *m,
		     struct grid *g, double *ctr, int dir){
  double q,Q_e,Q_i,vec[3],*ptr;
  int i,ii[3],j;
  double rvec[3],disp2,scale;

  scale=c->vol/(g->size[0]*g->size[1]*g->size[2]);
  Q_e=Q_i=0;
  q=0;

  ptr=g->data;
  if (!ptr) return 0;

  for(i=0;i<m->n;i++){
    for(j=0;j<3;j++){
      rvec[j]=m->atoms[i].frac[j]-ctr[j];
      rvec[j]=fmod(rvec[j]+0.5,1.0);
      if (rvec[j]<0) rvec[j]+=1;
      rvec[j]-=0.5;
    }
    for(j=0;j<3;j++)
      vec[j]=rvec[0]*c->basis[0][j]+
        rvec[1]*c->basis[1][j]+rvec[2]*c->basis[2][j];
    disp2=vec[dir]*vec[dir];
    Q_i+=m->atoms[i].chg*disp2;
  }

  if (debug) fprintf(stderr,"Ionic quadrupole_%c%c: Q %f eA^2\n",
		     dir+'a',dir+'a',Q_i);

  
  for(ii[0]=0;ii[0]<g->size[0];ii[0]++){
    for(ii[1]=0;ii[1]<g->size[1];ii[1]++){
      for(ii[2]=0;ii[2]<g->size[2];ii[2]++){
        for(i=0;i<3;i++){
          rvec[i]=(double)ii[i]/g->size[i]-ctr[i]; /* BUG!!! */
      /* force disp to range 0.5<=disp<=0.5 */
          rvec[i]=fmod(rvec[i]+0.5,1.0);
          if (rvec[i]<0) rvec[i]+=1;
          rvec[i]-=0.5;
        }
        for(i=0;i<3;i++)
          vec[i]=rvec[0]*c->basis[0][i]+
            rvec[1]*c->basis[1][i]+rvec[2]*c->basis[2][i];
        disp2=vec[dir]*vec[dir];
        q-=scale*(*ptr);
        Q_e-=scale*(*ptr)*disp2;
        ptr++;
      }
    }
  }

  if (debug) fprintf(stderr,"Electric quadrupole_%c%c: q=%fe, Q %f eA^2\n",
		     dir+'a',dir+'a',q,Q_e);

  if (debug) fprintf(stderr,"Total quadrupole_%c%c: %f eA^2\n",
		     dir+'a',dir+'a',Q_e+Q_i);
  return Q_e+Q_i;
  
}

void charge_corr(struct unit_cell *c, struct contents *m,
                 struct grid *g, struct es *elect){
  double charge,alpha,energy,quad,ctr[3];
  double i_charge,e_charge;
  double abc[6],a,dpole[3],dpole_frac[3];
  int i,j,n_grid_points,dir;

  if ((!g)||(!g->data)){
    fprintf(stderr,
	    "Unable to correct for charge as no electron density read\n");
    return;
  }

  dict_strcat(m->dict,"dipole_no_charge_warn","");
  
  if (elect->charge) charge=*elect->charge;
  else{ /* Need to calculate charge in cell */
    i_charge=0;
    for(i=0;i<m->n;i++) i_charge+=m->atoms[i].chg;
    if (!g->data){
      fprintf(stderr,"Cannot calculate charge correction as net charge "
              "not available\nDid you mean to read a density too?");
      return;
    }
    e_charge=0;
    n_grid_points=g->size[0]*g->size[1]*g->size[2];
    for(i=0;i<n_grid_points;i++) e_charge+=g->data[i];
    e_charge*=c->vol/n_grid_points;
    if (debug>1) fprintf(stderr,"Ionic charge: %f\nElectronic charge: %f\n",
                         i_charge,e_charge);
    charge=i_charge-e_charge;
  }

  if (*elect->dip_corr_dir=='m'){ /* 3D to 0D correction */
  
    alpha=madelung(c);

    energy=charge*charge*alpha/(8*M_PI*EPS0*sqrt(vmod2(c->basis[0])));

    if (debug) fprintf(stderr,"Total charge: %f\n",charge);
    fprintf(stderr,"Madelung energy correction: %12.6f eV\n",energy);
    if (elect->energy){
      fprintf(stderr,"Original energy:            %12.6f eV\n",*elect->energy);
      fprintf(stderr,"Corrected energy:           %12.6f eV\n",
	      *elect->energy+energy);
    }

    ctr[0]=ctr[1]=ctr[2]=0.5;
    cart2abc(c,NULL,abc,NULL,0);
    dipole_calc(c,m,g,ctr,dpole);
    for(i=0;i<3;i++){
      if (fabs(dpole[i])>fabs(charge*abc[i])){
	fprintf(stderr,"Charge too small / dipole too large for next term\n");
	return;
      }
      dpole_frac[i]=0;
      for(j=0;j<3;j++)
	dpole_frac[i]+=dpole[i]*c->recip[i][j];
      ctr[i]=ctr[i]+dpole_frac[i]/charge;
    }

    quad=quadrupole(c,m,g,ctr);

    fprintf(stderr,"Quadrupole correction:      %12.6f eV\n",
	    charge*quad/(6*EPS0*c->vol));
    fprintf(stderr,"Final corrected energy:     %12.6f eV\n",
	    *elect->energy+energy-charge*quad/(6*EPS0*c->vol));
  }
  else{  /* 3D to 2D slab correction */
    cart2abc(c,NULL,abc,NULL,0);
    dir=*elect->dip_corr_dir-'a';
    for(i=0;i<3;i++){
      if ((!aeq(abc[3+i],90))&&(dir!=i)){
	  fprintf(stderr,"Warning: charge correction axis not perpendicular "
		  "to other axes\n");
      }
    }
    a=abc[(dir+1)%3]*abc[(dir+2)%3]*sin(abc[3+dir]*M_PI/180);

    energy=-charge*charge*abc[dir]/(24*EPS0*a);

    if (debug) fprintf(stderr,"Total charge: %f\n",charge);
    fprintf(stderr,"Charged slab correction:  %12.6f eV\n",energy);
    if (elect->energy){
      fprintf(stderr,"Original energy:          %12.6f eV\n",*elect->energy);
      fprintf(stderr,"Corrected energy:         %12.6f eV\n",
	      *elect->energy+energy);

      ctr[0]=ctr[1]=ctr[2]=0.5;
      dipole_calc(c,m,g,ctr,dpole);
      if (fabs(dpole[dir])>fabs(charge*abc[dir])){
	fprintf(stderr,"Charge too small / dipole too large for next term\n");
	return;
      }
      ctr[dir]=ctr[dir]+dpole[dir]/(charge*abc[dir]);

      if (debug) fprintf(stderr,"Centre for quadrupole: (%f,%f,%f)\n",
			 ctr[0],ctr[1],ctr[2]);
      
      quad=quadrupole_ii(c,m,g,ctr,dir);
      quad=-quad*charge/(2*EPS0*c->vol);
      fprintf(stderr,"Quadrupole correction:    %12.6f eV\n",quad);
      fprintf(stderr,"Final corrected energy:   %12.6f eV\n",
	      *elect->energy+energy+quad);
    
    }
  }
  
}
