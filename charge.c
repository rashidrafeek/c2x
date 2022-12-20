/* Calculate Madelung constant and leading term in correction for
 * charged cell */


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "c2xsf.h"

// #undef EPS0
// #define EPS0 (1/180.952701)


double madelung(struct unit_cell *c){
  double M, M_real,M_recip,M_self,M_charged;
  double E_real,E_recip,E_self,E_charged;
  double dist,disp[3];
  double sigma;
  int i,j,k,ii,jj,kk,m;

  sigma=M_PI*pow(0.7/(c->vol*c->vol),1.0/3.0);
  fprintf(stderr,"sigma=%f\n",sigma);

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

void charge_corr(struct unit_cell *c, struct contents *m,
                 struct grid *g, struct es *elect){
  double charge,alpha,energy,quad,ctr[3];
  double i_charge,e_charge;
  int i,n_grid_points;

  alpha=madelung(c);

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

  energy=charge*charge*alpha/(8*M_PI*EPS0*sqrt(vmod2(c->basis[0])));

  if (debug) fprintf(stderr,"Total charge: %f\n",charge);
  fprintf(stderr,"Madelung energy correction: %.6f eV\n",energy);
  if (elect->energy){
    fprintf(stderr,"Original energy:            %.6f eV\n",*elect->energy);
    fprintf(stderr,"Corrected energy:           %.6f eV\n",
            *elect->energy+energy);
  }

  ctr[0]=ctr[1]=ctr[2]=0.5;

  quad=quadrupole(c,m,g,ctr);

  fprintf(stderr,"Fiddle=%f eV\n",charge*quad/(6*EPS0*c->vol));
  fprintf(stderr,"New corrected: %.6f eV\n",
          *elect->energy+energy-charge*quad/(6*EPS0*c->vol));

  
}
