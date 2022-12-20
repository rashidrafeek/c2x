/* Various functions for processing bands */


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
#include<stdlib.h> /* malloc */
#include<math.h>

#include "c2xsf.h"

void gcoeff_write(double *psi, int *pwgrid, int nplwv, int fft[3],
		  int gamma, int goff[3],
		  struct unit_cell *c, struct contents *m, struct kpts *kpt,
		  int ikpt, int isppol, int nb, double *eval, double occ,
		  struct es *e);

void wavecar_write(double *psi, int *pwgrid, int nplwv, int fft[3],
                   int gamma, int goff[3],
                   struct unit_cell *c, struct contents *m, struct kpts *kp,
                   int ikpt, int isppol, int nb,
                   double *eval, double occ, struct es *e);

void band_store(struct grid **gp, double *dptr, double occ, double wkpt,
		int nspr, int ns, int k, int b, struct es *elect,
		struct contents *m, int fft[3]);
void band2real(double *psi, double *out, int nfft[3], double kpoint[3]);

double *band2grid(double *dptr, int fft[3], int *pwgrid, int npw, int gamma){
  double *psi_g;
  int i,offset,n0,n1,n2,nfftpts,goff[3];

  nfftpts=fft[0]*fft[1]*fft[2];
  
  for(i=0;i<3;i++) goff[i]=0;
  if (gamma){
    gamma++;
    if (gamma&1) goff[0]=1;
    if (gamma>5) goff[1]=1;
    if ((gamma-2)&2) goff[2]=1;
  }
  
  psi_g=malloc(2*nfftpts*sizeof(double));
  if (!psi_g) error_exit("Malloc error for psi");

  for(i=0;i<2*nfftpts;i++) psi_g[i]=0;
  for(i=0;i<npw;i++){
    offset=pwgrid[3*i+2]+fft[2]*(pwgrid[3*i+1]+
                                 fft[1]*pwgrid[3*i]);
    if ((offset<0)||(offset>nfftpts)){
      fprintf(stderr,"Impossible offset in band2grid off=%d i=%d\n",
              offset,i);
      exit(1);
    }
    psi_g[2*offset]=dptr[2*i];
    psi_g[2*offset+1]=dptr[2*i+1];
  }

  if (gamma>1){ /* construct psi(-k)=conjg(psi(k)) */
    if (debug>1) fprintf(stderr,"Gamma point storage type %d\n",gamma);
    for(i=0;i<npw;i++){
      if ((gamma==2)&&(pwgrid[3*i]==0)&&
          (pwgrid[3*i+1]==0)&&(pwgrid[3*i+2]==0)) continue;
      n0=fft[2]-pwgrid[3*i+2]-goff[2];
      if (n0==fft[2]) n0=0;
      n1=fft[1]-pwgrid[3*i+1]-goff[1];
      if (n1==fft[1]) n1=0;
      n2=fft[0]-pwgrid[3*i]-goff[0];
      if (n2==fft[0]) n2=0;
      offset=n0+fft[2]*(n1+fft[1]*n2);
      if ((offset<0)||(offset>nfftpts)){
        fprintf(stderr,
                "Impossible -offset in band2grid off=%d i=%d\n",
                offset,i);
        exit(1);
      }
      psi_g[2*offset]=dptr[2*i];
      psi_g[2*offset+1]=-dptr[2*i+1];
    }
  }

  return psi_g;
  
}

/* Use abinit's definition of gamma, less one.
 *  gamma=0  -- not gamma
 *  gamma=1  -- gamma
 *  gamma>1  -- abinit's istwfk = gamma+1
 *
 * ikpt, ispinor, isppol and nb all start from 0, not 1
 *
 * isppol=-1 means both components of a spinor wavefunction being presented
 *   for gcoeff or wavecar output
 *
 */

void band_process(double *dptr, int fft[3], int *pwgrid, int npw, int gamma,
		  struct unit_cell *c,
		  struct grid **gp, struct es *elect, struct kpts *kp,
		  struct contents *m, int ikpt, int ispinor, int isppol,
		  int nb, int *i_grid){
  int i,off;
  int nfftpts,ffft[3],nfft[3],goff[3];
#if 0
  int n0,n1,n2,offset;
#endif
  double scale,occ,*psi,*kpoint,eval[2];

  off=elect->nbands*elect->nbspins*ikpt+isppol*elect->nbands;
  if (elect->occ)
    occ=elect->occ[off+nb];
  else
    occ=1;

  for(i=0;i<3;i++) nfft[i]=fft[i];
  nfftpts=fft[0]*fft[1]*fft[2];
  
  for(i=0;i<3;i++) goff[i]=0;
  if (gamma){
    gamma++;
    if (gamma&1) goff[0]=1;
    if (gamma>5) goff[1]=1;
    if ((gamma-2)&2) goff[2]=1;
    gamma--;
  }
  
  scale=1.0;
  kpoint=kp->kpts[ikpt].frac;

  if (flags&GCOEFF){
    if (elect->eval)
      eval[0]=elect->eval[off+nb];
    else
      eval[0]=1;
    eval[1]=0;
    if (dict_get(m->dict,"wavecar_output"))
      wavecar_write(dptr,pwgrid,npw,fft,gamma,goff,c,m,kp,ikpt,isppol,nb+1,
                    eval,occ,elect);
    else
      gcoeff_write(dptr,pwgrid,npw,fft,gamma,goff,c,m,kp,ikpt,isppol,nb+1,
                    eval,occ,elect);
    return;
  }

  if (isppol==-1){
    fprintf(stderr,"Unable to process spinor wavefunction\n");
    return;
  }
  
  psi=band2grid(dptr,fft,pwgrid,npw,gamma);
#if 0
  for(i=0;i<3;i++) goff[i]=0;
  if (gamma){
    gamma++;
    if (gamma&1) goff[0]=1;
    if (gamma>5) goff[1]=1;
    if ((gamma-2)&2) goff[2]=1;
  }

  psi=malloc(16*nfftpts);
  if (!psi) error_exit("Malloc error for psi");
  for(i=0;i<2*nfftpts;i++) psi[i]=0;
  for(i=0;i<npw;i++){
    offset=pwgrid[3*i+2]+fft[2]*(pwgrid[3*i+1]+
				 fft[1]*pwgrid[3*i]);
    if ((offset<0)||(offset>nfftpts)){
      fprintf(stderr,"Impossible offset in wave_read off=%d i=%d\n",
	      offset,i);
      exit(1);
    }
    psi[2*offset]=dptr[2*i];
    psi[2*offset+1]=dptr[2*i+1];
  }
  if (gamma>1){ /* construct psi(-k)=conjg(psi(k)) */
    if (debug>1) fprintf(stderr,"Gamma point storage type %d\n",gamma);
    for(i=0;i<npw;i++){
      if ((gamma==2)&&(pwgrid[3*i]==0)&&
	  (pwgrid[3*i+1]==0)&&(pwgrid[3*i+2]==0)) continue;
      n0=fft[2]-pwgrid[3*i+2]-goff[2];
      if (n0==fft[2]) n0=0;
      n1=fft[1]-pwgrid[3*i+1]-goff[1];
      if (n1==fft[1]) n1=0;
      n2=fft[0]-pwgrid[3*i]-goff[0];
      if (n2==fft[0]) n2=0;
      offset=n0+fft[2]*(n1+fft[1]*n2);
      if ((offset<0)||(offset>nfftpts)){
	fprintf(stderr,
		"Impossible -offset in wave_read off=%d i=%d\n",
		offset,i);
	exit(1);
      }
      psi[2*offset]=dptr[2*i];
      psi[2*offset+1]=-dptr[2*i+1];
    }
  }
#endif

  if ((aeq(kpoint[0],0)||aeq(fabs(kpoint[0]),0.5))&&
      (aeq(kpoint[1],0)||aeq(fabs(kpoint[1]),0.5))&&
      (aeq(kpoint[2],0)||aeq(fabs(kpoint[2]),0.5))&&
      (flags&BANDPARITY)) inv_parity(psi,fft,nb+1,kpoint);

  /* Was the parity all we were requested to report? */
  if (!(flags&BANDS))
    free(psi);
  else{

    /* Padding */
            
    if (i_grid){
      for(i=0;i<3;i++) nfft[i]=i_grid[i];
      if(debug>1)
	fprintf(stderr,"Padding wavefunction onto %dx%dx%d grid\n",
		nfft[0],nfft[1],nfft[2]);
      if ((fft[0]==nfft[0])&&(fft[1]==nfft[1])&&(fft[2]==nfft[2])){
	if (debug>1)
	  fprintf(stderr,"Skipping null padding operation\n");
      }
      else{
	pad_recip(psi,fft,&dptr,nfft);
	nfftpts=nfft[0]*nfft[1]*nfft[2];
	free(psi);
	psi=dptr;
      }
    }

          
    ffft[0]=nfft[2];
    ffft[1]=nfft[1];
    ffft[2]=nfft[0];
    fft3d(psi,ffft,1);
    dptr=malloc(nfftpts*sizeof(double));
    if(!dptr) error_exit("Malloc error for grid data");
    band2real(psi,dptr,nfft,kp->kpts[ikpt].frac);
    free(psi);

    /* Do we need to rescale? */
    if (((flags&RAW)==0)&&((flags&BANDPHASE)==0)){ /* Yes */
      if (flags&BANDDEN) scale=1/c->vol;
      else scale=1/sqrt(c->vol);
      if (debug>2) fprintf(stderr,"Scaling wavefun by %f\n",scale);
      for(i=0;i<nfftpts;i++) dptr[i]*=scale;
    }


    
    band_store(gp,dptr,occ,kp->kpts[ikpt].wt,
	       ispinor,isppol,ikpt+1,nb+1,elect,m,nfft);

  }
}

  
void band2real(double *psi, double *out, int nfft[3], double kpoint[3]){
  double phase_r,phase_r2,phase_i,phase_i2,phi,dtmp;
  double min,max,sum;
  int i,ii,jj,kk,ind,nfft_pts;

  nfft_pts=nfft[0]*nfft[1]*nfft[2];
  
  if ((!(flags&BANDDEN))&&
      ((kpoint[0]!=0)||(kpoint[1]!=0)||(kpoint[2]!=0))){ /* want psi,
							    but not at gamma! */
    if (debug)
      fprintf(stderr,"unwinding psi for non-gamma k-point...\n");
    for(ii=0;ii<nfft[0];ii++){
      for(jj=0;jj<nfft[1];jj++){
	for(kk=0;kk<nfft[2];kk++){
	  phi=2*M_PI*((ii*kpoint[0])/nfft[0]+
		      (jj*kpoint[1])/nfft[1]+
		      (kk*kpoint[2])/nfft[2]);
	  phase_r=cos(phi);
	  phase_i=sin(phi);
	  ind=2*(kk+nfft[2]*(jj+ii*nfft[1]));
	  dtmp=psi[ind];
	  psi[ind]=phase_r*psi[ind]-phase_i*psi[ind+1];
	  psi[ind+1]=phase_r*psi[ind+1]+phase_i*dtmp;
	}
      }
    }
  }
  phase_r=phase_i=phase_r2=phase_i2=0;
  for(i=0;i<nfft_pts;i++){
    if (psi[2*i]>0){
      phase_r+=psi[2*i];
      phase_i-=psi[2*i+1];
    }else{
      phase_r2-=psi[2*i];
      phase_i2+=psi[2*i+1];
    }
  } 
  phase_r+=phase_r2;
  phase_i+=phase_i2;
  dtmp=sqrt(phase_r*phase_r+phase_i*phase_i);
  phase_r/=dtmp;
  phase_i/=dtmp;
  ii=0;
  max=-1e300;
  min=1e300;
  sum=0;
  for (i=0;i<nfft_pts;i++){
    if (flags&BANDPHASE){
      out[i]=atan2(psi[2*i+1],psi[2*i]);
    }
    else if (flags&BANDREAL){
      out[i]=psi[2*i];
    }
    else if (flags&BANDIMAG){
      out[i]=psi[2*i+1];
    }
    else
      if (flags&BANDDEN)
	out[i]=psi[2*i]*psi[2*i]+psi[2*i+1]*psi[2*i+1];
      else{
	out[i]=psi[2*i]*phase_r-psi[2*i+1]*phase_i;
	dtmp=psi[2*i]*phase_i+psi[2*i+1]*phase_r;
	if((fabs(dtmp)>.05))ii++;
      }
    sum+=out[i];
    if(out[i]<min) min=out[i];
    if(out[i]>max) max=out[i];
  }
  if (debug>2) fprintf(stderr,"Min=%g Max=%g Sum=%g\n",min,max,sum);
  if ((debug>1)&&(ii>0)) fprintf(stderr,"Warning: %d components with "
			      " imaginary part >0.05\n",ii);

}

#define CBUFF 100
void band_store(struct grid **gp, double *dptr, double occ, double wkpt,
		int nspr, int ns, int k, int b, struct es *elect,
		struct contents *m, int fft[3]){
  double w;
  int i,nfft_pts;
  char cbuff[CBUFF+1];
  struct grid *g;

  g=*gp;
  nfft_pts=fft[0]*fft[1]*fft[2];
  w=1;

  /* Do we need to weight? */
  if ((flags&OCC_WEIGHT)||(flags&K_WEIGHT)){
    if (flags&OCC_WEIGHT) w=occ;
    if (flags&K_WEIGHT) w*=wkpt;
    /* If we want densities, and we do not have spins, each
     * band is doubly occupied */
    if ((elect->nspins==1)&&(elect->nspinors==1)&&
	(flags&BANDDEN)&&(flags&OCC_WEIGHT))
      w*=2;
    if ((w!=1)&&(!(flags&BANDDEN))) w=sqrt(w);
    if (debug)
      fprintf(stderr,"Using weight %f for ns=%d k=%d band=%d\n",
	      w,ns,k,b);
    if (debug>1)
      fprintf(stderr,"  kpt weight %f occupancy %f\n",wkpt,occ);
    
    if (w!=1)
      for(i=0;i<nfft_pts;i++) dptr[i]*=w;
  }

  if (debug) {
    fprintf(stderr,"Processing band %d kpoint %d",b,k);
    if (elect->nspinors==2)
      fprintf(stderr," spin %d\n",nspr);
    else if (elect->nspins==2)
      fprintf(stderr," spin %d\n",ns);
    else
      fprintf(stderr,"\n");
  }
  
  if (!(flags&ACCUMULATE)){
    g->data=dptr;
    for(i=0;i<3;i++) g->size[i]=fft[i];
    g->name=malloc(40);
    if (!g->name) error_exit("malloc error for name");
    if (elect->nspinors==2)
      sprintf(g->name,"band_vs%d_k%d_b%d",nspr,k,b);
    else if (elect->nspins==2)
      sprintf(g->name,"band_s%d_k%d_b%d",ns,k,b);
    else
      sprintf(g->name,"band_k%d_b%d",k,b);
    g->next=malloc(sizeof(struct grid));
    if (!g->next) error_exit("mallock error for grid");
    g=g->next;
    g->data=NULL;
    g->next=NULL;
    g->origin_abs=NULL;
    if ((flags&OCC_WEIGHT)||(flags&K_WEIGHT)){
      snprintf(cbuff,CBUFF,
	       "Weight %f used for spin %d kpt %d band %d",
	       w,(elect->nspinors==2)?nspr:ns,k,b);
      if (m) add_cmt(m->comment,cbuff);
    }
  }else{  /* Are accumulating */
    if (!g->data){  /* This is the first set */
      g->data=dptr;
      for(i=0;i<3;i++) g->size[i]=fft[i];
      g->name=malloc(40);
      if (!g->name) error_exit("mallock error for name");
      sprintf(g->name,"bands"); /* Don't move to a new grid */
    }else{
      for(i=0;i<nfft_pts;i++) g->data[i]+=dptr[i];
      free(dptr);
    }
  }
  *gp=g;
}


/* Find maximum |g|^2 in a pwgrid of gvectors */
double g2max(double recip[3][3], int *pwgrid, int nplwv, int fft[3],
             double *kpt){
  int i,j;
  double g2,g2m,gv[3],x[3];

  g2m=0;
  for(i=0;i<nplwv;i++){
    for(j=0;j<3;j++){
      x[j]=pwgrid[3*i+j]+kpt[j];
      if (x[j]>fft[j]/2) x[j]-=fft[j];
    }
    for(j=0;j<3;j++)
      gv[j]=x[0]*recip[0][j]+x[1]*recip[1][j]+x[2]*recip[2][j];
    g2=vmod2(gv);
    g2m=max(g2m,g2);
  }

  return g2m;

}
