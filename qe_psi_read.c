/* Read wavefunctions from PWscf assuming that they are
 * in the form of one file per k-point per spin, all in the same directory.
 * Non colinear spins not supported.
 */

/* Copyright (c) 2019 MJ Rutter 
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "c2xsf.h"

void fft3d(double *c, int *ngptar, int dir);
int inrange(int x, char *range); /* This and next two from check_read.c */
void band2real(double *psi, double *out, int nfft[3], double kpoint[3]);

void kcart2frac(double xk[3], double kpt[3], double recip[3][3]){
  int i;
  struct unit_cell ctmp;
  struct atom atmp;

  if (debug>1)
    fprintf(stderr,"Kpt given in abs terms as (%f,%f,%f)\n",
	    xk[0],xk[1],xk[2]);
  ctmp.basis=recip;
  real2rec(&ctmp);
  for(i=0;i<3;i++) atmp.abs[i]=xk[i];
  addfrac(&atmp,1,ctmp.recip);
  for(i=0;i<3;i++) kpt[i]=atmp.frac[i];
  if (debug>1)
    fprintf(stderr,"Kpt converted to frac coords: (%f,%f,%f)\n",
	    kpt[0],kpt[1],kpt[2]);
  
}

void qe_psi_read(char *dir, char *prefix, struct unit_cell *c,
		 struct contents *m,
                 struct kpts *k, struct symmetry *s, struct grid *g,
                 struct es *elect, int fft[3], int *i_grid){
  int i,j,ns,k2,nb,offset,tmp,ii;
  int n0,n1,n2;
  char *file,*sspin;
  int ik,ispin,gamma,ngw,igwx,npol,nbnd,ngpts;
  int *pwgrid,ffft[3],nfft[3];
  double xk[3],scalef,b[9],kpoint[3],*wtmp,*psi,*dptr;
  double occ,conv,lat_tmp[3][3];
  FILE *infile;

  if (i_grid){
    for(i=0;i<3;i++)
      if (i_grid[i]==0)
        nfft[i]=fft[i];
      else
        nfft[i]=i_grid[i];
  }
  else
    for(i=0;i<3;i++) nfft[i]=fft[i];
  
  if (dir==NULL) dir=".";
  file=malloc(strlen(dir)+30+(prefix?strlen(prefix):0));

  if ((elect->nspins!=1)&&(elect->nspins!=2)){
    fprintf(stderr,"Error: cannot read wavefunction for nspins=%d\n",
	    elect->nspins);
    return;
  }
  
  for(ns=0;ns<elect->nspins;ns++){
    
    if (elect->nspins==1)
      sspin="";
    else{
      if (!inrange(ns,elect->spin_range)) continue;
      if (ns==0) sspin="up";
      else sspin="dw";
    }
    
    for(k2=1;k2<=k->n;k2++){

      if (!inrange(k2,elect->kpt_range)) continue;
  
      sprintf(file,"%s/wfc%s%d.dat",dir,sspin,k2);

      infile=fopen(file,"r");
      if (!infile){
	if (prefix){
	  sprintf(file,"%s/%s.save/wfc%s%d.dat",dir,prefix,sspin,k2);
	  infile=fopen(file,"r");
	}
	if (!infile){
	  fprintf(stderr,"Failed to open %s/wfc%s%d.dat or %s for reading\n",
		  dir,sspin,k2,file);
	  return;
	}
      }

      if (debug>1) fprintf(stderr,"Reading %s\n",file);
      tmp=0;
      fread(&tmp,4,1,infile);
      if (tmp!=44){
	fprintf(stderr,"Unexpected first record in %s\n",file);
	return;
      }

      fread(&ik,4,1,infile);
      fread(xk,8,3,infile);
      fread(&ispin,4,1,infile);
      fread(&gamma,4,1,infile);
      fread(&scalef,8,1,infile);
      fread(&tmp,4,1,infile);

      if (ispin!=ns+1)
	fprintf(stderr,"Warning: unexpected value of ispin in %s. "
		"Expected %d, read %d\n",file,ns+1,ispin);

      
      fread(&tmp,4,1,infile);
      if (tmp!=16){
	fprintf(stderr,"Unexpected second record in %s\n",file);
	return;
      }

      fread(&ngw,4,1,infile);
      fread(&igwx,4,1,infile);
      if (debug>2) fprintf(stderr,"igwx=%d\n",igwx);
      fread(&npol,4,1,infile);
      elect->nspinors=npol;
      if (npol!=1) fprintf(stderr,"Warning: c2x's ability to read "
                           "non collinear spin wavefunctions untested\n");
      fread(&nbnd,4,1,infile);
      fread(&tmp,4,1,infile);

      fread(&tmp,4,1,infile);
      if (tmp!=72){
	fprintf(stderr,"Unexpected third record in %s\n",file);
	return;
      }

      fread(b,8,9,infile);
      fread(&tmp,4,1,infile);
    
      for(i=0;i<9;i++) lat_tmp[i/3][i%3]=b[i];
      kcart2frac(xk,kpoint,lat_tmp);
      if ((!aeq(kpoint[0],k->kpts[k2-1].frac[0]))||
          (!aeq(kpoint[1],k->kpts[k2-1].frac[1]))||
          (!aeq(kpoint[2],k->kpts[k2-1].frac[2]))){
        fprintf(stderr,"Warning: expected kpoint %d to be (%.6f,%.6f,%.6f)\n"
                "        found (%.6f,%.6f,%.6f)\n",k2,k->kpts[k2-1].frac[0],
                k->kpts[k2-1].frac[1],k->kpts[k2-1].frac[2],kpoint[0],
                kpoint[1],kpoint[2]);
      }
      
        
  
  
      pwgrid=malloc(igwx*3*sizeof(int));
      if (!pwgrid) error_exit("Malloc error for pwgrid");

      fread(&tmp,4,1,infile);
      if (tmp!=igwx*3*sizeof(int)) error_exit("Unexpected pwgrid record");
      fread(pwgrid,4,igwx*3,infile);
      fread(&tmp,4,1,infile);

      wtmp=malloc(igwx*2*npol*sizeof(double));
      if (!pwgrid) error_exit("Malloc error for wtmp");

      if (debug>1)
	fprintf(stderr,"FFT grid passed: %dx%dx%d\n",fft[0],fft[1],fft[2]);

      if (debug>2){
	fprintf(stderr,"FFT grid components in pwgrid: ");
	tmp=0;
	for(i=0;i<igwx;i++)
	  if (tmp>pwgrid[3*i]) tmp=pwgrid[3*i];
	fprintf(stderr,"%d:",tmp);
	tmp=0;
	for(i=0;i<igwx;i++)
	  if (tmp<pwgrid[3*i]) tmp=pwgrid[3*i];
	fprintf(stderr,"%dx",tmp);
	tmp=0;
	for(i=0;i<igwx;i++)
	  if (tmp>pwgrid[3*i+1]) tmp=pwgrid[3*i+1];
	fprintf(stderr,"%d:",tmp);
	tmp=0;
	for(i=0;i<igwx;i++)
	  if (tmp<pwgrid[3*i+1]) tmp=pwgrid[3*i+1];
	fprintf(stderr,"%dx",tmp);
	tmp=0;
	for(i=0;i<igwx;i++)
	  if (tmp>pwgrid[3*i+2]) tmp=pwgrid[3*i+2];
	fprintf(stderr,"%d:",tmp);
	tmp=0;
	for(i=0;i<igwx;i++)
	  if (tmp<pwgrid[3*i+2]) tmp=pwgrid[3*i+2];
	fprintf(stderr,"%d\n",tmp);
      }
    
      for(i=0;i<3;i++)
	for(j=0;j<igwx;j++)
	  if (pwgrid[3*j+i]<0) pwgrid[3*j+i]+=fft[i];
    
      ngpts=fft[0]*fft[1]*fft[2];

      psi=malloc(ngpts*2*sizeof(double));
      if (!psi) error_exit("Malloc error for psi");
    
      for(nb=1;nb<=nbnd;nb++){
	if (debug>2) fprintf(stderr,"Start band %d\n",nb);
	fread(&tmp,4,1,infile);
	if (tmp!=16*npol*igwx) error_exit("Error parsing wavefunction band");
	if (inrange(nb,elect->band_range)){
	  if (debug>2) fprintf(stderr,"Reading band %d\n",nb);
	  fread(wtmp,16,npol*igwx,infile);
	  fread(&tmp,4,1,infile);

          for(ii=0;ii<npol;ii++){
            if ((npol>1)&&(!inrange(ii+1,elect->spin_range))) continue;
            for(i=0;i<2*ngpts;i++) psi[i]=0;
            for(i=0;i<igwx;i++){
              offset=pwgrid[3*i+2]+fft[2]*(pwgrid[3*i+1]+
                                           fft[1]*pwgrid[3*i]);
              if ((offset<0)||(offset>ngpts)){
                fprintf(stderr,"Impossible offset in wave_read off=%d i=%d\n",
                        offset,i);
                exit(1);
              }
              
              psi[2*offset]=wtmp[2*ii*igwx+2*i];
              psi[2*offset+1]=wtmp[2*ii*igwx+2*i+1];
            }
            if (gamma){ /* construct psi(-k)=conjg(psi(k)) */
              for(i=0;i<igwx;i++){
                if ((pwgrid[3*i]==0)&&(pwgrid[3*i+1]==0)&&(pwgrid[3*i+2]==0))
                  continue;
                n0=fft[2]-pwgrid[3*i+2];
                if (n0==fft[2]) n0=0;
                n1=fft[1]-pwgrid[3*i+1];
                if (n1==fft[1]) n1=0;
                n2=fft[0]-pwgrid[3*i];
                if (n2==fft[0]) n2=0;
                offset=n0+fft[2]*(n1+fft[1]*n2);
                if ((offset<0)||(offset>ngpts)){
                  fprintf(stderr,
                          "Impossible -offset in wave_read off=%d i=%d\n",
                          offset,i);
                  exit(1);
                }
                psi[2*offset]=wtmp[ii*igwx+2*i];
                psi[2*offset+1]=-wtmp[ii*igwx+2*i+1];
              }
            }
            if (debug>2) fprintf(stderr,"Before FFT g=0 component is %g+%gi\n",
                                 psi[0],psi[1]);

            if (((kpoint[0]==0)||aeq(fabs(kpoint[0]),0.5))&&
                ((kpoint[1]==0)||aeq(fabs(kpoint[1]),0.5))&&
                ((kpoint[2]==0)||aeq(fabs(kpoint[2]),0.5))&&
                (flags&BANDPARITY)) inv_parity(psi,fft,nb,kpoint);

            /* Was the parity all we were requested to report? */
            if (!(flags&BANDS)) continue;

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
                ngpts=nfft[0]*nfft[1]*nfft[2];
                free(psi);
                psi=dptr;
              }
            }

	  
      
            ffft[0]=nfft[2];
            ffft[1]=nfft[1];
            ffft[2]=nfft[0];
            fft3d(psi,ffft,1);

            dptr=malloc(nfft[0]*nfft[1]*nfft[2]*sizeof(double));
            if (!dptr) error_exit("Malloc error for grid data");
            band2real(psi,dptr,nfft,kpoint);
	 
            /* Do we need to rescale? */
            if (((flags&RAW)==0)&&((flags&BANDPHASE)==0)){ /* Yes */
              if (flags&BANDDEN) conv=1/c->vol;
              else conv=1/sqrt(c->vol);
              if (debug>2) fprintf(stderr,"Scaling wavefun by %f\n",conv);
              for(i=0;i<nfft[0]*nfft[1]*nfft[2];i++) dptr[i]*=conv;
            }

            if (elect->occ)
              occ=elect->occ[elect->nspins*nbnd*(k2-1)+ns*nbnd+(nb-1)];
            else
              occ=1;
            band_store(&g,dptr,occ,k->kpts[k2-1].wt,ii+1,
                       ns,k2,nb,elect,m,nfft);
	  } /* End spinor loop */
	}
	else fseek(infile,4+16*npol*igwx,SEEK_CUR);
    

      } /* End loop over bands */
  
      free(wtmp);
      free(psi);
      free(pwgrid);
      fclose(infile);
    } /* End loop over kpoints */

  } /* End loop over spins */
  
  free(file);
  
}
