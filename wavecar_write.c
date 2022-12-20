/* Produce output in the style of
 * a VASP WAVECAR file
 *
 * Need to add re-indexing if it is to work with input from other codes
 *
 * (c) MJR 2/2020
 */

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#include "c2xsf.h"

int vasp_grid_fill(int *pwgrid, double recip[3][3], double kpt[3],
                   double ecut, double g2cut, int fft[3], int nplwv,
                   int gamma);

double g2max(double recip[3][3], int *pwgrid, int nplwv, int fft[3],
             double *kpt);

void wavecar_write(double *psi, int *pwgrid, int nplwv_in, int fft[3],
		  int gamma, int goff[3],
		  struct unit_cell *c, struct contents *m, struct kpts *kp,
		  int ikpt, int isppol, int nb,
		  double *eval, double occ, struct es *e){
  int i,j,nkpt,old_reclen,tmp,max_nplwv,*map,sizeofcomplex,trivial_map,hit;
  int map_inv[3];
  long fposn;
  FILE *outfile;
  double *kpt,c_double[2];
  float *psi_float,c_float[2];
  double junk,version,one,g2cut;
  static int firstcall=1,reclen=0,*remap=NULL,nplwv;

  kpt=kp->kpts[ikpt].frac;
  outfile=dict_get(m->dict,"out_file_handle");

  if (!outfile) error_exit("No output file in wavecar_write");

  if (isppol==-1) error_exit("Cannot write wavecar with spinor");
  
  one=1;
  g2cut=0;
  sizeofcomplex=8;
  if (flags&HIPREC) sizeofcomplex=16;

  
  if (firstcall==1){ /* write header */
    firstcall=0;
    /* Count kpoints */
    nkpt=0;
    for(i=0;i<kp->n;i++)
      if (inrange(i+1,e->kpt_range)) nkpt++;
    nplwv=nplwv_in;
    if (gamma){
      if (e->cut_off==0)
	g2cut=g2max(c->recip,pwgrid,nplwv_in,fft,kpt);
      nplwv=vasp_grid_fill(NULL,c->recip,kpt,e->cut_off,g2cut,fft,0,0);
    }
    if (nkpt==1)
      max_nplwv=nplwv;
    else
      max_nplwv=e->max_nplwv;

    if (max_nplwv==0)
      error_exit("WAVECAR output not possible as maximal nplwv not known");

    if (!e->eval)
      error_exit("WAVECAR output not possible as evals not read");
    
    reclen=((max_nplwv+1)/2)*16;
    if (flags&HIPREC) reclen=max_nplwv*16;
    old_reclen=((3*e->nbands+1)/2)*16+32;
    

    if (old_reclen<=reclen)
      version=45200;
    else{
      if (flags&ALT_OUT){
        version=45200;
        reclen=old_reclen;  /* As old_reclen>reclen */
      }
      else
        version=53300;
    }
    if (flags&HIPREC) version+=10;
    junk=reclen;    /* Assume single prec */
    fwrite(&junk,8,1,outfile);
    junk=e->nspins;
    fwrite(&junk,8,1,outfile);
    fwrite(&version,8,1,outfile);
    /* end first record */
    fseek(outfile,reclen,SEEK_SET);
    junk=kp->n;
    fwrite(&junk,8,1,outfile);
    junk=e->nbands;
    fwrite(&junk,8,1,outfile);
    junk=e->cut_off;
    fwrite(&junk,8,1,outfile);
    fwrite(c->basis,8,9,outfile);
    if (e->e_fermi)
      junk=*(e->e_fermi);
    else{
      junk=0;
      fprintf(stderr,"Warning: Fermi energy unknown. Written as zero\n");
    }
    fwrite(&junk,8,1,outfile);
    /* end second record */
    fseek(outfile,2*reclen,SEEK_SET);
  }

  if (nb==1){
    /* Do we need to remap? */
    if (remap) free(remap);
    remap=NULL;
    if (e->cut_off==0)
      g2cut=g2max(c->recip,pwgrid,nplwv_in,fft,kpt);
    nplwv=nplwv_in;
    if (gamma){
      nplwv=vasp_grid_fill(NULL,c->recip,kpt,e->cut_off,g2cut,fft,0,0);
      if (debug) fprintf(stderr,
			 "Expanding Gamma point from %d to %d plane waves\n",
			 nplwv_in,nplwv);
    }
    remap=malloc(nplwv*sizeof(int));
    map=malloc(3*nplwv*sizeof(int));
    if (!map) error_exit("Malloc error for map");
    if (!remap) error_exit("Malloc error for remap");
    vasp_grid_fill(map,c->recip,kpt,e->cut_off,g2cut,fft,nplwv,0);
    trivial_map=1;
    for(i=0;i<nplwv;i++){
      if((i<nplwv_in)&&(map[3*i]==pwgrid[3*i])&&(map[3*i+1]==pwgrid[3*i+1])
         &&(map[3*i+2]==pwgrid[3*i+2])) remap[i]=i;
      else{
        trivial_map=0;
        hit=0;
        for(j=0;j<nplwv_in;j++){
          if((map[3*i]==pwgrid[3*j])&&(map[3*i+1]==pwgrid[3*j+1])
             &&(map[3*i+2]==pwgrid[3*j+2])){
            remap[i]=j;
            hit=1;
            break;
          }
        }
        if (hit==0){
          if (gamma){
            for(j=0;j<3;j++)
              map_inv[j]=fft[j]-map[3*i+j]-goff[j];
            for(j=0;j<3;j++)
              if (map_inv[j]==fft[j]) map_inv[j]=0;
            for(j=0;j<nplwv_in;j++){
              if((map_inv[0]==pwgrid[3*j])&&(map_inv[1]==pwgrid[3*j+1])
                 &&(map_inv[2]==pwgrid[3*j+2])){
                remap[i]=-j;
                hit=1;
                break;
              }
            }
          }
          if (hit==0){
            fprintf(stderr,"Failed to find reverse mapping for %d %d %d\n",
                    map[3*i],map[3*i+1],map[3*i+2]);
            if (gamma)
              fprintf(stderr,"Also tried %d %d %d\n",map_inv[0],map_inv[1],
                      map_inv[2]);
            exit(1);
          }
        }
      }
    }
    free(map);
    if (trivial_map){
      free(remap);
      remap=NULL;
    }
  
    
    junk=nplwv;
    fwrite(&junk,8,1,outfile);
    fwrite(kp->kpts[ikpt].frac,8,3,outfile);
    tmp=(ikpt*e->nbspins+isppol)*e->nbands;
    junk=0;
    for(j=0;j<e->nbands;j++){
      fwrite(e->eval+tmp+j,8,1,outfile);
      fwrite(&junk,8,1,outfile); /* Imaginary part of energy?! */
      if (e->occ)
        fwrite(e->occ+tmp+j,8,1,outfile);
      else
        fwrite(&one,8,1,outfile);
    }
  }

  fposn=ftell(outfile);
  tmp=fposn%reclen;
  if (tmp){
    fposn=(fposn/reclen+1)*reclen;
    fseek(outfile,fposn,SEEK_SET);
  }

  if ((flags&HIPREC)==0){
    psi_float=malloc(2*nplwv_in*sizeof(float));
    if(!psi_float) error_exit("Malloc error for psi_float in wavecar_write");
    for(i=0;i<2*nplwv_in;i++) psi_float[i]=psi[i];
    if (remap==NULL)
      fwrite(psi_float,4,2*nplwv,outfile);
    else
      for(i=0;i<nplwv;i++){
        if (remap[i]>=0)
          fwrite(psi_float+2*remap[i],4,2,outfile);
        else{
          c_float[0]=*(psi_float-2*remap[i]);
          c_float[1]=-*(psi_float-2*remap[i]+1);
          fwrite(c_float,4,2,outfile);
        }
      }
    free(psi_float);
  }
  else{
    if (remap==NULL)
      fwrite(psi,8,2*nplwv,outfile);
    else
      for(i=0;i<nplwv;i++){
        if (remap[i]>=0)
          fwrite(psi+2*remap[i],8,2,outfile);
        else{
          c_double[0]=*(psi-2*remap[i]);
          c_double[1]=-*(psi-2*remap[i]+1);
          fwrite(c_double,8,2,outfile);
        }
      }
  }
    
  /* If this is the last call, which we don't know, need to get the filesize
   * correct, rounding up to multiple of reclen. ftruncate() does not suffice.
   */
  if (sizeofcomplex*nplwv<reclen){
    fseek(outfile,fposn+reclen-1,SEEK_SET);
    i=0;
    fwrite(&i,1,1,outfile);
  }
  
}
