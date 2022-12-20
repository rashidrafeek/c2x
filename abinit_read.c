/* Read data from an Abinit 8.x Fortran binary file.
 *
 * Unclear which versions of abinit this might work with
 * tested with 8.6.3
 *
 * Defaults to rescaling things recognised as densities from Bohr to A.
 * Does not if "raw" specified (-R flag)
 */

/* Copyright (c) 2018 MJ Rutter 
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

/* Note about abinit and variables:
 *
 * nspden is the number of spins in the density file
 * nsppol is the number of spins in the wavefunction file, and hence
 *   also number of occupations and evalues
 *
 * the occupation of a band is 2 if doubly occupied
 *
 * Note https://docs.abinit.org/guide/abinit/#64-the-header is
 * (currently) wrong. The format is
 *
 * Record 1: as URL
 * Record 2: Ditto
 * Record 3: Ditto
 * Record 4: Atomic positions and total energy
 * Record 5: Undocumented
 * Record 6 to 5+nspec: pseudopotentials
 * Two more records if usepaw==1
 */

#include<stdio.h>
#include<stdlib.h> /* malloc */
#include<math.h>

#include "c2xsf.h"

static void abinit_header_read(FILE* infile, struct unit_cell *c,
			       struct contents *m, struct kpts *kp,
			       struct es *elect,
			       int fft[3], int **gamma, int *fileform){
  int i,j,tmp,reclen;
  char version[6];
  int headform;
  int natoms,nsym,npsp,usepaw,ntypes,bantot,*nband,nbands;
  long offset;
  double *posn,*atnos;
  double x;
  int *atypes;

  nband=NULL;
  atnos=NULL;
  atypes=NULL;
  
  /* Read first record */
  tmp=0;
  fread(&tmp,4,1,infile);
  if(tmp!=14) error_exit("Unexpected first record length in abinit_read");
  fread(version,6,1,infile);
  fread(&headform,4,1,infile);
  fread(fileform,4,1,infile);
  tmp=0;
  fread(&tmp,4,1,infile);
  if (tmp!=14) error_exit("Unexpected first record length in abinit_read");

  if (debug) fprintf(stderr,"Abinit version %6s\n",version);

  if (headform<80)
    error_exit("Unable to read files from abinit <8.0");
  
  /* Read second record */

  fread(&reclen,4,1,infile);
  offset=ftell(infile);
  fread(&bantot,4,1,infile);
  fseek(infile,12,SEEK_CUR);
  fread(&natoms,4,1,infile);
  fread(fft,12,1,infile);
  fread(&kp->n,4,1,infile);
  fread(&elect->nspins,4,1,infile);
  fread(&elect->nspinors,4,1,infile);
  fread(&elect->nbspins,4,1,infile);
  fread(&nsym,4,1,infile);
  fread(&npsp,4,1,infile);
  fread(&ntypes,4,1,infile);
  fseek(infile,2*4,SEEK_CUR);
  fread(&usepaw,4,1,infile);
  fread(&elect->cut_off,8,1,infile);
  elect->cut_off*=H_eV;
  fseek(infile,6*8,SEEK_CUR);
  c->basis=malloc(72);
  if (!c->basis) error_exit("Malloc error for basis");
  fread(c->basis,9,8,infile);
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[i][j]*=BOHR;
  
  if (debug>1){
    fprintf(stderr,"natoms: %d\nnkpts: %d\nnspins: %d\nnspinors: %d\n",
            natoms,kp->n,elect->nspins,elect->nspinors);
    fprintf(stderr,"npsp: %d\nusepaw: %d\n",
            npsp,usepaw);
    fprintf(stderr,"FFT grid: %dx%dx%d\n",fft[0],fft[1],fft[2]);
  }
  fseek(infile,offset+reclen+4,SEEK_SET);

  m->n=natoms;
  m->atoms=malloc(natoms*sizeof(struct atom));
  if (!m->atoms) error_exit("Malloc error for atoms");
  init_atoms(m->atoms,m->n);
  
  /* Third record */

  reclen=0;
  fread(&reclen,4,1,infile);
  offset=ftell(infile);
  if (reclen==0) error_exit("Read error");
  nband=malloc(bantot*sizeof(int));
  if (!nband) error_exit("Malloc error for nband");
  *gamma=malloc(kp->n*sizeof(int));
  if (!(*gamma)) error_exit("Malloc error in abinit_header_read");
  fread(*gamma,4,kp->n,infile);
  fread(nband,4,kp->n*elect->nbspins,infile);
  nbands=nband[0];
  for(i=0;i<kp->n*elect->nbspins;i++)
    if (nbands<nband[i]) nbands=nband[i];
  free(nband);
  nband=NULL;
  elect->nbands=nbands;
  /* Need to find our atomic numbers */
  /* atypes gives the species number of each atom (starting from 1) */
  /* atnos (double!) gives the atomic number of each species */
  tmp=kp->n+npsp+nsym+9*nsym;
  fseek(infile,4*tmp,SEEK_CUR);
  atypes=malloc(4*natoms);
  fread(atypes,natoms*4,1,infile);
  atnos=malloc(8*ntypes);
  /* Now kpoints */
  kp->kpts=malloc(kp->n*sizeof(struct atom));
  if (!kp->kpts) error_exit("Malloc error for kpts");
  for(i=0;i<kp->n;i++)
    fread(kp->kpts[i].frac,8,3,infile);
  fseek(infile,offset+reclen-8*kp->n-8*ntypes,SEEK_SET);
  fread(atnos,ntypes*8,1,infile);
  for(i=0;i<natoms;i++)
    m->atoms[i].atno=atnos[atypes[i]-1];
  for(i=0;i<kp->n;i++)
    fread(&kp->kpts[i].wt,8,1,infile);
  fseek(infile,offset+reclen+4,SEEK_SET);

  /* Atom positions: fourth record */

  fread(&reclen,4,1,infile);
  offset=ftell(infile);
  fseek(infile,8,SEEK_CUR);
  posn=malloc(natoms*3*sizeof(double));
  if (!posn) error_exit("Malloc error for atomic positions");
  fread(posn,natoms*3*sizeof(double),1,infile);
  for(i=0;i<natoms;i++)
    for(j=0;j<3;j++)
      m->atoms[i].frac[j]=*(posn+3*i+j);
  free(posn);
  elect->energy=malloc(sizeof(double));
  fread(elect->energy,8,1,infile);
  *elect->energy*=H_eV;
  elect->e_fermi=malloc(sizeof(double));
  fread(elect->e_fermi,8,1,infile);
  *elect->e_fermi*=H_eV;
  fseek(infile,offset+reclen+4,SEEK_SET);

  real2rec(c);

  /* Reduce all atoms to unit cell, and deal with atoms which have
   * drifted outside of unit cell by tiny amounts so that -epsilon
   * ends up at zero, not one
   */
  x=tol;
  tol*=1e-10;
  reduce_cell_tol(m->atoms,m->n,c->basis,tol);
  tol=x;
  addabs(m->atoms,m->n,c->basis);

  /* Read another (undocumented?) record */

  tmp=0;
  fread(&tmp,4,1,infile);
  if (tmp==0) error_exit("Read error");
  fseek(infile,tmp+4,SEEK_CUR);
  
  /* Pseudopotentials */

  for(i=0;i<npsp;i++){
    reclen=tmp=0;
    fread(&reclen,4,1,infile);
    if (npsp==ntypes){  /* reuse atnos for pseudo charges */
      fseek(infile,132+8,SEEK_CUR);
      fread(atnos+i,8,1,infile);
      fseek(infile,reclen-132-16,SEEK_CUR);
    }
    else
      fseek(infile,reclen,SEEK_CUR);
    fread(&tmp,4,1,infile);
    if ((tmp==0)||(tmp!=reclen))
      error_exit("Read error in pseudopotential block");
  }

  if (npsp==ntypes){
  for(i=0;i<natoms;i++)
    m->atoms[i].chg=atnos[atypes[i]-1];
  }

  if (usepaw==1){
    if (headform>=56){
      for(i=0;i<2;i++){ /* Two records to read */
        reclen=tmp=0;
        fread(&reclen,4,1,infile);
        fseek(infile,reclen,SEEK_CUR);
        fread(&tmp,4,1,infile);
        if ((tmp==0)||(tmp!=reclen))
          error_exit("Read error in PAW block");
      }
    }
    else{
      error_exit("Unable to read PAW block from old version of abinit");
    }
  }

  if (atypes) free(atypes);
  if (atnos) free(atnos);

}

void abinit_charge_read(FILE* infile, struct unit_cell *c,
			struct contents *m, struct kpts *kp,
                        struct grid *gptr, struct es *elect){
  int i,j,k;
  int reclen,fft[3],fileform,*gamma;
  double scale,*ptr,*chgden;

  abinit_header_read(infile,c,m,kp,elect,fft,&gamma,&fileform);

  if ((!(flags&CHDEN))&&(!(flags&SPINDEN))) return;
  
  fread(&reclen,4,1,infile);
  if (reclen!=8*fft[0]*fft[1]*fft[2]){
    fprintf(stderr,"Volumetric data not understood, so not reading\n");
    fprintf(stderr,"Expected %d bytes, found %d\n",
            8*fft[0]*fft[1]*fft[2],reclen);
    return;
  }

  if (gptr->next) gptr=gptr->next;
  gptr->next=malloc(sizeof(struct grid));
  if (!gptr->next) error_exit("Malloc error for struct grid");
  gptr->next->next=NULL;
  gptr->next->data=NULL;

  for(i=0;i<3;i++) gptr->size[i]=fft[i];
  if(!(gptr->data=malloc(gptr->size[0]*gptr->size[1]*
                         gptr->size[2]*sizeof(double))))
    error_exit("Malloc error for grid data");

  
  if (fread(gptr->data,1,reclen,infile)!=reclen){
    fprintf(stderr,"Error reading grid data: ignoring grid\n");
    free(gptr->data);
    gptr->data=NULL;
    return;
  }

  fseek(infile,4,SEEK_CUR);
  
  /* Transpose data grid from Fortran to C ordering */

  if(!(ptr=malloc(gptr->size[0]*gptr->size[1]*
                         gptr->size[2]*sizeof(double))))
    error_exit("Malloc error for temporary data");

  for(i=0;i<fft[0];i++)
    for(j=0;j<fft[1];j++)
      for(k=0;k<fft[2];k++)
        ptr[k+j*fft[2]+i*fft[2]*fft[1]]=gptr->data[k*fft[0]*fft[1]+j*fft[0]+i];
  free(gptr->data);
  gptr->data=ptr;
  
  if ((fileform>=50)&&(fileform<100)){
    gptr->name="Density";
    if (!(flags&RAW)){
      if (debug) fprintf(stderr,"Rescaling data from Bohr^-3 to A^-3\n");
      scale=1/(BOHR*BOHR*BOHR);
      for(i=0;i<gptr->size[0]*gptr->size[1]*gptr->size[2];i++)
        gptr->data[i]*=scale;
    }
  }
  else if ((fileform>=100)&&(fileform<120)) gptr->name="Potential";
  else gptr->name="Unknown";

  if ((fileform!=52)||(elect->nspins!=2)) return;
  if (!(flags&SPINDEN)) return;

  /* We have some spin to read */

  fread(&reclen,4,1,infile);
  if (reclen!=8*fft[0]*fft[1]*fft[2]){
    fprintf(stderr,"Spin data not understood, so not reading\n");
    fprintf(stderr,"Expected %d bytes, found %d\n",
            8*fft[0]*fft[1]*fft[2],reclen);
    return;
  }
  if (!(ptr=malloc(reclen)))
    error_exit("Malloc error for temporary data");
  if (fread(ptr,1,reclen,infile)!=reclen){
    fprintf(stderr,"Error reading temporary spin data: ignoring\n");
    free(ptr);
    return;
  }

  chgden=gptr->data;
  if (flags&CHDEN){
    if (gptr->next) gptr=gptr->next;
    gptr->next=malloc(sizeof(struct grid));
    if (!gptr->next) error_exit("Malloc error for struct grid");
    gptr->next->next=NULL;
    gptr->next->data=NULL;
    for(i=0;i<3;i++) gptr->size[i]=fft[i];
    if(!(gptr->data=malloc(gptr->size[0]*gptr->size[1]*
                           gptr->size[2]*sizeof(double))))
      error_exit("Malloc error for spin data");
  }
  else{ /* Discard charge density data that was not requested 
         * by overwriting its grid entry */
    if(!(gptr->data=malloc(gptr->size[0]*gptr->size[1]*
                           gptr->size[2]*sizeof(double))))
      error_exit("Malloc error for spin data");
  }    

  /* Transpose data grid from Fortran to C ordering */

  for(i=0;i<fft[0];i++)
    for(j=0;j<fft[1];j++)
      for(k=0;k<fft[2];k++)
        gptr->data[k+j*fft[2]+i*fft[2]*fft[1]]=ptr[k*fft[0]*fft[1]+j*fft[0]+i];
  free(ptr);
  
  if (!(flags&RAW)){
    if (debug) fprintf(stderr,"Rescaling spin from Bohr^-3 to A^-3 "
                       "and calculating net spin\n");
    scale=2/(BOHR*BOHR*BOHR);
    for(i=0;i<gptr->size[0]*gptr->size[1]*gptr->size[2];i++)
      gptr->data[i]=gptr->data[i]*scale-chgden[i];
  }
  if (!(flags&CHDEN)) free(chgden);

  gptr->name="Spin";

}


void abinit_psi_read(FILE* infile, struct unit_cell *c,
                     struct contents *m, struct kpts *kp, struct grid *gptr,
                     struct es *elect, int *i_grid){
  int i,j,ikpt,nb,offset,isppol;
  int n0,n1,n2;
  int reclen,fft[3],fileform,*gamma,goff[3];
  int npw,nspinor,ispinor,nbnd,nfftpts,off,ffft[3],nfft[3];
  int *pwgrid;
  double scale,*psi,*dptr,*kpoint;

  gamma=NULL;
  abinit_header_read(infile,c,m,kp,elect,fft,&gamma,&fileform);

  if (debug>1){
    fprintf(stderr,"Gamma[] = [");
    for(i=0;i<kp->n;i++)
      fprintf(stderr,"%d ",gamma[i]);
    fprintf(stderr,"]\n");
  }
  
  if ((!(flags&BANDS))&&(!(flags&OCCUPANCIES))&&(!(flags&BANDPARITY)))
    return;
  
  nfftpts=fft[0]*fft[1]*fft[2];

  if (debug>1)
    fprintf(stderr,"nsppol=%d nspins=%d\n",elect->nbspins,elect->nspins);
  
  elect->occ=malloc(elect->nbands*elect->nbspins*kp->n*sizeof(double));
  if (!elect->occ) error_exit("Malloc error for occupancies");
  elect->eval=malloc(elect->nbands*elect->nbspins*kp->n*sizeof(double));
  if (!elect->occ) error_exit("Malloc error for eigenvalues");
  for(i=0;i<elect->nbands*elect->nbspins*kp->n;i++){
    elect->eval[i]=0;
    elect->occ[i]=0;
  }
  pwgrid=NULL;
  
  for(isppol=0;isppol<elect->nbspins;isppol++){
    for(ikpt=0;ikpt<kp->n;ikpt++){
      kpoint=kp->kpts[ikpt].frac;
      fread(&reclen,4,1,infile);
      if (reclen!=12) error_exit("Unexpected record length");
      fread(&npw,4,1,infile);
      if (debug>1) fprintf(stderr,"nplwv=%d\n",npw);
      if ((debug>1)&&(gamma[ikpt]>1))
        fprintf(stderr,"Gamma point storage for this kpt\n");
      if (gamma[ikpt]>1){
        for(i=0;i<3;i++) goff[i]=0;
        if (gamma[ikpt]&1) goff[0]=1;
        if (gamma[ikpt]>5) goff[1]=1;
        if ((gamma[ikpt]-2)&2) goff[2]=1;
      }
      fread(&nspinor,4,1,infile);
      fread(&nbnd,4,1,infile);
      fseek(infile,4,SEEK_CUR);
      if (nbnd>elect->nbands) error_exit("nbnd for kpt > nbands");
      pwgrid=realloc(pwgrid,3*npw*sizeof(int));
      if (!pwgrid) error_exit("Malloc error for pwgrid");
      fread(&reclen,4,1,infile);
      if (reclen!=12*npw) error_exit("Unexpected record length for npw");
      fread(pwgrid,4,3*npw,infile);
      fseek(infile,4,SEEK_CUR);
      for(i=0;i<3;i++)
        for(j=0;j<npw;j++)
          if (pwgrid[3*j+i]<0) pwgrid[3*j+i]+=fft[i];
      fread(&reclen,4,1,infile);
      if (reclen!=16*nbnd) error_exit("Unexpected record length for eval/occ");
      off=elect->nbands*elect->nbspins*ikpt+isppol*elect->nbands;
      fread(elect->eval+off,8,nbnd,infile);
      fread(elect->occ+off,8,nbnd,infile);
      /* Abinit has band occupancy of two if nbspins=1 */
      if (elect->nbspins==1)
        for(i=0;i<nbnd;i++) elect->occ[off+i]*=0.5;
      fseek(infile,4,SEEK_CUR);
      for(nb=0;nb<nbnd;nb++){
        i=fread(&reclen,4,1,infile);
        if (reclen!=16*npw*nspinor){
          error_exit("Unexpected record length for band");
        }
        if ((flags&BANDPARITY||((flags&BANDS)&&
				((inrange(nb+1,elect->band_range))&&
				 (inrange(ikpt+1,elect->kpt_range))&&
				 ((elect->nbspins==1)||
                                  (inrange(isppol,elect->spin_range))))))){
          if (debug>1) fprintf(stderr,"Reading band %d\n",nb+1);
          for(ispinor=1;ispinor<=nspinor;ispinor++){
            dptr=malloc(16*npw);
            if (!dptr) error_exit("Malloc error for band");
            if ((nspinor==1)||(inrange(ispinor,elect->spin_range))){
              fread(dptr,16,npw,infile);
              if (ispinor==nspinor) fseek(infile,4,SEEK_CUR);
            }
            else{
              fseek(infile,16*npw,SEEK_CUR);
              if (ispinor==nspinor) fseek(infile,4,SEEK_CUR);
              continue;
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
            if (gamma[ikpt]>1){ /* construct psi(-k)=conjg(psi(k)) */
              for(i=0;i<npw;i++){
                if ((gamma[ikpt]==2)&&(pwgrid[3*i]==0)&&
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
            free(dptr);


            if (((kpoint[0]==0)||aeq(fabs(kpoint[0]),0.5))&&
                ((kpoint[1]==0)||aeq(fabs(kpoint[1]),0.5))&&
                ((kpoint[2]==0)||aeq(fabs(kpoint[2]),0.5))&&
                (flags&BANDPARITY)) inv_parity(psi,fft,nb,kpoint);

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

	  
              ffft[0]=fft[2];
              ffft[1]=fft[1];
              ffft[2]=fft[0];
              fft3d(psi,ffft,1);
              dptr=malloc(nfftpts*sizeof(double));
              if(!dptr) error_exit("Malloc error for grid data");
              band2real(psi,dptr,fft,kp->kpts[ikpt].frac);
              free(psi);

              /* Do we need to rescale? */
              if (((flags&RAW)==0)&&((flags&BANDPHASE)==0)){ /* Yes */
                if (flags&BANDDEN) scale=1/c->vol;
                else scale=1/sqrt(c->vol);
                if (debug>2) fprintf(stderr,"Scaling wavefun by %f\n",scale);
                for(i=0;i<nfftpts;i++) dptr[i]*=scale;
              }
            
              band_store(&gptr,dptr,elect->occ[off+nb],kp->kpts[ikpt].wt,
                         ispinor,isppol,ikpt+1,nb+1,elect,m,fft);
            } /* end if (flags&BANDS) */
          } /* end spinor loop */
        }
        else fseek(infile,reclen+4,SEEK_CUR);
      } /* end for bands */
    }
  }
  if (pwgrid) free(pwgrid);
  if (elect->eval) /* We store these in eV */
    for(i=0;i<elect->nbands*elect->nbspins*kp->n;i++)
      elect->eval[i]*=H_eV;
}

