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



#include<stdio.h>
#include<stdlib.h> /* malloc */

#include "c2xsf.h"

void abinit_read(FILE* infile, struct unit_cell *c, struct contents *m,
               struct grid *gptr){
  int i,j,k,tmp,reclen;
  char version[6];
  int headform,fileform;
  int nbands,natoms,fft[3],nkpts,nspins,nspinors,nsym,npsp,
    usepaw,nsppol,ntypes;
  long offset;
  double *posn,*atnos,*ptr,*chgden;
  double scale,x;
  int *atypes;

  /* Read first record */
  tmp=0;
  fread(&tmp,4,1,infile);
  if(tmp!=14) error_exit("Unexpected first record length in abinit_read");
  fread(version,6,1,infile);
  fread(&headform,4,1,infile);
  fread(&fileform,4,1,infile);
  tmp=0;
  fread(&tmp,4,1,infile);
  if (tmp!=14) error_exit("Unexpected first record length in abinit_read");

  if (debug) fprintf(stderr,"Abinit version %6s\n",version);

  if (headform<80)
    error_exit("Unable to read files from abinit <8.0");
  
  /* Read second record */

  fread(&reclen,4,1,infile);
  offset=ftell(infile);
  fread(&nbands,4,1,infile);
  fseek(infile,12,SEEK_CUR);
  fread(&natoms,4,1,infile);
  fread(fft,12,1,infile);
  fread(&nkpts,4,1,infile);
  fread(&nspins,4,1,infile);
  fread(&nspinors,4,1,infile);
  fread(&nsppol,4,1,infile);
  fread(&nsym,4,1,infile);
  fread(&npsp,4,1,infile);
  fread(&ntypes,4,1,infile);
  fseek(infile,2*4,SEEK_CUR);
  fread(&usepaw,4,1,infile);
  fseek(infile,7*8,SEEK_CUR);
  c->basis=malloc(72);
  if (!c->basis) error_exit("Malloc error for basis");
  fread(c->basis,9,8,infile);
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[i][j]*=BOHR;
  
  if (debug>1){
    fprintf(stderr,"natoms: %d\nnkpts: %d\nnpins: %d\nnspinors: %d\n",
            natoms,nkpts,nspins,nspinors);
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
  /* Need to find our atomic numbers */
  /* atypes gives the species number of each atom (starting from 1) */
  /* atnos (double!) gives the atomic number of each species */
  tmp=nkpts+nkpts*nsppol+nkpts+npsp+nsym+9*nsym;
  fseek(infile,4*tmp,SEEK_CUR);
  atypes=malloc(4*natoms);
  fread(atypes,natoms*4,1,infile);
  atnos=malloc(8*ntypes);
  fseek(infile,offset+reclen-8*nkpts-8*ntypes,SEEK_SET);
  fread(atnos,ntypes*8,1,infile);
  for(i=0;i<natoms;i++)
    m->atoms[i].atno=atnos[atypes[i]-1];
  fseek(infile,offset+reclen+4,SEEK_SET);

  /* Atom positions */

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
  fseek(infile,offset+reclen+4,SEEK_SET);

  real2rec(c);

  /* Reduce all atoms to unit cell, and deal with atoms which have
   * drifted outside of unit cell by tiny amounts so that -epsilon
   * ends up at zero, not one
   */
  x=tol;
  tol*=1e-10;
  reduce_cell_tol(m->atoms,m->n,c->basis);
  tol=x;
  addabs(m->atoms,m->n,c->basis);

  /* Read another record */

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

  if ((fileform!=52)||(nspins!=2)) return;
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