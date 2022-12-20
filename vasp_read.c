/* Read a VASP CHGCAR or CHG file
 */

/* Copyright (c) 2017 MJ Rutter 
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

/* Format assumed to be:
 *
 * title (may be blank)
 * scale factor for axes (see below)
 * ax ay az
 * bx by bz
 * cx cy cz
 * atomic symbols of species present
 * number of atoms of each species present
 * A string starting with a [cCkK] (atom co-ords are cartesian),
 *           or not (atom co-ords are relative)
 * atomic co-ordinates, first species x number of atoms of that species,
 *           then second, ...
 * a blank line
 * three integers as grid size
 * charge density data, probably multiple items per line
 * if CHGCAR, augmentation occupancies
 * if CHGCAR and spin present, magnetic moments of atoms, muliple lines of
 *   multiple items per line
 * if spin present, three integers as grid size
 * if spin present, spin density data, probably multiple items per line
 * if vector spin, repeat two more times from mag moment section
 *
 *
 * the axis scale factor may be:
 *    i/ single float, >0: scale all axes linearly by this
 *   ii/ single float, <0: scale all axes to make cell volume abs(scale)
 *  iii/ three floats: linear scale factors for each of x, y and z
 *                     (not a, b and c)
 *
 */

#include<stdio.h>
#include<stdlib.h> /* malloc */
#include<string.h> /* memcpy */
#include<math.h>
#include<ctype.h>

#include "c2xsf.h"

#define LINE_SIZE 2049

static int formatted_read(double *dptr, unsigned int len, FILE *infile,
                          char *first);

void vasp_read(FILE* infile, struct unit_cell *c, struct contents *m,
                struct grid *gptr, struct es *elect){
  int i,j,k,nspec,natoms,frac,fft[3],ffts[3],isize;
  int *nionsp,*species;
  char buffer[LINE_SIZE+1], *ptr,*p2;
  double scale,ascale[3],vscale,*dptr,*dptr2,*magmom;

  dptr=dptr2=magmom=NULL;
  
  if (debug>2) fprintf(stderr,"vasp_read called\n");

  /* title */
  
  if(!fgets(buffer,LINE_SIZE,infile)){
    fprintf(stderr,"vasp_read() has failed to read the first line!\n");
    exit(1);
  }
  /* Was the title more than spaces and an end-of-line? */
  ptr=buffer;
  while((*ptr==' ')||(*ptr=='\n')) ptr++;
  if (*ptr){
    m->title=malloc(strlen(buffer)+1);
    if (!m->title){
      fprintf(stderr,"Malloc error for title in vasp_read()\n");
      exit(1);
    }
    strcpy(m->title,buffer);
    /* Remove newline */
    m->title[strlen(m->title)-1]=0;
  }

  /* Scale factor */

  fgets(buffer,LINE_SIZE,infile);
  vscale=0;
  ascale[0]=ascale[1]=ascale[2]=1;
  if (sscanf(buffer,"%lf %lf %lf",ascale,ascale+1,ascale+2)!=3){
    if (sscanf(buffer,"%lf",ascale)!=1)
      error_exit("Failed to parse scale factor in vasp_read()\n");
    if (ascale[0]<0){
      vscale=-ascale[0];
      ascale[0]=1;
    }
    else ascale[2]=ascale[1]=ascale[0];
  }
    
  /* Axes */

  if (!(c->basis=malloc(72)))
    error_exit("Malloc error in vasp_read for c->basis");
  
  for(i=0;i<3;i++){
    fgets(buffer,LINE_SIZE,infile);
    if(sscanf(buffer,"%lf %lf %lf",c->basis[i],c->basis[i]+1,c->basis[i]+2)!=3){
      fprintf(stderr,"Parse error in cell_read for axis %d\n",i);
      exit(1);
    }
  }

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[i][j]*=ascale[j];

  real2rec(c);

  if (vscale!=0){
    vscale=pow(vscale/fabs(c->vol),1./3.);
    if (debug) fprintf(stderr,
                       "Volume scaling results in linear scaling of %g\n",
                       vscale);
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
        c->basis[i][j]*=vscale;
    real2rec(c);
  }
  
  /* Species, symbols */
  fgets(buffer,LINE_SIZE,infile);
  nspec=0;

  ptr=buffer;
  while(*ptr){
    while(*ptr==' ') ptr++;
    if (isalpha(*ptr)){
      nspec++;
      while (isalpha(*ptr)) ptr++;
    }
    else if (*ptr=='\n') break;
    else{
      fprintf(stderr,"Unexpected character %c in species specification\n",*ptr);
      exit(1);
    }
  }
  if (nspec==0)
    error_exit("Failed to parse species specification\n");
  
  species=malloc(nspec*sizeof(int));
  nionsp=malloc(nspec*sizeof(int));

  ptr=buffer;
  i=0;
  while((*ptr)&&(i<nspec)){
    while(*ptr==' ') ptr++;
    if (isalpha(*ptr)){
      p2=ptr;
      while(isalpha(*p2)) p2++;
      *p2=0;
      species[i]=atsym2no(ptr);
      if (species[i]==0)
        fprintf(stderr,"Warning: failed to parse %s as atomic symbol\n",ptr);
      ptr=p2+1;
      i++;
    }
  }

  fgets(buffer,LINE_SIZE,infile);
  ptr=buffer;
  for(i=0;i<nspec;i++)
    nionsp[i]=strtol(ptr,&ptr,10);

  if(!nionsp[nspec-1])
    error_exit("Failed to parse number of each species\n");

  natoms=0;
  for(i=0;i<nspec;i++)
    natoms+=nionsp[i];
  
  fgets(buffer,LINE_SIZE,infile);
  frac=1;
  if ((buffer[0]=='c')||(buffer[0]=='C')||
      (buffer[0]=='k')||(buffer[0]=='K')) frac=0;

  if (debug>1)
    fprintf(stderr,"%d species, %d atoms, %s co-ordinates\n",nspec,natoms,
            frac?"fractional":"cartesian");

  m->n=natoms;
  m->atoms=malloc(natoms*sizeof(struct atom));
  if(!m->atoms) error_exit("Malloc error for m->atoms in vasp_read");

  k=0;
  for(i=0;i<nspec;i++){
    for(j=0;j<nionsp[i];j++){
      fgets(buffer,LINE_SIZE,infile);
      if (frac){
        if (sscanf(buffer,"%lf %lf %lf",m->atoms[k].frac,m->atoms[k].frac+1,
                   m->atoms[k].frac+2)!=3){
          fprintf(stderr,"Error parsing atom %d\n",k);
          exit(1);
        }
      }
      else{
        if (sscanf(buffer,"%lf %lf %lf",m->atoms[k].abs,m->atoms[k].abs+1,
                   m->atoms[k].abs+2)!=3){
          fprintf(stderr,"Error parsing atom %d\n",k);
          exit(1);
        }
      }
      m->atoms[k].atno=species[i];
      m->atoms[k].spin=0;
      m->atoms[k].chg=0;
      m->atoms[k].label=NULL;
      k++;
    }
  }

  /* Sort out structure */

  if (frac) addabs(m->atoms,m->n,c->basis);
  else addfrac(m->atoms,m->n,c->recip);

  /* Next have FFT grid if CHG/CHGCAR,
     or perhaps velocities if POSCAR/CONTCAR */

  fgets(buffer,LINE_SIZE,infile); /* Not sure what this is */
  if(!fgets(buffer,LINE_SIZE,infile)) return;
  
  if (sscanf(buffer,"%d %d %d",fft,fft+1,fft+2)!=3){
    fprintf(stderr,"Unable to find/parse grid data\n");
    return;
  }

  if (!(flags&(CHDEN|SPINDEN))) return;
  
  isize=fft[0]*fft[1]*fft[2];

  /* Data are not stored in the order that we want them */

  if(!(dptr=malloc(isize*sizeof(double))))
    error_exit("Malloc error in vasp_read for temp grid data");

  if (formatted_read(dptr,isize,infile,NULL)) exit(1);

  /* If we wanted the charge density, now rescue onto a grid */

  if (flags&CHDEN){
    if (gptr->next) gptr=gptr->next;
    gptr->next=malloc(sizeof(struct grid));
    if (!gptr->next) error_exit("Malloc error for struct grid");
    gptr->next->next=NULL;
    gptr->next->data=NULL;
    
    if(!(gptr->data=malloc(isize*sizeof(double))))
      error_exit("Malloc error in vasp_read for grid data");
    for(i=0;i<3;i++) gptr->size[i]=fft[i];


  /* Now reverse order, writing consecutively for efficiency */

    dptr2=gptr->data;
    for(i=0;i<fft[0];i++)
      for(j=0;j<fft[1];j++)
        for(k=0;k<fft[2];k++)
          *(dptr2++)=dptr[i+j*fft[0]+k*fft[0]*fft[1]];

  /* Rescale density */
  
    scale=1/fabs(c->vol);
    if (flags&AU) scale*=BOHR*BOHR*BOHR;
    
    if(!(flags&RAW)){
      for(i=0;i<isize;i++) gptr->data[i]*=scale;
      gptr->name="Density";
      if (flags&AU)
        add_cmt(m->comment,"Densities in e Bohr^-3");
      else
        add_cmt(m->comment,"Densities in e A^-3");
    }
    else{
      gptr->name="Density_raw";
      add_cmt(m->comment,"Densities are raw");
    }
  }
  free(dptr);
  dptr=NULL;
  
  /* Wonder what, if anything, else we have */

  /* The only thing we might have, and know how to read, is a spin density */
  
  if (!(flags&SPINDEN)) return;
  
  if (!fgets(buffer,LINE_SIZE,infile)) return;

  /* We might have some augmentation occupancies.
     For the moment, discard these. If norm-conserving pots used, this
     section will be absent. */

  while (!strncasecmp(buffer,"augmentation occupancies ",25)){
    if (sscanf(buffer+25,"%d %d",&i,&j)!=2){
      fprintf(stderr,"Failed to parse augmentation line. Ceasing reading\n");
      return;
    }
    if (debug>2) fprintf(stderr,"Skipping augmentation occupancy line\n");
    dptr=realloc(dptr,j*sizeof(double));
    formatted_read(dptr,j,infile,NULL);
    if (!fgets(buffer,LINE_SIZE,infile)){
      fprintf(stderr,"Spin density not found\n");
      return;
    }
  }

  /* End of augmentation occupancies */

  /* Might now have natom floats to read. These are the magnetic moments
     for the atoms */

  magmom=malloc(natoms*sizeof(double));
  formatted_read(magmom,natoms,infile,buffer);
  for(i=0;i<natoms;i++)
    m->atoms[i].spin=magmom[i];
  free(magmom);

  if (!fgets(buffer,LINE_SIZE,infile)) return;
  
  if (sscanf(buffer,"%d %d %d",ffts,ffts+1,ffts+2)==3){
    if ((ffts[0]==fft[0])&&
        (ffts[1]==fft[1])&&
        (ffts[2]==fft[2])){
      if (debug) fprintf(stderr,"Spin density found\n");
        if (gptr->next) gptr=gptr->next;
        gptr->next=malloc(sizeof(struct grid));
        if (!gptr->next) error_exit("Malloc error for struct grid");
        gptr->next->next=NULL;
        gptr->next->data=NULL;

        isize=fft[0]*fft[1]*fft[2];
        if(!(gptr->data=malloc(isize*sizeof(double))))
          error_exit("Malloc error in vasp_read for grid data");
        for(i=0;i<3;i++) gptr->size[i]=fft[i];

        /* Data are not stored in the order that we want them */

        if(!(dptr=malloc(isize*sizeof(double))))
          error_exit("Malloc error in vasp_read for temp grid data");

        if (formatted_read(dptr,isize,infile,NULL)) exit(1);

        /* Now reverse order, writing consecutively for efficiency */

        dptr2=gptr->data;
        for(i=0;i<fft[0];i++)
          for(j=0;j<fft[1];j++)
            for(k=0;k<fft[2];k++)
              *(dptr2++)=dptr[i+j*fft[0]+k*fft[0]*fft[1]];

        free(dptr);
        dptr=NULL;
  
        /* Rescale density */
  
        scale=1/fabs(c->vol);
        if (flags&AU) scale*=BOHR*BOHR*BOHR;
        
        if(!(flags&RAW)){
          for(i=0;i<isize;i++) gptr->data[i]*=scale;
          gptr->name="Spin";
          if (flags&AU)
            add_cmt(m->comment,"Spin in e Bohr^-3");
          else
            add_cmt(m->comment,"Spin in e A^-3");
        }
        else{
          gptr->name="Spin_raw";
        } 
        elect->nspins=2;
    }
  }
  else
    fprintf(stderr,"Spin density not found\n");
}

static int formatted_read(double *dptr,unsigned int len, FILE *infile,
                          char *first){
  char *ptr,buffer[LINE_SIZE];
  int i,j;
  double dummy;
  
  buffer[0]=0;
  if (len==0){ /* Read blank line */
    if (first) strncpy(buffer,first,LINE_SIZE);
    else fgets(buffer,LINE_SIZE,infile);
    if(sscanf(buffer,"%lf",&dummy)==1){
      fprintf(stderr,"Error reading formatted data in vasp_read."
              " Expected zero items, found one\n%s\n",buffer);
      return 1;
    }
    return 0;
  }
  i=0;
  if (first) strncpy(buffer,first,LINE_SIZE);
  else fgets(buffer,LINE_SIZE,infile);
  while(i<len){
    ptr=buffer;
    while (sscanf(ptr,"%lf%n",dptr+i,&j)==1){
      i++;
      ptr+=j;
    }
    if (i<len)
      if (!fgets(buffer,LINE_SIZE,infile)) break;
  }

  if (i!=len){
    fprintf(stderr,"Error reading formatted data in vasp_read."
            " Expected %d items, found %d\n",len,i);
    return 1;
  }
  return 0;
}
