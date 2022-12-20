/* Read a VASP CHGCAR or CHG, or also a POSCAR or a CONTCAR file
 */

/* Copyright (c) 2017-2019 MJ Rutter 
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

void to235(int *i);

static int formatted_read(double *dptr, unsigned int len, FILE *infile,
                          char *first);

static void vasp_kpoints_read(FILE *infile, struct kpts *k);

char *strrsubs(char *str, char *old, char *new);


void vasp_read(FILE* infile, char *filename, struct unit_cell *c,
               struct contents *m, struct kpts *kp, struct grid *gptr,
               struct es *elect){
  int i,j,k,nspec,natoms,frac,fft[3],ffts[3],isize;
  int *nionsp,*species;
  char buffer[LINE_SIZE+1], *ptr,*p2;
  char *names[]={"POSCAR","CHGCAR","CONTCAR","CHG",NULL};
  double scale,ascale[3],vscale,*dptr,*dptr2,*magmom;
  FILE *f;

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
  if (!strcasecmp(ptr,"comment\n")) *ptr=0;
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
  init_atoms(m->atoms,m->n);
  
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
      k++;
    }
  }

  /* Sort out structure */

  if (frac) addabs(m->atoms,m->n,c->basis);
  else addfrac(m->atoms,m->n,c->recip);

  /* Interlude to wonder about reading kpoints from IBZKPT */

  for(i=0;names[i];i++){
    ptr=names[i];
    if (strstr(filename,ptr)){
      p2=strrsubs(filename,ptr,"IBZKPT");
      if (p2){
        f=fopen(p2,"r");
        if (f) {
          vasp_kpoints_read(f,kp);
          fclose(f);
          if (kp->n){
            fprintf(stderr,"Additionally read %s\n",p2);
            free(p2);
            break;
          }
        }
        free(p2);
      }
    }
  }

  /* Next have FFT grid if CHG/CHGCAR,
     or perhaps velocities if POSCAR/CONTCAR */

  fgets(buffer,LINE_SIZE,infile); /* Not sure what this is */
  if(!fgets(buffer,LINE_SIZE,infile)) return;
  
  if (sscanf(buffer,"%d %d %d",fft,fft+1,fft+2)!=3){
    if (debug>1) fprintf(stderr,"Unable to find/parse grid data\n");
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
     section will be absent, as it also is in CHG rather than CHGCAR
     files. */

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
     for the atoms. But not if we have a CHG not CHGCAR file */

  if (sscanf(buffer,"%d %d %d",ffts,ffts+1,ffts+2)!=3){
    magmom=malloc(natoms*sizeof(double));
    formatted_read(magmom,natoms,infile,buffer);
    for(i=0;i<natoms;i++)
      m->atoms[i].spin=magmom[i];
    free(magmom);
    if (!fgets(buffer,LINE_SIZE,infile)) return;
  }
  
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

/* VASP does not record the plane wave component to grid mapping,
 * but the answer seems to be standard (and dodgy for involving
 * comparisons of doubles).
 *
 * Call with nplwv=0 and pwgrid=NULL for a count of the number of
 * plane waves within sphere.
 *
 * pwgrid must be malloced to 3*nplwv ints
 * gamma!=0 -> form half grid for gamma-point calculations
 */
static int grid_fill(int *pwgrid, double recip[3][3], double kpt[3],
                     double ecut, int fft[3], int nplwv, int gamma){
  int i,j,k,ii,count;
  int nx,ny,nz;
  double rx,ry,rz,v[3],e;
  const double vmagic=3.810019874080794559413826;

  count=0;
  for (i=0;i<fft[2];i++){
    nz=i;
    if (i>fft[2]/2) nz=i-fft[2];
    rz=nz+kpt[2];
    for (j=0;j<fft[1];j++){
      ny=j;
      if (j>fft[1]/2) ny=j-fft[1];
      ry=ny+kpt[1];
      for (k=0;k<fft[0];k++){
        nx=k;
        if (k>fft[0]/2) {
          if (gamma) break;
          else nx=k-fft[0];
        }
        if ((gamma)&&(k==0)&&(ny<0)) continue;
        if ((gamma)&&(k==0)&&(j==0)&&(nz<0)) continue;
        rx=nx+kpt[0];
        for(ii=0;ii<3;ii++)
          v[ii]=2*M_PI*(rx*recip[0][ii]+ry*recip[1][ii]+rz*recip[2][ii]);
        e=vmagic*vmod2(v);
        if (e<ecut){
          if (count<nplwv){
            pwgrid[3*count]=k;
            pwgrid[3*count+1]=j;
            pwgrid[3*count+2]=i;
          }
          count++;
        }
      }
    }
  }
  return count;
}

/* TO DO:
 *        Read grid from CHG/CHGCAR?
 *        Test lots!
 */

void vasp_psi_read(FILE* infile, char *filename, struct unit_cell *c,
                   struct contents *m, struct kpts *k, struct grid *gptr,
                   struct es *elect, int *i_grid){
  int i,j,ns,tmp,nb,ik,n0,n1,n2;
  int nspins,reclen,nkpts,nbands,nplwv,*pwgrid,fft[3],nfft[3],ffft[3];
  int ngpts,offset,gamma,wt_warn,khdr_len,okay;
  double junk,version,ecut,kpt[3],conv,occ,wkpt;
  float *psi_float;
  double *psi,*dptr,basis[3][3];
  const double vmagic=3.810019874080794559413826;
  char *newfile;
  char *files[]={"CHG","CHGCAR","CONTCAR","POSCAR",NULL}; /* Seven char max */
  FILE *f;
  
  fft[0]=fft[1]=fft[2]=0;
  gamma=0;
  wt_warn=1;
  khdr_len=1;
  
  fread(&junk,8,1,infile);
  if (junk!=(int)junk) error_exit("Unable to read WAVECAR");
  reclen=junk;
  if (debug>2) fprintf(stderr,"reclen=%d  (%d doubles)\n",reclen,reclen/8);
  
  fread(&junk,8,1,infile);
  nspins=junk;

  fread(&version,8,1,infile);
  fprintf(stderr,"WAVECAR version %f\n",version);
  if ((version!=45200)&&(version!=53300))
    error_exit("Unknown WAVECAR version");

  
  
  if (m->n==0){ /* We'd like to find some atoms too */
    newfile=malloc(strlen(filename)+4);
    if (!newfile) error_exit("Malloc error for filename");
    for (i=0;files[i];i++){
      newfile=strrsubs(filename,"WAVECAR",files[i]);
      if (!newfile) break;
      f=fopen(newfile,"r");
      if (f){
        fprintf(stderr,"Additionally reading %s\n",newfile);
        vasp_read(f,newfile,c,m,k,gptr,elect);
        fclose(f);
        free(newfile);
        break;
      }
      else if (debug>1) fprintf(stderr,"Tried *%s*\n",newfile);
      free(newfile);
    }
  }
  
  fseek(infile,reclen,SEEK_SET);

  fread(&junk,8,1,infile);
  nkpts=junk;
  fread(&junk,8,1,infile);
  nbands=junk;

  tmp=(3*nbands+4)*8;
  if (tmp>reclen){
    if (version==53300){
      khdr_len=(tmp+reclen-1)/reclen;
      if (debug>2) fprintf(stderr,"khdr_len=%d\n",khdr_len);
    }
    else
      error_exit("Confused by record length and nbands");
  }
  
  fread(&ecut,8,1,infile);
  if (debug>2) fprintf(stderr,"ecut=%f\n",ecut);

  fread(basis,8,9,infile);
  if (c->basis){
    okay=1;
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
        if (!(aeq(c->basis[i][j],basis[i][j]))) okay=0;
    if (!okay){
      fprintf(stderr,
              "Warning: basis in WAVECAR differs from that just read\n");
      fprintf(stderr,"Read:\n");
      print_basis(c->basis);
      fprintf(stderr,"WAVECAR contains\n");
      print_basis(basis);
    }
    free(c->basis); /* DP version in WAVECAR should be more accurate */
  }
  c->basis=malloc(9*sizeof(double));
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[i][j]=basis[i][j];

  real2rec(c);

  for(i=0;i<3;i++)
    if (fft[i]==0){
      fft[i]=4*sqrt(ecut/(4*M_PI*M_PI*vmagic*vmod2(c->recip[i])))+1;
      to235(fft+i);
    }

  if (debug) fprintf(stderr,"Grid chosen: %dx%dx%d\n",fft[0],fft[1],fft[2]);
  
  elect->e_fermi=malloc(sizeof(double));
  fread(elect->e_fermi,8,1,infile);
  
  if ((k->n!=0)&&(k->n!=nkpts)){
    fprintf(stderr,"%d kpts in %s, but %d in IBZKPT!",nkpts,filename,k->n);
    if (k->kpts) free(k->kpts);
    k->kpts=NULL;
    k->n=0;
  }
  if (k->n==0){
    k->n=nkpts;
    k->kpts=malloc(k->n*sizeof(struct atom));
    if (!k->kpts) error_exit("Malloc error for kpts");
    init_atoms(k->kpts,k->n);
  }
  else wt_warn=0;
  
  elect->nbands=nbands;
  elect->nbspins=elect->nspins=nspins;
  tmp=nspins*nkpts*nbands;
  elect->occ=malloc(tmp*sizeof(double));
  elect->eval=malloc(tmp*sizeof(double));
  if ((!elect->occ)||(!elect->eval))
    error_exit("Malloc error for occupations/evals");
  for(ns=0;ns<nspins;ns++){
    for(ik=0;ik<nkpts;ik++){
      fseek(infile,((ns*nkpts+ik)*(nbands+khdr_len)+2)*reclen,SEEK_SET);
      fread(&junk,8,1,infile);
      nplwv=junk;
      fread(kpt,8,3,infile);
      tmp=(ik*nspins+ns)*nbands;
      for(j=0;j<nbands;j++){
        fread(elect->eval+tmp+j,8,1,infile);
        fread(&junk,8,1,infile); /* Imaginary part of energy?! */
        fread(elect->occ+tmp+j,8,1,infile);
      }
      if (wt_warn==0){
        if ((!aeq(kpt[0],k->kpts[ik].frac[0]))||
            (!aeq(kpt[1],k->kpts[ik].frac[1]))||
            (!aeq(kpt[2],k->kpts[ik].frac[2]))){
          fprintf(stderr,"Warning: k-point coords inconsistent with IBZKPT");
          wt_warn=1;
        }
      }
      else k->kpts[ik].wt=1/(double)nkpts;
      for(j=0;j<3;j++) k->kpts[ik].frac[j]=kpt[j];
      pwgrid=malloc(3*nplwv*sizeof(int));
      if (((flags&BANDS)||(flags&BANDPARITY))&&
          (inrange(ik+1,elect->kpt_range))&&
          (inrange(ns,elect->spin_range))){
        tmp=grid_fill(pwgrid, c->recip, kpt, ecut, fft, nplwv,gamma);
        if (tmp!=nplwv){
          /* VASP fails to store the Gamma point kpoint as precisely zero */
          if ((tmp==2*nplwv-1)&&(vmod2(kpt)<1e-20)){
            gamma=1;
            tmp=grid_fill(pwgrid, c->recip, kpt, ecut, fft, nplwv,gamma);
            if (tmp!=nplwv){
              fprintf(stderr,
                      "Unexpected number of components in gamma mapping: "
                  "skipping kpt %d\n",ik+1);
              continue;
            }
            for(i=0;i<3;i++) k->kpts[ik].frac[i]=kpt[i]=0;
          }
          else{
            fprintf(stderr,"Unexpected number of plane wave components: "
                    "skipping kpt %d\n",ik+1);
            fprintf(stderr,"Expected %d, mapping finds %d\n",nplwv,tmp);
            continue;
          }
        }
        if (debug>2) fprintf(stderr,"nplwv=%d\n",nplwv);
        for(nb=0;nb<nbands;nb++){
          if (inrange(nb+1,elect->band_range)){
            if (debug>2) fprintf(stderr,"Reading band %d\n",nb+1);
            fseek(infile,
                  ((ns*nkpts+ik)*(nbands+khdr_len)+nb+2+khdr_len)*reclen,
                  SEEK_SET);
            psi_float=malloc(2*nplwv*sizeof(float));
            if (!psi_float) error_exit("Malloc error for psi_float");
            tmp=fread(psi_float,4,2*nplwv,infile);
            if (tmp!=2*nplwv) error_exit("Short read for psi_float");
            fprintf(stderr,"psi_float[0]=(%f,%f)\n",psi_float[0],psi_float[1]);
            ngpts=fft[0]*fft[1]*fft[2];
            psi=malloc(2*ngpts*sizeof(double));
            if (!psi) error_exit("Malloc error for psi");
            for(i=0;i<2*ngpts;i++) psi[i]=0;
            for(i=0;i<nplwv;i++){
              offset=pwgrid[3*i+2]+fft[2]*(pwgrid[3*i+1]+
                                            fft[1]*pwgrid[3*i]);
              if ((offset<0)||(offset>ngpts)){
                fprintf(stderr,"Impossible offset in wave_read off=%d i=%d\n",
                        offset,i);
                exit(1);
              }
              psi[2*offset]=psi_float[2*i];
              psi[2*offset+1]=psi_float[2*i+1];
            }
            if (gamma){ /* construct psi(-k)=conjg(psi(k)) */
              for(i=0;i<nplwv;i++){
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
                psi[2*offset]=psi_float[2*i];
                psi[2*offset+1]=-psi_float[2*i+1];
              }
            }
            free(psi_float);

            if (debug>2) fprintf(stderr,"Before FFT g=0 component is %g+%gi\n",
                                 psi[0],psi[1]);

            if (((kpt[0]==0)||aeq(fabs(kpt[0]),0.5))&&
                ((kpt[1]==0)||aeq(fabs(kpt[1]),0.5))&&
                ((kpt[2]==0)||aeq(fabs(kpt[2]),0.5))&&
                (flags&BANDPARITY)) inv_parity(psi,fft,nb,kpt);

            if (!(flags&BANDS)) {
              free(psi);
              continue;
            }

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
            else
             for(i=0;i<3;i++) nfft[i]=fft[i]; 
            
            ffft[0]=nfft[2];
            ffft[1]=nfft[1];
            ffft[2]=nfft[0];
            fft3d(psi,ffft,1);

            dptr=malloc(nfft[0]*nfft[1]*nfft[2]*sizeof(double));
            if (!dptr) error_exit("Malloc error for grid data");
            band2real(psi,dptr,nfft,kpt);
            free(psi);

            /* Do we need to rescale? */
            if (((flags&RAW)==0)&&((flags&BANDPHASE)==0)){ /* Yes */
              if (flags&BANDDEN) conv=1/c->vol;
              else conv=1/sqrt(c->vol);
              if (debug>2) fprintf(stderr,"Scaling wavefun by %f\n",conv);
              for(i=0;i<nfft[0]*nfft[1]*nfft[2];i++) dptr[i]*=conv;
            }

            if (elect->occ)
              occ=elect->occ[elect->nspins*nb*ik+ns*nbands+nb];
            else
              occ=1;
            wkpt=k->kpts[ik].wt;
            band_store(&gptr,dptr,occ,wkpt,1,ns,ik+1,nb+1,elect,m,nfft);

          } /* inrange(band) */

        } /* End band loop */
      } /* inrange(kpt) */
      free(pwgrid);
    } /* End kpt loop */

  } /* End spin loop */

  if ((flags&K_WEIGHT)&&(wt_warn))
    fprintf(stderr,"Warning: kpoint weights not read, "
            "weights likely to be wrong\n");

  
  if (debug){
    print_elect(elect);
  }

  if (debug>2)
    fprintf(stderr,"Final offset %d  reclen=%d\n",(int)ftell(infile),reclen);
  
}

/* Replace right-most occurance of old with new in str.
 * Return malloced pointer to new string
 */
char *strrsubs(char *str, char *old, char *new){
  char *cptr,*rtn,*cptr2;

  //  fprintf(stderr,"strrsubs called *%s* *%s* *%s*\n",str,old,new);
  
  if (!strstr(str,old)) return NULL;

  rtn=malloc(strlen(str)+strlen(new)-strlen(old)+1);
  if (!rtn) error_exit("Malloc error in strrsubs");

  cptr=str;
  while(strstr(cptr+1,old)) cptr=strstr(cptr+1,old); /* Find last occurance */

  strncpy(rtn,str,cptr-str); /* Copy part before old */
  cptr2=rtn+(cptr-str);

  strncpy(cptr2,new,strlen(new)); /* Copy new */
  cptr2+=strlen(new);
  cptr+=strlen(old);        /* No nulls copied */

  /* Copy tail of original string */
  while(*cptr) {
    *(cptr2++)=*(cptr++);
  }

  *cptr2=0; /* And final null */
  
  //  fprintf(stderr,"strrsubs returns *%s*\n",rtn);
  
  return rtn;
  
}

/* Read IBZKPT format only */
static void vasp_kpoints_read(FILE *infile, struct kpts *k){
  char buffer[LINE_SIZE+1];
  double total;
  int i;
  
  /* Skip first line */
  if (!fgets(buffer,LINE_SIZE,infile)) return;

  if (!fgets(buffer,LINE_SIZE,infile)) return;
  if (sscanf(buffer,"%d",&k->n)!=1) return;

  fgets(buffer,LINE_SIZE,infile);
  if (buffer[0]!='R'){
    fprintf(stderr,"Unsupported IBZKPT format");
    k->n=0;
    return;
  }

  k->kpts=malloc(k->n*sizeof(struct atom));
  if (!k->kpts) error_exit("Malloc error for kpts");
  init_atoms(k->kpts,k->n);
  for(i=0;i<k->n;i++){
    fgets(buffer,LINE_SIZE,infile);
    if (sscanf(buffer,"%lf %lf %lf %lf",k->kpts[i].frac,
               k->kpts[i].frac+1,k->kpts[i].frac+2,&k->kpts[i].wt)!=4){
      fprintf(stderr,"Parse error for kpt");
      free(k->kpts);
      k->kpts=NULL;
      k->n=0;
    }
  }

  /* Need to sort out unnormalised weights */

  total=0;
  for(i=0;i<k->n;i++) total+=k->kpts[i].wt;
  for(i=0;i<k->n;i++) k->kpts[i].wt/=total;
  
}
