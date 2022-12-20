/* Read some useful data from a CASTEP .check file
 *
 * Cope with either endianness
 *
 * Make various assumptions about the format...
 *
 * This version understands of spinors
 */


/* Copyright (c) 2007-2017 MJ Rutter 
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
#include<string.h> /* memcpy */
#include<errno.h>
#include<math.h>   /*sqrt*/

#include "c2xsf.h"

/* HEAD_LEN must be >=30 */
#define HEAD_LEN 30
#define CBUFF 80

/* scan for token x, read data of length len_x to targ, allocating
 * if required. Reverse endian if len_x==4
 */
#define SCAN(x,len_x,targ) \
  if (match(head,x)){ \
    fread(&tmp,4,1,infile); \
    if (endian) reverse4(&tmp); \
    if (tmp!=len_x) error_exit("Unexpected data in " #x); \
    if ((!targ)&&!(targ=malloc(len_x))) error_exit("Malloc error in " #x); \
    fread(targ,len_x,1,infile); \
    if (endian && (len_x==4)) reverse4(targ); \
    fseek(infile,4,SEEK_CUR); \
    if (debug>2) fprintf(stderr,"Found %s, %d bytes\n", #x, len_x); \
  }

/* scan for token x, read data of length len_x to targ, allocating
 * if required. Reverse endian assuming 8 byte tokens
 */
#define SCAN8R(x,len_x,targ) \
  if (match(head,x)){ \
    fread(&tmp,4,1,infile); \
    if (endian) reverse4(&tmp); \
    if (tmp!=len_x) error_exit("Unexpected data in " #x); \
    if ((!targ)&&!(targ=malloc(len_x))) error_exit("Malloc error in " #x); \
    fread(targ,len_x,1,infile); \
    if (endian) reverse8n(targ,len_x/8);		\
    fseek(infile,4,SEEK_CUR); \
    if (debug>2) fprintf(stderr,"Found %s, %d bytes\n", #x, len_x); \
  }


int match(char *s1, char *s2);
int inrange(int x, char *range);
void reverse4(void *data);
void reverse8(void *data);
void reverse4n(int *data,int n);
void reverse8n(double *data,int n);
static void wave_read(FILE *infile, int nkpts, double *kpts_pos,
                      double *wkpt, int gamma, struct grid *g,
                      struct es *elect, int *i_grid,struct contents *m);
static void occ_read(FILE *infile, double* occ, double *kpts_pos, int nkpts,
                     int nspins, int nbands);
static int find_kpt(double *kpt, double *kpts_pos, int nkpts);
void fft3d(double *c, int *ngptar, int dir);
void pad_recip(double *o, int fft[3], double **nptr, int nfft[3]);

/* variables shared with wavefunction reader */

static int endian;
static int nbands;
static double *cr_cell_vol;

void check_read(FILE* infile, struct unit_cell *c, struct contents *m,
                struct kpts *kp, struct symmetry *s, struct grid *gptr,
                struct es *elect, int *i_grid){
  int tmp;
  char head[HEAD_LEN+1];
  char cbuff[CBUFF+1];
  int i,j,k,density,na,fft[3],ion_sym_len,ps_pot_len,ver_maj,ver_min;
  int smisc_size;
  double *dptr1,*dptr2,*column,*castep_basis,conv;
  int *nsp,*nsp_max,*nionsp,*nsymops;
  char *ion_sym,*ps_pots;
  double *ion_pos,*ion_force,*mag_mom,*sym_mat,*sym_disp;
  double *species_charge;
  double *kpts_pos,*wkpts;
  int *mp_xyz;
  int *nkpts;
  double *mp_off;
  struct grid *gptr2;
  int section;
  char *spin_method;
  int nspins;
  double cut,x;
  char *sym,*lab,*cptr,*cptr2;
  char dipole_correction[21];
  char dipole_dir;
  
  ion_sym_len=8;
  endian=0;
  castep_basis=NULL;
  ion_pos=ion_force=mag_mom=NULL;
  nsp=nsp_max=nionsp=NULL;
  ps_pots=NULL;
  kpts_pos=NULL;
  nkpts=NULL;
  wkpts=mp_off=NULL;
  mp_xyz=NULL;
  s->tol=NULL;
  nsymops=NULL;
  sym_disp=NULL;
  sym_mat=NULL;
  ion_sym=NULL;
  species_charge=NULL;
  gptr2=NULL;
  section=1;
  ver_maj=0;
  nspins=1;

  if(debug>2) fprintf(stderr,"check_read called with flags=%d\n",flags);

  /* The first record is a string of length 30. Being Fortran, the first
   * item will therefore be an integer, 4 bytes, of value 30. If we
   * have an endian problem, we will see this as an integer of value
   * 30*(1<<24)
   *
   * Unless we have the new "castep_bin" file format, in which case
   * the first record has a length of 10
   */

  fread(&tmp,4,1,infile);

  if ((tmp!=30)&&(tmp!=10)){ /* Oh dear, wrong magic */
    if ((tmp==(30*(1<<24)))||(tmp==(10*(1<<24)))){
      endian=1;
    }else{
      fprintf(stderr,"Error: unable to read input as a CASTEP .check file.\n");
      exit(1);
    }
  }

  rewind(infile);
  density=0;

/* FORTRAN stores records as  int32     length
 *                           <length>   data
 *                            int32     length  (i.e. repeated of first)
 *
 * consider check file as five sections: start to 1st END_CELL_GLOBAL
 *                                       rest to the start of wavefuncts
 *                                       wavefunctions
 *                                       charge (& spin) density
 *                                       forces
 *
 * note that the cell appears twice, and we want the first copy
 * note that the wavefunctions will be absent in a .castep_bin file
 */

  while(fread(&tmp,4,1,infile)){
    if (endian) reverse4(&tmp);
    i=(tmp<HEAD_LEN)? tmp : HEAD_LEN;
    fread(head,i,1,infile);
    head[i]=0;
    if (debug>3) fprintf(stderr,"%d: %.4s\n",tmp,head);
    if (i<tmp) fseek(infile,tmp-i,SEEK_CUR);
    fseek(infile,4,SEEK_CUR); /* capricious */

    /* Need to find version of file so that we know how long the ionic
     * species symbol records will be. Version is stored in ASCII
     * as major.minor immediately after the BEGIN_PARAMETERS_DUMP
     * marker
     */

    switch (section){
    case(1):
      if(match(head,"BEGIN_PARAMETERS_DUMP")){
        fread(&tmp,4,1,infile);
        if (endian) reverse4(&tmp);
        i=(tmp<HEAD_LEN)? tmp : HEAD_LEN;
        fread(head,i,1,infile);
        fseek(infile,4,SEEK_CUR);
        head[i]=0;
        j=0;
        ver_maj=0;
        /* Swallow leading spaces */
        while(head[j]==' ') j++;
        while((head[j])&&(head[j]!='.')&&(j<i))
          ver_maj=10*ver_maj+head[j++]-'0';
        if (head[j]=='.') j++;
        ver_min=0;
        while((head[j])&&(head[j]!=' ')&&(j<i))
          ver_min=10*ver_min+head[j++]-'0';
        if(ver_maj<4) ion_sym_len=3;
        if (debug) fprintf(stderr,"Castep version %d.%d\n",ver_maj,ver_min);
      }

      SCAN("CELL%REAL_LATTICE",72,castep_basis);
      SCAN("CELL%NUM_SPECIES",4,nsp);
      SCAN("CELL%MAX_IONS_IN_SPECIES",4,nsp_max);
      SCAN("CELL%NUM_IONS_IN_SPECIES",(4*(*nsp)),nionsp);
      SCAN("CELL%IONIC_POSITIONS",(24*(*nsp)*(*nsp_max)),ion_pos);
      SCAN("CELL%SPECIES_SYMBOL",(ion_sym_len*(*nsp)),ion_sym);
      SCAN("CELL%INITIAL_MAGNETIC_MOMENT",(8*(*nsp)*(*nsp_max)),mag_mom);
      SCAN("CELL%IONIC_CHARGE_REAL",(8*(*nsp)),species_charge);
      SCAN8R("CELL%VOLUME",8,cr_cell_vol);
      SCAN("NKPTS",4,nkpts);
      SCAN8R("KPOINTS",24*(*nkpts),kpts_pos);
      SCAN8R("KPOINT_WEIGHTS",8*(*nkpts),wkpts);
      SCAN("SYMMETRY_GENERATE",4,s->gen);
      SCAN("SYMMETRY_TOL",8,s->tol);
      SCAN("KPOINT_MP_GRID",12,mp_xyz);
      SCAN("KPOINT_MP_OFF",24,mp_off);
      SCAN("NUM_SYMMETRY_OPERATIONS",4,nsymops);
      SCAN("SYMMETRY_OPERATIONS",72*(*nsymops),sym_mat);
      SCAN("SYMMETRY_DISPS",24*(*nsymops),sym_disp);
      if (match(head,"CELL%SPECIES_POT")){ /* We don't really know */
	fread(&tmp,4,1,infile);            /* the length of this record */
	if (endian) reverse4(&tmp);
	ps_pots=malloc(tmp+1);
	if (!ps_pots) error_exit("Malloc error for ps_pots");
	fread(ps_pots,tmp,1,infile);
	ps_pots[tmp]=0; /* .check file will space-pad, we want null term */
	fseek(infile,4,SEEK_CUR);
	if (debug>2) fprintf(stderr,"Found CELL%%SPECIES_POT, %d bytes\n",tmp);
      }
      if (match(head,"BEGIN_ELECTRONIC")){
        fseek(infile,4*12+4,SEEK_CUR);
        fread(&nbands,4,1,infile);
        if (endian) reverse4(&nbands);
        /* From Castep 16 nspins is always zero,
            so read spin_polarized instead */
        fseek(infile,4+16+2*12+4,SEEK_CUR);
        fread(&tmp,4,1,infile);  /* This is a logical */
        fseek(infile,4,SEEK_CUR);
        if (tmp==0) elect->nspins=1;
        else elect->nspins=2;
        if (debug>1) fprintf(stderr,"Spins=%d  nbands=%d\n",
			     elect->nspins,nbands);

        fseek(infile,18+5*16 ,SEEK_CUR);

        if (ver_maj>=9){
          fread(&tmp,4,1,infile);
          if (endian) reverse4(&tmp);
          if (tmp!=20) error_exit("Unexpected field length cr_sm\n");
          spin_method=malloc(21);
          elect->spin_method=spin_method;
          fread(spin_method,20,1,infile);
          spin_method[20]=0; /* null terminate string */
          if (debug>2) fprintf(stderr,"Spin method: %s\n",spin_method);
          fseek(infile,4,SEEK_CUR);
          if (!strncmp(spin_method,"NONE",4)){
	    if (elect->nspins!=1)
	      error_exit("Spin method none but nspins!=1\n");
	  }
          if (!strncmp(spin_method,"SCALAR",6)){
	    if (elect->nspins!=2)
	      error_exit("Spin method scalar but nspins!=2\n");
	  }
          if (!strncmp(spin_method,"VECTOR",6)){
	    if (elect->nspins!=1)
	      error_exit("Spin method vector but nspins!=1\n");
            elect->nspinors=2;
            if (flags&SPINDEN){
              fprintf(stderr,"Unable to output vector spin density\n"
		      " resetting flag requesting spin density output\n");
	      flags^=SPINDEN;
	    }
	  }

        }
	/* Use nspins for charge density
         * By setting as follows, the record length will be
         * nspins*16*fft[2] */
        nspins=1;
	if (elect->nspins==2) nspins=2;
	if (elect->nspinors==2) nspins=4;
      }
      if (match(head,"BEGIN_BASIS")){
        /* Skip first record */
        fread(&tmp,4,1,infile);
        if (endian) reverse4(&tmp);
        fseek(infile,4+tmp,SEEK_CUR);
        fread(&tmp,4,1,infile);
        if (endian) reverse4(&tmp);
        if (tmp==8){
          fread(&cut,8,1,infile);
          if (endian) reverse8(&cut);
          elect->cut_off=cut*H_eV;
          fseek(infile,4,SEEK_CUR);
        }
        else fseek(infile,4+tmp,SEEK_CUR);
        
      }
      if (match(head,"BEGIN_ELEC_MIN")){
        if ((ver_maj>6)||((ver_maj==6)&&(ver_min>=100))){
          /* We have at 18 entries to skip */
          for(i=0;i<18;i++){
            fread(&tmp,4,1,infile);
            if (endian) reverse4(&tmp);
	    if ((i==4)&&(tmp==8)){
	      fread(&x,8,1,infile);
	      if (endian) reverse8(&x);
	      elect->etol=x*H_eV;
	      tmp=0;
	    }
            fseek(infile,4+tmp,SEEK_CUR);
          }
          fread(&tmp,4,1,infile);
          if (endian) reverse4(&tmp);
          if (tmp==20){
            fread(dipole_correction,20,1,infile);
            dipole_correction[20]=0;
            elect->dip_corr=malloc(21);
            if (!elect->dip_corr) error_exit("Malloc err for char!");
            strncpy(elect->dip_corr,dipole_correction,21);
            fseek(infile,4,SEEK_CUR);
            fread(&tmp,4,1,infile);
            if (endian) reverse4(&tmp);
            if (tmp==1){
              fread(&dipole_dir,1,1,infile);
              fseek(infile,4,SEEK_CUR);
              /* Don't over-ride command-line specified dip_corr_dir */
              if ((elect->dip_corr_dir==NULL)&&(dipole_correction[0]!='N')){
                elect->dip_corr_dir=malloc(1);
                if (!elect->dip_corr_dir) error_exit("Malloc err for char!");
                *elect->dip_corr_dir=dipole_dir;
              }
              if ((debug&&(dipole_correction[0]!='N'))||(debug>1))
                fprintf(stderr,"Dipole correction: %s (direction %c)\n",
                        dipole_correction,dipole_dir);
            }
          }
        }
      }
      if (match(head,"END_CELL_GLOBAL")) {
        /* We have an urgent need for the volume so that densities can
         * be scaled as soon as they are read
         */
	if (!(flags&AU))
	  (*cr_cell_vol)*=BOHR*BOHR*BOHR; /* CASTEP uses Bohrs,
					     we want Angstroms */
        section=2;
      }
      break;
    case 2:      /* matches select(section) */
      if (match(head,"END_CELL_GLOBAL")){
        double energy;
        int wave_ok,den_ok;
        fseek(infile,4,SEEK_CUR);
        fread(&wave_ok,4,1,infile);
        fseek(infile,8,SEEK_CUR);
        fread(&den_ok,4,1,infile);
        fseek(infile,8,SEEK_CUR);
        fread(&energy,8,1,infile);
        if (endian) reverse8n(&energy,1);
        if (debug>1) fprintf(stderr,"Total Energy %.6f eV\n",energy*H_eV);
        elect->energy=malloc(sizeof(double));
        if (!elect->energy) error_exit("Malloc error for double!");
        *elect->energy=energy*H_eV;
        fseek(infile,8,SEEK_CUR);
        fread(&energy,8,1,infile);
        if (endian) reverse8n(&energy,1);
        if (debug>1) fprintf(stderr,"Fermi Energy %.6f eV\n",energy*H_eV);
        elect->e_fermi=malloc(sizeof(double));
        if (!elect->e_fermi) error_exit("Malloc error for double!");
        *elect->e_fermi=energy*H_eV;
        fseek(infile,4,SEEK_CUR);
        section=3;
        if(!nkpts) error_exit("Aborting: nkpts not set");
      }
      break;
    case 3: case 4: /* matches select(section) */
      if (match(head,"wave")||match(head,"Gpt")){
	if (match(head,"wave")){
	  if ((flags&BANDS)||(flags&OCCUPANCIES)||(flags&BANDPARITY))
            wave_read(infile, *nkpts, kpts_pos, wkpts, 0,
                      gptr, elect, i_grid, m);
	}
	if (match(head,"Gpt")){
	  if ((flags&BANDS)||(flags&OCCUPANCIES)||(flags&BANDPARITY))
            wave_read(infile, *nkpts, kpts_pos, wkpts, 1,
                      gptr, elect, i_grid, m);
	}
	section=4;
	if (flags&BANDS){
	  if (flags&RAW)
	    add_cmt(m->comment,"Wavefunction unscaled");
	  else{
            if (flags&AU){
              if ((flags&BANDDEN)==BANDDEN)
                add_cmt(m->comment,"Wavefunction as density in e Bohr^-3");
              else
                add_cmt(m->comment,"Wavefunction in Bohr^-1.5");
            }
            else{
              if ((flags&BANDDEN)==BANDDEN)
                add_cmt(m->comment,"Wavefunction as density in e A^-3");
              else
                add_cmt(m->comment,"Wavefunction in A^-1.5");
            }
            if ((flags&OCC_WEIGHT)||(flags&K_WEIGHT)){
              strcpy(cbuff,"Wavefunction weighted by ");
              if (flags&K_WEIGHT){
                strcat(cbuff,"kpt weight");
                if (flags&OCC_WEIGHT) strcat(cbuff," and ");
              }
              if (flags&OCC_WEIGHT) strcat(cbuff," occupancy");
              add_cmt(m->comment,cbuff);
            }
          }
          snprintf(cbuff,CBUFF,"kpt range: %s",elect->kpt_range);
          add_cmt(m->comment,cbuff);
          snprintf(cbuff,CBUFF,"band range: %s",elect->band_range);
          add_cmt(m->comment,cbuff);
          if (elect->nspins>1){
            snprintf(cbuff,CBUFF,"spin range: %s",elect->spin_range);
            add_cmt(m->comment,cbuff);
          }
          for (i=1;i<*nkpts;i++){
            if(inrange(i,elect->kpt_range)){
              snprintf(cbuff,CBUFF,"kpt %3d is (% f, % f, % f)",
                       i,kpts_pos[(i-1)*3],kpts_pos[(i-1)*3+1],
                       kpts_pos[(i-1)*3+2]);
              add_cmt(m->comment,cbuff);
            }
          }
	}
      }
      if (density==0){
        if (tmp==12){ /* we might have an FFT grid */
          memcpy(fft,head,12);
          fread(&tmp,4,1,infile);
          if (endian) {
            reverse4(&tmp);
            reverse4(fft);
            reverse4(fft+1);
            reverse4(fft+2);
          }
          if (debug>2) fprintf(stderr,"Potential FFT grid     %d %d %d\n",
                 fft[0],fft[1],fft[2]);
          if (tmp==nspins*16*fft[2]+8){ /* it might well be an FFT grid
               It will be stored as complex, we will read a complex column,
               then store as real */
            if (debug>2)
              fprintf(stderr,"We think we have found a charge density grid\n");

            if (flags&CHDEN){
              if (gptr->next) gptr=gptr->next;
	      gptr->next=malloc(sizeof(struct grid));
	      if (!gptr->next) error_exit("Malloc error for struct grid");
	      gptr->next->next=NULL;
	      gptr->next->data=NULL;
              gptr->name="Density";
              for(i=0;i<3;i++) gptr->size[i]=fft[i];
              gptr->data=malloc(8*fft[0]*fft[1]*fft[2]);
              if (!gptr->data){
                fprintf(stderr,"Error allocating density grid\n");
                exit(1);
              }
            } /* End if (flags&CHDEN) */
            if ((nspins==2)&&(flags&SPINDEN)){
              gptr2=gptr;
              if (gptr2->next) gptr2=gptr2->next;
	      gptr2->next=malloc(sizeof(struct grid));
	      if (!gptr2->next) error_exit("Malloc error for struct grid");
	      gptr2->next->next=NULL;
	      gptr2->next->data=NULL;
              gptr2->name="Spin";
              for(i=0;i<3;i++) gptr2->size[i]=fft[i];
              gptr2->data=malloc(8*fft[0]*fft[1]*fft[2]);
              if (!gptr2->data){
                fprintf(stderr,"Error allocating density grid\n");
                exit(1);
              }
            }
            if ((flags&CHDEN)||(flags&SPINDEN)){
              column=malloc(16*fft[2]);
              while(tmp==(nspins*16*fft[2]+8)){
                fread(&i,4,1,infile);
                if (endian) reverse4(&i);
                fread(&j,4,1,infile);
                if (endian) reverse4(&j);
                if ((i<1)||(i>fft[0])||(j<1)||(j>fft[1])){
		  fprintf(stderr,"i=%d j=%d\n",i,j);
		  error_exit("Unexpected FFT column in density\n");
		}
                if (flags&CHDEN){
                  fread(column,16*fft[2],1,infile);
                  dptr1=gptr->data+((i-1)*fft[1]+(j-1))*fft[2];
                  dptr2=column;
                  for(k=0;k<fft[2];k++){*dptr1++=*dptr2;dptr2+=2;}
                }
                else fseek(infile,16*fft[2],SEEK_CUR);
                if (nspins==2){
                  if (flags&SPINDEN){
                    fread(column,16*fft[2],1,infile);
                    dptr1=gptr2->data+((i-1)*fft[1]+(j-1))*fft[2];
                    dptr2=column;
                    for(k=0;k<fft[2];k++){*dptr1++=*dptr2;dptr2+=2;}
                  }
                  else fseek(infile,16*fft[2],SEEK_CUR);
                }
		else if (nspins==4)
		  fseek(infile,48*fft[2],SEEK_CUR);
                fseek(infile,4,SEEK_CUR);
                fread(&tmp,4,1,infile);
                if (endian) reverse4(&tmp);
              } /* end while(tmp==(spins... ) */
              free(column);
              fseek(infile,-4,SEEK_CUR);
              /* Correct endianness */
              if (endian){
                if (flags&CHDEN) reverse8n(gptr->data,fft[0]*fft[1]*fft[2]);
                if (flags&SPINDEN) reverse8n(gptr2->data,fft[0]*fft[1]*fft[2]);
              }
              /* Scale */
              if (((flags&RAW)==0)&&((gptr)||(gptr2))){
                conv=1/(*cr_cell_vol);
		if (flags&AU) add_cmt(m->comment,"Densities in e Bohr^-3");
		else add_cmt(m->comment,"Densities in e A^-3");
                if (debug>2)
                  fprintf(stderr,"Rescaling densities by %f\n",conv);
                if (flags&CHDEN)
                  for(i=0;i<fft[0]*fft[1]*fft[2];i++) gptr->data[i]*=conv;
                if (flags&SPINDEN)
                  for(i=0;i<fft[0]*fft[1]*fft[2];i++) gptr2->data[i]*=conv;
              }
	      if ((flags&RAW)&&((gptr)||(gptr2)))
		add_cmt(m->comment,"Densities unscaled, e per unit cell");
              if (gptr2) gptr=gptr2;
            } /* end if ((flags&CHDEN)||(flags&SPINDEN)) */
            else fseek(infile,(nspins*16*fft[2]+16)*fft[1]*fft[0]-4,
		       SEEK_CUR);
            density=1;
            section=5;
          }else{
            fseek(infile,4+tmp,SEEK_CUR);
	    /*  fprintf(stderr,"Skipped %d\n",tmp);  */
          }
        }
      }
      break;
    case 5:  
        SCAN("FORCES",24*(*nsp)*(*nsp_max),ion_force);
        break;
    } /* end select(section) */
  } /* end while(fread(&tmp,4...) */

/* Now sort out the endianness of everything we've read */

  if (endian){
    reverse8n(castep_basis,9);
/* remember single 4 byte objects will have been converted anyway */
    if (*nsp>1) reverse4n(nionsp,*nsp);
    reverse8n(ion_pos,3*(*nsp)*(*nsp_max));
    if (ion_force) reverse8n(ion_force,3*(*nsp)*(*nsp_max));
    if (mag_mom) reverse8n(mag_mom,(*nsp)*(*nsp_max));
    if (species_charge) reverse8n(species_charge,*nsp);
    /* The following two are both reversed when read in */
    //    if (kpts_pos) reverse8n(kpts_pos,3*(*nkpts));
    //    if (wkpts) reverse8n(wkpts,*nkpts);
    if (sym_mat) reverse8n(sym_mat,9*(*nsymops));
    if (sym_disp) reverse8n(sym_disp,3*(*nsymops));
  }

  /* Rescue k-points into sane structure */
  if ((kpts_pos)&&(wkpts)){
    if(!(kp->kpts=malloc((*nkpts)*sizeof(struct atom))))
      error_exit("Malloc error for kpts");
    for(i=0;i<*nkpts;i++){
      for(j=0;j<3;j++) kp->kpts[i].frac[j]=kpts_pos[3*i+j];
      kp->kpts[i].wt=wkpts[i];
    }
    kp->n=*nkpts;
    free(wkpts);
    free(kpts_pos);
  }

  if((mp_xyz)&&(mp_off)&&(mp_xyz[0])&&(mp_xyz[1])&&(mp_xyz[2])){
    if(!(kp->mp=malloc(sizeof(struct mp_grid))))
        error_exit("Malloc error for mp");
    /* Tests for zero are valid before endian reversal */
    if (endian) {reverse4n(mp_xyz,3);reverse8n(mp_off,3);}
    kp->mp->grid[0]=mp_xyz[0];
    kp->mp->grid[1]=mp_xyz[1];
    kp->mp->grid[2]=mp_xyz[2];
    kp->mp->disp[0]=mp_off[0];
    kp->mp->disp[1]=mp_off[1];
    kp->mp->disp[2]=mp_off[2];
    free(mp_off);
    free(mp_xyz);
  }


/* change basis from Fortran ordering to C ordering */

  if(!(c->basis=malloc(72))) error_exit("Error in malloc for basis");
  for(i=0;i<=2;i++)
    for(j=0;j<=2;j++) c->basis[i][j]=castep_basis[3*j+i];
  free(castep_basis);
  
  /* Convert basis to Angstoms */
  for(i=0;i<3;i++)
    for(j=0;j<3;j++) c->basis[i][j]=c->basis[i][j]*BOHR;

/* Add reciprocal basis and volume */

  real2rec(c);

  conv=*cr_cell_vol;
  if (flags&AU) conv*=BOHR*BOHR*BOHR;
  if ((fabs(c->vol-conv)/c->vol)>0.005*fabs(c->vol))
    fprintf(stderr,"Warning: calc cell vol %f A^3, CASTEP cell vol %f A^3\n",
            c->vol,conv);

  /* Rescue symmetry operations into sane structure */
  if ((nsymops)&&(*nsymops>1)){
    s->n=*nsymops;
    s->ops=malloc(s->n*sizeof(struct sym_op));
    if (!s->ops) error_exit("Malloc error for symops");
    for(i=0;i<s->n;i++){
      /* Translations are relative */
      s->ops[i].tr=malloc(3*sizeof(double));
      for(j=0;j<3;j++) s->ops[i].tr[j]=sym_disp[3*i]*c->basis[0][j]+
			 sym_disp[3*i+1]*c->basis[1][j]+
			 sym_disp[3*i+2]*c->basis[2][j];
      /* But rotations are absolute */
      for(j=0;j<3;j++)
	for(k=0;k<3;k++)
	  s->ops[i].mat[j][k]=sym_mat[9*i+3*j+k];
    }
    free(sym_disp);
    free(sym_mat);
  }

/* Pack ions */

  m->n=0;
  for(i=0;i<*nsp;i++) m->n+=nionsp[i];
  if (debug>1) fprintf(stderr,"%d atoms found\n",m->n);
  if (debug>2)
    for(i=0;i<*nsp;i++)
      fprintf(stderr,"Species %d atoms %d\n",i,nionsp[i]);

  if (!(m->atoms=malloc(m->n*sizeof(struct atom))))
    error_exit("Malloc error in check_read");

  init_atoms(m->atoms,m->n);

  na=0;
  for(i=0;i<*nsp;i++){
    sym=ion_sym+ion_sym_len*i;
    if ((sym[1]==':')||(sym[2]==':')||(sym[3]==':')){
      cptr=sym;
      while(*cptr!=':') cptr++;
      lab=malloc(ion_sym_len+1);
      j=0;
      while((cptr<ion_sym+ion_sym_len*(i+1))&&(*cptr!=' '))
        lab[j++]=*(cptr++);
      lab[j]=0;
    }
    else lab=NULL;
    for(j=0;j<nionsp[i];j++){
      if (mag_mom)
        m->atoms[na].spin=mag_mom[i*(*nsp_max)+j];
      else
        m->atoms[na].spin=0;
      if (species_charge)
        m->atoms[na].chg=species_charge[i];
      else
        m->atoms[na].chg=0;
      m->atoms[na].atno=atsym2no(sym);
      m->atoms[na].label=lab;
      if((!m->atoms[na].atno)&&(debug))
        fprintf(stderr,"Warning: atom symbol %s converted to 0\n",
                ion_sym+ion_sym_len*i);
      m->atoms[na].frac[0]=ion_pos[3*(i*(*nsp_max)+j)];  
      m->atoms[na].frac[1]=ion_pos[3*(i*(*nsp_max)+j)+1];  
      m->atoms[na].frac[2]=ion_pos[3*(i*(*nsp_max)+j)+2];  
      if (ion_force){
        m->atoms[na].force[0]=ion_force[3*(i*(*nsp_max)+j)];  
        m->atoms[na].force[1]=ion_force[3*(i*(*nsp_max)+j)+1];  
        m->atoms[na].force[2]=ion_force[3*(i*(*nsp_max)+j)+2];  
      }
      else for(k=0;k<3;k++)m->atoms[na].force[k]=0;
      na++;
    }
  }

  /* Recreate a pseudopot block */

  smisc_size=0;
  if (ps_pots){
    i=strlen(ps_pots);
    if (((i/(*nsp))*(*nsp))!=i)
      fprintf(stderr,"Unexpected length to pseudopot block -- ignoring\n");
    ps_pot_len=(i/(*nsp));
    m->species_misc=realloc(m->species_misc,smisc_size+
			    strlen("%block SPECIES POT\n")+1);
    if (!m->species_misc) error_exit("Realloc error in a species_ block");
    strcpy(m->species_misc+smisc_size,"%block SPECIES POT\n");
    smisc_size+=strlen("%block SPECIES_POT\n");
    for(j=0;j<(*nsp);j++){
      cptr=ps_pots+j*ps_pot_len;
      cptr2=cptr;
      while ((*cptr2!=' ')&&(cptr2-cptr)<ps_pot_len) cptr2++;
      m->species_misc=realloc(m->species_misc,smisc_size+ion_sym_len+1+
			      (cptr2-cptr)+3);
      if (!m->species_misc) error_exit("Realloc error in a species_ block");
      strncpy(m->species_misc+smisc_size,ion_sym+ion_sym_len*j,ion_sym_len);
      smisc_size+=ion_sym_len;
      strcpy(m->species_misc+smisc_size," ");
      smisc_size+=1;
      strncpy(m->species_misc+smisc_size,cptr,(cptr2-cptr)+1);
      smisc_size+=(cptr2-cptr)+1;
      strcpy(m->species_misc+smisc_size,"\n");
      smisc_size+=1;
    }      
    m->species_misc=realloc(m->species_misc,smisc_size+
			    strlen("%endblock SPECIES_POT\n")+1);
    if (!m->species_misc) error_exit("Realloc error in a species_ block");
    strcpy(m->species_misc+smisc_size,"%endblock SPECIES_POT\n");
    smisc_size+=strlen("%endblock SPECIES_POT\n");
  }
  
  if (ion_force) m->forces=1;
  if (ion_force) free(ion_force);
  if (ion_pos) free(ion_pos);
  if (ion_sym) free(ion_sym);
  if (species_charge) free (species_charge);
  if (mag_mom) free(mag_mom);
  if (nionsp) free(nionsp);
  if (nsp) free(nsp);
  if (nsp_max) free(nsp_max);
  if (nkpts) free(nkpts);
  if (nsymops) free(nsymops);
  if (ps_pots) free(ps_pots);
  
  addabs(m->atoms,m->n,c->basis);

  if(endian) reverse8n(s->tol,1);
/* Castep works in Bohr, we work in Angstroms */
  (*(s->tol))*=BOHR;

  if ((debug>1)&&nsymops){  /* Note nsymops is not a valid pointer here */
    fprintf(stderr,"%d symmetry operations\n",s->n);
    if (sym_mat) fprintf(stderr,"symmetry matrices read\n");
  }

  
  /* Worry about dipole corrections */

  /* .check dip_corr_dir is always upper case, and anything on the
     command line is always lower-case */

  if ((elect->dip_corr_dir)&&(*(elect->dip_corr_dir)<='Z')){
    if (*elect->dip_corr_dir=='A')  /* all */
      *elect->dip_corr_dir='m';
    else{
      i=*elect->dip_corr_dir-'X';
      if ((i<0)||(i>=3))
        fprintf(stderr,"Warning: unexpected value for dipole_dir\n");
      else{
        j=1;
        for(k=0;k<3;k++){
          if (i==k) continue;
          if (!aeq(c->basis[i][k],0)){
            fprintf(stderr,"Warning: dipole_dir=%c, "
                    "but %c does not lie along %c\n",
                    *elect->dip_corr_dir,'a'+i,*elect->dip_corr_dir);
            fprintf(stderr,"  unclear what this means in Castep\n");
            j=0;
          }
        }
        if (j) *elect->dip_corr_dir='a'+i;
      }
    }
  }
  
  if(!match(head,"END")){
    fprintf(stderr,"Warning: unexpected end to check file. "
                   "Continuing regardless.\n");
  }
  return;

}

int match(char *s1, char *s2){
/* Returns one if s2 is the initial substring of s1 */
/* Why not use "#define match(s1,s2) (!strncmp(s1,s2,strlen(s2)))"? */
  int i;

  i=0;
  while(s1[i] && s2[i] && s1[i]==s2[i]) i++;

  if(!s2[i]) return(1);
  else return(0);
}

void reverse4(void *data){
/* reverse endian a single 4 byte int */
   int out;
   char *p1,*p2;

   p1=(char*)data;
   p2=(char*)&out;

   p2=p2+3;

   *(p2--)=*(p1++);
   *(p2--)=*(p1++);
   *(p2--)=*(p1++);
   *(p2--)=*(p1++);

   *((int*)data)=out;
}

void reverse8(void *data){
/* reverse endian a single 8 byte int */
   long out;
   char *p1,*p2;

   p1=(char*)data;
   p2=(char*)&out;

   p2=p2+7;

   *(p2--)=*(p1++);
   *(p2--)=*(p1++);
   *(p2--)=*(p1++);
   *(p2--)=*(p1++);
   *(p2--)=*(p1++);
   *(p2--)=*(p1++);
   *(p2--)=*(p1++);
   *(p2--)=*(p1++);

   *((long*)data)=out;
}

void reverse4n(int *data, int n){
/* reverse endian n words of 4 byte data */
   int i,out;
   char *p1,*p2;

   for(i=0;i<n;i++){
     p1=(char*)(data+i);
     p2=(char*)&out;
  
     p2=p2+3;
  
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
   
     *(data+i)=out;
   }
}

void reverse8n(double *data, int n){
/* reverse endian n words of 8 byte data */
   int i;
   double out;
   char *p1,*p2;

   for(i=0;i<n;i++){
     p1=(char*)(data+i);
     p2=(char*)&out;
  
     p2=p2+7;
  
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
     *(p2--)=*(p1++);
   
     *(data+i)=out;
   }
}

static void wave_read(FILE *infile, int nkpts, double *kpts_pos,
                      double* wkpt, int gamma, struct grid *g,
                      struct es *elect, int *i_grid, struct contents *m){
  int tmp,ns,b,i,k,nplwv,fft[3],ffft[3],fft_pts,dummy[5],offset,n0,n1,n2;
  int nsp,nspr;
  int *pwgrid=NULL;
  long fpos;
  double *dptr1,*dptr2,*dptr3,sum,phase_r,phase_i,phase_r2,phase_i2,dtmp;
  double max,min,conv,w;
  double kpoint[3],phi,*occ;
  int ii,jj,kk,ind,k2;
  int nfft[3],nfft_pts;
  char cbuff[CBUFF+1];
  
  char *band_range,*kpt_range,*spin_range;

  band_range=elect->band_range;
  kpt_range=elect->kpt_range;
  spin_range=elect->spin_range;
  
  occ=NULL;
  dptr2=NULL;
  nfft_pts=0;
  if (i_grid){
    for(i=0;i<3;i++){
      nfft[i]=i_grid[i];
    }
  }

/* Assume we have just read the key "wave" or "Gpt"
 *
 * we are therefore faced with:
 * ngx,ngy,ngz
 * max_plane_waves,nspinors,nbands,nkpts,nspins
 * do over spins
 *   do over kpoints
 *     kpoint[3],num_plane_waves
 *     do i=1,3
 *       plane_wave_to_grid_point_component_i_mapping_array[num_plane_waves]
 *     do over bands
 *       do over spinors
 *         wavefunction(complex)[num_plane_waves]
 * do over kpoints
 *   kpoint[3]
 *   do over spins
 *     band_occupancies[nbands]
 *     band_evalues[nbands]
 *
 * I am not convinced that the two "do over kpoints" loops will always see
 * the kpoints in the same order. As they explicitly record the kpoint
 * co-ordinates, there is no reason for confusion...
 *
 * If Gpt has been read (gamma==1), only half the wavefunction is stored,
 *  and the rest needs reconstructing as psi(-k)=conjg(psi(k))
 */

  fread(&tmp,4,1,infile);
  if (endian) reverse4(&tmp);

  if (tmp!=12) error_exit("Error parsing wavefunction");

  fread(fft,12,1,infile);
  fseek(infile,4,SEEK_CUR);
  if (endian) {
     reverse4(fft);
     reverse4(fft+1);
     reverse4(fft+2);
  }
  fft_pts=fft[0]*fft[1]*fft[2];
  if (i_grid==NULL){
    nfft_pts=fft_pts;
    for(i=0;i<3;i++) nfft[i]=fft[i];
  }
  else{
    for(i=0;i<3;i++) nfft[i]=max(nfft[i],fft[i]);
    nfft_pts=nfft[0]*nfft[1]*nfft[2];
  }
  if (debug) fprintf(stderr,"Wavefunction grid     %d %d %d\n",
                 fft[0],fft[1],fft[2]);

  if (debug&&gamma) fprintf(stderr,"Wavefunction stored in gamma point form\n");

  fread(&tmp,4,1,infile);
  if (endian) reverse4(&tmp);

  if (tmp==16){
     fread(dummy,16,1,infile);
     if (endian) reverse4n(dummy,4);
  } else if (tmp==20){
     fread(dummy,20,1,infile);
     if (endian) reverse4n(dummy,5);
     if (dummy[1]!=elect->nspinors){
       fprintf(stderr,"Confused about number of spinors. Is it %d or %d?\n",
		  dummy[1],elect->nspinors);
       error_exit("Aborting");
     }
     dummy[1]=dummy[2];
     dummy[2]=dummy[3];
     dummy[3]=dummy[4];
     if (debug>2){
       if (elect->nspinors==1)
	 fprintf(stderr,"Wavefunction stored with single spinor\n");
       else
	 fprintf(stderr,"Wavefunction stored with two spinors\n");
     }
  }
  else error_exit("Error parsing wavefunction");

  if (debug>2) fprintf(stderr,"Read: %d %d %d %d\n",dummy[0],
                        dummy[1],dummy[2],dummy[3]);
  if (nbands!=dummy[1]){
    if (debug>2) fprintf(stderr,"Resetting nbands from %d to %d."
          " Reading a .orbitals file?\n",nbands,dummy[1]);
    nbands=dummy[1];
  }
  if (dummy[2]!=nkpts){
    fprintf(stderr,"Confused about number of k-points. Is it %d or %d?\n",
            dummy[2],nkpts);
    error_exit("Aborting");
  }
  if (dummy[3]!=elect->nspins){
    fprintf(stderr,"Confused about number of spins. Is it %d or %d?\n",
            dummy[3],elect->nspins);
    error_exit("Aborting");
  }

  fseek(infile,4,SEEK_CUR); 

  /* If we want weighted sums, we now need to skip forwards, read in
     weights, and then skip back to this point. Tedious. */

  if ((flags&OCC_WEIGHT)||(flags&K_WEIGHT)){
    fpos=ftell(infile);
    for(ns=0;ns<elect->nspins;ns++){
      for(k=1;k<=nkpts;k++){
	fread(&tmp,4,1,infile);
	if (endian) reverse4(&tmp);
	if (tmp!=28) error_exit("Error seeking in wavefunction");
       	fseek(infile,24,SEEK_CUR);
	fread(&nplwv,4,1,infile);
	if (endian) reverse4(&nplwv);
	fseek(infile,4,SEEK_CUR);
	/* Skip three plane wave to grid component arrays */
	fseek(infile,3*(8+4*(long)nplwv),SEEK_CUR);
	/* Skip nbands*nspinors complex wavefunction */
	fseek(infile,nbands*elect->nspinors*(8+16*(long)nplwv),SEEK_CUR);
      }
    }
    occ=malloc(8*nkpts*elect->nspins*nbands);
    if(!occ) error_exit("Malloc error for occs\n");
    occ_read(infile,occ,kpts_pos,nkpts,elect->nspins,nbands);

    /* Rest file pointer */
    fseek(infile,fpos,SEEK_SET);
  }
  
  for(ns=0;ns<elect->nspins;ns++){
    for(k2=1;k2<=nkpts;k2++){
/* record of kpoint[3],nplwv */
      fseek(infile,4,SEEK_CUR);
      fread(kpoint,3*8,1,infile);
      fread(&nplwv,4,1,infile);
      if (endian) reverse8n(kpoint,3);
      if (endian) reverse4(&nplwv);
      k=1+find_kpt(kpoint,kpts_pos,nkpts);
      if(debug>2) fprintf(stderr,"kpoint no %d (%f,%f,%f) nplwv=%d\n",
                          k,kpoint[0],kpoint[1],kpoint[2],nplwv);
      fseek(infile,4,SEEK_CUR);
/* Read component to FFT grid mapping */
      fread(&tmp,4,1,infile);
      if (endian) reverse4(&tmp);
      if (4*nplwv!=tmp) error_exit("Error parsing band");
      if (debug>2) fprintf(stderr,"Spin=%d kpt=%d nplwv=%d\n",ns,k,nplwv);
      if (inrange(k,kpt_range)&&
	  ((elect->nspinors==2)||(inrange(ns,spin_range)))){
        if (!(pwgrid=malloc(12*nplwv))) error_exit("Malloc error for pwgrid");
        fread(pwgrid,tmp,1,infile);
        fseek(infile,8,SEEK_CUR);
        fread(pwgrid+nplwv,tmp,1,infile);
        fseek(infile,8,SEEK_CUR);
        fread(pwgrid+2*nplwv,tmp,1,infile);
        fseek(infile,4,SEEK_CUR);
        if(endian) reverse4n(pwgrid,3*nplwv);
        if (debug>2) fprintf(stderr,"Read pwgrid\n");
/* CASTEP stores these reciprocal space coeffs as -n/2 to n/2, whereas
   we want 0 to n/2, -n/2 to -1, so .. */

        for(i=0;i<nplwv;i++){
          if(pwgrid[i]<0) pwgrid[i]+=fft[0];
          if(pwgrid[nplwv+i]<0) pwgrid[nplwv+i]+=fft[1];
          if(pwgrid[2*nplwv+i]<0) pwgrid[2*nplwv+i]+=fft[2];
        }
      }
      else fseek(infile,tmp*3+5*4,SEEK_CUR);
      for(b=1;b<=nbands;b++){
        if (debug>2) fprintf(stderr,"Start band %d\n",b);
        for(nspr=0;nspr<elect->nspinors;nspr++){
          if (elect->nspinors==2) nsp=nspr;
	  else nsp=ns;
	  fread(&tmp,4,1,infile);
	  if (endian) reverse4(&tmp);
	  if (tmp!=16*nplwv) error_exit("Error parsing wavefunction band");
	  if (((flags&BANDS)||(flags&BANDPARITY))&&
              inrange(k,kpt_range)&&inrange(nsp,spin_range)&&
	      inrange(b,band_range)){
	    if (debug>2) fprintf(stderr,"starting band read,"
				 " kpt %d, band %d, spin %d\n",k,b,nsp);
	    if (!(dptr1=malloc(16*nplwv)))
	      error_exit("Malloc error in band read");
	    fread(dptr1,16*nplwv,1,infile);
	    fseek(infile,4,SEEK_CUR);
	    if(endian) reverse8n(dptr1,2*nplwv);
	    if (!(dptr2=malloc(16*fft_pts)))
              error_exit("Malloc error in band read");
	    for(i=0;i<2*fft_pts;i++) dptr2[i]=0.0;
	    sum=0;
	    for(i=0;i<nplwv;i++){
	      /*!!*/
	      offset=pwgrid[i+2*nplwv]+fft[2]*(pwgrid[i+nplwv]+
					       fft[1]*pwgrid[i]);
	      if ((offset<0)||(offset>fft_pts)){
		fprintf(stderr,"Impossible offset in wave_read off=%d i=%d\n",
			offset,i);
		exit(1);
	      }
	      sum+=dptr1[2*i]*dptr1[2*i]+dptr1[2*i+1]*dptr1[2*i+1];
	      dptr2[2*offset]=dptr1[2*i];
	      dptr2[2*offset+1]=dptr1[2*i+1];
	      /* For the gamma point case */
	      if (gamma&&(offset!=0)){
		n0=fft[2]-pwgrid[i+2*nplwv];
		if (n0==fft[2]) n0=0;
		n1=fft[1]-pwgrid[i+nplwv];
		if (n1==fft[1]) n1=0;
		n2=fft[0]-pwgrid[i];
		if (n2==fft[0]) n2=0;
		offset=n0+fft[2]*(n1+fft[1]*n2);
		if ((offset<0)||(offset>fft_pts)){
		  fprintf(stderr,
			  "Impossible -offset in wave_read off=%d i=%d\n",
			  offset,i);
		  exit(1);
		}
		sum+=dptr1[2*i]*dptr1[2*i]+dptr1[2*i+1]*dptr1[2*i+1];
		dptr2[2*offset]=dptr1[2*i];
		dptr2[2*offset+1]=-dptr1[2*i+1];
	      }
	    }
	    free(dptr1);
	    if (debug>2) fprintf(stderr,"Before FFT g=0 component is %g+%gi\n",
				 dptr2[0],dptr2[1]);
	    if (debug>2) fprintf(stderr,"And normalisation is %g\n",sum);

            if (((kpoint[0]==0)||aeq(fabs(kpoint[0]),0.5))&&
                ((kpoint[1]==0)||aeq(fabs(kpoint[1]),0.5))&&
                ((kpoint[2]==0)||aeq(fabs(kpoint[2]),0.5))&&
                (flags&BANDPARITY)) inv_parity(dptr2,fft,b,kpoint);

            /* Was the parity all we were requested to report? */
            if (!(flags&BANDS)){
              free(dptr2);
              dptr2=NULL;
              continue;
            }

	    /* Padding */
	    
	    if (i_grid){
	      if(debug>1)
		fprintf(stderr,"Padding wavefunction onto %dx%dx%d grid\n",
			nfft[0],nfft[1],nfft[2]);
	      if ((fft[0]==nfft[0])&&(fft[1]==nfft[1])&&(fft[2]==nfft[2])){
		if (debug>1)
		  fprintf(stderr,"Skipping null padding operation\n");
	      }
	      else{
		pad_recip(dptr2,fft,&dptr3,nfft);
		nfft_pts=nfft[0]*nfft[1]*nfft[2];
		free(dptr2);
		dptr2=dptr3;
	      }
	    }
	    
	    /* A FORTRAN data order ... */
	    ffft[0]=nfft[2];
	    ffft[1]=nfft[1];
	    ffft[2]=nfft[0];
	    fft3d(dptr2,ffft,1);
	    
	    if (!(dptr3=malloc(8*nfft_pts)))
	      error_exit("Malloc error in band read");
	    
	    if (((flags&BANDDEN)==BANDS)&&
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
		    dtmp=dptr2[ind];
		    dptr2[ind]=phase_r*dptr2[ind]-phase_i*dptr2[ind+1];
		    dptr2[ind+1]=phase_r*dptr2[ind+1]+phase_i*dtmp;
		  }
		}
	      }
	    }
	    phase_r=phase_i=phase_r2=phase_i2=0;
	    for(i=0;i<nfft_pts;i++){
	      if (dptr2[2*i]>0){
		phase_r+=dptr2[2*i];
		phase_i-=dptr2[2*i+1];
	      }else{
		phase_r2-=dptr2[2*i];
		phase_i2+=dptr2[2*i+1];
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
	    for (i=0;i<nfft_pts;i++){
	      if (flags&BANDPHASE){
		dptr3[i]=atan2(dptr2[2*i+1],dptr2[2*i]);
	      }
	      else if (flags&BANDREAL){
		dptr3[i]=dptr2[2*i];
	      }
	      else if (flags&BANDIMAG){
		dptr3[i]=dptr2[2*i+1];
	      }
	      else
		if ((flags&BANDDEN)==BANDDEN)
		  dptr3[i]=dptr2[2*i]*dptr2[2*i]+dptr2[2*i+1]*dptr2[2*i+1];
		else{
		  dptr3[i]=dptr2[2*i]*phase_r-dptr2[2*i+1]*phase_i;
		  dtmp=dptr2[2*i]*phase_i+dptr2[2*i+1]*phase_r;
		  if((fabs(dtmp)>.05))ii++;
		}
	      sum+=dptr3[i];
	      if(dptr3[i]<min) min=dptr3[i];
	      if(dptr3[i]>max) max=dptr3[i];
	    }
	    if (debug>2) fprintf(stderr,"Min=%g Max=%g Sum=%g\n",min,max,sum);
	    if (ii>0) fprintf(stderr,"Warning: %d components with imaginary"
			      " part >0.05\n",ii);
	    free(dptr2);
	    /* Do we need to rescale? */
	    if (((flags&RAW)==0)&&((flags&BANDPHASE)==0)){ /* Yes */
	      if ((flags&BANDDEN)==BANDDEN) conv=1/(*cr_cell_vol);
	      else conv=1/sqrt(*cr_cell_vol);
	      if (debug>2) fprintf(stderr,"Scaling wavefun by %f\n",conv);
	      for(i=0;i<nfft_pts;i++) dptr3[i]*=conv;
	    }
	    /* Do we need to weight? */
	    if ((flags&OCC_WEIGHT)||(flags&K_WEIGHT)){
	      w=1;
	      /* Care: ns loop starts at 0, k and b loops at 1 */
	      if (flags&OCC_WEIGHT)
		w*=occ[elect->nspins*nbands*(k-1)+ns*nbands+(b-1)];
	      if (flags&K_WEIGHT)
		w*=wkpt[k-1];
	      /* If we want densities, and we do not have spins, each
	       * band is doubly occupied */
	      if ((elect->nspins==1)&&(elect->nspinors==1)&&(flags&BANDDEN))
		w*=2;
	      if ((w!=1)&&(!(flags&BANDDEN))) w=sqrt(w);
	      if (debug)
		fprintf(stderr,"Using weight %f for ns=%d k=%d band=%d\n",
			w,ns,k,b);
	      if (debug>1)
		fprintf(stderr,"  kpt weight %f occupancy %f\n",
			wkpt[k-1],
			occ[elect->nspins*nbands*(k-1)+ns*nbands+(b-1)]);
		
	      if (w!=1)
		for(i=0;i<nfft_pts;i++) dptr3[i]*=w;
	    }
	    if (!(flags&ACCUMULATE)){
	      g->data=dptr3;
	      for(i=0;i<3;i++) g->size[i]=nfft[i];
	      g->name=malloc(40);
	      if (elect->nspinors==2)
		sprintf(g->name,"band_vs%d_k%d_b%d",nspr,k,b);
	      else if (elect->nspins==2)
		sprintf(g->name,"band_s%d_k%d_b%d",ns,k,b);
	      else
		sprintf(g->name,"band_k%d_b%d",k,b);
	      g->next=malloc(sizeof(struct grid));
	      g=g->next;
	      g->data=NULL;
	      g->next=NULL;
              if ((flags&OCC_WEIGHT)||(flags&K_WEIGHT)){
                snprintf(cbuff,CBUFF,
                         "Weight %f used for spin %d kpt %d band %d",
                         w,(elect->nspinors==2)?nspr:ns,k,b);
                add_cmt(m->comment,cbuff);
              }
	    }else{  /* Are accumulating */
	      if (!g->data){  /* This is the first set */
		g->data=dptr3;
		for(i=0;i<3;i++) g->size[i]=nfft[i];
		g->name=malloc(40);
		sprintf(g->name,"bands"); /* Don't move to a new grid */
	      }else{
		for(i=0;i<nfft_pts;i++) g->data[i]+=dptr3[i];
		free(dptr3);
	      }
	    }
	  }  /* if (inrange) */
	  else fseek(infile,16*nplwv+4,SEEK_CUR);
	} /*for(nspr=...) */
      } /* for(b=...) */
      if (pwgrid) {free(pwgrid); pwgrid=NULL;}
    } /* for(k=...) */
  } /* for(ns=...) */
  if ((flags&ACCUMULATE)&&g->data){ /* Now move to a new grid */
    g->next=malloc(sizeof(struct grid));
    g=g->next;
    g->data=NULL;
    g->next=NULL;
  }    
  if (pwgrid) free(pwgrid);
  if (occ) free(occ); /* And we have printed occupancies */
  else{
    if (flags&OCCUPANCIES){
      dptr1=malloc(8*nkpts*elect->nspins*nbands);
      if(!dptr1) error_exit("Malloc error for occs\n");
      occ_read(infile,dptr1,kpts_pos,nkpts,elect->nspins,nbands);
      free(dptr1);
    }
  }
}

/* Note that the two spinor components share the same occupancy */
void occ_read(FILE *infile, double* occ, double *kpts_pos, int nkpts,
              int nspins, int nbands){
  int tmp,k,k2,ns,b,i;
  double *dptr1,*dptr2,*evals,total;
  double kpoint[3];

  if (!(evals=malloc(8*nbands*nspins*nkpts)))
    error_exit("Malloc error for evals\n");
  for (k2=0;k2<nkpts;k2++){
    fread(&tmp,4,1,infile);
    if (endian) reverse4(&tmp);
    if (tmp!=24) error_exit("Error parsing end of wavefunction");
    fread(kpoint,3*8,1,infile);
    if (endian) reverse8n(kpoint,3);
    fseek(infile,4,SEEK_CUR); 
    k=find_kpt(kpoint,kpts_pos,nkpts);
    for(ns=0;ns<nspins;ns++){
      dptr1=occ+nspins*nbands*k+ns*nbands;
      dptr2=evals+nspins*nbands*k+ns*nbands;
      fread(&tmp,4,1,infile);
      if (endian) reverse4(&tmp);
      if (tmp!=8*nbands) error_exit("Error parsing end of wavefunction");
      fread(dptr1,8*nbands,1,infile);
      if (endian) reverse8n(dptr1,nbands);
      fseek(infile,4,SEEK_CUR);
      fread(&tmp,4,1,infile);
      if (endian) reverse4(&tmp);
      if (tmp!=8*nbands) error_exit("Error parsing end of wavefunction");
      fread(dptr2,8*nbands,1,infile);
      if (endian) reverse8n(dptr2,nbands);
      fseek(infile,4,SEEK_CUR);
    }
  }
  if (flags&OCCUPANCIES){
    fprintf(stderr,"              kpoint              band spin "
	    " occupancy      evalue (eV)\n");
    i=0;
    total=0;
    for(k=0;k<nkpts;k++)
      for(ns=0;ns<nspins;ns++)
        for(b=0;b<nbands;b++){
	  fprintf(stderr,"( % 8f % 8f % 8f )  %3d  %d  %10f   %14f\n",
		  kpts_pos[3*k],kpts_pos[3*k+1],kpts_pos[3*k+2],b+1,ns,
                  occ[i],evals[i]*H_eV);
          total+=occ[i];
          i++;
        }
    fprintf(stderr,"                                  Total:  %11f\n",total);
  }
          
  free(evals);
}


int inrange(int x, char *range){
/* determine whether x is within a range given by comma-separated
 * whitespace-free list of ints and hyphen-separated ranges
 */
  char *cptr;
  int test1,test2;

  cptr=range;
  if (!*cptr) return(0);
  while(*cptr){
    test1=0;
    while ((*cptr>='0')&&(*cptr<='9'))
      test1=10*test1+(*(cptr++)-'0');
    if (*cptr=='-'){ /* a range */
      cptr++;
      if((*cptr>='0')&&(*cptr<='9')){
        test2=0;
        while ((*cptr>='0')&&(*cptr<='9'))
          test2=10*test2+(*(cptr++)-'0');
      } else test2=1<<30;
    }else test2=test1;

    if ((*cptr!=',')&&(*cptr!=0)) error_exit("Parse error in inrange");

    if ((x>=test1)&&(x<=test2)) return(1);

    if (*cptr==',') cptr++;
  }
  return(0);
}


static int find_kpt(double *kpt, double *kpts_pos, int nkpts){
  int i;
  int hit;
  double old_tol;
  
  hit=-1;
  /* Save old global tol, and reuse aeq() with a fixed tol */
  old_tol=tol;
  tol=1e-8;
  
  for(i=0;i<nkpts;i++){
    if((aeq(kpt[0],kpts_pos[3*i]))&&
       (aeq(kpt[1],kpts_pos[3*i+1]))&&
       (aeq(kpt[2],kpts_pos[3*i+2]))){
      if (hit==-1) hit=i;
      else{
        fprintf(stderr,"Confused: seem to have identical kpts at\n");
        fprintf(stderr,"%d (%f, %f, %f) and %d (%f, %f, %f)\n",
                hit,kpts_pos[3*hit],kpts_pos[3*hit+1],kpts_pos[3*hit+2],
                i,kpts_pos[3*i],kpts_pos[3*i+1],kpts_pos[3*i+2]);
      }
    }
  }

  if (hit==-1){
    fprintf(stderr,"Unexpected kpt (%f, %f, %f)\n",kpt[0],kpt[1],kpt[2]);
    exit(1);
  }

  tol=old_tol;
  return hit;
}
