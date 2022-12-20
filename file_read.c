/* Collect file reading selection into one place */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "c2xsf.h"


/* Abinit */

void abinit_charge_read(FILE* infile, struct unit_cell *c, struct contents *m,
                 struct kpts *kp, struct grid *gptr, struct es *elect);
void abinit_eig_read(FILE* infile, struct unit_cell *c, struct contents *m,
		     struct kpts *k, struct symmetry *s, struct es *e);
void abinit_in_read(FILE* infile, struct unit_cell *c, struct contents *m,
                 struct kpts *k, struct symmetry *s, struct es *e);
void abinit_psi_read(FILE* infile, struct unit_cell *c,
                     struct contents *m, struct kpts *kp, struct grid *gptr,
                     struct es *elect, int *i_grid);
/* Castep */

void bands_read(FILE* infile, struct unit_cell *c, struct contents *m,
		struct kpts *k, struct symmetry *s,struct es *e);
void cell_read(FILE* infile, struct unit_cell *c, struct contents *m,
               struct kpts *k, struct symmetry *s);
void chdiff_read(FILE* infile, struct grid *g);
void check_read(FILE* infile, struct unit_cell *c, struct contents *m,
                struct kpts *k, struct symmetry *s, struct grid *g,
                struct es *elect, int *i_grid);
void denfmt_read(FILE* infile, struct unit_cell *c, struct contents *m,
		 struct grid *gptr, struct es *elect, int rescale);
void esp_read(FILE* infile, struct contents *m, struct grid *g,
	      struct es *elect);
void geom_read(FILE *infile, struct unit_cell *c, struct contents *m,
	       struct time_series *ts);

/* Elk */

void elk_read(FILE* outfile, struct unit_cell *c, struct contents *m,
                  struct kpts *k, struct es *e);
void elk3d_read(FILE* infile, struct unit_cell *c, struct contents *m,
		struct kpts *k, struct grid *gptr, struct es *e);

/* QE */

void qe_rho_read(FILE* infile, struct unit_cell *c, struct contents *m,
                struct kpts *k, struct symmetry *s, struct grid *g,
                struct es *elect, int *i_grid);
void qe_xml_read(FILE* infile, char *filename, struct unit_cell *c,
		 struct contents *m,
                 struct kpts *k, struct symmetry *s, struct grid *g,
		 struct es *elect, struct time_series *ts, int *i_grid);

/* Siesta */

/* fdf_read() is prototyped in c2xsf.h */
void rho_read(FILE* infile, struct unit_cell *c,
              struct contents *m, struct grid *gptr, struct es *elect);

/* VASP */

void vasp_eigenval_read(FILE *infile, struct unit_cell *c, struct contents *m,
			struct kpts *k, struct es *e);
void vasp_psi_read(FILE* infile, char * filename, struct unit_cell *c,
		struct contents *m, struct kpts *k, struct grid *gptr,
		struct es *elect, int *i_grid);
void vasp_read(FILE* infile, char *filename,
		struct unit_cell *c, struct contents *m, struct kpts *k,
                struct grid *gptr, struct es *elect);

/* Generic */

void cif_read(FILE* infile, struct unit_cell *c, struct contents *m,
              struct symmetry *s);
void cube_read(FILE* infile, struct unit_cell *c, struct contents *m,
               struct grid *gptr);
void dx_read(FILE* infile, struct unit_cell *c, struct grid *g);
void npy_read(FILE* infile, struct grid *g);
void pdb_read(FILE* infile, struct unit_cell *c, struct contents *m);
void xsf_read(FILE* infile, struct unit_cell *c, struct contents *m,
	      struct grid *gptr);
void xyz_read(FILE* infile, struct unit_cell *c, struct contents *m);

/* Misc */

void crystal_read(FILE* infile, struct unit_cell *c, struct contents *m,
                  struct symmetry *s);
void fort34_read(FILE* infile, struct unit_cell *c, struct contents *m,
              struct symmetry *s);
void gcoeff_read(FILE *infile, struct unit_cell *c, struct contents *m,
                 struct kpts *k, struct grid *g, struct es *e,
		 int *i_grid);
void shelx_read(FILE* infile, struct unit_cell *c, struct contents *m);


void file_read(char *file1, FILE* infile, struct unit_cell *c,
	       struct contents *m, struct kpts *k, struct symmetry *s,
	       struct grid *g, struct es *elect, struct time_series *ts,
	       int *i_grid){
  int i,i_frame,*i_frame_ptr;
  char *cdfile,*cptr;

  i_frame_ptr=(int*)dict_get(m->dict,"i_frame");
  
  if (debug>1) fprintf(stderr,"Reading %s\n",file1);
  i=strlen(file1);

  if((i>4)&&(!strcmp(file1+i-4,".pdb")))
    pdb_read(infile,c,m);
  else if ((i==7)&&(!strcmp(file1,"fort.34")))
    fort34_read(infile,c,m,s);
  else if ((i==18)&&(!strcmp(file1,"charge-density.dat")))
    qe_rho_read(infile,c,m,k,s,g,elect,i_grid);
  else if ((i>=19)&&(!strcmp(file1+i-19,"/charge-density.dat")))
    qe_rho_read(infile,c,m,k,s,g,elect,i_grid);
  else if ((i==6)&&(!strcmp(file1+i-6,"elk.in")))
    elk_read(infile,c,m,k,elect);
  else if ((i>=7)&&(!strcmp(file1+i-7,"/elk.in")))
    elk_read(infile,c,m,k,elect);
  else if ((i==12)&&(!strcmp(file1+i-12,"GEOMETRY.OUT")))
    elk_read(infile,c,m,k,elect);
  else if ((i>=13)&&(!strcmp(file1+i-13,"/GEOMETRY.OUT")))
    elk_read(infile,c,m,k,elect);
  else if ((i>6)&&(!strcmp(file1+i-6,"3D.OUT")))
    elk3d_read(infile,c,m,k,g,elect);
  else if ((i>4)&&(!strcmp(file1+i-4,".xml")))
    qe_xml_read(infile,file1,c,m,k,s,g,elect,ts,i_grid);
  else if ((i>2)&&(!strcmp(file1+i-2,"12")))
    crystal_read(infile,c,m,s);
  else if ((i>4)&&(!strcmp(file1+i-4,".res")))
    shelx_read(infile,c,m);
  else if ((i>4)&&(!strcmp(file1+i-4,".cif")))
    cif_read(infile,c,m,s);
  else if ((i>6)&&(!strcmp(file1+i-6,".mmcif")))
    cif_read(infile,c,m,s);
  else if ((i>5)&&(!strcmp(file1+i-5,".pdbx")))
    cif_read(infile,c,m,s);
  else if ((i>4)&&(!strcmp(file1+i-4,".dat"))){
    dict_strcat(m->dict,"cell_is_onetep","");
    cell_read(infile,c,m,k,s);
  }
  else if ((i>4)&&(!strcmp(file1+i-4,".npy")))
    npy_read(infile,g);
  else if ((i>4)&&(!strcmp(file1+i-4,".cub")))
    cube_read(infile,c,m,g);
  else if ((i>5)&&(!strcmp(file1+i-5,".cube")))
    cube_read(infile,c,m,g);
  else if ((i>5)&&(!strcmp(file1+i-5,"_CUBE")))
    cube_read(infile,c,m,g);
  else if ((i>3)&&(!strcmp(file1+i-3,".dx")))
    dx_read(infile,c,g);
  else if ((i>4)&&(!strcmp(file1+i-4,".xyz")))
    xyz_read(infile, c, m);
  else if ((i>4)&&(!strcmp(file1+i-4,".xsf")))
    xsf_read(infile,c,m,g);
  else if ((i>=7)&&(!strcmp(file1+i-7,"CONTCAR")))
    vasp_read(infile,file1,c,m,k,g,elect);
  else if ((i>=6)&&(!strcmp(file1+i-6,"LOCPOT")))
    vasp_read(infile,file1,c,m,k,g,elect);
  else if ((i>=6)&&(!strcmp(file1+i-6,"POSCAR")))
    vasp_read(infile,file1,c,m,k,g,elect);
  else if ((i>=6)&&(!strcmp(file1+i-6,"CHGCAR")))
    vasp_read(infile,file1,c,m,k,g,elect);
  else if ((i>=3)&&(!strcmp(file1+i-3,"CHG")))
    vasp_read(infile,file1,c,m,k,g,elect);
  else if ((i>=7)&&(!strcmp(file1+i-7,"WAVECAR")))
    vasp_psi_read(infile,file1,c,m,k,g,elect,i_grid);
  else if ((i>=6)&&(!strcasecmp(file1+i-6,"GCOEFF")))
    gcoeff_read(infile,c,m,k,g,elect,i_grid);
  else if ((i>=10)&&(!strcmp(file1+i-10,"GCOEFF.txt")))
    gcoeff_read(infile,c,m,k,g,elect,i_grid);
  else if ((i>=8)&&(!strcmp(file1+i-8,"EIGENVAL")))
    vasp_eigenval_read(infile, c, m, k, elect);
  else if ((i>=3)&&(!strcmp(file1+i-3,"DEN")))
    abinit_charge_read(infile,c,m,k,g,elect);
  else if ((i>=3)&&(!strcmp(file1+i-3,"POT")))
    abinit_charge_read(infile,c,m,k,g,elect);
  else if ((i>=5)&&(!strcmp(file1+i-5,"VCLMB")))
    abinit_charge_read(infile,c,m,k,g,elect);
  else if ((i>=3)&&(!strcmp(file1+i-3,"WFK")))
    abinit_psi_read(infile,c,m,k,g,elect,i_grid);
  else if ((i>=3)&&(!strcmp(file1+i-3,"DDB")))
    abinit_in_read(infile,c,m,k,s,elect);
  else if ((i>=3)&&(!strcmp(file1+i-3,".in")))
    abinit_in_read(infile,c,m,k,s,elect);
  else if ((i>=4)&&(!strcmp(file1+i-4,".abi")))
    abinit_in_read(infile,c,m,k,s,elect);
  else if ((i>=3)&&(!strcmp(file1+i-3,"EIG")))
    abinit_eig_read(infile,c,m,k,s,elect);
  else if ((i>=7)&&(!strcmp(file1+i-7,"den_fmt")))
    denfmt_read(infile,c,m,g,elect,1);
  else if ((i>=7)&&(!strcmp(file1+i-7,"pot_fmt")))
    denfmt_read(infile,c,m,g,elect,2);
  else if ((i>=7)&&(!strcmp(file1+i-7,"elf_fmt")))
    denfmt_read(infile,c,m,g,elect,0);
  else if ((i>=7)&&(!strcmp(file1+i-7,"cst_esp")))
    esp_read(infile,m,g,elect);
  else if ((i>=4)&&(!strcmp(file1+i-4,".elf"))){
    flags|=RAW;
    dict_add(m->dict,"grid_name","ELF");
    esp_read(infile,m,g,elect);
  }
  else if ((i>=5)&&(!strcmp(file1+i-5,"bands")))
    bands_read(infile,c,m,k,s,elect);
  else if ((i>=5)&&(!strcmp(file1+i-5,".geom")))
    geom_read(infile,c,m,ts);
  else if ((i>=4)&&(!strcmp(file1+i-4,".fdf")))
    fdf_read(infile,c,m,k,elect);
  else if ((i>=3)&&(!strcasecmp(file1+i-3,".xv")))
    xv_read(infile,c,m);
  else if ((i>=4)&&(!strcasecmp(file1+i-4,".rho")))
    rho_read(infile,c,m,g,elect);
  else if ((i>=6)&&(!strcasecmp(file1+i-6,".rhoxc")))
    rho_read(infile,c,m,g,elect);
  else if ((i>=5)&&(!strcasecmp(file1+i-5,".toch")))
    rho_read(infile,c,m,g,elect);
  else if ((i>=3)&&(!strcasecmp(file1+i-3,".vh")))
    rho_read(infile,c,m,g,elect);
  else if ((i>=3)&&(!strcasecmp(file1+i-3,".vt")))
    rho_read(infile,c,m,g,elect);
  else if ((i>=3)&&(!strcasecmp(file1+i-3,".kp")))
    siesta_kp_read(infile,k);
  /* If not a Castep ending, check for a VASP beginning */
  else if ((!(((i>5)&&(!strcmp(file1+i-5,".cell")))||
              ((i>6)&&(!strcmp(file1+i-6,".check")))||
              ((i>4)&&(!strcmp(file1+i-4,".dat")))||
              ((i>11)&&(!strcmp(file1+i-11,".castep_bin")))||
              ((i>7)&&(!strcmp(file1+i-7,".chdiff")))))&&
           ((!strncmp(file1,"CONTCAR",7))||(!strncmp(file1,"POSCAR",6))||
            (!strncmp(file1,"CHGCAR",6))||(!strncmp(file1,"LOCPOT",6))||
	    (!strncmp(file1,"CHG",3))))
    vasp_read(infile,file1,c,m,k,g,elect);
  else{
    i=fgetc(infile);
    rewind(infile);
    if ((i==0)||(i==30)||(i==10)) /* possible first byte of .check or */
                                  /* .castep_bin */
      check_read(infile,c,m,k,s,g,elect,i_grid);
    else
      cell_read(infile,c,m,k,s);
  }
  fclose(infile);

  if ((i_frame_ptr)&&(ts->nsteps)){
    i_frame=*i_frame_ptr;
    if (i_frame<0) {
      i_frame=ts->nsteps+i_frame;
      if (debug) fprintf(stderr,"Setting frame to %d\n",i_frame);
    }
    if (i_frame<0)
      error_exit("Frame number cannot be negative!");
    if (i_frame>=ts->nsteps){
      fprintf(stderr,"Error: timeseries steps read numbered 0 to %d,"
	      " but frame %d requested.\n",ts->nsteps-1,i_frame);
      exit(1);
    }
    if (ts->nc>i_frame){
      if (c->basis) free(c->basis);
      *c=ts->cells[i_frame];
    }
    if (ts->nm>i_frame){
      if (m->atoms) free (m->atoms);
      m->atoms=ts->m[i_frame].atoms;
      m->n=ts->m[i_frame].n;
      m->forces=ts->m[i_frame].forces;
    }
    if (ts->nen>i_frame){
      if (elect->energy) free (elect->energy);
      elect->energy=ts->energies+i_frame;
    }
  }
  
  if (flags&CHDIFF) { /* We want to read a .chdiff file too */
    cdfile=malloc(strlen(file1)+8);
    if (!cdfile) error_exit("Malloc error for chdiff filename");
    strcpy(cdfile,file1);
    cptr=cdfile+strlen(cdfile);
    while((cptr>cdfile)&&(*cptr!='.')) cptr--;
    if (*cptr!='.') error_exit("Error constructing chdiff filename");
    strcpy(cptr+1,"chdiff");
    infile=fopen(cdfile,"rb");
    if(!infile){
      fprintf(stderr,"Error, unable to open %s for reading.\n",cdfile);
      exit(1);
    }
    chdiff_read(infile,g);
    fclose(infile);
    free(cdfile);
  }

  if (flags&CST_ESP) { /* We want to read a .cst_esp file too */
    cdfile=malloc(strlen(file1)+9);
    if (!cdfile) error_exit("Malloc error for cst_esp filename");
    strcpy(cdfile,file1);
    cptr=cdfile+strlen(cdfile);
    while((cptr>cdfile)&&(*cptr!='.')) cptr--;
    if (*cptr!='.') error_exit("Error constructing cst_esp filename");
    strcpy(cptr+1,"cst_esp");
    infile=fopen(cdfile,"rb");
    if(!infile){
      fprintf(stderr,"Error, unable to open %s for reading.\n",cdfile);
      exit(1);
    }
    esp_read(infile,m,g,elect);
    fclose(infile);
    free(cdfile);
  }
}
