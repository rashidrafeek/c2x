/* Convert CASTEP .check file to XCrySDen .xsf file
 * Extract other data from a .check file 
 * Form supercells, perform symmetry analysis, and
 * Convert between various crystallographic file formats.
 */


/* Copyright (c) 2007 -- 2022 MJ Rutter 
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
#include <unistd.h>  /* isatty */
#include <math.h> /* fabs */
#include <string.h>

#ifdef SPGLIB
#include <spglib.h>
#endif

#define C2X_MAIN
#include "c2xsf.h"
#undef C2X_MAIN

#define C2X_MU 4.2

void help(void);
void refs(void);
void formats(void);

void xsf_write(FILE* outfile, struct unit_cell *c, struct contents *m,
               int molecule, struct grid *g);
void cube_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                struct grid *g);
void denfmt_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                struct grid *g);
void xplor_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                 struct grid *g);
void xplor_fudge(struct unit_cell *c, struct contents *m, struct grid *g);
void pdb_write(FILE* outfile, struct unit_cell *c, struct contents *m);
void cif_write(FILE* outfile, struct unit_cell *c, struct contents *m, 
               struct symmetry *s, int mm);
void cell_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                struct kpts *k, struct symmetry *s);
void cell_write_abc(FILE* outfile, struct unit_cell *c, struct contents *m,
                    struct kpts *k, struct symmetry *s);
void cell_write_abs(FILE* outfile, struct unit_cell *c, struct contents *m,
                    struct kpts *k, struct symmetry *s);
void cell_write_abc_abs(FILE* outfile, struct unit_cell *c, struct contents *m,
                        struct kpts *k, struct symmetry *s);
void dx_write(FILE* outfile, struct unit_cell *c, struct grid *g);
void vasp_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                struct grid *g);
void xyz_write(FILE* outfile, struct unit_cell *c, struct contents *m,
		struct time_series *ts);
void cml_write(FILE* outfile, struct unit_cell *c, struct contents *m);
void fdf_write(FILE* fdf, char* filename, struct unit_cell *c,
               struct contents *m,struct kpts *kp, struct grid *g,
	       struct es *e);
void fdf_write_bp(FILE *outfile, struct kpts *k, struct contents *m);
void xv_write(FILE* infile, char *filename, struct unit_cell *c,
              struct contents *m, struct grid *g);
void ccp4_write(FILE* outfile,struct unit_cell *c, struct grid *g);
void line_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                struct grid *g, char *line_spec, int lflags);
void lscan(char **p, struct contents *m, double x[3]);
void shelx_write(FILE* outfile, struct unit_cell *c, struct contents *m);
void py_write(FILE* outfile, struct unit_cell *c, struct contents *m,
              char type);
void fbin_write(FILE* outfile, struct grid *g);
void f15_write(struct unit_cell *c, struct contents *m, struct kpts *k);
void abinit_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                  struct kpts *k, struct symmetry *s, struct es *e);
void qe_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                  struct kpts *k, struct es *e);
void elk_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                  struct kpts *k, struct es *e);
void npy_write(FILE* outfile, struct grid *g);

void file_read(char *file1, FILE* infile, struct unit_cell *c,
	       struct contents *m, struct kpts *k, struct symmetry *s,
	       struct grid *g, struct es *elect, struct time_series *ts,
	       int *i_grid);
void cell_read(FILE* infile, struct unit_cell *c, struct contents *m,
               struct kpts *k, struct symmetry *s);

void molecule_fix(int m_abc[3], double m_rel[3],struct unit_cell *c,
		  struct contents *m, struct grid *g);
void rotation(struct unit_cell *c, struct contents *m, double new_basis[3][3]);
void primitive(struct unit_cell *c, struct contents *m, double basis[3][3]);
void shorten(double basis[3][3]);
double interpolate0d(struct grid *gptr,double x[3]);
void fstar(struct kpts *ks, struct unit_cell *c, struct symmetry *rs);
void dipole(struct unit_cell *c, struct contents *m,
            struct grid *g, struct es *elect);
void charge_corr(struct unit_cell *c, struct contents *m,
                 struct grid *g, struct es *elect);
void es_pot(struct unit_cell *c, struct contents *m,
            struct grid *g, struct es *elect, double musq);

void sym_kpts(struct kpts *k_in, struct kpts *k_out, struct symmetry *s,
              double basis[3][3]);
void vacuum_adjust(struct unit_cell *c, struct contents *m, double new_abc[3]);
void bands_write(FILE* outfile, struct unit_cell *c,
                 struct kpts *k, struct es *e);
void geom_write(FILE *outfile, struct unit_cell *c, struct contents *m,
		struct es *elect, struct time_series *ts);
void bxsf_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                struct es *elect, struct kpts *kp,
                struct symmetry *rs);
void tube(struct unit_cell *c, struct contents *m, int rpt[3], double spacing);
void print_gap(struct es *elect, struct kpts *kpt, struct contents *m);
double lda(double dens);
void data_combine(int op, struct unit_cell *c,
		  struct contents *m, struct kpts *k, struct symmetry *s,
		  struct grid *g, struct es *elect, struct time_series *ts,
		  int *i_grid,struct unit_cell *c2,
		  struct contents *m2, struct kpts *k2, struct symmetry *s2,
		  struct grid *g2, struct es *elect2, struct time_series *ts2);

void check_read(FILE* infile, struct unit_cell *c, struct contents *m,
                struct kpts *k, struct symmetry *s, struct grid *g,
                struct es *elect, int *i_grid);

/* Global variables for system description */

int debug,flags;
double tol=1e-4;

void version(){
  printf("c2x version " C2XSF_VER "\n");
#ifdef SPGLIB
  printf("Linked with spglib version %d.%d.%d\n",
         spg_get_major_version(),
         spg_get_minor_version(),
         spg_get_micro_version());
#else
  printf("Not linked with spglib\n");
#endif

  if (debug){
    printf("\nInternal conversion factors:\n");
    printf("   One Hartree = %.11f eV\n",H_eV);
    printf("   One Rydberg = %.11f eV\n",0.5*H_eV);
    printf("      One Bohr = %.14f A\n",BOHR);
    printf("       1/eps_0 = %.11f e^-1 V A\n",1/EPS0);
    printf(" h_bar/Hartree = %.14f fs\n",1000*H_ps);
  }

  printf("\nFurther documentation at https://www.c2x.org.uk/\n");
  exit(0);
}

int main(int argc, char **argv)
{
  int i,j,k,n,opt=1,expand,rotate,half_shift=0,no_mp=0,prim=0,compact=0;
  int molecule=0,reduce=0,vexpand=0,lflags=0,pr_occ=0,pr_gap=0,xc=0;
  int sort_style=0,calc_ef=0,firstfile_arg=0,op=0,replace_rho=0;
  int no_sym=0;
  int circ[3];
  int spg_op;
  int gen_mp=0, failure=0,calc_esp=0,sym_k=0,sym_list=0;
  int keep_velocities=0;
  int format,preserve_c,mask;
  int *i_grid,*i_frame;
  int i_expand[3];
  char *optp,*file1=NULL,*outname=NULL,*line_spec=NULL,*pt_spec=NULL;
  char *cdfile,*cptr;
  FILE *infile,*outfile;
  double abc[6],new_cell[3][3],new_cell_rel[3][3],*new_mp,zpt[3],g_cut;
  double new_abc[3],dtmp,scale,rescale;
  double tolmin=1e-4,ionic_charge,musq,tube_spacing,*m_rel,chg;
  int *m_abc;
  struct grid *gptr;
  struct unit_cell cell,cell2,nc;
  struct contents motif,motif2;
  struct symmetry sym,sym2;
  struct kpts kp,kp2;
  struct grid grid1,grid2;
  struct es elect,elect2;
  struct time_series series,series2;

/* Initialise all pointers to NULL, etc. */

  init_elect(&elect);
  init_elect(&elect2);
  
  init_motif(&motif);
  init_motif(&motif2);
  
  init_cell(&cell);
  cell2=cell;

  init_sym(&sym);
  sym2=sym;

  init_kpts(&kp);
  kp2=kp;
  
  m_abc=i_grid=i_frame=NULL;
  m_rel=NULL;
  new_mp=NULL;
  sym.gen=NULL;

  init_grid(&grid1);
  grid2=grid1;

  flags=0;

  init_tseries(&series);

  rescale=1;
  
  opt=1;
  optp=argv[opt];
  debug=molecule=format=expand=rotate=preserve_c=0;
  spg_op=0;

  for(i=0;i<3;i++) circ[i]=0;
  tube_spacing=0;
  
  for(i=0;i<3;i++) new_abc[i]=0;

  if ((strlen(argv[0])>=7)&&
      (!strcmp("cellsym",&argv[0][strlen(argv[0])-7]))) format=CELL;
  
  while ((opt<argc)&&(!file1)){
    switch(*optp){
    case 0:
      opt++;
      optp=argv[opt];
      break;
    case '-':
      if(*(optp+1)=='-'){ /* Gnu-style "--" option */
        if (!strcmp(optp,"--xsf")) format=XSF;
        else if (!strcmp(optp,"--bxsf")) {format=BXSF; flags|=OCCUPANCIES;}
        else if (!strcmp(optp,"--cube")) format=CUBE;
        else if (!strcmp(optp,"--mocube")) {format=CUBE; flags|=ALT_OUT;}
        else if (!strcmp(optp,"--xplor")) format=XPLOR;
        else if (!strcmp(optp,"--pdb")) format=PDB;
        else if (!strcmp(optp,"--pdbn")) {format=PDB; flags|=PDBN;}
        else if (!strcmp(optp,"--pdbx"))
	  {format=CIF; dict_strcat(motif.dict,"CIF_is_PDBx","");}
        else if (!strcmp(optp,"--cell")) format=CELL;
        else if (!strcmp(optp,"--cell_abc")) format=CELL_ABC;
        else if (!strcmp(optp,"--cell_abs")) format=CELL_ABS;
        else if (!strcmp(optp,"--cell_abc_abs")) format=CELL_ABC_ABS;
        else if (!strcmp(optp,"--one"))
          {format=CELL; flags|=ALT_OUT;}
        else if (!strcmp(optp,"--one_abc"))
          {format=CELL_ABC; flags|=ALT_OUT;}
        else if (!strcmp(optp,"--one_abs")) 
          {format=CELL_ABS; flags|=ALT_OUT;}
        else if (!strcmp(optp,"--one_abc_abs"))
          {format=CELL_ABC_ABS; flags|=ALT_OUT;}
        else if (!strcmp(optp,"--dx")) format=DX;
        else if (!strcmp(optp,"--vasp")) format=VASP;
        else if (!strcmp(optp,"--chgcar"))
          {format=VASP; flags|=CHGCAR; flags|=CHDEN;}
        else if (!strcmp(optp,"--xyz")) format=XYZ;
        else if (!strcmp(optp,"--xyze")) {format=XYZ; flags|=ALT_OUT;}
        else if (!strcmp(optp,"--cml")) format=CML;
        else if (!strcmp(optp,"--fdf")) format=FDF;
        else if (!strcmp(optp,"--fdf_bp")) {format=FDF_BP; sym_k=1; no_mp=0;}
        else if (!strcmp(optp,"--gnu")) lflags|=1;
        else if (!strcmp(optp,"--shelx")) format=SHELX;
        else if (!strcmp(optp,"--airss")) {format=SHELX; flags|=SHELX_AIRSS;}
        else if (!strcmp(optp,"--fbin")) format=FBIN;
        else if (!strcmp(optp,"--nibf")) {format=FBIN; flags|=LITTLE_OUT;}
        else if (!strcmp(optp,"--fort15")) format=FORT15;
        else if (!strcmp(optp,"--cif")) format=CIF;
        else if (!strcmp(optp,"--mmcif")) format=MMCIF;
        else if (!strcmp(optp,"--py")) format=PY;
        else if (!strcmp(optp,"--pya")) format=PYA;
        else if (!strcmp(optp,"--denfmt")) format=DENFMT;
        else if (!strcmp(optp,"--den_fmt")) format=DENFMT;
        else if (!strcmp(optp,"--abinit")) format=ABINIT;
	else if (!strcmp(optp,"--abinit8")) {format=ABINIT; flags|=ALT_OUT;}
        else if (!strcmp(optp,"--npy")) format=NPY;
        else if (!strcmp(optp,"--qe")) format=QE;
        else if (!strcmp(optp,"--qef")) {format=QE; flags|=FRAC;}
        else if (!strcmp(optp,"--bands")) {
          format=CASTEP_BANDS;
          flags|=OCCUPANCIES;
        }
        else if (!strcmp(optp,"--geom")) format=CASTEP_GEOM;
        else if (!strcmp(optp,"--xv")) format=XV;
        else if (!strcmp(optp,"--ccp4")) format=CCP4;
        else if (!strcmp(optp,"--gcoeff")) {format=CNULL; flags|=GCOEFF;}
        else if (!strcmp(optp,"--gcoeff_sorted")) {
          format=CNULL;
          flags|=GCOEFF;
          dict_strcat(motif.dict,"gcoeff_sort","");
        }
        else if (!strcmp(optp,"--wavecar")) {
          format=CNULL;
          flags|=GCOEFF;
          dict_strcat(motif.dict,"wavecar_output","");
        }
        else if (!strcmp(optp,"--wavecar45")) {
          format=CNULL;
          flags|=GCOEFF;
          flags|=ALT_OUT;
          dict_strcat(motif.dict,"wavecar_output","");
        }
        else if (!strcmp(optp,"--elk")) format=ELK;
	else if (!strcmp(optp,"--null")) format=CNULL;

        else if (!strcmp(optp,"--primitive")) spg_op|=CSPG_PRIM;
        else if (!strcmp(optp,"--primitive_nr")) spg_op|=CSPG_PRIM_NR;
        else if (!strcmp(optp,"--refine")) spg_op|=CSPG_REF;
        else if (!strcmp(optp,"--standardise")) spg_op|=CSPG_STD;
        else if (!strcmp(optp,"--standardize")) spg_op|=CSPG_STD;
        else if (!strcmp(optp,"--std")) spg_op|=CSPG_STD;
        else if (!strcmp(optp,"--std_ideal")) spg_op|=CSPG_STD_IDEAL;
        else if (!strcmp(optp,"--snap")) spg_op|=CSPG_SNAP;
        else if (!strcmp(optp,"--snap_tr")) spg_op|=CSPG_SNAP+CSPG_SNAP_TR;
        else if (!strcmp(optp,"--int")) spg_op|=CSPG_INT;
        else if (!strcmp(optp,"--schoen")) spg_op|=CSPG_SCH;
        else if (!strcmp(optp,"--symmetry")) spg_op|=CSPG_SYM;
        else if (!strcmp(optp,"--sym")) spg_op|=CSPG_SYM;
        else if (!strcmp(optp,"--point")) spg_op|=CSPG_PNT;
        else if (!strcmp(optp,"--list")) spg_op|=(CSPG_LST|CSPG_SYM);

        else if (!strcmp(optp,"--help")) help();
        else if (!strcmp(optp,"--refs")) refs();
        else if (!strcmp(optp,"--constants")) {debug=1; version();}
        else if (!strcmp(optp,"--units")) {debug=1; version();}
        else if (!strcmp(optp,"--formats")) formats();
        else if (!strcmp(optp,"--version")) version();
        else if (!strcmp(optp,"--calc")){
          single_scan(argv[opt+1],&dtmp,NULL);
          if (debug) printf("%s = ",argv[opt+1]);
          printf("%.15g\n",dtmp);
          exit(0);
        }
        else if (!strcmp(optp,"--vel")) keep_velocities=1;
        else if (!strcmp(optp,"--gap")) {
	  flags|=OCCUPANCIES;
	  pr_gap=1;
	}
        else if (!strcmp(optp,"--efermi")) {
          flags|=OCCUPANCIES;
	  calc_ef=1;
	}
	else if ((!strcmp(optp,"--sym_list"))||(!strcmp(optp,"--symlist"))) {
           sym_list=1;
 	}
        else if (!strncmp(optp,"--frame=",8)) {
	  i_frame=malloc(sizeof(int));
	  if (!i_frame) error_exit("Malloc error for int!");
	  if (sscanf(optp+8,"%d",i_frame)!=1)
	    error_exit("Failed to parse --frame= argument");
	}
	else if (!strcmp(optp,"--merge")) op|=C2X_MERGE;
	else if (!strcmp(optp,"--add")) op=(op&~C2X_EXC_MASK)+C2X_ADD;
	else if (!strcmp(optp,"--sum")) op=(op&~C2X_EXC_MASK)+C2X_ADD;
	else if (!strcmp(optp,"--diff")) op=(op&~C2X_EXC_MASK)+C2X_DIFF;
	else if (!strcmp(optp,"--sub")) op=(op&~C2X_EXC_MASK)+C2X_DIFF;
	else if (!strcmp(optp,"--mask")) op=(op&~C2X_EXC_MASK)+C2X_MASK;
	else if (!strcmp(optp,"--mult")) op=(op&~C2X_EXC_MASK)+C2X_MASK;
	else if (!strcmp(optp,"--replace_rho")) replace_rho=1;
	else {
          fprintf(stderr,"Invalid option %s.\n%s -h for usage.\n",
                   optp,argv[0]);
          exit(1);
        }
        opt++;
        optp=argv[opt]; 
        break;
      }
      while(*(++optp)){
        switch(*optp){
        case 'v':
          debug++;
          break;
        case 'h':
          help();
          break;
        case 'V':
          version();
          break;
        case 'a':
          rotate=1;
          break;
        case 'A':
          flags|=ACCUMULATE;
          break;
        case 'b':
          flags|=BANDS;
          if(*(optp+1)=='='){
            elect.band_range=optp+2;
            while(*((++optp)+1));
          }
          break;
        case 'B':
          flags|=BANDS+BANDDEN;
          if(*(optp+1)=='='){
            elect.band_range=optp+2;
            while(*((++optp)+1));
          }
          break;
        case 'c':
          flags|=CHDEN;
          break;
        case 'C':
          compact=1;
          break;
        case 'd':
          flags|=CHDIFF;
          break;
        case 'D':
          elect.dip_ctr=malloc(3*sizeof(double));
	  if (!elect.dip_ctr) error_exit("malloc error for dip_ctr");
          if(((*(optp+1)>='a')&&(*(optp+1)<='c'))||(*(optp+1)=='m')){
            elect.dip_corr_dir=malloc(1);
            if (!elect.dip_corr_dir) error_exit("Malloc error for char!");
            *elect.dip_corr_dir=*(optp+1);
            optp++;
          }
          if(*(optp+1)=='='){
            i=sscanf(optp+2,"%lf,%lf,%lf",elect.dip_ctr,elect.dip_ctr+1,
                     elect.dip_ctr+2);
            if ((i==0)||(i==EOF))
              elect.dip_ctr[0]=elect.dip_ctr[1]=elect.dip_ctr[2]=0.5;
            else if (i!=3) error_exit("malformed option -D=");
            while(*((++optp)+1));
          }
          else error_exit("malformed option -D");
          break;
        case 'e':
          if(*(optp+1)=='='){
            i=sscanf(optp+2,"%lf-%lf",&tolmin,&tol);
            if (i==1) tol=tolmin;
            else if (i!=2) error_exit("malformed option -e=");
            while(*((++optp)+1));
          }
          else flags|=CST_ESP;
          break;
        case 'E':
          calc_esp=1;
          flags|=CHDEN;
          musq=0;
          j=0;
          if(*(optp+1)=='='){
            if(*(optp+2)=='-'){
              j=1;
              optp++;
            }
            sscanf(optp+2,"%lf",&musq);
            while(*((++optp)+1));
          }
          if (musq==0) musq=C2X_MU; /* Reasonable default? */
          musq=musq*musq;
          if (j) musq*=-1;
          break;
        case 'f':
          failure=1;
          break;
        case 'H':
          half_shift=1;
          break;
        case 'I':
          flags|=BANDPARITY;
          if(*(optp+1)=='='){
            elect.band_range=optp+2;
            while(*((++optp)+1));
          }
          break;
        case 'i':
          if(*(optp+1)=='='){
            i_grid=malloc(3*sizeof(int));
            if (!i_grid) error_exit("malloc error for three ints!");
            if ((sscanf(optp+2,"%d,%d,%d",i_grid,i_grid+1,i_grid+2)!=3)&&
		(sscanf(optp+2,"%dx%dx%d",i_grid,i_grid+1,i_grid+2)!=3))
              error_exit("malformed option -i=");
            while(*((++optp)+1));
          }
          else flags|=BANDIMAG;
          break;
        case 'K':
          sym_k=1;
          break;
        case 'k':
          if(*(optp+1)=='='){
            elect.kpt_range=optp+2;
            while(*((++optp)+1));
          }
          else error_exit("malformed option -k");
          break;
        case 'L':
          flags|=LHS_FUDGE;
          break;
        case 'l':
          no_mp=1;
	  dict_strcat(motif.dict,"CIF_site_label","");	  
          break;
        case 'M':
	  gen_mp=1;
          if (format!=FDF_BP) no_mp=1;
          if(*(optp+1)=='='){
	    new_mp=malloc(6*sizeof(double));
            if (!new_mp) error_exit("malloc error for six doubles!");
            i=sscanf(optp+2,"%lf,%lf,%lf,%lf,%lf,%lf",
		     new_mp,new_mp+1,new_mp+2,
		     new_mp+3,new_mp+4,new_mp+5);
	    if (i==3) new_mp[3]=new_mp[4]=new_mp[5]=0.0;
	    else if (i!=6) error_exit("malformed option -M=");
            while(*((++optp)+1));
          }
          break;
        case 'm':
          molecule=1;
          if(*(optp+1)=='='){
	    if (*(optp+2)=='('){
	      m_rel=malloc(3*sizeof(double));
	      if (!m_rel) error_exit("malloc error for three doubles!");
	      if (sscanf(optp+3,"%lf,%lf,%lf",m_rel,m_rel+1,m_rel+2)!=3)
		error_exit("malformed option -m=(");
              molecule=0;
	    }
	    else{
	      m_abc=malloc(3*sizeof(int));
	      if (!m_abc) error_exit("malloc error for three ints!");
	      if (sscanf(optp+2,"%d,%d,%d",m_abc,m_abc+1,m_abc+2)!=3)
		error_exit("malformed option -m=");
	    }
            while(*((++optp)+1));
          }
          break;
        case 'n':
          no_sym++;
          break;
        case 'N':
          reduce=1;
          break;
        case 'O':
          flags|=OCCUPANCIES;
          pr_occ=1;
          break;
        case 'p':
          flags|=BANDPHASE;
	  dict_strcat(motif.dict,"PDB_has_charge","");
          break;
        case 'P':
          if(*(optp+1)=='='){
            line_spec=malloc(strlen(optp+2)+1);
            if (!line_spec) error_exit("Malloc error for line_spec");
            strcpy(line_spec,optp+2);
            while(*((++optp)+1));
          }
          else prim++;
          break;
        case 'q':
	  /* All these corrections need a density read */
          flags|=CHDEN;
	  dict_strcat(motif.dict,"charge_correction","");
	  elect.dip_corr_dir=malloc(1);
	  if (!elect.dip_corr_dir) error_exit("Malloc error for char!");
	  *elect.dip_corr_dir='m';
	  if ((*(optp+1)>='a')&&(*(optp+1)<='c')){
	    *elect.dip_corr_dir=*(optp+1);
	    optp++;
	  }	    
	  break;	  
        case 'Q':
          sort_style=1;
          if ((*(optp+1)>='0')&&(*(optp+1)<='9')){
            optp++;
            sort_style=*optp-'0';
          }
          break;
        case 'r':
          flags|=BANDREAL;
	  dict_strcat(motif.dict,"CIF_symmetrise","");
          break;
        case 'R':
          flags|=RAW;
          if (*(optp+1)=='='){
	    cptr=optp+2;
            while(*((++cptr)+1));
            if (*cptr=='x') {
	      flags&=(~RAW);
	      *cptr=0;
	    }
            if (single_scan(optp+2,&rescale,NULL)!=1)
              fprintf(stderr,"Failed to parse -%s\n",optp);
	    while(*((++optp)+1));
          }            
          break;
        case 's':
          flags|=SPINDEN;
          break;
        case 'S':
          if(*(optp+1)=='='){
            elect.spin_range=optp+2;
            while(*((++optp)+1));
          }
          else error_exit("malformed option -S\n");
          break;
        case 't':
          new_cell_rel[2][0]=new_cell_rel[2][1]=new_cell_rel[2][2]=0.0;
          if(*(optp+1)=='='){
            if ((sscanf(optp+2,"(%lf,%lf,%lf)(%lf,%lf,%lf)(%lf,%lf,%lf)",
              &new_cell_rel[0][0],&new_cell_rel[0][1],&new_cell_rel[0][2],
	      &new_cell_rel[1][0],&new_cell_rel[1][1],&new_cell_rel[1][2],
              &new_cell_rel[2][0],&new_cell_rel[2][1],&new_cell_rel[2][2])!=9)
              &&(sscanf(optp+2,"(%lf,%lf,%lf)(%lf,%lf,%lf)",
              &new_cell_rel[0][0],&new_cell_rel[0][1],&new_cell_rel[0][2],
	      &new_cell_rel[1][0],&new_cell_rel[1][1],&new_cell_rel[1][2])!=6))
               error_exit("malformed option -t=");
             while(*((++optp)+1));
             expand=5;
             break;
          }
          else error_exit("malformed option -t=");
        case 'T':
          new_cell[2][0]=new_cell[2][1]=new_cell[2][2]=0.0;
          if(*(optp+1)=='='){
            if ((sscanf(optp+2,"(%lf,%lf,%lf)(%lf,%lf,%lf)(%lf,%lf,%lf)",
                       &new_cell[0][0],&new_cell[0][1],&new_cell[0][2],
		       &new_cell[1][0],&new_cell[1][1],&new_cell[1][2],
		       &new_cell[2][0],&new_cell[2][1],&new_cell[2][2])!=9)
                       &&(sscanf(optp+2,"(%lf,%lf,%lf)(%lf,%lf,%lf)",
                       &new_cell[0][0],&new_cell[0][1],&new_cell[0][2],
                       &new_cell[1][0],&new_cell[1][1],&new_cell[1][2])!=6))
               error_exit("malformed option -T=");
             while(*((++optp)+1));
             expand=4;
             break;
          }
          else error_exit("malformed option -T=");
        case 'u':
          flags|=AU;
          break;
        case 'U':
          flags|=DE_AU;
          break;
        case 'w':
	  if ((*(optp+1)=='=')&&(*(optp+2)=='k')){
	    optp+=2;
	    flags|=K_WEIGHT;
	  }
	  else
	    flags|=OCC_WEIGHT;
          break;
        case 'W':
          flags|=OCC_WEIGHT;
          flags|=K_WEIGHT;
          break;
        case 'x':
          if(*(optp+1)=='='){
            if (sscanf(optp+2,"(%lf,%lf,%lf)(%lf,%lf,%lf)(%lf,%lf,%lf)",
	         &new_cell_rel[0][0],&new_cell_rel[0][1],&new_cell_rel[0][2],
                 &new_cell_rel[1][0],&new_cell_rel[1][1],&new_cell_rel[1][2],
	         &new_cell_rel[2][0],&new_cell_rel[2][1],&new_cell_rel[2][2])
		==9){
	      while(*((++optp)+1));
	      expand=2;
	      break;
	    }
	    else if (sscanf(optp+2,"%dx%dx%d",
			    i_expand,i_expand+1,i_expand+2)==3){
	      while(*((++optp)+1));
	      expand=6;
	      break;
	    }	      
	    else error_exit("malformed option -x=");
          }
          else error_exit ("malformed option -x=");
        case 'X':
          if(*(optp+1)=='='){
            if (sscanf(optp+2,"(%lf,%lf,%lf)(%lf,%lf,%lf)(%lf,%lf,%lf)",
                       &new_cell[0][0],&new_cell[0][1],&new_cell[0][2],
                       &new_cell[1][0],&new_cell[1][1],&new_cell[1][2],
                       &new_cell[2][0],&new_cell[2][1],&new_cell[2][2])!=9)
               error_exit("malformed option -X=");
            while(*((++optp)+1));
            expand=3;
            break;
          }
          else if ((*(optp+1)>='a')&&(*(optp+1)<='c')){
            mask=0;
            while((*(optp+1)>='a')&&(*(optp+1)<='c')){
              if (*(optp+1)=='a') mask|=1;
              if (*(optp+1)=='b') mask|=2;
              if (*(optp+1)=='c') mask|=4;
              optp++;
            }
            if (debug>2) fprintf(stderr,"Mask=%d\n",mask);
            if(*(optp+1)!='=') error_exit("malformed option -X=");
            if (sscanf(optp+2,"%lf",&dtmp)!=1)
              error_exit("malformed option -X=");
            while(*((++optp)+1));
	    scale=1;
	    if (*optp=='B') scale=BOHR;
	    if ((*(optp-1)=='n')&&(*optp=='m')) scale=10;
            if (mask&1) new_abc[0]=dtmp*scale;
            if (mask&2) new_abc[1]=dtmp*scale;
            if (mask&4) new_abc[2]=dtmp*scale;
            vexpand=1;
            break;
          }
          else error_exit("malformed option -X=");
        case 'y':
          if(*(optp+1)=='='){
            if (sscanf(optp+2,"%d,%d:%lf",circ,circ+1,&tube_spacing)==3){
	      while(*((++optp)+1));
	      scale=1;
	      if (*optp=='B') scale=BOHR;
	      if ((*(optp-1)=='n')&&(*optp=='m')) scale=10;
	      tube_spacing*=scale;
	    }
	    else if (sscanf(optp+2,"%d,%d",circ,circ+1)!=2){
	      error_exit("malformed option -y=");
            }
	    else{
	      while(*((++optp)+1));
	    }
	    circ[2]=0;
            break;
          }
          else error_exit("malformed option -y=");
        case 'z':
        case 'Z':
          if (*optp=='Z') xc=1;
          if(*(optp+1)=='='){
	    pt_spec=malloc(strlen(optp+2)+1);
            if (!pt_spec) error_exit("Malloc error for pt_spec");
            strcpy(pt_spec,optp+2);
            while(*((++optp)+1));
            format=CNULL;
            break;
          }
          else error_exit("malformed option -z=");
        case '1':
          if (*(optp+1)=='5'){
	    flags|=HIPREC;
            optp++;
            break;
	  }
          else dict_strcat(motif.dict,"cell_is_onetep","");
        case '3':
          preserve_c+=2; /* Usually mangle axes 1 and 2. Add 2 (mod 3) */
          break;         /* to mangle 0 and 1, hence preserving c */
        default:
          fprintf(stderr,"Invalid option %c.\n%s -h for usage.\n",
                   *optp,argv[0]);
           exit(1);
         }
       }
       opt++;
       optp=argv[opt];
       break;
     default:
       file1=argv[opt];
       firstfile_arg=opt;
       break;
     }
  }

  dict_add(motif.dict,"i_frame",i_frame);
  
  if (!file1) error_exit("no input file specified.");

  if ((op==0)&&(opt<argc-1)) outname=argv[argc-1];

  if ((op!=0)&&(format!=CNULL)&&(opt<argc-2)) outname=argv[argc-1];
  
  flags|=((preserve_c%3)<<PC_SHIFT);
  if ((calc_esp)&&(flags&RAW)){
     flags&=(~RAW);
     if (rescale==1)
       fprintf(stderr,"-R ignored in combination with -E\n");
  }

  if ((!outname)&&(format==0)&&(spg_op&(CSPG_INT|CSPG_SCH|CSPG_PNT|CSPG_LST)))
    format=CNULL;

  n=argc-firstfile_arg;
  if (outname) n--; /* We have an output file */
  if ((n>1)&&(op==0)) /* If confused, exit before truncating output file */
    error_exit("confused: multiple input files, but no combining operation");

  if ((outname)&&(strcmp(outname,"-"))){
    if (replace_rho)
      outfile=fopen(outname,"r+b");
    else
      outfile=fopen(outname,"wb");
    if(!outfile){
      fprintf(stderr,"Error, unable to open %s for writing.\n",outname);
      exit(1);
    }
  }
  else
    outfile=stdout;

  dict_add(motif.dict,"out_file_handle",outfile);

  if ((outfile==stdout)&&(isatty(fileno(stdout)))&&
      (format!=CNULL)&&(format!=FORT15))
    error_exit("refusing to output to a terminal");

  /* n set to argc-firstfile_arg-(outname!=NULL) a few lines earlier */
  for(j=0;j<n;j++){
    if (j!=0) init_motif(&motif2);
    file1=argv[firstfile_arg+j];
    infile=fopen(file1,"rb");
    if(!infile){
      fprintf(stderr,"Error, unable to open %s for reading.\n",file1);
      exit(1);
    }
    cptr=malloc(strlen(file1)+1);
    if (!cptr) error_exit("Malloc error for filename");
    strcpy(cptr,file1);
    if (j==0)
      dict_add(motif.dict,"in_file",cptr);
    else
      dict_add(motif2.dict,"in_file",cptr);
    if (strchr(file1,'/')){
      i=strrchr(file1,'/')-file1;
      i=i+1;
      cptr=malloc(i+1);
      if (!cptr) error_exit("Malloc error for dirname");
      strncpy(cptr,file1,i);
      cptr[i]=0;
      if (j==0)
	dict_add(motif.dict,"in_dir",cptr);
      else
	dict_add(motif2.dict,"in_dir",cptr);
    }
    if (j==0)
      file_read(file1,infile,&cell,&motif,&kp,&sym,&grid1,&elect,&series,
		i_grid);
    else{
      init_cell(&cell2);
      init_kpts(&kp2);
      init_sym(&sym2);
      init_grid(&grid2);
      init_elect(&elect2);
      init_tseries(&series2);
      
      file_read(file1,infile,&cell2,&motif2,&kp2,&sym2,&grid2,&elect2,
		&series2,i_grid);
      if ((j==1)&&(op&C2X_MERGE)){
	data_combine(C2X_MERGE,&cell,&motif,&kp,&sym,&grid1,&elect,&series,
		     i_grid,&cell2,&motif2,&kp2,&sym2,&grid2,&elect2,&series2);
	if ((op&(~C2X_MERGE))) op=op&(~C2X_MERGE);
      }
      else
	data_combine(op,&cell,&motif,&kp,&sym,&grid1,&elect,&series,i_grid,
		     &cell2,&motif2,&kp2,&sym2,&grid2,&elect2,&series2);
    }
  }
  
  
  /* Should we also try to read a .cell file? */
  i=strlen(file1);
  if (((op&C2X_MERGE)==0)&&(((i>=7)&&(!strcmp(file1+i-7,"den_fmt")))||
			    ((i>=7)&&(!strcmp(file1+i-7,"pot_fmt")))||
			    ((i>=7)&&(!strcmp(file1+i-7,"elf_fmt")))||
			    ((i>=7)&&(!strcmp(file1+i-7,"cst_esp")))||
			    ((i>=4)&&(!strcmp(file1+i-4,".elf"))))){
    cdfile=malloc(i+5);
    if (!cdfile) error_exit("Malloc error for .cell filename");
    strcpy(cdfile,file1);
    cptr=cdfile+strlen(cdfile);
    while((cptr>cdfile)&&(*cptr!='.')) cptr--;
    if (*cptr!='.') error_exit("Error constructing .cell filename");
    strcpy(cptr+1,"cell");
    infile=fopen(cdfile,"rb");
    if (infile) {
      fprintf(stderr,"Also reading %s\n",cdfile);
      cell_read(infile,&cell,&motif,&kp,&sym);
      fclose(infile);
    }
  }

  /* Check that we have most of what we need */

  if (!cell.basis){
    if ((format!=NPY)&&(format!=CNULL)&&(!replace_rho))
      error_exit("no basis set found. Was format correctly detected?\n");
  }
  
  if (cell.vol<0) cell.vol=fabs(cell.vol);
  if (!keep_velocities) motif.velocities=0;

  cell_check(&cell,&motif);
  
  if (!motif.atoms){
    if ((format!=DENFMT)&&(format!=CUBE)&&(format!=DX)&&
        (format!=FBIN)&&(format!=BXSF)&&(format!=CASTEP_BANDS)&&
	(format!=NPY)&&
	(format!=CNULL)&&(format!=XSF))
      error_exit("no atoms found!");
    else
      if (format==XSF) fprintf(stderr,"Warning: no atoms found\n");
  }
  
  ionic_charge=0;
  if (motif.atoms)
    for(i=0;i<motif.n;i++) ionic_charge+=motif.atoms[i].chg;

  if (elect.nel==0){
    elect.nel=ionic_charge;
    if (elect.charge) elect.nel-=*elect.charge;
  }
  
  /* We define MP grid to include origin, and displacements to be in terms of
   * mesh cells. Castep defines displacements in reciprocal space cells, and
   * has a shift of half a mesh cell if the MP parameter is even.
   */

  if (new_mp){
    if (!kp.mp) kp.mp=malloc(sizeof(struct mp_grid));
    if (!kp.mp) error_exit("Malloc error\n");
    for(i=0;i<3;i++){
      kp.mp->grid[i]=new_mp[i];
    }
    for(i=0;i<3;i++){
      kp.mp->disp[i]=new_mp[i+3]/kp.mp->grid[i];
      if ((kp.mp->grid[i]&1)==0) /* Grid even */
	kp.mp->disp[i]-=1.0/(2*kp.mp->grid[i]);
    }
  }

  if (debug>1){
    fprintf(stderr,"Basis set\n");
    print_basis(cell.basis);
  }
  if (debug>2){
    fprintf(stderr,"Reciprocal basis set\n");
    print_basis(cell.recip);
  }
  if (debug){
    fprintf(stderr,"Cell volume %f\n",cell.vol);
    fprintf(stderr,"natoms      %d\n",motif.n);
    if (ionic_charge) fprintf(stderr,"Total ionic charge %f\n",ionic_charge);
    print_elect(&elect);

    if ((kp.mp)&&(kp.mp->grid[0]>0)){
      fprintf(stderr,"kpoint MP grid %d %d %d\n",kp.mp->grid[0],
              kp.mp->grid[1],kp.mp->grid[2]);
      fprintf(stderr,"        offset %lf %lf %lf\n",kp.mp->disp[0],
              kp.mp->disp[1],kp.mp->disp[2]);
    }
    if (sym.n) fprintf(stderr,"%d symmetry operations read\n",sym.n);
  }

  if (sym_list&&sym.n){
    for(i=0;i<sym.n;i++) ident_sym(&sym.ops[i],&cell,stderr);
  }
  
  if ((molecule)&&(m_abc)&&(i_grid)&&(!grid1.data)){
    m_rel=malloc(3*sizeof(double));
    if (!m_rel) error_exit("malloc error for three doubles!");
    for(i=0;i<3;i++) m_rel[i]=((double)m_abc[i])/i_grid[i];
    free(m_abc);
    m_abc=NULL;
  }
  
  if ((molecule)||(m_rel)) molecule_fix(m_abc,m_rel,&cell,&motif,&grid1);

  if (reduce) molecule=0;
  
  if (half_shift) xplor_fudge(&cell,&motif,&grid1);
  
  if ((elect.dip_ctr)&&(grid1.data)) dipole(&cell,&motif,&grid1,&elect);

  if (dict_get(motif.dict,"charge_correction"))
    charge_corr(&cell,&motif,&grid1,&elect);
  
  if (debug>2){
    fprintf(stderr,"Ionic positions, fractional\n");
    for(i=0;i<motif.n;i++){
      fprintf(stderr,"%d %f %f %f",motif.atoms[i].atno,motif.atoms[i].frac[0],
	      motif.atoms[i].frac[1],motif.atoms[i].frac[2]);
      if (motif.atoms[i].spin) fprintf(stderr,"  spin=%f",motif.atoms[i].spin);
      if (motif.atoms[i].label) fprintf(stderr,"  label=%s",
                                         motif.atoms[i].label);
      if (motif.atoms[i].chg) fprintf(stderr,"  chg=%f",motif.atoms[i].chg);
      fprintf(stderr,"\n");
    }
    fprintf(stderr,"Ionic positions, absolute\n");
    for(i=0;i<motif.n;i++)
      fprintf(stderr,"%d %f %f %f\n",motif.atoms[i].atno,motif.atoms[i].abs[0],
	      motif.atoms[i].abs[1],motif.atoms[i].abs[2]);
    if(motif.forces){
      fprintf(stderr,"Ionic forces\n");
      for(i=0;i<motif.n;i++)
        fprintf(stderr,"%d %f %f %f\n",motif.atoms[i].atno,
		motif.atoms[i].force[0],motif.atoms[i].force[1],
                motif.atoms[i].force[2]);
    }
  }
  if (debug){
    if (elect.cut_off>0){
      fprintf(stderr,"Requested cut-off energy %g eV\n",elect.cut_off);
      if (debug>1){
        g_cut=2*sqrt(2.0*elect.cut_off/H_eV)/BOHR;
        fprintf(stderr,"2|g_max| %g A^-1\n",g_cut);
        cart2abc(&cell,NULL,abc,NULL,0);
        fprintf(stderr,"Full FFT grid based on 2|g_max|: %ix%ix%i\n",
                2*(int)(ceil(g_cut*abc[0]/(2*M_PI)))+1,
                2*(int)(ceil(g_cut*abc[1]/(2*M_PI)))+1,
                2*(int)(ceil(g_cut*abc[2]/(2*M_PI)))+1);
      }
    }
    if ((debug>1)||(grid1.size[0])){
      fprintf(stderr,"First FFT grid     %d %d %d\n",
                      grid1.size[0],grid1.size[1],grid1.size[2]);
      fprintf(stderr,"spins=%d   spinors=%d\n",elect.nspins,elect.nspinors);
    }
  }
  gptr=&grid1;
  while((gptr)&&(gptr->data)){
    if ((calc_esp==0)&&(rescale!=1)){
      if (debug) fprintf(stderr,"Rescaling data by %lf\n",rescale);
      for(i=0;i<gptr->size[0]*gptr->size[1]*gptr->size[2];i++)
        gptr->data[i]*=rescale;
    }
    if (debug){
      double dmin,dmax,sum,sum_abs;
      dmin=1e20;
      dmax=-1e20;
      sum=sum_abs=0;

      fprintf(stderr,"Found 3D data for %s\n",gptr->name);
      for(i=0;i<gptr->size[0]*gptr->size[1]*gptr->size[2];i++){
	sum+=gptr->data[i];
        sum_abs+=fabs(gptr->data[i]);
	if (gptr->data[i]>dmax) dmax=gptr->data[i];
	if (gptr->data[i]<dmin) dmin=gptr->data[i];
      }
      fprintf(stderr,"  min=%g max=%g sum=%g int=%g",dmin,dmax,
	      sum,sum*cell.vol/(gptr->size[0]*gptr->size[1]*gptr->size[2]));
      if ((gptr->name)&&(strstr(gptr->name,"Spin")))
        fprintf(stderr," int|s|=%g",
                sum_abs*cell.vol/(gptr->size[0]*gptr->size[1]*gptr->size[2]));
      fprintf(stderr,
	      "\n  (integral is e per cell for charge and spin densities)\n");
    }

    if (i_grid){
      struct grid ng;

      /* If interpolation is to coarser grid, need to calc es pot first,
         but if it is to a finer grid, need to interpolate first. */
      if ((calc_esp)&&(gptr==&grid1)){
        for(i=0;i<3;i++) ng.size[i]=max(i_grid[i],gptr->size[i]);
        if ((ng.size[0]!=gptr->size[0])||
            (ng.size[1]!=gptr->size[1])||
            (ng.size[2]!=gptr->size[2])){
          interpolate3d(gptr,&ng);
          free(gptr->data);
          gptr->data=ng.data;
          for(i=0;i<3;i++) gptr->size[i]=ng.size[i];
        }
        es_pot(&cell,&motif,&grid1,&elect,musq);
        calc_esp=0;
	if (rescale!=1){
	  if (debug) fprintf(stderr,"Rescaling calculated potential by %lf\n",
			     rescale);
	  for(i=0;i<grid1.size[0]*grid1.size[1]*grid1.size[2];i++)
	    grid1.data[i]*=rescale;
	}
      }

      for(i=0;i<3;i++) ng.size[i]=i_grid[i];

      interpolate3d(gptr,&ng);

      free(gptr->data);
      gptr->data=ng.data;
      for(i=0;i<3;i++) gptr->size[i]=ng.size[i];

      if (debug){
	double min,max,sum;
	min=1e20;
	max=-1e20;
	sum=0;
	for(i=0;i<gptr->size[0]*gptr->size[1]*gptr->size[2];i++){
	  sum=sum+gptr->data[i];
	  if (gptr->data[i]>max) max=gptr->data[i];
	  if (gptr->data[i]<min) min=gptr->data[i];
	}
	fprintf(stderr,"After interpolation:\n");
	fprintf(stderr,"  min=%g  max=%g  sum=%g  int=%g\n",min,max,
                sum,sum*cell.vol/(gptr->size[0]*gptr->size[1]*gptr->size[2]));
      }
      
    }
    gptr=gptr->next;
  }

  if (calc_esp) {
    if ((i_grid)&&(!grid1.data)){
      for(i=0;i<3;i++) grid1.size[i]=i_grid[i];
    }
    es_pot(&cell,&motif,&grid1,&elect,musq);
    if (rescale!=1){
      if (debug) fprintf(stderr,"Rescaling calculated potential by %lf\n",
			 rescale);
      gptr=&grid1;
      for(i=0;i<gptr->size[0]*gptr->size[1]*gptr->size[2];i++)
        gptr->data[i]*=rescale;
    }
  }
    
  if (pt_spec){
    lscan(&pt_spec,&motif,zpt);
    gptr=&grid1;
    while(gptr&&(gptr->data)){
      dtmp=interpolate0d(gptr,zpt);
      printf("%s at (%f,%f,%f): %e\n",gptr->name,zpt[0],zpt[1],zpt[2],
             dtmp);
      if (xc) printf("PZ XC potential at %f eA^-3 = %f V\n",dtmp,lda(dtmp));
      gptr=gptr->next;
    }
  }

  if (gen_mp){
    if (kp.mp) mp_gen(&kp,&cell);
    else fprintf(stderr,"Warning: ignoring -M as no MP parameters found\n");
  }

  
  if (failure) fstar(&kp, &cell, &sym);

  if (reduce) reduce_cell(motif.atoms,motif.n,cell.basis);
  
  if (spg_op) cspg_op(&cell,&motif,&sym,&kp,spg_op,tolmin);

  if ((sym_k)&&(kp.n)){
    struct kpts kp2;
    struct symmetry ksym;

#ifdef SPGLIB
    if (sym.n==0){
      fprintf(stderr,"Calling spglib to determine symmetry\n");
      cspg_op(&cell,&motif,&sym,&kp,CSPG_SYM,tolmin);
    }
#endif
    
    ksym.n=0;
    sym2ksym(&sym,&ksym);
    
    kp2.n=0;
    addabs(kp.kpts,kp.n,cell.recip);
    sym_kpts(&kp,&kp2,&ksym,cell.basis);
    if (kp2.n){
      free(kp.kpts);
      kp.kpts=kp2.kpts;
      if (debug>1) fprintf(stderr,
                           "Before symmetrisation %d kpoints, after %d\n",
                           kp.n,kp2.n);
      kp.n=kp2.n;
    }
    free(ksym.ops);
  }

  
  if (prim) {
    double nb1[3][3],ob[3][3];
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
        ob[i][j]=cell.basis[i][j];
    primitive(&cell,&motif,nb1);
    if (debug>2){
      fprintf(stderr,"Primitive cell is:\n");
      for(i=0;i<3;i++){
        fprintf(stderr,"(%f,%f,%f)\n",nb1[i][0],nb1[i][1],nb1[i][2]);
      }
    }
    super(&cell,&motif,nb1,&kp,&sym,&grid1,0);
    if (debug){
      fprintf(stderr,"Old basis in terms of new: ");
      old_in_new(ob,cell.recip);
    }
  }

  if (compact){
    double nb1[3][3];
    for (i=0;i<3;i++)
      for(j=0;j<3;j++)
        nb1[i][j]=cell.basis[i][j];

    shorten(nb1);
    super(&cell,&motif,nb1,&kp,&sym,&grid1,0);
  }


  if (expand){
    switch(expand){
    case 6:
      simple_super(&cell,&motif,i_expand,&kp,&sym,&grid1);
      break;
    case 5: /* First vector given in relative terms */
      for(i=0;i<3;i++){
        new_cell[0][i]=new_cell_rel[0][0]*cell.basis[0][i]+
                       new_cell_rel[0][1]*cell.basis[1][i]+
                       new_cell_rel[0][2]*cell.basis[2][i];
        new_cell[1][i]=new_cell_rel[1][i];
        new_cell[2][i]=new_cell_rel[2][i];
      }
      rotation(&cell,&motif,new_cell);
      break;
    case 4:
      rotation(&cell,&motif,new_cell);
      break;
    case 3:  /* New basis given explicitly in absolute terms */
      nc.basis=new_cell;
      cart2abc(&nc,NULL,abc,&grid1,0);
      super(&cell,&motif,new_cell,&kp,&sym,&grid1,0);
      break;
    case 2:  /* New basis given explicitly in relative terms */
      for(i=0;i<3;i++)
        for(j=0;j<3;j++){
          new_cell[i][j]=0;
          for(k=0;k<3;k++)
            new_cell[i][j]+=new_cell_rel[i][k]*cell.basis[k][j];
        }
      nc.basis=new_cell;
      cart2abc(&nc,NULL,abc,&grid1,0);
      super(&cell,&motif,new_cell,&kp,&sym,&grid1,0);
      break;
    case 1:  /* We were expected to guess the appropriate new basis */
      fprintf(stderr,"-x: feature not present\n");
      break;
    default:
      error_exit("Internal error in expansion selection");
    } /* select(expand) ... */
  }  /* if (expand) ... */

  if (vexpand){
    vacuum_adjust(&cell,&motif,new_abc);
    gptr=&grid1;
    while (gptr->data){
      free(gptr->data);
      gptr->data=NULL;
      if (gptr->next){
	struct grid *g;
	g=gptr->next;
	gptr->next=NULL;
	gptr=g;
      }
    }
  }

  if ((circ[0])||(circ[1])||(circ[2]))
    tube(&cell,&motif,circ,tube_spacing);
  
  if (rotate) cart2abc(&cell,&motif,abc,&grid1,1);
  if ((no_mp)&&(kp.mp)) {free(kp.mp);kp.mp=NULL;}
  if (sort_style) sort_atoms(&motif,sort_style);
  if (no_sym){
    sym.n=0;
    if (sym.ops) free(sym.ops);
    if (format==XSF) motif.forces=0;
  }
  if (no_sym==2){
    kp.n=0;
    if (kp.mp){
      free(kp.mp);
      kp.mp=NULL;
    }
  }

  if (calc_ef){
    if (elect.eval){
      chg=elect.nel;
      if (chg==0){
	for(i=0;i<motif.n;i++)
	  chg+=motif.atoms[i].chg;
	if ((chg)&&(elect.charge)) chg-=*elect.charge;
	if (chg)
	  fprintf(stderr,"Assuming nel=%lf from total ionic charge\n",chg);
      }
      if (chg){
	if (!elect.e_fermi) elect.e_fermi=malloc(sizeof(double));
	*elect.e_fermi=calc_efermi(&elect,&kp,chg);
	fprintf(stderr,"Calculated E_Fermi (0K): ");
	print_energy(*elect.e_fermi);
	fprintf(stderr,"\n");
      }
      else
	fprintf(stderr,"Unable to calculate E_Fermi"
		" as unable to determine number of electrons\n");
    }
    else
	fprintf(stderr,"Unable to calculate E_Fermi"
		" as no eigenvalues read\n");
  }
  
  if (pr_gap) print_gap(&elect, &kp,&motif);
  if (pr_occ) print_occ(&elect, &kp);
  
  if (line_spec){
    line_write(outfile,&cell,&motif,&grid1,line_spec,lflags);
    exit(0);
  }

  if (replace_rho){
    flags|=(REPLACE_RHO+CHDEN);
    check_read(outfile,&cell,&motif,&kp,&sym,&grid1,&elect,NULL);
    exit(0);
  }
  
  switch(format){
    case XSF:
      xsf_write(outfile,&cell,&motif,molecule,&grid1);
      break;
    case BXSF:
      bxsf_write(outfile,&cell,&motif,&elect,&kp,&sym);
      break;
    case CUBE:
      cube_write(outfile,&cell,&motif,&grid1);
      break;
    case DENFMT:
      denfmt_write(outfile,&cell,&motif,&grid1);
      break;
    case XPLOR:
      xplor_write(outfile,&cell,&motif,&grid1);
      break;
    case PDB:
      pdb_write(outfile,&cell,&motif);
      break;
    case CELL:
      cell_write(outfile,&cell,&motif,&kp,&sym);
      break;
    case CELL_ABC:
      cell_write_abc(outfile,&cell,&motif,&kp,&sym);
      break;
    case CELL_ABS:
      cell_write_abs(outfile,&cell,&motif,&kp,&sym);
      break;
    case CELL_ABC_ABS:
      cell_write_abc_abs(outfile,&cell,&motif,&kp,&sym);
      break;
    case DX:
      dx_write(outfile,&cell,&grid1);
      break;
    case VASP:
      vasp_write(outfile,&cell,&motif,&grid1);
      break;
    case XYZ:
      xyz_write(outfile,&cell,&motif,&series);
      break;
    case CML:
      cml_write(outfile,&cell,&motif);
      break;
    case FDF:
      fdf_write(outfile,outname,&cell,&motif,&kp,&grid1,&elect);
      break;
    case FDF_BP:
      fdf_write_bp(outfile,&kp,&motif);
      break;
    case SHELX:
      shelx_write(outfile,&cell,&motif);
      break;
    case FBIN:
      fbin_write(outfile,&grid1);
      break;
    case FORT15:
      f15_write(&cell,&motif,&kp);
      break;
    case CIF:
      cif_write(outfile,&cell,&motif,&sym,0);
      break;
    case MMCIF:
      cif_write(outfile,&cell,&motif,&sym,1);
      break;
    case NPY:
      npy_write(outfile,&grid1);
      break;
    case PY:
      py_write(outfile,&cell,&motif,'d');
      break;
    case PYA:
      py_write(outfile,&cell,&motif,'a');
      break;
    case ABINIT:
      abinit_write(outfile,&cell,&motif,&kp,&sym,&elect);
      break;
    case QE:
      qe_write(outfile,&cell,&motif,&kp,&elect);
      break;
    case ELK:
      elk_write(outfile,&cell,&motif,&kp,&elect);
      break;
    case CASTEP_BANDS:
      bands_write(outfile,&cell,&kp,&elect);
      break;
    case CASTEP_GEOM:
      geom_write(outfile,&cell,&motif,&elect,&series);
      break;
    case XV:
      xv_write(outfile,outname,&cell,&motif,&grid1);
      break;
    case CCP4:
      ccp4_write(outfile,&cell,&grid1);
      break;
    case CNULL:
      break;

    default:
      fprintf(stderr,"This cannot happen. Sorry\n");
  }
  
  if (outfile!=stdout) fclose(outfile);
  exit(0); /* valgrind prefers this to return(0); */
}

void help(void){
  printf("Usage: c2x [-aAbBckmRsSvx] [--FORMAT] "
#ifdef SPGLIB
         "[--OPERATION] "
#endif
         "infile [outfile]\n\n"
         "-a           rotate as though outputing in abc format\n"
         "-A           accumulate (sum) requested bands\n"
         "-b[=range]   include bands (as psi, rescaled to A^-1.5)\n"
         "-B[=range]   include bands (as densities, rescaled to eA^-3)\n"
         "-c           include charge density (rescaled to eA^-3)\n"
         "--calc       evaluate next arg with arithmetic parser and exit\n"
         "--constants  report internal conversion constants and exit\n"
         "-C           find compact set of cell vectors\n"
         "-d           read also a .chdiff file, constructing its name from\n"
         "               the .cell or .check file given\n"
         "-D=[x,y,z]   calculate dipole moment about fractional co-ordinates\n"
         "               x,y,z (default .5,.5,.5) assuming -c also given\n"
         "-Da=[x,y,z]  ditto, and apply post-hoc slab energy correction with\n"
         "               perpendicular axis a. Valid values a, b and c.\n"
         "-Dm=[x,y,z]  ditto, but apply post-hoc correction for molecule in "
         "cubic box\n"
         "-e           read also a .cst_esp file, constructing its name from\n"
         "               the .cell or .check file given\n"
         "-e=eps       set tolerance for supercell operations. Default=%g\n"
         "-e=min-max   set tolerance range for spglib operations\n"
         "-E[=[-][mu]] calculate electrostatic potential. Default mu=%g\n"
         "-f           find first failure star of k-point mesh\n"
         "--formats    list supported file formats and exit\n"
	 "--frame=X    extract frame X (initial frame is 0) from timeseries\n"
	 "--gap        print band gap\n"
         "-H           shift atoms by half a grid cell\n",tol,C2X_MU);
  printf("-i           output imaginary part of band\n"
         "-i=n1,n2,n3  Fourier interpolate 3D grids to specified grid\n"
         "-I[=range]   report inversion symmetries of given bands\n"
         "-k=range     include given k-points (default 1) for bands\n"
         "-K           symmetrise k-point list\n"
         "-l           list k-points in .cell output, not MP parameters\n"
         "-L           produce (incorrect) left-handed abc output\n");
  printf("-m[=a,b,c]   assume input is molecule, not crystal, and move by\n"
         "               given nos of grid cells, or move automatically\n"
	 "-m=(x,y,z)   shift by given relative co-ords, rounded to nearest\n"
	 "               integer no of grid cells if grid present\n"
         "-M           generate Monkhurst-Pack k-point set from .cell values\n"
         "-M=nx,ny,nz  generate regular k-point mesh including origin\n"
         "-M=nx,ny,nx,dx,dy,dz\n"
         "               ditto displaced by fraction of mesh cell\n"
         "-n           discard symmetry information (and forces if XSF"
         " output)\n"
         "               give twice to discard k-points too\n"
         "-N           normalise by reducing fractional coords to 0<=x<1\n"
	 "               with -m, do write cell into .xsf file\n"
         "-O           print band occupancies and eigenvalues\n"
	 "-p           partial charges in cols 61-66 of PDB input\n");
  printf("-P           find primitive cell\n"
         "-P=X:Y:npts  output 1D data with npts points along line\n"
         "               X and Y either (x,y,z) as fractional co-ords\n"
         "               or atomic symbol followed by atom number\n"
         "               e.g. (0.5,0.5,0.5):C3:20 -- centre to 3rd C atom, "
         "20 points\n"
	 "-P=C:rl:npts ditto as centre, radius, points for cylindrical"
	 " symmetry with c\n"
	 "               cylindrical axis. E.g. (0.5,0.5,0.5):r10B:100\n"
         "-q           post hoc energy correction for charged cells (0D)\n"
         "-q[abc]      post hoc energy correction for charged cells, 2D\n"
	 "               non-periodic axis as specified\n"
         "-Q[n]        quicksort atoms in descending atomic order "
         "(n=1 or absent)\n"
         "                                ascending atomic order  (n=2)\n"
	 "-r           only reduced (symmetry inequivalent) atoms in cif"
	 " output\n"
         "-R           don't rescale grid data or bands, don't adjust"
	 " nanotube radius\n"
         "-R=x         rescale grid data by given factor, suffix with an x "
         "to include\n"
         "               c2x's usual rescaling\n"
         "-s           include spin densities\n"
	 "--sym_list   list symmetry operations read without calling spglib\n"
         "-S=range     include given spins/spinors (0 or 1) for bands\n");
  printf("-t=(x1,y1,z1)(x2,y2,z2)[(x3,y3,z3)]\n"
         "             rotate coords so 1st vector becomes 2nd, using third\n"
         "             vector as axis, else axis is perpendicular to others\n"
         "-T=(x1,y1,z1)(x2,y2,z2)[(x3,y3,z3)]\n"
         "             ditto, but first vec expressed in absolute coords\n"
         "-u           write .cell files and 1D axes in Bohr (atomic Units)\n"
         "             scale densities on writing .cube files from A^-3 to "
         "Bohr^-3\n"
	 "             write forces in .xsf files in Ha/A, not eV/A\n"
         "-U           scale densities on reading .cube files from "
         "Bohr^-3 to A^-3\n"
         "-v           be verbose (may be repeated)\n"
         "--vel        keep velocities\n"
	 "--version    report version and exit, also conversion factors if "
	 "preceded by -v\n"
         "-w           weight band by occupancy\n"
         "-w=k         weight band by kpt weight\n"
         "-W           weight band by kpt weight and occupancy\n");
  printf("-x=(x1,y1,z1)(x2,y2,z2)(x3,y3,z3)\n"
         "             re-express in new basis given in terms of old\n"
	 "-x=ixjxk     trivial tiling to make a supercell\n"
         "-X=(x1,y1,z1)(x2,y2,z2)(x3,y3,z3)\n"
         "             re-express in new basis given in absolute terms\n"
	 "-X[abc]=x    expand given axes to given length by adjusting vacuum\n"
	 "               length may be suffixed with B (Bohr) or nm, "
         "else A assumed\n"
	 "-y=i,j       form nanotube with circumferencial vector i*a+j*b.\n"
	 "-y=i,j:x     ditto, separate tubes by x (default A, else B or nm)\n"
         "-z=X         print volumetric data at point, "
	 "see -P= for point specification\n"
         "-Z=X         as -z, assume data is density (eA^-3) & report "
         "approx XC potential\n"
         "-3           if swapping axes to convert lhs to rhs, third axis\n"
         "               is special & first two swapped. Else first special.\n"
         "-15          use high precision in output.\n"
         "\n");
  printf("range specifies band and kpoint numbers as \"a,b-c,d\" starting"
	 " from 1\n"
         "-b and -B are mutually exclusive, as are -x and -t.\n\n");
  printf("OPERATION is used with multiple input files and can be one of:\n"
	 "  add or sum    sum corresponding datasets in input files\n"
	 "  sub or diff   difference corresponding datasets in two input"
	 " files\n"
	 "  mask or mult  multiply corresponding datasets in two input files\n"
	 "  merge    merge data into single output\n\n");
#ifdef SPGLIB
  printf("OPERATION is also used to call spglib and is one of:\n"
         "   primitive         call spg_find_primitive()\n"
         "   primitive_nr      call spg_standardize_cell(to_primitive=1, "
         "no_idealize=1)\n"
         "   refine            call spg_refine_cell()\n"
         "   int               call spg_get_dataset()"
         " and report international symbol\n"
         "   schoen            call spg_get_schoenflies()\n"
	 "   snap              call spg_standardize_cell() then expand back\n"
	 "                       to a snapped version of the original cell\n"
	 "   snap_tr           ditto, but also include any translation\n"
         " standardise or std  call spg_standardize_cell(no_idealize=1)\n"
	 "   std_ideal         call spg_standardize_cell(no_idealize=0)\n"
         "   symmetry or sym   call spg_get_dataset() and keep symmetry ops\n"
         "   list              call spg_get_dataset() and list symmetry ops\n"
         "   point             call spg_get_dataset() followed by "
         "spg_get_pointgroup()\n\n");
#endif
  printf("Valid values of FORMAT are listed by the --formats argument, as are"
	 "\nrecognised input formats.\n\n");
  printf("Version " C2XSF_VER ", (c) MJ Rutter 2007 - 2022"
         " licenced under the GPL v3.\n\n");
  printf("If useful to a published paper, please consider citing using the\n"
         "references shown with the --refs argument.\n\n");
  printf("This version is ");
#ifdef SPGLIB
  printf("linked with spglib version %d.%d.%d\n",
         spg_get_major_version(),
         spg_get_minor_version(),
         spg_get_micro_version());
#else
  printf("not linked with spglib\n");
#endif
   printf("\nFurther documentation at https://www.c2x.org.uk/\n");
  exit(0);
}

void formats(void){
  printf("Recognised output formats are:\n\n"
	 "            abinit    Abinit .in / .abi\n"
	 "            abinit8   Abinit version 8 .in\n"
	 "            bands     unsorted CASTEP .bands file\n"
	 "            bxsf      unsorted XCrysDen bands file\n"
	 "            ccp4      CCP4 (grid data only, no atoms)\n"
         "            cell      CASTEP .cell, cartesian and fractional\n"
         "            cell_abc                abc and fractional\n"
         "            cell_abs                cartesian and absolute\n"
         "            cell_abc_abs            abc and absolute\n"
         "            cif / mmcif  Basic maybe CIF compatible file\n"
         "            chgcar    VASP CHGCAR (implies -c)\n"
         "            cml       Chemical Markup Language\n"
         "            cube      Gaussian cube\n"
         "            denfmt    CASTEP formatted density\n"
         "            dx        OpenDX\n"
	 "            elk       Elk elk.in format\n"
         "            fdf       Flexible Data Format (Siesta),"
         " with .RHO if grid read\n"
	 "            gcoeff    ASCII wavefunction coefficients\n"
         "            gcoeff_sorted          ditto, sorted by |g|\n"
	 "            geom      CASTEP .geom format\n"
         "            gnu       Gnuplot (with -P=)\n"
         "            mocube    single dataset in Gaussian molecular"
         " orbital cube format\n"
	 "            npy       single dataset in NumPy array format, as"
	 " doubles with -15\n"
         "            null      Discard output\n"
         "            one       Onetep .dat, very similar to .cell,\n"
         "                              also one_abc, one_abs, one_abc_abs\n"
         "            pdb       PDB\n"
         "            pdbn      PDB with atoms numbered\n"
	 "            pdbx      Basic maybe PDBX compatible file\n"
         "            py        python dictionary\n"
         "            pya       python ASE Atoms\n"
         "            qe        Quantum Espresso .in\n"
         "            qef            ditto, fractional atomic coords\n"
         "            refs      List BibTeX references for this code "
         "and exit\n"
         "            shelx     SHELX97\n"
         "            vasp      VASP output (POSCAR/CHG)\n"
         "            wavecar   VASP WAVECAR output\n"
         "            wavecar45 ditto, VASP 4.5-compatible\n"
         "            xplor     Xplor\n"
         "            xsf       XCrySDen (default"
         " unless called as cellsym)\n"
         "            xv        Siesta's .XV format,"
         " .RHO also written if grid read\n"
         "            xyz       XYZ\n"
         "            xyze      Extended XYZ (adds lattice vectors)\n\n");
  printf("Input file type is detected by name. Files\n"
	 "            ending .abi are assumed to be Abinit input\n"
         "            ending .bands are assumed to be Castep .bands files\n"
         "            ending .cif or .mmcif are assumed to be in "
            "a cif format,\n"
	 "            ending .cst_esp are assumed to be in Castep potential"
	 " format,\n"
	 "            ending .cub, .cube or _CUBE are assumed to be in "
            "cube format,\n"
	 "            ending .dat are assumed to be in Onetep .dat format\n"
         "            ending .den_fmt are assumed to be in Castep format\n"
	 "            ending .EIG are assumed to be Abinit or Siesta"
	 " eigenvalues\n"
	 "            ending .elf or .elf_fmt are assumed to be in Castep ELF"
	 " format,\n"
	 "            ending .fdf are assumed to be in Siesta format,\n"
         "            ending gcoeff (any case) are assumed to be GCOEFF.txt"
         " format,\n"
	 "            ending .geom are assumed to be Castep .geom format,\n"
	 "            ending .in are assumed to be in either Abinit or\n"
         "              Quantum Espresso format (save for elk.in),\n"
	 "            ending .kp or .KP are assumed to be in Siesta kpoints,\n"
         "            ending .pdb are assumed to be in pdb format,\n"
         "            ending .pdbx are assumed to be in pdbx (mmCIF) format,\n"
         "            ending .pot_fmt are assumed to be in Castep format\n"
         "            ending .res are assumed to be in shelx97 format,\n"
         "            ending .rho or .RHO are assumed to be in Siesta's"
         " RHO format,\n"
	 "            ending .vh, .VH, .vt or .VT are assumed to be"
	 " Siesta potentials\n"
	 "            ending .xml are assumed to be QE PWscf output,\n"
         "            ending .xsf are assumed to be in xsf format,\n"
         "            ending .xv or .XV are assumed to be in Siesta's"
	 " XV format,\n"
	 "            ending .xyz are assumed to be in extended XYZ format,\n"
	 "            ending DEN, POT, VCLMB or WFK are assumed to be in\n"
	 "              Abinit binary format,\n");
  printf("            ending 3D.OUT are assumed to be in Elk 3D plot format\n");
  printf("            ending or beginning CHG, CHGCAR, POSCAR, CONTCAR, LOCPOT,"
         " EIGENVAL\n"
         "              or WAVECAR are assumed to be in VASP 5.x format,\n"
	 "            ending 12 are assumed to be in Crystal format,\n"
	 "            called fort.34 are assumed to be in Crystal binary "
	 "format.\n"
	 "            called GEOMETRY.OUT or elk.in are assumed to be in"
	 " Elk format\n\n"
         "Otherwise "
         "automatic detection of .cell or .check input. Compatible with\n"
         ".check files from CASTEP 3.0 to 20.11 (and perhaps beyond).\n\n");
  printf("\nFurther documentation at https://www.c2x.org.uk/\n");
  exit(0);
}

void refs(void){
  printf("@article{c2x,\n"
         "  title = {C2x: a tool for visualisation and input preparation for "
      "{C}astep and other electronic structure codes},\n"
         "  journal = {Computer Physics Communications},\n"
         "  author = {Rutter, M. J.},\n"
         "  year = {2018},\n"
         "  volume = {225},\n"
         "  pages = {174--179},\n"
         "  doi = \"10.1016/j.cpc.2017.12.008\"\n"
         "}\n");
#ifdef SPGLIB
  printf("\nIf making use of the spglib functionality,"
         " consider also citing\n\n");
  printf("@misc{spglib,\n"
         "  title = {Spglib, a library for finding and handling "
         "crystal symmetries},\n"
         "  author = {Togo, Atsushi},\n"
         "  url = {https://atztogo.github.io/spglib/}\n"
         "}\n");
  printf("or https://arxiv.org/abs/1808.01590\n");
#endif
  printf("\nThe post hoc dipole correction schemes are described by:\n\n");
  printf("Slabs:     https://doi.org/10.1103/PhysRevB.46.16067 and\n"
         "           https://doi.org/10.1103/PhysRevB.59.12301\n");
  printf("Cubes:     https://doi.org/10.1103/PhysRevB.51.4014 and\n"
         "           https://doi.org/10.1103/PhysRevB.60.15476\n");
  printf("Tetragons: https://doi.org/10.1088/1361-648X/ab20e1\n\n");
  printf("\nThe post hoc charge correction schemes are described by:\n\n");
  printf("0D:        https://link.aps.org/doi/10.1103/PhysRevB.51.4014\n");
  printf("2D:        https://doi.org/10.1088/2516-1075/abeda2\n");
  exit(0);
}

void error_exit(char *msg){
  fprintf(stderr,"Aborting: %s\n",msg);
  exit(1);
}

/* Append comment to linked list of comments */
void add_cmt(struct cmt *comment, char *txt){
  struct cmt *c;
  c=comment;
  while (c->next) c=c->next;
  c->next=malloc(sizeof(struct cmt));
  if (!c->next) error_exit("Malloc error in add_cmt");
  c->next->txt=NULL;
  c->next->next=NULL;
  c->txt=malloc(strlen(txt)+1);
  if (!c->txt) error_exit("Malloc error for txt in add_cmt");
  strcpy(c->txt,txt);
  if (c->txt[strlen(c->txt)-1]=='\n') c->txt[strlen(c->txt)-1]=0;
}
