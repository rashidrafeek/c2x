/* Convert CASTEP .check file to XCrySDen .xsf file
 * Extract other data from a .check file 
 * Form supercells, perform symmetry analysis, and
 * Convert between various crystallographic file formats.
 */


/* Copyright (c) 2007 -- 2018 MJ Rutter 
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
void xyz_write(FILE* outfile, struct unit_cell *c, struct contents *m);
void cml_write(FILE* outfile, struct unit_cell *c, struct contents *m);
void fdf_write(FILE* fdf, char* filename, struct unit_cell *c,
               struct contents *m,struct grid *g, struct es *e);
void line_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                struct grid *g, char *line_spec);
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

void cell_read(FILE* infile, struct unit_cell *c, struct contents *m,
               struct kpts *k, struct symmetry *s);
void check_read(FILE* infile, struct unit_cell *c, struct contents *m,
                struct kpts *k, struct symmetry *s, struct grid *g,
                struct es *elect, int *i_grid);
void chdiff_read(FILE* infile, struct grid *g);
void esp_read(FILE* infile, struct grid *g, struct es *elect);
void pdb_read(FILE* infile, struct unit_cell *c, struct contents *m);
void shelx_read(FILE* infile, struct unit_cell *c, struct contents *m);
void cif_read(FILE* infile, struct unit_cell *c, struct contents *m,
              struct symmetry *s);
void cube_read(FILE* infile, struct unit_cell *c, struct contents *m,
               struct grid *gptr);
void xsf_read(FILE* infile, struct unit_cell *c, struct contents *m,
	      struct grid *gptr);
void vasp_read(FILE* infile, struct unit_cell *c, struct contents *m,
                struct grid *gptr, struct es *elect);
void denfmt_read(FILE* infile, struct unit_cell *c, struct grid *gptr);
void fort34_read(FILE* infile, struct unit_cell *c, struct contents *m,
              struct symmetry *s);
void crystal_read(FILE* infile, struct unit_cell *c, struct contents *m,
                  struct symmetry *s);
void abinit_read(FILE* infile, struct unit_cell *c, struct contents *m,
                 struct grid *gptr);
void abinit_in_read(FILE* infile, struct unit_cell *c, struct contents *m,
                 struct kpts *k, struct symmetry *s, struct es *e);
void fdf_read(FILE* in, struct unit_cell *c, struct contents *m,
               struct kpts *kp, struct es *e);

void molecule_fix(int* m_abc, struct unit_cell *c, struct contents *m,
                  struct grid *g);
int super(struct unit_cell *c, struct contents *m,
           double new_basis[3][3], struct kpts *k, struct symmetry *s,
           struct grid *gptr, int rhs);
void rotation(struct unit_cell *c, struct contents *m, double new_basis[3][3]);
void primitive(struct unit_cell *c, struct contents *m, double basis[3][3]);
void shorten(double basis[3][3]);
void interpolate3d(struct grid *old_grid, struct grid *new_grid);
double interpolate0d(struct grid *gptr,double x[3]);
void fstar(struct kpts *ks, struct unit_cell *c, struct symmetry *rs);
int cspq_op(struct unit_cell *c, struct contents *m, struct symmetry *s,
            int op, double tolmin);
void dipole(struct unit_cell *c, struct contents *m,
            struct grid *g, double *dipole_ctr);
void es_pot(struct unit_cell *c, struct contents *m,
            struct grid *g, double musq);

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
  exit(0);
}

int main(int argc, char **argv)
{
  int i,j,k,opt=1,expand,rotate,half_shift=0,no_mp=0,prim=0,compact=0;
  int molecule=0;
  int sort_style=0;
  int no_sym=0;
  int spg_op;
  int gen_mp=0, failure=0,calc_esp=0;
  int format,preserve_c;
  int *i_grid;
  char *optp,*file1=NULL,*file2=NULL,*line_spec=NULL,*pt_spec=NULL;
  char *cdfile,*cptr;
  FILE *infile,*outfile;
  double abc[6],new_cell[3][3],new_cell_rel[3][3],*new_mp,zpt[3],g_cut;
  double tolmin=1e-4,ionic_charge,musq;
  double *dipole_ctr;
  int *m_abc;
  struct grid *gptr;
  struct unit_cell cell,nc;
  struct contents motif;
  struct symmetry sym;
  struct kpts kp;
  struct grid grid1;
  struct es elect;

/* Initialise all pointers to NULL, etc. */

  elect.band_range="-";
  elect.kpt_range="1";
  elect.spin_range="-";
  elect.nspins=1;
  elect.nspinors=1;
  elect.spin_method=NULL;
  elect.cut_off=0;
  elect.etol=0;

  motif.atoms=NULL;
  motif.n=motif.forces=0;
  motif.title=NULL;
  motif.species_misc=NULL;
  motif.block_species=NULL;
  motif.comment=malloc(sizeof(struct cmt));
  motif.comment->txt=NULL;
  motif.comment->next=NULL;
  cell.basis=NULL;
  cell.vol=0;
  sym.tol=NULL;
  sym.ops=NULL;
  sym.n=0;
  m_abc=i_grid=NULL;
  kp.n=0;
  kp.kpts=NULL;
  kp.mp=NULL;
  new_mp=NULL;
  sym.gen=NULL;
  for(i=0;i<3;i++){
    grid1.size[i]=0;
    for(j=0;j<3;j++) cell.recip[i][j]=0.0;
  }
  grid1.next=NULL;
  grid1.data=NULL;
  grid1.name=NULL;
  dipole_ctr=NULL;
  flags=0;

  opt=1;
  optp=argv[opt];
  debug=molecule=format=expand=rotate=preserve_c=0;
  spg_op=0;

  while (opt<argc){
    switch(*optp){
    case 0:
      opt++;
      optp=argv[opt];
      break;
    case '-':
      if(*(optp+1)=='-'){ /* Gnu-style "--" option */
        if (!strcmp(optp,"--xsf")) format=XSF;
        else if (!strcmp(optp,"--cube")) format=CUBE;
        else if (!strcmp(optp,"--xplor")) format=XPLOR;
        else if (!strcmp(optp,"--pdb")) format=PDB;
        else if (!strcmp(optp,"--pdbn")) {format=PDB; flags|=PDBN;}
        else if (!strcmp(optp,"--cell")) format=CELL;
        else if (!strcmp(optp,"--cell_abc")) format=CELL_ABC;
        else if (!strcmp(optp,"--cell_abs")) format=CELL_ABS;
        else if (!strcmp(optp,"--cell_abc_abs")) format=CELL_ABC_ABS;
        else if (!strcmp(optp,"--one"))
          {format=CELL; flags|=ONETEP_OUT;}
        else if (!strcmp(optp,"--one_abc"))
          {format=CELL_ABC; flags|=ONETEP_OUT;}
        else if (!strcmp(optp,"--one_abs")) 
          {format=CELL_ABS; flags|=ONETEP_OUT;}
        else if (!strcmp(optp,"--one_abc_abs"))
          {format=CELL_ABC_ABS; flags|=ONETEP_OUT;}
        else if (!strcmp(optp,"--dx")) format=DX;
        else if (!strcmp(optp,"--vasp")) format=VASP;
        else if (!strcmp(optp,"--chgcar"))
          {format=VASP; flags|=CHGCAR; flags|=CHDEN;}
        else if (!strcmp(optp,"--xyz")) format=XYZ;
        else if (!strcmp(optp,"--cml")) format=CML;
        else if (!strcmp(optp,"--fdf")) format=FDF;
        else if (!strcmp(optp,"--gnu")) flags|=GNUPLOT;
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
        else if (!strcmp(optp,"--qe")) format=QE;
        else if (!strcmp(optp,"--null")) format=CNULL;

        else if (!strcmp(optp,"--primitive")) spg_op|=CSPG_PRIM;
        else if (!strcmp(optp,"--refine")) spg_op|=CSPG_REF;
        else if (!strcmp(optp,"--int")) spg_op|=CSPG_INT;
        else if (!strcmp(optp,"--schoen")) spg_op|=CSPG_SCH;
        else if (!strcmp(optp,"--symmetry")) spg_op|=CSPG_SYM;
        else if (!strcmp(optp,"--point")) spg_op|=CSPG_PNT;
        else if (!strcmp(optp,"--list")) spg_op|=(CSPG_LST|CSPG_SYM);

        else if (!strcmp(optp,"--help")) help();
        else if (!strcmp(optp,"--refs")) refs();
        else if (!strcmp(optp,"--version")) version();
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
          flags|=BANDDEN;
          if(*(optp+1)=='='){
            elect.band_range=optp+2;
            while(*((++optp)+1));
          }
          break;
        case 'k':
          if(*(optp+1)=='='){
            elect.kpt_range=optp+2;
            while(*((++optp)+1));
          }
          else error_exit("malformed option -k");
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
          dipole_ctr=malloc(3*sizeof(double));
          if(*(optp+1)=='='){
            i=sscanf(optp+2,"%lf,%lf,%lf",dipole_ctr,dipole_ctr+1,
                     dipole_ctr+2);
            if ((i==0)||(i==EOF)) dipole_ctr[0]=dipole_ctr[1]=dipole_ctr[2]=0.5;
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
            if (sscanf(optp+2,"%d,%d,%d",i_grid,i_grid+1,i_grid+2)!=3)
              error_exit("malformed option -i=");
            while(*((++optp)+1));
          }
          else flags|=BANDIMAG;
          break;
        case 'L':
          flags|=LHS_FUDGE;
          break;
        case 'l':
          no_mp=1;
          break;
        case 'M':
	  gen_mp=1;
          no_mp=1;
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
          if(*(optp+1)=='='){
            m_abc=malloc(3*sizeof(int));
            if (!m_abc) error_exit("malloc error for three ints!");
            if (sscanf(optp+2,"%d,%d,%d",m_abc,m_abc+1,m_abc+2)!=3)
              error_exit("malformed option -m=");
            while(*((++optp)+1));
          }
          molecule=1;
          break;
        case 'n':
          no_sym++;
          break;
        case 'O':
          flags|=OCCUPANCIES;
          break;
        case 'p':
          flags|=BANDPHASE;
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
        case 'Q':
          sort_style=1;
          if ((*(optp+1)>='0')&&(*(optp+1)<='9')){
            optp++;
            sort_style=*optp-'0';
          }
          break;
        case 'r':
          flags|=BANDREAL;
          break;
        case 'R':
          flags|=RAW;
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
               &new_cell_rel[2][0],&new_cell_rel[2][1],&new_cell_rel[2][2])!=9)
               error_exit("malformed option -x=");
            while(*((++optp)+1));
            expand=2;
            break;
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
          else error_exit("malformed option -X=");
        case 'z':
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
          else flags|=ONETEP;
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
       if (!file1) file1=argv[opt];
       else if (!file2) file2=argv[opt];
       else{
         fprintf(stderr,"Unexpected argument %s\n%s -h for usage.\n",
                 argv[opt],argv[0]);
         exit(1);
       }
       opt++;
       optp=argv[opt];
       break;
     }
  }

  flags|=((preserve_c%3)<<PC_SHIFT);

  if (!file1) error_exit("no input file specified.");

  if ((!file2)&&(format==0)&&(spg_op&(CSPG_INT|CSPG_SCH|CSPG_PNT|CSPG_LST)))
    format=CNULL;

  if ((!file2)&&(isatty(fileno(stdout)))&&(format!=CNULL)&&(format!=FORT15))
    error_exit("refusing to output to a terminal");

  infile=fopen(file1,"rb");
  if(!infile){
    fprintf(stderr,"Error, unable to open %s for reading.\n",file1);
    exit(1);
  }

  if (file2){
    outfile=fopen(file2,"wb");
    if(!outfile){
      fprintf(stderr,"Error, unable to open %s for writing.\n",file2);
      exit(1);
    }
  }
  else
    outfile=stdout;

  i=strlen(file1);
  if((i>4)&&(!strcmp(file1+i-4,".pdb")))
    pdb_read(infile,&cell,&motif);
  else if ((i==7)&&(!strcmp(file1,"fort.34")))
    fort34_read(infile,&cell,&motif,&sym);
  else if ((i>2)&&(!strcmp(file1+i-2,"12")))
    crystal_read(infile,&cell,&motif,&sym);
  else if ((i>4)&&(!strcmp(file1+i-4,".res")))
    shelx_read(infile,&cell,&motif);
  else if ((i>4)&&(!strcmp(file1+i-4,".cif")))
    cif_read(infile,&cell,&motif,&sym);
  else if ((i>4)&&(!strcmp(file1+i-4,".cub")))
    cube_read(infile,&cell,&motif,&grid1);
  else if ((i>5)&&(!strcmp(file1+i-5,".cube")))
    cube_read(infile,&cell,&motif,&grid1);
  else if ((i>5)&&(!strcmp(file1+i-5,"_CUBE")))
    cube_read(infile,&cell,&motif,&grid1);
  else if ((i>4)&&(!strcmp(file1+i-4,".xsf")))
    xsf_read(infile,&cell,&motif,&grid1);
  else if ((i>=7)&&(!strcmp(file1+i-7,"CONTCAR")))
    vasp_read(infile,&cell,&motif,&grid1,&elect);
  else if ((i>=6)&&(!strcmp(file1+i-6,"POSCAR")))
    vasp_read(infile,&cell,&motif,&grid1,&elect);
  else if ((i>=6)&&(!strcmp(file1+i-6,"CHGCAR")))
    vasp_read(infile,&cell,&motif,&grid1,&elect);
  else if ((i>=3)&&(!strcmp(file1+i-3,"CHG")))
    vasp_read(infile,&cell,&motif,&grid1,&elect);
  else if ((i>=3)&&(!strcmp(file1+i-3,"DEN")))
    abinit_read(infile,&cell,&motif,&grid1);
  else if ((i>=3)&&(!strcmp(file1+i-3,"POT")))
    abinit_read(infile,&cell,&motif,&grid1);
  else if ((i>=3)&&(!strcmp(file1+i-3,".in")))
    abinit_in_read(infile,&cell,&motif,&kp,&sym,&elect);
  else if ((i>=7)&&(!strcmp(file1+i-7,"den_fmt")))
    denfmt_read(infile,&cell,&grid1);
  else if ((i>=4)&&(!strcmp(file1+i-4,".fdf")))
    fdf_read(infile,&cell,&motif,&kp,&elect);
  /* If not a Castep ending, check for a VASP beginning */
  else if ((!(((i>5)&&(!strcmp(file1+i-5,".cell")))||
              ((i>6)&&(!strcmp(file1+i-6,".check")))||
              ((i>4)&&(!strcmp(file1+i-4,".dat")))||
              ((i>11)&&(!strcmp(file1+i-11,".castep_bin")))||
              ((i>7)&&(!strcmp(file1+i-7,".chdiff")))))&&
           ((!strncmp(file1,"CONTCAR",7))||(!strncmp(file1,"POSCAR",6))||
            (!strncmp(file1,"CHGCAR",6))||(!strncmp(file1,"CHG",3))))
    vasp_read(infile,&cell,&motif,&grid1,&elect);
  else{
    i=fgetc(infile);
    rewind(infile);
    if ((i==0)||(i==30)||(i==10)) /* possible first byte of .check or */
                                  /* .castep_bin */
      check_read(infile,&cell,&motif,&kp,&sym,&grid1,&elect,i_grid);
    else
      cell_read(infile,&cell,&motif,&kp,&sym);
  }
  fclose(infile);

  if (!cell.basis){
    fprintf(stderr,
            "Error: no basis set found. Was format correctly detected?\n");
    exit(1);
  }
  
  if (cell.vol<0) cell.vol=fabs(cell.vol);

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
    chdiff_read(infile,&grid1);
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
    esp_read(infile,&grid1,&elect);
    fclose(infile);
    free(cdfile);
  }

  /* Check that we have most of what we need */

  if (!motif.atoms)
    if ((format!=DENFMT)&&(format!=CUBE)&&(format!=DX)&&
        (format!=FBIN)&&(format!=CNULL))
      error_exit("no atoms found!");

  ionic_charge=0;
  if (motif.atoms)
    for(i=0;i<motif.n;i++) ionic_charge+=motif.atoms[i].chg;
  
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
    for(i=0;i<=2;i++)
      fprintf(stderr,"% 11.7f % 11.7f % 11.7f   modulus % 11.7f\n",
              cell.basis[i][0],cell.basis[i][1],cell.basis[i][2],
              sqrt(cell.basis[i][0]*cell.basis[i][0]+
                   cell.basis[i][1]*cell.basis[i][1]+
                   cell.basis[i][2]*cell.basis[i][2]));
  }
  if (debug>2){
    fprintf(stderr,"Reciprocal basis set\n");
    for(i=0;i<=2;i++)
      fprintf(stderr,"% 11.7f % 11.7f % 11.7f   modulus % 11.7f\n",
              cell.recip[i][0],cell.recip[i][1],cell.recip[i][2],
              sqrt(cell.recip[i][0]*cell.recip[i][0]+
                   cell.recip[i][1]*cell.recip[i][1]+
                   cell.recip[i][2]*cell.recip[i][2]));
  }
  if (debug){
    fprintf(stderr,"Cell volume %f\n",cell.vol);
    fprintf(stderr,"natoms      %d\n",motif.n);
    if (ionic_charge) fprintf(stderr,"Total ionic charge %f\n",ionic_charge);

    if ((kp.mp)&&(kp.mp->grid[0]>0)){
      fprintf(stderr,"kpoint MP grid %d %d %d\n",kp.mp->grid[0],
              kp.mp->grid[1],kp.mp->grid[2]);
      fprintf(stderr,"        offset %lf %lf %lf\n",kp.mp->disp[0],
              kp.mp->disp[1],kp.mp->disp[2]);
    }

  }

  if ((dipole_ctr)&&(grid1.data)) dipole(&cell,&motif,&grid1,dipole_ctr);

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
    if (debug){
      double min,max,sum;
      min=1e20;
      max=-1e20;
      sum=0;

      fprintf(stderr,"Found 3D data for %s\n",gptr->name);
      for(i=0;i<gptr->size[0]*gptr->size[1]*gptr->size[2];i++){
	sum=sum+gptr->data[i];
	if (gptr->data[i]>max) max=gptr->data[i];
	if (gptr->data[i]<min) min=gptr->data[i];
      }
      fprintf(stderr,"  min=%g  max=%g  sum=%g  int=%g\n",min,max,
	      sum,sum*cell.vol/(gptr->size[0]*gptr->size[1]*gptr->size[2]));
      fprintf(stderr,
	      "  (integral is e per cell for charge and spin densities)\n");
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
        es_pot(&cell,&motif,&grid1,musq);
        calc_esp=0;
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

  if ((calc_esp)&&(grid1.data)) es_pot(&cell,&motif,&grid1,musq);
    
  if (pt_spec){
    lscan(&pt_spec,&motif,zpt);
    gptr=&grid1;
    while(gptr&&(gptr->data)){
      printf("%s at (%f,%f,%f): %e\n",gptr->name,zpt[0],zpt[1],zpt[2],
             interpolate0d(gptr,zpt));
      gptr=gptr->next;
    }
  }

  if (gen_mp){
    if (kp.mp) mp_gen(&kp,&cell);
    else fprintf(stderr,"Warning: ignoring -M as no MP parameters found\n");
  }

  if (failure) fstar(&kp, &cell, &sym);

  if (molecule) molecule_fix(m_abc,&cell,&motif,&grid1);


  if (half_shift) xplor_fudge(&cell,&motif,&grid1);

  if (spg_op) cspq_op(&cell,&motif,&sym,spg_op,tolmin);

  if (prim) {
    double nb1[3][3];
    primitive(&cell,&motif,nb1);
    if (debug>2){
      fprintf(stderr,"Primitive cell is:\n");
      for(i=0;i<3;i++){
        fprintf(stderr,"(%f,%f,%f)\n",nb1[i][0],nb1[i][1],nb1[i][2]);
      }
    }
    super(&cell,&motif,nb1,&kp,&sym,&grid1,0);
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
    } /* select(expand) ... */
  }  /* if (expand) ... */

  if (rotate) cart2abc(&cell,&motif,abc,&grid1,1);
  if ((no_mp)&&(kp.mp)) {free(kp.mp);kp.mp=NULL;}
  if (sort_style) sort_atoms(&motif,sort_style);
  if (no_sym){
    sym.n=0;
    if (sym.ops) free(sym.ops);
  }
  if (no_sym==2){
    kp.n=0;
    if (kp.mp){
      free(kp.mp);
      kp.mp=NULL;
    }
  }

  if (line_spec){
    line_write(outfile,&cell,&motif,&grid1,line_spec);
    exit(0);
  }

  switch(format){
    case XSF:
      xsf_write(outfile,&cell,&motif,molecule,&grid1);
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
      xyz_write(outfile,&cell,&motif);
      break;
    case CML:
      cml_write(outfile,&cell,&motif);
      break;
    case FDF:
      fdf_write(outfile,file2,&cell,&motif,&grid1,&elect);
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
    case CNULL:
      break;

    default:
      fprintf(stderr,"This cannot happen. Sorry\n");
  }

  return(0);
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
         "-C           find compact set of cell vectors\n"
         "-d           read also a .chdiff file, constructing its name from\n"
         "               the .cell or .check file given\n"
         "-D=[x,y,z]   calculate dipole moment about fractional co-ordinates\n"
         "               x,y,z (default .5,.5,.5) assuming -c also given\n"
         "-e           read also a .cst_esp file, constructing its name from\n"
         "               the .cell or .check file given\n"
         "-e=eps       set tolerance for supercell operations. Default=%g\n"
         "-e=min-max   set tolerance range for spglib operations\n"
         "-E[=[-][mu]] calculate electrostatic potential. Default mu=%g\n"
         "-f           find first failure star of k-point mesh\n"
         "-H           shift atoms by half a grid cell\n",tol,C2X_MU);
  printf("-i           output imaginary part of band\n"
         "-i=n1,n2,n3  Fourier interpolate 3D grids to specified grid\n"
         "-I[=range]   report inversion symmetries of given bands\n"
         "-k=range     include given k-points (default 1) for bands\n"
         "-l           list k-points in .cell output, not MP parameters\n"
         "-L           produce (incorrect) left-handed abc output\n");
  printf("-m[=a,b,c]   assume input is molecule, not crystal, and move by\n"
         "               given nos of grid cells, or move automatically\n"
         "-M           generate Monkhurst-Pack k-point set from .cell values\n"
         "-M=nx,ny,nz  generate regular k-point mesh including origin\n"
         "-M=nx,ny,nx,dx,dy,dz\n"
         "               ditto displaced by fraction of mesh cell\n"
         "-n           discard symmetry information\n"
         "               give twice to discard k-points too\n"
         "-O           print band occupancies and eigenvalues\n");
  printf("-P           find primitive cell\n"
         "-P=X:Y:npts  output 1D data with npts points along line\n"
         "               X and Y either (x,y,z) as fractional co-ords\n"
         "               or atomic symbol followed by atom number\n"
         "               e.g. (0.5,0.5,0.5):C3:20 -- centre to 3rd C atom, "
         "20 points\n"
         "-Q[n]        quicksort atoms in descending atomic order "
         "(n=1 or absent)\n"
         "                                ascending atomic order  (n=2)\n"
         "-R           don't rescale densities or bands\n"
         "-s           include spin densities\n"
         "-S=range     include given spins/spinors (0 or 1) for bands\n");
  printf("-t=(x1,y1,z1)(x2,y2,z2)[(x3,y3,z3)]\n"
         "             rotate coords so 1st vector becomes 2nd, using third\n"
         "             vector as axis, else axis is perpendicular to others\n"
         "-T=(x1,y1,z1)(x2,y2,z2)[(x3,y3,z3)]\n"
         "             ditto, but first vec expressed in absolute coords\n"
         "-u           write .cell files and 1D axes in Bohr (atomic Units)\n"
         "             scale densities on writing .cube files from A^-3 to "
         "Bohr^-3\n"
         "-U           scale densities on reading .cube files from "
         "Bohr^-3 to A^-3\n"
         "-v           be verbose (may be repeated)\n"
         "-w           weight band by occupancy\n"
         "-W           weight band by kpt weight and occupancy\n");
  printf("-x=(x1,y1,z1)(x2,y2,z2)(x3,y3,z3)\n"
         "             re-express in new basis given in terms of old\n"
         "-X=(x1,y1,z1)(x2,y2,z2)(x3,y3,z3)\n"
         "             re-express in new basis given in absolute terms\n"
         "-z=X         print volumetric data at point\n"
	 "               see -P= for point specification\n"
         "-1           assume input .cell file follows Onetep conventions\n"
         "-3           if swapping axes to convert lhs to rhs, third axis\n"
         "               is special & first two swapped. Else first special.\n"
         "-15          use high precision in output.\n"
         "\n");
  printf("FORMAT is one of: xsf       XCrySDen (default)\n"
	 "                  abinit    Abinit .in\n"
         "                  cell      CASTEP .cell, cartesian and fractional\n"
         "                  cell_abc                abc and fractional\n"
         "                  cell_abs                cartesian and absolute\n"
         "                  cell_abc_abs            abc and absolute\n"
         "                  cif / mmcif  Basic maybe CIF compatible file\n"
         "                  chgcar    VASP CHGCAR (implies -c)\n"
         "                  cml       Chemical Markup Language\n"
         "                  cube      Gaussian cube\n"
         "                  denfmt    CASTEP formatted density\n"
         "                  dx        OpenDX\n"
         "                  fdf       Flexible Data Format (Siesta) (beta)\n"
         "                  gnu       Gnuplot (with -P=)\n"
         "                  null      Discard output\n"
         "                  one       Onetep .dat, very similar to .cell\n"
         "                              also one_abc, one_abs, one_abc_abs\n"
         "                  pdb       PDB\n"
         "                  pdbn      PDB with atoms numbered\n"
         "                  py        python dictionary\n"
         "                  pya       python ASE Atoms\n"
         "                  qe        Quantum Espresso\n"
         "                  refs      List BibTeX references for this code "
         "and exit\n"
         "                  shelx     SHELX97\n"
         "                  vasp      VASP output (POSCAR/CHG)\n"
         "                  xplor     Xplor\n"
         "                  xyz       XYZ\n\n");
#ifdef SPGLIB
  printf("OPERATION is used to call spglib and is one of:\n"
         "   primitive      call spg_find_primitive()\n"
         "   refine         call spg_refine_cell()\n"
         "   int            call spg_get_dataset()"
         " and report international symbol\n"
         "   schoen         call spg_get_schoenflies()\n"
         "   symmetry       call spg_get_dataset() and keep symmetry ops\n"
         "   list           call spg_get_dataset() and list symmetry ops\n"
         "   point          call spg_get_dataset() followed by "
         "spg_get_pointgroup()\n\n");
#endif
  printf("range specifies band numbers as \"a,b-c,d\"\n"
         "-b and -B are mutually exclusive. Only one of the x and t"
         " options may be given\n\n");
  printf("Input files ending .cif are assumed to be in a cif format,\n"
	 "            ending .cub, .cube or _CUBE are assumed to be in "
            "cube format,\n"
	 "            ending .fdf are assumed to be in Siesta format,\n"
	 "            ending .in are assumed to be in either Abinit or\n"
         "              Quantum Espresso format,\n"
         "            ending .pdb are assumed to be in pdb format,\n"
         "            ending .res are assumed to be in shelx97 format,\n"
         "            ending .xsf are assumed to be in xsf format,\n"
	 "            ending POT or DEN are assumed to be in "
	 "Abinit binary format,\n");
  
  printf("            ending or beginning CHG, CHGCAR, POSCAR or CONTCAR "
         "are assumed \n"
         "              to be in VASP 5.x format,\n"
	 "            ending 12 are assumed to be in Crystal format,\n"
	 "            called fort.34 are assumed to be in Crystal binary "
	 "format.\n\n"
         "Otherwise "
         "automatic detection of .cell or .check input. Compatible with\n"
         ".check files from CASTEP 3.0 to 19.1 (and perhaps beyond).\n\n"
         "Version " C2XSF_VER ", (c) MJ Rutter 2007 - 2019"
         " licenced under the GPL v3.\n\n");
  printf("If useful to a published paper, please consider citing using the\n"
         "references shown with the --refs argument.\n\n");
  printf("This version is "
#ifndef SPGLIB
         "not "
#endif
         "linked against spglib.\n\n");
  exit(0);
}

void refs(void){
  printf("@article{c2x,\n"
         "  title = {C2x: a tool for visualisation and input preparation for "
      "{C}astep and other electronic structure codes},\n"
         "  journal = {Computer Physics Communications},\n"
         "  author = {Rutter, M. J.},\n"
         "  year = {2018},\n"
         "  volume = {225C},\n"
         "  pages = {152--157},\n"
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
#endif
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
}
