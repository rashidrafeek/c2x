#define C2XSF_VER "2.32"

/* Global variables for system description */

/* Note Castep 8, 16 and 17 use CODATA 2010
 *                18 and 19 use CODATA 2014
 *
 * See https://physics.nist.gov/cuu/Constants/index.html
 */

#define BOHR 0.52917721067   /* CODATA 2014 */
          /* 0.52917721092   is CODATA 2010 */
          /* 0.5291772101121114 is c2x <=2.20a */

#define H_eV 27.21138602  /* CODATA 2014 */
          /* 27.21138505  is CODATA 2010 */
          /* 27.21138342902473 is c2x <=2.20a */

/* epsilon_0 is 8.8541878e-12 Fm^-1 or C V^-1 m-1
 * We want A not m, and e not C, so
 * 8.8541878e-12/1.6021892e-19/1e10
 * Strangely the answer is generally expressed as its reciprocal */

#define EPS0 (1/180.952701)


#ifndef M_PI
#define M_PI 3.14159265358979324
#endif

/* constants for the file formats which can be written */
#define XSF 0
#define CUBE 1
#define XPLOR 2
#define PDB 3
#define CELL 4
#define CELL_ABC 5
#define CELL_ABS 6
#define CELL_ABC_ABS 7
#define DX 8
#define VASP 9
#define XYZ 10
#define CML 11
#define FDF 12
#define CNULL 13
#define SHELX 14
#define FBIN 15
#define FORT15 16
#define CIF 17
#define MMCIF 18
#define PYA 19
#define PY 20
#define DENFMT 21
#define ABINIT 22
#define QE 23
#define CASTEP_BANDS 24
#define CASTEP_GEOM 25
#define XV 26
#define CCP4 27

/* flags for reading and output */
#define CHDEN 1
#define SPINDEN 2
#define BANDS 4
#define BANDDEN 8
#define BANDPHASE 16
#define BANDREAL 32
#define BANDIMAG 64
#define ACCUMULATE 128
#define RAW 256
#define PDBN 512
#define CHDIFF 1024
#define SHELX_AIRSS 2048
#define LITTLE_OUT 4096
#define BANDPARITY 8192
#define FRAC 16384
#define CST_ESP 32768
#define AU 65536
#define HIPREC 262144
#define LHS_FUDGE 524288
#define OCC_WEIGHT 1048576
#define K_WEIGHT 2097152
#define OCCUPANCIES 4194304
#define ONETEP 8388608
/* Alternate output format: changes CELL to ONETEP,
 *                          changes CUBE to MO CUBE
 */
#define ALT_OUT 16777216
#define CHGCAR 33554432
#define DE_AU 67108864

/* Last valid number 2^(PC_SHIFT-1)=134217728 */

/* Other flags */

#define PC_SHIFT 28
/* 1<<28 + 1<<29 = 3<< PC_SHIFT = 805306368 */
#define PRESERVE_C 805306368

/* Flags for SPGLIB ops */

#define CSPG_PRIM 1
#define CSPG_REF 2
#define CSPG_INT 4
#define CSPG_SCH 8
#define CSPG_SYM 16
#define CSPG_PNT 32
#define CSPG_LST 64
#define CSPG_STD 128
#define CSPG_SNAP 256
#define CSPG_PRIM_NR 512
#define CSPG_NO_SORT 1024

struct dct {char *key; void *value; struct dct *next;};
struct cmt {char *txt; struct cmt *next;};
/* Currently atom.labels leak and should not be freed, as pointers for
   two different atoms may poin to the same location */
struct atom
   {unsigned int atno; double abs[3]; double frac[3]; double force[3];
     double wt; double spin; double chg; double site_chg; char *label;};
/* grid.next == NULL if grid unused */
/* grids storage order is size[0]=ngx, size[1]=ngy, size[2]=ngz,
     data[x*ngy*ngz+y*ngz+z] */
struct grid {char *name; int size[3]; double *data;
             struct grid *next;};
/* See ksym.c:mp_gen() for precise definition of mp_grid */
struct mp_grid {int grid[3]; double disp[3];};
struct vector {double v[3]; double mod2;};
/* We store symmetry matrix and translations in absolute co-ords */
struct sym_op {double mat[3][3]; double *tr;};
struct unit_cell {double (*basis)[3]; double recip[3][3]; double vol;
  double (*stress)[3];};
struct contents {int n; int forces; struct atom *atoms; char *title;
  struct cmt *comment; char *block_species; char *species_misc;
  struct dct *dict;};
struct kpts {int n; struct atom *kpts; struct mp_grid *mp; double *spacing;};
struct symmetry {int n; double *tol; int *gen; struct sym_op *ops;};
struct es {int nspins; int nspinors; char *spin_method; double cut_off;
  double etol; char *band_range; char *kpt_range; char *spin_range;
  char *dip_corr; char *dip_corr_dir; double *dip_ctr; double *charge;
  double *energy; double *e_fermi; int nbands; int nbspins; double nel;
  double nup_minus_down; double *occ; double *eval;};

struct time_series {int nsteps; int nc; struct unit_cell *cells;
  int nm; struct contents *m; int nen; double *energies;
  int nenth; double *enthalpies;};

void *dict_get(struct dct *dict, char *key);
void dict_add(struct dct *dict, char *key, void *value);
void dict_strcat(struct dct *dict, char *key, char *value);

void error_exit(char* msg);
void real2rec(struct unit_cell *cell);
void addfrac(struct atom *a,int natoms, double recip[3][3]);
void addabs(struct atom *a,int natoms, double basis[3][3]);
void reduce_cell(struct atom *a,int natoms, double basis[3][3]);
void reduce_cell_tol(struct atom *a,int natoms, double basis[3][3], double eps);
void abc2cart(double *abc, struct unit_cell *cell);
void basis2abc(double basis[3][3], double abc[6]);
void cart2abc(struct unit_cell *c, struct contents *m, double *abc, 
              struct grid *g, int fix);
void init_atoms(struct atom *a, int n);
void print_globals(int level);
void ident_sym(struct sym_op *s, struct unit_cell *c, FILE *out);
void equiv_sym(struct sym_op *s, struct unit_cell *c, FILE *out);
void mat_f2a(double m1[3][3], double m2[3][3], double basis[3][3],
             double recip[3][3]);
void mat_a2f(double m1[3][3], double m2[3][3], double basis[3][3],
             double recip[3][3]);
void equiv_read(struct sym_op *s, struct unit_cell *c, char *line);
int atom_in_list(struct atom *b, struct atom *a, int n, double basis[3][3]);
void sym_atom(struct atom *a, struct atom *b, struct sym_op *s,
              double recip[3][3]);
void sort_atoms(struct contents *mtf, int sort_style);
void add_cmt(struct cmt *comment, char *txt);
void sym_expand(struct unit_cell *c, struct contents *m, struct symmetry *s);
int tokenmatch(char **s1, const char *s2);
void inv_parity(double *d,int fft[3], int band, double kpt[3]);
void dipole_calc(struct unit_cell *c, struct contents *m,
                 struct grid *g, double *dipole_ctr, double *dpole);
void print_cell(struct unit_cell *c, struct contents *m);
void print_occ(struct es *elect, struct kpts *kp);
void print_elect(struct es *elect);
void print_basis(double basis[3][3]);
void old_in_new(double old_basis[3][3],double new_recip[3][3]);

unsigned int atsym2no(char* sym);
char* atno2sym(unsigned no);

void cspg_hall2sym(int hall, struct unit_cell *c, struct symmetry *s);
int spgr_is_double(int spgr);

int ascan(char *in, double *result);
int single_scan(char *buff, double *x, int *n);
int multi_scan(char *buff, double *x, int rep, int *n);

void mp_gen(struct kpts *ks, struct unit_cell *c);

int inrange(int x, char *range);
void fft3d(double *c, int *ngptar, int dir);
void band2real(double *psi, double *out, int nfft[3], double kpoint[3]);
void band_store(struct grid **g, double *dptr, double occ, double wkpt,
                int nspr, int ns, int k, int b, struct es *elect,
                struct contents *m, int fft[3]);
void pad_recip(double *o, int fft[3], double **nptr, int nfft[3]);


extern int igr2hall[];

int super(struct unit_cell *c, struct contents *m,
           double new_basis[3][3], struct kpts *k, struct symmetry *s,
           struct grid *gptr, int rhs);
void simple_super(struct unit_cell *c, struct contents *m,
           int expand[3], struct kpts *kp, struct symmetry *s,
           struct grid *gptr);

double dist(double a,double b);
double atom_dist(struct atom *a, struct atom *b, double basis[3][3]);
double vmod2(double v[3]);
void cell_check(struct unit_cell *c, struct contents *m);

#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif

#define aeq(a,b) (fabs((a)-(b))<((fabs(a)+0.5)*tol))

#ifndef C2X_MAIN
extern int debug,flags;
extern double tol;
#endif
