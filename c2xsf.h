#define C2XSF_VER "2.40c"

/* Global variables for system description */

/* Note Castep 8, 16 and 17 use CODATA 2010
 *                18 to at least 21 use CODATA 2014
 *
 * See https://physics.nist.gov/cuu/Constants/index.html
 */

#define BOHR 0.52917721067   /* CODATA 2014 */
          /* 0.52917721092   is CODATA 2010 */
          /* 0.5291772101121114 is c2x <=2.20a */

#define H_eV 27.21138602  /* CODATA 2014 */
          /* 27.21138505  is CODATA 2010 */
          /* 27.21138342902473 is c2x <=2.20a */

/* Atomic unit of time, in picoseconds */
#define H_ps 2.418884326585e-5 /* CODATA 2014 */

/* epsilon_0 is 8.85418782e-12 Fm^-1 or C V^-1 m-1
 * We want A not m, and e not C, so
 * 8.85418782e-12/1.60217662e-19/1e-10
 * Strangely the answer is generally expressed as its reciprocal */

#define EPS0 (1/180.95128)


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
#define BXSF 28
#define FDF_BP 29
#define ELK 30
#define NPY 31

/* Constants for data combining operations */
/* Exclusive operations */
#define C2X_ADD 1
#define C2X_DIFF 2
#define C2X_MASK 3
/* Non-exclusive operations */
#define C2X_MERGE 256
#define C2X_EXC_MASK 0xff

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
#define GCOEFF 131072
#define HIPREC 262144
#define LHS_FUDGE 524288
#define OCC_WEIGHT 1048576
#define K_WEIGHT 2097152
#define OCCUPANCIES 4194304
#define REPLACE_RHO 8388608
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

#define BANDREAD (BANDS+BANDPARITY+GCOEFF)

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
#define CSPG_SNAP_TR 2048
#define CSPG_STD_IDEAL 4096

struct dct {char *key; void *value; struct dct *next;};
struct cmt {char *txt; struct cmt *next;};
/* Currently atom.labels leak and should not be freed, as pointers for
   two different atoms may point to the same location */
struct atom
   {unsigned int atno; double abs[3]; double frac[3]; double force[3];
     double v[3]; double wt; double spin; double chg; double site_chg;
     char *label;};
/* This structure will probably be extended... */
struct species {unsigned int atno;};
/* grid.next == NULL if grid unused */
/* grids storage order is size[0]=ngx, size[1]=ngy, size[2]=ngz,
     data[x*ngy*ngz+y*ngz+z] */
struct grid {char *name; int size[3]; double *origin_abs; double *data;
             struct grid *next;};
/* See ksym.c:mp_gen() for precise definition of mp_grid */
struct mp_grid {int grid[3]; double disp[3];};
struct vector {double v[3]; double mod2;};
/* We store symmetry matrix and translations in absolute co-ords
   And the matrix is the transpose of Castep's convention, or
   identical to SPGlib's convention */
struct sym_op {double mat[3][3]; double *tr;};
struct unit_cell {double (*basis)[3]; double recip[3][3]; double vol;
  double (*stress)[3];};
struct contents {int n; int forces; int velocities; struct atom *atoms;
  char *title; struct cmt *comment; char *block_species; char *species_misc;
  struct dct *dict; int nspec; struct species *spec;};
struct kpts {int n; struct atom *kpts; struct mp_grid *mp; double *spacing;
             int bs_n; struct atom *bs_kpts; struct mp_grid *bs_mp;
             double *bs_spacing;
             int path_n; struct atom *path; double *path_spacing;};
struct symmetry {int n; double *tol; int *gen; struct sym_op *ops;};
struct es {int nspins; int nspinors; char *spin_method; double cut_off;
  double etol; char *band_range; char *kpt_range; char *spin_range;
  char *dip_corr; char *dip_corr_dir; double *dip_ctr; double *charge;
  double *energy; double *e_fermi; int nbands; int nbspins; double nel;
  double nup_minus_down; double *occ; double *eval; int max_nplwv;
  double *path_eval; double *path_occ; struct kpts *path_kpt;
  int path_nbands;};

struct time_series {int nsteps; int nc; struct unit_cell *cells;
  int nm; struct contents *m; int nen; double *energies;
  int nenth; double *enthalpies;};

struct infiles { FILE* f; struct infiles *next; struct infiles *last;
  char *name; int line; int include; struct infiles *ret; int count;};

void *dict_get(struct dct *dict, char *key);
void dict_add(struct dct *dict, char *key, void *value);
void dict_strcat(struct dct *dict, char *key, char *value);

#ifdef __GNUC__
void error_exit(char* msg) __attribute__((__noreturn__));
#else
void error_exit(char* msg);
#endif

void real2rec(struct unit_cell *cell);
void addfrac(struct atom *a,int natoms, double recip[3][3]);
void addabs(struct atom *a,int natoms, double basis[3][3]);
void addspec(struct contents *m);
void reduce_cell(struct atom *a,int natoms, double basis[3][3]);
void reduce_cell_tol(struct atom *a,int natoms, double basis[3][3], double eps);
void abc2cart(double *abc, struct unit_cell *cell);
void basis2abc(double basis[3][3], double abc[6]);
void cart2abc(struct unit_cell *c, struct contents *m, double *abc, 
              struct grid *g, int fix);
void cart2abc_sym(struct unit_cell *c, struct contents *m, double *abc, 
                  struct grid *g, int fix, struct symmetry *s);
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
void print_bandwidths(struct es *elect, struct kpts *kp);
void print_energy(double e);
double calc_efermi(struct es *elect, struct kpts *kpt, double nel);
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
void pad_recip(double *o, int fft[3], double **nptr, int nfft[3]);

void band_process(double *dptr, int fft[3], int *pwgrid, int npw, int gamma,
		  struct unit_cell *c,
		  struct grid **gp, struct es *elect, struct kpts *kp,
		  struct contents *m, int ikpt, int ispinor, int isppol,
		  int nb, int *i_grid);
double *band2grid(double *dptr, int fft[3], int *pwgrid, int npw, int gamma);
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
int is_identity(double m[3][3]);
void make_rhs(struct unit_cell *c, struct contents *m, double *abc,
              struct grid *gptr);
void cell_check(struct unit_cell *c, struct contents *m);
void include_file(struct infiles **file, char *dir, char *ptr);
void sym2ksym(struct symmetry *rs, struct symmetry *ks);
char *strrsubs(char *str, char *old, char *new);
int cspg_op(struct unit_cell *c, struct contents *m, struct symmetry *s,
            struct kpts *k, int op, double tolmin);
void xv_read(FILE* infile, struct unit_cell *c, struct contents *m);
void fdf_read(FILE* in, struct unit_cell *c,
              struct contents *m, struct kpts *kp, struct es *e);
void siesta_kp_read(FILE* infile, struct kpts *k);
int sscanfmsn(char *buffer, char **str, int *n);
struct grid *grid_new(struct grid *gptr);
void init_grid(struct grid *gptr);
void init_cell(struct unit_cell *c);
void init_motif(struct contents *m);
void init_elect(struct es *e);
void init_kpts(struct kpts *k);
void init_sym(struct symmetry *s);
void init_tseries(struct time_series *ts);
void interpolate3d(struct grid *old_grid, struct grid *new_grid);

void reverse4n(int *data,int n);
void reverse8n(double *data,int n);
int self_little_endian(void);

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
