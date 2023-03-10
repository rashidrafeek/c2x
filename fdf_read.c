/* Read some useful data from a SIESTA .fdf file
 * 
 */

/* 
 * This version does not parse quotation marks properly.
 */

/* Note that Siesta ignores any of the three characters -_. appearing
 * in its keywords.
 */


/* Note too Siesta's strange reading of include files
 *
 *  include foo.dat  -- entirely conventional, insert foo.dat here
 *
 * label < foo.dat -- parse foo.dat as fdf file, update only label,
 * which may be a block...
 *
 * %block label < foo.dat
 *
 * parses as
 *
 * %block label
 * [contents of foo.dat]
 * [implicit %endblock label]
 *
 */

/* Copyright (c) 2018-2021 MJ Rutter 
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
#include<string.h>
#include<ctype.h>
#include<errno.h>

#include "c2xsf.h"

#define LINE_SIZE 200

#define PARSE_ERROR  fprintf(stderr,"Parse error on line %d %s\n",files->line,(files->name)?files->name:"")

static int fdfreadline(char *buffer, int len, char *dir,
                       struct infiles **filesp);
static int fdfreadlength(char *buff, double *x);
static int fdfreadenergy(char *buff, double *x);
static int fdfstrcasecmp(char **s1_in, char *s2);
static int fdftokenmatch(char **s1_in, char *s2);
static int fdfscanfile(char *filename, char *dirname, char *buffer,
                       char *token, struct infiles **infile);
static void fdf_unwind(struct infiles **infile);
double siesta_pot_charge(char *prefix, char *label); /* Used by xv_read.c */
void dequote(char **str);

static char *ptr;

void fdf_read(FILE* in, struct unit_cell *c, struct contents *m,
              struct kpts *kp, struct es *e){
  int i,j,k,n,*atomspecies,mp_grid[3][3],okay;
  int nspec,*species_to_atno;
  double scale,shift[3],mp_off[3];
  char *ptr2,*format;
  static char buffer[LINE_SIZE+1];
  double *vec,*abc,*atomcoords,*dmetol,latticeconstant,*atomspins;
  double *species_to_charge;
  char *siesta_misc[]={"ElectronicTemperature","MeshCutoff",
                       "MaxSCFIterations","MinSCFIterations",
                       "OccupationFunction","SaveDeltaRho",
                       "SaveElectrostaticPotential",
                       "SaveRho","SaveRhoXC",
                       "SaveTotalCharge","SaveTotalPotential",
                       "SlabDipoleCorrection","SolutionMethod",
                       "WriteDM","XCauthors","XCfunctional",
                       ""};
  char *spaces="                              ";  /* 30 spaces */
  char *prefix;
  char *dir;
  struct infiles *files;
 
  dir=dict_get(m->dict,"in_dir");
  if (dir)
    prefix=dir;
  else
    prefix="";
  
  files=malloc(sizeof(struct infiles));
  if (!files) error_exit("Malloc error for files in fdf_read");
  files->f=in;
  files->next=files->last=NULL;
  files->name=NULL;
  files->include=0;
  files->line=0;
  files->ret=NULL;
  files->count=0;
  
  m->n=0;
  nspec=0;
  vec=abc=NULL;
  latticeconstant=0;
  format=NULL;
  species_to_atno=NULL;
  species_to_charge=NULL;
  atomspecies=NULL;
  shift[0]=shift[1]=shift[2]=0;
  dmetol=NULL;
  atomspins=NULL;
  atomcoords=NULL;
  
  if (debug>2) fprintf(stderr,"fdf_read called\n");

  if (!(c->basis=malloc(72)))
    error_exit("Malloc error in fdfread for c->basis");

  while(fdfreadline(buffer,LINE_SIZE,dir,&files)){
    ptr=buffer;
    
    if (!fdftokenmatch(&ptr,"SystemName")){
      m->title=malloc(strlen(ptr)+1);
      strcpy(m->title,ptr);
      continue;
    }

    if (!fdftokenmatch(&ptr,"NumberOfAtoms")){
      if (sscanf(ptr,"%d",&(m->n))!=1)
        error_exit("Error parsing number of atoms");
      continue;
    }

    if (!fdftokenmatch(&ptr,"NumberOfSpecies")){
      if (sscanf(ptr,"%d",&nspec)!=1)
        error_exit("Error parsing number of species");
      continue;
    }

    if (!fdftokenmatch(&ptr,"LatticeConstant")){
      if (!fdfreadlength(ptr,&latticeconstant))
        error_exit("Error parsing lattice constant");
      continue;
    }

    if (!fdftokenmatch(&ptr,"AtomicCoordinatesFormat")){
      while (*ptr==' ') ptr++;
      ptr2=ptr;
      while (*ptr2&&(*ptr2!=' ')) ptr2++;
      format=malloc(ptr2-ptr+1);
      strncpy(format,ptr,(ptr2-ptr)+1);
      continue;
    }

    if (!fdftokenmatch(&ptr,"AtomicCoordinatesOrigin")){
      if (sscanf(ptr,"%lf %lf %lf",shift,shift+1,shift+2)!=3)
        error_exit("Error parsing atomic coordinates origin");
      continue;
    }

    if (!fdftokenmatch(&ptr,"DMEnergyTolerance")){
      dmetol=malloc(sizeof(double));
      if (!fdfreadenergy(ptr,dmetol)){
        fprintf(stderr,"Error parsing DM.Energy.Tolerance\n");
        free(dmetol);
        dmetol=NULL;
      }
      continue;
    }
        
    if (!fdftokenmatch(&ptr,"NetCharge")){
      e->charge=malloc(sizeof(double));
      if (sscanf(ptr,"%lf",e->charge)!=1){
        fprintf(stderr,"Warning: error parsing NetCharge\n");
        free(e->charge);
        e->charge=NULL;
      }
      continue;
    }

    if (!fdftokenmatch(&ptr,"SystemLabel")){
      ptr2=malloc(strlen(ptr)+1);
      strcpy(ptr2,ptr);
      dict_add(m->dict,"Siesta_SystemLabel",ptr2);
    }

    for(i=0;siesta_misc[i][0];i++){
      if (!fdftokenmatch(&ptr,siesta_misc[i])){
        dict_strcat(m->dict,"Siesta_misc",siesta_misc[i]);
        while ((*ptr)&&isspace(*ptr)) ptr++;
        if (strlen(ptr)){
          if (strlen(siesta_misc[i])<28) /* 28 < length spaces array */
            dict_strcat(m->dict,"Siesta_misc",spaces+strlen(siesta_misc[i]));
          else
            dict_strcat(m->dict,"Siesta_misc","  ");
          dict_strcat(m->dict,"Siesta_misc",ptr);
        }
        dict_strcat(m->dict,"Siesta_misc","\n");
      }
    }
    

/* What remains ought to be a block or a keyword.  */
    if (strncasecmp(ptr,"%block",6)) continue;
    ptr+=6;
/* Eat leading spaces */
    while(isspace(*ptr)) ptr++;
/* Kill trailing spaces and things after # */
    ptr2=ptr;
    while((*ptr2)&&(*ptr2!=' ')&&(*ptr2!='#')&&(*ptr2!='!')&&(*ptr2!=';'))
      ptr2++;
    *ptr2=0; /* Either it was a null, or it should be one... */

    if (debug>2) fprintf(stderr,"Found a %s block\n",ptr);

    if (!fdfstrcasecmp(&ptr,"latticevectors")){
      fdfreadline(buffer,LINE_SIZE,dir,&files);
      vec=malloc(9*sizeof(double));
      if (!vec) error_exit("malloc error for vectors");
      for(i=0;i<3;i++){
        if (sscanf(buffer,"%lf %lf %lf",vec+3*i,vec+3*i+1,vec+3*i+2)!=3)
          error_exit("Error reading lattice vectors");
        fdfreadline(buffer,LINE_SIZE,dir,&files);
      }

    }else if (!fdfstrcasecmp(&ptr,"latticeparameters")){
      fdfreadline(buffer,LINE_SIZE,dir,&files);
      abc=malloc(6*sizeof(double));
      if (!abc) error_exit("malloc error for abc");
      if (sscanf(buffer,"%lf %lf %lf %lf %lf %lf",
                 abc,abc+1,abc+2,abc+3,abc+4,abc+5)!=6)
          error_exit("Error reading lattice parameters");
      fdfreadline(buffer,LINE_SIZE,dir,&files);

    }else if (!fdfstrcasecmp(&ptr,"AtomicCoordinatesAndAtomicSpecies")){
      if (m->n==0)
        error_exit("AtomicCoordinatesAndAtomicSpecies preceeds NumberOfAtoms");
      atomspecies=malloc(m->n*sizeof(int));
      atomcoords=malloc(3*m->n*sizeof(double));
      if ((!atomspecies)||(!atomcoords)) error_exit("malloc error for atoms");
      for(i=0;i<m->n;i++){
        fdfreadline(buffer,LINE_SIZE,dir,&files);
        if (sscanf(buffer,"%lf %lf %lf %d",
                   atomcoords+3*i,atomcoords+3*i+1,atomcoords+3*i+2,
                   atomspecies+i)!=4)
          error_exit("Error reading atomic positions");
      }
      fdfreadline(buffer,LINE_SIZE,dir,&files);
    }
    else if (!fdfstrcasecmp(&ptr,"ChemicalSpeciesLabel")){
      if (!nspec)
        error_exit("ChemicalSpeciesLabel preceeds NumberOfSpecies");
      species_to_atno=malloc((nspec+1)*sizeof(int));
      species_to_charge=malloc((nspec+1)*sizeof(double));
      if ((!species_to_atno)||(!species_to_charge))
	error_exit("malloc error for species");
      for(i=0;i<nspec;i++){
        fdfreadline(buffer,LINE_SIZE,dir,&files);
        if (sscanf(buffer,"%d %d%n",&j,&k,&n)!=2)
          error_exit("Error reading ChemicalSpeciesLabel");
        if ((j<1)||(j>nspec)){
          error_exit("Invalid species number in ChemicalSpeciesLabel");
        }
        species_to_atno[j]=k;
        if (flags&CHDEN){ /* No point in attempting to find pseudo */
          ptr=buffer+n;   /* charges unless also processing elect charge */
          while((*ptr)&&(isspace(*ptr))) ptr++;
          ptr2=ptr;
          while((*ptr2)&&(!isspace(*ptr2))) ptr2++;
          *ptr2=0;
          species_to_charge[j]=siesta_pot_charge(prefix,ptr);
        }
      }
      fdfreadline(buffer,LINE_SIZE,dir,&files);
    }
    else if (!fdfstrcasecmp(&ptr,"DMInitSpin")){
      if (m->n==0)
        error_exit("DMInitSpin preceeds NumberOfAtoms");
      atomspins=malloc(m->n*sizeof(double));
      for(i=0;i<m->n;i++) atomspins[i]=0;
      while(1){
        fdfreadline(buffer,LINE_SIZE,dir,&files);
        i=sscanf(buffer,"%d%n",&j,&k);
        if (i==0) break;
        if ((j<1)||(j>m->n))
          error_exit("Atom index out of range in spin block");
        ptr=buffer+k;
        while (*ptr==' ') ptr++;
        if (!tokenmatch(&ptr,"+")){
          atomspins[j-1]=1;
        }
        else if (!tokenmatch(&ptr,"-")){
          atomspins[j-1]=-1;
        }
        else{
          if (sscanf(ptr,"%lf",atomspins+j-1)!=1)
            error_exit("Error parsing spin");
        }
      }
    }
    else if ((kp)&&(!fdfstrcasecmp(&ptr,"kgridMonkhorstPack"))){
      for(i=0;i<3;i++){
        fdfreadline(buffer,LINE_SIZE,dir,&files);
        j=sscanf(buffer,"%d %d %d %lf",mp_grid[i],mp_grid[i]+1,mp_grid[i]+2,
                 mp_off+i);
        if (j!=4){
          PARSE_ERROR;
          break;
        }
      }
      if (j!=4)
        fprintf(stderr,"Ignoring MP block\n");
      else{
        okay=1;
        /* Need grid to be diagonal */
        for(i=0;i<3;i++)
          for(j=0;j<3;j++)
            if (i==j) continue;
            else if (mp_grid[i][j]!=0) okay=0;
        if (okay){
          kp->mp=malloc(sizeof(struct mp_grid));
          if (!kp->mp) error_exit("Malloc error for MP grid");
          for(i=0;i<3;i++)
            kp->mp->grid[i]=mp_grid[i][i];
          for(i=0;i<3;i++){
            if ((kp->mp->grid[i]&1)==0) mp_off[i]-=0.5;
            mp_off[i]/=kp->mp->grid[i];
            kp->mp->disp[i]=mp_off[i];
          }
        }
        else{
          fprintf(stderr,"c2x is unable to process non-diagonal MP grids\n");
        }
      }
    }
    else{
      fprintf(stderr,"Ignoring %s\n",ptr);
      while(fdfreadline(buffer,LINE_SIZE,dir,&files)){
        ptr=buffer;
        if (!strncasecmp(ptr,"%endblock",9)) break;
      }
    }
  }

  /* Extract unit cell */
  if (vec){
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
        c->basis[i][j]=vec[3*i+j]*latticeconstant;
    free(vec);
    vec=NULL;
  }
  else if (abc){
    for(i=0;i<3;i++)
      abc[i]*=latticeconstant;
    abc2cart(abc,c);
    free(abc);
    abc=NULL;
  }
  else{
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
        c->basis[i][j]=(i==j)?latticeconstant:0;
  }
  real2rec(c);

  /* Extract atomic positions */

  if (!atomcoords)
    error_exit("No atomic coordinates found!");
  
  m->atoms=malloc(m->n*sizeof(struct atom));
  if (!m->atoms) error_exit("malloc error for atoms");
  init_atoms(m->atoms,m->n);
  
  scale=-1;
  ptr=format;
  if (!format)  /* default */
    scale=BOHR;
  else if ((!fdftokenmatch(&ptr,"Bohr"))||
           (!fdftokenmatch(&ptr,"NotScaledCartesianBohr")))
    scale=BOHR;
  else if ((!fdftokenmatch(&ptr,"Ang"))||
           (!fdftokenmatch(&ptr,"NotScaledCartesianAng")))
    scale=1;
  else if (!fdftokenmatch(&ptr,"ScaledCartesian"))
    scale=latticeconstant;
  else if ((!fdftokenmatch(&ptr,"Fractional"))||
           (!fdftokenmatch(&ptr,"ScaledByLatticeVectors")))
    scale=0;

  if (scale==-1){
    fprintf(stderr,"Unrecognised atomic coord format: %s\n",format);
    exit(1);
  }
  if (format) free(format);

  if (scale){
    for(i=0;i<m->n;i++)
      for(j=0;j<3;j++)
        m->atoms[i].abs[j]=(atomcoords[3*i+j]+shift[j])*scale;
    addfrac(m->atoms,m->n,c->recip);
  }
  else{
    for(i=0;i<m->n;i++)
      for(j=0;j<3;j++)
        m->atoms[i].frac[j]=atomcoords[3*i+j]+shift[j];
    addabs(m->atoms,m->n,c->basis);
  }    
  free(atomcoords);
  
  /* Add atomic numbers */

  if (!atomspecies)
    error_exit("chemical_species_label block not found");
  
  for(i=0;i<m->n;i++){
    m->atoms[i].atno=species_to_atno[atomspecies[i]];
    m->atoms[i].chg=species_to_charge[atomspecies[i]];
  }

  m->nspec=nspec;
  m->spec=malloc(nspec*sizeof(struct species));
  if (!m->spec) error_exit("malloc error for m->spec in fdf_read");
  for(i=0;i<nspec;i++) m->spec[i].atno=species_to_atno[i+1];
  
  free(atomspecies);
  free(species_to_atno);
  dict_add(m->dict,"Siesta_species_to_charge",species_to_charge);

  /* Add spins */
  
  if (atomspins){
    for(i=0;i<m->n;i++)
      m->atoms[i].spin=atomspins[i];
    free(atomspins);
  }
  
  /* Miscellaneous */

  if (dmetol) {
    e->etol=(*dmetol)/m->n;
    free(dmetol);
  }

}

int fdfreadline(char *buffer, int len, char *dir, struct infiles **filesp){
  int i,off,success2,inc;
  char *ptr,*success,*p2,token[LINE_SIZE];
  
  while((success=fgets(buffer,len,(*filesp)->f))){ /* fgets() always
                                                    null terminates,
                                                    gcc likes extra brackets */
    (*filesp)->line++;

/* Kill trailing spaces and newlines / carriage returns */
    ptr=buffer+strlen(buffer)-1;
    while((ptr>=buffer)&&((*ptr==' ')||(*ptr=='\n')||(*ptr=='\r'))) ptr--;
    *(ptr+1)=0;

/* Eat leading spaces */
    ptr=buffer;
    while(isspace(*ptr)) ptr++;
/* Skip comments and blank lines */
    if ((*ptr=='#')||(*ptr=='!')||(*ptr==';')||(*ptr==0)) continue;
    break;
  }

  if (!success){ /* Need to recurse out of nested include files */
    success2=0;
    while ((*filesp)->last){
      fclose((*filesp)->f);
      free((*filesp)->name);
      inc=(*filesp)->include;
      (*filesp)=(*filesp)->last;
      free((*filesp)->next);
      (*filesp)->next=NULL;
      if (!inc){
        strcpy(buffer,"%endblock\n");
        return 1;
      }
      success2=fdfreadline(buffer,LINE_SIZE,dir,filesp);
      if (success2) return success2;
    }
    if (!success2) return(0);
  }

  /* Shift buffer to remove leading whitespace */
  off=ptr-buffer;
  if (off){
    while(*ptr) {*(ptr-off)=*ptr;ptr++;}
    *(ptr-off)=0;
  }

  ptr=buffer;
  
  /* return read %endblock line, but reset files to old value */
  if (((*filesp)->include==2)&&(!strncasecmp(ptr,"%endblock",9))){
    fdf_unwind(filesp);
    return 1;
  }

  if (!tokenmatch(&ptr,"%include")){
    while (*ptr==' ') ptr++;
    p2=ptr;
    while ((*p2)&&(*p2!=' ')) p2++;
    *p2=0;
    include_file(filesp,dir,ptr);
    return fdfreadline(buffer,LINE_SIZE,dir,filesp);
  }

  /* Might have "label < include" syntax */

  /* consume label */
  p2=ptr;
  /* ought to add $%&@/^~ */
  while ((isalnum(*p2))||(*p2=='.')||(*p2=='_')||(*p2=='-')||(*p2=='+')) p2++;

  /* and following whitespace */
  while (isspace(*p2)) p2++;

  /* do we now have a < ? */
  if (*p2=='<'){
    for(i=0;!(isspace(buffer[i]));i++)
      token[i]=buffer[i];
    token[i]=0;
    ptr=p2+1;
    while (*ptr==' ') ptr++;
    p2=ptr;
    while ((*p2)&&(*p2!=' ')) p2++;
    *p2=0;

    i=fdfscanfile(ptr,dir,buffer,token,filesp);
    if (!i){
      fprintf(stderr,"Warning: %s not found in < included file %s\n",
	      token,ptr);
      return fdfreadline(buffer,LINE_SIZE,dir,filesp);
    }
    else return 1;
  }
    
  
  if (!tokenmatch(&ptr,"%block")){ /* Need to scan for < includes */
    while(*(++ptr)){
      if (*ptr=='<'){
        *ptr=0; /* stop rest of parser seeing the < part */
        ptr++;
        while (*ptr==' ') ptr++;
        p2=ptr;
        while ((*p2)&&(*p2!=' ')) p2++;
        *p2=0;
        include_file(filesp,dir,ptr);
        (*filesp)->include=0;
        return 1; /* returns buffer from old file with %block line
                   * next read is from new file
                   */
      }
    }
  }
  
  return (1);
}

void include_file(struct infiles **file, char *dir, char *ptr){
  struct infiles *fs;
  char *name,*infile;

  name=NULL;
  infile=ptr;
  dequote(&infile);
  ptr=infile;
  
  fs=*file;
  fs->next=malloc(sizeof(struct infiles));
  if (!fs->next) error_exit("Malloc error for files_next in cell_read");
  fs->next->last=fs;
  fs->next->next=NULL;
  fs=fs->next;
  fs->f=NULL;
  fs->count=(*file)->count+1;
  if (fs->count>30)
    error_exit("maximum include depth of 30 exceeded -- infinite recursion?");
  
  if (dir){
    name=malloc(strlen(dir)+strlen(ptr)+1);
    if (!name) error_exit("malloc error for filename");
    strcpy(name,dir);
    strcat(name,ptr);
    fs->f=fopen(name,"r");
    fs->name=name;
  }
  if (!fs->f) {
    fs->f=fopen(ptr,"r");
    fs->name=malloc(strlen(ptr)+1);
    if (!fs->name) error_exit("malloc error for filename");
    strcpy(fs->name,ptr);
  }
  if (!fs->f){
    fprintf(stderr,"Error: failed to open include file ");
    if (name) fprintf(stderr,"%s or",name);
    fprintf(stderr,"%s\n",ptr);
    
  }
  if (debug)
    fprintf(stderr,"Opened included file %s\n",fs->name);

  fs->line=0;
  fs->include=1;  /* This gets set to zero if a block include which
                     * needs an %endblock added */   
  if ((*file)->include)
  fs->include=(*file)->include;
  fs->count=(*file)->count+1;
  fs->ret=(*file)->ret;
  *file=fs;
}

static int fdfreadlength(char *buff, double *x){
  char *p,*ptr;
  int n,i;  

  while((*buff)&&((*buff==' '))) buff++;
  /*
   * n=sscanf(buff,"%lf %ms",x,&ptr);
   */
  n=sscanf(buff,"%lf%n",x,&i);
  n+=sscanfmsn(buff+i,&ptr,NULL);
  p=ptr; /* tokenmatch may modify ptr, but we wish to free it */
  
  if (n==0) return 0;

  if (n==1) {
    error_exit("No units specified for length");
  }

  if (!tokenmatch(&p,"bohr")){
    (*x)*=BOHR;
    free(ptr);
    return 1;
  }

  if (!tokenmatch(&p,"ang")){
    free(ptr);
    return 1;
  }

  if (!tokenmatch(&p,"nm")){
    (*x)*=10;
    free(ptr);
    return 1;
  }

  fprintf(stderr,"Unexpected length unit: %s\n",ptr);
  free(ptr);
  return 0;

}


static int fdfreadenergy(char *buff, double *x){
  char *p,*ptr;
  int n,i;  

  while((*buff)&&((*buff==' '))) buff++;
  /*
   * n=sscanf(buff,"%lf %ms",x,&ptr);
   */
  n=sscanf(buff,"%lf%n",x,&i);
  n+=sscanfmsn(buff+i,&ptr,NULL);
  p=ptr; /* tokenmatch may modify ptr, but we wish to free it */
  
  if (n==0) return 0;

  if (n==1) {
    error_exit("No units specified for energy");
  }

  if (!tokenmatch(&p,"Ry")){
    (*x)*=0.5*H_eV;
    free(ptr);
    return 1;
  }

  if (!tokenmatch(&p,"mRy")){
    (*x)*=0.5*0.001*H_eV;
    free(ptr);
    return 1;
  }

  if (!tokenmatch(&p,"eV")){
    free(ptr);
    return 1;
  }

  if (!tokenmatch(&p,"meV")){
    (*x)*=0.001;
    free(ptr);
    return 1;
  }

  if (!tokenmatch(&p,"Hartree")){
    (*x)*=H_eV;
    free(ptr);
    return 1;
  }

  if (!tokenmatch(&p,"mHartree")){
    (*x)*=0.001*H_eV;
    free(ptr);
    return 1;
  }

  fprintf(stderr,"Unexpected length unit: %s\n",ptr);
  free(ptr);
  return 0;

}


/* Compare two strings up to the length of the second.
 * If they compare equal, move the first pointer to the end of
 * the equal section.
 * Ignore the three characters -_. in s1 when comparing
 */
static int fdfstrcasecmp(char **s1_in, char *s2){
  char *s1;

  s1=*s1_in;
  while((*s1)&&(*s2)){
    if ((*s1=='_')||(*s1=='-')||(*s1=='.'))
      s1++;
    else if ((*s2=='_')||(*s2=='-')||(*s2=='.'))
      s2++;
    else if (toupper(*s1)==toupper(*s2)){
      s1++;
      s2++;
    }
    else break;
  }

  if (*s2==0){
    *s1_in=s1;
    return 0;
  }
  
  if ((*s1)>(*s2)) return 1;
  return -1;

}

static int fdftokenmatch(char **s1_in, char *s2){
  char *s1;
  int tmp;

  s1=*s1_in;
  tmp=fdfstrcasecmp(&s1,s2);
  if (tmp!=0) return tmp;

  /* s1 will have been moved to end of the equal section */

  if ((*s1==0)||(isspace(*s1))){
    *s1_in=s1;
    return 0;
  }
  return 1;
}

/* Read charge from .psf file */
double siesta_pot_charge(char *prefix, char *label){
  char *name;
  int i,idummy;
  double dummy,charge;
  FILE *pseudo;
  char buffer[LINE_SIZE+1];

  if (prefix)
    name=malloc(strlen(prefix)+strlen(label)+10);
  else
    name=malloc(strlen(label)+10);

  name[0]=0;
  if (prefix) strcpy(name,prefix);
  strcat(name,label);
  strcat(name,".psf");

  pseudo=fopen(name,"rb");

  if (!pseudo) {
    name[0]=0;
    if (prefix) strcpy(name,prefix);
    strcat(name,label);
    strcat(name,".out.psf");
    pseudo=fopen(name,"rb");
  }

  if (!pseudo){
    fprintf(stderr,"Failed to open %s\n",name);
    free(name);
    return 0;
  }

  for(i=0;i<4;i++)
    fgets(buffer,LINE_SIZE,pseudo);

  i=sscanf(buffer,"%d %d %d %lf %lf %lf",&idummy,&idummy,&idummy,
	   &dummy,&dummy,&charge);

  if (i!=6) {
    fprintf(stderr,"Failed to read charge from %s\n",name);
    charge=0;
  }
  else if (debug)
    fprintf(stderr,"Label %s charge %lf read from %s\n",
	    label,charge,name);
    
  fclose(pseudo);
  free(name);

  return charge;
}

/* scan file(s) for a given token in response to "<" syntax */

static int fdfscanfile(char *filename, char *dirname, char *buffer,
                       char *token, struct infiles **infile){
  char *ptr;
  int found,block;
  struct infiles *files;

  if (debug>2)
    fprintf(stderr,"fdfscanfile for %s in %s called\n",
	    token,filename);

  block=found=0;

  files=*infile;
  include_file(infile,dirname,filename);
  (*infile)->ret=files;
  if ((*infile)->f==NULL){
    free((*infile)->name);
    free(*infile);
    *infile=files;
    return 0;
  }
  
  while(fdfreadline(buffer,LINE_SIZE,dirname,infile)){
    ptr=buffer;
    block=0;
    if (!strncasecmp(ptr,"%block",6)){
      block=1;
      ptr+=6;
      while(isspace(*ptr)) ptr++;
    }
    if (!fdftokenmatch(&ptr,token)){
      found=1;
      break;
    }
  }

  if (debug>1)
    fprintf(stderr,"fdfscanfile for %s in %s returns %d block=%d %s\n",
	    token,(*infile)->last->name,found,block,buffer);

  if (block&&found){
    (*infile)->include=2;
  }
  else{
    fdf_unwind(infile);
  }
    
  return found;
}

/* Deal with potentially multiple back-outs caused by Siesta's "<" syntax */
static void fdf_unwind(struct infiles **infile){
  struct infiles *target;

  target=(*infile)->ret;

  if (!target) error_exit("asked to unwind to null in fdf_unwind");
  
  while(*infile!=target){
    fclose((*infile)->f);
    free((*infile)->name);
    *infile=(*infile)->last;
    free((*infile)->next);
  }
}

void dequote(char **str){
  char *p;

  p=*str;
  
  if (((*p)=='"')&&(strchr(p+1,'"'))){
    p++;
    *(strchr(p,'"'))=0;
    *str=p;
    return;
  }

  if (((*p)=='\'')&&(strchr(p+1,'\''))){
    p++;
    *(strchr(p,'\''))=0;
    *str=p;
    return;
  }
}
