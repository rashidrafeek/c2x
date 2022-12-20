/* Read some useful data from a SIESTA .fdf file
 * 
 */

/* To do: add DM.InitSpin block
 *        use PARSE_ERROR
 *        test!!
 *
 * This version does not parse quotation marks properly.
 */

/* Note that Siesta ignores any of the three characters -_. appearing
 * in its keywords.
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
#include<string.h>
#include<ctype.h>
#include<errno.h>

#include "c2xsf.h"

#define LINE_SIZE 200

#define PARSE_ERROR  fprintf(stderr,"Parse error on line %d %s\n",files->line,(files->name)?files->name:"")

static int fdfreadline(char *buffer, int len);
static int fdfreadlength(char *buff, double *x);
static int fdfreadenergy(char *buff, double *x);
static int fdfstrcasecmp(char **s1_in, char *s2);
static int fdftokenmatch(char **s1_in, char *s2);
static void fdfinclude(char *ptr);

static char *ptr;

static struct infiles {
  FILE* f; struct infiles *next; struct infiles *last; char *name;
  int line; int include;
} *files;
static FILE *infile;


void fdf_read(FILE* in, struct unit_cell *c, struct contents *m,
              struct kpts *kp, struct es *e){
  int i,j,k,*atomspecies;
  int nspec,*species_to_atno;
  double scale,shift[3];
  char *ptr2,*format;
  static char buffer[LINE_SIZE+1];
  double *vec,*abc,*atomcoords,*dmetol,latticeconstant,*atomspins;
  
  files=malloc(sizeof(struct infiles));
  if (!files) error_exit("Malloc error for files in cell_read");
  files->f=in;
  files->next=files->last=NULL;
  files->name=NULL;
  files->line=0;
  infile=in;
  
  m->n=0;
  nspec=0;
  vec=abc=NULL;
  latticeconstant=0;
  format=NULL;
  species_to_atno=NULL;
  atomspecies=NULL;
  shift[0]=shift[1]=shift[2]=0;
  dmetol=NULL;
  atomspins=NULL;
  atomcoords=NULL;
  
  if (debug>2) fprintf(stderr,"fdf_read called\n");

  if (!(c->basis=malloc(72)))
    error_exit("Malloc error in fdfread for c->basis");

  while(fdfreadline(buffer,LINE_SIZE)){
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
    

/* What remains ought to be a block or a keyword.  */
    if (strncasecmp(ptr,"%block",6)) continue;
    ptr+=6;
/* Eat leading spaces */
    while(*ptr==' ') ptr++;
/* Kill trailing spaces and things after # */
    ptr2=ptr;
    while((*ptr2)&&(*ptr2!=' ')&&(*ptr2!='#')) ptr2++;
    *ptr2=0; /* Either it was a null, or it should be one... */

    if (debug>2) fprintf(stderr,"Found a %s block\n",ptr);

    if (!fdfstrcasecmp(&ptr,"latticevectors")){
      fdfreadline(buffer,LINE_SIZE);
      vec=malloc(9*sizeof(double));
      for(i=0;i<3;i++){
        if (sscanf(buffer,"%lf %lf %lf",vec+3*i,vec+3*i+1,vec+3*i+2)!=3)
          error_exit("Error reading lattice vectors");
        fdfreadline(buffer,LINE_SIZE);
      }

    }else if (!fdfstrcasecmp(&ptr,"latticeparameters")){
      fdfreadline(buffer,LINE_SIZE);
      abc=malloc(6*sizeof(double));
      if (sscanf(buffer,"%lf %lf %lf %lf %lf %lf",
                 abc,abc+1,abc+2,abc+3,abc+4,abc+5)!=6)
          error_exit("Error reading lattice parameters");
      fdfreadline(buffer,LINE_SIZE);

    }else if (!fdfstrcasecmp(&ptr,"AtomicCoordinatesAndAtomicSpecies")){
      if (m->n==0)
        error_exit("AtomicCoordinatesAndAtomicSpecies preceeds NumberOfAtoms");
      atomspecies=malloc(m->n*sizeof(int));
      atomcoords=malloc(3*m->n*sizeof(double));
      for(i=0;i<m->n;i++){
        fdfreadline(buffer,LINE_SIZE);
        if (sscanf(buffer,"%lf %lf %lf %d",
                   atomcoords+3*i,atomcoords+3*i+1,atomcoords+3*i+2,
                   atomspecies+i)!=4)
          error_exit("Error reading atomic positions");
      }
      fdfreadline(buffer,LINE_SIZE);
    }
    else if (!fdfstrcasecmp(&ptr,"ChemicalSpeciesLabel")){
      if (!nspec)
        error_exit("ChemicalSpeciesLabel preceeds NumberOfSpecies");
      species_to_atno=malloc((nspec+1)*sizeof(int));
      for(i=0;i<nspec;i++){
        fdfreadline(buffer,LINE_SIZE);
        if (sscanf(buffer,"%d %d",&j,&k)!=2)
          error_exit("Error reading ChemicalSpeciesLabel");
        if ((j<1)||(j>nspec)){
          error_exit("Invalid species number in ChemicalSpeciesLabel");
        }
        species_to_atno[j]=k;
      }
      fdfreadline(buffer,LINE_SIZE);
    }
    else if (!fdfstrcasecmp(&ptr,"DMInitSpin")){
      if (m->n==0)
        error_exit("DMInitSpin preceeds NumberOfAtoms");
      atomspins=malloc(m->n*sizeof(double));
      for(i=0;i<m->n;i++) atomspins[i]=0;
      while(1){
        fdfreadline(buffer,LINE_SIZE);
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
    else{
      while(fdfreadline(buffer,LINE_SIZE)){
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
  
  for(i=0;i<m->n;i++)
    m->atoms[i].atno=species_to_atno[atomspecies[i]];
  free(atomspecies);
  free(species_to_atno);

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

int fdfreadline(char *buffer, int len){
  int off,success2,inc;
  char *ptr,*success,*p2;

  while((success=fgets(buffer,len,infile))){ /* fgets() always
                                                    null terminates,
                                                    gcc likes extra brackets */
    files->line++;

/* Kill trailing spaces and newlines / carriage returns */
    ptr=buffer+strlen(buffer)-1;
    while((ptr>=buffer)&&((*ptr==' ')||(*ptr=='\n')||(*ptr=='\r'))) ptr--;
    *(ptr+1)=0;

/* Eat leading spaces */
    ptr=buffer;
    while(*ptr==' ') ptr++;
/* Skip comments and blank lines */
    if ((*ptr=='#')||(*ptr==0)) continue;
    break;
  }

  if (!success){ /* Need to recurse out of nested include files */
    success2=0;
    while (files->last){
      fclose(infile);
      free(files->name);
      inc=files->include;
      files=files->last;
      free(files->next);
      files->next=NULL;
      infile=files->f;
      if (!inc){
        strcpy(buffer,"%endblock\n");
        return 1;
      }
      success2=fdfreadline(buffer,LINE_SIZE);
      if (success2) return success2;
    }
    if (!success2) return(0);
  }

  off=ptr-buffer;
  if (off){
    while(*ptr) {*(ptr-off)=*ptr;ptr++;}
    *(ptr-off)=0;
  }

  if (!tokenmatch(&ptr,"%include")){
    fdfinclude(ptr);
    return fdfreadline(buffer,LINE_SIZE);
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
    ptr=p2+1;
    fdfinclude(ptr);
    return fdfreadline(buffer,LINE_SIZE);
  }
    
  
  if (!tokenmatch(&ptr,"%block")){ /* Need to scan for < includes */
    while(*(++ptr)){
      if (*ptr=='<'){
        *ptr=0; /* stop rest of parser seeing the < part */
        ptr++;
        fdfinclude(ptr);
        files->include=0;
        return 1; /* returns buffer from old file with %block line
                   * next read is from new file
                   */
      }
    }
  }
  
  return (1);
}

static void fdfinclude(char *ptr){
  char *p2;

  while (*ptr==' ') ptr++;
  p2=ptr;
  while ((*p2)&&(*p2!=' ')) p2++;
  *p2=0;
  files->next=malloc(sizeof(struct infiles));
  if (!files) error_exit("Malloc error for files_next in cell_read");
  files->next->last=files;
  files->next->next=NULL;
  files->next->f=fopen(ptr,"r");
  files=files->next;
  files->name=malloc(strlen(ptr)+1);
  strcpy(files->name,ptr);
  files->line=0;
  files->include=1; /* This gets set to zero if a block include which
                     * needs an %endblock added */
  infile=files->f;
  if (!files->f){
    fprintf(stderr,"Error, unable to open %s in fdf_read\n",ptr);
    exit(1);
  }
  if (debug) fprintf(stderr,"Opened included file %s\n",ptr);
}

static int fdfreadlength(char *buff, double *x){
  char *p,*ptr;
  int n;  

  while((*buff)&&((*buff==' '))) buff++;
  
  n=sscanf(buff,"%lf %ms",x,&ptr);
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
  int n;  

  while((*buff)&&((*buff==' '))) buff++;
  
  n=sscanf(buff,"%lf %ms",x,&ptr);
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
 * Ignore the three characters -_. when comparing
 */
static int fdfstrcasecmp(char **s1_in, char *s2){
  char *s1;

  s1=*s1_in;
  while((*s1)&&(*s2)){
    if ((*s1=='_')||(*s1=='-')||(*s1=='.'))
      s1++;
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
