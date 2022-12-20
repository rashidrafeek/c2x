/* Read an abinit .in file */

/* Abinit does not regard end of line as special -- there may be
 * multiple keywords on a single line, or data for a single keyword
 * spread across multiple lines.
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
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<math.h>

#include "c2xsf.h"

#define LINE_SIZE 132  /* Thus says abinit */

static int ab_len_read(char *line, int off, double *data, int n,
		       FILE* infile);
static int ab_energy_read(char *line, int off, double *data, int n,
			  FILE* infile);
static int ab_float_read(char *line, int off, double *data, int n,
			  FILE* infile);
static int ab_int_read(char *line, int off, int *data, int n,
		       FILE* infile);
static int ab_readline(char *buffer, char **p, FILE* infile);

/* from cif_read */
void sym_expand(struct unit_cell *c, struct contents *m, struct symmetry *s);

void qe_read(FILE* infile, struct unit_cell *c, struct contents *m,
             struct kpts *k, struct symmetry *s, struct es *e);

void abinit_in_read(FILE* infile, struct unit_cell *c, struct contents *m,
		    struct kpts *k, struct symmetry *s, struct es *e){
  
  double *acell,*rprim,*scalecart,*znucl,*xcart,*xred,*kpt,*wtk,*shiftk;
  double *ecut,*toldfe,*spinat,*tnons,*kptnrm;
  double x,xa,xc,angles[3],mat[3][3];
  char line[LINE_SIZE+1];
  char *ptr;
  int i,j,jj;
  int *natom,*typat,*ntypat,*nkpt,*ngkpt,*nshiftk,*kptopt,*ndtset;
  int *nspden,*nsppol,*nsym,*symrel,*natrd,*spgroup,*spgaxor,*spgorig;
  int hall;


  acell=rprim=scalecart=znucl=xcart=xred=kpt=wtk=shiftk=NULL;
  ecut=toldfe=spinat=tnons=kptnrm=NULL;
  natom=typat=ntypat=nkpt=ngkpt=nshiftk=kptopt=ndtset=NULL;
  nspden=nsppol=nsym=symrel=natrd=spgroup=spgaxor=spgorig=NULL;
  ptr=line;
  *ptr=0;
  
  /* Read file in two passes. First pass scalars and fixed-size arrays,
   * second pass arrays whose size depends on scalars read from the file
   */
  
  while(1){
    /* remove leading spaces */
    while(*ptr&&(isspace(*ptr))) ptr++;
    if (*ptr==0){
      if (!ab_readline(line,&ptr,infile)) break;
    }

    if (!tokenmatch(&ptr,"&control")){
      rewind(infile);
      qe_read(infile,c,m,k,s,e);
      return;
    }
    
    if (!tokenmatch(&ptr,"acell")){
      if (!acell) acell=malloc(3*sizeof(double));
      ptr=line+ab_len_read(line,ptr-line,acell,3,infile);
    }
    else if (!tokenmatch(&ptr,"rprim")){
      if (!rprim) rprim=malloc(9*sizeof(double));
      ptr=line+ab_float_read(line,ptr-line,rprim,9,infile);
    }
    else if (!tokenmatch(&ptr,"scalecart")){
      if (!scalecart) scalecart=malloc(3*sizeof(double));
      ptr=line+ab_float_read(line,ptr-line,scalecart,3,infile);
    }
    else if (!tokenmatch(&ptr,"angdeg")){
      ptr=line+ab_float_read(line,ptr-line,angles,3,infile);
      if (!rprim) rprim=malloc(9*sizeof(double));
      if ((angles[0]==angles[1])&&(angles[0]==angles[2])){
	if (angles[0]==90){
	  rprim[0]=1;
	  rprim[1]=0;
	  rprim[2]=0;
	  rprim[3]=0;
	  rprim[4]=1;
	  rprim[5]=0;
	  rprim[6]=0;
	  rprim[7]=0;
	  rprim[8]=1;
	}
	else{
	  xc=sqrt((2*cos(M_PI*angles[0]/180)+1)/3);
	  xa=sqrt(1-xc*xc);
	  rprim[0]=xa;
	  rprim[1]=0;
	  rprim[2]=xc;
	  rprim[3]=-xa/2;
	  rprim[4]=xa*sqrt(3.0)/2;
	  rprim[5]=xc;
	  rprim[6]=-xa/2;
	  rprim[7]=-rprim[4];
	  rprim[8]=xc;
	}
      }
      else{
	for(i=0;i<3;i++)angles[i]*=M_PI/180;
	rprim[0]=1;
	rprim[1]=0;
	rprim[2]=0;
	rprim[3]=cos(angles[2]);
	rprim[4]=sin(angles[2]);
	rprim[5]=0;
	rprim[6]=cos(angles[1]);
	x=cos(angles[0])-rprim[0]*rprim[3];
	rprim[7]=0;
	if (fabs(rprim[7])>1e-20){
	  if (fabs(rprim[4])>1e-30)
	    rprim[7]=x/rprim[4];
	  else
	    error_exit("Unable to expand angdeg to rprim");
	}
	rprim[8]=sqrt(1-rprim[6]*rprim[6]-rprim[7]*rprim[7]);
      }	  
    }
    else if (!tokenmatch(&ptr,"ndtset")){
      if (!ndtset) ndtset=malloc(sizeof(int));
      ptr=line+ab_int_read(line,ptr-line,ndtset,1,infile);
      if (*ndtset!=1)
	error_exit("Cannot parse abinit .in files with multiple datasets");
    }
    else if (!tokenmatch(&ptr,"natom")){
      natom=malloc(sizeof(int));
      ptr=line+ab_int_read(line,ptr-line,natom,1,infile);
    }
    else if (!tokenmatch(&ptr,"natrd")){
      natrd=malloc(sizeof(int));
      ptr=line+ab_int_read(line,ptr-line,natrd,1,infile);
    }
    else if (!tokenmatch(&ptr,"ntypat")){
      ntypat=malloc(sizeof(int));
      ptr=line+ab_int_read(line,ptr-line,ntypat,1,infile);
    }
    else if (!tokenmatch(&ptr,"nsppol")){
      nsppol=malloc(sizeof(int));
      ptr=line+ab_int_read(line,ptr-line,nsppol,1,infile);
    }
    else if (!tokenmatch(&ptr,"nspden")){
      nspden=malloc(sizeof(int));
      ptr=line+ab_int_read(line,ptr-line,nspden,1,infile);
    }
    else if (!tokenmatch(&ptr,"nkpt")){
      nkpt=malloc(sizeof(int));
      ptr=line+ab_int_read(line,ptr-line,nkpt,1,infile);
    }
    else if (!tokenmatch(&ptr,"ngkpt")){
      ngkpt=malloc(3*sizeof(int));
      ptr=line+ab_int_read(line,ptr-line,ngkpt,3,infile);
    }
    else if (!tokenmatch(&ptr,"kptopt")){
      kptopt=malloc(sizeof(int));
      ptr=line+ab_int_read(line,ptr-line,kptopt,1,infile);
    }
    else if (!tokenmatch(&ptr,"kptnrm")){
      kptnrm=malloc(sizeof(double));
      ptr=line+ab_float_read(line,ptr-line,kptnrm,1,infile);
    }
    else if (!tokenmatch(&ptr,"nshiftk")){
      nshiftk=malloc(sizeof(int));
      ptr=line+ab_int_read(line,ptr-line,nshiftk,1,infile);
    }
    else if (!tokenmatch(&ptr,"nsym")){
      nsym=malloc(sizeof(int));
      ptr=line+ab_int_read(line,ptr-line,nsym,1,infile);
    }
    else if (!tokenmatch(&ptr,"spgroup")){
      spgroup=malloc(sizeof(int));
      ptr=line+ab_int_read(line,ptr-line,spgroup,1,infile);
    }
    else if (!tokenmatch(&ptr,"spgaxor")){
      spgaxor=malloc(sizeof(int));
      ptr=line+ab_int_read(line,ptr-line,spgaxor,1,infile);
    }
    else if (!tokenmatch(&ptr,"spgorig")){
      spgorig=malloc(sizeof(int));
      ptr=line+ab_int_read(line,ptr-line,spgorig,1,infile);
    }
    else if (!tokenmatch(&ptr,"ecut")){
      if (!ecut) ecut=malloc(sizeof(double));
      else
	fprintf(stderr,
		"Warning: ecut appears more than once. Last one used");
      ptr=line+ab_energy_read(line,ptr-line,ecut,1,infile);
    }
    else if (!tokenmatch(&ptr,"toldfe")){
      toldfe=malloc(sizeof(double));
      ptr=line+ab_energy_read(line,ptr-line,toldfe,1,infile);
    }
    else if (!tokenmatch(&ptr,"charge")){
      e->charge=malloc(sizeof(double));
      ptr=line+ab_float_read(line,ptr-line,e->charge,1,infile);
    }
    else{
      i=0;
      while (*(ptr+i)&&(!isspace(*(ptr+i)))) i++;
      if ((debug)&&(isalpha(*ptr)))
	fprintf(stderr,"Warning, ignoring token %.*s\n",i,ptr);
      ptr+=i;
    }
  }

  /* Now read variable-length arrays */

  if (!natrd) natrd=natom;
  if (!natrd) error_exit("natom/natrd not found");

  rewind(infile);
  
  while(1){
    /* remove leading spaces */
    while(*ptr&&(isspace(*ptr))) ptr++;
    if (*ptr==0){
      if (!ab_readline(line,&ptr,infile)) break;
    }

    if (!tokenmatch(&ptr,"typat")){
      typat=malloc(*natrd*sizeof(int));
      ptr=line+ab_int_read(line,ptr-line,typat,*natrd,infile);
    }
    else if (!tokenmatch(&ptr,"xangst")){
      xcart=malloc(*natrd*3*sizeof(double));
      ptr=line+ab_float_read(line,ptr-line,xcart,3*(*natrd),infile);
    }
    else if (!tokenmatch(&ptr,"xcart")){
      xcart=malloc(*natrd*3*sizeof(double));
      ptr=line+ab_len_read(line,ptr-line,xcart,3*(*natrd),infile);
    }
    else if (!tokenmatch(&ptr,"xred")){
      xred=malloc(*natrd*3*sizeof(double));
      ptr=line+ab_float_read(line,ptr-line,xred,3*(*natrd),infile);
    }
    else if (!tokenmatch(&ptr,"znucl")){
      znucl=malloc(*ntypat*sizeof(double));
      ptr=line+ab_float_read(line,ptr-line,znucl,*ntypat,infile);
    }
    else if (!tokenmatch(&ptr,"spinat")){
      spinat=malloc(3*(*natrd)*sizeof(double));
      ptr=line+ab_float_read(line,ptr-line,spinat,*natrd*3,infile);
    }
    else if (!tokenmatch(&ptr,"kpt")){
      if (!nkpt) error_exit("kpt found without nkpt");
      kpt=malloc(*nkpt*3*sizeof(double));
      ptr=line+ab_float_read(line,ptr-line,kpt,3*(*nkpt),infile);
    }
    else if (!tokenmatch(&ptr,"wtk")){
      if (!nkpt) error_exit("wtk found without nkpt");
      wtk=malloc(*nkpt*sizeof(double));
      ptr=line+ab_float_read(line,ptr-line,wtk,*nkpt,infile);
    }
    else if (!tokenmatch(&ptr,"shiftk")){
      shiftk=malloc(*nshiftk*3*sizeof(double));
      ptr=line+ab_float_read(line,ptr-line,shiftk,3*(*nshiftk),infile);
    }
    else if (!tokenmatch(&ptr,"symrel")){
      if (!nsym) error_exit("symrel found without nsym");
      symrel=malloc(*nsym*9*sizeof(int));
      ptr=line+ab_int_read(line,ptr-line,symrel,9*(*nsym),infile);
    }
    else if (!tokenmatch(&ptr,"tnons")){
      if (!nsym) error_exit("tnons found without nsym");
      tnons=malloc(*nsym*3*sizeof(double));
      ptr=line+ab_float_read(line,ptr-line,tnons,3*(*nsym),infile);
    }
    else{ /* skip token */
      i=0;
      while (*(ptr+i)&&(!isspace(*(ptr+i)))) i++;
      ptr+=i;
    }
  }
  
  /* Now convert abinit variables to ours */

  /* basis */
  if (!(c->basis=malloc(9*sizeof(double))))
    error_exit("Malloc error for basis");

  if (rprim)
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	c->basis[i][j]=rprim[3*i+j];
  else  /* default value of rprim is unit axes */
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	c->basis[i][j]=(i==j)?1:0;
  
  if (acell)
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	c->basis[i][j]*=acell[i];
  else  /* default value for acell is 1 in Bohr */
    for(i=0;i<3;i++)
      for(j=0;j<3;j++) c->basis[i][j]*=BOHR;

  if (scalecart)
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	c->basis[i][j]*=scalecart[j];
  
  real2rec(c);

  if (rprim) free(rprim);
  if (acell) free(acell);
  if (scalecart) free(scalecart);
  
  /* atomic co-ordinates */
  
  if ((xred)&&(xcart)) error_exit("Both xred and xcart specified");
  if ((xred)||(xcart)){
    m->atoms=malloc(*natrd*sizeof(struct atom));
    if (!m->atoms) error_exit("malloc error for atoms");
    m->n=*natrd;
    init_atoms(m->atoms,m->n);
    
    if (xred){
      for(i=0;i<m->n;i++){
	for (j=0;j<3;j++)
	  m->atoms[i].frac[j]=xred[3*i+j];
	m->atoms[i].atno=(int)znucl[typat[i]-1];
      }
      addabs(m->atoms,m->n,c->basis);
      free(xred);
      xred=NULL;
    }
    if (xcart){
      for(i=0;i<m->n;i++){
	for (j=0;j<3;j++)
	  m->atoms[i].abs[j]=xcart[3*i+j];
	m->atoms[i].atno=(int)znucl[typat[i]-1];
      }
      addfrac(m->atoms,m->n,c->recip);
      free(xcart);
      xcart=NULL;
    }
    free(typat);
    typat=NULL;
    free(znucl);
    znucl=NULL;
  }

  /* spin */

  if (spinat){
    if (!nsppol){ /* set default value */
      nsppol=malloc(sizeof(int));
      *nsppol=1;
    }
    if ((*nsppol==2)||((*nsppol==1)&&(nspden)&&(*nspden==2))){
      for(i=0;i<m->n;i++)
	m->atoms[i].spin=spinat[3*i+2];
      free(spinat);
    }
    else fprintf(stderr,"Discarding spins\n");
  }

  /* symmetry */

  if ((spgroup)&&(*spgroup)){
    hall=igr2hall[*spgroup];
    if (spgaxor){
      i=*spgroup;
      if ((i==146)||(i==148)||(i==155)||(i==160)||(i==161)||
          (i==166)||(i==167)){
        if (*spgaxor==2) hall+=1;
        else if (*spgaxor!=1) error_exit("Invalid value for spgaxor");
      }
      else if (*spgaxor!=1)
        error_exit("Value of spgaxor not supported by c2x");
    }
    if (spgorig){
      if ((*spgorig!=1)&&(*spgorig!=2))
        error_exit("Invalid value for spgorig");
      if (!spgr_is_double(*spgroup))
        if (*spgorig!=1)
          error_exit("Second origin requested in single origin system");
      if (*spgorig==2) hall++;
    }
    
    if (hall>1){

      if (debug) fprintf(stderr,"Using Hall number of %d\n",hall);
      cspg_hall2sym(hall,c,s);

      if (debug) fprintf(stderr,"%d symops returned\n",s->n);
      sym_expand(c,m,s);
      if ((natom)&&(m->n!=*natom))
        fprintf(stderr,"WARNING: found %d atoms after symmetrisation, "
                "expected %d\n",m->n,*natom);
    }
  }
  else if ((nsym)&&(*nsym>1)){
    if (!symrel) error_exit("nsym>1 but symrel not found");
    s->n=*nsym;
    s->ops=malloc(s->n*sizeof(struct sym_op));
    if (!s->ops) error_exit("Malloc error for symmetry ops");
    for(i=0;i<s->n;i++){
      for(j=0;j<3;j++)
	for(jj=0;jj<3;jj++)
	  mat[j][jj]=symrel[9*i+3*j+jj];
      mat_f2a(mat,s->ops[i].mat,c->basis,c->recip);
      if (tnons){
	if ((tnons[3*i]!=0)||(tnons[3*i+1]!=0)||(tnons[3*i+2]!=0)){
	  s->ops[i].tr=malloc(3*sizeof(double));
	  for(j=0;j<3;j++)
	    s->ops[i].tr[j]=tnons[3*i+j];
	}
      }
    }
    sym_expand(c,m,s);
    if ((natom)&&(m->n!=*natom))
      fprintf(stderr,"WARNING: found %d atoms after symmetrisation, "
	      "expected %d\n",m->n,*natom);
  }
  
  /* k-points */

  if ((kptopt)&&(*kptopt==0)){
    if (!nkpt){
      nkpt=malloc(sizeof(int));
      *nkpt=1;
    }
  }
      
  
  if (kpt){
    k->n=*nkpt;
    k->kpts=malloc(*nkpt*sizeof(struct atom));
    if (!k->kpts) error_exit("malloc error for kpts");
    if (wtk){ /* may need to normalise */
      x=0;
      for(i=0;i<*nkpt;i++)
	x+=wtk[i];
      for(i=0;i<*nkpt;i++)
	wtk[i]/=x;
    }
    /* Also may need to scale positions */
    if (kptnrm) x=1/(*kptnrm);
    else x=1;
    
    for(i=0;i<*nkpt;i++){
      for(j=0;j<3;j++)
	k->kpts[i].frac[j]=kpt[3*i+j]*x;
      if (wtk)
	k->kpts[i].wt=wtk[i];
      else
	k->kpts[i].wt=1.0/(*nkpt);
    }
    free(wtk);
  }
  else if ((ngkpt)&&((nshiftk==NULL)||(*nshiftk==1))){
    k->mp=malloc(sizeof(struct mp_grid));
    if (!k->mp) error_exit("malloc error for MP grid");
    for(i=0;i<3;i++){
      k->mp->grid[i]=ngkpt[i];
      /* Convert to Castep's convention for shifts */
      if (nshiftk)
        k->mp->disp[i]=shiftk[i];
      else
        k->mp->disp[i]=0;
      if ((ngkpt[i]&1)==0) /* is grid even? */
        k->mp->disp[i]-=0.5;
      k->mp->disp[i]/=ngkpt[i];
    }
    if (nshiftk) {free(shiftk); shiftk=NULL;}
    if (ngkpt) {free(ngkpt); ngkpt=NULL;}
  }
  else if ((nkpt)&&(*nkpt==1)){
    k->n=1;
    k->kpts=malloc(sizeof(struct atom));
    if (!k->kpts) error_exit("malloc error for kpts");
    for(i=0;i<3;i++)
      k->kpts[0].frac[i]=0;
    k->kpts[0].wt=1.0;
  }
  else if ((ngkpt)&&(nshiftk)&&(*nshiftk!=1)){
    fprintf(stderr,"Warning, discarding k point grid:\n"
	    "  c2x does not currently support multiple shifted grids\n");
  }

  if (kpt) free(kpt);
  if (nkpt) free(nkpt);
  
  /* misc */

  if (ecut) {
    e->cut_off=*ecut;
    free(ecut);
  }
  if (toldfe) {
    e->etol=*toldfe/m->n;
    free(toldfe);
  }
  
}

static int ab_len_read(char *line, int off, double *data, int n,
		       FILE* infile){
  int i;
  off=ab_float_read(line,off,data,n,infile);
  /* Do we have some units specified? */
  while (line[off]==' ') off++;
  if (!strncasecmp(line+off,"angstr",6)){
    off+=6;
    while(isalpha(line[off])) off++;
    return off;
  }
  /* Default for lengths is Bohrs */
  for(i=0;i<n;i++)
    data[i]*=BOHR;

  return off;
  
}

static int ab_energy_read(char *line, int off, double *data, int n,
		       FILE* infile){
  int i;
  off=ab_float_read(line,off,data,n,infile);
  /* Do we have some units specified? */
  while (line[off]==' ') off++;
  if (!strncasecmp(line+off,"eV",2)){
    off+=2;
    return off;
  }
  else if (!strncasecmp(line+off,"Ry",2)){
    off+=2;
    for(i=0;i<n;i++)
      data[i]*=0.5;  /* convert to Hartrees, then to eV below */
  }
  /* Default for energy is Hartrees */
  for(i=0;i<n;i++)
    data[i]*=H_eV;

  return off;
  
}


/* This routine will read multiple floats, potentially spread across
 * multiple lines. It understands the "rep*data" syntax, and an
 * implicit missing rep. It can parse arithmetic expressions, understanding
 * more syntax than abinit */

static int ab_float_read(char *line, int off, double *data, int n,
			 FILE* infile){
  int i,j,tmp,rep,t2;
  double result;
  char *ptr;
  
  if (debug>3) fprintf(stderr,"ab_float_read for %d items\n",n);
  if (!data) error_exit("malloc error in abinit_in_read()");
  i=0;
  
  while(i<n){
    while((line[off])&&(isspace(line[off]))) off++;
    if (line[off]==0){
      ab_readline(line,&ptr,infile); /* will remove leading spaces */
      off=ptr-line;
    }
    /* check for repeat syntax */
    rep=1;
    if (line[off]=='*'){
      rep=n-i;
      off++;
    }
    else{
      sscanf(line+off,"%d%n",&tmp,&t2);
      if (line[off+t2]=='*'){
	rep=tmp;
	off+=t2+1;
      }
    }
    /* Find end of current token */
    tmp=off;
    while(line[tmp]&&(!isspace(line[tmp]))) tmp++;
    line[tmp]=0;

    if(!ascan(line+off,&result)){
      fprintf(stderr,"Failed to parse \"%s\" as a float\n",line+off);
      exit(1);
    }

    for(j=0;j<rep;j++) data[i++]=result;

    off=tmp+1;
	
  }

  if (debug>3){
    fprintf(stderr,"Read:");
    for(i=0;i<n;i++) fprintf(stderr," %f",data[i]);
    fprintf(stderr,"\n");
  }
  return off;
}


/* This routine will read multiple integers, potentially spread across
 * multiple lines. It understands the "rep*data" syntax, and an
 * implicit missing rep */

static int ab_int_read(char *line, int off, int *data, int n,
			 FILE* infile){
  int i,j,tmp,rep,dat,tmp1,tmp2;
  char *ptr;

  if (debug>3) fprintf(stderr,"ab_int_read for %d items\n",n);
  if (!data) error_exit("malloc error in abinit_in_read()");
  i=0;
  while(i<n){
    while ((j=sscanf(line+off,"%d%n*%d%n",&rep,&tmp1,&dat,&tmp2))){
      if (j==EOF) break;
      if (j==1){
	data[i++]=rep;
	off+=tmp1;
      }
      else{
	if (debug>1) fprintf(stderr,"Read %d * %d\n",rep,dat);
	if (rep>n-i) rep=n-i;
	for(j=0;j<rep;j++)
	  data[i++]=dat;
	off+=tmp2;
      }
      if (i==n) break;
    }
    if (i==n) break;
    while((*(line+off))&&(isspace(*(line+off)))) off++;
    if (*(line+off)==0){
      ab_readline(line,&ptr,infile);
      off=ptr-line;
    }
    else { /* pick up the missing rep "*dat" syntax here */
      if (sscanf(line+off,"*%d%n",&dat,&tmp)==1){
	for(;i<n;i++) data[i]=dat;
	off+=tmp;
	break;
      }
      else{
	fprintf(stderr,"Unexpected content in .in: %s\n",line+off);
	exit(1);
      }
    }
  }

  if (debug>3){
    fprintf(stderr,"Read:");
    for(i=0;i<n;i++) fprintf(stderr," %d",data[i]);
    fprintf(stderr,"\n");
  }
  return off;
}

static int ab_readline(char *buffer, char **p, FILE* infile){
  char *ptr;
  
  if (!fgets(buffer,LINE_SIZE+1,infile)) return 0;

  /* anything after a # or ! is a comment and should be removed */

  ptr=buffer;
  for(;*ptr;ptr++)
    if ((*ptr=='#')||(*ptr=='!')) {*ptr=0; break;}

  /* any = is treated as a space */

  ptr=buffer;
  for(;*ptr;ptr++)
    if (*ptr=='=') *ptr=' ';

  /* eat leading spaces */
  *p=buffer;
  while((**p)&&(isspace(**p))) (*p)++;

  if ((**p)==0) return ab_readline(buffer,p,infile);

  return 1;
  
}
