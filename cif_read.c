/* A very simplistic reader which, whilst not an (mm)cif reader,
* sometimes reads very basic versions of those files */

/* Copyright (c) 2014-2019 MJ Rutter 
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
#include<errno.h>
#include<ctype.h>

#include "c2xsf.h"

#define LINE_SIZE 2049

static char* cif_loop(FILE *infile, struct unit_cell *c, struct contents *m,
		      struct symmetry *s, char ***sxyz, int *need_abs);

int scmp(char *s1, char *s2){
  return (strncasecmp(s1,s2,strlen(s2))==0)?1:0;
}

void cif_read(FILE* infile, struct unit_cell *c, struct contents *m,
	      struct symmetry *s){
  int have_data=0,have_basis=0,hit,i,no_read,sym_warn=0,need_abs=0;
  double abc[6];
  char buffer[LINE_SIZE+1], *buff;
  char **sym_as_xyz;

  for(i=0;i<6;i++) abc[i]=123456;
  m->n=0;
  no_read=0;
  sym_as_xyz=NULL;
  
  if (debug>2) fprintf(stderr,"cif_read called\n");

  if (!(c->basis=malloc(72)))
    error_exit("Malloc error in cif_read for c->basis");

  while(1){
    if (no_read) no_read=0;
    else{
      for(i=0;i<LINE_SIZE-1;i++) buffer[i]=0;
      if (!fgets(buffer,LINE_SIZE,infile)) break;
    }

    buff=buffer;
    while (*buff==' ') buff++;
    if(buff[0]=='\n') continue;
    if(buff[0]=='\r') continue;
    if(buff[0]=='#') continue;

    if (!have_data){
      if (scmp(buff,"data_")) {have_data=1; continue;}
      else error_exit("data_ not found");
    }

    if (scmp(buff,"_struct_title ")||scmp(buff,"_struct.title ")){
      m->title=malloc(strlen(buff+strlen("_struct_title")));
      strcpy(m->title,buff+strlen("_struct_title")+1);
    }
    else if (scmp(buff,"_cell_length_a ")||scmp(buff,"_cell.length_a "))
      sscanf(buff+strlen("_cell_length_a "),"%lf",abc);
    else if (scmp(buff,"_cell_length_b ")||scmp(buff,"_cell.length_b "))
      sscanf(buff+strlen("_cell_length_b "),"%lf",abc+1);
    else if (scmp(buff,"_cell_length_c ")||scmp(buff,"_cell.length_c "))
      sscanf(buff+strlen("_cell_length_c "),"%lf",abc+2);
    else if (scmp(buff,"_cell_angle_alpha ")||
             scmp(buff,"_cell.angle_alpha "))
      sscanf(buff+strlen("_cell_angle_alpha "),"%lf",abc+3);
    else if (scmp(buff,"_cell_angle_beta ")||
             scmp(buff,"_cell.angle_beta "))
      sscanf(buff+strlen("_cell_angle_beta "),"%lf",abc+4);
    else if (scmp(buff,"_cell_angle_gamma ")||
             scmp(buff,"_cell.angle_gamma "))
      sscanf(buff+strlen("_cell_angle_gamma "),"%lf",abc+5);
    else if (scmp(buff,"_symmetry"))
      sym_warn=1;
    else if (scmp(buff,"loop_")){
      buff=cif_loop(infile,c,m,s,&sym_as_xyz,&need_abs);
      if (buff){
	strncpy(buffer,buff,LINE_SIZE);
	no_read=1;
      }
    }
    else if (debug>1)
      fprintf(stderr,"Ignoring %s",buff);
    /* fill in basis as soon as possible */
    if (!have_basis){
      hit=1;
      for(i=0;i<6;i++)
	if (abc[i]==123456) hit=0;
      if (hit){
	abc2cart(abc,c);
	have_basis=1;
        if(debug>2) fprintf(stderr,"Basis found in cif_read\n");
      }
    }
  }

  if (!have_basis) error_exit("Incomplete basis in cif_read");

  if (m->n==0) error_exit("No atoms found");

  real2rec(c);
  if (need_abs==1)
    addabs(m->atoms,m->n,c->basis);
  if (need_abs==-1)
    addfrac(m->atoms,m->n,c->recip);
  
  if ((sym_warn)&&(s->n==0)){
    fprintf(stderr,"Warning: CIF file contains symmetry information, but\n"
	    "no list of symmetry operations. Symmetry information ignored.\n");
  }

  /* Need to wait until basis read for this */

  if (sym_as_xyz){
    s->ops=malloc((s->n+1)*sizeof(struct sym_op));
    if(!s->ops) error_exit("realloc error in cif_read");
    for(i=0;i<s->n;i++){
      s->ops[i].tr=NULL;
      equiv_read(s->ops+i,c,sym_as_xyz[i]);
      free(sym_as_xyz[i]);
    }
    free(sym_as_xyz);
  }

  if (s->n>0){ /* Need to symmetrise atoms */
    sym_expand(c,m,s);
  }
}

void sym_expand(struct unit_cell *c, struct contents *m, struct symmetry *s){
  int i,j,k,sn,hit;
  double *wt;
  struct atom b,*satoms;
  
  if (debug>1) fprintf(stderr,"%d atoms before symmetrisation\n",m->n);
  reduce_cell_tol(m->atoms,m->n,c->basis,tol);
  wt=malloc(m->n*s->n*sizeof(double)); /* We cannot end up with more atoms
                                          than this */
  if (!wt) error_exit("malloc error in sym_expand");
  satoms=NULL;
  sn=0;

  for(i=0;i<m->n;i++){
    for(j=0;j<s->n;j++){
      sym_atom(m->atoms+i,&b,s->ops+j,c->recip);
      hit=atom_in_list(&b,satoms,sn,c->basis);
      if (hit==-1){
        satoms=realloc(satoms,(sn+1)*sizeof(struct atom));
        if(!satoms) error_exit("realloc error in sym_expand");
        addabs(&b,1,c->basis);
        satoms[sn]=b;
        wt[sn]=1;
        sn++;
      }
      else{ /* This deals, imperfectly, with low accuracy initial coords */
        for(k=0;k<3;k++){
          satoms[hit].frac[k]=(wt[hit]*satoms[hit].frac[k]+b.frac[k])/
            (wt[hit]+1);
        }
        wt[hit]+=1;
        addabs(satoms+hit,1,c->basis);
      }
    }
  }
  free(m->atoms);
  m->atoms=satoms;
  m->n=sn;
  if (debug>1) fprintf(stderr,"%d atoms after symmetrisation\n",m->n);
  free(wt);
  sort_atoms(m,1);
}


static char* cif_loop(FILE *infile, struct unit_cell *c, struct contents *m,
		      struct symmetry *s, char ***sxyz, int *need_abs){
  int i,j,maxfield,first;
  int fx,fy,fz,fsym,flab,fsym_as_xyz,fst_chg,*f;
  int cx,cy,cz;
  static char buffer[LINE_SIZE];
  char *p1,sym[4],*buff;
  char **sym_as_xyz;
  
  fx=fy=fz=fsym=flab=fsym_as_xyz=fst_chg=cx=cy=cz=-1;
  sym_as_xyz=*sxyz;
  f=NULL;
  i=0;
  if (debug>2) fprintf(stderr,"loop_ being parsed\n");

  while(1){
    for(j=0;j<LINE_SIZE-1;j++) buffer[j]=0;
    if (!fgets(buffer,LINE_SIZE,infile)) break;

    buff=buffer;
    while (*buff==' ') buff++;
    if(buff[0]=='\n') continue;
    if(buff[0]=='\r') continue;
    if(buff[0]=='#') continue;

    if (scmp(buff,"_atom_site_type_symbol")||
	scmp(buff,"_atom_site.type_symbol"))
      fsym=i;
    else if (scmp(buff,"_atom_site_label_atom_id")||
             scmp(buff,"_atom_site.label_atom_id"))
      flab=i;
    else if (scmp(buff,"_atom_site_charge")||
             scmp(buff,"_atom_site.charge")){
      fst_chg=i;
      dict_add(m->dict,"site_charge",(void*)1);
    }
    else if (scmp(buff,"_atom_site_label")||
	scmp(buff,"_atom_site.label"))
      {if (isspace(buff[strlen("_atom_site.label")])) flab=i;}
    else if (scmp(buff,"_atom_site_fract_x")||
             scmp(buff,"_atom_site.fract_x"))
      fx=i;
    else if (scmp(buff,"_atom_site_fract_y")||
	     scmp(buff,"_atom_site.fract_y"))
      fy=i;
    else if (scmp(buff,"_atom_site_fract_z")||
	     scmp(buff,"_atom_site.fract_z"))
      fz=i;
    else if (scmp(buff,"_atom_site_cartn_x")||
             scmp(buff,"_atom_site.cartn_x"))
      cx=i;
    else if (scmp(buff,"_atom_site_cartn_y")||
	     scmp(buff,"_atom_site.cartn_y"))
      cy=i;
    else if (scmp(buff,"_atom_site_cartn_z")||
	     scmp(buff,"_atom_site.cartn_z"))
      cz=i;
    else if (scmp(buff,"_symmetry_equiv_pos_as_xyz")||
	     scmp(buff,"_symmetry_equiv.pos_as_xyz"))
      fsym_as_xyz=i;
    else if (buff[0]!='_') break;
    else if (debug>1)
      fprintf(stderr,"In loop_ ignoring %s",buff);

    i++;
  }

  if (debug>2) fprintf(stderr,"fx=%d fy=%d fz=%d fsym=%d flab=%d "
		       "fsym_as_xyz=%d fst_chg=%d i=%d\n",
		       fx,fy,fz,fsym,flab,fsym_as_xyz,fst_chg,i);

  if (fsym==-1) fsym=flab;

  maxfield=-1;

  if (fsym_as_xyz!=-1){
    maxfield=fsym_as_xyz;
  }
  else{
    if (fx>maxfield) maxfield=fx;
    if (fy>maxfield) maxfield=fy;
    if (fz>maxfield) maxfield=fz;
    if (cx>maxfield) maxfield=cx;
    if (cy>maxfield) maxfield=cy;
    if (cz>maxfield) maxfield=cz;
    if (fsym>maxfield) maxfield=fsym;
    if (flab>maxfield) maxfield=flab;
    if (fst_chg>maxfield) maxfield=fst_chg;
    /* We must either have all of the below, or none */
    if (maxfield!=-1){
      if ((fx==-1)&&(cx==-1))
        error_exit("neither _atom_site frac_x nor cartn_x found");
      if ((fy==-1)&&(cy==-1))
        error_exit("neither _atom_site frac_y nor cartn_y found");
      if ((fz==-1)&&(cz==-1))
        error_exit("neither _atom_site frac_z nor cartn_z found");
      if (fsym==-1)
	error_exit("Neither _atom_site type_symbol nor label found");
      if ((cx==-1)||(cy==-1)||(cz==-1)){
        if ((fx!=-1)&&(fy!=-1)&&(fz!=-1))
          *need_abs=1;
        else
          error_exit("Mixed fractional and cartesian atomic co-ords");
      }
      if ((fx==-1)||(fy==-1)||(fz==-1)){
        if ((cx!=-1)&&(cy!=-1)&&(cz!=-1))
          *need_abs=-1;
        else
          error_exit("Mixed fractional and cartesian atomic co-ords");
      }
    }
  }

  f=malloc((maxfield+1)*sizeof(int));
  first=1;

  while(1){
    if (first) first=0;
    else{
      for(i=0;i<LINE_SIZE-1;i++) buffer[i]=0;
      if (!fgets(buffer,LINE_SIZE,infile)){
	free(f);
	return(NULL);
      }
    }

    buff=buffer;
    while (*buff==' ') buff++;
    if (*buff=='\n') continue;
    if (*buff=='\r') continue;
    if (*buff=='#') continue;
    if (*buff=='_') break;
    if (scmp(buff,"loop_")) break;

    if (maxfield==-1) continue;

    p1=buff;
    for(i=0;i<=maxfield;i++){
      while(*p1==' ') p1++;
      f[i]=p1-buff;
      if (*p1=='\''){
	p1++;
	while((*p1!='\'')&&(*p1!=0)) p1++;
      }
      if (*p1=='"'){
	p1++;
	while((*p1!='"')&&(*p1!=0)) p1++;
      }
      else while((*p1!=' ')&&(*p1!=0)) p1++;
      if (*p1==0) break;
    }
    if (i<maxfield) {
      fprintf(stderr,"short of fields in _loop, %d found, %d expected\n",
              i,maxfield);
      exit(1);
    }

    if (fsym_as_xyz!=-1){
      sym_as_xyz=realloc(sym_as_xyz,(s->n+1)*sizeof(char*));
      if(!sym_as_xyz) error_exit("realloc error in cif_read");
      *sxyz=sym_as_xyz;
      sym_as_xyz[s->n]=malloc(strlen(buff+f[fsym_as_xyz])+1);
      if(!sym_as_xyz[s->n]) error_exit("realloc error in cif_read");
      strcpy(sym_as_xyz[s->n],buff+f[fsym_as_xyz]);
      s->n++;
      continue;
    }    
    
    m->atoms=realloc(m->atoms,(m->n+1)*sizeof(struct atom));
    if(!m->atoms) error_exit("realloc error in cif_read");
    init_atoms(m->atoms+m->n,1);
    
    if (fx!=-1) sscanf(buff+f[fx],"%lf",m->atoms[m->n].frac);
    if (fy!=-1) sscanf(buff+f[fy],"%lf",m->atoms[m->n].frac+1);
    if (fz!=-1) sscanf(buff+f[fz],"%lf",m->atoms[m->n].frac+2);
    if (cx!=-1) sscanf(buff+f[cx],"%lf",m->atoms[m->n].abs);
    if (cy!=-1) sscanf(buff+f[cy],"%lf",m->atoms[m->n].abs+1);
    if (cz!=-1) sscanf(buff+f[cz],"%lf",m->atoms[m->n].abs+2);

    if (fst_chg!=-1) sscanf(buff+f[fst_chg],"%lf",&(m->atoms[m->n].site_chg));

    
    /* The symbol field ends on a space, the symbol part of the label
       field ends on a non-alpha, so we can read the symbol in the same
       manner from either */
    i=0;
    while(isalpha(buff[f[fsym]+i])&&(i<3)){
      sym[i]=buff[f[fsym]+i];
      i++;
    }
    sym[i]=0;
    m->atoms[m->n].atno=atsym2no(sym);

    m->n++;
  }

  /* Cannot add fractional / absolute co-ords, as might not have basis yet */

  *sxyz=sym_as_xyz;
  
  free(f);
  return(buffer);

}
