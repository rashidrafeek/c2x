/* Read some useful data from a CASTEP .cell file
 * 
 * Should read lattice_cart, lattice_abc, positions_frac and positions_abs
 * blocks. Should skip blanks lines and comments, and should cope with
 * UNIX, DOS or Mac line-endings.
 */


/* Copyright (c) 2007, 2017 MJ Rutter 
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

#include "c2xsf.h"

#define LINE_SIZE 100

#define PARSE_ERROR  fprintf(stderr,"Parse error on line %d %s\n",files->line,(files->name)?files->name:"")


static int cellreadline(char *buffer, int len);
static int cell_read_length(char *buff, double *x);

static char *ptr;
static char **title;

void spin_read(char *buff, double *spin){
  while (*buff==' ') buff++;
  if (!strncasecmp(buff,"spin",4)){
    buff+=4;
    while ((*buff==' ')||(*buff=='=')||(*buff==':')) buff++;
    sscanf(buff,"%lf",spin);
  }
}

static struct infiles {
  FILE* f; struct infiles *next; struct infiles *last; char *name;
  int line;
} *files;
static FILE *infile;


void cell_read(FILE* in, struct unit_cell *c, struct contents *m,
               struct kpts *kp, struct symmetry *s){
  int units,i,j,k,n,frac,nsym,success,default_units,nlabels;
  int smisc_size=0,sblock_size=0;
  double lat_abc[6],*dptr;
  char sym[4],*ptr2;
  static char buffer[LINE_SIZE+1];
  double *sym_mat,*sym_disp;
  struct label {char *l; char *sym;} *labels;
  
  files=malloc(sizeof(struct infiles));
  if (!files) error_exit("Malloc error for files in cell_read");
  files->f=in;
  files->next=files->last=NULL;
  files->name=NULL;
  files->line=0;
  infile=in;
  
  default_units=0;
  if (flags&ONETEP) default_units=1;
  units=default_units;
  m->n=0;
  nsym=0;
  frac=0;
  sym_mat=sym_disp=NULL;
  title=&(m->title);
  labels=NULL;
  nlabels=0;

  if (debug>2) fprintf(stderr,"Cell read called\n");

  if (!(c->basis=malloc(72))) error_exit("Malloc error in cellread for c->basis");

  while(cellreadline(buffer,LINE_SIZE)){
    ptr=buffer;
    
    if (!strncasecmp(ptr,"symmetry_tol",12)){
      if (!s->tol){
        s->tol=malloc(sizeof(double));
        if (!s->tol) error_exit("Malloc error for single double!");
        if (!cell_read_length(ptr+12,s->tol)){
          fprintf(stderr,"Error parsing:\n%s\n",buffer);
          free(s->tol);
          s->tol=NULL;
        }
      }
      continue;
    }
    if (!strncasecmp(ptr,"symmetry_generate",17)){
      s->gen=malloc(sizeof(int));
      *(s->gen)=1;
      continue;
    }

    if ((!strncasecmp(ptr,"kpoint_mp_grid",14))||
        (!strncasecmp(ptr,"kpoint_mp_offset",16))){
      /* Set ptr2 to char after string matched */
      ptr2=ptr+14;
      if (*ptr2=='e') ptr2+=2;
      /* Increment ptr2 over spaces, = and : */
      while ((*ptr2)&&((*ptr2==' ')||(*ptr2==':')||(*ptr2=='='))) ptr2++;
      if (!kp->mp){
        kp->mp=malloc(sizeof(struct mp_grid));
        if (!kp->mp) error_exit("Malloc error for struct mp_grid!");
        for(i=0;i<3;i++) kp->mp->grid[i]=0;
        for(i=0;i<3;i++) kp->mp->disp[i]=0;
      }
      if(!strncasecmp(ptr,"kpoint_mp_grid",14)){
        if (sscanf(ptr2,"%d %d %d",kp->mp->grid,
                   kp->mp->grid+1,kp->mp->grid+2)!=3){
          fprintf(stderr,"Error parsing:\n%s\n",buffer);
          exit(1);
        }
      }else{
        if (sscanf(ptr2,"%lf %lf %lf",kp->mp->disp,kp->mp->disp+1,
                   kp->mp->disp+2)!=3){
          fprintf(stderr,"Error parsing:\n%s\n",buffer);
          exit(1);
        }
      }
    }

/* What remains ought to be a block or a keyword.  */
    if (strncasecmp(ptr,"%block",6)) continue;
    ptr+=6;
/* Eat leading spaces */
    while(*ptr==' ') ptr++;
/* Kill trailing spaces and things after ! */
    ptr2=ptr;
    while((*ptr2)&&(*ptr2!=' ')&&(*ptr2!='!')) ptr2++;
    *ptr2=0; /* Either it was a null, or it should be one... */

    if (debug>2) fprintf(stderr,"Found a %s block\n",ptr);

    if (!strcasecmp(ptr,"lattice_cart")){
      units=default_units;
      cellreadline(buffer,LINE_SIZE);
      ptr=buffer;
      if (!strncasecmp(ptr,"ang",3)){
        units=0;
        cellreadline(buffer,LINE_SIZE);
      }else if(!strncasecmp(ptr,"bohr",4)){
        units=1;
        cellreadline(buffer,LINE_SIZE);
      }
      for(i=0;i<3;i++){
        if(sscanf(buffer,"%lf %lf %lf",c->basis[i],c->basis[i]+1,c->basis[i]+2)!=3){
          PARSE_ERROR;
          if (debug) fprintf(stderr,"%s\n",buffer);
        }
        cellreadline(buffer,LINE_SIZE);
      }
      if (units==1)
        for(i=0;i<3;i++)
          for(j=0;j<3;j++)
            c->basis[i][j]*=BOHR;

    }else if(!strcasecmp(ptr,"lattice_abc")){
      units=default_units;
      cellreadline(buffer,LINE_SIZE);
      ptr=buffer;
      if (!strncasecmp(ptr,"ang",3)){
        units=0;
        cellreadline(buffer,LINE_SIZE);
      }else if(!strncasecmp(ptr,"bohr",4)){
        units=1;
        cellreadline(buffer,LINE_SIZE);
      }
      for(i=0;i<2;i++){
        if(sscanf(buffer,"%lf %lf %lf",lat_abc+3*i,lat_abc+3*i+1,
                  lat_abc+3*i+2)!=3){
          PARSE_ERROR;
          exit(1);
        }
        cellreadline(buffer,LINE_SIZE);
      }
      if (units==1)
        for(i=0;i<3;i++)
          *(lat_abc+i)*=BOHR;
      abc2cart(lat_abc,c);

    }else if((!strcasecmp(ptr,"positions_frac"))||
            (!strcasecmp(ptr,"positions_abs"))){
      frac=1;
      if(!strcasecmp(ptr,"positions_abs")) frac=0;
      cellreadline(buffer,LINE_SIZE);
      ptr=buffer;
      units=default_units;
      if(!frac){
        if (!strncasecmp(ptr,"ang",3)){
          units=0;
          cellreadline(buffer,LINE_SIZE);
        }else if(!strncasecmp(ptr,"bohr",4)){
          units=1;
          cellreadline(buffer,LINE_SIZE);
        }
      }
      if ((flags&ONETEP)&&(!nlabels)){
        fprintf(stderr,"Warning: ONETEP input specified, but %%block species "
                "does not preceed atoms\n");
      }
      for(i=0;;i++){
	m->atoms=realloc(m->atoms,(i+1)*sizeof(struct atom));
	if (!m->atoms) error_exit("realloc error in cell_read");
        m->atoms[i].spin=0;
        m->atoms[i].label=NULL;
        if (frac) dptr=m->atoms[i].frac;
        else dptr=m->atoms[i].abs;
        if ((flags&ONETEP)&&(nlabels)){
          if(!strncasecmp(ptr,"%endblock",9)) break;
          success=0;
          for(j=0;j<nlabels;j++){
            if((!strncasecmp(buffer,labels[j].l,strlen(labels[j].l)))&&
               (buffer[strlen(labels[j].l)]==' ')){
              success=1;
              m->atoms[i].atno=atsym2no(labels[j].sym);
              /* Record label only if is not an atomic symbol */
              if (strcmp(labels[j].l,atno2sym(m->atoms[i].atno))){
                m->atoms[i].label=malloc(strlen(labels[j].l)+1);
                if(!m->atoms[i].label) error_exit("Malloc error for label");
                strcpy(m->atoms[i].label,labels[j].l);
              }
              break;
            }
          }
          if (!success){
            fprintf(stderr,"Failed to find symbol for label\n");
            PARSE_ERROR;
            if (debug) fprintf(stderr,"%s\n",buffer);
            exit(1);
          }
          if(sscanf(buffer,"%*s %lf %lf %lf%n",dptr,
                    dptr+1,dptr+2,&n)>=3)
            spin_read(buffer+n,&m->atoms[i].spin);
          else{
            PARSE_ERROR;
            if (debug) fprintf(stderr,"%s\n",ptr); 
            exit(1);
          }            
        }
        else{ /* Not ONETEP */
          if(sscanf(buffer,"%d %lf %lf %lf%n",&m->atoms[i].atno,dptr,
                    dptr+1,dptr+2,&n)>=4)
            spin_read(buffer+n,&m->atoms[i].spin);
          else if(sscanf(buffer,"%3s %lf %lf %lf%n",sym,dptr,
                         dptr+1,dptr+2,&n)>=4){
            m->atoms[i].atno=atsym2no(sym);
            m->atoms[i].label=NULL;
            spin_read(buffer+n,&m->atoms[i].spin);
          }
          else if(!strncasecmp(ptr,"%endblock",9)) break;
          else{
            PARSE_ERROR;
            if (debug) fprintf(stderr,"%s\n",ptr); 
            exit(1);
          }
        }
        cellreadline(buffer,LINE_SIZE);
      }
      m->n=i;
      if (debug>2) fprintf(stderr,"%d m->atoms read\n",m->n);
      if (units==1)
        for(i=0;i<m->n;i++)
          for(j=0;j<3;j++)
            m->atoms[i].abs[j]*=BOHR;
    }else if(!strcasecmp(ptr,"kpoints_list")){
      i=0;
      while(1){
        cellreadline(buffer,LINE_SIZE);
        if (!strncasecmp(buffer,"%endblock",9)) break;
        kp->kpts=realloc(kp->kpts,(i+1)*sizeof(struct atom));
        if (!kp->kpts) error_exit("realloc error for kpts");
        if (sscanf(buffer,"%lf %lf %lf %lf",kp->kpts[i].frac,kp->kpts[i].frac+1,
                   kp->kpts[i].frac+2,&(kp->kpts[i].wt))!=4){
          PARSE_ERROR;
          exit(1);
        }
        i++;
      }
      kp->n=i;
    }else if(!strcasecmp(ptr,"symmetry_ops")){
      sym_mat=NULL;
      sym_disp=NULL;
      for(nsym=0;;nsym++){
        cellreadline(buffer,LINE_SIZE);
        if (!strncasecmp(buffer,"%endblock",9)) break;
        sym_mat=realloc(sym_mat,(nsym+1)*9*sizeof(double));
        sym_disp=realloc(sym_disp,(nsym+1)*3*sizeof(double));
	if ((!sym_mat)||(!sym_disp)) error_exit("realloc error in cell_read");
        for(i=0;i<3;i++){
          if(sscanf(buffer,"%lf %lf %lf",sym_mat+9*nsym+3*i,
                   sym_mat+9*nsym+3*i+1,sym_mat+9*nsym+3*i+2)!=3){
            PARSE_ERROR;
            exit(1);
          }
          cellreadline(buffer,LINE_SIZE);
        }
        if(sscanf(buffer,"%lf %lf %lf",sym_disp+3*nsym,
                 sym_disp+3*nsym+1,sym_disp+3*nsym+2)!=3){
          PARSE_ERROR;
          exit(1);
        }
      }
    }else if(!strcasecmp(ptr,"species")){
      if (!(flags&ONETEP)) fprintf(stderr,"Warning: %%block species found, "
                                   "but not parsing as Onetep input.\n"
                                   "Was -1 forgotten as commandline flag?\n");

      if (!cellreadline(buffer,LINE_SIZE))
        error_exit("Read error in species block");
      while(strncasecmp(buffer,"%endblock",9)){
        /* If we use cellreadline here, we will strip comments */
        /* If we don't we will fail to parse includefile lines */
        m->block_species=realloc(m->block_species,sblock_size+strlen(buffer)+2);
        if (!m->block_species) error_exit("Realloc error in species block");
        strcpy(m->block_species+sblock_size,buffer);
        sblock_size+=strlen(buffer);
        strcpy(m->block_species+sblock_size,"\n");
        sblock_size++;
        labels=realloc(labels,(nlabels+1)*sizeof(struct label));
        if (!labels) error_exit("Realloc error for labels");
        sscanf(buffer,"%ms %ms",&(labels[nlabels].l),&(labels[nlabels].sym));
        nlabels++;
	if (!cellreadline(buffer,LINE_SIZE))
	  error_exit("Read error in species block");
      }
    }
    else if((!strcasecmp(ptr,"species_pot"))||
	     (!strcasecmp(ptr,"species_lcao_states"))||
	     (!strcasecmp(ptr,"species_q"))||
	     (!strcasecmp(ptr,"species_mass"))||
             (!strcasecmp(ptr,"species_atomic_set"))||
             (!strcasecmp(ptr,"species_gamma"))||
             (!strcasecmp(ptr,"species_ngwf_plot"))){
      m->species_misc=realloc(m->species_misc,smisc_size+strlen("%BLOCK ")+1);
      if (!m->species_misc) error_exit("Realloc error in a species_ block");
      strcpy(m->species_misc+smisc_size,"%BLOCK ");
      smisc_size+=strlen("%BLOCK ");

      m->species_misc=realloc(m->species_misc,smisc_size+strlen(ptr)+2);
      if (!m->species_misc) error_exit("Realloc error in a species_ block");
      strcpy(m->species_misc+smisc_size,ptr);
      smisc_size+=strlen(ptr);
      strcpy(m->species_misc+smisc_size,"\n");
      smisc_size+=1;

      while(strncasecmp(buffer,"%endblock",9)){
        /* If we use cellreadline here, we will strip comments */
        /* If we don't we will fail to parse includefile lines */
	if (!fgets(buffer,LINE_SIZE,infile))
	  error_exit("Read error in a species_ block");
	files->line++;
        m->species_misc=realloc(m->species_misc,smisc_size+strlen(buffer)+1);
        if (!m->species_misc) error_exit("Realloc error in a species_ block");
        strcpy(m->species_misc+smisc_size,buffer);
        smisc_size+=strlen(buffer);
      }
      m->species_misc=realloc(m->species_misc,smisc_size+2);
      if (!m->species_misc) error_exit("Realloc error in a species_ block");
      strcpy(m->species_misc+smisc_size,"\n");
      smisc_size+=1;
    }else{
      while(cellreadline(buffer,LINE_SIZE)){
        ptr=buffer;
        if (!strncasecmp(ptr,"%endblock",9)) break;
      }
    }
  }

  real2rec(c);
  if (frac) addabs(m->atoms,m->n,c->basis);
  else addfrac(m->atoms,m->n,c->recip);

  /* Rescue symmetry operations into sane structure */
  if ((sym_mat)&&(nsym>1)){
    s->n=nsym;
    s->ops=malloc(nsym*sizeof(struct sym_op));
    if (!s->ops) error_exit("Malloc error for symops");
    for(i=0;i<nsym;i++){
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
    if (debug) fprintf(stderr,"%d symmetry operations read\n",nsym);
  }

  if ((sym_mat)&&(nsym<=1)){
    free(sym_disp);
    free(sym_mat);
    nsym=0;
  }

  free(labels);

}

int cellreadline(char *buffer, int len){
  int off,success2;
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
/* Parse a very special comment */
    if (!strncasecmp(ptr,"#titl ",6)){ /* We have a shelx-style title */
      ptr+=6;
      p2=ptr;
      while((*p2)&&(*p2!='\n')&(*p2!='\r')) p2++;
      *title=malloc(p2-ptr+1);
      if(!*title) error_exit("Malloc error in cell_read()");
      strncpy(*title,ptr,p2-ptr);
      *title[p2-ptr]=0;
      continue;
    }
/* Skip comments and blank lines */
    if (!strncasecmp(ptr,"comment",7)) continue;
    if ((*ptr=='#')||(*ptr=='!')||(*ptr==';')||(*ptr==0)) continue;
    break;
  }

  if (!success){ /* Need to recurse out of nested include files for Onetep */
    success2=0;
    while (files->last){
      fclose(infile);
      free(files->name);
      files=files->last;
      free(files->next);
      files->next=NULL;
      infile=files->f;
      success2=cellreadline(buffer,LINE_SIZE);
      if (success2) return success2;
    }
    if (!success2) return(0);
  }

  off=ptr-buffer;
  if (off){
    while(*ptr) {*(ptr-off)=*ptr;ptr++;}
    *(ptr-off)=0;
  }

  if (!strncasecmp(ptr,"includefile",11)){
    ptr=buffer+11;
    while ((*ptr==' ')||(*ptr==':')) ptr++;
    if (*ptr=='"'){
      ptr++;
      p2=ptr;
      while ((*p2)&&(*p2!='"')) p2++;
      if (*p2!='"') error_exit("Closing quote missing in cell_read");
      *p2=0;
    }
    files->next=malloc(sizeof(struct infiles));
    if (!files) error_exit("Malloc error for files_next in cell_read");
    files->next->last=files;
    files->next->next=NULL;
    files->next->f=fopen(ptr,"r");
    files=files->next;
    files->name=malloc(strlen(ptr+1));
    strcpy(files->name,ptr);
    files->line=0;
    infile=files->f;
    if (!files->f){
      fprintf(stderr,"Error, unable to open %s in cell_read\n",ptr);
      exit(1);
    }
    if (debug) fprintf(stderr,"Opened included file %s\n",ptr);
    return cellreadline(buffer,LINE_SIZE);
  }
  
  return (1);
}

static int cell_read_length(char *buff, double *x){
  char *p;
  int n;  

  while((*buff)&&((*buff==' ')||(*buff==':')||(*buff=='='))) buff++;
  
  n=sscanf(buff,"%lf %ms",x,&p);
  
  if (n==0) return 0;

  if (n==1) {
    if (flags&ONETEP) (*x)*=BOHR;
    return 1;
  }

  if (!strncasecmp(p,"bohr",4)){
    (*x)*=BOHR;
    free(p);
    return 1;
  }

  if (!strncasecmp(p,"ang",3)){
    free(p);
    return 1;
  }

  fprintf(stderr,"Unexpected length unit: %s\n",p);
  free(p);
  return 0;

}


