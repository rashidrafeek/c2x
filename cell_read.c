/* Read some useful data from a CASTEP or ONETEP .cell file
 * 
 * Should read lattice_cart, lattice_abc, positions_frac and positions_abs
 * blocks. Should skip blanks lines and comments, and should cope with
 * UNIX, DOS or Mac line-endings.
 */


/* Copyright (c) 2007-2020 MJ Rutter 
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

static int cellreadline(char *buffer, int len, char *dir,
                       struct infiles **filesp);
static int cell_read_length(char *buff, double *x);

static char *ptr;
static char **title;
static int is_onetep;
  
static void spin_read(char *buff, double *spin){
  while (*buff==' ') buff++;
  if (!strncasecmp(buff,"spin",4)){
    buff+=4;
    while ((*buff==' ')||(*buff=='=')||(*buff==':')) buff++;
    sscanf(buff,"%lf",spin);
  }
}


static int set_length_unit(char *ptr, double *unit){
  if (!strncasecmp(ptr,"ang",3)){
    *unit=1;
    return 1;
  }
  else if(!strncasecmp(ptr,"bohr",4)){
    *unit=BOHR;
    return 1;
  }
  else if(!strncasecmp(ptr,"a0",2)){
    *unit=BOHR;
    return 1;
  }
  else if(!strncasecmp(ptr,"nm",2)){
    *unit=10;
    return 1;
  }
  return 0;
}

void cell_read(FILE* in, struct unit_cell *c, struct contents *m,
               struct kpts *kp, struct symmetry *s){
  int i,j,k,n,n2,frac,nsym,success,nlabels;
  int smisc_size=0,sblock_size=0;
  double scale,default_scale;
  double lat_abc[6],*dptr;
  char sym[14],*ptr2;
  static char buffer[LINE_SIZE+1];
  double *sym_mat,*sym_disp;
  struct label {char *l; char *sym;} *labels;
  struct infiles *files;
  char *dir;

  dir=dict_get(m->dict,"in_dir");
  
  files=malloc(sizeof(struct infiles));
  if (!files) error_exit("Malloc error for files in cell_read");
  files->f=in;
  files->next=files->last=NULL;
  files->name=NULL;
  files->line=0;

  is_onetep=0;
  if (dict_get(m->dict,"cell_is_onetep")) is_onetep=1;
  default_scale=1;
  if (dict_get(m->dict,"cell_is_onetep")) default_scale=BOHR;
  scale=default_scale;
  m->n=0;
  nsym=0;
  frac=0;
  sym_mat=sym_disp=NULL;
  title=&(m->title);
  labels=NULL;
  nlabels=0;

  if (debug>2) fprintf(stderr,"Cell read called\n");

  if (!(c->basis=malloc(72))) error_exit("Malloc error in cellread for c->basis");

  while(cellreadline(buffer,LINE_SIZE,dir,&files)){
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
        (!strncasecmp(ptr,"kpoints_mp_grid",15))||
        (!strncasecmp(ptr,"kpoint_mp_offset",16))||
        (!strncasecmp(ptr,"kpoints_mp_offset",17))){
      /* Set ptr2 to char after string matched */
      ptr2=ptr+14;
      if ((ptr[6]=='s')||(ptr[6]=='S')) ptr2++;
      if ((*ptr2=='e')||(*ptr2=='E')) ptr2+=2;
      /* Increment ptr2 over spaces, = and : */
      while ((*ptr2)&&((*ptr2==' ')||(*ptr2==':')||(*ptr2=='='))) ptr2++;
      if (!kp->mp){
        kp->mp=malloc(sizeof(struct mp_grid));
        if (!kp->mp) error_exit("Malloc error for struct mp_grid!");
        for(i=0;i<3;i++) kp->mp->grid[i]=0;
        for(i=0;i<3;i++) kp->mp->disp[i]=0;
      }
      if((!strncasecmp(ptr,"kpoint_mp_grid",14))||
         (!strncasecmp(ptr,"kpoints_mp_grid",15))){
        if (sscanf(ptr2,"%d %d %d",kp->mp->grid,
                   kp->mp->grid+1,kp->mp->grid+2)!=3){
          fprintf(stderr,"Error parsing:\n%s\n",buffer);
          exit(1);
        }
      }else{
        if (multi_scan(ptr2,kp->mp->disp,3,NULL)!=3){
          fprintf(stderr,"Error parsing:\n%s\n",buffer);
          exit(1);
        }
      }
    }
    if ((!strncasecmp(ptr,"kpoint_mp_spacing",17))||
        (!strncasecmp(ptr,"kpoints_mp_spacing",18))){
      ptr2=ptr+17;
      if ((ptr[6]=='s')||(ptr[6]=='S')) ptr2++;
      /* Increment ptr2 over spaces, = and : */
      while ((*ptr2)&&((*ptr2==' ')||(*ptr2==':')||(*ptr2=='='))) ptr2++;
      dptr=malloc(sizeof(double));
      if (sscanf(ptr2,"%lf",dptr)!=1){
        fprintf(stderr,"Warning: error parsing kpoint_mp_spacing\n");
        free(dptr);
      }
      else kp->spacing=dptr;
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
      scale=default_scale;
      cellreadline(buffer,LINE_SIZE,dir,&files);
      ptr=buffer;
      if (set_length_unit(ptr,&scale))
        cellreadline(buffer,LINE_SIZE,dir,&files);
      for(i=0;i<3;i++){
        if(multi_scan(buffer,c->basis[i],3,NULL)!=3){
          PARSE_ERROR;
          if (debug) fprintf(stderr,"%s\n",buffer);
        }
        cellreadline(buffer,LINE_SIZE,dir,&files);
      }
      if (scale!=1)
        for(i=0;i<3;i++)
          for(j=0;j<3;j++)
            c->basis[i][j]*=scale;

    }else if(!strcasecmp(ptr,"lattice_abc")){
      scale=default_scale;
      cellreadline(buffer,LINE_SIZE,dir,&files);
      ptr=buffer;
      if (set_length_unit(ptr,&scale))
        cellreadline(buffer,LINE_SIZE,dir,&files);
      for(i=0;i<2;i++){
        if(multi_scan(buffer,lat_abc+3*i,3,NULL)!=3){
          PARSE_ERROR;
          exit(1);
        }
        cellreadline(buffer,LINE_SIZE,dir,&files);
      }
      if (scale!=1)
        for(i=0;i<3;i++)
          *(lat_abc+i)*=scale;
      abc2cart(lat_abc,c);

    }else if((!strcasecmp(ptr,"positions_frac"))||
            (!strcasecmp(ptr,"positions_abs"))){
      frac=1;
      if(!strcasecmp(ptr,"positions_abs")) frac=0;
      cellreadline(buffer,LINE_SIZE,dir,&files);
      ptr=buffer;
      scale=default_scale;
      if(!frac){
        if (set_length_unit(ptr,&scale))
          cellreadline(buffer,LINE_SIZE,dir,&files);
      }
      if (dict_get(m->dict,"cell_is_onetep")&&(!nlabels)){
        fprintf(stderr,"Warning: ONETEP input specified, but %%block species "
                "does not preceed atoms\n");
      }
      for(i=0;;i++){
	m->atoms=realloc(m->atoms,(i+1)*sizeof(struct atom));
	if (!m->atoms) error_exit("realloc error in cell_read");
	init_atoms(m->atoms+i,1);
        if (frac) dptr=m->atoms[i].frac;
        else dptr=m->atoms[i].abs;
        if (dict_get(m->dict,"cell_is_onetep")&&(nlabels)){
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
          if(!strncasecmp(ptr,"%endblock",9)) break;
          /* Try to read symbol / atomic number */
          sym[1]=sym[2]=sym[3]=0;
          if(sscanf(buffer,"%d%n",&m->atoms[i].atno,&n)>=1){
            ;
          }
          else if (sscanf(buffer,"%12s%n",sym,&n)>=1){
            /* Does symbol contain a colon? */
            /* If so, colon and everything after it is stored in the label,
             * and the part before should be an atomic symbol */
            if ((sym[1]==':')||(sym[2]==':')||(sym[3]==':')){
              ptr2=sym;
              while(*ptr2!=':') ptr2++;
              m->atoms[i].label=malloc(strlen(ptr2)+1);
              strcpy(m->atoms[i].label,ptr2);
              *ptr2=0;
            }
            else{
              m->atoms[i].label=NULL;
            }
            m->atoms[i].atno=atsym2no(sym);
          }
          else{
            PARSE_ERROR;
            if (debug) fprintf(stderr,"%s\n",ptr); 
            exit(1);
          }
          /* Now need to read co-ordinates */
          if(multi_scan(buffer+n,dptr,3,&n2)==3){
            n+=n2;
            spin_read(buffer+n,&m->atoms[i].spin);
          }
          else{
            PARSE_ERROR;
            if (debug) fprintf(stderr,"%s\n",ptr); 
            exit(1);
          }
        }
        cellreadline(buffer,LINE_SIZE,dir,&files);
      }
      m->n=i;
      if (debug>2) fprintf(stderr,"%d m->atoms read\n",m->n);
      if ((!frac)&&(scale!=1))
        for(i=0;i<m->n;i++)
          for(j=0;j<3;j++)
            m->atoms[i].abs[j]*=scale;
    }else if(!strcasecmp(ptr,"ionic_velocities")){
      if (!m->n) fprintf(stderr,
                         "ionic velocities preceed positions, ignoring\n");
      else{
        cellreadline(buffer,LINE_SIZE,dir,&files);
        ptr=buffer;
        scale=default_scale;
        if (set_length_unit(ptr,&scale)){
          while((*ptr)&&(*ptr!='/')) ptr++;
          if (*ptr!='/'){
            fprintf(stderr,"Warning: failed to parse units in %s\n",buffer);
            scale=default_scale;
          }
          else{
            ptr++;
            if (strcasecmp(ptr,"ps")==0)
              ;
            else if (strcasecmp(ptr,"ns")==0)
              scale*=1000;
            else if (strcasecmp(ptr,"fs")==0)
              scale/=1000;
            else {
              fprintf(stderr,"Warning: failed to parse units in %s\n",buffer);
              scale=default_scale;
            }
          }
          cellreadline(buffer,LINE_SIZE,dir,&files);
        }
        for(i=0;i<m->n;i++){
          if (sscanf(buffer,"%lf %lf %lf",m->atoms[i].v,m->atoms[i].v+1,
                     m->atoms[i].v+2)!=3)
            error_exit("Parse error for velocities");
          for(j=0;j<3;j++) m->atoms[i].v[j]*=scale;
          cellreadline(buffer,LINE_SIZE,dir,&files);
        }
        ptr=buffer;
        while(*ptr==' ') ptr++;
        if (strncasecmp(ptr,"%endblock",9))
          error_exit("%endblock ionic_velocities expected but not found");
        m->velocities=1;
      }
    }else if(!strcasecmp(ptr,"kpoints_list")){
      i=0;
      while(1){
        cellreadline(buffer,LINE_SIZE,dir,&files);
        if (!strncasecmp(buffer,"%endblock",9)) break;
        kp->kpts=realloc(kp->kpts,(i+1)*sizeof(struct atom));
        if (!kp->kpts) error_exit("realloc error for kpts");
        if (multi_scan(buffer,kp->kpts[i].frac,3,&j)!=3){
          PARSE_ERROR;
          exit(1);
        }
        if (single_scan(buffer+j,&(kp->kpts[i].wt),NULL)!=1){
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
        cellreadline(buffer,LINE_SIZE,dir,&files);
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
          cellreadline(buffer,LINE_SIZE,dir,&files);
        }
        if(sscanf(buffer,"%lf %lf %lf",sym_disp+3*nsym,
                 sym_disp+3*nsym+1,sym_disp+3*nsym+2)!=3){
          PARSE_ERROR;
          exit(1);
        }
      }
    }else if(!strcasecmp(ptr,"species")){
      if (!dict_get(m->dict,"cell_is_onetep"))
        fprintf(stderr,"Warning: %%block species found, "
                                   "but not parsing as Onetep input.\n"
                                   "Was -1 forgotten as commandline flag?\n");

      if (!cellreadline(buffer,LINE_SIZE,dir,&files))
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
     /*
      * sscanf(buffer,"%ms %ms",&(labels[nlabels].l),&(labels[nlabels].sym));
      */
        sscanfmsn(buffer,&(labels[nlabels].l),&i);
        sscanfmsn(buffer+i,&(labels[nlabels].sym),&i);
        nlabels++;
	if (!cellreadline(buffer,LINE_SIZE,dir,&files))
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
      m->species_misc=realloc(m->species_misc,smisc_size+strlen("%block ")+1);
      if (!m->species_misc) error_exit("Realloc error in a species_ block");
      strcpy(m->species_misc+smisc_size,"%block ");
      smisc_size+=strlen("%block ");

      m->species_misc=realloc(m->species_misc,smisc_size+strlen(ptr)+2);
      if (!m->species_misc) error_exit("Realloc error in a species_ block");
      strcpy(m->species_misc+smisc_size,ptr);
      smisc_size+=strlen(ptr);
      strcpy(m->species_misc+smisc_size,"\n");
      smisc_size+=1;

      while(strncasecmp(buffer,"%endblock",9)){
        /* If we use cellreadline here, we will strip comments */
        /* If we don't we will fail to parse includefile lines */
	if (!fgets(buffer,LINE_SIZE,files->f))
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
      while(cellreadline(buffer,LINE_SIZE,dir,&files)){
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

int cellreadline(char *buffer, int len, char *dir, struct infiles **filesp){
  int off,success2;
  char *ptr,*success,*p2;

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
    while(*ptr==' ') ptr++;
/* Parse a very special comment */
    if (!strncasecmp(ptr,"#titl ",6)){ /* We have a shelx-style title */
      ptr+=6;
      p2=ptr;
      while((*p2)&&(*p2!='\n')&(*p2!='\r')) p2++;
      *title=malloc(p2-ptr+1);
      if(!*title) error_exit("Malloc error in cell_read()");
      strncpy(*title,ptr,p2-ptr);
      (*title)[p2-ptr]=0;
      continue;
    }
/* Skip comments and blank lines */
    if (!strncasecmp(ptr,"comment",7)) continue;
    if ((*ptr=='#')||(*ptr=='!')||(*ptr==';')||(*ptr==0)) continue;
    break;
  }

  if (!success){ /* Need to recurse out of nested include files for Onetep */
    success2=0;
    while ((*filesp)->last){
      fclose((*filesp)->f);
      free((*filesp)->name);
      (*filesp)=(*filesp)->last;
      free((*filesp)->next);
      (*filesp)->next=NULL;
      success2=cellreadline(buffer,LINE_SIZE,dir,filesp);
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
    include_file(filesp,dir,ptr);
    return cellreadline(buffer,LINE_SIZE,dir,filesp);
  }
  
  return (1);
}

static int cell_read_length(char *buff, double *x){
  char *p;
  int n,i;  

  while((*buff)&&((*buff==' ')||(*buff==':')||(*buff=='='))) buff++;

  /*
   *  n=sscanf(buff,"%lf %ms",x,&p);
   */

  n=sscanf(buff,"%lf%n",x,&i);
  n+=sscanfmsn(buff+i,&p,NULL);
  
  if (n==0) return 0;

  if (n==1) {
    if (is_onetep) (*x)*=BOHR;
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

  if (!strncasecmp(p,"nm",2)){
    (*x)*=10;
    free(p);
    return 1;
  }

  
  fprintf(stderr,"Unexpected length unit: %s\n",p);
  free(p);
  return 0;

}


