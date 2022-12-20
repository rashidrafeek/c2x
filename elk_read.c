/* Copyright (c) 2020 MJ Rutter 
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
#include<math.h>
#include<ctype.h>

#include "c2xsf.h"

#define LINE_SIZE 132 

static int elk_readline(char *buffer, int len, FILE* infile);

void elk_read(FILE* infile, struct unit_cell *c, struct contents *m,
              struct kpts *k, struct es *e){
  char line[LINE_SIZE+1], *ptr, *p2;
  double scale[3],x;
  int i,j,n,nspec,atno,natoms;
  int *ngridk;
  double *vkloff;

  ngridk=NULL;
  vkloff=NULL;
  
  /* Elk's natural units are Bohr */
  for(i=0;i<3;i++) scale[i]=BOHR;

  ptr=line;
  *ptr=0;
  
  while(1){
    if (!elk_readline(line,LINE_SIZE,infile)) break;

    if (!strcasecmp(line,"avec")){
      c->basis=malloc(9*sizeof(double));
      if (!c->basis) error_exit("malloc error for basis");
      for(i=0;i<3;i++){
	if (!elk_readline(line,LINE_SIZE,infile))
	  error_exit("EOF reading basis");
	n=sscanf(line,"%lf %lf %lf",c->basis[i],c->basis[i]+1,c->basis[i]+2);
	if (n!=3) error_exit("error parsing basis");
      }
    }
    else if (!strcasecmp(line,"scale")){
      if (!elk_readline(line,LINE_SIZE,infile)) error_exit("Unexpected EOF");
      n=sscanf(line,"%lf",&x);
      if (n!=1) error_exit("error parsing scale");
      for(i=0;i<3;i++)
	scale[i]*=x;
    }
    else if (!strcasecmp(line,"scale1")){
      if (!elk_readline(line,LINE_SIZE,infile)) error_exit("Unexpected EOF");
      n=sscanf(line,"%lf",&x);
      if (n!=1) error_exit("error parsing scale");
      scale[0]*=x;
    }
    else if (!strcasecmp(line,"scale2")){
      if (!elk_readline(line,LINE_SIZE,infile)) error_exit("Unexpected EOF");
      n=sscanf(line,"%lf",&x);
      if (n!=1) error_exit("error parsing scale");
      scale[1]*=x;
    }
    else if (!strcasecmp(line,"scale3")){
      if (!elk_readline(line,LINE_SIZE,infile)) error_exit("Unexpected EOF");
      n=sscanf(line,"%lf",&x);
      if (n!=1) error_exit("error parsing scale");
      scale[2]*=x;
    }
    else if (!strcasecmp(line,"atoms")){
      if (!elk_readline(line,LINE_SIZE,infile)) error_exit("Unexpected EOF");
      n=sscanf(line,"%d",&nspec);
      if (n!=1) error_exit("error parsing nspec");
      for(i=0;i<nspec;i++){
	if (!elk_readline(line,LINE_SIZE,infile)) error_exit("Unexpected EOF");
	p2=strrchr(line,'\'');
	if (!p2) error_exit("failed to find end quote for spfname");
	if (strncmp(p2-3,".in",3)) error_exit("unexpected spfname");
	p2-=3;
	*p2=0;
	p2=strrchr(line,'/');
	if (!p2) p2=strrchr(line,'\'');
	if (!p2) error_exit("failed to find quote for spfname");
	p2++;
	atno=atsym2no(p2);
	if (!elk_readline(line,LINE_SIZE,infile)) error_exit("Unexpected EOF");
	n=sscanf(line,"%d",&natoms);
	if (n!=1) error_exit("error parsing natoms");
	for(j=0;j<natoms;j++){
	  m->atoms=realloc(m->atoms,(m->n+1)*sizeof(struct atom));
	  if (!m->atoms) error_exit("realloc error for atoms");
	  init_atoms(m->atoms+m->n,1);
	  if (!elk_readline(line,LINE_SIZE,infile))
	    error_exit("Unexpected EOF");
	  n=sscanf(line,"%lf %lf %lf",m->atoms[m->n].frac,
		   m->atoms[m->n].frac+1,m->atoms[m->n].frac+2);
	  if (n!=3) error_exit("error parsing atomic coords");
	  m->atoms[m->n].atno=atno;
	  m->n++;
	}
      }
    }
    else if (!strcasecmp(line,"ngridk")){
      ngridk=malloc(3*sizeof(int));
      if (!ngridk) error_exit("malloc error for ngridk");
      if (!elk_readline(line,LINE_SIZE,infile))
	error_exit("EOF reading ngridk");
      n=sscanf(line,"%d %d %d",ngridk,ngridk+1,ngridk+2);
      if (n!=3) error_exit("error parsing ngridk");
    }
    else if (!strcasecmp(line,"vkloff")){
      vkloff=malloc(3*sizeof(double));
      if (!vkloff) error_exit("malloc error for vkloff");
      if (!elk_readline(line,LINE_SIZE,infile))
	error_exit("EOF reading vkloff");
      n=sscanf(line,"%lf %lf %lf",vkloff,vkloff+1,vkloff+2);
      if (n!=3) error_exit("error parsing vkloff");
    }
    else if (!strcasecmp(line,"sppath")){
      if (!elk_readline(line,LINE_SIZE,infile))
	error_exit("EOF reading sppath");
      dict_strcat(m->dict,"Elk_sppath",line);
    }
  }

  if (!c->basis) error_exit("No basis found");
  
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[i][j]*=scale[i];

  real2rec(c);
  addabs(m->atoms,m->n,c->basis);

  if (ngridk){
    k->mp=malloc(sizeof(struct mp_grid));
    for(i=0;i<3;i++)
      k->mp->grid[i]=ngridk[i];
    for(i=0;i<3;i++)
      if ((ngridk[i]&1)==1)
	k->mp->disp[i]=0;
      else
	k->mp->disp[i]=0.5/ngridk[i];
    if (vkloff){
      for(i=0;i<3;i++){
	k->mp->disp[i]=fmod(k->mp->disp[i]+vkloff[i]/ngridk[i],1.0/ngridk[i]);
	if (aeq(k->mp->disp[i],1.0/ngridk[i])) k->mp->disp[i]=0;
      }
      free(vkloff);
    }
    free(ngridk);
  }

}

void elk3d_read(FILE* infile, struct unit_cell *c, struct contents *m,
		struct kpts *kp, struct grid *gptr, struct es *e){
  char line[LINE_SIZE+1],*filename,*ptr;
  int i,j,k,ii,okay;
  int ngx,ngy,ngz;
  double x,v[3],basis[3][3],scale;
  FILE *f;

  elk_readline(line,LINE_SIZE,infile);

  if (sscanf(line,"%d %d %d",&ngx,&ngy,&ngz)!=3)
    error_exit("failed to read grid from Elk 3D file");

  gptr=grid_new(gptr);

  gptr->size[0]=ngx;
  gptr->size[1]=ngy;
  gptr->size[2]=ngz;
  if (!(gptr->data=malloc(ngx*ngy*ngz*sizeof(double))))
    error_exit("Malloc error for grid data");

  if (!c->basis){ /* We'd like a basis and some positions too */
    if (dict_get(m->dict,"in_dir")){
      filename=malloc(strlen(dict_get(m->dict,"in_dir"))+7);
      if (!filename) error_exit("malloc error for filename");
      strcpy(filename,dict_get(m->dict,"in_dir"));
      strcat(filename,"GEOMETRY.OUT");
    }
    else
      filename="GEOMETRY.OUT";

    f=fopen(filename,"r");
    if (f){
      fprintf(stderr,"Additionally reading %s\n",filename);
      elk_read(f,c,m,kp,e);
      fclose(f);
      if (dict_get(m->dict,"in_dir")) free(filename);
    }
    else
      fprintf(stderr,"Failed to open %s, so atomic positions unavailable.\n",
              filename);
  }
  
  elk_readline(line,LINE_SIZE,infile);
  if (sscanf(line,"%lf %lf %lf %lf",v,v+1,v+2,&x)!=4)
    error_exit("failed to read first data line from Elk 3D file");

  if ((v[0]!=0)||(v[1]!=0)||(v[2]!=0))
    error_exit("confused: grid appears not to start at origin");

  *gptr->data=x;
  
  for(k=0;k<ngz;k++){
    for(j=0;j<ngy;j++){
      for(i=0;i<ngx;i++){
	if (i+j+k==0) continue;
	if (!fgets(line,LINE_SIZE,infile)) error_exit("Unexpected EOF");
	if (sscanf(line,"%lf %lf %lf %lf",v,v+1,v+2,
		   gptr->data+i*ngy*ngz+j*ngz+k)!=4)
	  error_exit("failed to parse data line from Elk 3D file");
	if ((k==0)&&(j==0)&&(i==ngx-1))
	  for(ii=0;ii<3;ii++)
	    basis[0][ii]=v[ii];
	if ((k==0)&&(j==ngy-1)&&(i==0))
	  for(ii=0;ii<3;ii++)
	    basis[1][ii]=v[ii];
	if ((k==ngz-1)&&(j==0)&&(i==0))
	  for(ii=0;ii<3;ii++)
	    basis[2][ii]=v[ii];
      }
    }
  }

  for(i=0;i<3;i++)
    basis[0][i]*=ngx*BOHR/(ngx-1);
  for(i=0;i<3;i++)
    basis[1][i]*=ngy*BOHR/(ngy-1);
  for(i=0;i<3;i++)
    basis[2][i]*=ngz*BOHR/(ngz-1);

  if(!c->basis){
    c->basis=malloc(9*sizeof(double));
    if (!c->basis) error_exit("malloc error for basis");
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	c->basis[i][j]=basis[i][j];
    real2rec(c);
  }
  else{
    okay=1;
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	if (!aeq(basis[i][j],c->basis[i][j])) okay=0;
    if (!okay){
      fprintf(stderr,"Warning: bases differ. From GEOMETRY.OUT:\n");
      for(i=0;i<3;i++)
	fprintf(stderr,"%.8f %.8f %.8f\n",c->basis[i][0],c->basis[i][1],
		c->basis[i][2]);
      fprintf(stderr,"From this file:\n");
      for(i=0;i<3;i++)
	fprintf(stderr,"%.8f %.8f %.8f\n",basis[i][0],basis[i][1],
		basis[i][2]);
    }
  }

  scale=1;
  filename=dict_get(m->dict,"in_file");
  ptr=strrchr(filename,'/');
  if (!ptr) ptr=filename;
  else ptr++;
  filename=malloc(strlen(ptr)+1);
  strcpy(filename,ptr);
  ptr=strstr(filename,"3D.OUT");
  if (!ptr) {
    fprintf(stderr,"Confused by name: not rescaling");
    return;
  }
  *ptr=0;
  if (!(strcmp(filename,"RHO"))){
    if (flags&RAW){
      gptr->name="Density_raw";
      add_cmt(m->comment,"Densities are raw");
    }
    else if (!(flags&AU)) {
      scale=1/(BOHR*BOHR*BOHR);
      if (debug>1) fprintf(stderr,"Rescaling density from Bohrs to A\n");
      gptr->name="Density";
      add_cmt(m->comment,"Densities in e A^-3");
    }
    else{
      gptr->name="Density";
      add_cmt(m->comment,"Densities in e Bohr^-3");
    }
  }
  else if ((!(strcmp(filename,"VCL")))||(!(strcmp(filename,"VXC")))){
    if (flags&RAW){
      gptr->name="Potential_raw";
      add_cmt(m->comment,"Potentials are raw");
    }
    else if (!(flags&AU)){
      scale=H_eV;
      if (debug>1)
	fprintf(stderr,"Rescaling potential from Hartrees to eV\n");
      gptr->name="Potential";
      add_cmt(m->comment,"Potentials in eV");
    }
    else{
      gptr->name="Potential";
      add_cmt(m->comment,"Potentials in Ha");
    }
  }
  
  if (scale!=1)
    for(i=0;i<ngx*ngy*ngz;i++)
      gptr->data[i]*=scale;
}


static int elk_readline(char *buffer, int len, FILE* infile){
  int off;
  char *ptr,*success;

  while((success=fgets(buffer,len,infile))){ /* fgets() always
                                                    null terminates,
                                                    gcc likes extra brackets */

/* Kill trailing spaces and newlines / carriage returns */
    ptr=buffer+strlen(buffer)-1;
    while((ptr>=buffer)&&((*ptr==' ')||(*ptr=='\n')||(*ptr=='\r'))) ptr--;
    *(ptr+1)=0;

/* Eat leading spaces */
    ptr=buffer;
    while(*ptr==' ') ptr++;
/* Skip comments and blank lines */
    if ((*ptr=='!')||(*ptr==0)) continue;
    break;
  }

  if (!success) return(0);

  off=ptr-buffer;
  if (off){
    while(*ptr) {*(ptr-off)=*ptr;ptr++;}
    *(ptr-off)=0;
  }

  return (1);
}
