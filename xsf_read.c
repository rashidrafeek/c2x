/* Read data set from an XSF file */

/* Copyright (c) 2017 MJ Rutter 
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
#include<math.h> /* fabs from aeq */

#include "c2xsf.h"

int super(struct unit_cell *c, struct contents *m,
           double new_basis[3][3], struct kpts *k, struct symmetry *s,
           struct grid *gptr, int rhs);

#define LINE_SIZE 100

static int xsfreadline(char *buffer, int len, FILE* infile);
static int xsfreadatom(struct contents *m, int i, FILE* infile);
static int line_count;


void xsf_read(FILE* infile, struct unit_cell *c, struct contents *m,
               struct grid *gptr){
  char buffer[LINE_SIZE+1];
  int i,j,k,n,tmp,nn[3],result,have_convvec=0,v2g[3],g2v[3];
  int is,js,ks;
  double x,y,z,b[3][3],convvec[3][3];
  
  if (debug>2) fprintf(stderr,"xsf_read called\n");

  xsfreadline(buffer,LINE_SIZE,infile);

  if(!strncasecmp(buffer,"atoms",5)){
    n=0;
    m->atoms=malloc(sizeof(struct atom));
    while(xsfreadatom(m,n,infile)){
      n++;
      m->atoms=realloc(m->atoms,n*sizeof(struct atom));
    }
    m->n=n;
    if (debug) fprintf(stderr,"%d atoms read from xsf molecule file\n",n);
    fprintf(stderr,"Warning, no unit cell in xsf file\n"
	    "Creating dummy 10A box\n");
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	c->basis[i][j]=0;
    for(i=0;i<3;i++) c->basis[i][i]=10;
    return;
  }
  else if(!strncasecmp(buffer,"crystal",7)){
    while(xsfreadline(buffer,LINE_SIZE,infile)){
      if(!strncasecmp(buffer,"primvec",7)){
	if (!(c->basis=malloc(72)))
	  error_exit("Malloc error in xsf_read for c->basis");
	for(i=0;i<3;i++){
	  xsfreadline(buffer,LINE_SIZE,infile);
	  if (sscanf(buffer,"%lf %lf %lf",&(c->basis[i][0]),
		     &(c->basis[i][1]),&(c->basis[i][2]))!=3)
	    error_exit("error parsing primvec");
	}
      }
      else if(!strncasecmp(buffer,"convvec",7)){
        for(i=0;i<3;i++){
	  xsfreadline(buffer,LINE_SIZE,infile);
	  if (sscanf(buffer,"%lf %lf %lf",&(convvec[i][0]),
		     &(convvec[i][1]),&(convvec[i][2]))!=3)
	    error_exit("error parsing convvec");
	}
	have_convvec=1;
      }
      else if(!strncasecmp(buffer,"primcoord",9)){
	xsfreadline(buffer,LINE_SIZE,infile);
	if (sscanf(buffer,"%d %d",&n,&tmp)!=2)
	  error_exit("error parsing first line after primccord");
	if (tmp!=1)
	  error_exit("error parsing first line after primcoord");
	if(!(m->atoms=malloc(n*sizeof(struct atom))))
	  error_exit("Malloc error in xsf_read for atoms");
	
	for(i=0;i<n;i++){
	  xsfreadatom(m,i,infile);
	}
	m->n=n;
	real2rec(c);
        c->vol=fabs(c->vol);
	addfrac(m->atoms,m->n,c->recip);
	if (debug>1) fprintf(stderr,"%d atoms read in primcoord\n",n);
	if (have_convvec){
	  result=1;
	  for(i=0;i<3;i++)
	    for(j=0;j<3;j++)
	      if (!aeq(convvec[i][j],c->basis[i][j])) result=0;
	  if (result==0){
	    if (debug)
	      fprintf(stderr,"convvec and primvec differ, calling super\n");
	    if (super(c,m,convvec,NULL,NULL,NULL,0))
	      error_exit("Cell conversion from primvec to convvec failed");
	  }
	}
      }
      else if(!strncasecmp(buffer,"convcoord",9)){
	xsfreadline(buffer,LINE_SIZE,infile);
	if (sscanf(buffer,"%d %*d",&tmp)!=2)
	  error_exit("error parsing first line after convcoord");
	if (m->atoms){
	  if (debug) fprintf(stderr,"Discarding convcoord section\n");
	  for(i=0;i<tmp;i++)
	    xsfreadline(buffer,LINE_SIZE,infile);
	}
	else{
	  n=tmp;
	  if(!(m->atoms=malloc(n*sizeof(struct atom))))
	    error_exit("Malloc error in xsf_read for atoms");
	  
	  for(i=0;i<n;i++){
	    xsfreadatom(m,i,infile);
	  }
	  m->n=n;
	  real2rec(c);
          c->vol=fabs(c->vol);
	  addfrac(m->atoms,m->n,c->recip);
	}
      }
      else if(!strncasecmp(buffer,"begin_block_datagrid_3d",23)){
      /* read title */
	xsfreadline(buffer,LINE_SIZE,infile);
	xsfreadline(buffer,LINE_SIZE,infile);
	do{
          /* The spec says begin_datagrid_3d, but XCrysDen supports
           * simply datagrid_3d, and some programs use this... */
	  if((!strncasecmp(buffer,"begin_datagrid_3d",17))||
             (!strncasecmp(buffer,"datagrid_3d",11))){
	    if (gptr->next) gptr=gptr->next;
	    gptr->next=malloc(sizeof(struct grid));
	    if (!gptr->next) error_exit("Malloc error for struct grid");
	    gptr->next->next=NULL;
	    gptr->next->data=NULL;
	    gptr->next->name=NULL;
	    if ((!strncasecmp(buffer,"begin_datagrid_3d",17))&&
                (strlen(buffer)>18)){
	      gptr->name=malloc(strlen(buffer)-17);
	      strcpy(gptr->name,buffer+18);
	    }
	    if ((!strncasecmp(buffer,"datagrid_3d",11))&&
                (strlen(buffer)>12)){
	      gptr->name=malloc(strlen(buffer)-11);
	      strcpy(gptr->name,buffer+12);
	    }
	    if (debug>1) fprintf(stderr,"Start of data for %s\n",
				 gptr->name?gptr->name:"(unknown)");
	    
	    xsfreadline(buffer,LINE_SIZE,infile);
	    if(sscanf(buffer,"%d %d %d",nn,nn+1,nn+2)!=3)
	      error_exit("Error parsing size in datagrid_3d");
	    /* XSF grids are periodic */
	    for(i=0;i<3;i++) nn[i]--;
	    
	    xsfreadline(buffer,LINE_SIZE,infile);
	    if(sscanf(buffer,"%lf %lf %lf",&x,&y,&z)!=3)
	      error_exit("Error parsing displacement in datagrid_3d");
	    if ((x!=0)||(y!=0)||(z!=0))
	      error_exit("Non-zero displacement in datagrid_3d not supported");
	    
	    for(i=0;i<3;i++){
	      xsfreadline(buffer,LINE_SIZE,infile);
	      if (sscanf(buffer,"%lf %lf %lf",&(b[i][0]),
			 &(b[i][1]),&(b[i][2]))!=3)
		error_exit("error parsing datagrid_3d axes");
	    }
	    
	    if(!(gptr->data=malloc(nn[0]*nn[1]*nn[2]*sizeof(double))))
	      error_exit("Malloc error for 3D grid data");

	    for(i=0;i<3;i++) v2g[i]=g2v[i]=-1;
	    for(j=0;j<3;j++){
	      for(i=0;i<3;i++){
		if((aeq(b[i][0],c->basis[j][0]))&&
		   (aeq(b[i][1],c->basis[j][1]))&&
		   (aeq(b[i][2],c->basis[j][2]))){
		  g2v[i]=j;
		  v2g[j]=i;
		}
	      }
	    }

	    if((g2v[0]==-1)||(g2v[1]==-1)||(g2v[2]==-1))
	      error_exit("Datagrid_3d axes are not convvecs");

	    if (debug&&((g2v[0]!=0)||(g2v[1]!=1)||(g2v[2]!=2)))
	      fprintf(stderr,"Permuting axes in data grid\n");

	    for(i=0;i<3;i++)
	      gptr->size[i]=nn[v2g[i]];

	    if (g2v[2]==2) ks=1;
	    else if (g2v[2]==1) ks=gptr->size[2];
	    else ks=gptr->size[2]*gptr->size[1];
	    if (g2v[1]==2) js=1;
	    else if (g2v[1]==1) js=gptr->size[2];
	    else js=gptr->size[2]*gptr->size[1];
	    if (g2v[0]==2) is=1;
	    else if (g2v[0]==1) is=gptr->size[2];
	    else is=gptr->size[2]*gptr->size[1];

	    for(k=0;k<nn[2];k++){
	      for(j=0;j<nn[1];j++){
		for(i=0;i<nn[0];i++)
		  fscanf(infile,"%lf",gptr->data+k*ks+j*js+i*is);
		fscanf(infile,"%*f");
	      }
	      for(i=0;i<=nn[0];i++)
		fscanf(infile,"%*f");
	    }
	    for(i=0;i<(nn[0]+1)*(nn[1]+1);i++)
	      fscanf(infile,"%*f");
    
	    xsfreadline(buffer,LINE_SIZE,infile);
	    if (strncasecmp(buffer,"end_datagrid_3d",15))
	      error_exit("end_datagrid_3d not found");
	  } /* if ! strncmp begin_datagrid_3d */
	  if (!xsfreadline(buffer,LINE_SIZE,infile)) break;
	}while (strncasecmp(buffer,"end_block_datagrid_3d",21));
	if (strncasecmp(buffer,"end_block_datagrid_3d",21))
	  fprintf(stderr,"Warning, end_block_datagrid_3d not found\n"
		  "Ended on %s\n",buffer);
      }   /* if ! strncmp block_datagrid */
    } /* while */
  } /* if ! strncmp crystal */
  else
    error_exit("XSF file not understood");
  
}

static int xsfreadline(char *buffer, int len, FILE* infile){
  int off;
  char *ptr,*success;

  while((success=fgets(buffer,len,infile))){ /* fgets() always
                                                    null terminates,
                                                    gcc likes extra brackets */
    line_count++;

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

  if (!success) return(0);

  off=ptr-buffer;
  if (off){
    while(*ptr) {*(ptr-off)=*ptr;ptr++;}
    *(ptr-off)=0;
  }

  return (1);
}

static int xsfreadatom(struct contents *m, int i, FILE* infile){
  int j;
  char buffer[LINE_SIZE+1],sym[4];
  
  if (!xsfreadline(buffer,LINE_SIZE,infile)) return 0;
  if (sscanf(buffer,"%d %lf %lf %lf %lf %lf %lf",
	     &(m->atoms[i].atno),
	     &(m->atoms[i].abs[0]),
	     &(m->atoms[i].abs[1]),
	     &(m->atoms[i].abs[2]),
	     &(m->atoms[i].force[0]),
	     &(m->atoms[i].force[1]),
	     &(m->atoms[i].force[2]))==7) m->forces=1;
  else if (sscanf(buffer,"%3s %lf %lf %lf %lf %lf %lf",
		  sym,
		  &(m->atoms[i].abs[0]),
		  &(m->atoms[i].abs[1]),
		  &(m->atoms[i].abs[2]),
		  &(m->atoms[i].force[0]),
		  &(m->atoms[i].force[1]),
		  &(m->atoms[i].force[2]))==7){
    m->forces=1;
    m->atoms[i].atno=atsym2no(sym);
  }
  else if (sscanf(buffer,"%d %lf %lf %lf",
		  &(m->atoms[i].atno),
		  &(m->atoms[i].abs[0]),
		  &(m->atoms[i].abs[1]),
		  &(m->atoms[i].abs[2]))==4)
    for(j=0;j<3;j++) m->atoms[i].force[j]=0.0;
  else if (sscanf(buffer,"%3s %lf %lf %lf",
		  sym,
		  &(m->atoms[i].abs[0]),
		  &(m->atoms[i].abs[1]),
		  &(m->atoms[i].abs[2]))==4){
    for(j=0;j<3;j++) m->atoms[i].force[j]=0.0;
    m->atoms[i].atno=atsym2no(sym);
  }
  else	  
    error_exit("error parsing atoms line in xsf_read");
  m->atoms[i].spin=0;
  m->atoms[i].chg=0;
  m->atoms[i].label=NULL;
  return 1;
}
