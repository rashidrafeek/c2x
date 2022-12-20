/* Read a tiny subset of .dx files, hoping to read c2x's output,
 * and APBS/VMD .dx files */

/* MJR 12/2020 */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>

#include "c2xsf.h"

#define LINE_SIZE 132

void dx_read(FILE* infile,struct unit_cell *c, struct grid *gptr){
  double origin[3],delta[3][3];
  int i,j,ndata,have_axes;
  char line[LINE_SIZE+1],*cptr,*cptr2;

  gptr=grid_new(gptr);

  ndata=have_axes=0;
  
  while((ndata==0)||(have_axes==0)){

    if (!fgets(line,LINE_SIZE+1,infile))
      error_exit("read error");
    cptr=line;
    while(isspace(*cptr)) cptr++;
    if (*cptr=='#') continue;
    if (*cptr==0) continue;

    if (strncmp(cptr,"object",5)){
      if (debug>1) fprintf(stderr,"Ignoring %s\n",line);
      continue;
    }
    
    cptr=strstr(cptr,"class");
    if (!cptr) error_exit("object with no class");

    cptr+=5;
    while(isspace(*cptr)) cptr++;

    /* Now we need to process different objects */

    if (!(strncmp(cptr,"array",5))){ /* the data */
      cptr2=strstr(cptr,"items");
      if (!cptr2) error_exit("Failed to find item count for array");
      cptr2+=5;
      sscanf(cptr2,"%d",&ndata);

      cptr2=strstr(cptr2,"data follows");
      if (!cptr2) error_exit("Failed to find array data");
      
      gptr->data=malloc(ndata*sizeof(double));
      if (!gptr->data) error_exit("malloc error for data");
      for(i=0;i<ndata;i++)
	if (fscanf(infile,"%lf",gptr->data+i)!=1)
	  error_exit("error reading array data");
    }
    else if (!(strncmp(cptr,"gridconnections",15))){
      continue;
    }
    else if (!(strncmp(cptr,"gridpositions",13))){
      cptr+=13;
      cptr=strstr(cptr,"counts");
      if (!cptr) error_exit("failed to find gridposition counts");
      cptr+=6;
      if (sscanf(cptr,"%d %d %d",gptr->size,gptr->size+1,gptr->size+2)!=3)
	error_exit("failed to parse counts");
      fgets(line,LINE_SIZE+1,infile);
      cptr=line;
      while(isspace(*cptr)) cptr++;
      if (strncmp(cptr,"origin",6))
	error_exit("Failed to find origin");
      cptr+=6;
      if (sscanf(cptr,"%lf %lf %lf",origin,origin+1,origin+2)!=3)
	error_exit("Failed to parse origin");
      for(i=0;i<3;i++){
	fgets(line,LINE_SIZE+1,infile);
	cptr=line;
	while(isspace(*cptr)) cptr++;
	if (strncmp(cptr,"delta",5))
	  error_exit("Failed to find deltas");
	cptr+=5;
	if (sscanf(cptr,"%lf %lf %lf",
		   &(delta[i][0]),&(delta[i][1]),&(delta[i][2]))!=3)
	  error_exit("Failed to parse deltas");
      }
      have_axes=1;
    }
  }

  /* This follows from the bizarre fixed format of VMD/APBS */

  if (fgets(line,LINE_SIZE+1,infile)){
    cptr=line;
    while(isspace(*cptr)) cptr++;
    while (*cptr==0) {
      if (!fgets(line,LINE_SIZE+1,infile)) break;
      cptr=line;
      while(isspace(*cptr)) cptr++;
    }
    if (!strncmp(cptr,"object ",7)){
      cptr+=7;
      while(isspace(*cptr)) cptr++;
      if ((strstr(cptr,"class"))&&(strstr(cptr,"field"))&&(*cptr=='"')){
	cptr++;
	cptr2=index(cptr,'"');
	if (cptr2){
	  *cptr2=0;
	  gptr->name=malloc((cptr2-cptr)+1);
	  if (!gptr->name) error_exit("malloc error for grid name");
	  strcpy(gptr->name,cptr);
	}
      }
    }
  }
  
  if (ndata!=gptr->size[0]*gptr->size[1]*gptr->size[2]){
    fprintf(stderr,"Error: expected grid of %dx%dx%d=%d points,"
	    " found %d points\n",
	    gptr->size[0],gptr->size[1],gptr->size[2],
	    gptr->size[0]*gptr->size[1]*gptr->size[2],ndata);
    exit(1);
  }

  c->basis=malloc(9*sizeof(double));
  if(!c->basis) error_exit("malloc error for basis");
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[i][j]=delta[i][j]*gptr->size[i];

  real2rec(c);

  if ((origin[0]!=0)||(origin[1]!=0)||(origin[2]!=0)){
    gptr->origin_abs=malloc(3*sizeof(double));
    if (!gptr->origin_abs) error_exit("mallock error");
    for(i=0;i<3;i++)
      gptr->origin_abs[i]=origin[i];
  }
  
}
