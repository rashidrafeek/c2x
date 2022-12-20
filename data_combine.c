#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "c2xsf.h"

static int cells_equal(struct unit_cell *c, struct unit_cell *c2){
  double abc[6],abc2[6];
  int i,ok;

  if ((!c->basis)||(!c2->basis)) return 1;

  basis2abc(c->basis,abc);
  basis2abc(c2->basis,abc2);

  ok=1;
  for(i=0;i<6;i++)
    if (!aeq(abc[i],abc2[i])) ok=0;

  if ((debug>1)&&(ok==0)){
    fprintf(stderr,"Cells %f %f %f %f %f %f\n"
	    "  and %f %f %f %f %f %f compare different\n",
	    abc[0],abc[1],abc[2],abc[3],abc[4],abc[5],
	    abc2[0],abc2[1],abc2[2],abc2[3],abc2[4],abc2[5]);
  }
  
  return ok;
  
}

void data_combine(int op, struct unit_cell *c,
		  struct contents *m, struct kpts *k, struct symmetry *s,
		  struct grid *g, struct es *elect, struct time_series *ts,
		  int *i_grid,struct unit_cell *c2,
		  struct contents *m2, struct kpts *k2, struct symmetry *s2,
		  struct grid *g2, struct es *elect2, struct time_series *ts2){

  struct grid *gptr,*gptr2,ng;
  int i,ngpts,ok;
  char *name2;

  name2=(char*)dict_get(m2->dict,"in_file");
  if (!name2) name2="(unknown)";
  
  if (op==C2X_MERGE){
    if ((!c->basis)&&(c2->basis)){
      *c=*c2;
      if (debug) fprintf(stderr,"Taking cell from %s\n",name2);
    }
    if ((!m->n)&&(m2->n)){
      *m=*m2;
      if (debug) fprintf(stderr,"Taking atoms from %s\n",name2);
    }
    if ((!k->n)&&(k2->n)){
      *k=*k2;
      if (debug) fprintf(stderr,"Taking kpoints from %s\n",name2);
    }
    if ((!s->n)&&(s->n)){
      *s=*s2;
      if (debug) fprintf(stderr,"Taking symmetry from %s\n",name2);
    }
    if (g2->data){
      if (!cells_equal(c,c2)){
	fprintf(stderr,"Unable to merge data from different cells\n");
	return;
      }
      gptr=g;
      while(gptr->next) gptr=gptr->next;
      while(g2->data){
	gptr=grid_new(gptr);
	gptr->data=g2->data;
	for(i=0;i<3;i++) gptr->size[i]=g2->size[i];
	gptr->name=g2->name;
	gptr->origin_abs=g2->origin_abs;
	if (debug)
	  fprintf(stderr,"Merging %s from %s\n",
		  (g2->name)?g2->name:"unknown data",name2);
	g2=g2->next;
      }
    }
    return;
  }

  gptr=g;
  gptr2=g2;
  while((gptr->data)&&(gptr2->data)){
    if (i_grid){
      if ((i_grid[0]==0)&&(i_grid[1]==0)&&(i_grid[2]==0))
	for(i=0;i<3;i++)
	  i_grid[i]=max(gptr->size[i],gptr2->size[i]);
      ok=1;
      for(i=0;i<3;i++)
	if (i_grid[i]!=gptr->size[i]) ok=0;
      if (!ok){
	fprintf(stderr,"Interpolating %s from %s\n",gptr->name,
		(char*)dict_get(m->dict,"in_file"));
	for(i=0;i<3;i++) ng.size[i]=i_grid[i];

	interpolate3d(gptr,&ng);

	free(gptr->data);
	gptr->data=ng.data;
	for(i=0;i<3;i++) gptr->size[i]=ng.size[i];
      }
      
      ok=1;
      for(i=0;i<3;i++)
	if (i_grid[i]!=gptr2->size[i]) ok=0;
      if (!ok){
	fprintf(stderr,"Interpolating %s from %s\n",gptr2->name,
		(char*)dict_get(m2->dict,"in_file"));
	for(i=0;i<3;i++) ng.size[i]=i_grid[i];

	interpolate3d(gptr2,&ng);

	free(gptr2->data);
	gptr2->data=ng.data;
	for(i=0;i<3;i++) gptr2->size[i]=ng.size[i];
      }
    }
    if ((gptr->size[0]!=gptr2->size[0])||
	(gptr->size[1]!=gptr2->size[1])||
	(gptr->size[2]!=gptr2->size[2])){
      fprintf(stderr,"Unable to work with different grid sizes,"
	      " have %dx%dx%d and %dx%dx%d\n",
	      gptr->size[0],gptr->size[1],gptr->size[2],
	      gptr2->size[0],gptr2->size[1],gptr2->size[2]);
      exit(1);
    }
    if (!cells_equal(c,c2)){
      fprintf(stderr,"Unable to combine data from different cells\n");
      exit(1);
    }
    ngpts=gptr->size[0]*gptr->size[1]*gptr->size[2];
    if (op==C2X_ADD){
      if (debug)
	fprintf(stderr,"Adding %s from %s to %s\n",
		(gptr2->name)?gptr2->name:"unknown data",name2,
		(gptr->name)?gptr->name:"unknown data");
      for(i=0;i<ngpts;i++)
	gptr->data[i]+=gptr2->data[i];
      free(gptr2->data);
      gptr2->data=NULL;
      gptr2=gptr2->next;
    }
    else if (op==C2X_DIFF){
      if (debug)
	fprintf(stderr,"Subtracting %s from %s from %s\n",
		(gptr2->name)?gptr2->name:"unknown data",name2,
		(gptr->name)?gptr->name:"unknown data");
      for(i=0;i<ngpts;i++)
	gptr->data[i]-=gptr2->data[i];
      free(gptr2->data);
      gptr2->data=NULL;
      gptr2=gptr2->next;
    }
    else if (op==C2X_MASK){
      if (debug)
	fprintf(stderr,"Masking with %s from %s\n",
		(gptr2->name)?gptr2->name:"unknown data",name2);

      for(i=0;i<ngpts;i++)
	gptr->data[i]*=gptr2->data[i];
    }
    gptr=gptr->next;
  }

  if (op==C2X_MASK){
    free(gptr2->data);
    gptr2->data=NULL;
  }

}
