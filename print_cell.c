#include<stdio.h>
#include<string.h>
#include<math.h>

#include "c2xsf.h"

void print_cell(struct unit_cell *c, struct contents *m){
  int i;
  char *fmt;
  double abc[6];

  if (flags&HIPREC)
    fmt="% 19.15f % 19.15f % 19.15f\n";
  else
    fmt="% 11.7f % 11.7f % 11.7f\n";
  fprintf(stderr,"Lattice:\n");
  for(i=0;i<3;i++)
      fprintf(stderr,fmt,
                    c->basis[i][0],c->basis[i][1],c->basis[i][2]);
  fprintf(stderr,"\n");

  if (flags&HIPREC)
    fmt="%c=% 19.15f ";
  else
    fmt="%c=% 11.7f ";
  cart2abc(c,NULL,abc,NULL,0);
  for(i=0;i<3;i++)
  fprintf(stderr,fmt,'a'+i,abc[i]);
  if (flags&HIPREC)
    fmt="\nalpha=% 19.15f beta=% 19.15f gamma=% 19.15f\n";
  else
    fmt="\nalpha=% 11.7f beta=% 11.7f gamma=% 11.7f\n";
  fprintf(stderr,fmt,abc[3],abc[4],abc[5]);
    
  if (m){
    fprintf(stderr,"\nAtomic positions:\n");
    if (flags&HIPREC)
      fmt="%3s % 19.15f % 19.15f % 19.15f\n";
    else
      fmt="%3s % 11.7f % 11.7f % 11.7f\n";
    for(i=0;i<m->n;i++)
      fprintf(stderr,fmt,atno2sym(m->atoms[i].atno),m->atoms[i].frac[0],
	    m->atoms[i].frac[1],m->atoms[i].frac[2]);
  }
  
}
