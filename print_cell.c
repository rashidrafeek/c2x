#include<stdio.h>
#include<string.h>
#include<math.h>

#include "c2xsf.h"

void print_basis(double basis[3][3]){
  int i;
  char *fmt;
  
  if (flags&HIPREC)
    fmt="% 19.15f % 19.15f % 19.15f   modulus %19.15f\n";
  else
    fmt="% 11.7f % 11.7f % 11.7f   modulus %11.7f\n";
  for(i=0;i<3;i++)
      fprintf(stderr,fmt,basis[i][0],basis[i][1],basis[i][2],
              sqrt(vmod2(basis[i])));
}

void print_cell(struct unit_cell *c, struct contents *m){
  int i;
  char *fmt;
  double abc[6];

  if (flags&HIPREC)
    fmt="% 19.15f % 19.15f % 19.15f\n";
  else
    fmt="% 11.7f % 11.7f % 11.7f\n";
  fprintf(stderr,"Lattice:\n");
  print_basis(c->basis);
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

void print_occ(struct es *elect, struct kpts *kp){
  double total,wtotal,scale;
  int i,k,ns,b;

  scale=1;
  if (flags&AU) scale=1/H_eV;
  
  if (!elect->occ){
    fprintf(stderr,"No occupancies found to report\n");
    return;
  }
  fprintf(stderr,"                   kpoint              band spin "
          " occupancy      evalue (%s)\n",(scale==1)?"eV":"Ha");
  i=0;
  total=wtotal=0;
  for(k=0;k<kp->n;k++)
    for(ns=0;ns<elect->nbspins;ns++)
      for(b=0;b<elect->nbands;b++){
        fprintf(stderr,"%3d: ( % 8f % 8f % 8f )  %3d  %d  %10f   %14f\n",k+1,
                kp->kpts[k].frac[0],kp->kpts[k].frac[1],kp->kpts[k].frac[2],
                b+1,ns,elect->occ[i],
                elect->eval?elect->eval[i]*scale:0);
        total+=elect->occ[i];
        wtotal+=kp->kpts[k].wt*elect->occ[i];
        i++;
      }
  fprintf(stderr,"                                       Total:  %11f\n",
		  total);
  
  fprintf(stderr,"                              Weighted total:  %11f\n",
		  wtotal);
}

void print_elect(struct es *elect){
  double scale;
  scale=1;
  if (flags&AU) scale=1/H_eV;
  
  if ((elect->nspins!=1)||(elect->nbspins!=1)){
    fprintf(stderr,"Spin components (density): %d\n",elect->nspins);
    fprintf(stderr,"Spin components (bands):   %d\n",elect->nbspins);
  }
  if (elect->nspinors!=1)
    fprintf(stderr,"Spinors: %d\n",elect->nspinors);
  if (elect->nbands) fprintf(stderr,"Bands: %d\n",elect->nbands);
  if (elect->charge) fprintf(stderr,"Total charge: %f e\n",*elect->charge);
  if (elect->energy) fprintf(stderr,"Total energy: %.6f %s\n",
                             *elect->energy*scale,(scale==1)?"eV":"Ha");
  if (elect->e_fermi) fprintf(stderr,"Fermi energy: %.6f %s\n",
                              *elect->e_fermi*scale,(scale==1)?"eV":"Ha");
}
