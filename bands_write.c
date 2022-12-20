/* Write a .bands file (CASTEP-style) */
/* (c) MJR 2019-2020 */

#include<stdio.h>

#include "c2xsf.h"

void bands_write(FILE* outfile, struct unit_cell *c,
                 struct kpts *k, struct es *e){
  int i,ik,nb,ns;

  if (!e->eval) {
    fprintf(stderr,"No eigenvalues to write!\n");
    return;
  }
  
  /* Print header */
  fprintf(outfile,"Number of k-points %d\n",k->n);
  fprintf(outfile,"Number of spin components %d\n",e->nspins);
  if (e->nspins==1){
    fprintf(outfile,"Number of electrons %.4g\n",e->nel);
    fprintf(outfile,"Number of eigenvalues %d\n",e->nbands);
    fprintf(outfile,"Fermi energy (in atomic units) %12f\n",
            (e->e_fermi)?(*e->e_fermi)/H_eV:0.0);
  }
  else{
    fprintf(outfile,"Number of electrons %.4g %.4g\n",
            0.5*(e->nel+e->nup_minus_down),0.5*(e->nel-e->nup_minus_down));
    fprintf(outfile,"Number of eigenvalues %d %d\n",e->nbands,e->nbands);
    fprintf(outfile,"Fermi energies (in atomic units) %12f %12f\n",
            (e->e_fermi)?(*e->e_fermi)/H_eV:0.0,
            (e->e_fermi)?(*e->e_fermi)/H_eV:0.0);
  }
  fprintf(outfile,"Unit cell vectors\n"); /* In Bohr... */
  for(i=0;i<3;i++)
    fprintf(outfile,"%12.6f %12.6f %12.6f\n",c->basis[i][0]/BOHR,
            c->basis[i][1]/BOHR,c->basis[i][2]/BOHR);
  
  /* Loop over bands etc */
  i=0;
  for(ik=0;ik<k->n;ik++){
    fprintf(outfile,"K-point %4d %12.8f %12.8f %12.8f %12.8f\n",
            ik+1,k->kpts[ik].frac[0],k->kpts[ik].frac[1],k->kpts[ik].frac[2],
            k->kpts[ik].wt);

    for(ns=0;ns<e->nspins;ns++){
      fprintf(outfile,"Spin component %d\n",ns+1);

      for(nb=0;nb<e->nbands;nb++)
        fprintf(outfile,"%14.8f\n",e->eval[i++]/H_eV);

    }
  }
}
