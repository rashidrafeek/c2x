/* Write an abinit .in file */


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include "c2xsf.h"

void abinit_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                  struct kpts *k, struct symmetry *s, struct es *e){
  double abc[6],acell[3],akpt_disp[3],tot_spin;
  int i,j,hit,nspec,nspin;
  int *natomsp,*spatno;
  char *fmt;
  
  if (m->title) fprintf(outfile,"# %s\n",m->title);

  /* The cell */
  
  cart2abc(c,NULL,abc,NULL,0);
  for(i=0;i<3;i++) acell[i]=abc[i];

  /* Abinit likes to make rprim into rationals by modifying acell.
   * Deal with one special case:
   *  convert (sqrt(0.5),sqrt(0.5),0) and its permutations to (0.5,0.5,0)
   */
  for(i=0;i<3;i++){
    hit=0;
    for(j=0;j<3;j++)
      if ((fabs(c->basis[i][j])<tol)||
	  (aeq(fabs(c->basis[i][j]/acell[i]),sqrt(0.5)))) hit++;
    if (hit==3) acell[i]*=sqrt(2);
  }
  /* and one more */
  for(i=0;i<3;i++){
    if ((aeq(fabs(c->basis[i][0]/acell[i]),sqrt(1./3.)))&&
	(aeq(fabs(c->basis[i][1]/acell[i]),sqrt(1./3.)))&&
	(aeq(fabs(c->basis[i][2]/acell[i]),sqrt(1./3.))))
      acell[i]*=2*sqrt(1./3.);
  }

  if (flags&HIPREC)
    fmt="acell %.15g %.15g %.15g\n";
  else
    fmt="acell %.12g %.12g %.12g\n";
  fprintf(outfile,fmt,acell[0]/BOHR,
          acell[1]/BOHR,acell[2]/BOHR);

  fprintf(outfile,"rprim\n");
  if (flags&HIPREC)
    fmt="     %18.15f %18.15f %18.15f\n";
  else
    fmt="     %14.11f %14.11f %14.11f\n";
  for(i=0;i<3;i++)
    fprintf(outfile,fmt,
            c->basis[i][0]/acell[i],c->basis[i][1]/acell[i],
	    c->basis[i][2]/acell[i]);

  /* Atom types */

  /* Now we need to know the number of species.
   * It must be fewer than the number of atoms...
   */

  nspec=0;
  natomsp=malloc(m->n*sizeof(int));
  if (!natomsp) error_exit("Malloc error in abinit_write");
  spatno=malloc(m->n*sizeof(int));
  if (!spatno) error_exit("Malloc error in abinit_write");

  for(i=0;i<m->n;i++){
    for(j=0;j<nspec;j++) if (m->atoms[i].atno==spatno[j]) break;
    if (j==nspec){  /* new species */
      spatno[j]=m->atoms[i].atno;
      natomsp[j]=1;
      nspec++;
    }else{          /* existing species */
      natomsp[j]++;
    }
  }

  fprintf(outfile,"\nntypat %d\n",nspec);
  fprintf(outfile,"znucl ");
  for(i=0;i<nspec;i++) fprintf(outfile," %d",spatno[i]);
  fprintf(outfile,"\n\n");

  /* Atoms */

  fprintf(outfile,"natom %d\n",m->n);
  fprintf(outfile,"typat");
  for(i=0;i<m->n;i++){
    for(j=0;j<nspec;j++)
      if(m->atoms[i].atno==spatno[j]){
        fprintf(outfile," %2d",j+1);
        break;
      }
    if ((i%16)==15) fprintf(outfile,"\n     ");
  }
  fprintf(outfile,"\n");
  fprintf(outfile,"xred\n");
  if (flags&HIPREC)
    fmt="  %.15g";
  else
    fmt="  %.11g";
  for(i=0;i<m->n;i++){
    for(j=0;j<3;j++) fprintf(outfile,fmt,m->atoms[i].frac[j]);
    fprintf(outfile,"\n");
  }
  fprintf(outfile,"\n");

  /* Spin */

  nspin=0;
  tot_spin=0;
  for(i=0;i<m->n;i++){
    if (m->atoms[i].spin) nspin=1;
    tot_spin+=m->atoms[i].spin;
  }
    
  if (nspin){
    if (tot_spin==0){
      fprintf(outfile,"\n# antiferromagnetic system assumed\n");
      fprintf(outfile,"nsppol 1\nnspden 2\n");
    }
    else{
      fprintf(outfile,"\nnsppol 2\n");
      fprintf(outfile,"occopt 7  # specify some form of metallic smearing\n"
	      "                  # for ferromagnetic spin\n");
    }
    fprintf(outfile,"spinat");
    for(i=0;i<m->n;i++){
      fprintf(outfile,"  0 0 %f",m->atoms[i].spin);
      if ((i%4)==3) fprintf(outfile,"\n      ");
    }
    fprintf(outfile,"\n");
  }
  
  /* k-points */

  if ((k->mp)&&(k->mp->grid[0]>0)){
    fprintf(outfile,"kptopt 1\n");
    fprintf(outfile,"ngkpt %d %d %d\n",k->mp->grid[0],
            k->mp->grid[1],k->mp->grid[2]);
    fprintf(outfile,"nshiftk 1\n");
    /* Convert from Castep's convention to Abinit's */
    for(i=0;i<3;i++){
      akpt_disp[i]=k->mp->disp[i]*k->mp->grid[i];
      if ((k->mp->grid[i]&1)==0) akpt_disp[i]+=0.5; /* is grid even? */
      akpt_disp[i]=fmod(akpt_disp[i],1.0);
    }
    if (flags&HIPREC)
      fmt="shiftk %.15g %.15g %.15g\n";
    else
      fmt="shiftk %.11g %.11g %.11g\n";
    fprintf(outfile,fmt,akpt_disp[0],akpt_disp[1],akpt_disp[2]);
  }
  else if(k->n) {
    fprintf(outfile,"kptopt 0\n");
    fprintf(outfile,"nkpt %d\n",k->n);
    fprintf(outfile,"kpt\n");
    if (flags&HIPREC)
      fmt="  %.15f %.15f %.15f\n";
    else
      fmt="  %.11f %.11f %.11f\n";
    for(i=0;i<k->n;i++)
      fprintf(outfile,fmt,k->kpts[i].frac[0],
              k->kpts[i].frac[1],k->kpts[i].frac[2]);
    fprintf(outfile,"wtk\n");
    if (flags&HIPREC)
      fmt="%.15f\n";
    else
      fmt="%.11f\n";
    for(i=0;i<k->n;i++)
      fprintf(outfile,fmt,k->kpts[i].wt);
  }


  if (e->cut_off)
    if (flags&AU)
      fprintf(outfile,"\necut %f\n",e->cut_off/H_eV);
    else
      fprintf(outfile,"\necut %f eV\n",e->cut_off);
  else
    fprintf(outfile,"\n#ecut 100 eV\n");

  if (e->etol)
    if (flags&AU)
      fprintf(outfile,"\ntoldfe %g\n",e->etol/H_eV);
    else
      fprintf(outfile,"\ntoldfe %g eV\n",e->etol*m->n);
  else
    if (flags&AU)
      fprintf(outfile,"\ntoldfe %g\n",3.675e-7*m->n);
    else
      fprintf(outfile,"\ntoldfe %g eV\n",1e-5*m->n);

  fprintf(outfile,"\nchksymbreak 0\nchkprim 0\n");

  if ((e->charge)&&(*e->charge!=0))
    fprintf(outfile,"charge %f\n",*e->charge);
  
}
