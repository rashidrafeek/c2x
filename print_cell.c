#include<stdio.h>
#include<stdlib.h>
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
  double total,wtotal,scale,wt;
  int i,k,ns,b;

  scale=1;
  if (flags&AU) scale=1/H_eV;
  
  if ((!elect->occ)&&(!elect->eval)){
    fprintf(stderr,"No occupancies or evals found to report\n");
    return;
  }
  fprintf(stderr,"                   kpoint              band spin "
          " occupancy%s   evalue (%s)\n",(flags&K_WEIGHT)?"*wt":"   ",
	  (scale==1)?"eV":"Ha");
  i=0;
  total=wtotal=0;
  for(k=0;k<kp->n;k++){
    wt=1;
    if (flags&K_WEIGHT) wt=kp->kpts[k].wt;
    if ((flags&OCC_WEIGHT)&&(elect->nspins==1)&&(elect->nspinors==1)) wt*=2;
    for(ns=0;ns<elect->nbspins;ns++){
      for(b=0;b<elect->nbands;b++){
        fprintf(stderr,"%3d: ( % 8f % 8f % 8f )  %3d  %d  %11.7f  %14f\n",k+1,
                kp->kpts[k].frac[0],kp->kpts[k].frac[1],kp->kpts[k].frac[2],
                b+1,ns,elect->occ?elect->occ[i]*wt:0,
                elect->eval?elect->eval[i]*scale:0);
        if (elect->occ) {
	  total+=elect->occ[i];
	  wtotal+=kp->kpts[k].wt*elect->occ[i];
	}
        i++;
      }
    }
  }
  fprintf(stderr,"                                       Total:  %11f\n",
		  total);
  
  if ((flags&OCC_WEIGHT)&&(elect->nspins==1)&&(elect->nspinors==1)) wtotal*=2;
  fprintf(stderr,"                              Weighted total:  %11f\n",
		  wtotal);
}

void print_bandwidths(struct es *elect, struct kpts *kp){
  int i,off,b,ns,cross;
  double emin,emax,scale,occ;
  
  scale=1;
  if (flags&AU) scale=1/H_eV;
  
  if (!elect->eval){
    fprintf(stderr,"No evals found to report\n");
    return;
  }

  if (debug==0) fprintf(stderr,"Bands crossing Fermi level:\n");
  
  fprintf(stderr,
	  "   band spin    occupancy       min eval     max eval (%s)\n",
          (scale==1)?"eV":"Ha");
  for(b=0;b<elect->nbands;b++){
    for(ns=0;ns<elect->nbspins;ns++){
      emin=1e100;
      emax=-1e100;
      occ=0;
      for(i=0;i<kp->n;i++){
        off=i*elect->nbspins*elect->nbands+ns*elect->nbands+b;
        emin=min(emin,elect->eval[off]);
        emax=max(emax,elect->eval[off]);
        if (elect->occ) occ+=elect->occ[off];
      }
      cross=0;
      if ((elect->e_fermi)&&(emin<*elect->e_fermi)&&(emax>*elect->e_fermi))
	cross=1;
      if (cross)
	fprintf(stderr,"*  %4d  %1d %14f  %14f  %14f\n",b+1,ns,occ/kp->n,
		emin*scale,emax*scale);
      else
	if (debug)
	  fprintf(stderr,"   %4d  %1d %14f  %14f  %14f\n",b+1,ns,occ/kp->n,
		  emin*scale,emax*scale);
    }
  }
}

int band_cmp(const void *ptr1, const void *ptr2){
  struct band { double energy; double weight;};

  if (((struct band*)ptr1)->energy<((struct band*)ptr2)->energy) return -1;
  if (((struct band*)ptr1)->energy>((struct band*)ptr2)->energy) return 1;
  return 0;
  
}

double calc_efermi(struct es *elect, struct kpts *kpt, double nel){
  int k,ns,b,nb,i;
  double wt,total,fill,ef;
  struct band { double energy; double weight;} *bands;

  if (!elect->eval){
    fprintf(stderr,"Unable to calculate Fermi level without eigenvalues\n");
    return NAN;
  }
  
  ef=0;
  nb=kpt->n*elect->nbspins*elect->nbands;
  bands=malloc(nb*sizeof(struct band));
  if (!bands) error_exit("malloc error in calc_fermi");

  i=0;
  for(k=0;k<kpt->n;k++){
    wt=kpt->kpts[k].wt;
    if (elect->nbspins==1) wt*=2;
    for(ns=0;ns<elect->nbspins;ns++){
      for(b=0;b<elect->nbands;b++){
	bands[i].energy=elect->eval[i];
	bands[i].weight=wt;
	i++;
      }
    }
  }

  qsort(bands,nb,sizeof(struct band),band_cmp);

  i=0;
  total=0;
  while((i<nb)&&(total+bands[i].weight<nel)) total+=bands[i++].weight;

  if (i==nb){
    fprintf(stderr,"No empty bands\n");
    ef=bands[nb-1].energy;
  }
  else{
    fill=(nel-total)/bands[i].weight;
    if (fill<sqrt(tol)){
      ef=0.5*(bands[i-1].energy+bands[i].energy);
    }
    else if (fill>1-sqrt(tol)){
      i++;
      if (i==nb){
	fprintf(stderr,"No empty bands\n");
	ef=bands[nb-1].energy;
      }
      else
	ef=0.5*(bands[i-1].energy+bands[i].energy);
    }
    else{
      ef=bands[i-1].energy;
      if (debug) fprintf(stderr,"Band at Fermi level partially filled\n");
      if (debug>1) fprintf(stderr,"   occupancy=%lf\n",fill);
    }
  }

  return ef;
  
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

void print_energy(double e){
  if (flags&AU)
    fprintf(stderr,"%lf Ha",e/H_eV);
  else
    fprintf(stderr,"%lf eV",e);
}

void print_gap(struct es *elect, struct kpts *kpt, struct contents *m){
  int ns,i,k,nel,ivb,ind;
  int vbm_ind,cbm_ind,dgap_ind;
  double chg,vbm,cbm,dgap;

  chg=0;
  vbm_ind=cbm_ind=dgap_ind=-1;
  
  if (!elect->eval){
    fprintf(stderr,"Cannot print band gap without eigenvalues\n");
    return;
  }

  if (elect->nel==0){
    for(i=0;i<m->n;i++)
      chg+=m->atoms[i].chg;
    if ((chg)&&(elect->charge)) chg-=*elect->charge;
    if (chg>0)
      fprintf(stderr,"Assuming nel=%lf from total ionic charge\n",chg);
  }
  else chg=elect->nel;
  
  nel=(int)(chg+0.5);

  if (nel==0){
    fprintf(stderr,"Unable to determine number of electrons\n");
    return;
  }

  if (fabs(chg-nel)>tol){
    fprintf(stderr,"Non integer number of electrons\n");
    return;
  }

  if (nel&1){
    fprintf(stderr,"Odd number of electrons, so metallic\n");
    return;
  }

  if (debug) fprintf(stderr,"Nel=%d\n",nel);

  ivb=nel/2;

  if (ivb+1>elect->nbands){
    fprintf(stderr,"No conduction bands calculated\n");
    return;
  }
  
  for(ns=0;ns<elect->nbspins;ns++){
  
    /* Set to huge values */
    vbm=-1e20;
    cbm=1e20;
    dgap=1e20;
    
    for(k=0;k<kpt->n;k++){
      ind=k*(elect->nbspins*elect->nbands)+ns*elect->nbands+ivb-1;
      if (elect->eval[ind]>vbm){
	vbm=elect->eval[ind];
	vbm_ind=k;
      }
      if (elect->eval[ind+1]<cbm){
	cbm=elect->eval[ind+1];
	cbm_ind=k;
      }
      if (elect->eval[ind+1]-elect->eval[ind]<dgap){
	dgap=elect->eval[ind+1]-elect->eval[ind];
	dgap_ind=k;
      }
    }

    if (elect->nbspins==2){
      if (ns==0) fprintf(stderr,"\nSpin up:\n");
      else fprintf(stderr,"\nSpin down:\n");
    }

    fprintf(stderr,"\nValence band maximum:    ");
    print_energy(vbm);
    fprintf(stderr," at k=(%f,%f,%f)\n",kpt->kpts[vbm_ind].frac[0],
	    kpt->kpts[vbm_ind].frac[1],kpt->kpts[vbm_ind].frac[2]);


    if ((elect->e_fermi)&&(vbm>*elect->e_fermi))
      fprintf(stderr,"*** Problem: VMB > E_Fermi\n");
    
    fprintf(stderr,"Conduction band minimum: ");
    print_energy(cbm);
    fprintf(stderr," at k=(%f,%f,%f)\n",kpt->kpts[cbm_ind].frac[0],
	    kpt->kpts[cbm_ind].frac[1],kpt->kpts[cbm_ind].frac[2]);

    if ((elect->e_fermi)&&(cbm<*elect->e_fermi))
      fprintf(stderr,"*** Problem: CBM < E_Fermi\n");

    if (cbm-vbm>0){
      fprintf(stderr,"Direct gap:              ");
      print_energy(dgap);
    
      fprintf(stderr," at k=(%f,%f,%f)\n",kpt->kpts[dgap_ind].frac[0],
	      kpt->kpts[dgap_ind].frac[1],kpt->kpts[dgap_ind].frac[2]);

      if (vbm_ind!=cbm_ind){
	fprintf(stderr,"Indirect gap:            ");
	print_energy(cbm-vbm);
	fprintf(stderr,"\n");
      }
    }
    else fprintf(stderr,"System is metallic\n");
  }
}

