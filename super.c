/* Conversion to arbitrary, specified supercells, with grid interpolation */

/* Copyright (c) 2007 MJ Rutter 
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
#include<math.h>

#include "c2xsf.h"

void vcross(double a[3],double b[3],double c[3]);


void to235(int *i);
void grid_interp(struct unit_cell *c, struct grid *grid,int old_fft[3],
                    double old_basis[3][3],double old_recip[3][3]);
/* From basis.c */
int is_rhs(double b[3][3]);
/* From ksym.c */
void sym_basis(struct symmetry *s, struct unit_cell *c);
/* From sort_atoms.c */
int atom_sort(const void *a, const void *b);

int super(struct unit_cell *c, struct contents *mtf,
           double new_basis[3][3], struct kpts *kp, struct symmetry *s, 
           struct grid *gptr, int sflags){
  int i,j,k,l,m,na,at,old_fft[3],same,rhs,quiet,kpts_only;
  double old_basis[3][3],old_recip[3][3],old_vol,dtmp;
  double new_in_old[3][3],abc[6],vtmp[3];
  struct atom *old_atoms;
  int *atom_ctr;
  int old_natoms;
  int scan_min[3],scan_max[3];
  double corner, fscan_min[3],fscan_max[3],disp[3],new_abs[3],new_frac[3];
  double fft_res;

  rhs=sflags&1;
  quiet=sflags&2;
  kpts_only=sflags&4;

  /* Correct new basis if not rhs and rhs wanted */
  /* Correction is exchange of 2nd and 3rd vectors */
  if (rhs){
    if (!is_rhs(new_basis)){
      for(i=0;i<3;i++){
        dtmp=new_basis[1][i];
        new_basis[1][i]=new_basis[2][i];
        new_basis[2][i]=dtmp;
      }
    }
  }

  
  same=1;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      if (!aeq(c->basis[i][j],new_basis[i][j])) same=0;

  if (same) return 0;

  /* Try to cope with k-points, else throw them away */

  if (kp && kp->mp){ /* Attempt conversion if axes simply scaled */
    int new_mp_grid[3],success;
    double new_mp_disp[3];
    double vtmp[3];
    
    for(i=0;i<3;i++){
      new_mp_grid[i]=0;
      vcross(c->basis[i],new_basis[i],vtmp);
      if (aeq(vmod2(vtmp),0)){ /* New and old vector are parallel */
        dtmp=sqrt(vmod2(c->basis[i])/vmod2(new_basis[i]));
        new_mp_grid[i]=dtmp*kp->mp->grid[i];
        if ((int)((new_mp_grid[i]/dtmp)+0.5)!=kp->mp->grid[i])
          new_mp_grid[i]=0;
      
        new_mp_disp[i]=kp->mp->disp[i];
        if ((kp->mp->grid[i]&1)==0)
          new_mp_disp[i]+=1.0/(2*kp->mp->grid[i]); /* For thus says Castep */
        new_mp_disp[i]/=dtmp;
        if ((new_mp_grid[i]&1)==0)
          new_mp_disp[i]+=1.0/(2*new_mp_grid[i]); /* For thus says Castep */
        new_mp_disp[i]=fmod(new_mp_disp[i],1.0/new_mp_grid[i]);
      }
    }
    success=1;
    for(i=0;i<3;i++) if (new_mp_grid[i]==0) success=0;

    if (success){
      for(i=0;i<3;i++)
        kp->mp->grid[i]=new_mp_grid[i];
      for(i=0;i<3;i++)
        kp->mp->disp[i]=new_mp_disp[i];
      if (debug>1) fprintf(stderr,"New MP grid %dx%dx%d offset (%f,%f,%f)\n",
                           kp->mp->grid[0],kp->mp->grid[1],kp->mp->grid[2],
                           kp->mp->disp[0],kp->mp->disp[1],kp->mp->disp[2]);
    }
    else{
      if (debug>1) fprintf(stderr,"Discarding MP grid\n");
      free(kp->mp);
      kp->mp=NULL;
    }
  }
  
  if (kp && kp->n){
    struct unit_cell rcell,new_cell;
    struct contents rmtf;

    rcell.basis=malloc(9*sizeof(double));
    for(i=0;i<3;i++)
      for(j=0;j<3;j++){
	rcell.basis[i][j]=c->recip[i][j];
        rcell.recip[i][j]=c->basis[i][j];
      }
    rcell.vol=1.0/c->vol;

    rmtf.n=kp->n;
    rmtf.atoms=kp->kpts;
    addabs(rmtf.atoms,rmtf.n,rcell.basis);

    new_cell.basis=malloc(9*sizeof(double));
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	new_cell.basis[i][j]=new_basis[i][j];
    real2rec(&new_cell);
    
    if (super(&rcell,&rmtf,new_cell.recip,NULL,NULL,NULL,2)==0){

      if (debug>2){
        fprintf(stderr,"New kpoints:\n");
        for(i=0;i<rmtf.n;i++)
          fprintf(stderr,"%f %f %f\n",rmtf.atoms[i].frac[0],
    	      rmtf.atoms[i].frac[1],rmtf.atoms[i].frac[2]);
      }

      kp->n=rmtf.n;
      kp->kpts=rmtf.atoms;
      for(i=0;i<rmtf.n;i++) kp->kpts[i].wt*=fabs(new_cell.vol/c->vol);
    }
    else{
      if (debug)
	fprintf(stderr,"Warning: unable to convert k points to new cell\n");
      kp->n=0;
      if (kp->kpts)  {free(kp->kpts);kp->kpts=NULL;}
    }
    if (kp->mp)    {free(kp->mp);kp->mp=NULL;}
    free(new_cell.basis);
    free(rcell.basis);
  }

  if (kpts_only) return 0;
  
  /* Save old cell */

  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      old_basis[i][j]=c->basis[i][j];
      old_recip[i][j]=c->recip[i][j];
    }
  }

  old_atoms=mtf->atoms;
  old_natoms=mtf->n;
  old_vol=c->vol;

  /* Make new cell */

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[i][j]=new_basis[i][j];

  real2rec(c); /* This will also update c->vol */

  c->vol=fabs(c->vol);


  if (!quiet){
    if (debug) fprintf(stderr,"New cell volume %f (%g times old)\n",
                       c->vol,c->vol/old_vol);
    if (debug>1){
      fprintf(stderr,"New basis set\n");
      for(i=0;i<=2;i++)
        fprintf(stderr,"%f %f %f\n",c->basis[i][0],
		c->basis[i][1],c->basis[i][2]);
    }

    if (debug>2){
      fprintf(stderr,"New reciprocal basis set\n");
      for(i=0;i<=2;i++)
        fprintf(stderr,"%f %f %f\n",c->recip[i][0],
		c->recip[i][1],c->recip[i][2]);
    }
  }

  /* Deal with symmetry operations */

  if (s) sym_basis(s,c);

  /* Worry about atoms */

  mtf->n=old_natoms*c->vol/old_vol+0.5;
  if (fabs(mtf->n-old_natoms*c->vol/old_vol)>min(0.05,tol*mtf->n)){
    if (quiet) return 1;
    fprintf(stderr,"Impossible cell transformation leads to %f atoms\n",
            old_natoms*c->vol/old_vol);
    exit(1);
  }
  
  if(!(mtf->atoms=malloc(mtf->n*sizeof(struct atom))))
    error_exit("Malloc error in super");
  if(!(atom_ctr=malloc(mtf->n*sizeof(int))))
    error_exit("Malloc error in super");

  /* Find extent of new corners in old cell */

  for(i=0;i<3;i++) fscan_min[i]=fscan_max[i]=0;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      new_in_old[i][j]=c->basis[i][0]*old_recip[j][0]+
                       c->basis[i][1]*old_recip[j][1]+
                       c->basis[i][2]*old_recip[j][2];

  for(i=0;i<2;i++)
    for(j=0;j<2;j++)
      for(k=0;k<2;k++)
        for(l=0;l<3;l++){
          corner=i*new_in_old[0][l]+j*new_in_old[1][l]+k*new_in_old[2][l];
          if (corner<fscan_min[l]) fscan_min[l]=corner;
          if (corner>fscan_max[l]) fscan_max[l]=corner;
        }

  for(i=0;i<3;i++){
    scan_max[i]=fscan_max[i]+0.01;
    scan_min[i]=fscan_min[i]-0.99;
    if (debug>2) fprintf(stderr,"i=%d scan_min=%d scan_max=%d\n",i,
                         scan_min[i],scan_max[i]);
  }

  /* Now loop over required array of cells */

  na=0;
  for(i=scan_min[0];i<=scan_max[0];i++){
    for(j=scan_min[1];j<=scan_max[1];j++){
      for(k=scan_min[2];k<=scan_max[2];k++){
        for(l=0;l<3;l++)
          disp[l]=i*old_basis[0][l]+j*old_basis[1][l]+k*old_basis[2][l];
        for(at=0;at<old_natoms;at++){
          for(l=0;l<3;l++) new_abs[l]=old_atoms[at].abs[l]+disp[l];
  /* reduce this atom to new basis set */
          for(l=0;l<3;l++){
            new_frac[l]=new_abs[0]*c->recip[l][0]+
                        new_abs[1]*c->recip[l][1]+
                        new_abs[2]*c->recip[l][2];
            new_frac[l]=fmod(new_frac[l],1.0);
	    if (new_frac[l]<0.0) new_frac[l]+=1.0;
            if (new_frac[l]>1.0-1e-9) new_frac[l]=0.0;
          }

  /* add to our list if we have not yet seen it */
          for(l=0;l<na;l++)
            if((dist(mtf->atoms[l].frac[0],new_frac[0])*
                dist(mtf->atoms[l].frac[0],new_frac[0])+
                dist(mtf->atoms[l].frac[1],new_frac[1])*
                dist(mtf->atoms[l].frac[1],new_frac[1])+
                dist(mtf->atoms[l].frac[2],new_frac[2])*
                dist(mtf->atoms[l].frac[2],new_frac[2]))<tol*tol) break;
          if (l==na){
            if (na>=mtf->n){
	      if (quiet) {free(atom_ctr); return 1;}
              if (debug<2) error_exit("Too many new atoms found in super.c");
              else{
                fprintf(stderr,"Too many atoms found in super.c\n");
                for (l=0;l<na;l++)
                  fprintf(stderr,"%3d %3d %f %f %f\n",l,mtf->atoms[l].atno,
                          mtf->atoms[l].frac[0],mtf->atoms[l].frac[1],
			  mtf->atoms[l].frac[2]);
                fprintf(stderr,"%d %f %f %f\n",old_atoms[at].atno,
                    new_frac[0],new_frac[1],new_frac[2]);
                exit(1);
              }
            }
	    init_atoms(mtf->atoms+na,1);
            for(l=0;l<3;l++){
              mtf->atoms[na].frac[l]=new_frac[l];
              mtf->atoms[na].abs[l]=new_frac[0]*c->basis[0][l]+
		                     new_frac[1]*c->basis[1][l]+
		                     new_frac[2]*c->basis[2][l];
            }           
            mtf->atoms[na].atno=old_atoms[at].atno;
            mtf->atoms[na].wt=old_atoms[at].wt;
            mtf->atoms[na].spin=old_atoms[at].spin;
            mtf->atoms[na].chg=old_atoms[at].chg;
            mtf->atoms[na].site_chg=old_atoms[at].site_chg;
            mtf->atoms[na].label=old_atoms[at].label;
	    for(l=0;l<3;l++)
	      mtf->atoms[na].force[l]=old_atoms[at].force[l];
	    atom_ctr[na]=1;
            na++;
          }
	  else{ /* (l!=na) -- we have seen this atom before, so average */
	    for(m=0;m<3;m++){
              /* Need new fractional co-ords to be same as old,
               * i.e. not -0.2.. and 0.8..., but 0.8... and 0.8...
               */
              if (new_frac[m]-mtf->atoms[l].frac[m]>0.5) new_frac[m]-=1.0;
              if (new_frac[m]-mtf->atoms[l].frac[m]<-0.5) new_frac[m]+=1.0;
	      mtf->atoms[l].frac[m]=(atom_ctr[l]*mtf->atoms[l].frac[m]+
				      new_frac[m])/(atom_ctr[l]+1);
	      mtf->atoms[l].abs[m]=new_frac[0]*c->basis[0][m]+
		              new_frac[1]*c->basis[1][m]+
                              new_frac[2]*c->basis[2][m];
	    }
	    atom_ctr[l]++;
	  }
        }
      }
    }
  }
  free(atom_ctr);

  sort_atoms(mtf,1);
  
  if ((!quiet)&&(debug>1))
    fprintf(stderr,"New cell: na=%d, mtf->n=%d\n",na,mtf->n);

  if (na!=mtf->n){
    if (quiet) return 1;
    fprintf(stderr,"Surprise in super.c. Expected %d atoms, found %d\n",
            mtf->n,na);
    mtf->n=na;
  }

  free(old_atoms);
  
  /* Worry about grids */

  while((gptr)&&(gptr->data)){
    fft_res=0;
    /* Find maximal grid point density */
    for(i=0;i<3;i++){
      dtmp=gptr->size[i]/sqrt(old_basis[i][0]*old_basis[i][0]+
                              old_basis[i][1]*old_basis[i][1]+
                              old_basis[i][2]*old_basis[i][2]);
      if (dtmp>fft_res) fft_res=dtmp;
    }
    cart2abc(c,NULL,abc,NULL,0);
    for(i=0;i<3;i++){
      old_fft[i]=gptr->size[i];
      gptr->size[i]=abc[i]*fft_res;
      /* But if new axis simply a multiple of old, make new grid size
       * that multiple of old too */
      vcross(old_basis[i],c->basis[i],vtmp);
      if (vmod2(vtmp)<1e-6*vmod2(old_basis[i])){
        dtmp=sqrt(vmod2(c->basis[i])/vmod2(old_basis[i]));
        gptr->size[i]=old_fft[i]*dtmp+0.5;
      }
      to235(gptr->size+i);
    }
    if (debug) fprintf(stderr,"New FFT grid is %d %d %d\n",
                       gptr->size[0],gptr->size[1],gptr->size[2]);
    grid_interp(c,gptr,old_fft,old_basis,old_recip);
    gptr=gptr->next;
  }
  return 0;
}

void to235(int *i){
/* Force i to next highest int whose only factors are 2, 3 and 5 */
  int tmp;

  if ((*i)<1) *i=2;
  tmp=*i;
  while (tmp%2==0) tmp/=2;
  while (tmp%3==0) tmp/=3;
  while (tmp%5==0) tmp/=5;

  if (tmp==1) return;

  (*i)++;
  to235(i);
}

void grid_interp(struct unit_cell *c, struct grid *grid, int old_fft[3],
                    double old_basis[3][3], double old_recip[3][3]){
  int i,j,k,l,ii,jj,kk,pt[3][2];
  double pabs[3],pfrac[3],ifrac[3][2],*gnew;
  double di,dj,dk,x,sum;

  sum=0;

  if (debug>1) fprintf(stderr,"Moving from %dx%dx%d grid to %dx%dx%d\n",
		       old_fft[0],old_fft[1],old_fft[2],
		       grid->size[0],grid->size[1],grid->size[2]);

  if (!(gnew=malloc(sizeof(double)*grid->size[0]*grid->size[1]*grid->size[2])))
     error_exit("Malloc error in grid_interp");

  for(i=0;i<grid->size[0];i++){
    di=(double)i/grid->size[0];
    for(j=0;j<grid->size[1];j++){
      dj=(double)j/grid->size[1];
      for(k=0;k<grid->size[2];k++){
        dk=(double)k/grid->size[2];
        for(l=0;l<3;l++)
          pabs[l]=di*c->basis[0][l]+dj*c->basis[1][l]+dk*c->basis[2][l];
        /* convert to old cell.basis */
        for(l=0;l<3;l++){
          pfrac[l]=pabs[0]*old_recip[l][0]+pabs[1]*old_recip[l][1]+
                   pabs[2]*old_recip[l][2];
        /* reduce to unit cell */
          pfrac[l]=fmod(pfrac[l],1.0);
	  if (pfrac[l]<0.0) pfrac[l]+=1.0;
        /* convert to old grid */
          pfrac[l]=pfrac[l]*old_fft[l];
        /* and to int and fractional parts */
          pt[l][0]=pfrac[l];
          ifrac[l][1]=pfrac[l]-pt[l][0];
	  if (pt[l][0]==old_fft[l]) pt[l][0]=0;
	  if (pt[l][0]>old_fft[l]) error_exit("Impossible in grid_interp");
        /* and the other side of the grid "cube" */
          pt[l][1]=(pt[l][0]+1)%old_fft[l];
          ifrac[l][0]=1-ifrac[l][1];
        }
        /* Trilinear interpolation */
#define OX(x,y,z) grid->data[((x)*old_fft[1]+(y))*old_fft[2]+(z)]
        x=0;
        for(ii=0;ii<2;ii++)
          for(jj=0;jj<2;jj++)
            for(kk=0;kk<2;kk++)
              x+=ifrac[0][ii]*ifrac[1][jj]*ifrac[2][kk]*
                   OX(pt[0][ii],pt[1][jj],pt[2][kk]);
        gnew[((i*grid->size[1])+j)*grid->size[2]+k]=x;
	if (fabs(x)>99) {
	  fprintf(stderr,"%lf\n",x);
	  for(ii=0;ii<3;ii++) fprintf(stderr,"ifrac[0][%d]=%lf\n",ii,ifrac[0][i]);
	  for(ii=0;ii<3;ii++)
	    fprintf(stderr,"pt[i][0]=%5d  pt[i][1]=%5d\n",pt[ii][0],pt[ii][1]);  
	}
        sum+=x;
      }
    }
  }

  free(grid->data);
  grid->data=gnew;

  if (debug>1) fprintf(stderr,"On new grid sum=%g int=%g\n",sum,
                     sum*c->vol/(grid->size[0]*grid->size[1]*grid->size[2]));
}

/* Independent code for the simple case of new vectors being integer
 * multiples of old
 */

void simple_super(struct unit_cell *c, struct contents *m,
           int expand[3], struct kpts *kp, struct symmetry *s, 
           struct grid *gptr){
  double new_basis[3][3],old_basis[3][3],*dptr;
  int new_fft[3],old_fft[3];
  int i,j,k,p,q,old_n;
  int ii,jj,kk;
  int old_off,new_off;
  struct atom *old_atoms;
  
  if ((expand[0]==1)&&(expand[1]==1)&&(expand[2]==1)) return;

  for(i=0;i<3;i++)
    if(expand[i]<1) error_exit("Invalid cell tranformation");

  if (debug>1)
    fprintf(stderr,"Cell expansion by %dx%dx%d\n",
            expand[0],expand[1],expand[2]);

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      old_basis[i][j]=c->basis[i][j];
  old_atoms=m->atoms;
  old_n=m->n;

  /* k-points */

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      new_basis[i][j]=expand[i]*c->basis[i][j];

  if (debug){
    fprintf(stderr,"New basis:\n");
    print_basis(new_basis);
  }
  
  super(c,m,new_basis,kp,s,gptr,4);

  /* New basis */

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[i][j]=new_basis[i][j];
  real2rec(c);
  c->vol=fabs(c->vol);
  
  /* Deal with symmetry operations */

  if (s) sym_basis(s,c);

  /* Atoms */
  
  m->n*=expand[0]*expand[1]*expand[2];
  m->atoms=malloc(m->n*sizeof(struct atom));
  if (!m->atoms) error_exit("Malloc error for expanded atoms");

  q=0;
  for(i=0;i<expand[0];i++){
    for(j=0;j<expand[1];j++){
      for(k=0;k<expand[2];k++){
	for(p=0;p<old_n;p++){
	  m->atoms[q]=old_atoms[p];
	  for(ii=0;ii<3;ii++)
	    m->atoms[q].abs[ii]+=i*old_basis[0][ii];
	  for(ii=0;ii<3;ii++)
	    m->atoms[q].abs[ii]+=j*old_basis[1][ii];
	  for(ii=0;ii<3;ii++)
	    m->atoms[q].abs[ii]+=k*old_basis[2][ii];
	  q++;
	}
      }
    }
  }

  addfrac(m->atoms,m->n,c->recip);

  free(old_atoms);

  /* Worry about grids */

  while((gptr)&&(gptr->data)){  
    for(i=0;i<3;i++)
      old_fft[i]=gptr->size[i];
    for(i=0;i<3;i++)
      new_fft[i]=expand[i]*old_fft[i];
    dptr=malloc(new_fft[0]*new_fft[1]*new_fft[2]*sizeof(double));
    if (!dptr) error_exit("Malloc error for new grid");
    for(i=0;i<new_fft[0];i++){
      ii=i%old_fft[0];
      for(j=0;j<new_fft[1];j++){
	jj=j%old_fft[1];
	for(k=0;k<new_fft[2];k++){
	  kk=k%old_fft[2];
	  old_off=kk+old_fft[2]*(jj+ii*old_fft[1]);
	  new_off=k+new_fft[2]*(j+i*new_fft[1]);
	  dptr[new_off]=gptr->data[old_off];
	}
      }
    }
    free(gptr->data);
    gptr->data=dptr;
    for(i=0;i<3;i++)
      gptr->size[i]=new_fft[i];
    gptr=gptr->next;
  }
}
