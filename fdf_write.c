/* Various Siesta writers:
 * fdf_write -- write fdf file, also .RHO file if grid exists
 * xv_write  -- ditto, XV file
 * rho_write -- write RHO file, single component, or two if
 *                 CHDEN and SPINDEN both set in flags
 */


/* Copyright (c) 2007-2019 MJ Rutter 
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

#include "c2xsf.h"

void rho_write(char *oldname, struct unit_cell *c, struct grid *g);

void fdf_write(FILE* fdf, char* filename, struct unit_cell *c,
               struct contents *m, struct kpts *kp,
               struct grid *g, struct es *e){
  int i,j,nspec,*natomsp,*spatno;
  double dtmp;

  fprintf(fdf,"# fdf file created by c2x\n");

  if (m->title) fprintf(fdf,"\nSystemName %s\n\n",m->title);

  if (dict_get(m->dict,"Siesta_SystemLabel"))
    fprintf(fdf,"SystemLabel %s\n\n",
            (char*)dict_get(m->dict,"Siesta_SystemLabel"));
  
  fprintf(fdf,"LatticeConstant  1.0 Ang\n");

  fprintf(fdf,"%%block LatticeVectors\n");
  for(i=0;i<3;i++)fprintf(fdf,"%f %f %f\n",c->basis[i][0],
        c->basis[i][1],c->basis[i][2]);
  fprintf(fdf,"%%endblock LatticeVectors\n\n");

  /* Now we need to know the number of species.
   * It must be fewer than the number of atoms...
   */

  nspec=0;
  natomsp=malloc(m->n*sizeof(int));
  if (!natomsp) error_exit("Malloc error in fdf_write");
  spatno=malloc(m->n*sizeof(int));
  if (!spatno) error_exit("Malloc error in fdf_write");


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

  fprintf(fdf,"NumberOfAtoms    %d\n",m->n);
  fprintf(fdf,"NumberOfSpecies  %d\n",nspec);

  fprintf(fdf,"%%block ChemicalSpeciesLabel\n");
  for(i=0;i<nspec;i++)
    fprintf(fdf,"%d  %d   %s\n",i+1,spatno[i],atno2sym(spatno[i]));
  fprintf(fdf,"%%endblock ChemicalSpeciesLabel\n\n");

  fprintf(fdf,"AtomicCoordinatesFormat Fractional\n");
    /* Write out atoms, sorted by species */
  fprintf(fdf,"%%block AtomicCoordinatesAndAtomicSpecies\n");
  for(i=0;i<nspec;i++)
    for(j=0;j<m->n;j++)
      if (m->atoms[j].atno==spatno[i])
        fprintf(fdf," %f %f %f %d\n",m->atoms[j].frac[0],
                m->atoms[j].frac[1],m->atoms[j].frac[2],i+1);
  fprintf(fdf,"%%endblock AtomicCoordinatesAndAtomicSpecies\n");

  /* Optional colinear spins */
  
  j=0;
  for(i=0;i<m->n;i++)
    if (m->atoms[i].spin) {j=1;break;}

  if (j){
    fprintf(fdf,"\nSpinPolarized .true.\n");
    fprintf(fdf,"%%block DMInitSpin\n");
    for(i=0;i<m->n;i++)
      if (m->atoms[i].spin) fprintf(fdf,"  %d  %f\n",i+1,m->atoms[i].spin);
    fprintf(fdf,"%%endblock DMInitSpin\n");
  }                           
      
  
  if (e->etol){
    fprintf(fdf,"\nDM.Require.Energy.Convergence .true.\n");
    fprintf(fdf,"DM.Energy.Tolerance           %g eV\n",e->etol*m->n);
  }

  if ((e->charge)&&(*e->charge!=0.0))
    fprintf(fdf,"NetCharge %f\n",*e->charge);

  /* MP k-points */
  if ((kp)&&(kp->mp)){
    fprintf(fdf,"\n%%block kgrid_Monkhorst_Pack\n");
    fprintf(fdf,"%d 0 0 ",kp->mp->grid[0]);
    dtmp=kp->mp->disp[0]*kp->mp->grid[0];
    if ((kp->mp->grid[0]&1)==0) dtmp+=0.5;
    dtmp=fmod(dtmp,1.0);
    fprintf(fdf,"%f\n",dtmp);
    fprintf(fdf,"0 %d 0 ",kp->mp->grid[1]);
    dtmp=kp->mp->disp[1]*kp->mp->grid[1];
    if ((kp->mp->grid[1]&1)==0) dtmp+=0.5;
    dtmp=fmod(dtmp,1.0);
    fprintf(fdf,"%f\n",dtmp);
    fprintf(fdf,"0 0 %d ",kp->mp->grid[2]);
    dtmp=kp->mp->disp[2]*kp->mp->grid[2];
    if ((kp->mp->grid[2]&1)==0) dtmp+=0.5;
    dtmp=fmod(dtmp,1.0);
    fprintf(fdf,"%f\n",dtmp);
    fprintf(fdf,"%%endblock kgrid_Monkhorst_Pack\n");
  }

  if (dict_get(m->dict,"Siesta_misc"))
    fprintf(fdf,"\n%s",(char*)dict_get(m->dict,"Siesta_misc"));
  
  /* And now write density */

  if ((g==NULL)||(g->data==0)) return;

  rho_write(filename,c,g);
}


/* XV files are in Bohr
   3 x (cell vector, cell vector velocity)
   natoms
   natoms x (spec no, atomic no, abs coords, velocity)
 */

void xv_write(FILE* outfile, char *filename, struct unit_cell *c,
              struct contents *m, struct grid *g){
  int i,j;
  int *atomspec,*specatno,nspec;
  char *fmt;

  if (flags&HIPREC)
    fmt="  %19.14f %19.14f %19.14f   0.00 0.00 0.00\n";
  else
    fmt="    %17.9f %17.9f %17.9f"
      "          0.000000000       0.000000000       0.000000000\n";

  
  for(i=0;i<3;i++)
    fprintf(outfile,fmt,
            c->basis[i][0]/BOHR,c->basis[i][1]/BOHR,c->basis[i][2]/BOHR);

  fprintf(outfile,"    %8d\n",m->n);

  /* Create species number for each atom */

  atomspec=malloc(m->n*sizeof(int));
  specatno=malloc(m->n*sizeof(int));
  if ((!atomspec)||(!specatno)) error_exit("Malloc error in xv_write");

  nspec=0;
  for(i=0;i<m->n;i++){
    for(j=0;j<nspec;j++)
      if (m->atoms[i].atno==specatno[j]) break;
    if (j<nspec)
      atomspec[i]=j+1;
    else{
      atomspec[i]=nspec+1;
      specatno[nspec]=m->atoms[i].atno;
      nspec++;
    }
  }

  if (flags&HIPREC)
    fmt=" %2d %5d %19.14f %19.14f %19.14f   0.00 0.00 0.00\n";
  else
    fmt=" %2d %5d %17.9f %17.9f %17.9f"
      "          0.000000000       0.000000000       0.000000000\n";
  
  for(i=0;i<m->n;i++)
    fprintf(outfile,fmt,
            atomspec[i],m->atoms[i].atno,
            m->atoms[i].abs[0]/BOHR,m->atoms[i].abs[1]/BOHR,
            m->atoms[i].abs[2]/BOHR);

  free(specatno);
  free(atomspec);

  /* And now write density */

  if ((g==NULL)||(g->data==0)) return;

  rho_write(filename,c,g);
  
}

/* RHO files are in Fortran binary format
 *
 * Record 1: basis (Bohr)
 * Record 2: grid (3*int) and nspins (int)
 * do spin=1,nspins
 *  do iz=1,grid(3)
 *    do iy=1,grid(2)
 *     Record: data(1:grid(1)) as single precision float
 *
 * The default unit is e per cubic Bohr
 *
 * If spin is present, the spin up density and spin down density are
 *   written, not the sum and difference (which is what c2x stores)
 */

void rho_write(char *oldname, struct unit_cell *c, struct grid *g){
  char *cptr,*filename;
  int i,j,k,len;
  float *fptr,scale;
  double *dptr1,*dptr2,basis[3][3];
  FILE *bin;

  if (oldname==NULL) return;

  filename=malloc(strlen(oldname)+5);
  if (!filename) error_exit("Malloc error for filename");
  strcpy(filename,oldname);

  cptr=filename+strlen(filename);
  while ((cptr>filename)&&(*cptr!='.')) cptr--;
  if (*cptr!='.') return;

  cptr[1]='R';
  cptr[2]='H';
  cptr[3]='O';

  bin=fopen(filename,"wb");
  if(!bin){
    fprintf(stderr,"Error: unable to open %s for writing\n",filename);
    return;
  }

  scale=1;
  if ((!(flags&RAW))&&(!(flags&AU))) scale=BOHR*BOHR*BOHR;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      basis[i][j]=c->basis[i][j]/BOHR;
  
  len=9*sizeof(double);
  fwrite(&len,sizeof(int),1,bin);
  fwrite(basis,sizeof(double),9,bin);
  fwrite(&len,sizeof(int),1,bin);

  len=4*sizeof(int);
  fwrite(&len,sizeof(int),1,bin);
  fwrite(g->size,sizeof(int),3,bin);
  i=1;
  if ((flags&CHDEN)&&(flags&SPINDEN)&&(g->next->data)) i=2;
  fwrite(&i,sizeof(int),1,bin);
  fwrite(&len,sizeof(int),1,bin);

  if(!(fptr=malloc(g->size[0]*sizeof(float))))
      error_exit("Malloc error in fdf_write\n");

  len=g->size[0]*sizeof(float);

  if ((flags&CHDEN)&&(flags&SPINDEN)&&(g->next->data)){
    /* Need to write spin up as 0.5*(sum+difference) and
     * spin down as 0.5*(sum-difference)
     */
    if (debug) fprintf(stderr,"Also writing spin\n");
    scale*=0.5;
    for (k=0;k<g->size[2];k++){
      for (j=0;j<g->size[1];j++){
        dptr1=g->data+k+j*g->size[2];
        dptr2=g->next->data+k+j*g->size[2];
        for (i=0;i<g->size[0];i++)
          fptr[i]=(*(dptr1+i*g->size[1]*g->size[2])
                   +*(dptr2+i*g->size[1]*g->size[2]))*scale;
        fwrite(&len,sizeof(int),1,bin);
        fwrite(fptr,sizeof(float),g->size[0],bin);
        fwrite(&len,sizeof(int),1,bin);
      }
    }
    for (k=0;k<g->size[2];k++){
      for (j=0;j<g->size[1];j++){
        dptr1=g->data+k+j*g->size[2];
        dptr2=g->next->data+k+j*g->size[2];
        for (i=0;i<g->size[0];i++)
          fptr[i]=(*(dptr1+i*g->size[1]*g->size[2])
                   -*(dptr2+i*g->size[1]*g->size[2]))*scale;
        fwrite(&len,sizeof(int),1,bin);
        fwrite(fptr,sizeof(float),g->size[0],bin);
        fwrite(&len,sizeof(int),1,bin);
      }
    }
  }
  else{
    for (k=0;k<g->size[2];k++){
      for (j=0;j<g->size[1];j++){
        dptr1=g->data+k+j*g->size[2];
        for (i=0;i<g->size[0];i++)
          fptr[i]=*(dptr1+i*g->size[1]*g->size[2])*scale;
        fwrite(&len,sizeof(int),1,bin);
        fwrite(fptr,sizeof(float),g->size[0],bin);
        fwrite(&len,sizeof(int),1,bin);
      }
    }
  }
  
  fprintf(stderr,"Also wrote %s\n",filename);

  free(fptr);
  free(filename);
  
}



