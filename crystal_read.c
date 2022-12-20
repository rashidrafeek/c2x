/* A reader for Crystal's geometry files */

/* Copyright (c) 2018 MJ Rutter 
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
#include<math.h>

#include "c2xsf.h"
#ifdef SPGLIB
#include "spglib.h"
#endif

struct sym_op *sym_frac2abs(int spg_rot[][3][3],double spg_tr[][3],
                            struct unit_cell *c,int nsym);


#define LINE_SIZE 2049

void crystal_read(FILE* infile, struct unit_cell *c, struct contents *m,
              struct symmetry *s){
  int i,j,natoms;
  int iflag,ifhr,ifso,igr,hall;
  int ishift[3];
  char buffer[LINE_SIZE+1];
  double *ptr,abc[6],junk;
  
  buffer[0]=0;
  i=0;
  
  while(fgets(buffer,LINE_SIZE,infile)){
    if (!strncmp(buffer,"CRYSTAL",7)) break;
    if (!strncmp(buffer,"SLAB",4))
      error_exit("Reading of 2D Crystal systems not supported");
    if (!strncmp(buffer,"POLYMER",7))
      error_exit("Reading of 1D Crystal systems not supported");
    if (!strncmp(buffer,"HELIX",5))
      error_exit("Reading of 1D Crystal systems not supported");
    if (!strncmp(buffer,"MOLECULE",8))
      error_exit("Reading of 0D Crystal systems not supported");
    if (!strncmp(buffer,"EXTERNAL",8))
      error_exit("EXTERNAL keyword not supported. Read fort.34 directly");
    if (!strncmp(buffer,"DLVINPUT",8))
      error_exit("DLVINPUT keyword not supported");
  }

  fprintf(stderr,"Warning: using little-tested Crystal fort.12 reader - "
          "check output!\n");
  
  fgets(buffer,LINE_SIZE,infile);
  if (sscanf(buffer,"%d %d %d",&iflag,&ifhr,&ifso)!=3)
    error_exit("Error reading flags in CRYSTAL file");
  if ((iflag<0)||(iflag>1)||(ifhr<0)||(ifhr>1)||(ifso<0)||(ifso>24))
    error_exit("Invalid value for flags in CRYSTAL file");

  if (iflag==1)
    error_exit("H-M codes unsupported by c2x");

  fgets(buffer,LINE_SIZE,infile);
  if (sscanf(buffer,"%d",&igr)!=1)
    error_exit("Error reading spacegroup in CRYSTAL file");
  if ((igr<1)||(igr>230))
    error_exit("Invalid value for space group in CRYSTAL file");

  hall=igr2hall[igr];

  if (ifhr==1){ /*ifhr=0 is hexagonal, 1 is rhombohedral */
    if ((igr==146)||(igr==148)||(igr==155)||(igr==160)||(igr==161)||
        (igr==166)||(igr==167))
      hall+=1;
    else
      fprintf(stderr,"Unexpected value of ifhr ignored\n");
  }

  if (spgr_is_double(igr)){
    if (debug>1) {
      fprintf(stderr,"Two origins for space group, ");
      if (ifso==1) fprintf(stderr,"first used\n");
      else if (ifso==0) fprintf(stderr,"second used\n");
    }
    if (ifso==0) hall+=1;
  }
  else if (ifso==1){
    fprintf(stderr,"Problem: input file says two origins possible, "
        "but only one known to c2x\n"
        "Result may be wrong\n");
  }

  if (ifso>1){
    fgets(buffer,LINE_SIZE,infile);
    if (sscanf(buffer,"%d %d %d",ishift,ishift+1,ishift+2)!=3)
      error_exit("Error reading non-standard shift in CRYSTAL file");
  }
  
  if(hall>=489){ /* Cubic */
    fgets(buffer,LINE_SIZE,infile);
    if (sscanf(buffer,"%lf %lf",abc,&junk)!=1)
      error_exit("Error reading abc in CRYSTAL file");
    abc[1]=abc[0];
    abc[2]=abc[0];
    abc[3]=abc[4]=abc[5]=90;
  }
  else if (hall>=430){              /* Hexagonal or trigonal */
    fgets(buffer,LINE_SIZE,infile);
    if (ifhr==0){
      if (sscanf(buffer,"%lf %lf %lf",abc,abc+2,&junk)!=2)
        error_exit("Error reading abc in CRYSTAL file");
      abc[1]=abc[0];
      abc[3]=abc[4]=90;
      abc[5]=120;
    }
    else{
      if (sscanf(buffer,"%lf %lf %lf",abc,abc+3,&junk)!=2)
        error_exit("Error reading abc in CRYSTAL file");
      abc[1]=abc[0];
      abc[2]=abc[0];
      abc[4]=abc[3];
      abc[5]=abc[3];
    }
  }
  else if (hall>=349){ /* Tetragonal */
    if (sscanf(buffer,"%lf %lf %lf",abc,abc+2,&junk)!=2)
      error_exit("Error reading abc in CRYSTAL file");
    abc[1]=abc[0];
    abc[3]=abc[4]=abc[5]=90;
  }
  else if (hall>=108){ /* Orthorhombic */
    if (sscanf(buffer,"%lf %lf %lf %lf",abc,abc+1,abc+2,&junk)!=3)
      error_exit("Error reading abc in CRYSTAL file");
    abc[3]=abc[4]=abc[5]=90;
  }
  else if (hall>=3){ /* Monoclinic */  
    if (sscanf(buffer,"%lf %lf %lf %lf %lf",abc,abc+1,abc+2,abc+4,&junk)!=4)
      error_exit("Error reading abc in CRYSTAL file");
    abc[3]=abc[5]=90;
  }
  else { /* Triclinic */
    if (sscanf(buffer,"%lf %lf %lf %lf %lf %lf %lf",abc,abc+1,abc+2,
               abc+3,abc+4,abc+4,&junk)!=6)
      error_exit("Error reading abc in CRYSTAL file");
    abc[3]=abc[5]=90;
  }

  if (!(c->basis=malloc(9*sizeof(double))))
    error_exit("Malloc error in crystal_read for c->basis");

  abc2cart(abc,c);

  if (debug) fprintf(stderr,"Using Hall number of %d from IGR=%d\n",hall,igr);

  cspg_hall2sym(hall,c,s);

  if (debug) fprintf(stderr,"%d symops returned\n",s->n);
  
  natoms=0;
  fgets(buffer,LINE_SIZE,infile);
  if (sscanf(buffer,"%d",&natoms)!=1)
    error_exit("Error reading number of atoms in CRYSTAL file");

  m->atoms=malloc(natoms*sizeof(struct atom));
  if (!m->atoms) error_exit("Malloc error for atoms");
  m->n=natoms;
  init_atoms(m->atoms,m->n);
  
  for(i=0;i<natoms;i++){
    fgets(buffer,LINE_SIZE,infile);
    ptr=m->atoms[i].frac;
    if (sscanf(buffer,"%d %lf %lf %lf",&m->atoms[i].atno,ptr,ptr+1,ptr+2)!=4)
      error_exit("Error parsing atom in CRYSTAL file");
    m->atoms[i].atno=m->atoms[i].atno%100;
  }
  addabs(m->atoms,natoms,c->basis);

#ifdef SPGLIB
  if (ifso>1){
    for(i=0;i<m->n;i++)
      for(j=0;j<3;j++){
        m->atoms[i].frac[j]+=ishift[j]/24.0;
        m->atoms[i].frac[j]=fmod(m->atoms[i].frac[j],1.0);
        if (m->atoms[i].frac[j]<0) m->atoms[i].frac[j]+=1.0;
      }
    addabs(m->atoms,m->n,c->basis);
  }
  sym_expand(c,m,s);
  if (ifso>1){
    for(i=0;i<m->n;i++)
      for(j=0;j<3;j++){
        m->atoms[i].frac[j]-=ishift[j]/24.0;
        m->atoms[i].frac[j]=fmod(m->atoms[i].frac[j],1.0);
        if (m->atoms[i].frac[j]<0) m->atoms[i].frac[j]+=1.0;
      }
    addabs(m->atoms,m->n,c->basis);
  }
#endif
  
}
