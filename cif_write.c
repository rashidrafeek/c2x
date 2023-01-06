/* Write a CIF file */

/* Copyright (c) 2014-2021 MJ Rutter 
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

/* Write a file in format similar to that which Castep would use 
 * for CIF output. Though this format may be compatible with CIF / mmCIF,
 * it is very far from making use of all the features of that format. */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

#include "c2xsf.h"

/* From c2x2spg.c */

struct sym_op *sym_frac2abs(int spg_rot[][3][3],double spg_tr[][3],
			    struct unit_cell *c,int nsym);

extern int periodic_max_el;

struct contents *reduce_atoms(struct contents *m, struct symmetry *s,
			      struct unit_cell *c);

void cif_write(FILE* outfile, struct unit_cell *c, struct contents *m, 
	       struct symmetry *s, int mm){
  int i,j,prec,site_charge,sym,*n_in_el,label,pdbx;
  unsigned char sep;
  double abc[6];

  pdbx=0;
  
  if (dict_get(m->dict,"CIF_is_PDBx")){
    mm=1;
    pdbx=1;
  }

  sep='_';
  if (mm) sep='.';

  prec=7;
  if (flags&HIPREC) prec=15;

  sym=0;
  if (dict_get(m->dict,"CIF_symmetrise")) sym=1;

  label=0;
  if (dict_get(m->dict,"CIF_site_label")) label=1;
  
  n_in_el=calloc(periodic_max_el+1,sizeof(int));
  if (!n_in_el) error_exit("Calloc error in cif_write");

  
  fprintf(outfile,"# written by c2x\n");

  fprintf(outfile,"data_c2x\n");

  if (m->title) fprintf(outfile,"_struct%ctitle %s\n",sep,m->title);

  cart2abc_sym(c,m,abc,NULL,s);
  //  make_rhs(c,m,NULL,NULL);
  //  cart2abc(c,NULL,abc,NULL,1);

  fprintf(outfile,"\n");
  fprintf(outfile,"_cell%clength_a      %.*f\n",sep,prec,abc[0]);
  fprintf(outfile,"_cell%clength_b      %.*f\n",sep,prec,abc[1]);
  fprintf(outfile,"_cell%clength_c      %.*f\n",sep,prec,abc[2]);
  fprintf(outfile,"_cell%cangle_alpha   %.*f\n",sep,prec,abc[3]);
  fprintf(outfile,"_cell%cangle_beta    %.*f\n",sep,prec,abc[4]);
  fprintf(outfile,"_cell%cangle_gamma   %.*f\n",sep,prec,abc[5]);
  fprintf(outfile,"\n");

  site_charge=0;
  if (dict_get(m->dict,"site_charge")) site_charge=1;

  if (!pdbx){
    fprintf(outfile,"\n");
    fprintf(outfile,"loop_\n"
	    "_atom_site%ctype_symbol\n",sep);
    if (label) fprintf(outfile,"_atom_site%clabel\n",sep);
    fprintf(outfile,"_atom_site%cfract_x\n"
	    "_atom_site%cfract_y\n"
	    "_atom_site%cfract_z\n"
	    "_atom_site%cU_iso_or_equiv\n"
	    "_atom_site%coccupancy\n",sep,sep,sep,sep,sep);
    if (site_charge)
      fprintf(outfile,"_atom_site%ccharge\n",sep);

    if (prec<10) prec=10;

    if (sym) m=reduce_atoms(m,s,c);
  
    for(i=0;i<m->n;i++){
      fprintf(outfile,"%s ",atno2sym(m->atoms[i].atno));
      if (label){
	fprintf(outfile,"%s%-3d ",atno2sym(m->atoms[i].atno),
		++n_in_el[m->atoms[i].atno]);
      }
      for(j=0;j<3;j++)
	fprintf(outfile,"%.*f ",prec,m->atoms[i].frac[j]);
      fprintf(outfile,"0.01 1.00");
      if (site_charge) fprintf(outfile," %.4f",m->atoms[i].site_chg);
      fprintf(outfile,"\n");
    }
    if (sym){
      free(m->atoms);
      free(m);
    }
  }
  else{
    fprintf(outfile,"\n");
    fprintf(outfile,"loop_\n"
	    "_atom_site.group_PDB\n"
	    "_atom_site.id\n"
	    "_atom_site.type_symbol\n"
	    "_atom_site.label_atom_id\n"
	    "_atom_site.label_alt_id\n"
	    "_atom_site.label_comp_id\n"
	    "_atom_site.label_asym_id\n"
	    "_atom_site.label_entity_id\n"
	    "_atom_site.label_seq_id\n"
	    "_atom_site.pdbx_PDB_ins_code\n"
	    "_atom_site.Cartn_x\n"
	    "_atom_site.Cartn_y\n"
	    "_atom_site.Cartn_z\n"
	    "_atom_site.occupancy\n"
	    "_atom_site.B_iso_or_equiv\n"
	    "_atom_site.Cartn_x_esd\n"
	    "_atom_site.Cartn_y_esd\n"
	    "_atom_site.Cartn_z_esd\n"
	    "_atom_site.occupancy_esd\n"
	    "_atom_site.B_iso_or_equiv_esd\n"
	    "_atom_site.pdbx_formal_charge\n"
	    "_atom_site.auth_seq_id\n"
	    "_atom_site.auth_comp_id\n"
	    "_atom_site.auth_asym_id\n"
	    "_atom_site.auth_atom_id\n"
	    "_atom_site.pdbx_PDB_model_num\n");
    for(i=0;i<m->n;i++){
      fprintf(outfile,"ATOM %4d ",i);
      fprintf(outfile,"%s ",atno2sym(m->atoms[i].atno));
      fprintf(outfile,"%s%-3d ",atno2sym(m->atoms[i].atno),
	      ++n_in_el[m->atoms[i].atno]);
      fprintf(outfile,". . . . . . ");
      for(j=0;j<3;j++)
	fprintf(outfile,"%.*f ",prec,m->atoms[i].abs[j]);
      fprintf(outfile,"1.0 . . . . . . . . . . ");
      fprintf(outfile,"%s%-3d 1\n",atno2sym(m->atoms[i].atno),
	      n_in_el[m->atoms[i].atno]);
    }
  }
  if (s->n){ /* Must print identity as first operation */
    fprintf(outfile,"\nloop_\n");
    if (mm==0){
      fprintf(outfile,"_symmetry_equiv_pos_as_xyz\n");
      fprintf(outfile,"x,y,z\n");
      for(i=0;i<s->n;i++)
	if ((s->ops[i].tr)||(!is_identity(s->ops[i].mat)))
	  equiv_sym(s->ops+i,c,outfile);
    }
    else{
      fprintf(outfile,"_space_group_symop.id\n");
      fprintf(outfile,"_space_group_symop.operation_xyz\n");
      fprintf(outfile,"1  x,y,z\n");
      j=2;
      for(i=0;i<s->n;i++){
	if ((s->ops[i].tr)||(!is_identity(s->ops[i].mat))){
	  fprintf(outfile,"%-2d ",j);
	  equiv_sym(s->ops+i,c,outfile);
	  j++;
	}
      }
    }
  }

  free(n_in_el);
  
}

struct contents *reduce_atoms(struct contents *m, struct symmetry *s,
			      struct unit_cell *c){
  int i,j,k,hit,del;
  struct contents *m2;
  struct atom new_atom;

  m2=malloc(sizeof(struct contents));
  if (!m2) error_exit("Malloc error for struct contents");
  *m2=*m;

  m2->atoms=malloc(m->n*sizeof(struct atom));
  if (!m2->atoms) error_exit("Malloc error for struct atoms");
  memcpy(m2->atoms,m->atoms,m->n*sizeof(struct atom));

  for(i=0;i<m2->n;i++){
    for(j=0;j<s->n;j++){
      sym_atom(m2->atoms+i,&new_atom,s->ops+j,c->recip);
      hit=atom_in_list(&new_atom,m2->atoms,m2->n,c->basis);
      if (hit==-1){
	hit=atom_in_list(&new_atom,m->atoms,m->n,c->basis);
	if (hit!=-1) continue;
	fprintf(stderr,"Symmetry error\n");
	fprintf(stderr,"Atom atno=%d (%lf,%lf,%lf)\n",m2->atoms[i].atno,
		m2->atoms[i].frac[0],m2->atoms[i].frac[1],m2->atoms[i].frac[2]);
	fprintf(stderr,"Sym op ");
	ident_sym(s->ops+j,c,m,stderr);
	fprintf(stderr,"New atom at (%lf,%lf,%lf)\n",new_atom.frac[0],
		new_atom.frac[1],new_atom.frac[2]);
	exit(1);
      }
      if (hit==i) continue;
      if (hit>i){
	del=hit;
      }
      else{
	del=i;
      }
      for(k=del;k<m2->n-1;k++){
	m2->atoms[k]=m2->atoms[k+1];
      }
      m2->n--;
      if (del==i) {i--; break;}
    }
  }
      
  m2->atoms=realloc(m2->atoms,m2->n*sizeof(struct atom));

  if (debug)
    fprintf(stderr,"%d atoms before symmetrisation, %d after\n",
	    m->n,m2->n);

  return m2;

}
