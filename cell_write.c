/* Write a .cell file */


/* Copyright (c) 2007,2014 MJ Rutter 
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
#include<string.h>
#include<math.h>

#include "c2xsf.h"


static void cell_write_common(FILE* outfile, struct unit_cell *c,
			      struct contents *m, struct kpts *k,
			      struct symmetry *s);

void cell_write_species_block(FILE* outfile, struct contents *m){
  int i,j,hit;
  
  if (m->block_species){
    fprintf(outfile,"%%block SPECIES\n");
    fprintf(outfile,"%s",m->block_species);
    fprintf(outfile,"%%endblock SPECIES\n\n");   
  }
  else{ /* create basic species block */
    fprintf(outfile,"# Autogenerated species block by c2x\n");
    fprintf(outfile,"%%block SPECIES\n");
    for(i=0;i<m->n;i++){
      if (m->atoms[i].label){
        hit=0;
        for(j=0;j<i;j++){
          if ((m->atoms[j].label)&&
              (!strcmp(m->atoms[j].label,m->atoms[i].label))){
            hit=1;
            break;
          }
        }
        if (!hit) fprintf(outfile,"%s %3s 0 0 0.0\n",m->atoms[i].label,
                          atno2sym(m->atoms[i].atno));
      }
      else{  /* No label */
        hit=0;
        for(j=0;j<i;j++){
          if (m->atoms[j].atno==m->atoms[i].atno){
            hit=1;
            break;
          }
        }
        if (!hit) fprintf(outfile,"%3s %3s 0 0 0.0\n",
                          atno2sym(m->atoms[i].atno),
                          atno2sym(m->atoms[i].atno));
      }
    }
    fprintf(outfile,"%%endblock SPECIES\n\n");   
  }
}

static void cell_write_lattice_cart(FILE* outfile, struct unit_cell *c){
  int i;
  char *fmt;

  if (flags&HIPREC)
    fmt="% 19.15f % 19.15f % 19.15f\n";
  else
    fmt="% 11.7f % 11.7f % 11.7f\n";

  fprintf(outfile,"%%block LATTICE_CART\n");
  if (flags&AU){
    fprintf(outfile,"bohr\n");
    for(i=0;i<3;i++)
      fprintf(outfile,fmt,
                  c->basis[i][0]/BOHR,c->basis[i][1]/BOHR,c->basis[i][2]/BOHR);
  }
  else{
    fprintf(outfile,"ang\n");
    for(i=0;i<3;i++)
      fprintf(outfile,fmt,
                    c->basis[i][0],c->basis[i][1],c->basis[i][2]);
  }
  fprintf(outfile,"%%endblock LATTICE_CART\n\n");
}

static void cell_write_lattice_abc(FILE* outfile, double* abc){
  char *fmt;

  if (flags&HIPREC)
    fmt="% .15f % .15f % .15f\n% .15f % .15f % .15f\n";
  else
    fmt="% .7f % .7f % .7f\n% .7f % .7f % .7f\n";

  fprintf(outfile,"%%block LATTICE_ABC\n");
  if (flags&AU){
    fprintf(outfile,"bohr\n");
    fprintf(outfile,fmt,
                    abc[0]/BOHR,abc[1]/BOHR,abc[2]/BOHR,abc[3],abc[4],abc[5]);
  }
  else{
    fprintf(outfile,"ang\n");
    fprintf(outfile,fmt,
	    abc[0],abc[1],abc[2],abc[3],abc[4],abc[5]);
  }
  fprintf(outfile,"%%endblock LATTICE_ABC\n\n");


}

static void cell_write_atoms_abs(FILE* outfile, struct contents *m){
  int i;
  double scale;
  char *fmt,*fmt2;

  if (flags&HIPREC){
    fmt="%3s % .15f % .15f % .15f";
    fmt2="%3s%s % .15f % .15f % .15f";
  }
  else{
    fmt="%3s % .9f % .9f % .9f";
    fmt2="%3s%s % .9f % .9f % .9f";
  }

  if (flags&ONETEP_OUT) cell_write_species_block(outfile,m);

  fprintf(outfile,"%%block POSITIONS_ABS\n");
  if (flags&AU){
    fprintf(outfile,"bohr\n");
    scale=1.0/BOHR;
  }
  else{
    fprintf(outfile,"ang\n");
    scale=1;
  }
    
  for(i=0;i<m->n;i++){
    if ((flags&ONETEP_OUT)&&(m->atoms[i].label))
      fprintf(outfile,fmt,
              m->atoms[i].label,m->atoms[i].abs[0]*scale,
              m->atoms[i].abs[1]*scale,m->atoms[i].abs[2]*scale);
    else if ((m->atoms[i].label)&&(m->atoms[i].label[0]==':'))
      fprintf(outfile,fmt2,
              atno2sym(m->atoms[i].atno),m->atoms[i].label,
              m->atoms[i].abs[0],m->atoms[i].abs[1],m->atoms[i].abs[2]);
    else
      fprintf(outfile,fmt,
              atno2sym(m->atoms[i].atno),m->atoms[i].abs[0]*scale,
              m->atoms[i].abs[1]*scale,m->atoms[i].abs[2]*scale);
    if (m->atoms[i].spin)
      fprintf(outfile," SPIN=%f\n",m->atoms[i].spin);
    else
      fprintf(outfile,"\n");
    
  }
  fprintf(outfile,"%%endblock POSITIONS_ABS\n");
}

static void cell_write_atoms(FILE *outfile, struct contents *m){
  int i;
  char *fmt,*fmt2;

  if (flags&HIPREC){
    fmt="%3s % .15f % .15f % .15f";
    fmt2="%3s%s % .15f % .15f % .15f";
  }
  else{
    fmt="%3s % .9f % .9f % .9f";
    fmt2="%3s%s % .9f % .9f % .9f";
  }
  
  if (flags&ONETEP_OUT) cell_write_species_block(outfile,m);

  fprintf(outfile,"%%block POSITIONS_FRAC\n");
  for(i=0;i<m->n;i++){
    if ((flags&ONETEP_OUT)&&(m->atoms[i].label))
      fprintf(outfile,fmt,
              m->atoms[i].label,m->atoms[i].frac[0],
              m->atoms[i].frac[1],m->atoms[i].frac[2]);
    else if ((m->atoms[i].label)&&(m->atoms[i].label[0]==':'))
      fprintf(outfile,fmt2,
              atno2sym(m->atoms[i].atno),m->atoms[i].label,
              m->atoms[i].frac[0],m->atoms[i].frac[1],m->atoms[i].frac[2]);
    else
      fprintf(outfile,fmt,
              atno2sym(m->atoms[i].atno),m->atoms[i].frac[0],
              m->atoms[i].frac[1],m->atoms[i].frac[2]);
    if (m->atoms[i].spin)
      fprintf(outfile," SPIN=%f\n",m->atoms[i].spin);
    else
      fprintf(outfile,"\n");
  }
  fprintf(outfile,"%%endblock POSITIONS_FRAC\n");
}


void cell_write_abc(FILE* outfile, struct unit_cell *c, struct contents *m,
                    struct kpts *k, struct symmetry *s){
  double abc[6];

  cart2abc(c,m,abc,NULL,1);

  if (m->title)
    fprintf(outfile,"#TITL %s\n",m->title);

  cell_write_lattice_abc(outfile,abc);
  cell_write_atoms(outfile,m);
  cell_write_common(outfile,c,m,k,s);
}

void cell_write(FILE* outfile, struct unit_cell *c, struct contents *m,
                struct kpts *k, struct symmetry *s){
  if (m->title)
    fprintf(outfile,"#TITL %s\n",m->title);

  cell_write_lattice_cart(outfile,c);
  cell_write_atoms(outfile,m);
  cell_write_common(outfile,c,m,k,s);
}

void cell_write_abs(FILE* outfile, struct unit_cell *c, struct contents *m,
                    struct kpts *k, struct symmetry *s){

  if (m->title)
    fprintf(outfile,"#TITL %s\n",m->title);

  cell_write_lattice_cart(outfile,c);
  cell_write_atoms_abs(outfile,m);
  cell_write_common(outfile,c,m,k,s);
}


void cell_write_abc_abs(FILE* outfile, struct unit_cell *c, struct contents *m,
                        struct kpts *k, struct symmetry *s){
  double abc[6];

  cart2abc(c,m,abc,NULL,1);

  if (m->title)
    fprintf(outfile,"#TITL %s\n",m->title);

  cell_write_lattice_abc(outfile,abc);
  cell_write_atoms_abs(outfile,m);
  cell_write_common(outfile,c,m,k,s);
}


void cell_write_common(FILE* outfile, struct unit_cell *c, struct contents *m,
                     struct kpts *k, struct symmetry *s){
  int i,j;
  double v[3];
  char *fmt1,*fmt2,*fmt3;

  fmt1="% 19.15f % 19.15f % 19.15f\n";
  if (flags&HIPREC){
    fmt2="% 19.15f % 19.15f % 19.15f     %19.15f\n";
    fmt3="KPOINT_MP_OFFSET %.9f %.9f %.9f\n";
  }
  else{
    fmt2="% 13.9f % 13.9f % 13.9f     %12.9f\n";
    fmt3="KPOINT_MP_OFFSET %.9f %.9f %.9f\n";
  }

  if (m->species_misc)
    fprintf(outfile,"\n%s",m->species_misc);
  
  if (((s->ops)&&(s->n>1))||((s->gen)&&(*s->gen))){
    if (s->tol){
      if (flags&AU)
        fprintf(outfile,"\nSYMMETRY_TOL %.8f bohr\n",(*s->tol)/BOHR);
      else
        fprintf(outfile,"\nSYMMETRY_TOL %.8f ang\n",*s->tol);
    }
  }

  if ((s->ops)&&(s->n>1)){
    fprintf(outfile,"\n%%block SYMMETRY_OPS");
    for(i=0;i<s->n;i++){
      fprintf(outfile,"\n");
      if (debug){
        fprintf(outfile,"! ");
        equiv_sym(s->ops+i,c,outfile);
        fprintf(outfile,"! ");
        ident_sym(s->ops+i,c,outfile);
      }
      for(j=0;j<3;j++)
        fprintf(outfile,fmt1,s->ops[i].mat[j][0],
                s->ops[i].mat[j][1],s->ops[i].mat[j][2]);
      /* Though the matrix is in cartesian co-ords, the associated
       * translation is expected in fractional co-ords... */
      if (s->ops[i].tr){
        for(j=0;j<3;j++)
          v[j]=s->ops[i].tr[0]*c->recip[j][0]+s->ops[i].tr[1]*c->recip[j][1]+
            s->ops[i].tr[2]*c->recip[j][2];
        fprintf(outfile,fmt1,v[0],v[1],v[2]);
      }
      else
        fprintf(outfile,fmt1,0.0,0.0,0.0);
    }
    fprintf(outfile,"%%endblock SYMMETRY_OPS\n");
  }
  else if ((s->gen)&&(*s->gen))
    fprintf(outfile,"\nSYMMETRY_GENERATE\n");

  if ((k->mp)&&(k->mp->grid[0]>0)){
    fprintf(outfile,"\nKPOINT_MP_GRID %d %d %d\n",k->mp->grid[0],
            k->mp->grid[1],k->mp->grid[2]);
    /* A shift of +/- 1/12 with a grid size of 6 is much better written
     * as a shift of 1/4, which is equivalent as 1/12+1/6=1/4
     * This generalises.
     */
    for(i=0;i<3;i++){
      v[i]=k->mp->disp[i];
      if ((k->mp->grid[i])&&((k->mp->grid[i]&1)==0)&&
          (aeq(fabs(v[i]),0.5/k->mp->grid[i]))){
        v[i]=0.5;
        j=1;
        while((j&k->mp->grid[i])==0){
          j=j<<1;
          v[i]*=0.5;
        }
      }
    }

    fprintf(outfile,fmt3,v[0],v[1],v[2]);
  } else if(k->n) {
    fprintf(outfile,"\n%%block KPOINTS_LIST\n");
    for(i=0;i<k->n;i++)
      fprintf(outfile,fmt2,k->kpts[i].frac[0],
              k->kpts[i].frac[1],k->kpts[i].frac[2],k->kpts[i].wt);
    fprintf(outfile,"%%endblock KPOINTS_LIST\n");
  }

}
