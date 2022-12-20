/* Read a PWscf input file */

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

/*
 * Like abinit, multiple keywords can appear on a single line, so
 * abinit's readline is reused
 *
 */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<math.h>

#include "c2xsf.h"

#define LINE_SIZE 132

#define QECONTROL 1
#define QESYSTEM 2
#define QEELECTRONS 3
#define QEIONS 4
#define QECELL 5

static int qe_readline(char **p, FILE* infile);

static int qe_int_read(char **ptr,int *n){
  int i;
  
  if ((sscanf(*ptr,"%d%n",n,&i)!=1)&&(sscanf(*ptr," = %d%n",n,&i)!=1))
    return 0;

  (*ptr)+=i;
  while((*ptr)&&(isspace(**ptr))) (*ptr)++;
  if (**ptr==',') (*ptr)++;

  return 1;
}

static int qe_float_read(char **ptr,double *x){
  int i;
  
  if ((sscanf(*ptr,"%lf%n",x,&i)!=1)&&(sscanf(*ptr," = %lf%n",x,&i)!=1))
    return 0;

  *ptr+=i;
  while((*ptr)&&(isspace(**ptr))) (*ptr)++;
  if (**ptr==',') (*ptr)++;

  return 1;
}

void qe_read(FILE* infile, struct unit_cell *c, struct contents *m,
             struct kpts *k, struct symmetry *s, struct es *e){
  int i,j,namelist;
  int ibrav,ntyp;
  double celldm[7];
  double *abc,x,y,z,c2,alpha,beta,gamma,scale;
  double tx,ty,tz,u,v;
  char sym[4],*ptr,*p2;
  int *aspec;
  int atoms_are_abs,off[3];

  ntyp=0;
  ibrav=-1;
  namelist=0;
  abc=NULL;
  aspec=NULL;
  atoms_are_abs=0;
  for(i=0;i<7;i++)celldm[i]=0;
  ptr="";
  
  if (!(c->basis=malloc(9*sizeof(double))))
    error_exit("Malloc error for basis");
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[i][j]=0;

  while(1){
    /* remove leading spaces */
    while(*ptr&&(isspace(*ptr))) ptr++;
    if ((*ptr=='\n')||(*ptr==0)){
      if (!qe_readline(&ptr,infile)) break;
    }

#if 0
    fprintf(stderr,"Parsing %s\n",ptr);
    sleep(1);
#endif
    
    if ((namelist)&&(!tokenmatch(&ptr,"/"))){
      namelist=0;
      continue;
    }
    if ((!namelist)&&(!tokenmatch(&ptr,"&control"))){
      namelist=QECONTROL;
      continue;
    }
    if ((!namelist)&&(!tokenmatch(&ptr,"&system"))){
      namelist=QESYSTEM;
      continue;
    }
    if ((!namelist)&&(!tokenmatch(&ptr,"&electrons"))){
      namelist=QEELECTRONS;
      continue;
    }
    if ((!namelist)&&(!tokenmatch(&ptr,"&ions"))){
      namelist=QEIONS;
      continue;
    }
    if ((!namelist)&&(!tokenmatch(&ptr,"&cell"))){
      namelist=QECELL;
      continue;
    }

    if (namelist==QECONTROL){
      if (!tokenmatch(&ptr,"title")){
        while ((*ptr)&&(isspace(*ptr))) ptr++;
        if (*ptr=='=') ptr++;
        while ((*ptr)&&(isspace(*ptr))) ptr++;
        if (*ptr!='\'')
          error_exit("Syntax error for title");
        ptr++;
        p2=ptr;
        while ((*p2)&&(*p2!='\'')) p2++;
        *p2=0;
        i=strlen(ptr)+1;
        m->title=malloc(i);
        strncpy(m->title,ptr,i);
        ptr=p2;
      }
      else{
        while ((*ptr)&&(!isspace(*ptr))) ptr++;  /* consume a token */
      }
    }      
    else if (namelist==QESYSTEM){
      if (!tokenmatch(&ptr,"nat")){
        if (!qe_int_read(&ptr,&(m->n)))
          error_exit("Error parsing nat");
        m->atoms=malloc(m->n*sizeof(struct atom));
        if (!m->atoms)
          error_exit("Error in malloc for atoms");
        init_atoms(m->atoms,m->n);
      }
      else if (!tokenmatch(&ptr,"ntyp")){
        if (!qe_int_read(&ptr,&ntyp))
          error_exit("Error parsing ntyp");
      }
      else if (!tokenmatch(&ptr,"ibrav")){
        if (!qe_int_read(&ptr,&ibrav))
          error_exit("Error parsing ibrav");
      }
      else if (!tokenmatch(&ptr,"ecutwfc")){
        if (!qe_float_read(&ptr,&(e->cut_off)))
          error_exit("Error parsing ecutwfc");
        e->cut_off*=H_eV;
      }
      else if (!strncasecmp(ptr,"celldm(",7)){
        ptr+=7;
        if (sscanf(ptr,"%d",&i)!=1)
          error_exit("Error parsing celldm");
        ptr+=2;
        if (!qe_float_read(&ptr,celldm+i))
          error_exit("Error parsing celldm");
      }
      else if (!tokenmatch(&ptr,"A")){
        if (!abc) abc=malloc(6*sizeof(double));
        if (!qe_float_read(&ptr,abc))
          error_exit("Error parsing A");
      }
      else if (!tokenmatch(&ptr,"B")){
        if (!abc) abc=malloc(6*sizeof(double));
        if (!qe_float_read(&ptr,abc+1))
          error_exit("Error parsing B");
      }
      else if (!tokenmatch(&ptr,"C")){
        if (!abc) abc=malloc(6*sizeof(double));
        if (!qe_float_read(&ptr,abc+2))
          error_exit("Error parsing C");
      }
      else if (!tokenmatch(&ptr,"cosAB")){
        if (!abc) abc=malloc(6*sizeof(double));
        if (!qe_float_read(&ptr,&x))
          error_exit("Error parsing cosAB");
        abc[5]=acos(x)*180/M_PI;
      }
      else if (!tokenmatch(&ptr,"cosAC")){
        if (!abc) abc=malloc(6*sizeof(double));
        if (!qe_float_read(&ptr,&x))
          error_exit("Error parsing cosAC");
        abc[4]=acos(x)*180/M_PI;
      }
      else if (!tokenmatch(&ptr,"cosBC")){
        if (!abc) abc=malloc(6*sizeof(double));
        if (!qe_float_read(&ptr,&x))
          error_exit("Error parsing cosBC");
        abc[3]=acos(x)*180/M_PI;
      }
      else{
        while ((*ptr)&&(!isspace(*ptr))) ptr++;
      }
    }
    else if (namelist==0){
      if (!tokenmatch(&ptr,"atomic_species")){
        if (ntyp==0) error_exit("atomic_species found before ntyp");
        aspec=malloc(ntyp*sizeof(int));
        for(i=0;i<ntyp;i++){
          qe_readline(&ptr,infile);
          sym[3]=0;
          for(j=0;j<3;j++){
            sym[j]=0;
            if (!isalpha(*ptr)) break;
            sym[j]=*(ptr++);
          }
          for(j=3;j>0;j++){
            sym[j]=0;
            aspec[i]=atsym2no(sym);
            if (aspec[i]) break;
          }
        }
      }
      else if (!tokenmatch(&ptr,"atomic_positions")){
        while (*ptr==' ') ptr++;
        if ((*ptr=='(')||(*ptr=='{')) ptr++;
        scale=0;
        if (!strncasecmp(ptr,"angstrom",8)){
          scale=1;
          atoms_are_abs=1;
        }
        else if (!strncasecmp(ptr,"bohr",4)){
          scale=BOHR;
          atoms_are_abs=1;
        }
        else if ((!strncasecmp(ptr,"alat",4))||
                 (*ptr=='!')||(*ptr=='#')||(*ptr==0)||(*ptr=='\n')){
          scale=celldm[1]*BOHR;
          if (abc) scale=abc[0];
          atoms_are_abs=1;
        }
        else if (!strncasecmp(ptr,"crystal",7)){
          atoms_are_abs=0;
        }
        else fprintf(stderr,
                     "Warning: confused by %s on atomic_positions line\n",
                     ptr);
        if (m->n==0) error_exit("atomic_positions found before nat");
        for(i=0;i<m->n;i++){
          qe_readline(&ptr,infile);
          /* process label */
          sym[3]=0;
          for(j=0;j<3;j++){
            sym[j]=0;
            if (!isalpha(*ptr)) break;
            sym[j]=*(ptr++);
          }
          for(j=3;j>0;j++){
            sym[j]=0;
            m->atoms[i].atno=atsym2no(sym);
            if (m->atoms[i].atno) break;
          }
          /* Eat any trailing non-alpha on symbol */
          if ((*ptr)&&(!isspace(*ptr))) ptr++;
          if (atoms_are_abs){
            if (multi_scan(ptr,m->atoms[i].abs,3,NULL)!=3)
              error_exit("Error parsing atomic positions");
            for(j=0;j<3;j++)
              m->atoms[i].abs[j]*=scale;
          }
          else{
            if (multi_scan(ptr,m->atoms[i].frac,3,NULL)!=3)
              error_exit("Error parsing atomic positions");
          }
        }
        /* Move to next line */
        *ptr=0;
      }
      else if (!tokenmatch(&ptr,"k_points")){
        while (*ptr==' ') ptr++;
        if ((*ptr=='(')||(*ptr=='{')) ptr++;
        if (!strncasecmp(ptr,"gamma",5)){
          k->n=1;
          k->kpts=malloc(sizeof(struct atom));
          k->kpts[0].wt=1;
          for(i=0;i<3;i++)
            k->kpts[0].frac[i]=0;
          /* Move to next line */
          *ptr=0;
        }
        else if (!strncasecmp(ptr,"automatic",9)){
          k->mp=malloc(sizeof(struct mp_grid));
          qe_readline(&ptr,infile);
          if (sscanf(ptr,"%d %d %d %d %d %d",k->mp->grid,
                     k->mp->grid+1,k->mp->grid+2,
                     off,off+1,off+2)!=6)
            error_exit("Error parsing k-point MP grid");
          for(i=0;i<3;i++)
            /* QE's convention is that all grids include origin,
             * ours that only odd grids include the origin.
             * So shift if off[i]==0 and even grid, or if
             * off[i]==1 and odd grid.
             */
            if (((off[i]==0)&&((k->mp->grid[i]&1)==0))||
                ((off[i]==1)&&((k->mp->grid[i]&1)==1)))
              k->mp->disp[i]=0.5/k->mp->grid[i];
            else
              k->mp->disp[i]=0;
          /* Move to next line */
          *ptr=0;
        }
        else if ((!strncasecmp(ptr,"crystal",7))&&(ptr[7]!='_')){
          qe_readline(&ptr,infile);
          if (sscanf(ptr,"%d",&k->n)!=1)
            error_exit("Error parsing kpoints");
          k->kpts=malloc(k->n*sizeof(struct atom));
          if (!k->kpts) error_exit("Malloc error for kpoints");
          for(i=0;i<k->n;i++){
            qe_readline(&ptr,infile);
            if (sscanf(ptr,"%lf %lf %lf %lf",k->kpts[i].frac,
                       k->kpts[i].frac+1,k->kpts[i].frac+2,&(k->kpts[i].wt))!=4)
              error_exit("Error parsing kpoint");
          }
          /* Move to next line */
          *ptr=0;
        }
        else if ((!strncasecmp(ptr,"tpiba",5))||
                 (!strncasecmp(ptr,"crystal_b",9))||
                 (!strncasecmp(ptr,"crystal_c",9))){
          fprintf(stderr,"Unable to parse kpoints card of type %s\n",ptr);
          qe_readline(&ptr,infile);
          if (sscanf(ptr,"%d",&i)!=1)
            error_exit("Error parsing number of kpoints");
          for(j=0;j<i;j++) qe_readline(&ptr,infile);
          *ptr=0;
        }
        else{
          fprintf(stderr,"Unable to parse kpoints card of this type\n");
          *ptr=0;
        }
      }
      else if (!tokenmatch(&ptr,"cell_parameters")){
        while (*ptr==' ') ptr++;
        if ((*ptr=='(')||(*ptr=='{')) ptr++;
        scale=0;
        if (!strncasecmp(ptr,"angstrom",8)){
          scale=1;
        }
        else if (!strncasecmp(ptr,"bohr",4)){
          scale=BOHR;
        }
        else if ((!strcasecmp(ptr,"alat"))||
                 (*ptr=='!')||(*ptr=='#')||(*ptr==0)){
          scale=celldm[1]*BOHR;
          if (abc) scale=abc[0];
        }
        else fprintf(stderr,
                     "Warning: confused by %s on cell_parameters line\n",
                     ptr);
        for(i=0;i<3;i++){
          if (!qe_readline(&ptr,infile))
            error_exit("Read error for cell parameters");
          if (sscanf(ptr,"%lf %lf %lf",c->basis[i],
                     c->basis[i]+1,c->basis[i]+2)!=3)
            error_exit("Parse error for cell_parameters");
        }
        for(i=0;i<3;i++)
          for(j=0;j<3;j++)
            c->basis[i][j]*=scale;
        *ptr=0;
      }
      else{
        while ((*ptr)&&(!isspace(*ptr))) ptr++;
      }      
    }
    else{
      while ((*ptr)&&(!isspace(*ptr))) ptr++;
    }
  }

  /* Sort out lattice */

  if (ibrav==1){ /* cubic */
    x=celldm[1]*BOHR;
    if (abc) x=abc[0];
    c->basis[0][0]=c->basis[1][1]=c->basis[2][2]=x;
  }
  else if (ibrav==2){ /* fcc */
    x=0.5*celldm[1]*BOHR;
    if (abc) x=0.5*abc[0];
    c->basis[0][2]=c->basis[1][1]=c->basis[1][2]=c->basis[2][1]=x;
    c->basis[0][0]=c->basis[2][0]=-x;
  }
  else if (ibrav==3){ /* bcc */
    x=0.5*celldm[1]*BOHR;
    if (abc) x=0.5*abc[0];
    c->basis[0][0]=c->basis[0][1]=c->basis[0][2]=
      c->basis[1][1]=c->basis[1][2]=c->basis[2][2]=x;
    c->basis[1][0]=c->basis[2][0]=c->basis[2][1]=-x;
  }
  else if (ibrav==4){ /* hexagonal */
    x=celldm[1]*BOHR;
    if (abc) x=abc[0];
    y=celldm[3]*x;
    if (abc) y=abc[2];
    c->basis[0][0]=x;
    c->basis[1][0]=-0.5*x;
    c->basis[1][1]=0.5*sqrt(3)*x;
    c->basis[2][2]=y;
  }
  else if (ibrav==5){ /* Trigonal R, 3-axis c */
    x=celldm[1]*BOHR;
    if (abc) x=abc[0];
    c2=celldm[4];
    if (abc) c2=cos(M_PI*abc[3]/180);
    tx=sqrt((1-c2)/2);
    ty=sqrt((1-c2)/6);
    tz=sqrt((1+2*c2)/3);
    c->basis[0][0]=x*tx;
    c->basis[0][1]=-x*ty;
    c->basis[0][2]=x*tz;
    c->basis[1][1]=2*x*ty;
    c->basis[1][2]=x*tz;
    c->basis[2][0]=-x*tx;
    c->basis[2][1]=-x*ty;
    c->basis[2][2]=x*tz;
  }
  else if (ibrav==-5){ /* Trigonal R, 3-axis 111 */
    x=celldm[1]*BOHR;
    if (abc) x=abc[0];
    c2=celldm[4];
    if (abc) c2=cos(M_PI*abc[3]/180);
    tx=sqrt((1-c2)/2);
    ty=sqrt((1-c2)/6);
    tz=sqrt((1+2*c2)/3);
    x=x/sqrt(3);
    u=tz-2*sqrt(2)*ty;
    v=tz+sqrt(2)*ty;
    
    c->basis[0][0]=x*u;
    c->basis[0][1]=x*v;
    c->basis[0][2]=x*v;
    c->basis[1][0]=x*v;
    c->basis[1][1]=x*u;
    c->basis[1][2]=x*v;
    c->basis[2][0]=x*v;
    c->basis[2][1]=x*v;
    c->basis[2][2]=x*u;
  }
  else if (ibrav==6){ /* Tetragonal P */   
    x=celldm[1]*BOHR;
    if (abc) x=abc[0];
    c->basis[0][0]=x;
    c->basis[1][1]=x;
    y=x*celldm[3];
    if (abc) y=abc[2];
    c->basis[2][2]=y;
  }
  else if (ibrav==7){ /* Tetragonal P */   
    x=celldm[1]*BOHR;
    if (abc) x=abc[0];
    y=x*celldm[3];
    if (abc) y=abc[2];
    c->basis[0][0]=0.5*x;
    c->basis[0][1]=-0.5*x;
    c->basis[0][2]=0.5*y;
    c->basis[1][0]=0.5*x;
    c->basis[1][1]=0.5*x;
    c->basis[2][2]=0.5*y;
    c->basis[2][2]=y;
  }
  else if (ibrav==8){ /* Orthorhombic P */   
    x=celldm[1]*BOHR;
    if (abc) x=abc[0];
    c->basis[0][0]=x;
    y=x*celldm[2];
    if (abc) y=abc[1];
    c->basis[1][1]=y;
    y=x*celldm[3];
    if (abc) y=abc[2];
    c->basis[2][2]=y;
  }
  else if (ibrav==9){ /* Orthorhombic base centered */   
    x=celldm[1]*BOHR;
    if (abc) x=abc[0];
    y=x*celldm[2];
    if (abc) y=abc[1];
    c->basis[0][0]=x/2;
    c->basis[0][1]=y/2;
    c->basis[1][0]=-x/2;
    c->basis[1][1]=y/2;
    y=x*celldm[3];
    if (abc) y=abc[2];
    c->basis[2][2]=y;
  }
  else if (ibrav==-9){ /* Orthorhombic base centered */   
    x=celldm[1]*BOHR;
    if (abc) x=abc[0];
    y=x*celldm[2];
    if (abc) y=abc[1];
    c->basis[0][0]=x/2;
    c->basis[0][1]=-y/2;
    c->basis[1][0]=x/2;
    c->basis[1][1]=y/2;
    y=x*celldm[3];
    if (abc) y=abc[2];
    c->basis[2][2]=y;
  }
  else if (ibrav==10){ /* Orthorhombic face centered */   
    x=celldm[1]*BOHR;
    if (abc) x=abc[0];
    y=x*celldm[2];
    if (abc) y=abc[1];
    z=x*celldm[3];
    if (abc) z=abc[2];
    c->basis[0][0]=x/2;
    c->basis[0][2]=z/2;
    c->basis[1][0]=x/2;
    c->basis[1][1]=y/2;
    c->basis[2][1]=y/2;
    c->basis[2][2]=z/2;
  }
  else if (ibrav==11){ /* Orthorhombic body centered */   
    x=celldm[1]*BOHR;
    if (abc) x=abc[0];
    y=x*celldm[2];
    if (abc) y=abc[1];
    z=x*celldm[3];
    if (abc) z=abc[2];
    c->basis[0][0]=x/2;
    c->basis[0][1]=y/2;
    c->basis[0][2]=z/2;
    c->basis[1][0]=-x/2;
    c->basis[1][1]=y/2;
    c->basis[1][2]=z/2;
    c->basis[2][0]=-x/2;
    c->basis[2][1]=-y/2;
    c->basis[2][2]=z/2;
  }
  else if (ibrav==12){ /* Monoclinic P, unique c */   
    x=celldm[1]*BOHR;
    if (abc) x=abc[0];
    y=x*celldm[2];
    if (abc) y=abc[1];
    z=x*celldm[3];
    if (abc) z=abc[2];
    gamma=acos(celldm[4]);
    if (abc) gamma=M_PI*abc[5]/180;
    c->basis[0][0]=x;
    c->basis[1][0]=y*cos(gamma);
    c->basis[1][1]=y*sin(gamma);
    c->basis[2][2]=z;
  }
  else if (ibrav==-12){ /* Monoclinic P, unique b */   
    x=celldm[1]*BOHR;
    if (abc) x=abc[0];
    y=x*celldm[2];
    if (abc) y=abc[1];
    z=x*celldm[3];
    if (abc) z=abc[2];
    beta=acos(celldm[5]);
    if (abc) beta=M_PI*abc[4]/180;
    c->basis[0][0]=x;
    c->basis[1][0]=y;
    c->basis[2][0]=z*cos(beta);
    c->basis[2][2]=z*sin(beta);
  }
  else if (ibrav==13){ /* Monoclinic base-centred */   
    x=celldm[1]*BOHR;
    if (abc) x=abc[0];
    y=x*celldm[2];
    if (abc) y=abc[1];
    z=x*celldm[3];
    if (abc) z=abc[2];
    gamma=acos(celldm[4]);
    if (abc) gamma=M_PI*abc[5]/180;
    c->basis[0][0]=x/2;
    c->basis[0][2]=-z/2;
    c->basis[1][0]=y*cos(gamma);
    c->basis[1][1]=y*sin(gamma);
    c->basis[2][0]=x/2;
    c->basis[2][2]=z/2;
  }
  else if ((ibrav==14)||(abc)){  /* Triclinic */   
    x=celldm[1]*BOHR;
    if (abc) x=abc[0];
    y=x*celldm[2];
    if (abc) y=abc[1];
    z=x*celldm[3];
    if (abc) z=abc[2];
    alpha=acos(celldm[4]);
    beta=acos(celldm[5]);
    gamma=acos(celldm[6]);
    if (abc){
      alpha=M_PI*abc[3]/180;
      beta=M_PI*abc[4]/180;
      gamma=M_PI*abc[5]/180;
    }
    c->basis[0][0]=x;
    c->basis[1][0]=y*cos(gamma);
    c->basis[1][1]=y*sin(gamma);
    c->basis[2][0]=z*cos(beta);
    c->basis[2][1]=z*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);
    c->basis[2][2]=z*sqrt(1+2*cos(alpha)*cos(beta)*cos(gamma)-
                          cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-
                          cos(gamma)*cos(gamma))/sin(gamma);
                      
  }
  else if (ibrav!=0)
    error_exit("Invalid value of ibrav");

  real2rec(c);
  
  if (atoms_are_abs)
    addfrac(m->atoms,m->n,c->recip);
  else
    addabs(m->atoms,m->n,c->basis);

}



static int qe_readline(char **p, FILE* infile){
  static char buffer[LINE_SIZE+1];
  char *ptr;
  
  if (!fgets(buffer,LINE_SIZE+1,infile)) return 0;

  /* anything after a # or ! is a comment and should be removed */

  ptr=buffer;
  for(;*ptr;ptr++)
    if ((*ptr=='#')||(*ptr=='!')||(*ptr=='\n')) {*ptr=0; break;}

  /* eat leading spaces */
  *p=buffer;
  while((**p)&&(isspace(**p))) (*p)++;

  if ((**p)==0) return qe_readline(p,infile);

  return 1;

}
