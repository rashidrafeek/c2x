/* A very simplistic SHELX reader */


/* Copyright (c) 2013 MJ Rutter 
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


/* Does not handle continuation lines.
 *
 * Parses CELL and SFAC "commands"
 *
 * Recognises all(?) other commands, so what remains must be an atom
 * label, and is parsed as same.
 */

#include<stdio.h>
#include<stdlib.h> /* malloc */
#include<string.h>
#include<errno.h>
#include<ctype.h>

#include "c2xsf.h"

#define LINE_SIZE 200

void shelx_read(FILE* infile, struct unit_cell *c, struct contents *m){
  double abc[6],dummy;
  int have_basis=0,i,cnt,is_cmd;
  int spec[102],nspec=1;
  char buffer[LINE_SIZE+1],*cptr,*cptr2;

  char *cmd[]={"TITL","CELL","ZERR","LATT","SFAC","UNIT",
	       "SYMM","DISP","LAUE","REM","MORE","TIME","END","HKLF",
	       "OMIT","SHEL","BASF","TWIN","EXTI","SWAT","HOPE","MERG",
	       "SPEC","RESI","MOVE","ANIS","AFIX","HFIX","FRAG","FEND",
	       "EXYZ","EADP","EQIV","CONN","PART","BIND","DFIX","DANG",
	       "BUMP","SAME","SADI","CHIV","FLAT","DELU","SIMU","DEFS",
	       "ISOR","NCSY","SUMP","L.S.","CGLS","BLOC","DAMP","STIR",
	       "WGHT","FVAR","BOND","CONF","MPLA","RTAB","HTAB","LIST",
	       "ACTA","SIZE","TEMP","WPDB","FMAP","GRID","PLAN","MOLE",0};

  m->n=0;

  if (debug>2) fprintf(stderr,"shelx_read called\n");

  if (!(c->basis=malloc(72))) error_exit("Malloc error in shelx_read for c->basis");

  while(1){
    for(i=0;i<LINE_SIZE;i++) buffer[i]=0;
    if (!fgets(buffer,LINE_SIZE-1,infile)) break;

    cptr=buffer;

    if (*cptr==0) continue;   /* Skip blank lines */
    if (*cptr=='!') continue; /* Skip comments */
    if (*cptr==' ') continue; /* Lines starting with spaces are comments */

    cptr2=cptr;
    while((*cptr2)&&(*cptr2!=' ')&&(*cptr2!='_')&&(*cptr2!='\n')
	  &&(*cptr2!='\r')&&((cptr2-cptr)<4)) cptr2++;
    *cptr2=0;

    if (!strcasecmp(cptr,"REM")) continue;
    if (!strcasecmp(cptr,"TITL")){
      cptr2++;
      cptr=cptr2;
      while((*cptr2)&&(*cptr2!='\n')&(*cptr2!='\r')) cptr2++;
      m->title=malloc(cptr2-cptr+1);
      if(!m->title) error_exit("Malloc error in shelx_read()");
      strncpy(m->title,cptr,cptr2-cptr);
      m->title[cptr2-cptr]=0;
      continue;
    }
    if (!strcasecmp(cptr,"CELL")){
      cnt=sscanf(cptr2+1,"%lf %lf %lf %lf %lf %lf %lf",&dummy,
	       abc,abc+1,abc+2,abc+3,abc+4,abc+5);
      if (cnt!=7){
        fprintf(stderr,"%s",cptr+4);
        error_exit("error parsing CELL line");
      }
      have_basis=1;
      continue;
    }
    if (!strcasecmp(cptr,"SFAC")){ /* We now have an unknown number of
				    * atomic species */
      cptr=cptr2+1;
      while(*cptr!=0){
	while(*cptr==' ') cptr++;
	cptr2=cptr;
	while(isalpha(*cptr2)) cptr2++;
	*cptr2=0;
	spec[nspec]=atsym2no(cptr);
	cptr=cptr2+1;
	nspec++;
      }
      if (debug) fprintf(stderr,"%d species in shelx file\n",nspec-1);
      continue;
    }

    /* Do we have another command? */
    is_cmd=0;
    for(i=0;cmd[i];i++){
      if (!strcasecmp(cptr,cmd[i])){
	is_cmd=1;
	break;
      }
    }

    if(!is_cmd){ /* We have an atom */
      m->atoms=realloc(m->atoms,(m->n+1)*sizeof(struct atom));
      if (!m->atoms) error_exit("realloc error in shelx_read");
      init_atoms(m->atoms+m->n,1);
      cnt=sscanf(cptr2+1,"%d %lf %lf %lf",&i,m->atoms[m->n].frac,
	       m->atoms[m->n].frac+1,m->atoms[m->n].frac+2);
      if (cnt!=4){
	fprintf(stderr,"Error parsing following line as atom.\n%s\n",buffer);
	fprintf(stderr,"Blundering on regardless\n");
      }
      else{
	if (i>=nspec) error_exit("Species number too big in shelx_read");
	m->atoms[m->n].atno=spec[i];
	m->n++;
      }
    }
  }

  if (!have_basis) error_exit("No basis in shelx file");

  abc2cart(abc,c);
 
  if (debug>1) fprintf(stderr,"%d atoms read\n",m->n);

  real2rec(c);
  addabs(m->atoms,m->n,c->basis);

}
