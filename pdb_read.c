/* A very simplistic .pdb reader */

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
#include<stdlib.h> /* malloc */
#include<string.h>
#include<errno.h>

#include "c2xsf.h"

#define LINE_SIZE 100

static int strnncpy(char *dest, char *src, int n);

void pdb_read(FILE* infile, struct unit_cell *c, struct contents *m){
  double abc[6],*dptr;
  int have_basis=0,i,j;
  char buffer[LINE_SIZE+1],buff2[LINE_SIZE],*cptr,*cptr2;

  m->n=0;

  if (debug>2) fprintf(stderr,"pdb_read called\n");

  if (!(c->basis=malloc(72)))
    error_exit("Malloc error in pdb_read for c->basis");

  while(1){
   for(i=0;i<LINE_SIZE-1;i++) buffer[i]=0;
    if (!fgets(buffer,LINE_SIZE,infile)) break;
/* First six characters are record name */
    strnncpy(buff2,buffer,6);

    if (!strcasecmp(buff2,"REMARK")) continue;
    if (!strcasecmp(buff2,"TER")) continue;
    if (!strcasecmp(buff2,"END")) break;
    if (!strcasecmp(buff2,"MTRIX1")) continue; /* MTRXn are symmetry */
    if (!strcasecmp(buff2,"MTRIX2")) continue; /* operations which can */
    if (!strcasecmp(buff2,"MTRIX3")) continue; /* be safely ignored */

    if ((!strcasecmp(buff2,"TITLE"))||(!strcasecmp(buff2,"HEADER"))){
      if (!m->title){
        m->title=malloc(strlen(buffer)-9);
        cptr=buffer+10;
        cptr2=m->title;
        while((*cptr)&&(*cptr!='\n')) *(cptr2++)=*(cptr++);
        *cptr2=0;
      }
      else add_cmt(m->comment,buffer+10);
    }
    
    if (!strcasecmp(buff2,"CRYST1")){
      sscanf(buffer+7,"%lf %lf %lf %lf %lf %lf",abc,abc+1,abc+2,
                                              abc+3,abc+4,abc+5);
/* By convention, a=b=c=1, alpha=beta=gamma=90 is a dummy entry here... */
      if ((abc[0]==1)&&(abc[1]==1)&&(abc[2]==1)&&
          (abc[3]==90)&&(abc[4]==90)&&(abc[5]==90)) continue;
      abc2cart(abc,c);
      have_basis=1;
      continue;
    }

    if ((strcasecmp(buff2,"ATOM")==0)||(strcasecmp(buff2,"HETATM")==0)){
      m->atoms=realloc(m->atoms,(m->n+1)*sizeof(struct atom));
      if (!m->atoms) error_exit("realloc error in pdb_read");
      init_atoms(m->atoms+m->n,1);
      strnncpy(buff2,buffer+30,24); /* grab co-ords section */
      dptr=m->atoms[m->n].abs;
      if (sscanf(buff2," %lf %lf %lf",dptr,dptr+1,dptr+2)!=3){
	fprintf(stderr,"Error parsing line\n%s\n",buffer);
	fprintf(stderr,"buff2 was '%s'\n",buff2);
	exit(1);
      }
/* The atomic symbol ought to be in columns 77 to 78, right justified */
      strnncpy(buff2,buffer+76,2);
      if (buff2[0])
	m->atoms[m->n].atno=atsym2no(buff2);
/* Now we have a problem. The "atom name" is in cols 13-16, but 
 * there is little convention for how this is formatted. A standard
 * says that cols 13-14 are the symbol, right justified, col 15 is
 * a letter from ABGDEZH (think Greek), and col 16 is a digit, except
 * if the atom is H, in which case it may be in column 13 if its suffix
 * is 3 characters. However, other people use symbol followed by a multidigit
 * number. So " CA1" ought to be carbon, and "CA1 " ought to be calcium,
 * but not everyone obeys this strictly. As for "HG12", it could be
 * hydrogen or mercury...
 */
#if 0
/* This was check2xsf's original code. It is not very good... */
      else {
	strnncpy(buff2,buffer+12,4);
	cptr=buff2;
	while(*cptr==' ') cptr++;
	cptr[2]=0;
	m->atoms[m->n].atno=atsym2no(cptr);
      }
#else
/* This should be better, on average... */
      else {
	strnncpy(buff2,buffer+12,2);
	cptr=buff2;
	cptr[2]=0;
        /* I have seen things (mostly H) prefixed by a digit */
	if ((*cptr>='0')&&(*cptr<='9')) *cptr=' ';
        /* If starts with an H, is H, unless He, Hf, Hg or Ho */
	if ((*cptr=='H')&&((((cptr[1]^'E')&95)!=0)||
			   (((cptr[1]^'F')&95)!=0)||
			   (((cptr[1]^'G')&95)!=0)||
			   (((cptr[1]^'O')&95)!=0))) cptr[1]=0;
        if ((cptr[1]>='0')&&(cptr[1]<='9')) cptr[1]=0;
	if ((cptr[0]==' ')&&(cptr[1]==' ')){
	  strnncpy(buff2,buffer+14,2);
	  cptr=buff2;
	}
	m->atoms[m->n].atno=atsym2no(cptr);
      }
#endif
      m->n++;
      continue;
    }

    if (debug)
      fprintf(stderr,"Warning, ignoring PDB line entitled '%s'\n",buff2);
  }

  if(have_basis==0){
    if (debug) fprintf(stderr,"Warning, no unit cell in pdb file\n"
                               "Creating dummy 10A box\n");
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
        c->basis[i][j]=0;
    for(i=0;i<3;i++) c->basis[i][i]=10;
  }

  if (debug>1) fprintf(stderr,"%d atoms read\n",m->n);

  real2rec(c);
  addfrac(m->atoms,m->n,c->recip);

}


/* Function to copy n characters, to null-terminate, and to kill trailing
 * whitespace. Destination must be n+1 characters long, unless source was
 * null terminated */

static int strnncpy(char *dest, char *src, int n){
  int i;
  for (i=0;src[i]&&i<n;i++) dest[i]=src[i];
  dest[i--]=0;
  while((i>=0)&&((dest[i]==' ')||(dest[i]=='\n')||(dest[i]=='\r'))) dest[i--]=0;
  return(i+2);
}
