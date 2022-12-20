/* Write a 1D file. */

/* Copyright (c) 2014-2018 MJ Rutter 
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
#include<ctype.h>
#include<math.h>
#include<string.h>

#include "c2xsf.h"

/* line_spec must be AtN:AtN:n
 * with At - atomic symbol of element present
 *      N  - index of required atom within species (1 assumed if omitted)
 *      n  - number of points
 */

void interpolate1d(struct grid *gptr, double st[3], double end[3],
                   int n, double *pts);

/* lscan: parse "(x,y,z)" or "AtN" or "At"
 *        expects termination by : or null
 *        used by main() for point value, and this file for line
 */
 
void lscan(char **p, struct contents *m, double x[3]){
  char *ptr,sym[5];
  int i,j,atno,ns,found;

  ptr=*p;

  if (*ptr=='('){
    ptr++;
    if (sscanf(ptr,"%lf,%lf,%lf",x,x+1,x+2)!=3){
      fprintf(stderr,"Malformed line_spec: %s\n",ptr);
      exit(1);
    }
    while(*ptr&&(*ptr!=')')) ptr++;
    if (*ptr!=')'){
      fprintf(stderr,"Malformed line_spec: %s\n",ptr);
      exit(1);
    }
    ptr++;
  }
  else if ((*ptr=='0')&&((*(ptr+1)==':')||(*(ptr+1)==0))){
    ptr++;
    x[0]=x[1]=x[2]=0;
  }
  else{
    for(i=0;i<4;i++){
      if(isalpha(*ptr)) sym[i]=*ptr;
      else break;
      ptr++;
    }

    sym[i]=0;
    if ((!isdigit(*ptr))&&(*ptr!=':')&&(*ptr!=0)){
      fprintf(stderr,"Malformed line_spec: %s\n",ptr);
      exit(1);
    }
    atno=atsym2no(sym);
    if (atno==0) error_exit("Invalid atom in line_spec");
    if ((*ptr==':')||(*ptr==0))
      ns=1;
    else
      if(sscanf(ptr,"%d",&ns)!=1)
	error_exit("Invalid number in atom spec in line_spec");

    found=0;
    for(i=0;i<m->n;i++){
      if (m->atoms[i].atno==atno){
        found++;
        if (found==ns){
          for(j=0;j<3;j++) x[j]=m->atoms[i].frac[j];
          break;
        }
      }
    }
    if (found!=ns) error_exit("Failed to find atom in line_spec");
    while(*ptr&&(*ptr!=':')) ptr++;
  }

  *p=ptr;

}

void line_write(FILE* outfile, struct unit_cell *c,
                struct contents *m, struct grid *gptr, char *line_spec){
  char *ptr;
  int i,j,n;
  double start[3],end[3],v[3],*points,len;
  struct grid *g;
  char *fmt;
  struct cmt *comment;

  if (!gptr||(!gptr->data))
    error_exit("No data found to plot");

  if (flags&HIPREC)
    fmt="%.12g %.12g\n";
  else
    fmt="%8g %g\n";

  /* parse line spec */


  ptr=line_spec;

  if (!strcmp(line_spec,"a")){
    start[0]=start[1]=start[2]=0;
    end[0]=1;
    end[1]=end[2]=0;
    n=gptr->size[0]+1;
  }
  else if (!strcmp(line_spec,"b")){
    start[0]=start[1]=start[2]=0;
    end[1]=1;
    end[0]=end[2]=0;
    n=gptr->size[1]+1;
  }
  else if (!strcmp(line_spec,"c")){
    start[0]=start[1]=start[2]=0;
    end[2]=1;
    end[0]=end[1]=0;
    n=gptr->size[2]+1;
  }
  else{
  
    lscan(&ptr,m,start);

    if (*ptr!=':') error_exit("Failed to find first colon in line_spec");
    ptr++;

    lscan(&ptr,m,end);

    if (*ptr!=':') error_exit("Failed to find second colon in line_spec");
    ptr++;

    if(sscanf(ptr,"%d",&n)!=1)
      error_exit("Invalid number of points in line_spec");
  }
  
  if (debug)
    fprintf(stderr,"Requested line (%f,%f,%f) to (%f,%f,%f) with %d points.\n",
            start[0],start[1],start[2],end[0],end[1],end[2],n);

  /* Find length of line */

  /* Convert to absolute co-ords */

  for(i=0;i<3;i++){
    v[i]=0;
    for(j=0;j<3;j++)
      v[i]+=(end[j]-start[j])*c->basis[j][i];
  }

  len=0;
  for(i=0;i<3;i++) len+=v[i]*v[i];
  len=sqrt(len);
  if (flags&AU) len=len/BOHR;

  points=malloc(n*sizeof(double));
  if(!points) error_exit("Malloc error for points in line_write");

  if (m->title) fprintf(outfile,"# %s\n\n",m->title);

  fprintf(outfile,"# %s\n",line_spec);
  fprintf(outfile,"# (%f,%f,%f) to (%f,%f,%f) with %d points\n",
          start[0],start[1],start[2],end[0],end[1],end[2],n);
  fprintf(outfile,"# distance in %s\n",(flags&AU)?"Bohr":"Angstrom");
  if (m->comment->txt){
    fprintf(outfile,"\n");
    comment=m->comment;
    while((comment)&&(comment->txt)){
      fprintf(outfile,"# %s\n",comment->txt);
      comment=comment->next;
    }
    fprintf(outfile,"\n");
  }


  if (flags&GNUPLOT){
    fprintf(outfile,"set xlabel \"%s\"\n",(flags&AU)?"Bohr":"Angstrom");
    fprintf(outfile,"set title \"%s (length %g %s)\"\n",line_spec,len,
            (flags&AU)?"Bohr":"A");
    fprintf(outfile,"plot [0:%g] ",len);
    g=gptr;
    while(g&&(g->data)){
      fprintf(outfile,"'-' w lp title \"%s\"",g->name);
      g=g->next;
      if (g&&(g->data)) fprintf(outfile,",");
    }
    fprintf(outfile,"\n");
  }

  while(gptr&&(gptr->data)){
    interpolate1d(gptr,start,end,n,points);

    fprintf(outfile,"# %s\n",gptr->name);

    for(i=0;i<n;i++) fprintf(outfile,fmt,i*len/(n-1),points[i]);

    if (flags&GNUPLOT) fprintf(outfile,"e\n");

    gptr=gptr->next;

    if (gptr&&(gptr->data)) fprintf(outfile,"\n");
  }

  if (flags&GNUPLOT) fprintf(outfile,"pause -1 \"Press return to exit\"\n");

}
