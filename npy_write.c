/* Write a numpy array file, single density only, as floats or doubles */

/* MJR 11/2020 */

#include<stdio.h>
#include<stdlib.h>

#include "c2xsf.h"

#define MAX_HDR_LEN 256

void npy_write(FILE* outfile, struct grid *g){
  char magic[6]={0x93,'N','U','M','P','Y'};
  char *hdr,c;
  float *x;
  int i,ndata;

  if ((!g)||(!g->data)) error_exit("No grid data to write");
  
  if (fwrite(magic,1,6,outfile)!=6)
    error_exit("Error writing magic number");

  c=1;
  fwrite(&c,1,1,outfile);
  c=0;
  fwrite(&c,1,1,outfile);

  hdr=malloc(MAX_HDR_LEN);
  for(i=0;i<MAX_HDR_LEN;i++) hdr[i]=' ';

  i=snprintf(hdr,MAX_HDR_LEN,"{'descr': '%sf%d', 'fortran_order': False,"
	     " 'shape': (%d, %d, %d) }",(self_little_endian())?"<":">",
	     (flags&HIPREC)?8:4,
	     g->size[0],g->size[1],g->size[2]);

  /* Remove terminating null */
  hdr[i]=' ';

  /* Some docs say round to multiple of 16, some to multiple of 64.
   * Here 64 is used, so round i+10 up to multiple of 64 */
  i=i+10;
  i=i+(64-(i&63));
  i=i-10;
  hdr[i-1]='\n';

  c=i&0xff;
  fwrite(&c,1,1,outfile);
  c=(i&0xff00)>>8;
  fwrite(&c,1,1,outfile);
  
  fwrite(hdr,1,i,outfile);
  free(hdr);
  
  ndata=g->size[0]*g->size[1]*g->size[2];
  if (flags&HIPREC)
    i=fwrite(g->data,sizeof(double),ndata,outfile);
  else{
    x=malloc(ndata*sizeof(float));
    if (!x) error_exit("Malloc error in npy_write");
    for(i=0;i<ndata;i++) x[i]=g->data[i];
    i=fwrite(x,sizeof(float),ndata,outfile);
    free(x);
  }
    

  if (i!=ndata)
    error_exit("Error writing array");

}
