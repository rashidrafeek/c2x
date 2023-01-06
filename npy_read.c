/* Read a numpy array file, single density only, as floats or doubles 
* Cope with Fortran or C ordering, and either endianness */

/* MJR 12/2020 */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "c2xsf.h"



void npy_read(FILE* infile, struct grid *gptr){
/* NB "fortran" is an optionally-reserved word in C99, so fortrn here */
  int i,j,k,hdr_len,version,float_len,ndata,fortrn,data_little_endian;
  char magic[6]={0x93,'N','U','M','P','Y'};
  unsigned char buff[6];
  char *hdr,*cptr;
  void *data;
  double *dptr;

  data=NULL;
  
  if (fread(buff,1,6,infile)!=6)
    error_exit("Failed to read magic number");

  if (strncmp((char*)buff,magic,6))
    error_exit("Wrong magic number for .npy file");

  fread(buff,1,2,infile);
  version=buff[0];

  if (version==1){
    fread(buff,1,2,infile);
    hdr_len=buff[0]+(buff[1]<<8);
  }
  else{
    fread(buff,1,4,infile);
    hdr_len=buff[0]+(buff[1]<<8)+(buff[2]<<16)+(buff[3]<<24);
  }

  hdr=malloc(hdr_len+1);
  if (!hdr) error_exit("Malloc error for header");
  fread(hdr,1,hdr_len,infile);
  hdr[hdr_len]=0;

  cptr=strstr(hdr,"'shape'");
  if (!cptr) error_exit("Failed to find array shape");
  cptr=strchr(cptr,'(');
  if (!cptr) error_exit("Failed to parse array shape");
  cptr++;
  if (sscanf(cptr,"%d, %d, %d",gptr->size,gptr->size+1,gptr->size+2)!=3){
    fprintf(stderr,"%s\n",cptr);
    error_exit("Failed to parse dimensions from array shape");
  }
  ndata=gptr->size[0]*gptr->size[1]*gptr->size[2];
  
  cptr=strstr(hdr,"'descr'");
  if (!cptr) error_exit("Failed to find array datatype");
  cptr=strchr(cptr,':');
  if (!cptr) error_exit("Failed to parse array datatype");
  cptr++;
  while(*cptr==' ') cptr++;
  if ((*cptr!='\'')&&(*cptr!='"')){
    fprintf(stderr,"%s\n",cptr);
    error_exit("Failed to parse datatype");
  }
  cptr++;
  /* Ignore endianness for now */
  data_little_endian=self_little_endian(); /* start by assuming native */
  if (*cptr=='=')
    cptr++;
  else if (*cptr=='<') {
    data_little_endian=1;
    cptr++;
  }
  else if (*cptr=='>') {
    data_little_endian=0;
    cptr++;
  }
  
  if (*cptr!='f') error_exit("Unexpected datatype in descr");
  cptr++;
  float_len=4;
  sscanf(cptr,"%d",&float_len);
  if ((float_len!=4)&&(float_len!=8))
    error_exit("Unexpected float len in descr");

  fortrn=0;
  cptr=strstr(hdr,"'fortrn_order'");
  if (cptr){
    cptr=strchr(cptr,':');
    if (cptr){
      cptr++;
      while(*cptr==' ') cptr++;
      if (!strncmp(cptr,"True",4)) fortrn=1;
      else if (!strncmp(cptr,"False",5)) fortrn=0;
      else error_exit("Failed to parse fortrn_order");
    }
  }
  
  if (debug>1)
    fprintf(stderr,"Reading numpy array of %s size %dx%dx%d\n",
	    (float_len==4)?"floats":"doubles",
	    gptr->size[0],gptr->size[1],gptr->size[2]);
  
  data=malloc(float_len*ndata);
  if (!data) error_exit("Malloc error for data");
  
  if (fread(data,float_len,ndata,infile)!=ndata)
    error_exit("Read error reading array data");

  if (data_little_endian!=self_little_endian()){
    if (debug) fprintf(stderr,"Reversing endianness\n");
    if (float_len==4) reverse4n(data,ndata);
    else reverse8n(data,ndata);
  }
  
  gptr=grid_new(gptr);

  if ((!fortrn)&&(float_len==8)){
    gptr->data=(double*)data;
    return;
  }

  gptr->data=malloc(ndata*sizeof(double));
  if (!gptr->data) error_exit("Malloc error for data");

  
  if (float_len==8){
    dptr=gptr->data;
    for(i=0;i<gptr->size[0];i++)
      for(j=0;j<gptr->size[1];j++)
	for(k=0;k<gptr->size[2];k++)
	  *(dptr++)=((double*)data)[i+j*gptr->size[0]+
				    k*gptr->size[0]*gptr->size[1]];
  }
  else{
    if (!fortrn)
      for(i=0;i<ndata;i++)
	gptr->data[i]=((float*)data)[i];
    else{
      dptr=gptr->data;
      for(i=0;i<gptr->size[0];i++)
	for(j=0;j<gptr->size[1];j++)
	  for(k=0;k<gptr->size[2];k++)
	    *(dptr++)=((float*)data)[i+j*gptr->size[0]+
				     k*gptr->size[0]*gptr->size[1]];

    }
  }
  
  free(data);
  
}

int self_little_endian(){
  const float x=1;
  unsigned char *ptr;

  ptr=(unsigned char*)(&x);
  if (*ptr==0) return 1;
  if (*ptr==63) return 0;

  fprintf(stderr,"Endianness check failed. Assuming little\n");
  return 1;
  
  return 1;
}
