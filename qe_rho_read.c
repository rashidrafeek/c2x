

#include <stdio.h>
#include <stdlib.h>
#include "c2xsf.h"

void fft3d(double *c, int *ngptar, int dir);
void to235(int *i); /* from super.c */

void qe_rho_read(FILE* infile, struct unit_cell *c, struct contents *m,
                 struct kpts *k, struct symmetry *s, struct grid *g,
		 struct es *elect, int *i_grid){
  int i,j,tmp;
  int gamma_only,ngm_g,nspin;
  int *mill_g;
  double *rho_g,recip[3][3],*data,*data2;
  int ngx_min,ngx_max,ngy_min,ngy_max,ngz_min,ngz_max;
  int ngx,ngy,ngz,offset,fft[3];
  int nx,ny,nz;
  
  if (debug) fprintf(stderr,"QE binary reader called\n");

  fread(&tmp,4,1,infile);
  if (tmp!=12) error_exit("Unexpected start to QE binary file");

  fread(&gamma_only,4,1,infile);
  fread(&ngm_g,4,1,infile);
  fread(&nspin,4,1,infile);

  if (nspin>2) error_exit("nspin>2 not supported");
  if ((nspin==1)&&(flags&SPINDEN))
    error_exit("Spin density requested, but nspin=1 in density file");

  if (debug) fprintf(stderr,"QE density file with %d spins\n",nspin);

  elect->nspins=nspin;
  
  fread(&tmp,4,1,infile);
  if (tmp!=12) error_exit("Unexpected data in QE binary file");

  fread(&tmp,4,1,infile);
  if (tmp!=72) error_exit("Unexpected data in QE binary file");

  fread(recip,72,1,infile);

  fread(&tmp,4,1,infile);
  if (tmp!=72) error_exit("Unexpected data in QE binary file");

  mill_g=malloc(3*ngm_g*sizeof(int));
  if (!mill_g) error_exit("Malloc error for mill_g in qe_bin_read()");

  fread(&tmp,4,1,infile);
  if (tmp!=3*ngm_g*4) error_exit("Unexpected data in QE binary file");
  fread(mill_g,3*ngm_g*4,1,infile);
  fread(&tmp,4,1,infile);
  if (tmp!=3*ngm_g*4) error_exit("Unexpected data in QE binary file");

  rho_g=malloc(ngm_g*2*sizeof(double));  /* rho_g is complex */
  if (!rho_g) error_exit("Malloc error for rho_g in qe_bin_read()");

  fread(&tmp,4,1,infile);
  if (tmp!=ngm_g*2*sizeof(double)) error_exit("Unexpected data in QE binary file");
  fread(rho_g,ngm_g*16,1,infile);
  fread(&tmp,4,1,infile);
  if (tmp!=ngm_g*2*sizeof(double)) error_exit("Unexpected data in QE binary file");

  if ((nspin==2)&&(flags&SPINDEN)){
    fread(&tmp,4,1,infile);
    if (tmp!=ngm_g*2*sizeof(double))
      error_exit("Unexpected data in QE binary file");
    fread(rho_g,ngm_g*16,1,infile);
    fread(&tmp,4,1,infile);
    if (tmp!=ngm_g*2*sizeof(double))
      error_exit("Unexpected data in QE binary file");
  }
  
  if (debug>1) fprintf(stderr,"Data file read\n");

  /* We have reciprocal axes, in odd units. First convert the units */

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      recip[i][j]/=2*M_PI;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      recip[i][j]/=BOHR;

  /* Given that a real2rec routine and a rec2real routine would be identical,
     abuse existing real2rec... */
  
  c->basis=malloc(72);
  if (!c->basis) error_exit("Malloc error for basis in qe_bin_read()");
  
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[i][j]=recip[i][j];

  real2rec(c);

  /* Now we need to exchange real and recip, and correct the cell volume */
  
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      c->basis[i][j]=c->recip[i][j];
  
  real2rec(c);

  /* Next work out our FFT grid size */

  ngx_min=ngx_max=ngy_min=ngy_max=ngz_min=ngz_max=0;

  for(i=0;i<ngm_g;i++){
    if (mill_g[3*i]<ngx_min) ngx_min=mill_g[3*i];
    if (mill_g[3*i]>ngx_max) ngx_max=mill_g[3*i];
    if (mill_g[3*i+1]<ngy_min) ngy_min=mill_g[3*i+1];
    if (mill_g[3*i+1]>ngy_max) ngy_max=mill_g[3*i+1];
    if (mill_g[3*i+2]<ngz_min) ngz_min=mill_g[3*i+2];
    if (mill_g[3*i+2]>ngz_max) ngz_max=mill_g[3*i+2];
  }

  if (debug) fprintf(stderr,"grid extent in file %d:%d,%d:%d,%d:%d\n",
	  ngx_min,ngx_max,ngy_min,ngy_max,ngz_min,ngz_max);  

  ngx=ngx_max-ngx_min+1;
  to235(&ngx);
  ngy=ngy_max-ngy_min+1;
  to235(&ngy);
  ngz=ngz_max-ngz_min+1;
  to235(&ngz);
  
  if (i_grid){
    if (i_grid[0]!=0) ngx=i_grid[0];
    if (i_grid[1]!=0) ngy=i_grid[1];
    if (i_grid[2]!=0) ngy=i_grid[2];
  }
  else{
    if (debug) fprintf(stderr,"FFT grid used %dx%dx%d\n",ngx,ngy,ngz);
  }
  ngx_max=ngx/2;
  ngx_min=-(ngx-1)/2;
  ngy_max=ngy/2;
  ngy_min=-(ngy-1)/2;
  ngz_max=ngz/2;
  ngz_min=-(ngz-1)/2;
    
  
  data=malloc(ngx*ngy*ngz*16);

  if (!data) error_exit("Error in malloc for data in qe_bin_read");

  for(i=0;i<ngx*ngy*ngz*2;i++) data[i]=0.0;

  for(i=0;i<ngm_g;i++){
    nx=mill_g[3*i];
    if (nx>ngx_max) continue;
    if (nx<ngx_min) continue;
    if (nx<0) nx+=ngx;
    ny=mill_g[3*i+1];
    if (ny>ngy_max) continue;
    if (ny<ngy_min) continue;
    if (ny<0) ny+=ngy;
    nz=mill_g[3*i+2];
    if (nz>ngz_max) continue;
    if (nz<ngz_min) continue;
    if (nz<0) nz+=ngz;
    offset=nz+ngz*(ny+ngy*nx);
    data[2*offset]=rho_g[2*i];
    data[2*offset+1]=rho_g[2*i+1];
  }

  fft[0]=ngz;
  fft[1]=ngy;
  fft[2]=ngx;

  fft3d(data,fft,1);

  data2=malloc(ngx*ngy*ngz*8);

  if (!data2) error_exit("Error in malloc for data2 in qe_bin_read");

  for(i=0;i<ngx*ngy*ngz;i++)
    data2[i]=data[2*i]/(BOHR*BOHR*BOHR);

  free(data);

  g->size[0]=ngx;
  g->size[1]=ngy;
  g->size[2]=ngz;
  g->data=data2;
  if (flags&SPINDEN)
    g->name="Spin";
  else
    g->name="Charge";
  g->next=malloc(sizeof(struct grid));
  if (!g->next) error_exit("malloc failure for struct grid");
  g=g->next;
  g->data=NULL;
  g->next=NULL; 
  
}
