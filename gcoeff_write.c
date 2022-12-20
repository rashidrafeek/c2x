/* Produce output in the style of
 * https://www.andrew.cmu.edu/user/feenstra/wavetrans/ 
 *
 * (c) 2020 MJ Rutter
 */

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>

#include "c2xsf.h"

static int gcmp(const void *a, const void *b){
  struct gcoeff { double g2; int gpt[3]; double psi_r; double psi_i;} *x,*y;

  x=(struct gcoeff*)a;
  y=(struct gcoeff*)b;

  if (x->g2>y->g2) return 1;
  if (x->g2<y->g2) return -1;
  return 0;
}

void gcoeff_write(double *psi, int *pwgrid, int nplwv, int fft[3],
		  int gamma, int goff[3],
		  struct unit_cell *c, struct contents *m, struct kpts *kp,
		  int ikpt, int isppol, int nb,
		  double *eval, double occ, struct es *e){
  int i,j,ii,nkpt,nspins,gpt[3],gpt2[3],sort,okay;
  FILE *outfile;
  double *kpt,gvec[3];
  char *order,required_order[5];
  static int firstcall=1;
  struct gcoeff { double g2; int gpt[3]; double psi_r; double psi_i;} *band;

  kpt=kp->kpts[ikpt].frac;
  outfile=dict_get(m->dict,"out_file_handle");

  sort=0;
  if (dict_get(m->dict,"gcoeff_sort")) sort=1;
  
  if (!outfile) error_exit("No output file in gcoeff_write");
  
  if (firstcall==1){ /* write header */
    firstcall=0;
    /* count spins */
    nspins=0;
    for(i=0;i<e->nspins;i++)
      if (inrange(i,e->spin_range)) nspins++;
    /* Count kpoints */
    nkpt=0;
    for(i=0;i<kp->n;i++)
      if (inrange(i+1,e->kpt_range)) nkpt++;
    /* Check input will be in correct order */
    if ((order=dict_get(m->dict,"band_read_order"))){
      required_order[0]=0;
      if (nspins>1)
        strcat(required_order,"s");
      if (nkpt>1)
        strcat(required_order,"k");
      strcat(required_order,"b");
      okay=1;
      for(i=0;i<strlen(required_order)-1;i++){
        if (strchr(order,required_order[i])==NULL) okay=0;
        for(j=i+1;j<strlen(required_order);j++){
          if(strchr(order,required_order[i])>strchr(order,required_order[j]))
            okay=0;
        }
      }
      if (!okay){
        fprintf(stderr,"Need input in order %s, have %s, exiting\n",
                required_order,order);
        exit(1);
      }
    }
    fprintf(outfile,"      %d\n",nspins);
    fprintf(outfile," %6d\n",nkpt);
    fprintf(outfile," %6d",e->nbands);
    if (e->cut_off){
      fprintf(outfile,"   %14.8f",e->cut_off);
      if (e->e_fermi)
        fprintf(outfile,"   %14.8f",*(e->e_fermi));
    }
    fprintf(outfile,"\n");
    for(i=0;i<3;i++)
      fprintf(outfile,"%14.8f  %14.8f  %14.8f\n",c->basis[i][0],
              c->basis[i][1],c->basis[i][2]);
    for(i=0;i<3;i++)
      fprintf(outfile,"%14.8f  %14.8f  %14.8f\n",2*M_PI*c->recip[i][0],
              2*M_PI*c->recip[i][1],2*M_PI*c->recip[i][2]);
  }

  if (nb==1)
    fprintf(outfile,"%14.8f %14.8f %14.8f\n",kpt[0],kpt[1],kpt[2]);

  /* Need to count nplwv carefully if expanding gamma storage */

  if (gamma==0)
    fprintf(outfile,"        %4d      %6d\n",nb,nplwv);
  else{
    ii=nplwv;
    for(i=0;i<nplwv;i++){
      for(j=0;j<3;j++){
        gpt[j]=fft[j]-pwgrid[3*i+j]-goff[j];;
        if (gpt[j]>fft[j]/2) gpt[j]-=fft[j];
        gpt2[j]=pwgrid[3*i+j];
        if (gpt2[j]>fft[j]/2) gpt2[j]-=fft[j];
      }
      if ((gpt[0]==gpt2[0])&&(gpt[1]==gpt2[1])&&(gpt[2]==gpt2[2])) continue;
      ii++;
    }
    fprintf(outfile,"        %4d      %6d\n",nb,ii);
  }
    
  fprintf(outfile,"( %8g  , %8g  )  %8g\n",eval[0],eval[1],occ);

  if (!sort){
    for(i=0;i<nplwv;i++){
      for(j=0;j<3;j++){
        gpt[j]=pwgrid[3*i+j];
        if (gpt[j]>fft[j]/2) gpt[j]-=fft[j];
      }
      fprintf(outfile,"%6d %6d %6d  (  % .6e , % .6e )\n",gpt[0],
              gpt[1],gpt[2],psi[2*i],psi[2*i+1]);
      if (gamma){
        for(j=0;j<3;j++){
          gpt2[j]=fft[j]-pwgrid[3*i+j]-goff[j];
          if (gpt2[j]>fft[j]/2) gpt2[j]-=fft[j];
        }
        if (!((gpt[0]==gpt2[0])&&(gpt[1]==gpt2[1])&&(gpt[2]==gpt2[2])))
          fprintf(outfile,"%6d %6d %6d  (  % .6e , % .6e )\n",gpt2[0],
                  gpt2[1],gpt2[2],psi[2*i],-psi[2*i+1]);
      }
    }
  }
  else{
    band=malloc(nplwv*sizeof(struct gcoeff));
    if (!band) error_exit("malloc error for band in gcoeff_write");
    for(i=0;i<nplwv;i++){
      for(j=0;j<3;j++){
        gpt[j]=pwgrid[3*i+j];
        if (gpt[j]>fft[j]/2) gpt[j]-=fft[j];
      }
      for(j=0;j<3;j++)
        gvec[j]=(gpt[0]+kpt[0])*c->recip[0][j]+
                (gpt[1]+kpt[1])*c->recip[1][j]+
                (gpt[2]+kpt[2])*c->recip[2][j];
      band[i].g2=vmod2(gvec);
      for(j=0;j<3;j++)
        band[i].gpt[j]=gpt[j];
      band[i].psi_r=psi[2*i];
      band[i].psi_i=psi[2*i+1];
    }
    if (gamma){
      band=realloc(band,2*nplwv*sizeof(struct gcoeff));
      if (!band) error_exit("realloc error for band in gcoeff_write");
      ii=nplwv;
      for(i=0;i<nplwv;i++){
        for(j=0;j<3;j++){
          gpt[j]=fft[j]-pwgrid[3*i+j]-goff[j];;
          if (gpt[j]>fft[j]/2) gpt[j]-=fft[j];
          gpt2[j]=pwgrid[3*i+j];
          if (gpt2[j]>fft[j]/2) gpt2[j]-=fft[j];
        }
        if ((gpt[0]==gpt2[0])&&(gpt[1]==gpt2[1])&&(gpt[2]==gpt2[2])) continue;
        for(j=0;j<3;j++)
          gvec[j]=(gpt[0]+kpt[0])*c->recip[0][j]+
            (gpt[1]+kpt[1])*c->recip[1][j]+
            (gpt[2]+kpt[2])*c->recip[2][j];
        band[ii].g2=vmod2(gvec);
        for(j=0;j<3;j++)
          band[ii].gpt[j]=gpt[j];
        band[ii].psi_r=psi[2*i];
        band[ii].psi_i=psi[2*i+1];
        ii++;
      }
      nplwv=ii;
    }
        
      
    qsort(band,nplwv,sizeof(struct gcoeff),gcmp);
    for(i=0;i<nplwv;i++)
      fprintf(outfile,"%6d %6d %6d  (  % .6e , % .6e )\n",band[i].gpt[0],
              band[i].gpt[1],band[i].gpt[2],band[i].psi_r,band[i].psi_i);
      
    free(band);
  }
}
