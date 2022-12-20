
/* Analyse a complex wavefunction */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "c2xsf.h"

/* When comparing angles, the natural scale is pi */
#define aeqa(a,b) (fabs((a)-(b))<(3.14*tol))

int pfit(double *out, double *phase, double *a, int fft[3], int dir){
  int i,imin,imax;
  int j,jmin,jmax;
  int k,kmin,kmax;
  int ok,ind,ind2,off;
  double wph1,wph2,wph3;
  double w,wt,ph,t;

  imin=jmin=kmin=0;
  imax=fft[0];
  jmax=fft[1];
  kmax=fft[2];
  if (dir==0){
    imin=1;
    imax=fft[0]/2+1;
    off=fft[2]*fft[1];
  }
  if (dir==1){
    jmin=1;
    jmax=fft[1]/2+1;
    off=fft[2];
  }
  if (dir==2){
    kmin=1;
    kmax=fft[2]/2+1;
    off=1;
  }
  
  wt=0;
  wph1=wph2=wph3=0;
  for(i=imin;i<imax;i++){
    for(j=jmin;j<jmax;j++){
      for(k=kmin;k<kmax;k++){
        ind=k+fft[2]*(j+i*fft[1]);
        ind2=ind-off;
        w=a[ind]*a[ind2];
        wt+=w;
        ph=phase[ind]-phase[ind2];
        
        if (ph>M_PI) ph-=2*M_PI;
        if (ph<-M_PI) ph+=2*M_PI;
        wph1+=w*ph;
        if (ph<-M_PI/2) ph+=2*M_PI;
        wph2+=w*ph;
        if (ph<0) ph+=2*M_PI;
        wph3+=w*ph;
      }
    }
  }

  wph1/=wt;
  wph2/=wt;
  wph3/=wt;

  ok=1;
  if (aeqa(wph1,wph2)) ph=0.5*(wph1+wph2);
  else{
    if (aeqa(wph2,wph3)) ph=0.5*(wph2+wph3);
    else{
      t=wph3;
      if (t>M_PI) t-=2*M_PI;
      if (aeqa(wph1,t)) ph=0.5*(wph1+t);
      else{
        ok=0;
        if (debug>2)
          fprintf(stderr,"Pfit failing: %d %lf %lf %lf\n",
                  dir,wph1,wph2,wph3);
      }
    }
  }
  if (ph>M_PI) ph-=2*M_PI;
  *out=ph;
  
  return ok;

}
  

void inv_parity(double *d,int fft[3],int band,double kpt[3]){
  double w,wt;
  int size,ind,ind_inv;
  int i,imin,imax,ii,ii_inv;
  int j,jmin,jmax,jj,jj_inv;
  int k,kmin,kmax,kk,kk_inv;
  int ishift,jshift,kshift;
  int ok;
  double *phi,*a;
  double x,y,z,ph,ph1,ph2,ph3;
  double old_tol,residual;
  static double old_kpt[3]={99,99,99};
  
  if (debug>2) fprintf(stderr,"Band parity called\n");

  if ((kpt[0]!=old_kpt[0])||(kpt[1]!=old_kpt[1])||(kpt[2]!=old_kpt[2])){
    fprintf(stderr,"Band parity called for k=(%f,%f,%f)\n",
            kpt[0],kpt[1],kpt[2]);
    for(i=0;i<3;i++) old_kpt[i]=kpt[i];
  }
  
  old_tol=tol;
  tol=sqrt(tol);
  
  size=fft[0]*fft[1]*fft[2];

  phi=malloc(size*sizeof(double));
  a=malloc(size*sizeof(double));
  if (!phi) error_exit("Error allocating phi");
  if (!a) error_exit("Error allocating a");

  for(i=0;i<size;i++) phi[i]=0;
  for(i=0;i<size;i++) a[i]=0;

  /* k=0, n odd: g=0 maps to itself, each other point has its own inverse
   * k=0, n even: g=0 and g=n/2 maps to self, other points have own inverse
   * k=0.5, n odd: each point has own inverse save g=n/2 which has none
   * k=-0.5, n odd: each point has own inverse save g=-n/2 which has none
   * k=+/-0.5, n even: each point has own inverse
   */ 

  
  imax=fft[0]/2;
  imin=(1-fft[0])/2;
  jmax=fft[1]/2;
  jmin=(1-fft[1])/2;
  kmax=fft[2]/2;
  kmin=(1-fft[2])/2;
  ok=1;

  ishift=jshift=kshift=0;
  if (aeq(kpt[0],0.5)) ishift=1;
  if (aeq(kpt[0],-0.5)) ishift=-1;
  if (aeq(kpt[1],0.5)) jshift=1;
  if (aeq(kpt[1],-0.5)) jshift=-1;
  if (aeq(kpt[2],0.5)) kshift=1;
  if (aeq(kpt[2],-0.5)) kshift=-1;
  
  if (debug>3) fprintf(stderr,"Shifts: %d %d %d\n",ishift,jshift,kshift);
  
  for(i=imin;i<=imax;i++){
    if (i>=0)
      ii=i;
    else
      ii=fft[0]+i;
    ii_inv=-i-ishift;
    if ((ii_inv<imin)||(ii_inv>imax)) continue;
    if (ii_inv<0)
      ii_inv+=fft[0];
    for(j=jmin;j<=jmax;j++){
      if (j>=0)
        jj=j;
      else
        jj=fft[1]+j;
      jj_inv=-j-jshift;
      if ((jj_inv<jmin)||(jj_inv>jmax)) continue;
      if (jj_inv<0)
        jj_inv+=fft[1];
      for(k=kmin;k<=kmax;k++){
        if (k>=0)
          kk=k;
        else
          kk=fft[2]+k;
        kk_inv=-k-kshift;
        if ((kk_inv<kmin)||(kk_inv>kmax)) continue;
        if (kk_inv<0)
          kk_inv+=fft[2];
        ind=2*(kk+fft[2]*(jj+ii*fft[1]));
        ind_inv=2*(kk_inv+fft[2]*(jj_inv+ii_inv*fft[1]));

        phi[ind/2]=atan2(d[ind+1]*d[ind_inv]-d[ind]*d[ind_inv+1],
                         d[ind]*d[ind_inv]+d[ind+1]*d[ind_inv+1]);

        a[ind/2]=sqrt(d[ind]*d[ind]+d[ind+1]*d[ind+1]);

        if (!aeq(sqrt(d[ind]*d[ind]+d[ind+1]*d[ind+1]),
                 sqrt(d[ind_inv]*d[ind_inv]+d[ind_inv+1]*d[ind_inv+1])))
          ok=0;
      }
    }
  }

#if 0
  for(i=0;i<fft[0];i++){
    fprintf(stderr,"%d %f %f (%f,%f)\n",i,a[i],phi[i],d[2*i],d[2*i+1]);
  }
#endif
  
  if (!ok){
    if (debug>1) fprintf(stderr,"No parity: modulus check failed\n");
    free(phi);
    free(a);
    return;
  }
  
  if (!pfit(&x,phi,a,fft,0)) ok=0;
  if (!pfit(&y,phi,a,fft,1)) ok=0;
  if (!pfit(&z,phi,a,fft,2)) ok=0;

  if (!ok){
    if (debug>1) fprintf(stderr,"No parity: one of x, y or z failed\n");
    free(phi);
    free(a);
    return;
  }

  /* Find constant */

  wt=0;
  ph1=ph2=ph3=0;
  for(i=imin;i<=imax;i++){
    if (i>=0)
      ii=i;
    else
      ii=fft[0]+i;
    for(j=jmin;j<=jmax;j++){
      if (j>=0)
        jj=j;
      else
        jj=fft[1]+j;
      for(k=kmin;k<=kmax;k++){
        if (k>=0)
          kk=k;
        else
          kk=fft[2]+k;
        ind=kk+fft[2]*(jj+ii*fft[1]);
        ph=fmod(phi[ind]-x*(i+kpt[0])-y*(j+kpt[1])-z*(k+kpt[2]),2*M_PI);
        if (ph>M_PI) ph-=2*M_PI;
        if (ph<-M_PI) ph+=2*M_PI;
        w=a[ind]*a[ind];
        wt+=w;
        ph1+=ph*w;
        if (ph<-M_PI/2) ph+=2*M_PI;
        ph2+=ph*w;
        if (ph<0) ph+=2*M_PI;
        ph3+=ph*w;
      }
    }
  }
  ph1/=wt;
  ph2/=wt;
  ph3/=wt;
  ok=1;
  if (aeqa(ph1,ph2)) ph=0.5*(ph1+ph2);
  else{
    if (aeqa(ph2,ph3)) ph=0.5*(ph2+ph3);
    else{
      if (aeqa(ph1,ph3)) ph=0.5*(ph1+ph3);
      else{
        ok=0;
        if (debug>2)
          fprintf(stderr,"Const failing: %lf %lf %lf\n",ph1,ph2,ph3);
      }
    }
  }

  if (!ok){
    if (debug>1) fprintf(stderr,"No parity: constant failed\n");
    free(phi);
    free(a);
    return;
  }

  if (debug>3) fprintf(stderr,"Const: %lf %lf %lf\n",ph1,ph2,ph3);
  
  residual=0;
  wt=0;
  for(i=imin;i<=imax;i++){
    if (i>=0)
      ii=i;
    else
      ii=fft[0]+i;
    for(j=jmin;j<=jmax;j++){
      if (j>=0)
        jj=j;
      else
        jj=fft[1]+j;
      for(k=kmin;k<=kmax;k++){
        if (k>=0)
          kk=k;
        else
          kk=fft[2]+k;
        ind=kk+fft[2]*(jj+ii*fft[1]);
        ph1=fmod(x*(i+kpt[0])+y*(j+kpt[1])+z*(k+kpt[2])+ph-phi[ind],2*M_PI);
        if (ph1>M_PI) ph1-=2*M_PI;
        if (ph1<-M_PI) ph1+=2*M_PI;
        w=a[ind]*a[ind];
        wt+=w;
        residual+=ph1*ph1*w;
      }
    }
  }
  residual=sqrt(residual/wt);

  if (residual<tol){
    fprintf(stderr,"Band %3d: ",band);
    if (aeqa(ph,0)){
      fprintf(stderr,"inversion with even parity");
    }
    else if (aeqa(fabs(ph),M_PI)){
      fprintf(stderr,"inversion with odd parity");
    }
    else
      fprintf(stderr,"unexpected mixed parity of %f",ph);
    fprintf(stderr," about (%f,%f,%f)\n",
            -0.25*x/M_PI,-0.25*y/M_PI,-0.25*z/M_PI);
    if (debug>1)
      fprintf(stderr,"(residual %f)\n",residual);
  }
  else{
    if (debug>1)
      fprintf(stderr,
              "Band %3d: residual of %f too large to identify symmetry\n",
              band,residual);
  }

  free(phi);
  free(a);
  tol=old_tol;
  
}
