/* Need to check symmetry reading for six-fold axes etc. */

/* This should be a proper XML parser. It isn't */

/* Energies here are in Hartrees, but in a .in file are in Rydbergs?! */
/* Lengths are in Bohr */

#include<stdio.h>
#include<stdlib.h> /* malloc */
#include<string.h>
#include<ctype.h>
#include<math.h>

#include "c2xsf.h"

/* The header line in a pseudo_pot can be rather long */
#define LINE_SIZE 2000

static void read_step(FILE *infile, char *buffer, struct time_series *ts);
static void read_basis(FILE *infile, struct unit_cell *c);
static void read_atoms(FILE *infile, struct contents *m);

void qe_rho_read(FILE* infile, struct unit_cell *c, struct contents *m,
		 struct kpts *k, struct symmetry *s, struct grid *g,
		 struct es *elect, int *i_grid);
void qe_psi_read(char *dir, char *prefix, struct unit_cell *c,
		 struct contents *m,
                 struct kpts *k, struct symmetry *s, struct grid *g,
                 struct es *elect, int fft[3], int *i_grid);
FILE *qe_pspot_file(char *name, char *prefix, char *cwd);
double qe_pot_charge(FILE *infile);
void kcart2frac(double xk[3], double kpt[3], double recip[3][3]);

/* If ts==NULL, don't read time series */

void qe_xml_read(FILE* infile, char *filename, struct unit_cell *c,
		 struct contents *m, struct kpts *k, struct symmetry *s,
		 struct grid *g, struct es *e, struct time_series *ts,
                 int *i_grid){
  char buffer[LINE_SIZE+1];
  char key[6];
  char *cptr,*cptr2,*cptr3,*den_name,*den_name2,*cwd;
  int i,j,jj,tmp;
  FILE *ch,*pot;
  struct sp {char *name; double mass; char *p_pot;
    double chg; double mag;} *species;
  int nspec,nsym,fft[3],*iptr,nkpt,nbands;
  double *symrel,*symtr,mat[3][3],*dptr,alat,sum;
  struct kpts *k_in;
  
  if (debug>2) fprintf(stderr,"QE xml read called\n");

  fft[0]=fft[1]=fft[2]=0;
  alat=0;
  k_in=malloc(sizeof(struct kpts));
  if (!k_in) error_exit("malloc error for struct kpts");
  k_in->n=0;
  symrel=NULL;
  symtr=NULL;
  species=NULL;
  nbands=0;
  nkpt=0;

  /* Get directory part of filename */
  cwd=NULL;
  cptr=filename+strlen(filename);
  while((cptr>filename)&&(*cptr!='/')) cptr--;
  if (*cptr=='/'){
    cwd=malloc((cptr-filename)+2);
    if (!cwd) error_exit("malloc error for struct dirname");
    strncpy(cwd,filename,(cptr-filename)+1);
    cwd[(cptr-filename)+1]=0;
    /* We need the trailing / only if the directory part is simply "/" */
    if (cptr!=filename) cwd[cptr-filename]=0;
  }
  
  nspec=0;
  
  /* Skip first line */

  while((i=fgetc(infile))!=EOF)
    if (i=='\n') break;

  if (i==EOF) error_exit("Input file contains no lines!");

  fgets(buffer,LINE_SIZE,infile);

  if (!strstr(buffer,"quantum-espresso"))
    error_exit("'quantum-espresso' not found on 2nd line of XML file");

  /* First scan for input section */
  
  while(fgets(buffer,LINE_SIZE,infile))
    if (strstr(buffer,"<input>")) break;

  while(fgets(buffer,LINE_SIZE,infile)){
    if (strstr(buffer,"</input>")) break;

    if (strstr(buffer,"<title>")){
      cptr=strstr(buffer,"<title>")+7;
      cptr2=strstr(cptr,"</title>");
      if (cptr2&&(cptr2!=cptr)){
        m->title=malloc((cptr2-cptr)+2);
	if (!m->title) error_exit("Malloc error for title");
	strncpy(m->title,cptr,(cptr2-cptr)+1);
	m->title[(cptr2-cptr)+1]=0;
      }
    }
    else if (strstr(buffer,"<prefix>")){
      cptr=strstr(buffer,"<prefix>")+8;
      cptr2=strstr(cptr,"</prefix>");
      if (cptr2&&(cptr2!=cptr)){
        cptr3=malloc((cptr2-cptr)+1);
	if (!cptr3) error_exit("Malloc error for prefix");
	strncpy(cptr3,cptr,cptr2-cptr);
	cptr3[cptr2-cptr]=0;
	dict_add(m->dict,"QE_prefix",cptr3);
      }
    }
    else if (strstr(buffer,"<pseudo_dir>")){
      cptr=strstr(buffer,"<pseudo_dir>")+12;
      cptr2=strstr(cptr,"</pseudo_dir>");
      if (cptr2&&(cptr2!=cptr)){
        cptr3=malloc((cptr2-cptr)+1);
	if (!cptr3) error_exit("Malloc error for pseudo_dir");
	strncpy(cptr3,cptr,cptr2-cptr);
	cptr3[cptr2-cptr]=0;
	dict_add(m->dict,"QE_pseudo_dir",cptr3);
      }
    }
    else if (strstr(buffer,"<calculation>")){
      cptr=strstr(buffer,"<calculation>")+13;
      cptr2=strstr(cptr,"</calculation>");
      if (cptr2&&(cptr2!=cptr)){
        cptr3=malloc((cptr2-cptr)+1);
	if (!cptr3) error_exit("Malloc error for calculation");
	strncpy(cptr3,cptr,cptr2-cptr);
	cptr3[cptr2-cptr]=0;
	dict_add(m->dict,"QE_calculation",cptr3);
      }
    }
    else if (strstr(buffer,"<forc_conv_thr>")){
      cptr=strstr(buffer,"<forc_conv_thr>")+15;
      dptr=malloc(sizeof(double));
      i=sscanf(cptr,"%lf",dptr);
      if (i) dict_add(m->dict,"QE_forc_conv_thr",dptr);
    }
    else if (strstr(buffer,"<tot_charge>")){
      cptr=strstr(buffer,"<tot_charge>")+12;
      e->charge=malloc(sizeof(double));
      i=sscanf(cptr,"%lf",e->charge);
      if (!i) {
        fprintf(stderr,"Warning: failed to parse <tot_charge>\n");
        free(e->charge);
        e->charge=NULL;
      }
    }
    else if (strstr(buffer,"<conv_thr>")){
      cptr=strstr(buffer,"<conv_thr>")+10;
      i=sscanf(cptr,"%lf",&e->etol);
      if (!i) {
        fprintf(stderr,"Warning: failed to parse <conv_thr>\n");
        e->etol=0;
      }
      else e->etol*=H_eV; /* units were Hartrees per cell */
    }
    else if (strstr(buffer,"<nosym>")){
      cptr=strstr(buffer,"<nosym>")+7;
      if (!strncmp(cptr,"true",4)){ /* record only if not the default */
	iptr=malloc(sizeof(int));
	*iptr=1;
	dict_add(m->dict,"QE_nosym",iptr);
      }
    }
    else if ((ts)&&(strstr(buffer,"<cell>"))){ /* Need to read initial cell */
      ts->cells=realloc(ts->cells,(ts->nc+1)*sizeof(struct unit_cell));
      if (!ts->cells) error_exit("Realloc error for ts->cells");
      ts->cells[ts->nc].basis=NULL;
      ts->cells[ts->nc].stress=NULL;
      read_basis(infile,ts->cells+ts->nc);
      ts->nc++;
    }
    else if ((ts)&&(strstr(buffer,"<atomic_positions>"))){
      ts->m=realloc(ts->m,(ts->nm+1)*sizeof(struct contents));
      if (!ts->m) error_exit("Realloc error for ts->m");
      ts->m[ts->nm].atoms=NULL;
      ts->m[ts->nm].title=NULL;
      ts->m[ts->nm].comment=NULL;
      ts->m[ts->nm].block_species=NULL;
      ts->m[ts->nm].species_misc=NULL;
      ts->m[ts->nm].dict=NULL;
      ts->m[ts->nm].forces=0;
      ts->m[ts->nm].n=0;
      read_atoms(infile,ts->m+ts->nm);
      ts->nm++;
    }
  }
  
  /* Now scan for output section, also reading timesteps if req */

  while(fgets(buffer,LINE_SIZE,infile)){
    if ((ts)&&(strstr(buffer,"<step n_step=")))
      read_step(infile,buffer,ts);
    if (strstr(buffer,"<output>")) break;
  }
    
  if (!strstr(buffer,"<output>")) error_exit("output section not found");


  /* Don't make any assumptions about the order of the items */
  while(fgets(buffer,LINE_SIZE,infile)){
    if (strstr(buffer,"</output>")) break;
    if(strstr(buffer,"<cell>")){
      read_basis(infile,c);
    }
    else if(strstr(buffer,"<atomic_structure ")){
      cptr=strstr(buffer,"<atomic_structure>");
      cptr=strstr(buffer,"alat=");
      cptr+=5;
      if (*cptr=='"') cptr++;
      sscanf(cptr,"%lf",&alat);
    }
    else if(strstr(buffer,"<atomic_positions>")){
      read_atoms(infile,m);
    }
    else if(strstr(buffer,"<starting_k_points>")){
      fgets(buffer,LINE_SIZE,infile);
      cptr=strstr(buffer,"<nk>");
      if (cptr) {
	sscanf(cptr+4,"%d",&k_in->n);
	k_in->kpts=malloc(k_in->n*sizeof(struct atom));
	if (!k_in->kpts) error_exit("Malloc error for kpts");
	for(i=0;i<k_in->n;i++){
	  fgets(buffer,LINE_SIZE,infile);
	  cptr=strstr(buffer,"<k_point");
	  if (!cptr) error_exit("Unexpected entry in kpts");
	  cptr=strstr(cptr,"weight=");
	  if (!cptr) error_exit("No weight entry in kpts");
	  cptr+=strlen("weight=");
	  cptr++;
	  j=sscanf(cptr,"%lf",&k_in->kpts[i].wt);
	  if (j!=1) error_exit("Error parsing kpt weight");
	  while((*cptr)&&(*cptr!='>')) cptr++;
	  if (!*cptr) error_exit("Kpt misparse 1");
	  cptr++;
	  j=sscanf(cptr,"%lf %lf %lf",k_in->kpts[i].abs,
		   k_in->kpts[i].abs+1,k_in->kpts[i].abs+2);
	  if (j!=3) error_exit("Kpt misparse 2");
	}
      }
      else{
	cptr=strstr(buffer,"<monkhorst_pack");
	if (cptr){
	  k->mp=malloc(sizeof(struct mp_grid));
	  if (!k->mp) error_exit("Malloc error for struct mp_grid!");
	  for(i=0;i<3;i++) k->mp->disp[i]=0;

	  for(i=0;i<3;i++){
	    sprintf(key,"nk%1d=\"",i+1);
	    cptr=strstr(buffer,key);
	    if(!cptr){
	      fprintf(stderr,"Failed to find %s\n",key);
	      exit(1);
	    }
	    j=sscanf(cptr+5,"%d",k->mp->grid+i);
	    if (j!=1) error_exit("Error parsing MP grid");
	  }
	  for(i=0;i<3;i++){
	    sprintf(key," k%1d=\"",i+1);  /* NB leading space */
	    cptr=strstr(buffer,key);
	    if(!cptr){
	      fprintf(stderr,"Failed to find %s\n",key);
	      exit(1);
	    }
	    j=sscanf(cptr+5,"%d",&tmp);
	    if (j!=1) error_exit("Error parsing MP disp");
	    /* QE's convention is that all grids include origin,
             * ours that only odd grids include the origin.
             * So shift if off[i]==0 and even grid, or if
             * off[i]==1 and odd grid.
             */
	    if (((tmp==0)&&((k->mp->grid[i]&1)==0))||
                ((tmp==1)&&((k->mp->grid[i]&1)==1)))
              k->mp->disp[i]=0.5/k->mp->grid[i];
            else
              k->mp->disp[i]=0;
	  }
	}
      }  /* end MP */
    }
    else if(strstr(buffer,"<nks>")){
      cptr=strstr(buffer,"<nks>");
      sscanf(cptr+5,"%d",&k->n);
      k->kpts=malloc(k->n*sizeof(struct atom));
      if (!k->kpts) error_exit("Malloc error for kpts");
      nkpt=0;
    }
    else if (strstr(buffer,"<fermi_energy>")){
      e->e_fermi=malloc(sizeof(double));
      if (!e->e_fermi) error_exit("Malloc error");
      cptr=strstr(buffer,"<fermi_energy>");
      if (sscanf(cptr+14,"%lf",e->e_fermi)==0){
        free(e->e_fermi);
        e->e_fermi=NULL;
      }
      else
        *e->e_fermi*=H_eV;
    }
    else if(strstr(buffer,"<ks_energies>")){
      if (k->n==0) error_exit("ks_energies before nks");
      while(!strstr(buffer,"</ks_energies>")){
	fgets(buffer,LINE_SIZE,infile);
	if (strstr(buffer,"<k_point")){
	  if (nkpt>=k->n) error_exit("Too many kpoints found!");
	  cptr=strstr(buffer,"weight=");
	  if (!cptr) error_exit("No weight entry in kpts");
	  cptr+=strlen("weight=");
	  cptr++;
	  j=sscanf(cptr,"%lf",&k->kpts[nkpt].wt);
	  if (j!=1) error_exit("Error parsing kpt weight");
	  while((*cptr)&&(*cptr!='>')) cptr++;
	  if (!*cptr) error_exit("Kpt misparse 1");
	  cptr++;
	  j=sscanf(cptr,"%lf %lf %lf",k->kpts[nkpt].abs,
		   k->kpts[nkpt].abs+1,k->kpts[nkpt].abs+2);
	  if (j!=3) error_exit("Kpt misparse 2");
	}
	if (strstr(buffer,"<occupations ")){
	  cptr=strstr(buffer,"size=");
	  if (!cptr) error_exit("No size entry in occupations");
	  cptr+=strlen("size=");
	  cptr++;
	  tmp=0;
	  j=sscanf(cptr,"%d",&tmp);
	  if (!j) fprintf(stderr,"Error reading size in occupations\n"); 
	  if (nbands==0) nbands=tmp/e->nspins;
	  if (tmp==nbands*e->nspins){
	    if (!e->occ)
              e->occ=malloc(nbands*e->nspins*k->n*sizeof(double));
	    if (!e->occ) error_exit("Malloc error for occupations");
	    j=nbands*e->nspins*nkpt;
	    while((*cptr)&&(*cptr!='>')) cptr++;
	    i=0;
	    while(i<nbands*e->nspins){
	      while(sscanf(cptr,"%lf%n",e->occ+j+i,&tmp)==1){
		cptr+=tmp;
		i++;
	      }
	      if (i<nbands*e->nspins){
		fgets(buffer,LINE_SIZE,infile);
		cptr=buffer;
	      }
	    }
	  }
	  else
	    fprintf(stderr,"Unexpected size in occupations. Ignoring\n");
	}
	if (strstr(buffer,"<eigenvalues ")){
	  cptr=strstr(buffer,"size=");
	  if (!cptr) error_exit("No size entry in eigenvalues");
	  cptr+=strlen("size=");
	  cptr++;
	  tmp=0;
	  j=sscanf(cptr,"%d",&tmp);
	  if (!j) fprintf(stderr,"Error reading size in eigenvalues\n"); 
	  if (nbands==0) nbands=tmp/e->nspins;
	  if (tmp==nbands*e->nspins){
	    if (!e->eval) e->eval=malloc(nbands*e->nspins*k->n*sizeof(double));
	    if (!e->eval) error_exit("Malloc error for evals");
	    j=nbands*e->nspins*nkpt;
	    while((*cptr)&&(*cptr!='>')) cptr++;
	    i=0;
	    while(i<nbands*e->nspins){
	      while(sscanf(cptr,"%lf%n",e->eval+j+i,&tmp)==1){
		cptr+=tmp;
		i++;
	      }
	      if (i<nbands*e->nspins){
		fgets(buffer,LINE_SIZE,infile);
		cptr=buffer;
	      }
	    }
	  }
	  else
	    fprintf(stderr,"Unexpected size in eigenvalues. Ignoring\n");
	}
      }
      nkpt++;

    }
    else if(strstr(buffer,"<lsda>true</lsda>")){
      e->nspins=e->nbspins=2;
    }
    else if(strstr(buffer,"<noncolin>true</noncolin>")){
      e->nspinors=2;
    }
    else if(strstr(buffer,"<forces ")){
      for(i=0;i<m->n;i++){
        fgets(buffer,LINE_SIZE,infile);
        j=sscanf(buffer,"%lf %lf %lf",m->atoms[i].force,m->atoms[i].force+1,
                 m->atoms[i].force+2);
        if (j!=3) fprintf(stderr,"Warning: error parsing forces\n");
        else
          for(j=0;j<3;j++)  /* Units assumed to be Ha/alat */
            m->atoms[i].force[j]*=H_eV/(alat*BOHR);
      }
      m->forces=1;
    }
    else if(strstr(buffer,"<atomic_species ")){
      cptr=strstr(buffer,"<atomic_species ")+strlen("<atomic_species ");
      cptr2=strstr(cptr,"ntyp=\"");
      if (!cptr2) error_exit("Error parsing atomic_species");
      cptr2+=6;
      i=sscanf(cptr2,"%d",&nspec);
      if (i!=1) error_exit("Error parsing atomic_species nspec\n%s");
      if (debug>2) fprintf(stderr,"nspec=%d\n",nspec);
      species=malloc(nspec*sizeof(struct sp));
      for(i=0;i<nspec;i++){
        species[i].mass=0;
        species[i].mag=0;
        species[i].name=NULL;
        species[i].p_pot=NULL;
        while(fgets(buffer,LINE_SIZE,infile))
          if (strstr(buffer,"<species ")) break;
        cptr=strstr(buffer,"name=\"");
        if (!cptr) error_exit("Error parsing species");
        cptr+=6;
        cptr2=cptr;
        while(*cptr2&&(*cptr2!='"')) cptr2++;
        if (*cptr2!='"') error_exit("Error parsing species name");
        *cptr2=0;
        species[i].name=malloc((cptr2-cptr)+1);
        strcpy(species[i].name,cptr);
        while(fgets(buffer,LINE_SIZE,infile)){
          if (strstr(buffer,"</species>")) break;
          if (strstr(buffer,"<mass>")){
            cptr=strstr(buffer,"<mass>")+6;
            j=sscanf(cptr,"%lf",&species[i].mass);
            if (j!=1) fprintf(stderr,"Warning, error parsing species mass\n");
          }
          else if (strstr(buffer,"<pseudo_file>")){
            cptr=strstr(buffer,"<pseudo_file>")+13;
            cptr2=strstr(cptr,"</pseudo_file>");
            if ((!cptr)||(!cptr2))
              fprintf(stderr,"Warning, error parsing pseudo_file\n");
            else{
              species[i].p_pot=malloc((cptr2-cptr)+1);
              strncpy(species[i].p_pot,cptr,cptr2-cptr);
              species[i].p_pot[cptr2-cptr]=0;
            }
	  }
	  else if (strstr(buffer,"<starting_magnetization>")){
	    cptr=strstr(buffer,"<starting_magnetization>")+
	      strlen("<starting_magnetization>");
	    j=sscanf(cptr,"%lf",&species[i].mag);
            if (j!=1)
	      fprintf(stderr,
		      "Warning, error parsing species magnetic moment\n");
          }
        }
      }
    }
    else if(strstr(buffer,"<symmetries>")){
      while(fgets(buffer,LINE_SIZE,infile))
        if(strstr(buffer,"<nsym>")) break;
      cptr=strstr(buffer,"<nsym>");
      if (cptr){
        cptr+=6;
        i=sscanf(cptr,"%d",&nsym);
        if ((i!=1)||(nsym<1)) fprintf(stderr,"Error parsing nsym\n");
        else{
          s->n=nsym;
          if (debug>2) fprintf(stderr,"Reading %d sym ops\n",nsym);
          s->ops=malloc(s->n*sizeof(struct sym_op));
          symrel=malloc(s->n*9*sizeof(double));
          symtr=malloc(s->n*3*sizeof(double));
          if ((!s->ops)||(!symrel)||(!symtr))
            error_exit("Malloc error for symmetry ops");
          for(i=0;i<s->n;i++){
            while(fgets(buffer,LINE_SIZE,infile))
              if(strstr(buffer,"<symmetry>")) break;
            while(fgets(buffer,LINE_SIZE,infile))
              if(strstr(buffer,"<info ")) break;
            if(!strstr(buffer,"crystal_symmetry")){i--;continue;}
            while(fgets(buffer,LINE_SIZE,infile))
              if(strstr(buffer,"<rotation ")) break;
            if (!fgets(buffer,LINE_SIZE,infile))
              error_exit("Unexpected EOF reading symmetries");
            cptr=buffer;
            j=0;
            do{
              while (sscanf(cptr,"%lf%n",symrel+9*i+j,&jj)>0){
                cptr+=jj;
                j++;
              }
              if (!fgets(buffer,LINE_SIZE,infile))
                error_exit("Error reading symmetry operations");
              cptr=buffer;
            } while (j<9);
            while(fgets(buffer,LINE_SIZE,infile))
              if(strstr(buffer,"<fractional_translation>")) break;
            cptr=strstr(buffer,"<fractional_translation>");
            if (!cptr) error_exit("Error reading symmetry operations");
            cptr+=strlen("<fractional_translation>");
            s->ops[i].tr=malloc(3*sizeof(double));
            if (!s->ops[i].tr)
              error_exit("Malloc error for symmetry translation");
            j=0;
            do{
              while (sscanf(cptr,"%lf%n",symtr+3*i+j,&jj)>0){
                cptr+=jj;
                j++;
              }
              if (!fgets(buffer,LINE_SIZE,infile))
                error_exit("Error reading symmetry operations");
              cptr=buffer;
            } while (j<3);
          } /* end loop of sym ops */
        }
      }
    } /* end symmetry */
    else if(strstr(buffer,"<total_energy>")){
      while(fgets(buffer,LINE_SIZE,infile)){
        if(strstr(buffer,"</total_energy>")) break;
        cptr=strstr(buffer,"<etot>");
        if(cptr){
          cptr+=6;
          e->energy=malloc(sizeof(double));
          i=sscanf(cptr,"%lf",e->energy);
          if (i!=1){
            fprintf(stderr,"Warning: error parsing energy\n");
            free(e->energy);
            e->energy=NULL;
          }
          else *e->energy*=H_eV; /* Its units were Ha */
        }
      }
    }
    else if(strstr(buffer,"<basis_set>")){
      while(fgets(buffer,LINE_SIZE,infile)){
        if(strstr(buffer,"</basis_set>")) break;
        cptr=strstr(buffer,"<ecutwfc>");
        if(cptr){
          cptr+=9;
          i=sscanf(cptr,"%lf",&e->cut_off);
          if (i!=1){
            fprintf(stderr,"Warning: error parsing cut-off\n");
            e->cut_off=0;
          }
          else e->cut_off*=H_eV; /* Its units were Ha */
        }
	if ((cptr=strstr(buffer,"<fft_grid "))){
	  cptr+=10;
	  i=sscanf(cptr,"nr1=\"%d\" nr2=\"%d\" nr3=\"%d\"",
		   fft,fft+1,fft+2);
	  if (i!=3){
	    fprintf(stderr,"Warning: error parsing FFT grid size\n");
	    fft[0]=fft[1]=fft[2]=0;
	  }
	}
      }
    }
  }

  if(!strstr(buffer,"</output>"))
    fprintf(stderr,"Warning, unexpected end to XML file\n");

  if (!c->basis) error_exit("No basis found in qe_xml_read");

  if(nspec){
    for(i=0;i<nspec;i++){
      pot=qe_pspot_file(species[i].p_pot,dict_get(m->dict,"QE_pseudo_dir"),cwd);
      if (pot){
        species[i].chg=qe_pot_charge(pot);
        fclose(pot);
      }
      if (debug>2) fprintf(stderr,"%d %s %lf %s z=%lf\n",i,species[i].name,
			   species[i].mass,species[i].p_pot,species[i].chg);
    }
  }
  
  /* Fix all units etc */

  if (nkpt!=k->n) {
    fprintf(stderr,"Warning: expected %d kpoints, found %d, deleting\n",
	    k->n,nkpt);
    k->n=0;
  }
  
  if (alat){
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
        c->basis[i][j]/=BOHR;
    real2rec(c);
    for(i=0;i<k_in->n;i++){
      if (debug>2)
	fprintf(stderr,"Unscaled kpt (%f,%f,%f)\n",k_in->kpts[i].abs[0],
		k_in->kpts[i].abs[1],k_in->kpts[i].abs[2]);
      for(j=0;j<3;j++) k_in->kpts[i].abs[j]*=1/alat;
      kcart2frac(k_in->kpts[i].abs,k_in->kpts[i].frac,c->recip);
      for(j=0;j<3;j++) k_in->kpts[i].abs[j]=0;
    }
    sum=0;
    for(i=0;i<k_in->n;i++) sum+=k_in->kpts[i].wt;
    for(i=0;i<k_in->n;i++) k_in->kpts[i].wt/=sum;
    
    for(i=0;i<k->n;i++){
      if (debug>2)
	fprintf(stderr,"Unscaled kpt (%f,%f,%f)\n",k->kpts[i].abs[0],
		k->kpts[i].abs[1],k->kpts[i].abs[2]);
      for(j=0;j<3;j++) k->kpts[i].abs[j]*=1/alat;
      kcart2frac(k->kpts[i].abs,k->kpts[i].frac,c->recip);
      for(j=0;j<3;j++) k->kpts[i].abs[j]=0;
    }
    sum=0;
    for(i=0;i<k->n;i++) sum+=k->kpts[i].wt;
    for(i=0;i<k->n;i++) k->kpts[i].wt/=sum;
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
        c->basis[i][j]*=BOHR;
  }
  else{
    fprintf(stderr,"alat not set -- kpoints cannot be retained\n");
    k->n=0;
    k_in->n=0;
  }

  if (k_in->n) dict_add(m->dict,"QE_k_in",(void*)k_in);
  

  real2rec(c);

  
  addfrac(m->atoms,m->n,c->recip);

  if (nspec){
    for(i=0;i<m->n;i++){
      cptr=m->atoms[i].label;
      if (!cptr) cptr=atno2sym(m->atoms[i].atno);
      if (cptr){
	for(j=0;j<nspec;j++){
	  if (!strcasecmp(cptr,species[j].name)){
	    m->atoms[i].chg=species[j].chg;
	    m->atoms[i].spin=species[j].mag;
	    break;
	  }
	}
      }
    }
  }

  if (nspec){
    cptr=NULL;
    cptr2=cptr;
    j=0;
    tmp=0;
    for(i=0;i<nspec;i++){
      j+=snprintf(NULL,0," %3s  %lf  %s\n",species[i].name,
		  species[i].mass,species[i].p_pot);
      cptr=realloc(cptr,j+1);
      sprintf(cptr+tmp," %3s  %lf  %s\n",species[i].name,
	      species[i].mass,species[i].p_pot);
      tmp=j;
    }
    dict_add(m->dict,"QE_atomic_species",cptr);
  }


  /* Normalise etol (will be zero if not read) */

  e->etol/=m->n;
  
  if (s->n){
    for(i=0;i<s->n;i++){
      for(j=0;j<3;j++)
        for(jj=0;jj<3;jj++)
          mat[j][jj]=symrel[9*i+3*j+jj];

      mat_f2a(mat,s->ops[i].mat,c->basis,c->recip);

      for(j=0;j<3;j++){
        s->ops[i].tr[j]=0;
        for(jj=0;jj<3;jj++)
          s->ops[i].tr[j]+=symtr[3*i+jj]*c->basis[jj][j];
      }
    }
  }

  if (e->eval) /* We store these in eV */
    for(i=0;i<nbands*e->nspins*k->n;i++)
      e->eval[i]*=H_eV;
  e->nbands=nbands;
  
  if ((flags&CHDEN)||(flags&SPINDEN)){
    i=strlen(filename);
    cptr2=filename+i;
    while((cptr2>filename)&&(*cptr2!='/')) cptr2--;
    if (*cptr2=='/') i=cptr2-filename+1;
    den_name=malloc(1);
    *den_name=0;
    if (cptr2!=filename){
      den_name=realloc(den_name,strlen(filename)+1);
      strncpy(den_name,filename,i);
      den_name[i]=0;
    }
    den_name=realloc(den_name,strlen(den_name)+strlen("charge-density.dat")+1);
    strcat(den_name,"charge-density.dat");
    ch=fopen(den_name,"r");
    if (!ch){
      den_name2=NULL;
      if (!strcmp(filename+strlen(filename)-4,".xml")){
	den_name2=malloc(strlen(filename)-3);
	strncpy(den_name2,filename,strlen(filename)-4);
	den_name2[strlen(filename)-4]=0;
	den_name2=realloc(den_name2,strlen(den_name2)+
			  strlen(".save/charge-density.dat")+1);
	strcat(den_name2,".save/charge-density.dat");
	ch=fopen(den_name2,"r");
      }
      if (!ch){
	fprintf(stderr,
		"Warning: density requested, tried to open %s ",den_name);
	if (den_name2) fprintf(stderr,"and %s ",den_name2);
	fprintf(stderr,"but failed\n");
      }
    }
    if (ch) qe_rho_read(ch, c, m, k, s, g, e, i_grid);
  }

  if (flags&BANDREAD)
    qe_psi_read(cwd,(char*)dict_get(m->dict,"QE_prefix"),c, m, k, s, g,
		e, fft, i_grid);
  
  
}

/* See if we can find a pspot file */
FILE *qe_pspot_file(char *name, char *prefix, char *cwd){
  FILE *p;
  char *path,*ptr;

  if (debug>3)
    fprintf(stderr,"qe_pspot_file called with prefix=%s cwd=%s\n",prefix,cwd);
  
  if (prefix){
    if ((prefix[0]!='/')&&(cwd)){
      path=malloc(strlen(cwd)+strlen(prefix)+strlen(name)+3);
      if (!path) error_exit("malloc error for filename");
      strcpy(path,cwd);
      strcat(path,"/");
      strcat(path,prefix);
      strcat(path,"/");
      strcat(path,name);
      if (debug>2) fprintf(stderr,"Trying to open %s\n",path);
      p=fopen(path,"r");
      if (p){
        if (debug>1) fprintf(stderr,"Found pseudopot file %s\n",path);
        free(path);
        return p;
      }
      free(path);
    }
    
    path=malloc(strlen(prefix)+strlen(name)+2);
    if (!path) error_exit("malloc error for filename");
    strcpy(path,prefix);
    strcat(path,"/");
    strcat(path,name);
    if (debug>2) fprintf(stderr,"Trying to open %s\n",path);
    p=fopen(path,"r");

    if (p){
      if (debug>1) fprintf(stderr,"Found pseudopot file %s\n",path);
      free(path);
      return p;
    }
    free(path);
  }

  ptr=getenv("PSEUDO_DIR");
  if (ptr){
    if ((ptr[0]!='/')&&(cwd)){
      path=malloc(strlen(cwd)+strlen(ptr)+strlen(name)+3);
      if (!path) error_exit("malloc error for filename");
      strcpy(path,cwd);
      strcat(path,"/");
      strcat(path,ptr);
      strcat(path,"/");
      strcat(path,name);
      if (debug>2) fprintf(stderr,"Trying to open %s\n",path);
      p=fopen(path,"r");
      if (p){
        if (debug>1) fprintf(stderr,"Found pseudopot file %s\n",path);
        free(path);
        return p;
      }
      free(path);
    }
    path=malloc(strlen(ptr)+strlen(name)+2);
    if (!path) error_exit("malloc error for filename");
    strcpy(path,ptr);
    strcat(path,"/");
    strcat(path,name);

    if (debug>2) fprintf(stderr,"Trying to open %s\n",path);
    p=fopen(path,"r");
    if (p){
      if (debug>1) fprintf(stderr,"Found pseudopot file %s\n",path);
      free(path);
      return p;
    }
    free(path);
  }

  if (cwd){
    path=malloc(strlen(cwd)+strlen(name)+2);
    if (!path) error_exit("malloc error for filename");
    strcpy(path,cwd);
    strcat(path,"/");
    strcat(path,name);
    if (debug>2) fprintf(stderr,"Trying to open %s\n",path);
    p=fopen(path,"r");
    if (p){
      if (debug>1) fprintf(stderr,"Found pseudopot file %s\n",path);
      free(path);
      return p;
    }
    free(path);
  }
  
  if (debug>2) fprintf(stderr,"Trying to open %s\n",name);
  p=fopen(name,"r");
  if (debug&&(!p))
    fprintf(stderr,"Failed to find pseudopot file for %s\n",name);
  if ((p)&&(debug>1)) fprintf(stderr,"Found pseudopot file %s\n",name);
    
  return p;

}

double qe_pot_charge(FILE *infile){
  char buffer[LINE_SIZE+1],*ptr;
  double chg;
  int i,hit;

  hit=0;
  while(fgets(buffer,LINE_SIZE,infile)){
    ptr=strstr(buffer,"<PP_HEADER");
    if ((ptr)&&(*(ptr+10)=='>')) {hit=1;break;}
    if ((ptr)&&((*(ptr+10)==' ')||(*(ptr+10)=='\n'))) {hit=2;break;}
  }
  
  if (hit==2){  /* Version 2 UPF file */
    ptr=strstr(ptr,"z_valence=\"");
    if (!ptr) { /* Has <PP_HEADER been split over multiple lines? */
      while (!strstr(buffer,"/>")){
        if (!fgets(buffer,LINE_SIZE,infile)) break;
        ptr=strstr(buffer,"z_valence=\"");
        if (ptr) break;
      }
      if (!ptr){
        fprintf(stderr,"Warning: ionic pseudo charge not found\n");
        return 0;
      }
    }
    ptr+=strlen("z_valence=\"");
    i=sscanf(ptr,"%lf",&chg);
    if (i==1)
      return chg;
    else
      fprintf(stderr,"Warning: ionic pseudo charge not parsed\n");
  } /* Pre version 2 UPF file */
  else if (hit==1){
    /* Skip five lines */
    for(i=0;i<5;i++) fgets(buffer,LINE_SIZE,infile);
    if (!fgets(buffer,LINE_SIZE,infile)){
      fprintf(stderr,"Warning: unexpected end to pseudopot file\n");
      return 0;
    }
    i=sscanf(buffer,"%lf",&chg);
    if (i==1)
      return chg;
    else
      fprintf(stderr,"Warning: ionic pseudo charge not parsed\n");
  }
  else
    fprintf(stderr,"Warning: PP_HEADER not found in pseudopot file\n");  
  
  return 0;
}

static void read_step(FILE *infile, char *buffer,struct time_series *ts){
  int i,j,n;
  double alat;
  char *ptr,*cptr;
  struct contents *m;
  struct atom atmp;
  
  ptr=strstr(buffer,"<step n_step=");
  if (!ptr) return;

  ptr+=strlen("<step n_step=")+1;
  i=sscanf(ptr,"%d",&n);
  if (n!=ts->nsteps+1){
    fprintf(stderr,"Ignoring step %d as %d expected\n",n,ts->nsteps+1);
    return;
  }

  ts->nsteps=n;

  while(fgets(buffer,LINE_SIZE,infile)){
    if (strstr(buffer,"</step>")) break;
    if(strstr(buffer,"<cell>")){
      ts->cells=realloc(ts->cells,(ts->nc+1)*sizeof(struct unit_cell));
      if (!ts->cells) error_exit("Realloc error for ts->cells");
      ts->cells[ts->nc].basis=NULL;
      ts->cells[ts->nc].stress=NULL;
      read_basis(infile,ts->cells+ts->nc);
      ts->nc++;
    }
    else if(strstr(buffer,"<atomic_structure ")){
      cptr=strstr(buffer,"<atomic_structure>");
      cptr=strstr(buffer,"alat=");
      cptr+=5;
      if (*cptr=='"') cptr++;
      sscanf(cptr,"%lf",&alat);
    }
    else if(strstr(buffer,"<atomic_positions>")){
      ts->m=realloc(ts->m,(ts->nm+1)*sizeof(struct contents));
      if (!ts->m) error_exit("Realloc error for ts->m");
      ts->m[ts->nm].atoms=NULL;
      ts->m[ts->nm].title=NULL;
      ts->m[ts->nm].comment=NULL;
      ts->m[ts->nm].block_species=NULL;
      ts->m[ts->nm].species_misc=NULL;
      ts->m[ts->nm].dict=NULL;
      ts->m[ts->nm].forces=0;
      ts->m[ts->nm].n=0;
      read_atoms(infile,ts->m+ts->nm);
      ts->nm++;
    }
    else if(strstr(buffer,"<forces ")){
      if (ts->nm<=1) error_exit("Impossible: forces before atoms");
      if (debug) fprintf(stderr,"%d parsing forces with alat=%f\n",
                         n,alat);
      m=ts->m+(ts->nm-2);
      for(i=0;i<m->n;i++){
        fgets(buffer,LINE_SIZE,infile);
        j=sscanf(buffer,"%lf %lf %lf",m->atoms[i].force,m->atoms[i].force+1,
                 m->atoms[i].force+2);
        if (j!=3) fprintf(stderr,"Warning: error parsing forces\n");
        else
          if (alat){
            for(j=0;j<3;j++){  /* Units assumed to be Ha/alat */
              m->atoms[i].force[j]*=H_eV/(alat*BOHR);
            }
            /* The following seems to be necessary... */
            for(j=0;j<3;j++)
              atmp.abs[j]=m->atoms[i].force[j];
            addfrac(&atmp,1,ts->cells[ts->nc-2].recip);
            addabs(&atmp,1,ts->cells[ts->nc-1].basis);
            for(j=0;j<3;j++)
              m->atoms[i].force[j]=atmp.abs[j];
          }
          else
            fprintf(stderr,"Warning: forces without alat\n");
      }
      m->forces=1;
    }
    else if (strstr(buffer,"<etot>")){
      ts->energies=realloc(ts->energies,(ts->nen+1)*sizeof(double));
      cptr=strstr(buffer,"<etot>")+6;
      i=sscanf(cptr,"%lf",ts->energies+ts->nen);
      if (i!=1){
        fprintf(stderr,"Warning: error parsing energy\n");
        ts->energies[ts->nen]=-999;
      }
      else ts->energies[ts->nen]*=H_eV; /* Its units were Ha */
      ts->nen++;
    }
  }

}

static void read_basis(FILE *infile, struct unit_cell *c){
  int i,j;
  char *cptr,key[6];
  char buffer[LINE_SIZE+1];
  
  c->basis=malloc(72);
  if (!c->basis) error_exit("Malloc error for basis");
  for(i=0;i<3;i++){
    sprintf(key,"<a%1d>",i+1);
    fgets(buffer,LINE_SIZE,infile);
    cptr=strstr(buffer,key);
    if(!cptr){
      fprintf(stderr,"Failed to find %s\n",key);
      exit(1);
    }
    j=sscanf(cptr+4,"%lf %lf %lf",c->basis[i],c->basis[i]+1,c->basis[i]+2);
    if (j!=3){
      fprintf(stderr,"Failed to scan %s\n",key);
      exit(1);
    }
    for(j=0;j<3;j++) c->basis[i][j]*=BOHR;
  }
  real2rec(c);
}


static void read_atoms(FILE *infile, struct contents *m){
  int i;
  char *cptr,*cptr2,*cptr3;
  char buffer[LINE_SIZE+1];
  
  while(fgets(buffer,LINE_SIZE,infile)){
    if (strstr(buffer,"</atomic_positions>")) break;
    cptr=strstr(buffer,"name=");
    if (!cptr){
      fprintf(stderr,"Warning: ignoring %s\n",buffer);
      continue;
    }
    m->n++;
    m->atoms=realloc(m->atoms,m->n*sizeof(struct atom));
    if (!m->atoms) error_exit("realloc error");
    init_atoms(m->atoms+(m->n-1),1);
    cptr+=5;
    while((*cptr)&&(*cptr!='"')) cptr++;
    if (!*cptr) error_exit("Atom misparse 1");
    cptr++;
    cptr2=cptr;
    while(isalpha(*cptr2)) cptr2++;
    if (!*cptr2) error_exit("Atom misparse 2");
    if (*cptr!='"'){ /* We have a label which is not a simple atomic sym */
      cptr3=cptr2;
      while(*cptr3!='"') cptr3++; /* cptr3 points one past end */
      if (!*cptr3) error_exit("Atom misparse 2.5");
      m->atoms[m->n-1].label=malloc((cptr3-cptr)+1);
      strncpy(m->atoms[m->n-1].label,cptr,cptr3-cptr);
      m->atoms[m->n-1].label[cptr3-cptr]=0;
    }
    *cptr2=0;
    m->atoms[m->n-1].atno=atsym2no(cptr);
    cptr=cptr2+1;
    while((*cptr)&&(*cptr!='>')) cptr++;
    if (!*cptr) error_exit("Atom misparse 3");
    cptr++;
    i=sscanf(cptr,"%lf %lf %lf",m->atoms[m->n-1].abs,
		 m->atoms[m->n-1].abs+1,m->atoms[m->n-1].abs+2);
    if (i!=3) error_exit("Atom misparse 4");
    for(i=0;i<3;i++) m->atoms[m->n-1].abs[i]*=BOHR;
  }

}
