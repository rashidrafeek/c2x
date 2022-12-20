#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<math.h>
#include<string.h>
#include "c2xsf.h" /* for debug flag */

/* One function for helping to tokenise input files, used
   by the abinit reader and others */

/* Act like strncasecmp but ensure that token stops after
 * match, i.e. "sin" matches "sin" and "sin .*", but not "since"
 *
 * Accept whitespace or = as terminators to token
 *
 * Also, if there is a match, increment s1 to first char beyond
 * match.
 */

int tokenmatch(char **s1, const char *s2){
  int tmp,len;

  len=strlen(s2);
  tmp=strncasecmp(*s1,s2,len);
  if (tmp!=0) return tmp;
  if (((*s1)[len]==0)||(isspace((*s1)[len]))||((*s1)[len]=='=')){
    *s1+=len;
    return 0;
  }
  return 1;
}

/* And a function to act like sscanf("%ms%n",buffer,str,n) for
 *  MacOS X does not support %ms in sscanf (which is not in C99, but
 *  is in POSIX.1-2008)
 */

int sscanfmsn(char *buffer, char **str, int *n){
  char *ptr1,*ptr2,*s;
  int i;

  ptr1=buffer;
  if (n) *n=0;
  while ((*ptr1)&&(isspace(*ptr1))) ptr1++;

  ptr2=ptr1;
  while ((*ptr2)&&(!isspace(*ptr2))) ptr2++;

  if (ptr1==ptr2) return 0;

  s=malloc(ptr2-ptr1+1);
  if (!s) error_exit("malloc error in sscanfms");
  *str=s;
  
  i=0;
  while(ptr1<ptr2) s[i++]=*(ptr1++);
  s[i]=0;
  if (n) *n=ptr2-buffer;
  return 1;
}


/* The rest of this file contains routines for parsing general
 * arithmetic expressions, as required by the Castep cell file
 * format since verison 18. C2x also uses this to parse certain
 * parts of Abinit files, although Abinit has a much more restricted
 * syntax.
 */


int parse(char* in);
int evaluate_st(double *result);
static struct op {char *op; int prec;} *op_stack;
static struct a {void *data; int type;} *stack;
/* pointers point to free element at end of stack */
static int op_ptr,st_ptr;

/* Return zero on failure, one on success, like sscan */
int ascan(char *in, double *result){
  int i;
  if (parse(in)) return 0;
  i=evaluate_st(result);
  if (stack) {free(stack); stack=NULL;}
  if (op_stack) {free (op_stack); op_stack=NULL;}
  return (!i);
}

/* act like sscanf(buff,"%lf%n",x,n)
 * but evaluate arithmetic expressions
 * and ignore n if NULL
 */
int single_scan(char *buff, double *x, int *n){
  char s[81];
  int i,dummy;

  if (!n) n=&dummy;
  
  i=sscanf(buff,"%lf%n",x,n);
  if ((i>=1)&&((*(buff+(*n))==0)||(*(buff+(*n))==' '))) return 1;
  
  i=sscanf(buff,"%80s%n",s,n);
  if (i==0) return 0;
  if (ascan(s,x)==0) return 0;
  return 1;
}

/* as single_scan, but have repeat count and read a vector */
int multi_scan(char *buff, double *x, int rep, int *n){
  int i,n2,dummy;
  n2=0;

  if (!n) n=&dummy;
  *n=0;
  
  for(i=0;i<rep;i++){
    if (single_scan(buff+(*n),x+i,&n2)==0) return i;
    *n+=n2;
  }

  return rep;
}


static void pop_op(void){
  stack=realloc(stack,(st_ptr+1)*sizeof(struct a));
  stack[st_ptr].type='O';
  stack[st_ptr].data=op_stack[op_ptr-1].op;
  st_ptr++;
  op_ptr--;
}

static void push_op(char* op, int pr){
  op_stack=realloc(op_stack,(op_ptr+1)*sizeof(struct op));
  op_stack[op_ptr].prec=pr;
  op_stack[op_ptr].op=op;
  op_ptr++;
}

/* Take infix expression, parse according to shunting yard algorithm */
/* On success returns zero, on failure returns one */
int parse(char* in){
  char *p1,*p2,*op,*c;
  int i,unitary,l_assoc,is_num,pr,implicit_mul;
  double *x;
  
  p1=in;
  op_stack=NULL;
  stack=NULL;
  op_ptr=0;
  st_ptr=0;
  
  unitary=1;
  is_num=0;

  if (debug>2) fprintf(stderr,"Parsing '%s'\n",in);

  while (isspace(*p1)) p1++; /* consume whitespace */
  if ((*p1=='!')||(*p1=='#')) return 1; /* all we have is a comment */
 
  while(*p1){ /* Find end of token */
    implicit_mul=0;
    if (isdigit(*p1)||(*p1=='.')||(unitary&&((*p1=='+')||(*p1=='-')))){
      /* We have a number */
      /* (If this is a unitary + or -, we detect this at the next stage) */
      p2=p1+1;
      while (isdigit(*p2)||(*p2=='.')) p2++;
      /* we might have an exponential format */
      if (((*p2=='e')||(*p2=='E')||(*p2=='d')||(*p2=='D'))&&
	  ((*(p2+1)=='+')||(*(p2+1)=='-')||(isdigit(*(p2+1))))){
	*p2='e'; /* sscanf, used later, does not recognise d */
	p2++;
	if ((*p2=='+')||(*p2=='-')) p2++;
	while (isdigit(*p2)) p2++;
      }
      
      if (debug>4) fprintf(stderr,"Token: %.*s\n",(int)(p2-p1),p1);
      unitary=0;
      is_num=1;
    }
    else if ((*p1=='+')||(*p1=='-')||(*p1=='*')||(*p1=='/')||(*p1=='^')||
	(*p1=='(')||(*p1==')')){
      if ((*p1=='(')&&(unitary==0)) implicit_mul=1;
      if (*p1==')') unitary=0;
      else unitary=1;
      is_num=0;
      if (debug>4) fprintf(stderr,"Token: %.1s\n",p1);
      p2=p1+1;
    }
    else if (isalpha(*p1)){
      p2=p1+1;
      while (isalpha(*p2)) p2++;
      if (debug>4) fprintf(stderr,"Token: %.*s\n",(int)(p2-p1),p1);
      is_num=0;
      if (unitary==0) implicit_mul=1;
      unitary=1;
    }
    else {
      fprintf(stderr,"Error tokenising %s\n",in);
      fprintf(stderr,"                 ");
      for (i=0;i<p1-in;i++) fprintf(stderr," ");
      fprintf(stderr,"^ error\n");
      return(1);
    }
    /* Now have token */
    /* Add token to stack following shunting algorithm */
    if (is_num){
      if(((p2-p1)==1)&&(*p1=='-')){ /* We have a unitary - */
        while((op_ptr)&&
              ((op_stack[op_ptr-1].prec>=50)
               &&(*op_stack[op_ptr-1].op)!='(')){
          pop_op();
        }
        /* ops must be free()able */
        c=malloc(4);
        if (!c) error_exit("malloc error for 4 bytes");
        strncpy(c,"neg",4);
        push_op(c,50);
        is_num=0;
        unitary=1;
      }
      else if (((p2-p1)==1)&&(*p1=='+')){ /* We have a unitary + */
        is_num=0;
        unitary=1;
      }
      else{
        x=malloc(sizeof(double));
        sscanf(p1,"%lf",x);
        stack=realloc(stack,(st_ptr+1)*sizeof(struct a));
        stack[st_ptr].type='N';
        stack[st_ptr].data=x;
        st_ptr++;
      }
    }
    else{ /* p2 points one char beyond the end of token */
      op=malloc(p2-p1+1);
      for(i=0;i<p2-p1;i++) op[i]=*(p1+i);
      op[p2-p1]=0;
      pr=100;
      l_assoc=0;
      if ((*op=='+')||(*op=='-')||
	  (*op=='*')||(*op=='/')||(*op=='/')) l_assoc=1;
      if ((*op=='+')||(*op=='-')) pr=10;
      if ((*op=='*')||(*op=='/')) pr=20;
      if (*op=='^') pr=100;
      if (*op=='(') pr=1000;
      if (*op==')'){
	while((op_ptr)&&(*(op_stack[op_ptr-1].op)!='(')){
	  pop_op();
	}
        if ((op_ptr==0)||(*(op_stack[op_ptr-1].op)!='(')){
	  fprintf(stderr,"Parenthesis problem parsing %s\n",in);
	  return(1);
	}
	free(op_stack[op_ptr-1].op); /* Pop the ( */
	op_ptr--;
	free(op); /* We don't store the ) */
      }
      else{
	if (implicit_mul){ /* Need to add an implict *, which we assume
			      has the same precidence as an 
			      explicit multiplication. This means that
                              1/2sqrt(3) is 0.5sqrt(3) which is probably
                              what would be meant in .cell files */
	  while((op_ptr)&&
		((op_stack[op_ptr-1].prec>=20)
		 &&(*op_stack[op_ptr-1].op)!='(')){
	    pop_op();
	  }
	  /* ops must be free()able */
	  c=malloc(2);
	  c[0]='*';
	  c[1]=0;
	  push_op(c,20);
	  implicit_mul=0;
	}
	if (*op=='('){
	  push_op(op,pr);
	}
	else if ((strcmp(op,"pi")==0)||(strcmp(op,"Ha")==0)||
		 (strcmp(op,"Ry")==0)||(strcmp(op,"B")==0)){
	  stack=realloc(stack,(st_ptr+1)*sizeof(struct a));
	  stack[st_ptr].type='N';
	  x=malloc(sizeof(double));
	  if (strcmp(op,"pi")==0)
	    *x=M_PI;
	  else if (strcmp(op,"Ha")==0)
	    *x=H_eV;
	  else if (strcmp(op,"Ry")==0)
	    *x=0.5*H_eV;
	  else if (strcmp(op,"B")==0)
	    *x=BOHR;
	  else
	    error_exit("Impossible in parser");
	  free(op);
	  stack[st_ptr].data=x;
	  st_ptr++;
	  is_num=1;
	  unitary=0;
	}
	else{
	  while((op_ptr)&&
		(((op_stack[op_ptr-1].prec>pr)||
		  ((l_assoc)&&(op_stack[op_ptr-1].prec==pr)))
		 &&(*op_stack[op_ptr-1].op)!='(')){
	    pop_op();
	  }
	  push_op(op,pr);
	}
      }
    }
    /* Token added to stack, increment our way through infix string */
    p1=p2;	
    while (isspace(*p1)) p1++; /* consume whitespace */
  }

  /* Move remaining operators from operator stack to output stack */
  while(op_ptr>0){
    if ((*op_stack[op_ptr-1].op)=='('){
      fprintf(stderr,"Parenthesis problem parsing %s\n",in);
      return(1);
    }
    pop_op();
  }

  if (debug>2){
    fprintf(stderr,"Parsed to: ");
    for(i=0;i<st_ptr;i++){
      if (stack[i].type=='O') fprintf(stderr,"%s ",(char*)stack[i].data);
      else fprintf(stderr,"%lf ",*((double*)stack[i].data));
    }
    fprintf(stderr,"\n");
  }
  return(0);
}

#define stack_underflow() {fprintf(stderr,"Stack underflow on evaluation\n"); \
    return(1);}

/* On success returns zero, on failure returns one */
int evaluate_st(double *result){
  double *num_st;
  char *op;
  int num_ptr,i,j;

  num_st=NULL;
  /* num_ptr points to first free item */
  num_ptr=0;

  for(i=0;i<st_ptr;i++){
    if (stack[i].type=='N'){
      num_st=realloc(num_st,(num_ptr+1)*sizeof(double));
      num_st[num_ptr]=*(double*)stack[i].data;
      num_ptr++;
    }
    else{
      op=stack[i].data;
      if (!strcmp(op,"+")){
	if (num_ptr<2) stack_underflow();
	num_st[num_ptr-2]+=num_st[num_ptr-1];
	num_ptr--;
      }
      else if (!strcmp(op,"-")){
	if (num_ptr<2) stack_underflow();
	num_st[num_ptr-2]-=num_st[num_ptr-1];
	num_ptr--;
      }
      else if (!strcmp(op,"*")){
	if (num_ptr<2) stack_underflow();
	num_st[num_ptr-2]*=num_st[num_ptr-1];
	num_ptr--;
      }
      else if (!strcmp(op,"/")){
	if (num_ptr<2) stack_underflow();
	num_st[num_ptr-2]/=num_st[num_ptr-1];
	num_ptr--;
      }
      else if (!strcmp(op,"^")){
	if (num_ptr<2) stack_underflow();
	num_st[num_ptr-2]=pow(num_st[num_ptr-2],num_st[num_ptr-1]);
	num_ptr--;
      }
      else if (!strcmp(op,"neg")){
	if (num_ptr<1) stack_underflow();
	num_st[num_ptr-1]=-num_st[num_ptr-1];
      }
      else if (!strcmp(op,"sqrt")){
	if (num_ptr<1) stack_underflow();
	num_st[num_ptr-1]=sqrt(num_st[num_ptr-1]);
      }
      else if (!strcmp(op,"cos")){
	if (num_ptr<1) stack_underflow();
	num_st[num_ptr-1]=cos(M_PI*num_st[num_ptr-1]/180);
      }
      else if (!strcmp(op,"sin")){
	if (num_ptr<1) stack_underflow();
	num_st[num_ptr-1]=sin(M_PI*num_st[num_ptr-1]/180);
      }
      else if (!strcmp(op,"tan")){
	if (num_ptr<1) stack_underflow();
	num_st[num_ptr-1]=tan(M_PI*num_st[num_ptr-1]/180);
      }
      else {
	fprintf(stderr,"Unknown operator %s\n",(char*)stack[i].data);
	return(1);
      }
    }
    free(stack[i].data);
    if (debug>4){
      fprintf(stderr,"Iteration %d: ",i);
      for(j=0;j<num_ptr;j++) fprintf(stderr,"%lf ",num_st[j]);
      fprintf(stderr,"\n");
    }
  }
  if (num_ptr!=1){
    fprintf(stderr,"Stack imbalance at end of evaluation\n");
    return(1);
  }
  *result=num_st[0];
  if (debug>2) fprintf(stderr,"Result returned: %lf\n",*result);
  free(num_st);
  return (0);
}
