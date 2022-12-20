/* Command-line routine to check arithmetic parser in c2x
 *
 * compile with "gcc test_parser.c parse.c -lm"
 */

#include<stdio.h>
#include<stdlib.h>

int ascan(char *in, double *result);
int debug;

int main(int argc, char **argv){
  double x;
  char *test;
  size_t i,n;
  
  test=NULL; n=0;
  debug=6;
  
  printf("Please input line to parse and evaluate: ");
  i=getline(&test,&n,stdin);
  test[i-1]=0; /* Remove newline */
  
  ascan(test,&x);

  printf("Final result: %f\n",x);
}

void error_exit(char *msg){
  fprintf(stderr,"Aborting: %s\n",msg);
  exit(1);
}


