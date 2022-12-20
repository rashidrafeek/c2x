/* A poor man's dictionary, with order(N) search */

#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include "c2xsf.h"

void *dict_get(struct dct *dict, char *key){

  while(dict){
    if ((dict->key)&&(!strcmp(dict->key,key))) return dict->value;
    dict=dict->next;
  }
  return NULL;
}

void dict_add(struct dct *dict, char *key, void *value){

  if (dict_get(dict,key)){
    while(dict){
      if ((dict->key)&&(!strcmp(dict->key,key))){
        dict->value=value;
        return;
      }
    }
    error_exit("Confusion in dict_add");
  }

  if (!dict) error_exit("Null pointer passed to dict_add");

  while((dict->key)&&(dict->next)) dict=dict->next;

  if (dict->key){
    dict->next=malloc(sizeof(struct dct));
    if (!dict->next) error_exit("Malloc error for dict");
    dict=dict->next;
    dict->next=NULL;
  }
  dict->key=malloc(strlen(key)+1);
  if (!dict->key) error_exit("Malloc error for dict key");
  strcpy(dict->key,key);
  dict->value=value;
}

/* Concatenate a string to an existing entry, or make a new entry
 * if no existing entry. Use malloc() always, so can use realloc,
 * even if called with string constants.
 */
void dict_strcat(struct dct *dict, char *key, char *value){
  char *p;

  if (!dict) error_exit("Null pointer passed to dict_strcat");

  if (!dict_get(dict,key)){
    p=malloc(strlen(value)+1);
    if (!p) error_exit("Malloc failed in dict_strcat");
    strcpy(p,value);
    dict_add(dict,key,p);
    return;
  }

  while(dict){
    if ((dict->key)&&(!strcmp(dict->key,key))) break;
    dict=dict->next;
  }

  if (strcmp(dict->key,key))
    error_exit("Confusion in dict_strcat");

  dict->value=realloc(dict->value,strlen(dict->value)+strlen(value)+1);

  if (!dict->value)
    error_exit("Realloc() failed in dict_strcat");

  strcat(dict->value,value);
  
}
