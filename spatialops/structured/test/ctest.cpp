#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cu_lib.h"

int main()
{
  float *a,*b;
  int j;
  int i;

  for ( j=10; j<1000000; j*=2){
    a =( float *) malloc(sizeof(float)*j);
    b =( float *) malloc(sizeof(float)*j);

    for (i=0; i<j; i++){
      a[i] = 2;
      b[i] = 3;
    }

    if(!func_add(a,b,j)){ printf("Failed\n"); return 1; }

    if(!func_mul(a,b,j)){ printf("Failed\n"); return 1; }

    free(a);
    free(b);
  }

  printf("Success\n");

  return 0;
}
