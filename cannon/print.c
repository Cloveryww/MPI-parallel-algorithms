// print.cxx
// print the double precision matrix from the input file

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

int main(int argc, char** argv)
{
  FILE* fh;
  int i,j;
  if (argc < 2) {
    printf("Invalid arguments!\n");
    printf("Run the program as ./print filename\n");
    exit(-1);
  } 

  if(!(fh = fopen(argv[1], "r"))) {
    printf("Can't open file %s\n", argv[1]);
    exit(-1);
  }

  struct stat fstat;
  int n1, n2, fsize;
  char* fstream;

  stat(argv[1], &fstat);
  fsize = fstat.st_size;
  fstream = (char *)malloc(fsize);
  fread(fstream, sizeof(char), fsize, fh);

  n1 = ((int*)fstream)[0]; 
  n2 = ((int*)fstream)[1]; 

  if (n1 <=0 || n2 <=0) {
    printf("Matrix size error, %dx%d\n", n1, n2);
    exit(-1);
  } 
 
  if (fsize < (sizeof(int)*2 + sizeof(double)*n1*n2)) {
    printf("Actual size mismatches with stated size\n");
    exit(-1);
  }

  double* A = (double*)(fstream + sizeof(int)*2);

  printf("       ---- %s: %d*%d Matrix -----\n", argv[1], n1, n2);
  for (i = 0; i < n1; i++) {
    for (j = 0; j < n2; j++) {
      printf("%.4f  ", *(A+i*n2+j)); // A[i,j]
    }
    printf("\n");
  }

  free(fstream);
  fclose(fh);
  return 0;
}
