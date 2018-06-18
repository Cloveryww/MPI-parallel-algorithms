// matgen.cxx
// Generate a m*n double precision matrix file

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(int argc, char** argv)
{
  int m, n;
  char* filename;
  double* a;
  int i,j;
  if (argc < 4) {
    printf("Invalid arguments!\n");
    printf("Run the program as ./matgen m n filename\n");
    printf("m, n are row/col number of the matrix, filename is the file to write\n");
    return 0;
  } else {
    m = atoi(argv[1]);
    n = atoi(argv[2]);
    filename = argv[3];
  }

  int bufsize = sizeof(int)*2 + sizeof(double)*m*n;
  a = (double*) malloc(bufsize);

  ((int*)a)[0] = m;
  ((int*)a)[1] = n;

  double *ptr = (double*)((int*)a + 2);

  srand48(time(NULL)); // Use time as a seed
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      *(ptr + i*n + j) = drand48();
    }
  }

  FILE* file;
  if(!(file = fopen(filename, "w"))) {
    printf("Can't open file %s\n", filename);
  }

  fwrite(a, sizeof(char), bufsize, file);
  fclose(file);

  return 0;
}
