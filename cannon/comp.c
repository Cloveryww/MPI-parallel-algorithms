// comp.cxx
//
// Compare two matrixes if equal.


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>

// Preprocess the command line, read in matrix A, B from input files, allocate
// memory for buffers, i.e, fstreama, fstreamb to cache them.
int setup(int argc, char** argv, char* &fstreama, char* &fstreamb)
{
  int error = 0;

  if (argc < 3) {
    printf("Invalid arguments!\n"); 
    printf("Usage: ./serial filea fileb\n");
    printf("filea, fileb are file names for matrix A, B to be compared\n");
    return 1;
  }

  FILE *fha, *fhb;
  int fsizea, fsizeb;

  if (!(fha = fopen(argv[1], "r"))) {
    printf("Can't open matrix file %s, Errno=%d\n", argv[1], errno);
    return 1;
  }
  
  if (!(fhb = fopen(argv[2], "r"))) {
    printf("Can't open file %s, Errno=%d\n", argv[2], errno);
    return 1;
  }

  struct stat fstata, fstatb;
  stat(argv[1], &fstata);
  stat(argv[2], &fstatb);
  fsizea = fstata.st_size;
  fsizeb = fstatb.st_size;

  fstreama = (char *)malloc(fsizea);
  fstreamb = (char *)malloc(fsizeb);
  fread(fstreama, sizeof(char), fsizea, fha);
  fread(fstreamb, sizeof(char), fsizeb, fhb);

  int n1, n2, n3, n4;
  n1 = ((int*)fstreama)[0]; 
  n2 = ((int*)fstreama)[1]; 
  n3 = ((int*)fstreamb)[0]; 
  n4 = ((int*)fstreamb)[1]; 

  if (n1 <=0 || n2 <=0 || n3 <= 0 || n4 <=0) {
    printf("Matrix size error, %dx%d with %dx%d\n", n1, n2, n3, n4);
    return 1;
  } 
 
  if (fsizea < (sizeof(int)*2 + sizeof(double)*n1*n2)) {
    printf("Actual size of A mismatches with stated size\n");
    return 1;
  }

  if (fsizeb < (sizeof(int)*2 + sizeof(double)*n3*n4)) {
    printf("Actual size of B mismatches with stated size\n");
    return 1;
  }

  fclose(fha); 
  fclose(fhb); 
  return 0;
}

void comp(char* fstreama, char* fstreamb)
{
  int n1, n2, n3, n4;
  int i,j;
  n1 = ((int*)fstreama)[0]; 
  n2 = ((int*)fstreama)[1]; 
  n3 = ((int*)fstreamb)[0]; 
  n4 = ((int*)fstreamb)[1]; 
  
  if (n1 != n3 || n2 != n4) {
    printf("Matrix size mismatch, %dx%d with %dx%d\n", n1, n2, n3, n4);
    return;
  }
  
  double *A = (double*)(fstreama + sizeof(int)*2);
  double *B = (double*)(fstreamb + sizeof(int)*2);

#define A(i,j)  *(A + i*n2 + j)
#define B(i,j)  *(B + i*n2 + j)

  double norm = 0.0;
  for (i = 0; i < n1; i++) {
    for (j = 0; j < n2; j++) {
      double diff = A(i,j)-B(i,j);
      norm += diff*diff;
    }
  }
  norm = sqrt(norm);

  if (norm > 0.000001) {
    printf("Matrix compare failed, with norm = %.8f\n", norm);
  } else {
    printf("Matrix compare succeeded, with norm = %.8f\n", norm);
  } 

}

int main(int argc, char** argv)
{
  // Suppose A's size is n1xn2, B's is n3xn4. n1~n4 are read from input files 
  int n1, n2, n3, n4;

  // Buffers to cache matrix files of A, B
  char *fstreama, *fstreamb;

  // preprocess the command line, read in files for A, B.
  if (setup(argc, argv, fstreama, fstreamb)) {
    exit(-1); // Something error during preprocessing 
  }

  comp(fstreama, fstreamb);

  free(fstreama);
  free(fstreamb);

  return 0;
}
