// serial.cxx
// Serial verions of matrix multiplicaton on matrixes read from input files.
// Run as ./serial A B C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/time.h>
#include <sys/times.h>

double wtime()
{
  struct timeval tv; 
  gettimeofday(&tv, NULL);
  return (tv.tv_sec + 1E-6 * tv.tv_usec);
}

// Preprocess the command line, read in matrix A, B from input files, allocate
// memory for buffers, i.e, fstreama, fstreamb to cache them. Suppose A, B's 
// size are n1*n2, n2*n3, then n1~n3 will be stored at dim[0~2]
// Return value 0 means no error occurred during preprocessing, otherwise a
// non-zero returns.
int setup(int argc, char** argv, char* &fstreama, char* &fstreamb, int* dim)
{
  int error = 0;

  if (argc < 4) {
    printf("Invalid arguments!\n"); 
    printf("Usage: ./serial filea fileb filec\n");
    printf("filea, fileb and filec are file names for matrix A, B and C\n");
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

  if (n1 <=0 || n2 <=0 || n3 <= 0 || n4 <=0 || n2 != n3) {
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

  dim[0] = n1;
  dim[1] = n2;
  dim[2] = n4;
 
  fclose(fha); 
  fclose(fhb); 
  return 0;
}

// Compute C = A*B. A is a n1*n2 matrix. B is a n2*n3 matrix.
void matmul(double* A, double* B, double* C, int n1, int n2, int n3)
{
#define A(i,j)  *(A + i*n2 + j)
#define B(i,j)  *(B + i*n3 + j)
#define C(i,j)  *(C + i*n3 + j)

  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n3; j++) {
      C(i,j) = 0.0;
      for (int k = 0; k < n2; k++) {
        C(i,j) += A(i,k)*B(k,j); 
      } 
    }
  } 
}


int main(int argc, char** argv)
{
  // Suppose A's size is n1xn2, B's is n2xn3. n1~n3 are read from input files 
  int n1, n2, n3;

  // Buffers for matrix A, B, C. Because A, B will be shifted, so they
  // each have two buffers
  double *A, *B, *C, *bufA, *bufB;

  // On proc 0, buffers to cache matrix files of A, B and C
  char *fstreama, *fstreamb, *fstreamc;

  double elapsed_time;

  // On proc 0, preprocess the command line, read in files for A, B and
  // put their sizes in dim[].
  int dim[3]; 
  if (setup(argc, argv, fstreama, fstreamb, dim)) {
    exit(-1); // Something error during preprocessing 
  }

  n1 = dim[0];
  n2 = dim[1];
  n3 = dim[2];

  FILE* fhc;
  int fsizec = sizeof(int)*2 + sizeof(double)*n1*n3;
  if (!(fhc = fopen(argv[3], "w"))) {
    printf("Can't open file %s, Errno=%d\n", argv[3], errno);
    exit(-1);
  }
  fstreamc = (char *)malloc(fsizec);
  ((int*)fstreamc)[0] = n1; 
  ((int*)fstreamc)[1] = n3; 


  elapsed_time = wtime();

  matmul((double*)(fstreama + sizeof(int)*2),
    (double*)(fstreamb + sizeof(int)*2),
    (double*)(fstreamc + sizeof(int)*2),
    n1, n2, n3);

  elapsed_time = wtime()- elapsed_time;

  printf("Serial algrithm: multiply a %dx%d with a %dx%d, use %.2f(s)\n", 
    n1, n2, n2, n3, elapsed_time);

  fwrite(fstreamc, sizeof(char), fsizec, fhc);
  fclose(fhc);
  free(fstreama);
  free(fstreamb);
  free(fstreamc);

  return 0;
}
