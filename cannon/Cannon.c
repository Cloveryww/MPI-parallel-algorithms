#include<stdio.h>
#include<malloc.h>
#include<stdlib.h>
#include<mpi.h>
#include<pthread.h>
#include<math.h>
#include<cstring>
int myrank, p;

// Compute C = A*B. A is a n1*n2 matrix. B is a n2*n3 matrix.
void matmul(double* A, double* B, double* C, int n1, int n2, int n3)//做矩阵乘法，结果累加到C矩阵中(需要保证C矩阵初始化过)
{
  int i,j,k;
  //简单的串行矩阵乘法
  for (i = 0; i < n1; i++) {
    for (j = 0; j < n3; j++) {
      for (k = 0; k < n2; k++) {
        C[i*n3+j]+=A[i*n2+k]*B[k*n3+j];
      } 
    }
  } 
}


int setup(int argc,char**argv,double** fstreama,double** fstreamb,int* dim)
{
	FILE* fha;
	FILE* fhb;
	int n1,n2,n3;
	int re=1;
	if (!(fha = fopen(argv[1], "r"))) //打开存储A矩阵的文件
	{
		printf("Can't open file %s, Errno=%d\n", argv[1], 1);//打开失败输出信息
		return -1;
	}
	if(fread(&n1,sizeof(int),1,fha)==0)//读取矩阵的行数
	{
		printf("fread error1!\n");
		return -1;
	}
	if(fread(&n2,sizeof(int),1,fha)==0)//读取矩阵的列数
	{
		printf("fread error2!\n");
		return -1;
	}
	*fstreama = (double *) malloc (n1*n2*sizeof(double));//为矩阵申请内存
	if(fread(*fstreama,sizeof(double),n1*n2,fha)==0)//读取矩阵内容
	{
		printf("fread error3!\n");
		return -1;
	}
	
	fclose(fha);//关闭矩阵文件
	
	if (!(fhb = fopen(argv[2], "r"))) //打开存储A矩阵的文件
	{
		printf("Can't open file %s, Errno=%d\n", argv[2], 2);//打开失败输出信息
		return -1;
	}
	if(fread(&n2,sizeof(int),1,fhb)==0)//读取矩阵的行数
	{
		printf("fread error4!\n");
		return -1;
	}
	if(fread(&n3,sizeof(int),1,fhb)==0)//读取矩阵的列数
	{
		printf("fread error5!\n");
		return -1;
	}
	*fstreamb = (double *) malloc (n2*n3*sizeof(double));//为矩阵申请内存
	if(fread(*fstreamb,sizeof(double),n2*n3,fhb)==0)//读取矩阵内容
	{
		printf("fread error6!\n");
		return -1;
	}
	
	fclose(fhb);//关闭矩阵文件	
	dim[0] = n1;//返回矩阵的大小参数
	dim[1] = n2;//返回矩阵的大小参数
	dim[2] = n3;//返回矩阵的大小参数
	return 0;
}

void scatter_matrix(double* matrixbuf, int rows, int cols, double* local_matrix, int rootp)//将矩阵划分为小矩阵块并分发到各个节点
{
	int row, column, i, j, count;
	int maxrows_block = (rows + rootp - 1)/rootp;//小A矩阵块行数的最大值
	int maxcols_block = (cols + rootp - 1)/rootp;//小矩阵块列数的最大值
	double * matrixbuf2 = NULL;//用来格式化原矩阵的缓冲区
	MPI_Status status;//返回通信的状态
	
	if(myrank == 0)//0号线程
	{
		if(!(matrixbuf2 = (double *)malloc(maxcols_block*maxrows_block*rootp*rootp*sizeof(double))))//为缓冲区申请内存
		{
			printf("Memory allocation failed\n");
		}
		//将矩阵转化为按块连续存放的形式，方便分发每块小矩阵，同时对于边界没有对齐的小矩阵，补零对齐，方便计算
		count = 0;
		for (i = 0; i < rootp; i++){
			for (j = 0; j < rootp; j++){
				if(i!=(rootp-1)&&j==(rootp-1))//特殊处理除了最后一行以外的最后一列
				{
					for (row = 0; row < maxrows_block; row++){
						for (column = 0; column < maxcols_block; column++){
							if((j * maxcols_block + column)>=cols)//补零对齐
							{
								matrixbuf2[count] = 0;
							}else{
								matrixbuf2[count] = matrixbuf[(i * maxrows_block + row ) * cols +j * maxcols_block + column];
							}
							count++;
						}
					}
				}else if(i==(rootp-1)&&j!=(rootp-1))//特殊处理除了最后一列以外的最后一行
				{
					for (row = 0; row < maxrows_block; row++){
						for (column = 0; column < maxcols_block; column++){
							if((i * maxrows_block + row)>=rows)//补零对齐
							{
								matrixbuf2[count] = 0;
							}else{
								matrixbuf2[count] = matrixbuf[(i * maxrows_block + row)*cols + j * maxcols_block + column];
							}
							count++;
						}
					}
				}else if(i==(rootp-1)&&j==(rootp-1))//特殊处理最后一列最后一行的那个块
				{
					for (row = 0; row < maxrows_block; row++){
						for (column = 0; column < maxcols_block; column++){
							if(((j * maxcols_block + column)>=cols) || ((i * maxrows_block + row)>=rows))//补零对齐
							{
								matrixbuf2[count] = 0;
							}else{
								matrixbuf2[count] = matrixbuf[(i * maxrows_block + row) * cols + j * maxcols_block + column];
							}
							count++;
						}
					}
				}else{//普通的块
					for (row = 0; row < maxrows_block; row++){
						for (column = 0; column < maxcols_block; column++){
							matrixbuf2[count] = matrixbuf[(i * maxrows_block + row)*cols + j * maxcols_block + column];
							count++;
						}
					}
				}
			}
		}
		if(count!=maxcols_block*maxrows_block*rootp*rootp)//检查是否出错
		{
			printf("scatter_matrix error!\n");
			return ;
		}
		//将属于本地的那个块留下来
		for(i = 0; i < maxrows_block*maxcols_block; i++)
		{
			local_matrix[i] = matrixbuf2[i];
		}
		//分发其他块到对应的线程
		for(i = 1; i < rootp*rootp; i++)
		{
			MPI_Send((matrixbuf2 + (i * maxcols_block * maxrows_block)), maxcols_block * maxrows_block, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}		
	} else {//非0号线程
		MPI_Recv(local_matrix, maxcols_block * maxrows_block , MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);//非零线程接受0线程发送的小矩阵块
	}
	if(matrixbuf2!=NULL){//释放缓冲区
		free(matrixbuf2);
	}
	return;
}
//A and bufA : n1_block*n2_block;  B and bufB : n2_block*n3_block
//进行cannon算法，各个节点在本地进行矩阵乘法，并交换矩阵块，进行循环，直到计算完毕
void cannon(double* A, double* bufA, double* B, double* bufB, double* C, int n1_block, int n2_block, int n3_block, int rootp)
{
	MPI_Request send_a_req, send_b_req;
	MPI_Status send_a_status, send_b_status, recv_a_status, recv_b_status;
	int cycle_count;
	int rank_next_a,rank_next_b;
	int rank_last_a,rank_last_b;
	int curRowP,curColP,i,j;
	int tag=0;//表示当前正确数据再A中还是bufA中，0表示在A中，1表示在bufA中
	//先初始化各个块，即A_ij循环左移i步，B_ij循环上移j步，C_ij初始化为零
	
	//初始化矩阵C，值全部设为0
	for(i=0;i<n1_block;i++)
	{
		for(j=0;j<n3_block;j++)
		{
			C[i*n3_block+j] = 0;
		}
	}

	//循环传播小矩阵
	curRowP = myrank/rootp;//当前节点所在行
	curColP = myrank%rootp;//当前节点所在列
	
	//获得左移i步后的节点号
	if((curColP-curRowP)<0)
	{
		rank_next_a = myrank+rootp-curRowP;
	}else{
		rank_next_a = myrank-curRowP;
	}
	
	
	//获得上移j步后的节点号
	if((curRowP-curColP)<0)
	{
		rank_next_b = myrank -curColP*rootp + rootp*rootp;
	}else{
		rank_next_b = myrank -curColP*rootp;
	}
	
	//获得接受左移i步的节点号
	if((curColP+curRowP)>=rootp)
	{
		rank_last_a = myrank+curRowP-rootp;
	}else{
		rank_last_a = myrank+curRowP;
	}
	
	
	//获得接受上移j步的节点号
	if((curRowP+curColP)>=rootp)
	{
		rank_last_b = myrank + curColP*rootp - rootp*rootp;
	}else{
		rank_last_b = myrank + curColP*rootp;
	}
	
	//非阻塞发送矩阵，如果不需要移动，则直接本地memcpy
	if(rank_next_a!=myrank)
	{
		MPI_Isend(A, n1_block*n2_block, MPI_DOUBLE, rank_next_a, 0, MPI_COMM_WORLD, &send_a_req);//非阻塞发送矩阵A，避免死锁
	}else
	{
		memcpy(bufA, A, n1_block*n2_block*sizeof(double));//本地直接memcpy
	}
	if(rank_next_b!=myrank)
	{
		MPI_Isend(B, n2_block*n3_block, MPI_DOUBLE, rank_next_b, 0, MPI_COMM_WORLD, &send_b_req);//非阻塞发送矩阵B，避免死锁
	}else
	{
		memcpy(bufB, B, n2_block*n3_block*sizeof(double));//本地直接memcpy
	}
	
	//阻塞接受矩阵
	if(rank_last_a!=myrank)
	{
		MPI_Recv(bufA, n1_block*n2_block, MPI_DOUBLE, rank_last_a, 0, MPI_COMM_WORLD, &recv_a_status);//阻塞接受矩阵A
	}
	if(rank_last_b!=myrank)
	{
		MPI_Recv(bufB, n2_block*n3_block, MPI_DOUBLE, rank_last_b, 0, MPI_COMM_WORLD, &recv_b_status);//阻塞接受矩阵B
	}

	//阻塞等待发送矩阵结束
	if(rank_next_a!=myrank)
	{
		MPI_Wait(&send_a_req, &send_a_status);//阻塞发送矩阵A到结束
	}
	if(rank_next_b!=myrank)
	{
		MPI_Wait(&send_b_req, &send_b_status);//阻塞发送矩阵B到结束
	}
	
	MPI_Barrier(MPI_COMM_WORLD);//同步
	tag=1;
	if(myrank%rootp==0)//第一列的节点
	{
		rank_next_a = myrank+rootp-1;
	}else{
		rank_next_a = myrank-1;
	}
	
	if(myrank/rootp==0)//第一行的节点
	{
		rank_next_b = myrank+rootp*(rootp-1);
	}else{
		rank_next_b = myrank - rootp;
	}
	
	if(myrank%rootp==(rootp-1))//最后一列的节点
	{
		rank_last_a = myrank-rootp+1;
	}else{
		rank_last_a = myrank+1;
	}
	
	if(myrank/rootp==(rootp-1))//最后一行的节点
	{
		rank_last_b = myrank-rootp*(rootp-1);
	}else{
		rank_last_b = myrank + rootp;
	}
	//循环，每次做当前块的乘加运算，并使得A_ij循环左移1步，B_ij循环上移1步
	for(cycle_count = 0; cycle_count < rootp; cycle_count++)
	{
		if(tag==1)//数据在bufA中
		{
			matmul(bufA, bufB, C, n1_block, n2_block, n3_block);//做当前节点的矩阵乘法
			//循环传播小矩阵
			MPI_Isend(bufA, n1_block*n2_block, MPI_DOUBLE, rank_next_a, 0, MPI_COMM_WORLD, &send_a_req);//非阻塞发送矩阵A，避免死锁
            		MPI_Isend(bufB, n2_block*n3_block, MPI_DOUBLE, rank_next_b, 0, MPI_COMM_WORLD, &send_b_req);//非阻塞发送矩阵B，避免死锁

            		MPI_Recv(A, n1_block*n2_block, MPI_DOUBLE, rank_last_a, 0, MPI_COMM_WORLD, &recv_a_status);//阻塞接受矩阵A
            		MPI_Recv(B, n2_block*n3_block, MPI_DOUBLE, rank_last_b, 0, MPI_COMM_WORLD, &recv_b_status);//阻塞接受矩阵B

            		MPI_Wait(&send_a_req, &send_a_status);//阻塞发送矩阵A到结束
            		MPI_Wait(&send_b_req, &send_b_status);//阻塞发送矩阵B到结束
			tag = 0;
		}else{//数据在A中
			matmul(A, B, C, n1_block, n2_block, n3_block);//做当前节点的矩阵乘法
			//循环传播小矩阵
			MPI_Isend(A, n1_block*n2_block, MPI_DOUBLE, rank_next_a, 0, MPI_COMM_WORLD, &send_a_req);//非阻塞发送矩阵A，避免死锁
            		MPI_Isend(B, n2_block*n3_block, MPI_DOUBLE, rank_next_b, 0, MPI_COMM_WORLD, &send_b_req);//非阻塞发送矩阵B，避免死锁

            		MPI_Recv(bufA, n1_block*n2_block, MPI_DOUBLE, rank_last_a, 0, MPI_COMM_WORLD, &recv_a_status);//阻塞接受矩阵A
            		MPI_Recv(bufB, n2_block*n3_block, MPI_DOUBLE, rank_last_b, 0, MPI_COMM_WORLD, &recv_b_status);//阻塞接受矩阵B

            		MPI_Wait(&send_a_req, &send_a_status);//阻塞发送矩阵A到结束
            		MPI_Wait(&send_b_req, &send_b_status);//阻塞发送矩阵B到结束
			tag = 1;
		}
		MPI_Barrier(MPI_COMM_WORLD);//同步
	}
	return;
}

//gather_matrix((double*)(fstreamc + sizeof(int)*2), n1, n3, C, rootp);
//将各个节点的小矩阵C收集到0号节点
void gather_matrix(double* matrixCbuf, int rows, int cols, double* local_C, int rootp, int rows_block_pad, int cols_block_pad)
{
	int curRow, curCol, i, j, curP;
	MPI_Status status;
	double * matrixC_pad = NULL;//有零填充的矩阵C
	if(myrank == 0) {//0号线程 
		if(!(matrixC_pad = (double *)malloc(rows_block_pad*cols_block_pad*rootp*rootp*sizeof(double))))//为缓冲区申请内存
		{
			printf("Memory allocation failed\n");
		}
		//将本地计算结果直接复制过来
		for(i = 0; i < rows_block_pad * cols_block_pad; i++){
			matrixC_pad[i] = local_C[i];
		}
		//接受其他非0线程的计算结果
		for(i = 1; i < rootp*rootp; i++){
			MPI_Recv(matrixC_pad + (i * rows_block_pad * cols_block_pad), rows_block_pad * cols_block_pad, MPI_DOUBLE, i, 0,MPI_COMM_WORLD, &status);
		}
		
		//重新整理矩阵C，除去零填充，并且重新整理顺序
		for(i=0;i<rows;i++)
		{
			for(j=0;j<cols;j++)
			{
				curP = (i/rows_block_pad)*rootp+(j/cols_block_pad);//属于第几个节点，从0开始
				curRow = i%rows_block_pad;//属于小矩阵的第几行
				curCol = j%cols_block_pad;//属于小矩真的第几列
				matrixCbuf[i * cols + j] = matrixC_pad[curP * rows_block_pad * cols_block_pad +curRow*cols_block_pad+curCol];
			}
		}		
	} else {//非0号线程
		MPI_Send(local_C,rows_block_pad * cols_block_pad, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);//给0号线程发送计算结果
	}
	if(matrixC_pad!=NULL)
	{
		free(matrixC_pad);//释放缓冲区
	}
	return ;
}
int main(int argc, char** argv)
{
	double elapsed_time;
	// Suppose A:n1xn2, B:n2xn3. n1~n3 are read from input files
	int n1, n2, n3,rootp;
	// Buffers for matrix A, B, C. Because A, B will be shifted, so they each have two buffers
	double *A, *B, *C, *bufA, *bufB;
	// On proc 0, buffers to cache matrix files of A, B and C
	double *fstreama, *fstreamb;
	char *fstreamc;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	rootp = sqrt(p);
	if (p != rootp*rootp) {
		printf("Processor number must be a square!\n");
	}
	// On proc 0, preprocess the command line, read in files for A, B and
	// put their sizes in dim[].
	int dim[3];
	if (myrank == 0) {//0号线程负责从文件中读取矩阵A和B以及他们的大小信息
		if (setup(argc, argv, &fstreama, &fstreamb, dim)!=0) {
			MPI_Finalize(); // Something error during preprocessing
			exit(-1);
		}
	}
	MPI_Bcast(dim, 3, MPI_INT, 0, MPI_COMM_WORLD);//0号线程将A和B矩阵的size广播给所有线程
	n1 = dim[0];//A： n1*n2
	n2 = dim[1];//B:  n2*n3
	n3 = dim[2];

	// Allocate memories for A, B, C, bufA and bufB.
	// Suppose an m*n matrix is 2D block-distributed on a rootp*rootp processor grid.
	// If rootp doesn't divide m or n, then submatrixes won't have the same size.
	// Because we will shift A, B, so we allocate memories according to the max
	// rows and cols of A and B.
	//因为有可能rootp不能整除n1,n2,n3,所以在申请内存的时候考虑最大的块的大小
	int maxrows_a = (n1 + rootp - 1)/rootp;//A矩阵块行数的最大值
	int maxcols_a = (n2 + rootp - 1)/rootp;//A矩阵块列数的最大值
	int maxrows_b = maxcols_a;//B矩阵块行数的最大值
	int maxcols_b = (n3 + rootp - 1)/rootp;//B矩阵块列数的最大值
	int bufA_size = sizeof(double)*maxrows_a*maxcols_a;//大小为一个A矩阵块的大小
	int bufB_size = sizeof(double)*maxrows_b*maxcols_b;//大小为一个B矩阵块的大小
	int bufC_size = sizeof(double)*maxrows_a*maxcols_b;//大小为一个C矩阵块的大小
	char* buf;
	int i;
	if(!(buf = (char *)malloc(bufA_size*2 + bufB_size*2 + bufC_size)))//申请两个A矩阵块，两个B矩阵块，和一个C矩阵块
	{
		printf("Memory allocation failed\n");
	}
	//或者以下4个缓存区的指针位置
	A = (double*)buf;
	bufA = (double*) (buf + bufA_size);
	B = (double*) (buf + bufA_size*2);
	bufB = (double*) (buf + bufA_size*2 + bufB_size);
	C = (double*) (buf + bufA_size*2 + bufB_size*2);
	// Proc 0 scatters A, B to other procs in a 2D block distribution fashion
	scatter_matrix((double*)fstreama, n1, n2, A, rootp);//0号线程分发A矩阵块到各个线程
	MPI_Barrier(MPI_COMM_WORLD);//同步
	scatter_matrix((double*)fstreamb, n2, n3, B, rootp);//0号线程分发B矩阵块到各个线程
	MPI_Barrier(MPI_COMM_WORLD);//同步
	elapsed_time = MPI_Wtime();//记录计算开始的时间戳
	// Compute C=A*B by Cannon algorithm
	cannon(A, bufA, B, bufB, C, maxrows_a,maxcols_a,maxcols_b, rootp);
	MPI_Barrier(MPI_COMM_WORLD);//同步
	elapsed_time = MPI_Wtime() - elapsed_time;//记录计算所用的时间
	// Proc 0 gathers C from other procs and write it out
	FILE* fhc;
	int fsizec = sizeof(int)*2 + sizeof(double)*n1*n3;//存储C矩阵以及两个大小参数的空间大小
	if(myrank == 0)
	{
		if (!(fhc = fopen(argv[3], "w"))) //打开输出C矩阵的文件
		{
			printf("Can't open file %s, Errno=%d\n", argv[3], 3);//打开失败输出信息
			MPI_Finalize();
		}
		fstreamc = (char *)malloc(fsizec);//申请存储矩阵C的内存空间
		((int*)fstreamc)[0] = n1;//记录矩阵C的行数
		((int*)fstreamc)[1] = n3;//记录矩阵C的列数
	}
	gather_matrix((double*)(fstreamc + sizeof(int)*2), n1, n3, C, rootp, maxrows_a, maxcols_b);//聚集计算结果，其他线程将自己的C矩阵块发送给线程0
	MPI_Barrier(MPI_COMM_WORLD); // Make sure proc 0 read all it needs
	if(myrank == 0)
	{
		printf("Cannon algrithm: multiply a %dx%d with a %dx%d, use %.2f(s)\n",n1, n2, n2, n3, elapsed_time);
		fwrite(fstreamc, sizeof(char), fsizec, fhc);//线程0将矩阵C写入文件
		fclose(fhc);//关闭文件
		free(fstreama);//释放内存
		free(fstreamb);//释放内存
		free(fstreamc);//释放内存
	}
	free(buf);//释放存储小矩阵块的内存空间
	MPI_Finalize();
	return 0;
}
