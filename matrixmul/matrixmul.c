#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "mpi.h"
#define M 5
#define N 5
#define P 5
typedef struct threadArg
{
	int tid;
	float (*B)[P];
	float* A_row;
	float* C_row;
	int numthreads;
}threadArg;
//计算线程平均分配B的所有列B中的一列与A行相乘，计算结果存入C一行中的对应位置
void* worker(void* arg) {
	int i, j;
	threadArg* myarg = (threadArg*)arg;//获取当前线程计算参数
	for (i = myarg->tid; i < P; i += myarg->numthreads) {//计算当前线程所分配的任务
		myarg->C_row[i] = 0.0;//清零
		for (j = 0; j < N; j++) {
			myarg->C_row[i] += myarg->A_row[j] * myarg->B[j][i];//两个向量相乘
		}
	}
	return NULL;
}

int main(int argc, char *argv[]) {
	float A[M][N], B[N][P], C[M][P]; //变量声明
	int myid, numprocs, numthreads;
	int i, j, numsend;
	int sender;
	pthread_t *tids;
	float *A_row;
	float *C_row;
	threadArg *targs;
	MPI_Status status;
	MPI_Init(&argc, &argv);//MPI初始化
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);//获取当前进程id
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);//获取通信域进程数

	if (!myid) {//0号进程初始化矩阵A和矩阵B
		for (i = 0; i < M; i++)//初始化矩阵A
			for (j = 0; j < N; j++)
				A[i][j] = i * j + 1;
		for (i = 0; i < N; i++)//初始化矩阵B
			for (j = 0; j < P; j++)
				B[i][j] = i * j + 1;
		//打印矩阵A 和矩阵B
		printf("matrix A = \n");
		for (i = 0; i < M; i++)//初始化矩阵A
		{
			for (j = 0; j < N; j++)
				printf("%f	",A[i][j]);
			printf("\n");
		}
		printf("matrix B = \n");
		for (i = 0; i < N; i++)//初始化矩阵B
		{
			for (j = 0; j < P; j++)
				printf("%f	",B[i][j]);
			printf("\n");
		}
	}
	//广播矩阵B给非0进程
	MPI_Bcast(B[0], N * P, MPI_FLOAT, 0, MPI_COMM_WORLD);
	if (!myid) { /*0进程：分配任务和回收结果*/
		j = (numprocs - 1) < M ? (numprocs - 1) : M;//取计算进程数和矩阵A行数的较小值
		for (i = 1; i <= j; i++) {//发送A矩阵的前j行给j个线程
			MPI_Send(A[i - 1], N, MPI_FLOAT, i, 99, MPI_COMM_WORLD);
			printf("envid 0 sent no.%d row\n",i-1);
		}
		numsend = j;//记录已经发送的A矩阵的行数

		for (i = 1; i <= M; i++) {//回收计算结果
			sender = (i - 1) % (numprocs - 1) + 1;//获取计算A矩阵第i行的进程的id
			MPI_Recv(C[i - 1], P, MPI_FLOAT, sender, 100, MPI_COMM_WORLD, &status);//回收计算结果

			if (numsend < M) {//如果矩阵A中还有行没有分配计算，则分配给刚刚完成计算任务的进程进行计算
							  //   i ==》numsend
				MPI_Send(A[numsend], N, MPI_FLOAT, sender, 99, MPI_COMM_WORLD);//发送任务
				printf("envid 0 sent no.%d row\n",numsend);
				numsend++;//已经发送的A矩阵的行数加一
			}
			else {
				MPI_Send(&j, 0, MPI_INT, sender, 0, MPI_COMM_WORLD);//A矩阵全部计算结束，终止计算进程
			}
		}
		if((numprocs-1)>M)
		{
			for(i=M+1;i<numprocs;i++)
			{
				MPI_Send(&j, 0, MPI_INT, i, 0, MPI_COMM_WORLD);//A矩阵全部计算结束，终止计算进程
			}
		}
		//打印计算结果
		printf("\nmatrix C = \n");
		for (i = 0; i < M; i++)//初始化矩阵A
		{
			for (j = 0; j < P; j++)
				printf("%f	",C[i][j]);
			printf("\n");
		}
	}
	else {//从进程
		  /*从进程(myid > 0)：接收0进程发来的任务，计算完毕发回主进程*/
		numthreads = get_nprocs();//获取进程所在节点的CPU数
		tids = (pthread_t*)malloc(numthreads * sizeof(pthread_t));//用来存放线程ID
		A_row = (float*)malloc(N * sizeof(float));//用来存放A矩阵的一行
		C_row = (float*)malloc(P * sizeof(float));//用来存放计算结果
		targs = (threadArg*)malloc(numthreads * sizeof(threadArg));//用来存放线程计算参数

		for (i = 0; i < numthreads; i++) {//初始化线程计算参数
			targs[i].tid = i;//线程id
			targs[i].B = B;//矩阵B的地址
			targs[i].A_row = A_row;//矩阵A某行的开始地址
			targs[i].C_row = C_row;//保存计算结果的的开始地址
			targs[i].numthreads = numthreads;//当前节点的线程数
		}

		while (1)
		{
			MPI_Recv(A_row, N, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);//接受第0号进程发送的A矩阵的某行
			if (status.MPI_TAG == 0)//发送的是终止消息，则跳出死循环，结束进程
				break;
			printf("env %d get one new row\n",myid);
			for (i = 0; i < numthreads; i++) //创建计算线程
			{
				pthread_create(&tids[i], NULL, worker, &targs[i]);
			}
			for (i = 0; i < numthreads; i++) //同步计算线程
			{
				pthread_join(tids[i], NULL);
			}
			MPI_Send(C_row, P, MPI_FLOAT, 0, 100, MPI_COMM_WORLD);//发送计算结果给0号进程
			printf("Over env %d return one row\n",myid);
		}
		printf("envid %d killed\n",myid);
	}
	MPI_Finalize(); 
	return 0;
}

