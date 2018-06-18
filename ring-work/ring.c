#include<stdio.h>
#include "mpi.h"
int main(int argc,char *argv[])
{
	int myid,numprocs;
	char message[100];
	char message2[100];
	double begintime,endtime;	
	int nums = 1;
	int i;
	double alltime = 0;
	MPI_Status status;
	MPI_Init(&argc,&argv);//初始化
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);//获取当前进程的id
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);//获取当前通信域的进程数
	if(myid==0)//如果是第0个进程
	{
		for(i=0;i<nums;i++)
		{
			strcpy(message,"This is the ring message!");//初始化消息
			begintime = MPI_Wtime();//记录开始时间戳
			MPI_Send(message,strlen(message)+1,MPI_CHAR,1,58,MPI_COMM_WORLD);//发送消息
			printf("envid 0: Send message\n");
			MPI_Recv(message2,100,MPI_CHAR,numprocs-1,58,MPI_COMM_WORLD,&status);//接受消息
			endtime = MPI_Wtime();//记录结束时间戳
			alltime = alltime + endtime - begintime;//记录总时间
			printf("envid 0: Recv message\n");
		//	printf("envid 0: Finish one ring! time is %f\n",endtime-begintime);
		//	printf("envid 0: message : %s\n",message2);
		}
		printf("envid 0: Finish 100 ring! Average time is %f\n",alltime/nums);
	}
	else//如果是其他进程
	{
		for(i=0;i<nums;i++)
		{
			MPI_Recv(message2,100,MPI_CHAR,myid-1,58,MPI_COMM_WORLD,&status);//接收消息
			printf("envid %d: Recv message\n",myid);
			printf("envid %d: Send message\n",myid);
			if(myid!=numprocs-1)//如果不是最后一个进程
			{
				MPI_Send(message2,strlen(message2)+1,MPI_CHAR,myid+1,58,MPI_COMM_WORLD);//发给第0个进程
			}
			else{//是最后一个进程
				MPI_Send(message2,strlen(message2)+1,MPI_CHAR,0,58,MPI_COMM_WORLD);//发给下一个进程
			}
		}
	}
	MPI_Finalize();//结束
	return 0;
}
