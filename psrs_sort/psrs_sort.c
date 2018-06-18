#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define INIT_TYPE 10
#define ALLTOONE_TYPE 100
#define ONETOALL_TYPE 200
#define MULTI_TYPE 300
#define RESULT_TYPE 400
#define RESULT_LEN 500
#define MULTI_LEN 600

int Spt;
long DataSize;
int *arr,*arr1;
int mylength;
int *index;
int *temp1;

main(int argc,char* argv[])
{
    long BaseNum = 1;
    int PlusNum;
    int MyID;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&MyID);

    PlusNum=200;
    DataSize = BaseNum*PlusNum;//初始化数组总长度

    if (MyID==0)
        printf("The DataSize is : %lu\n",PlusNum);
    Psrs_Main();

    if (MyID==0)
        printf("\n");

    MPI_Finalize();
}


Psrs_Main( )
{
    int i,j;
    int MyID,SumID;
    int n,c1,c2,c3,c4,k,l;
    FILE * fp;
    int ready;
    MPI_Status status[32*32*2];
    MPI_Request request[32*32*2];

    MPI_Comm_rank(MPI_COMM_WORLD,&MyID);//获得当前进程的id
    MPI_Comm_size(MPI_COMM_WORLD,&SumID);//获得总的进程数

    Spt = SumID-1;

	/*初始化参数*/
    arr = (int *)malloc(2*DataSize*sizeof(int));
    if (arr==0) merror("malloc memory for arr error!");//为存储数组申请内存空间
    arr1 = &arr[DataSize];//为归并排序申请内存空间

    if (SumID>1)
    {
        temp1 = (int *)malloc(sizeof(int)*SumID*Spt);//为正则抽样和主元申请内存空间
        if (temp1==0) merror("malloc memory for temp1 error!");
        index = (int *)malloc(sizeof(int)*2*SumID);//为根据主元划分段后的分界下标数组申请内存空间
        if (index==0) merror("malloc memory for index error!");
    }

    MPI_Barrier( MPI_COMM_WORLD);

    mylength = DataSize / SumID;//每个进程排序数据的长度
    srand(MyID);//初始化随机数生成种子

    printf("This is node %d \n",MyID);
    printf("On node %d the input data is:\n",MyID);
    for (i=0;i<mylength;i++)//生成数据
    {
        arr[i] = (int)rand();
        printf("%d : ",arr[i]);//显示数据
    }
    printf("\n");

	/*每个处理器将自己的n/P个数据用串行快速排序(Quicksort)，得到一个排好序的序列，对应于算法13.5步骤（1）*/
    MPI_Barrier( MPI_COMM_WORLD);
    quicksort(arr,0,mylength - 1);
    MPI_Barrier( MPI_COMM_WORLD);

	/*每个处理器从排好序的序列中选取第w，2w，3w，…，(P-1)w个共P-1个数据作为代表元素，其中w=n/P*P，对应于算法13.5步骤（2）*/
    if (SumID>1)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        n = (int)(mylength/(Spt+1));
        for (i=0;i<Spt;i++)//正则采样p-1个元素
            temp1[i] = arr[(i+1)*n-1];

        MPI_Barrier(MPI_COMM_WORLD);

        if (MyID==0)
        {
			/*每个处理器将选好的代表元素送到处理器P0中，对应于算法13.5步骤（3） */
            j = 0;
            for (i=1;i<SumID;i++)
                MPI_Irecv(&temp1[i*Spt], sizeof(int)*Spt, MPI_CHAR, i,ALLTOONE_TYPE+i, MPI_COMM_WORLD, &request[j++]);
            MPI_Waitall(SumID-1, request, status);//等待所有通信操作完成之后才返回，否则将一直等待

			/* 处理器P0将上一步送来的P段有序的数据序列做P路归并，再选择排序后的第P-1，2(P-1)，…，(P-1)(P-1)个共P-1个主元，，对应于算法13.5步骤（3）*/
            MPI_Barrier(MPI_COMM_WORLD);
            quicksort(temp1,0,SumID*Spt-1);//对正则采样数据进行快排(没有使用p路归并)
            MPI_Barrier( MPI_COMM_WORLD);

            for (i=1;i<Spt+1;i++)
                temp1[i] = temp1[i*Spt-1];//等间隔选取主元
			/*处理器P0将这P-1个主元播送到所有处理器中，对应于算法13.5步骤（4）*/
            MPI_Bcast(temp1, sizeof(int)*(1+Spt), MPI_CHAR, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        }
        else
        {
            MPI_Send(temp1,sizeof(int)*Spt,MPI_CHAR,0,ALLTOONE_TYPE+MyID, MPI_COMM_WORLD);//将自己选择的正则采样数据发送个0号进程
            MPI_Barrier( MPI_COMM_WORLD);
            MPI_Barrier( MPI_COMM_WORLD);
            MPI_Bcast(temp1, sizeof(int)*(1+Spt), MPI_CHAR, 0, MPI_COMM_WORLD);//接受0号进程返回的p-1个主元
            MPI_Barrier(MPI_COMM_WORLD);
        }

		/*每个处理器根据上步送来的P-1个主元把自己的n/P个数据分成P段，记为处理器Pi的第j+1段，其中i=0,…,P-1，j=0,…,P-1，对应于算法13.5步骤（5）*/
        n = mylength;
        index[0] = 0;
        i = 1;
        while ((arr[0]>=temp1[i])&&(i<SumID))//对每个进程数据中长度为0的段的上下界下标设为0
        {
            index[2*i-1] = 0;
            index[2*i] = 0;
            i++;
        }
        if (i==SumID) index[2*i-1] = n;//处理特殊情况，该进程的数据实际只分了1段
        c1 = 0;
        while (i<SumID)
        {
            c4 = temp1[i];//2分查找第i个分界点
            c3 = n;
            c2 = (int)((c1+c3)/2);
            while ((arr[c2]!=c4)&&(c1<c3))
            {
                if (arr[c2]>c4)//在左半段
                {
                    c3 = c2-1;
                    c2 = (int)((c1+c3)/2);
                }
                else//在右半段
                {
                    c1 = c2+1;
                    c2 = (int)((c1+c3)/2);
                }
            }
            while ((arr[c2]<=c4)&&(c2<n)) c2++;//特殊处理连续几个数都等于分界点的情况
            if (c2==n)//到了最后一个数据，即后面的分段的长度都为0
            {
                index[2*i-1] = n;
                for (k=i;k<SumID;k++)//后面的分段的长度都设为0
                {
                    index[2*k] = 0;
                    index[2*k+1] = 0;
                }
                i = SumID;//分段完毕
            }
            else//找到第i段的分界点
            {
                index[2*i] = c2;
                index[2*i-1] = c2;
            }
            c1 = c2;//继续分下
            c2 = (int)((c1+c3)/2);//？？？似乎没用
            i++;
        }
        if (i==SumID) index[2*i-1] = n;//给最后一个分界赋值

        MPI_Barrier( MPI_COMM_WORLD);//分段完毕

		/*每个处理器送它的第i+1段给处理器Pi，从而使得第i个处理器含有所有处理器的第i段数据(i=0,…,P-1)，，对应于算法13.5步骤（6）*/

        j = 0;
        for (i=0;i<SumID;i++)
        {
            if (i==MyID)//当前段就不用
            {
                temp1[i] = index[2*i+1]-index[2*i];//计算当前进程对应的分段的长度，其不用发送这一段
                for (n=0;n<SumID;n++)//发送其他分段的长度给对应的进程
                    if (n!=MyID)
                {
                    k = index[2*n+1]-index[2*n];
                    MPI_Send(&k, sizeof(int), MPI_CHAR, n, MULTI_LEN+MyID,MPI_COMM_WORLD);
                }

            }
            else
            {
                MPI_Recv(&temp1[i], sizeof(int), MPI_CHAR, i,MULTI_LEN+i, MPI_COMM_WORLD, &status[j++]);//接受其他进程发过来的分段长度，记录长度
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        j = 0;
        k = 0;
        l = 0;
		//开始真正发送数据分段
        for (i=0;i<SumID;i++)
        {
            MPI_Barrier(MPI_COMM_WORLD);

            if (i==MyID)
            {
                for (n=index[2*i];n<index[2*i+1];n++)//复制本进程对应的分段，不需要通信就可以获得
                    arr1[k++] = arr[n];
            }

            MPI_Barrier(MPI_COMM_WORLD);

            if (i==MyID)
            {
                for (n=0;n<SumID;n++)//向其他进程发送对应的分段数据
                    if (n!=MyID)
                {
                    MPI_Send(&arr[index[2*n]], sizeof(int)*(index[2*n+1]-index[2*n]),MPI_CHAR, n, MULTI_TYPE+MyID, MPI_COMM_WORLD);
                }

            }
            else//接受对应i进程发送过来的分段数据
            {
                l = temp1[i];//取处分段长度
                MPI_Recv(&arr1[k], l*sizeof(int), MPI_CHAR, i ,MULTI_TYPE+i, MPI_COMM_WORLD, &status[j++]);//接受分段数据
                k=k+l;
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }
        mylength = k;
        MPI_Barrier(MPI_COMM_WORLD);

		/*每个处理器再通过P路归并排序将上一步的到的数据排序；从而这n个数据便是有序的，，对应于算法13.5步骤（7） */
        k = 0;
        multimerge(arr1,temp1,arr,&k,SumID);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    printf("On node %d the sorted data is : \n",MyID);
	if(k%1==1)
	{
		for (i=0;i<mylength;i++)
			printf("%d : ",arr[i]);
	}else{
		for (i=0;i<mylength;i++)
			printf("%d : ",arr1[i]);
	}
    printf("\n");
}


/*输出错误信息*/
merror(char* ch)
{
    printf("%s\n",ch);
    exit(1);
}


/*串行快速排序算法*/
quicksort(int *datas,int bb,int ee)
{
    int tt,i,j;
    tt = datas[bb];
    i = bb;
    j = ee;

    if (i<j)
    {
		//先扫描一遍，把数组划分成两部分
        while(i<j)
        {
            while ((i<j)&&(tt<=datas[j])) j--;//移动j指针
            if (i<j)
            {
                datas[i] = datas[j];//交换
                i++;//移动i指针
                while ((i<j)&&(tt>datas[i])) i++;//移动i指针
                if (i<j)
                {
                    datas[j] = datas[i];//交换
                    j--;//移动j指针
                    if (i==j) datas[i] = tt;//放回mid
                }
                else datas[j] = tt;//放回mid
            } else datas[i] = tt;//放回mid
        }
		
        quicksort(datas,bb,i-1);//对前一半数组递归调用quicksort
        quicksort(datas,i+1,ee);//对后一半数组递归调用quicksort
    }
}


/*串行多路归并算法*/
multimerge(int *data1,int *ind,int *data,int *iter,int SumID)//data1是原数据，data是排序后的数据，ind是每分段的长度
{
    int i,j,n;

    j = 0;
    for (i=0;i<SumID;i++)//删除长度为0的分段，之后ind数组中存的是各个分段的长度(>0)
        if (ind[i]>0)
    {
        ind[j++] = ind[i];
        if (j<i+1) ind[i] = 0;
    }
    if ( j>1 )//非零分段数大于1才需要归并
    {
        n = 0;
        for (i=0;i<j,i+1<j;i=i+2)
        {
            merge(&(data1[n]),ind[i],ind[i+1],&(data[n]));//相邻两路归并
            ind[i] += ind[i+1];
            ind[i+1] = 0;//清零
            n += ind[i];//指针右移
        }
        if (j%2==1 )//特殊处理最后一段单独的
		{
            for (i=0;i<ind[j-1];i++) 
			{
				data[n]=data1[n];
				n++;
			}
		}
        (*iter)++;
        multimerge(data,ind,data1,iter,SumID);//递归并归
    }
}


merge(int *data1,int s1,int s2,int *data2)//并归两路数据
{
    int i,l,m;

    l = 0;
    m = s1;
    for (i=0;i<s1+s2;i++)
    {
        if (l==s1)
		{
            data2[i]=data1[m++];
        }else if(m==s2+s1)
		{
            data2[i]=data1[l++];
		}
        else if (data1[l]>data1[m])
		{
			data2[i]=data1[m++];
        }else
		{
			data2[i]=data1[l++];
		}
	}
}
