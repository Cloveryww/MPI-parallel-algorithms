#include<stdlib.h>
#include <malloc.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#define  MAX(m,n)    (m>n?m:n)

typedef struct{
  	int pedlen;
  	int psuffixlen;
  	int pednum;
}pntype;


/*对模式串进行周期分析，并计算相应的next和nextval值*/
void Next(char *W,int patlen,int *nextval,pntype *pped)
{
	int i,j,plen;
   	int *next;

   	if((next=(int *)malloc((patlen+1)*sizeof(int)))==NULL){
     	printf("no enough memory\n");
       	exit(1);
   	}

  	/*计算next和nextval*/    
    next[0]=nextval[0]=-1;
  	j=1;
   	while(j<=patlen){//计算next值
   		i=next[j-1];
      	while(i!=(-1) && W[i]!=W[j-1]) i=next[i];//正常计算next值
      	next[j]=i+1;

      	if(j!=patlen){
          	if( W[j]!=W[i+1])
              	nextval[j]=i+1;//优化计算nextval的值
           	else
    			nextval[j]=nextval[i+1]; //优化计算nextval的值
      	}
        j++;
  	} 

    pped->pedlen=patlen-next[patlen]; //模式串最小周期长度
 	pped->pednum=(int)(patlen/pped->pedlen); //模式串周期数
   	pped->psuffixlen=patlen%pped->pedlen; //最后一段长度 

  	free(next);
}

/*改进的KMP算法*/
void kmp(char *T,char*W,int textlen,int patlen,int *nextval,pntype *pped,int prefix_flag,int matched_num,int *match,int *prefixlen)
{
	int i,j;

  	i=matched_num;              
  	j=matched_num;           

    while(i<textlen)
    {
    	if((prefix_flag==1)&&((patlen-j)>(textlen-i))) {//前面匹配成功但是长度不足的情况
          	break;
		}

        while(j!=(-1) && W[j]!=T[i])  j=nextval[j];//根据next滑动窗口

        if(j==(patlen-1)) {    
			match[i-(patlen-1)]=1;
        	if(pped->pednum+pped->psuffixlen==1)
            j = -1 ; 
   		else                                   
          	j=patlen-1-pped->pedlen; 
        }
   		j++;
      	i++; 
   	}
   	(*prefixlen)=j;//返回是否跨段匹配
}

/*重构模式串以及next函数*/
void Rebuild_info(int patlen,pntype *pped,int *nextval,char *W) 
{ 
	int i; 
   	if (pped->pednum == 1) 
   		memcpy(W+pped->pedlen,W,pped->psuffixlen); 
	else {  
       	memcpy(W+pped->pedlen,W,pped->pedlen);
       	for (i=3; i<=pped->pednum; i++){ 
        	memcpy(W+(i-1)*pped->pedlen,W,pped->pedlen);
           	memcpy(nextval+(i-1)*pped->pedlen,nextval+pped->pedlen,pped->pedlen*sizeof(int));
      	} 

       	if(pped->psuffixlen!=0){
       		memcpy(W+(i-1)*pped->pedlen,W,pped->psuffixlen);
           	memcpy(nextval+(i-1)*pped->pedlen,nextval+pped->pedlen,pped->psuffixlen*sizeof(int));            
       	}
 	} 
} 
 
/*生成文本串*/
void gen_string(int strlen,int pedlen,char *string,int seed)
{
	int suffixlen,num,i,j;

   	srand(seed);//设置随机数seed
   	for(i=0;i<pedlen;i++){//生成随机字母
    	num=rand()%26;
        string[i]='a'+num;
   	}
   	for(j=1;j<(int)(strlen/pedlen);j++)
    	strncpy(string+j*pedlen,string,pedlen);//周期性复制字符
   	if((suffixlen=strlen%pedlen)!=0)//特殊处理最后一段
    	strncpy(string+j*pedlen,string,suffixlen);
	string[strlen]='\0';
}  

/*从文件读入模式串信息*/ 
void GetFile(char *filename,char **place, int *number) 
{ 
	FILE *fp;
    struct stat statbuf;  

    if ((fp=fopen(filename,"rb"))==NULL) {//打开模式串文件
		printf ("Error open file %s\n",filename); 
        exit(0); 
	} 

    fstat(fileno(fp),&statbuf);  
    if(((*place)=(char *)malloc(sizeof(char)*statbuf.st_size)) == NULL){//为模式串申请内存空间
		printf ("Error alloc memory\n");
        exit(1);
 	}
     
   	if (fread((*place),1,statbuf.st_size,fp)!=statbuf.st_size){//读取模式串的长度
		printf ("Error in reading num\n"); 
        exit(0); 
	} 
    fclose (fp); //关闭文件
    (*number)=statbuf.st_size; //返回模式串长度
} 

/*打印运行参数信息*/
void PrintFile_info(char *filename,char *T,int id)
{ 
	FILE *fp; 
	int i;

	if ((fp=fopen(filename,"a"))==NULL){//打开文件
		printf ("Error open file %s\n",filename); 
        exit(0); 
	} 

 	fprintf (fp,"The Text on node %d is %s .\n",id,T); //打印在节点id上的字符串内容
   	
	fclose (fp); 
} 


/*打印匹配结果*/
void PrintFile_res(char *filename,int *t,int len,int init,int id)
{ 
	FILE *fp; 
	int i;

	if ((fp=fopen(filename,"a"))==NULL){//打开输出信息文件
		printf ("Error open file %s\n",filename); 
        exit(0); 
	} 

 	fprintf (fp,"This is the match result on node %d \n",id); 
   	for (i=0; i<=len-1; i++) //输出当前节点的匹配结果
    	if(t[i]==1)
 			fprintf (fp,"(%d)  +\n",i+init); //匹配成功
    	else
  			fprintf (fp,"(%d)  -\n",i+init);//匹配失败
	fclose (fp); //关闭文件
} 

void main(int argc,char *argv[]) 
{ 
	char *T,*W; 
	int	*nextval,*match; 
  	int	textlen,patlen,pedlen,nextlen_send; 
   	pntype pped; 
 	int	i,myid,numprocs,prefixlen,ready; 
  	MPI_Status  status; 

	MPI_Init(&argc,&argv); 
   	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); //获取全局进程数
   	MPI_Comm_rank(MPI_COMM_WORLD,&myid); //获取当前进程id

   	nextlen_send=0;
   	ready=1;//设置标志位
   	/*读如文本串长度*/
   	textlen=atoi(argv[1]);//获取字符串长度参数
   	textlen=textlen/numprocs;//获取每个进程负责的字符串长度
  	pedlen=atoi(argv[2]);//获取字符串周期
   	if((T=(char *)malloc(textlen*sizeof(char)))==NULL){//为字符串申请内存
     	printf("no enough memory\n");
       	exit(1);
   	}
   	gen_string(textlen,pedlen,T,myid);//生成字符串

	if(myid==0){
		PrintFile_info("match_result",T,myid);//0号进程打印运行参数信息
		if(numprocs>1)
			MPI_Send(&ready,1,MPI_INT,1,0,MPI_COMM_WORLD);//向1号进程发送准备好了的消息
	}
   	else{
  		MPI_Recv(&ready,1,MPI_INT,myid-1,myid-1,MPI_COMM_WORLD,&status);//接受前一个进程的消息
		PrintFile_info("match_result",T,myid);//打印运行参数信息
		if(myid!=numprocs-1)
			MPI_Send(&ready,1,MPI_INT,myid+1,myid,MPI_COMM_WORLD);//给下一个进程发送消息
	}

	printf("\n");
   	
	if((match=(int *)malloc(textlen*sizeof(int)))==NULL){//为是否匹配成功记录数组申请内存
		printf("no enough memory\n");
		exit(1);
	}
  
  	/*处理器0读入模式串，并记录运行参数*/
   	if(myid==0){ 
  		printf("processor num = %d \n",numprocs);//输出进程数量
    	printf("textlen = %d\n",textlen*numprocs); //输出字符串总长度

        GetFile("pattern.dat",&W,&patlen); //读取模式串
    	printf("patlen= %d\n",patlen); //输出模式串长度

    	if((nextval=(int *)malloc(patlen*sizeof(int)))==NULL){//为next数组分配空间
        	printf("no enough memory\n");
           	exit(1);
        }
		/*对模式串进行分析，对应于算法14.6步骤（1）*/
      	Next(W,patlen,nextval,&pped);//计算Next数组
        if(numprocs>1){
        	if (pped.pednum==1) 
           		nextlen_send = patlen;
            else 
        		nextlen_send = pped.pedlen*2;
        }
    }

	/*向各个处理器播送模式串的信息，对应于算法14.6步骤（2）*/
  	if(numprocs>1){
     	MPI_Bcast(&patlen, 1, MPI_INT, 0, MPI_COMM_WORLD);  //广播模式串长度
  		if(myid!=0)
    		if(((nextval=(int *)malloc(patlen*sizeof(int)))==NULL)//为模式串和next数组申请内存
				||((W=(char *)malloc(patlen*sizeof(char)))==NULL)){
           		printf("no enough memory\n");
            	exit(1);
            }

 	 	MPI_Barrier(MPI_COMM_WORLD);//同步
    	MPI_Bcast(&pped,3,MPI_INT,0,MPI_COMM_WORLD);  //广播其他进程所需的变量
    	MPI_Bcast(&nextlen_send,1,MPI_INT,0,MPI_COMM_WORLD);
    	MPI_Bcast(nextval,nextlen_send,MPI_INT,0,MPI_COMM_WORLD); 
    	MPI_Bcast(W,pped.pedlen,MPI_CHAR,0,MPI_COMM_WORLD);
   	}

    MPI_Barrier(MPI_COMM_WORLD);

   	/*调用修改过的KMP算法进行局部串匹配，对应于算法14.6步骤（3）*/
  	if(numprocs==1) {
  		kmp(T,W,textlen,patlen,nextval,&pped,1,0,match+patlen-1,&prefixlen);
   	}
   	else { 
    	if(myid!=0)
    		/*各个处理器分别根据部分串数据以及周期信息重构模式串*/
        	Rebuild_info(patlen,&pped,nextval,W); //重构模式串
    	if(myid!=numprocs-1)
  			kmp(T,W,textlen,patlen,nextval,&pped,0,0,match+patlen-1,&prefixlen);
		else
  			kmp(T,W,textlen,patlen,nextval,&pped,1,0,match+patlen-1,&prefixlen);

   		MPI_Barrier(MPI_COMM_WORLD);

		/*各个处理器进行段间匹配，对应于算法14.6步骤（4）*/
    	if(myid<numprocs-1) 
        	MPI_Send(&prefixlen,1,MPI_INT,myid+1,99,MPI_COMM_WORLD); //给当前进程的下一个进程发送前段长度

    	if(myid>0) 
    		MPI_Recv(&prefixlen,1,MPI_INT,myid-1,99,MPI_COMM_WORLD,&status); //除了第一个进程外接受前一个进程发送的前端信息

    	MPI_Barrier(MPI_COMM_WORLD);//同步

    	if((myid>0)&&(prefixlen!=0))  //段间匹配
   			kmp(T-prefixlen,W,prefixlen+patlen-1,patlen,nextval,&pped,1,prefixlen,match+patlen-1-prefixlen,&prefixlen); 

   		MPI_Barrier(MPI_COMM_WORLD);
   	}

	/*输出匹配结果*/
	if(myid==0){
		PrintFile_res("match_result",match+patlen-1,textlen-patlen+1,0,myid);//输出匹配结果
		if(numprocs>1)
			MPI_Send(&ready,1,MPI_INT,1,0,MPI_COMM_WORLD);
	}
   	else{
  		MPI_Recv(&ready,1,MPI_INT,myid-1,myid-1,MPI_COMM_WORLD,&status);
		PrintFile_res("match_result",match,textlen,myid*textlen-patlen+1,myid);//输出匹配结果
		if(myid!=numprocs-1)
			MPI_Send(&ready,1,MPI_INT,myid+1,myid,MPI_COMM_WORLD);
	}

	free(T);
    free(W);
    free(nextval);
    MPI_Finalize(); 
 } 
