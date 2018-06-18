#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

int main(int argc,char *argv[])
{
	int strlen,pedlen,suffixlen,num,i,j;
   	char *string;
   	FILE *fp;

   	strlen=atoi(argv[1]);
   	pedlen=atoi(argv[2]);
   	srand(atoi(argv[3]));

   	string=(char*)malloc(strlen*sizeof(char));
   	if(string==NULL){
      	printf("malloc error\n");
      	exit(1);
   	}

   	for(i=0;i<pedlen;i++){
        num=rand()%26;
        string[i]='a'+num;
  	}
 
  	for(j=1;j<(int)(strlen/pedlen);j++)
       	strncpy(string+j*pedlen,string,pedlen);

   	if((suffixlen=strlen%pedlen)!=0)
   		strncpy(string+j*pedlen,string,suffixlen);

  	if((fp=fopen(argv[4],"w"))!=NULL){
     	fprintf(fp,"%s",string); 
     	fclose(fp);
   	}
   	else{
     	printf("file open error\n");
     	exit(1);
   	}

	return;
}
