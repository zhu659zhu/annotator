#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<omp.h>
#define TOTAL 486881

char *name[TOTAL+1];
int st[TOTAL+1];
int ed[TOTAL+1];
int st1[22][TOTAL+1];
int ed1[22][TOTAL+1];
int count[22];
double result[TOTAL];
char filename[22][20] = {"KN9204.1A","KN9204.1B","KN9204.1D","KN9204.2A","KN9204.2B","KN9204.2D","KN9204.3A","KN9204.3B","KN9204.3D","KN9204.4A","KN9204.4B","KN9204.4D","KN9204.5A","KN9204.5B","KN9204.5D","KN9204.6A","KN9204.6B","KN9204.6D","KN9204.7A","KN9204.7B","KN9204.7D","KN9204.Un"};

int readbed(){
    FILE* fp;
    fp = fopen("KN_9204.bed","r");
    int i,index =0 ;
    for(i=0;i<TOTAL+1;i++)
        name[i] = (char *)malloc(30*sizeof(char));
    while(feof(fp) == 0){
        fscanf(fp,"%s %d %d",name[index],&st[index],&ed[index]);
        index++;
    }
    fclose(fp);
}

void readKN(int i){
    int index = 0;
    char a[30];
    sprintf(a,"%s%s%s","TE_",filename[i],".txt");
    FILE *fp1;
    fp1 = fopen(a,"r");
    while(feof(fp1)==0){
        fscanf(fp1,"%d,%d",&st1[i][index],&ed1[i][index]);
        index++;
        //printf("index:%d\n",index);
    }
    count[i] = index - 1;
    fclose(fp1);
    printf("%d:(%s,%d)\n",i,a,count[i]);
}

double check(int i){
    int j;
    for(j=0;j<22;j++)
        if(strcmp(name[i],filename[j])==0)
            break;
    int findex = j;
    for(j=0;j<count[findex];j++){
        int mst = (st[i]>=st1[findex][j])?st[i]:st1[findex][j];
        int med = (ed[i]<=ed1[findex][j])?ed[i]:ed1[findex][j];
        if(med - mst >= 0){
            if(i==0)
            printf("%d %d %d %d %d\n",findex,st[i],ed[i],st1[findex][j],ed1[findex][j]);
            int allrange = ed[i] - st[i]+1;
            int truerange = med - mst+1;
            return (truerange*1.0/allrange);
        }
    }
    return 0;

}

int main(){
    readbed();
    int i;
    for(i=0;i<22;i++)
        readKN(i);
    for(i=0;i<TOTAL;i++){
        printf("%d\r",i);
        fflush(stdout);
        result[i] = check(i);
    }
    FILE* fpo = fopen("result.txt","w");
    for(i=0;i<TOTAL;i++)
        fprintf(fpo,"%f\n",result[i]);
    fclose(fpo);
    return 0;
}
