#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "interfaceCluster.h"

int hammingDistance(int n1, int n2){
    int x = n1 ^ n2;
    int setBits = 0;
 
    while (x > 0) {
        setBits += x & 1;
        x >>= 1;
    }
 
    return setBits;
}

void countFile(char path[256], int *vectorsCounter, int *coordsCounter){
    FILE* fp = fopen(path, "r");
	(*coordsCounter)=0;
	(*vectorsCounter)=0;
    char line[8000];
    while (fgets(line, 8000, fp)){
		///(*vectorsCounter)++;
		char* tmp = strdup(line);
		const char* tok;
		tok = strtok(line,"\t\n");
		if(tok != NULL) (*vectorsCounter)++;
		while((tok = strtok(NULL,"\t\n")) != NULL){
			if((*vectorsCounter) == 1){
				(*coordsCounter)++;
			}
		}
	}
	printf("v = %d c = %d\n",*vectorsCounter,*coordsCounter);

/*	char ch;
	FILE *fp;
	fp = fopen(path,"r");
	(*coordsCounter)=-1;
	(*vectorsCounter)=0;
	while(1){
		ch = fgetc(fp);
		if(ch==EOF){
			break;
		}
        else if(ch=='\n'){
			(*vectorsCounter)++;
		}
        else if((ch=='\t') && ((*vectorsCounter)==0)){
			(*coordsCounter)++;
		}
	}
	printf("vectorsCounter = %d, coordsCounter = %d\n",*vectorsCounter, *coordsCounter);
	fclose(fp);*/
}

void saveFile(char path[256], vectorCluster *vectors){
	char ch;
	FILE *fp;
	fp = fopen(path,"r");
	int digitNumber=0;
	int coordsCounter;
	int vectorsCounter=0;
	int flag = 0;

    char line[8000];
    while (fgets(line, 8000, fp)){
		coordsCounter = 0;
		char* tmp = strdup(line);
		const char* tok;
		tok = strtok(line,"\t\n");
		if(tok != NULL){
			vectors[vectorsCounter].id= malloc(sizeof(char)*(strlen(tok)+1));
			strcpy(vectors[vectorsCounter].id,tok);
			while((tok = strtok(NULL,"\t\n")) != NULL){
				vectors[vectorsCounter].coords[coordsCounter]=atof(tok);
				coordsCounter++;
			}
		}
		free(tmp);
		vectorsCounter++;
	}

/*	while(1){
		ch = fgetc(fp);
		if(ch==EOF){
			break;
		}
        else if(ch=='\n'){
			printf("changed line %c\n",ch);
			vectorsCounter++;
			coordsCounter=-1;
		}
        else if(ch=='\t'){
			digitNumber=0;
			if(coordsCounter!=0){
				//printf("%s  ",num);
				vectors[vectorsCounter].coords[coordsCounter]=atof(num);
				coordsCounter++;
			}
            else{
				//printf("id  =      %s \n",num);
				//printf("strlen = %ld\n",strlen(num));
				vectors[vectorsCounter].id= malloc(sizeof(char)*(strlen(num)+1));
				strcpy(vectors[vectorsCounter].id,num);
				coordsCounter++;
			}
		}
        else{
			num[digitNumber++]=ch;
			num[digitNumber]='\0';
		}	
	}*/
	fclose(fp);
}

void writeResults(char path[256], int method, int numOfClusters, vectorCluster *vectors, int vectorsCounter, int coordsCounter, vectorCluster* centroids, double clusterTime, int complete, float* silhouettes, float* stotal){
	int i;
	int j = 0;
	FILE *fp;

	fp = fopen(path, "w");
	if(method == 0)
		fprintf(fp, "Algorithm: Lloyds\n");
	else if(method == 1)
		fprintf(fp, "Algorithm: Range Search LSH\n");
	else if(method == 2)
		fprintf(fp, "Algorithm:  Range Search Hypercube\n");
	else
		fprintf(fp, "Algorithm:  Range Search Hypercube\n");

	int vectorsInCluster;
	for(i=0; i<numOfClusters; i++){
		vectorsInCluster = 0;
		if(complete == 0){
			for(j = 0 ; j < vectorsCounter ; j++){
				if(vectors[j].nearest == i)	vectorsInCluster++;
			}
			fprintf(fp, "CLUSTER-%d {size: %d, centroid:", i+1,vectorsInCluster);
			for(j = 0 ; j < coordsCounter ; j++){
				fprintf(fp, " %f", centroids[i].coords[j]);
			}
			fprintf(fp, "}\n");
		}
		else{
			for(j = 0 ; j < vectorsCounter ; j++){
				if(vectors[j].nearest == i)	vectorsInCluster++;
			}
			fprintf(fp, "CLUSTER-%d {size: %d, centroid:", i+1,vectorsInCluster);
			for(j = 0 ; j < coordsCounter ; j++){
				fprintf(fp, " %f", centroids[i].coords[j]);
			}
			for(j = 0 ; j < vectorsCounter ; j++){
				if(vectors[j].nearest == i)	fprintf(fp, ", %s", vectors[j].id);
			}			
			fprintf(fp, "}\n");	
		}
	}
	fprintf(fp, "clustering_time: %f\n",clusterTime);
	float averageSilhouette = 0.0;
	if (silhouettes != NULL){
		fprintf(fp, "Silhouette: [");
		for(j = 0 ; j < numOfClusters ; j++){
			averageSilhouette += silhouettes[j];
			fprintf(fp,"%f,",silhouettes[j]);
		}
		fprintf(fp," %f]\n",*stotal);
	}
	fclose(fp);
}

int eucmod(int a, int base){    
  return (a < 0 ? ((a % base) + base) % base : a % base);
}


float inner(float *vertexP, float *vertexV, int n)
{
    float product = 0;
    int i;

    //printf("n = %d\n", n);
    for (i = 0; i < n; i++)
    {
		//printf("%d %f\n",vertexP[i], vertexV[i]);
        //printf("P = %2d, A[%d] = %d, B[%d] = %d\n", product, i, A[i], i, B[i]);
        product = product + vertexP[i] * vertexV[i];
		
    }

	
    return product;
}

void uniformDistr(float* vertex, int d, double a, double b){

    int i;
    for(i=0; i < d; i++){
        vertex[i] = rand() / (RAND_MAX + 1.0) * (b - a) + a;
    }
	return;
}

void normalDistr(float* vertex, int d){
    double sigma = 1;
    double m = 0;
    
    double v1, v2, temp;

    int i;
    for(i=0; i < d; i++){
        v1 = ((double)(rand()) + 1. )/( (double)(RAND_MAX) + 1);
        v2 = ((double)(rand()) + 1. )/( (double)(RAND_MAX) + 1);
        temp = cos(2*3.14*v2)*sqrt(-2.*log(v1));
        vertex[i] = temp * (sigma + m);
        //printf("%f", vertex[i]);
    }
	return;
}

float euclideanDistance(int coordsCounter,float* coords1, float* coords2){
	float dist = 0;
	//printf("Coords %d ", coordsCounter);
	for(int i=0 ; i<coordsCounter ; i++){
		
		dist+=((coords1[i])-(coords2[i])) * ((coords1[i])-(coords2[i]));
	}
	//printf("%f\n",dist);
	return sqrtf((dist));
}




void swap1(double* xp, double* yp)
{
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void swap2(int *xp, int *yp){

	int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void sortArrays(double *dist, int *result, int n){
	//int n = sizeof(dist) / sizeof(dist[0]);
	int i, j, minDistPos;
	

	for (i = 0; i < n - 1; i++){

		minDistPos = i;
		for(j = i + 1; j < n; j++){
			if(dist[j] < dist[minDistPos]){
				minDistPos = j;
			}
		}

		swap1(&dist[minDistPos], &dist[i]);
		swap2(&result[minDistPos], &result[i]);
	}

	
}

void hash(list_node **Hash, int pos, int ID, int i){
	list_node *cur;
	

	cur=Hash[pos];
	
	if(cur==NULL){
		Hash[pos]=malloc(sizeof(list_node));
		Hash[pos]->next=NULL;
		Hash[pos]->g=ID;
		Hash[pos]->vec_pos=i;
	}else{
		while(cur->next!=NULL){
			cur= cur->next;
		}
		cur->next= malloc(sizeof(list_node));
		cur->next->g=ID;
		cur->next->next=NULL;
		cur->next->vec_pos=i;
	}
}