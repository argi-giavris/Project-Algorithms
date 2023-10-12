#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "interface.h"

//counts number of vectors in file and number of coordinates
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

//save vectors to structs
void saveFile(char path[256], vector *vectors){
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


//write results to ouput file
void writeResults(char path[256], int queriesCounter, vector *queries, vector *vectors, int **lsh_results, double **distanceLSH, double **distanceTrue, double *tLSH, double *tTrue, int** radius_results, int N, int flag){
	int i;
	int j = 0;
	FILE *fp;

	fp = fopen(path, "w");
	
	printf("%d\n",queriesCounter);
	for(i=0; i<queriesCounter; i++){
		fprintf(fp, "Query: %s\n", queries[i].id);
		for(j = 0 ; j < N ; j++){
			if(lsh_results[i][j] != -1){
				fprintf(fp, "Nearest neighbor-%d: %s\n",j+1, vectors[lsh_results[i][j]].id);
				if(flag == 1){ //LSH
					fprintf(fp, "distanceLSH: %f\n", distanceLSH[i][j]);
				}else{
					fprintf(fp, "distanceCube: %f\n", distanceLSH[i][j]);
				}
				fprintf(fp, "distanceTrue: %f\n", distanceTrue[i][j]);
			}

			
		}

		if(flag == 1){
			fprintf(fp, "tLSH:  %f\n", tLSH[i]);
		}else{
			fprintf(fp, "tCube:  %f\n", tLSH[i]);
		}
		
		fprintf(fp, "tTrue: %f\n", tTrue[i]);
		fprintf(fp, "R-near neighbors\n");
		///fprintf(fp, "toso einai%d\n",sizeof(radius_results[i])/sizeof(radius_results[i][0]));
		j = 0;
		while(radius_results[i][j] != -1){
			fprintf(fp, "%s\n",vectors[radius_results[i][j]].id);
			j++;
		}
	}

	fclose(fp);

}


//find nearest neighbors of query given
void query_knn(int vectorsCounter, int coordsCounter, vector *vectors, vector query, double *distanceTrue, int N, int* trueResults){
	int i, j, min, min_pos;
	double distance;

	distance=0;
	int  maxPosition;
	double max;
	for(i=0; i<vectorsCounter; i++){
		distance = euclideanDistance(coordsCounter, vectors[i].coords, query.coords);
		if(i < N){//save N first distances to array
			distanceTrue[i] = distance;
			trueResults[i] = i;
			}
		else{
			max = distanceTrue[0];
			maxPosition = 0;
			for(j = 1; j < N; j++){
				if(max < distanceTrue[j]){
					max = distanceTrue[j];
					maxPosition = j;
				}
			}
			if(distance < max){
				distanceTrue[maxPosition] = distance;
				trueResults[maxPosition] = i;
			}

		}
	}
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