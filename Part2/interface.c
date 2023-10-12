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
		free(tmp);
	}
	printf("v = %d c = %d\n",*vectorsCounter,*coordsCounter);
	fclose(fp);

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
			vectors[vectorsCounter].id = malloc(sizeof(char)*(strlen(tok)+1));
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
void writeResults(char path[256], int queriesCounter, vector *queries, vector *vectors, int **lsh_results, double **distanceLSH, double **distanceTrue, double *tLSH, double *tTrue, int** trueResults, int N, int flag){
	int i;
	int j = 0;
	FILE *fp;
	double tTrueAverage = 0.0;
	double tApproximateAverage = 0.0;
	double MAF = 0.0;
	fp = fopen(path, "w");
	
	printf("%d\n",queriesCounter);
	for(i=0; i<queriesCounter; i++){
		fprintf(fp, "Query: %s\n", queries[i].id);
		if(flag == 0){ //LSH
			fprintf(fp, "Algorithm: LSH_Vector\n");
		}
		else if(flag == 1){ //HyperCube
			fprintf(fp, "Algorithm: HyperCube\n");
		}
		else if(flag == 2){ //LSH_Frechet_Discrete
			fprintf(fp, "Algorithm: LSH_Frechet_Discrete\n");
		}
		else if(flag == 3){ //LSH_Frechet_Continuous
			fprintf(fp, "Algorithm: LSH_Frechet_Continuous\n");
		}			
		for(j = 0 ; j < N ; j++){
			//if(lsh_results[i][j] != -1){
				if(distanceLSH[i][j] != -1.0)
					fprintf(fp, "Aproximate Nearest neighbor: %s\n", vectors[lsh_results[i][j]].id);
				else
					fprintf(fp, "Aproximate Nearest neighbor: -\n");
				//printf("na ta %d\n",trueResults[i][j]);
				fprintf(fp, "True Nearest neighbor: %s\n", vectors[trueResults[i][j]].id);
				if(distanceLSH[i][j] != -1.0)
					fprintf(fp, "distanceApproximate: %f\n", distanceLSH[i][j]);
				else
					fprintf(fp, "distanceApproximate: -\n");
				fprintf(fp, "distanceTrue: %f\n", distanceTrue[i][j]);

				tApproximateAverage = tApproximateAverage + tLSH[i];
				tTrueAverage = tTrueAverage + tTrue[i];
				if((distanceLSH[i][j] > 0.0) && (distanceTrue[i][j] == 0.0)){
					MAF = -1.0; //INFINITY
				}
				else if((MAF != -1.0) && (distanceLSH[i][j] == 0.0) && (distanceTrue[i][j] == 0.0)){
					if(MAF < 1.0) MAF = 1.0;
				}
				else if((MAF != -1.0) && (distanceLSH[i][j]/distanceTrue[i][j] > MAF)){
						MAF = distanceLSH[i][j]/distanceTrue[i][j];
				}

			//}

			
		}
		
	}
	fprintf(fp, "tApproximateAverage: %f\n", tApproximateAverage/queriesCounter);
	fprintf(fp, "tTrueAverage: %f\n", tTrueAverage/queriesCounter);
	if(MAF == -1.0)
		fprintf(fp, "MAF: INF\n");
	else
		fprintf(fp, "MAF: %f\n", MAF);

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
			//printf("true results = %d\n",trueResults[i]);
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

void uniformDistr(float* vertex, int d, double a, double b){

    int i;
    for(i=0; i < d; i++){
        vertex[i] = rand() / (RAND_MAX + 1.0) * (b - a) + a;
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