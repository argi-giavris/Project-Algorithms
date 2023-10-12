#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "interface.h"

int hammingDistance(int n1, int n2){ // calculate hamming distance between two integers
    int x = n1 ^ n2;
    int setBits = 0;
 
    while (x > 0) {
        setBits += x & 1;
        x >>= 1;
    }
 
    return setBits;
}

void cubeTrain(list_node** cube, vector * vectors, struct h_func * h, list_node** f, int vectorsCounter, int k, int coordsCounter,int w){
	int i = 0, z = 0, t = 0;
	srand(time(NULL));
	int cubePos = 0;
	list_node* temp;

	for(i=0; i<vectorsCounter; i++){
		cubePos = 0;
		for(t=0; t<k; t++){	//calculate h for euclidean space
			h[t].h_floor = 0;
			h[t].v = malloc(sizeof(float)* coordsCounter);
			normalDistr(h[t].v, coordsCounter);
			float pv = inner(vectors[i].coords,h[t].v,coordsCounter);
			h[t].h_floor =floor((pv+h[t].t)/w);
			
			temp = f[t];
			if(temp == NULL){
				f[t] = malloc(sizeof(list_node));
				f[t]->g = h[t].h_floor;
				f[t]->vec_pos = 2*(rand() / (RAND_MAX +1.0)); // toss coin
				f[t]->next = NULL;
				cubePos += f[t]->vec_pos*pow(2,k-1 -t);
				free(h[t].v);
				continue;
			}
			int flag = 0; // value doesnt exist
			list_node* prevTemp;
			while (temp != NULL){	// check if h_floor already exists
				if(temp->g == h[t].h_floor){
					flag = 1; // value exists
					cubePos += temp->vec_pos*pow(2,k-1 -t);
					break;
				}
				prevTemp = temp;
				temp = temp->next;
			
			}
			if(flag == 0){ // create a new node
				list_node* new_node = malloc(sizeof(list_node));
				new_node->next = NULL;
				new_node->g = h[t].h_floor;
				new_node->vec_pos = 2*(rand() / (RAND_MAX +1.0)); // toss coin
				prevTemp->next = new_node;
				cubePos += new_node->vec_pos*pow(2,k-1 -t);
			}

			free(h[t].v);
		}
		if(cube[cubePos] == NULL){	//insert vector into the cube
			cube[cubePos] = malloc(sizeof(list_node));
			cube[cubePos]->next = NULL;
			cube[cubePos]->vec_pos = i;
		}
		else{
			temp = cube[cubePos];
			while(temp->next!=NULL){ //find end of the list
				temp = temp->next;
			}
			temp->next = malloc(sizeof(list_node)); // add new vector at the end of the list
			temp->next->next = NULL;
			temp->next->vec_pos = i;
		}
		

	}

}

void cubeSearch(vector* vectors, vector query, struct h_func *h, list_node **cube, list_node** f, int vectorsCounter, int coordsCounter, int M, int k, int probes, int w, double* distanceCube, int *cubeResults, int N, int* radius_results, int radius, int* radiusCounter){
	//printf("Searching...\n");
	
	
	int j, h_floor,t;
	int *binaryString = malloc(k*sizeof(int));
	int* cubePos = malloc(probes*sizeof(int));
	list_node* temp;

	srand(time(0));
	for(t=0; t<k; t++){								
		h[t].h_floor = 0;
		h[t].v = malloc(sizeof(float)* coordsCounter);
		normalDistr(h[t].v, coordsCounter);
		float pv = inner(query.coords,h[t].v,coordsCounter);
		h[t].h_floor =floor((pv+h[t].t)/w);

		temp = f[t];
		int flag = 0;
		if(temp == NULL){
			binaryString[t] = 2*(rand() / (RAND_MAX +1.0));
		}
		while (temp != NULL){	// check if h_floor already exists
			if(temp->g == h[t].h_floor){
				flag = 1;
				binaryString[t] = temp->vec_pos;
				break;
			}
				temp = temp->next;
		}
		if(flag == 0){
			binaryString[t] = 2*(rand() / (RAND_MAX +1.0)); // toss coin
		}
		free(h[t].v);
	}

	int i, probesCounter = 1;
	for(i=0; i<probes; i++){
		cubePos[i]=0;
	}
	
	for(j=0; j<k; j++){
		cubePos[0] += binaryString[j]*pow(2,k-1 -j); 
	}
	
	int min = 1;
	int cubeSize = pow(2,k);
	int tempHamming;
	i = 0;
	j = 1;
	while(probesCounter < probes){
		if(cube[i] != NULL){
			tempHamming = hammingDistance(cubePos[0],i);
			if(tempHamming == min){
				cubePos[probesCounter] = i;
				probesCounter++;
			}
		}
		i++;
		if(i == cubeSize){
			i = 0;
			min++;
		}
	}
	
	
	float distance;
	double max;
	int flag = 0, counter = 0, maxPosition, z = 0;

	for(z = 0 ; z < probes ; z++){
		temp = cube[cubePos[z]];
		while(temp != NULL) {
			if(counter > M) break;
			//calculate distance between query and bucket vectors
			//store N nearest neighbors in sorted array
			
			distance = euclideanDistance(coordsCounter, vectors[temp->vec_pos].coords, query.coords);
			if(distance < radius){
				//check if already exists
				for(i= 0; i < (*radiusCounter) - 1; i++){
					if(radius_results[i] == temp->vec_pos){
						flag = 1;
						break;
					}
				}
				if(flag == 0){
					radius_results[(*radiusCounter)] = temp->vec_pos;
					(*radiusCounter)++;
				}
				
			}
			if(counter < N){
				distanceCube[counter] = distance;
				cubeResults[counter] = temp->vec_pos;
			}
			else{
				max = distanceCube[0];
				maxPosition = 0;
				for(i = 1; i < N; i++){
					if(max < distanceCube[i]){
						max = distanceCube[i];
						maxPosition = i;
					}
				}
				if(distance < max){
					distanceCube[maxPosition] = distance;
					cubeResults[maxPosition] = temp->vec_pos;
				}

			}
			counter++;
			temp = temp->next;
		}
		
	}

	double *newDistance;
	int *newResults;
	if(counter < N){//we didnt find N Nearest neighbors
		newDistance = malloc(counter *sizeof(double));
		newResults = malloc(counter * sizeof(int));

		for(int i = 0; i < counter; i++){
			newDistance[i] = distanceCube[i];
			newResults[i] = cubeResults[i];
		}
		sortArrays(newDistance, newResults, counter);

		for(i = 0; i < N; i++){
			if(i < counter){
				distanceCube[i] = newDistance[i];
				cubeResults[i] = newResults[i];
			}else{
				distanceCube[i] = -1.0;
				cubeResults[i] = -1;
			}
		}
	}else{
		sortArrays(distanceCube, cubeResults, N);
	}
	
	



	//frees

	free(binaryString);
	free(cubePos);
}


void cube(vector *vectors, vector *queries, int vectorsCounter, int queriesCounter, int coordsCounter, int k, int M, int w, int N, int radius, int probes, char* queryFile, int gotQueries, char* outputFile, int gotOutput){
	int i, j, z, m, hash_pos, TableSize;
	struct h_func *h; 
	///list_node ***HashTables;
	clock_t start, stop;


	h = malloc(k*sizeof(struct h_func));

	srand(time(0));
	for(j=0; j<k; j++){
		h[j].t = w*(rand() / (RAND_MAX +1.0));
	}
	
	list_node ** cube;
	int cubeSize = pow(2,k);
	cube = malloc(cubeSize*sizeof(list_node *));
	for(i=0; i<cubeSize; i++){
		cube[i] = NULL;
	}

	list_node **f = malloc(k*sizeof(list_node*));

	for(i = 0 ; i < k ; i++){
		f[i] = NULL;
	}


	cubeTrain(cube, vectors,h, f, vectorsCounter, k,coordsCounter,w);
	//done hashing vertexes

	//printf("Ready for queries\n");
	char answer[10];
	strcpy(answer,"yes");

	while((strcmp(answer,"yes") == 0) || (strcmp(answer,"y") == 0) || (strcmp(answer,"YES") == 0)){

		if(gotOutput == 0){
			printf("Insert query file path:\n");
			scanf("%s",queryFile);
		}
		gotOutput = 0;

		countFile(queryFile, &queriesCounter, &coordsCounter);
		queries = malloc(queriesCounter*sizeof(vector));
		i=0;
		for(i=0; i<queriesCounter; i++){
			queries[i].coords = malloc(coordsCounter*sizeof(int));
		}
		saveFile(queryFile, queries);

		double **distanceCube = malloc(queriesCounter*sizeof(double*));
		for(i=0; i<queriesCounter; i++){
			distanceCube[i] = malloc(N * sizeof(double));
		}

		double **distanceTrue = malloc(queriesCounter*sizeof(double*));
		for(i=0; i<queriesCounter; i++){
			distanceTrue[i] = malloc(N * sizeof(double));
		}

		int **cubeResults = malloc(queriesCounter*sizeof(int*));

		for(i=0; i<queriesCounter; i++){
			cubeResults[i] = malloc(N * sizeof(int));
		}

		int **trueResults = malloc(queriesCounter*sizeof(int*));

		for(i=0; i<queriesCounter; i++){
			trueResults[i] = malloc(N * sizeof(int));
		}

		double *tCube = malloc(queriesCounter*sizeof(double));

		double *tTrue = malloc(queriesCounter*sizeof(double));


		int **radius_results = malloc(queriesCounter*sizeof(int*));

		clock_t start, stop;
		for(i=0; i<queriesCounter; i++){
			start = clock();
			query_knn(vectorsCounter, coordsCounter, vectors, queries[i], distanceTrue[i], N, trueResults[i]);
			sortArrays(distanceTrue[i],trueResults[i], N);
			stop = clock();
			tTrue[i] = (double)(stop-start) / CLOCKS_PER_SEC;
		}


		int *temp_radius_results = malloc(vectorsCounter*sizeof(int));
		int radiusCounter = 0;
		for(i = 0; i < queriesCounter; i++){
			radiusCounter = 0;
			start = clock();
			cubeSearch(vectors, queries[i], h, cube, f, vectorsCounter, coordsCounter, M, k, probes, w, distanceCube[i], cubeResults[i], N, temp_radius_results, radius, &radiusCounter);
			stop = clock();
			tCube[i] = (double)(stop-start) / CLOCKS_PER_SEC;

			radius_results[i] = malloc((radiusCounter+1)*sizeof(int));
			for(j = 0 ; j < radiusCounter ; j++){
				radius_results[i][j] = temp_radius_results[j];
			}
			radius_results[i][radiusCounter]=-1;
		}

		int flag = 0; // 0 for cube and 1 for LSH
		writeResults(outputFile, queriesCounter,queries, vectors, cubeResults, distanceCube, distanceTrue, tCube, tTrue,radius_results, N, flag);

		for(i = 0; i < queriesCounter; i++){
			free(queries[i].coords);
			free(distanceCube[i]);
			free(distanceTrue[i]);
			free(cubeResults[i]);
			free(trueResults[i]);
			free(radius_results[i]);
		}

		free(queries);
		free(distanceCube);
		free(distanceTrue);
		free(cubeResults);
		free(trueResults);
		free(radius_results);
		free(temp_radius_results);

		free(tCube);
		free(tTrue);

		printf("Do you want to insert a new query file:y\\n\n");
		scanf("%s",answer);

	}

	list_node *temp;
	for(i=0; i<cubeSize; i++){
		
		while(cube[i] != NULL){
			temp = cube[i];
			cube[i] = cube[i]->next;
			free(temp);
		}
	}
	free(cube);

	for(i=0; i<k; i++){
		
		while(f[i] != NULL){
			temp = f[i];
			f[i] = f[i]->next;
			free(temp);
		}
	}
	
	free(f);
	free(h);

	
	
}


void main (int argc, char *argv[]){
	int i, k = 14 ,M = 10, probes = 2,  vectorsCounter, queriesCounter, coordsCounter, N = 1,  gotInput = 0, gotQueries = 0, gotOutput = 0;
	char ch, inputFile[256], queryFile[256], outputFile[256];
	double radius = 10000;
	vector *vectors, *queries;
	int argCount = 0;

    for(argCount = 0 ; argCount < argc ; argCount++){ //reading arguments from commadn line
        if(strcmp(argv[argCount],"-i") == 0){
            strcpy(inputFile, argv[argCount+1]);
			gotInput = 1;
        }
        else if (strcmp(argv[argCount],"-q") == 0){
            strcpy(queryFile, argv[argCount+1]);
			gotQueries = 1;
        }
        else if (strcmp(argv[argCount],"-k") == 0){
		    k = atoi(argv[argCount+1]);
        }
        else if (strcmp(argv[argCount],"-M") == 0){
		    M = atoi(argv[argCount+1]);
        }  
        else if (strcmp(argv[argCount],"-o") == 0){
		    strcpy(outputFile, argv[argCount+1]);
			gotOutput = 1;
        }
        else if (strcmp(argv[argCount],"-N") == 0){
		    N = atoi(argv[argCount+1]);
        }
        else if (strcmp(argv[argCount],"-R") == 0){
		    radius = atof(argv[argCount+1]);
        }
        else if (strcmp(argv[argCount],"-probes") == 0){
		    probes = atoi(argv[argCount+1]);
        }
    }
	countFile(inputFile, &vectorsCounter, &coordsCounter);
	vectors = malloc(vectorsCounter*sizeof(vector));
	i=0;
	for(i=0; i<vectorsCounter; i++){
		vectors[i].coords = malloc(coordsCounter*sizeof(int));
	}
	saveFile(inputFile, vectors);

	int w = 1000; //w is averagedistance of vectors * 4;

	cube(vectors,queries, vectorsCounter, queriesCounter, coordsCounter, k, M, w, N, radius, probes, queryFile,gotQueries,outputFile,gotOutput);

	//frees

	for(i = 0; i < vectorsCounter; i++){
		free(vectors[i].coords);
	}
	free(vectors);


	return; 
}