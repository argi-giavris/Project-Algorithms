#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "interface.h"

void lshTrain(vector *vectors, struct h_func **h, list_node ***HashTables, int vectorsCounter, int coordsCounter, int M, int k, int L, int w, int TableSize){
	int i, j, z, t, hash_pos, *a;
	int g;
	float f;
	srand(time(NULL));

	//save vectors to L hash tables
	for(i=0; i<vectorsCounter; i++){
		for(z=0; z<L; z++){
			g=0;
			for(t=0; t<k; t++){								
				h[z][t].h_floor = 0;
				h[z][t].v = malloc(sizeof(float)* coordsCounter);
				normalDistr(h[z][t].v, coordsCounter);
				//printf("TEST %d ",vectors[i].coordsCounter);
				float pv = inner(vectors[i].coords,h[z][t].v,coordsCounter);//esoteriko ginomeno
				h[z][t].h_floor =floor((pv+h[z][t].t)/w);
				int r;
				r = rand();
				g = eucmod(eucmod(g,M)+eucmod(h[z][t].r*h[z][t].h_floor,M),M);
				//g = g+h[z][t].r*h[z][t].h_floor;
				free(h[z][t].v);
			}
			//int ID = eucmod(g,M);
			int ID = g;
			g = eucmod(ID,TableSize);
			
			//printf("G = %d\n", g);
			hash(HashTables[z], g, ID, i);
		}
		//printf("Done hashing vector %d", i);
	}
}

void lshSearch(vector* vectors, vector query, struct h_func **h, list_node ***HashTables, int vectorsCounter, int coordsCounter, int M, int k, int L, int w, int TableSize, double* distanceLSH, int *lsh_results, int N, int* radius_results, int radius, int* radiusCounter){

	//printf("\nSearching..\n");
	int i, j, z, t, hash_pos, *a;
	int g;
	float f;
	srand(time(NULL));

	int flag = 0; 
	float distance;
	float minDistance;
	int minPosition;
	*radiusCounter = 0;
	int counter = 0;

	
	for(z=0; z<L; z++){//for every hash table
		g=0;
		for(t=0; t<k; t++){								
			h[z][t].h_floor = 0;
			h[z][t].v = malloc(sizeof(float)* coordsCounter);
			normalDistr(h[z][t].v, coordsCounter);
			//printf("TEST %d ",vectors[i].coords);
			float pv = inner(query.coords,h[z][t].v,coordsCounter);
			h[z][t].h_floor =floor((pv+h[z][t].t)/w);
			int r;
			r = rand();
			//g = g+h[z][t].r*h[z][t].h_floor;
			g = eucmod(eucmod(g,M)+eucmod(h[z][t].r*h[z][t].h_floor,M),M);
			free(h[z][t].v);
		}
		int ID = g;
		g = eucmod(ID,TableSize);
		//hash(HashTables[z], g, g, i);
		list_node *temp;
		temp = HashTables[z][g];
		double max;
		int maxPosition;
		

		while(temp != NULL) {//traverse list bucket
			if(temp->g == ID){
				//calculate distance between query and bucket vectors
				//store N nearest neighbors in sorted array
				flag = 0;
				distance = euclideanDistance(coordsCounter, vectors[temp->vec_pos].coords, query.coords);
				if(distance <= radius){
					//check if already exists
					for(i= 0; i < (*radiusCounter); i++){
						if(radius_results[i] == temp->vec_pos){
							flag = 1;
						}
					}
					if(flag == 0){//if doesnt exists add it to array
						radius_results[(*radiusCounter)] = temp->vec_pos;
						(*radiusCounter)++;
					}
					
				}

				if(counter < N){
					distanceLSH[counter] = distance;
					lsh_results[counter] = temp->vec_pos;
				}
				else{
					max = distanceLSH[0];
					maxPosition = 0;
					for(i = 1; i < N; i++){
						if(max <= distanceLSH[i]){
							max = distanceLSH[i];
							maxPosition = i;
						}
					}
					if(distance <= max){
						distanceLSH[maxPosition] = distance;
						lsh_results[maxPosition] = temp->vec_pos;
					}

				}
				counter++;
			}
			
			temp = temp->next;
		}
		
	}
	
	sortArrays(distanceLSH, lsh_results, N);

}

void lsh(vector *vectors, vector *queries, int vectorsCounter, int queriesCounter, int coordsCounter, int k, int L, int w, int N, int radius, char* queryFile, int gotQueries, char* outputFile, int gotOutput){
	//first we save data set vertexes to hash tables and then we search with queries
	
	
	int i, j, z, M, hash_pos, TableSize;
	struct h_func **h; 
	list_node ***HashTables;
	clock_t start, stop;

	h = malloc(L*sizeof(struct h_func *));
	for(i=0; i<L; i++){	
		h[i] = malloc(k*sizeof(struct h_func));
	}

	srand(time(0));
	for(i=0; i<L; i++){
		for(j=0; j<k; j++){
			h[i][j].t = w*(rand() / (RAND_MAX +1.0));
		}
	}

	for(j=0; j<k; j++){
		h[0][j].r = rand();
		for(i = 1 ; i < L ; i++){
			h[i][j].r = h[0][j].r;
		}
	}
	
	TableSize = vectorsCounter/4;	
	HashTables = malloc(L*sizeof(struct list_node **));
	for(i=0; i<L; i++){
		HashTables[i] = malloc(TableSize*sizeof(struct list_node *));
	}
	for(i=0; i<L; i++){	
		for(j=0; j<TableSize; j++){
			HashTables[i][j] = NULL;
		}
	}

	M = pow(2, 32/k);

	lshTrain(vectors, h, HashTables, vectorsCounter, coordsCounter, M, k, L, w, TableSize);
	//done hashing vertexes

	printf("Ready for queries\n");


	char answer[10];
	strcpy(answer,"yes");


	//lsh_search queries until user exit
	while((strcmp(answer,"yes") == 0) || (strcmp(answer,"y") == 0) || (strcmp(answer,"YES") == 0)){

		//if not given from terminal
		if(gotOutput == 0){
			printf("Insert query file path:\n");
			scanf("%s",queryFile);
		}
		gotOutput = 0;

		countFile(queryFile, &queriesCounter, &coordsCounter);

		//queries array vector
		queries = malloc(queriesCounter*sizeof(vector));
		i=0;
		for(i=0; i<queriesCounter; i++){
			queries[i].coords = malloc(coordsCounter*sizeof(int));
		}
		printf("saving queries\n");
		saveFile(queryFile, queries);
		
		double **distanceLSH = malloc(queriesCounter*sizeof(double*));
		for(i=0; i<queriesCounter; i++){
			distanceLSH[i] = malloc(N * sizeof(double));
			distanceLSH[i][0] = -1.0;
		}

		double **distanceTrue = malloc(queriesCounter*sizeof(double*));
		for(i=0; i<queriesCounter; i++){
			distanceTrue[i] = malloc(N * sizeof(double));
		}

		//positions of estimated NN
		int **lsh_results = malloc(queriesCounter*sizeof(int*));

		for(i=0; i<queriesCounter; i++){
			lsh_results[i] = malloc(N * sizeof(int));
		}

		//positions of true NN
		int **trueResults = malloc(queriesCounter*sizeof(int*));

		for(i=0; i<queriesCounter; i++){
			trueResults[i] = malloc(N * sizeof(int));
		}

		double *tLSH = malloc(queriesCounter*sizeof(double));

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

		//search NN for each query
		for(i = 0; i < queriesCounter; i++){
			radiusCounter = 0;
			start = clock();
			lshSearch(vectors, queries[i], h, HashTables, vectorsCounter, coordsCounter, M, k, L, w, TableSize, distanceLSH[i], lsh_results[i], N, temp_radius_results, radius, &radiusCounter);
			stop = clock();
			tLSH[i] = (double)(stop-start) / CLOCKS_PER_SEC;

			radius_results[i] = malloc((radiusCounter+1)*sizeof(int));
			for(j = 0 ; j < radiusCounter ; j++){
				radius_results[i][j] = temp_radius_results[j];
			}
			radius_results[i][radiusCounter]=-1;
			//printf("Query with id %d-> Nearest neighbor %d with distance %f\n", queries[i].id, vectors[*lsh_results[i]].id, distanceLSH[i]);
			
		}
		int flag = 1;

		//write results to output file
		writeResults(outputFile, queriesCounter,queries, vectors, lsh_results, distanceLSH, distanceTrue, tLSH, tTrue,trueResults, N, 0);

		for(i=0; i< queriesCounter; i++){
			free(queries[i].coords);
			free(queries[i].id);
			free(distanceLSH[i]);
			free(distanceTrue[i]);
			free(lsh_results[i]);
			free(trueResults[i]);
			free(radius_results[i]);
		}
		free(queries);
		free(distanceLSH);
		free(distanceTrue);
		free(lsh_results);
		free(trueResults);
		free(radius_results);

		free(tLSH);
		free(tTrue);

		free(temp_radius_results);
		printf("Do you want to insert a new query file:y\\n\n");
		scanf("%s",answer);
	}
	



	//free data set vector hashtable
	list_node *temp;

	for(i = 0; i <L; i++){
		free(h[i]);
	}
	free(h);

	for(i = 0; i < L; i++){
		for(j = 0; j < TableSize; j++){
			while(HashTables[i][j] != NULL){
				temp = HashTables[i][j];
				HashTables[i][j] = HashTables[i][j]->next;
				free(temp);
			}
			
		}
		free(HashTables[i]);

	}
	free(HashTables);


	
}