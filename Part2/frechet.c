#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "interface.h"

vector* snap(vector* vectors, double delta,int vectorsCounter, int coordsCounter, float* t){
	int i, j;
	vector* vectorsSnapped = malloc(vectorsCounter*sizeof(vector));
	for(i=0; i<vectorsCounter; i++){
		vectorsSnapped[i].coords = malloc((2*coordsCounter)*sizeof(float));
	}
	for(i = 0 ; i < vectorsCounter ; i++){
		for(j = 0 ; j < 2*coordsCounter ; j+=2){
			vectorsSnapped[i].coords[j] = (j/2);
		}
	}

	for(i = 0 ; i < vectorsCounter ; i++){
		vectorsSnapped[i].id =  malloc(sizeof(char)*(strlen(vectors[i].id)+1));
		strcpy(vectorsSnapped[i].id,vectors[i].id);
		float* x = vectors[i].coords;
		for(j = 0 ; j < 2*coordsCounter ; j+=2){
			//printf("t = %f\n",t[j/2]);
			vectorsSnapped[i].coords[j] = floor((vectorsSnapped[i].coords[j]-t[j/2])/delta +1/2) * delta + t[j/2]; //x
			vectorsSnapped[i].coords[j+1] = floor((x[j/2]-t[j/2])/delta +1/2) * delta + t[j/2]; //y
		}
	}
	return vectorsSnapped;
}

vector* snapContinuous(vector* vectors, double delta,int vectorsCounter, int coordsCounter, float* t, int* filteringCounter){
	int i, j;
	vector* vectorsSnapped = malloc(vectorsCounter*sizeof(vector));
	for(i=0; i<vectorsCounter; i++){
		vectorsSnapped[i].coords = malloc((coordsCounter)*sizeof(float));
	}

	for(i = 0 ; i < vectorsCounter ; i++){
		vectorsSnapped[i].id =  malloc(sizeof(char)*(strlen(vectors[i].id)+1));
		strcpy(vectorsSnapped[i].id,vectors[i].id);
		float* x = vectors[i].coords;
		for(j = 0 ; j < coordsCounter - filteringCounter[i] ; j++){
			//printf("t = %f\n",t[j/2]);
			vectorsSnapped[i].coords[j] = floor((x[j]+t[j])/delta) * delta; //x
		}
	}
	return vectorsSnapped;
}

void frechetTrain(vector *vectors, struct h_func **h, list_node ***HashTables, int vectorsCounter, int coordsCounter, int M, int k, int L, int w, int TableSize, double delta){
	int i, j, z, t, hash_pos, *a;
	int g;
	float f;
	srand(time(NULL));

	vector** vectorsSnapped = malloc(L*sizeof(vector*));
	for(z = 0 ; z < L ; z++){
		vectorsSnapped[z] = snap(vectors, delta, vectorsCounter, coordsCounter,h[z][0].tFrechet);
		int j;
		for(i = 0 ; i < vectorsCounter ; i++){
			for(j = 2 ; j < coordsCounter*2 ; j+=2){
				if((vectorsSnapped[z][i].coords[j] == vectorsSnapped[z][i].coords[j-2]) && (vectorsSnapped[z][i].coords[j+1] == vectorsSnapped[z][i].coords[j-1])){
					for(int x = j ; x < coordsCounter*2 ; x+=2){
						vectorsSnapped[z][i].coords[x] = vectorsSnapped[z][i].coords[x+2];
						vectorsSnapped[z][i].coords[x+1] = vectorsSnapped[z][i].coords[x+3];
					}
					vectorsSnapped[z][i].coords[(coordsCounter*2)-2] = 2147483647; //INT_MAX
					vectorsSnapped[z][i].coords[(coordsCounter*2)-1] = 2147483647; //INT_MAX
				}
			}
		}
	}

	//save vectors to L hash tables
	for(i=0; i<vectorsCounter; i++){
		for(z=0; z<L; z++){
			g=0;
			for(t=0; t<k; t++){								
				h[z][t].h_floor = 0;
				h[z][t].v = malloc(sizeof(float)* coordsCounter*2);
				normalDistr(h[z][t].v, coordsCounter*2);
				//printf("TEST %d ",vectors[i].coordsCounter);
				float pv = inner(vectorsSnapped[z][i].coords,h[z][t].v,coordsCounter*2);//esoteriko ginomeno
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
	for(i = 0 ; i < L ; i++){
		for(j = 0 ; j < vectorsCounter ; j++){
			free(vectorsSnapped[i][j].coords);
			free(vectorsSnapped[i][j].id);
		}
		free(vectorsSnapped[i]);
	}
	free(vectorsSnapped);
}

void frechetContinuousTrain(vector *vectors, struct h_func *h, list_node **HashTables, int vectorsCounter, int coordsCounter, int M, int k, int L, int w, int TableSize, double delta){
	int i, j, z, t, hash_pos, *a, *filteringCounter;
	int g;
	float f, epsilon = 0.5;
	srand(time(NULL));

	filteringCounter = malloc(vectorsCounter*sizeof(int));
	for(i = 0 ; i < vectorsCounter ; i++){ //Filtering
		filteringCounter[i] = 0;
		for(j = 0 ; j < coordsCounter ; j+=3){
			if((fabs(vectors[i].coords[j] - vectors[i].coords[j+1]) <= epsilon) && (fabs(vectors[i].coords[j+1] - vectors[i].coords[j+2]) <= epsilon)){
				vectors[i].coords[j+1] = -1.0;
				filteringCounter[i]++;
			}
		}
		for(j = 0 ; j < coordsCounter ; j++){
			if(vectors[i].coords[j] == -1.0){
				for(int x = j ; x < coordsCounter ; x++){
					vectors[i].coords[j] = vectors[i].coords[j+1];
				}
				vectors[i].coords[coordsCounter-1] = -1.0;
			}
		}
	}

	vector* vectorsSnapped;
	vectorsSnapped = snapContinuous(vectors, delta, vectorsCounter, coordsCounter,h[0].tFrechet,filteringCounter);
	for(i = 0 ; i < vectorsCounter ; i++){ //Min / Max
		for(j = 0 ; j < coordsCounter ; j+=3){
			if(((vectorsSnapped[i].coords[j+1] > vectorsSnapped[i].coords[j]) && (vectorsSnapped[i].coords[j+1] < vectorsSnapped[i].coords[j+2])) || (vectorsSnapped[i].coords[j+1] == vectorsSnapped[i].coords[j]) || (vectorsSnapped[i].coords[j+1] == vectorsSnapped[i].coords[j+2])){
				vectorsSnapped[i].coords[j+1] = -1.0;
				filteringCounter[i]++;
			}
		}
		for(j = 0 ; j < coordsCounter ; j++){
			if(vectorsSnapped[i].coords[j] == -1.0){
				for(int x = j ; x < coordsCounter ; x++){
					vectorsSnapped[i].coords[j] = vectorsSnapped[i].coords[j+1];
				}
				vectorsSnapped[i].coords[coordsCounter-1] = -1.0;
			}
		}
	}

	for(i = 0 ; i < vectorsCounter ; i++){
		for(j = coordsCounter - filteringCounter[i] ; j < coordsCounter ; j++){
			vectorsSnapped[i].coords[j] = 2147483647; //INT_MAX
		}
	}

	//save vectors to L hash tables
	for(i=0; i<vectorsCounter; i++){
		g=0;
		for(t=0; t<k; t++){								
			h[t].h_floor = 0;
			h[t].v = malloc(sizeof(float)* coordsCounter*2);
			normalDistr(h[t].v, coordsCounter*2);
			//printf("TEST %d ",vectors[i].coordsCounter);
			float pv = inner(vectorsSnapped[i].coords,h[t].v,coordsCounter);//esoteriko ginomeno
			h[t].h_floor =floor((pv+h[t].t)/w);
			int r;
			r = rand();
			g = eucmod(eucmod(g,M)+eucmod(h[t].r*h[t].h_floor,M),M);
			//g = g+h[z][t].r*h[z][t].h_floor;
			free(h[t].v);
		}
		//int ID = eucmod(g,M);
		int ID = g;
		g = eucmod(ID,TableSize);
		
		//printf("G = %d\n", g);
		hash(HashTables, g, ID, i);
		//printf("Done hashing vector %d", i);
	}
}

float maxFloat(float num1, float num2){
    return (num1 > num2 ) ? num1 : num2;
}

float minFloat(float num1, float num2, float num3) {
	if(num1 < num2 && num1 < num3)
	{
		return num1;
	}
	else if(num2 < num3)
	{
		return num2;
	}
	else
	{
		return num3;
	}	
	return 0;
}

float recursiveC(float** ca,vector a, vector b, int i, int j){
	//printf("%d %d\n",i,j);
	if(ca[i/2][j/2] > -1){
		return ca[i/2][j/2];
	}
	else if((i == 0) && (j == 0)){
		float tempA[2];
		tempA[0] = a.coords[i];
		tempA[1] = a.coords[i+1];
		float tempB[2];
		tempB[0] = b.coords[j];
		tempB[1] = b.coords[j+1];
		ca[i/2][j/2] = euclideanDistance(2,tempA,tempB);
	}
	else if((i > 0) && (j == 0)){
		float tempA[2];
		tempA[0] = a.coords[i];
		tempA[1] = a.coords[i+1];
		float tempB[2];
		tempB[0] = b.coords[j];
		tempB[1] = b.coords[j+1];
		ca[i/2][j/2] = maxFloat(recursiveC(ca, a, b, i-2, 0), euclideanDistance(2,tempA,tempB));
	}
	else if((i == 0) && (j > 0)){
		float tempA[2];
		tempA[0] = a.coords[i];
		tempA[1] = a.coords[i+1];
		float tempB[2];
		tempB[0] = b.coords[j];
		tempB[1] = b.coords[j+1];
		ca[i/2][i/2] = maxFloat(recursiveC(ca, a, b, 0, j-2), euclideanDistance(2,tempA,tempB));
	}
    else if((i > 0) && (j > 0)){
		float tempA[2];
		tempA[0] = a.coords[i];
		tempA[1] = a.coords[i+1];
		float tempB[2];
		tempB[0] = b.coords[j];
		tempB[1] = b.coords[j+1];
		//printf("%f, %f %f %f\n", tempA[0], tempA[1],tempB[0],tempB[1]);
        ca[i/2][j/2] = maxFloat(minFloat(recursiveC(ca, a, b, i-2, j),recursiveC(ca, a, b, i-2, j-2),recursiveC(ca, a, b, i, j-2)),euclideanDistance(2,tempA,tempB));
	}
}

float discreteFrectet(vector a, vector b, int coordsCounter){
	float **ca;
	int i;
	ca = malloc(coordsCounter*sizeof(float*));
	for(i = 0 ; i < coordsCounter ; i++){
		ca[i] = malloc(coordsCounter*sizeof(float));
		for(int j = 0 ; j < coordsCounter ; j++)
			ca[i][j] = -1.0;
	}
	float result = recursiveC(ca,a,b,(coordsCounter*2)-2,(coordsCounter*2)-2);
	for(i = 0 ; i < coordsCounter ; i++){
		free(ca[i]);
	}
	free(ca);
    return result;	

}

void frechetSearch(vector* vectors, vector query, struct h_func **h, list_node ***HashTables, int vectorsCounter, int coordsCounter, int M, int k, int L, int w, int TableSize, double* distanceLSH, int *lsh_results, int N, int* radius_results, int radius, int* radiusCounter, double delta){

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
		vector querySnapped, queryCurved;
		querySnapped.coords = malloc(2*coordsCounter*sizeof(float));
		queryCurved.coords = malloc(2*coordsCounter*sizeof(float));
		for(j = 0 ; j < 2*coordsCounter ; j+=2){
			querySnapped.coords[j] = (j/2);
			queryCurved.coords[j] = (j/2);
			queryCurved.coords[j+1] = query.coords[j/2];
		}
		querySnapped.id =  malloc(sizeof(char)*(strlen(query.id)+1));
		strcpy(querySnapped.id,query.id);
		for(j = 0 ; j < 2*coordsCounter ; j+=2){
			float* x = query.coords;
			querySnapped.coords[j] = floor((querySnapped.coords[j]-h[z][0].tFrechet[j/2])/delta +1/2) * delta + h[z][0].tFrechet[j/2]; //x
			querySnapped.coords[j+1] = floor((x[j/2]-h[z][0].tFrechet[j/2])/delta +1/2) * delta + h[z][0].tFrechet[j/2]; //y
		}
		for(j = 2 ; j < coordsCounter*2 ; j+=2){
			if((querySnapped.coords[j] == querySnapped.coords[j-2]) && (querySnapped.coords[j+1] == querySnapped.coords[j-1])){
				for(int z = j ; z < coordsCounter*2 ; z+=2){
					querySnapped.coords[z] = querySnapped.coords[z+2];
					querySnapped.coords[z+1] = querySnapped.coords[z+3];
				}
				querySnapped.coords[(coordsCounter*2)-2] = 2147483647; //INT_MAX
				querySnapped.coords[(coordsCounter*2)-1] = 2147483647; //INT_MAX
			}
		}


		for(t=0; t<k; t++){								
			h[z][t].h_floor = 0;
			h[z][t].v = malloc(sizeof(float)* coordsCounter*2);
			normalDistr(h[z][t].v, coordsCounter*2);
			//printf("TEST %d ",vectors[i].coords);
			float pv = inner(querySnapped.coords,h[z][t].v,coordsCounter*2);
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
				//distance = euclideanDistance(coordsCounter, vectors[temp->vec_pos].coords, query.coords);
				
				vector vectorCurved;
				vectorCurved.coords = malloc((2*coordsCounter)*sizeof(float));
				for(j = 0 ; j < 2*coordsCounter ; j+=2){
					vectorCurved.coords[j] = (j/2);
					vectorCurved.coords[j+1] = vectors[temp->vec_pos].coords[j/2];
				}

				distance = discreteFrectet(queryCurved, vectorCurved, coordsCounter);
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
				free(vectorCurved.coords);
			}
			
			temp = temp->next;
		}
		free(querySnapped.coords);
		free(querySnapped.id);
		free(queryCurved.coords);
	}
	
	sortArrays(distanceLSH, lsh_results, N);

}

void frechet(vector *vectors, vector *queries, int vectorsCounter, int queriesCounter, int coordsCounter, int k, int L, int w, int N, int radius, char* queryFile, int gotQueries, char* outputFile, int gotOutput, double delta){
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
		h[i][0].tFrechet = malloc(sizeof(float)* coordsCounter);
		uniformDistr(h[i][0].tFrechet, coordsCounter, 0.0, delta); // uniform
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

	frechetTrain(vectors, h, HashTables, vectorsCounter, coordsCounter, M, k, L, w, TableSize, delta);
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
			frechetSearch(vectors, queries[i], h, HashTables, vectorsCounter, coordsCounter, M, k, L, w, TableSize, distanceLSH[i], lsh_results[i], N, temp_radius_results, radius, &radiusCounter, delta);
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
		writeResults(outputFile, queriesCounter,queries, vectors, lsh_results, distanceLSH, distanceTrue, tLSH, tTrue,trueResults, N, 2);

		for(i=0; i< queriesCounter; i++){
			free(queries[i].id);
			free(queries[i].coords);
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
		free(h[i][0].tFrechet);
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

void frechetContinuousSearch(vector* vectors, vector query, struct h_func *h, list_node **HashTables, int vectorsCounter, int coordsCounter, int M, int k, int L, int w, int TableSize, double* distanceLSH, int *lsh_results, int N, double delta){

	//printf("\nSearching..\n");
	int i, j, z, t, hash_pos, *a;
	int g;
	float f;
	srand(time(NULL));

	int flag = 0; 
	float distance;
	float minDistance;
	int minPosition;
	int counter = 0;

	g=0;
	vector querySnapped;
	querySnapped.coords = malloc(coordsCounter*sizeof(float));
	querySnapped.id =  malloc(sizeof(char)*(strlen(query.id)+1));
	strcpy(querySnapped.id,query.id);

	int filteringCounter = 0;
	double epsilon = 0.5;
	for(j = 0 ; j < coordsCounter ; j+=3){ //Filtering
		if((fabs(query.coords[j] - query.coords[j+1]) <= epsilon) && (fabs(query.coords[j+1] - query.coords[j+2]) <= epsilon)){
			query.coords[j+1] = -1.0;
			filteringCounter++;
		}
	}
	for(j = 0 ; j < coordsCounter ; j++){
		if(query.coords[j] == -1.0){
			for(int x = j ; x < coordsCounter ; x++){
				query.coords[j] = query.coords[j+1];
			}
			query.coords[coordsCounter-1] = -1.0;
		}
	}

	for(j = 0 ; j < coordsCounter - filteringCounter ; j++){ // snapping
		float* x = query.coords;
		querySnapped.coords[j] = floor((x[j]+h[0].tFrechet[j])/delta) * delta; //x
	}
	 //Min / Max
	for(j = 0 ; j < coordsCounter ; j+=3){
		if(((querySnapped.coords[j+1] > querySnapped.coords[j]) && (querySnapped.coords[j+1] < querySnapped.coords[j+2])) || (querySnapped.coords[j+1] == querySnapped.coords[j]) || (querySnapped.coords[j+1] == querySnapped.coords[j+2])){
			querySnapped.coords[j+1] = -1.0;
			filteringCounter++;
		}
	}
	for(j = 0 ; j < coordsCounter ; j++){
		if(querySnapped.coords[j] == -1.0){
			for(int x = j ; x < coordsCounter ; x++){
				querySnapped.coords[j] = querySnapped.coords[j+1];
			}
			querySnapped.coords[coordsCounter-1] = -1.0;
		}
	}

	for(j = coordsCounter-filteringCounter ; j < coordsCounter ; j++){
		querySnapped.coords[j] = 2147483647; //INT_MAX
	}


	for(t=0; t<k; t++){								
		h[t].h_floor = 0;
		h[t].v = malloc(sizeof(float)* coordsCounter);
		normalDistr(h[t].v, coordsCounter*2);
		//printf("TEST %d ",vectors[i].coords);
		float pv = inner(querySnapped.coords,h[t].v,coordsCounter*2);
		h[t].h_floor = floor((pv+h[t].t)/w);
		int r;
		r = rand();
		//g = g+h[z][t].r*h[z][t].h_floor;
		g = eucmod(eucmod(g,M)+eucmod(h[t].r*h[t].h_floor,M),M);
		free(h[t].v);
	}
	int ID = g;
	g = eucmod(ID,TableSize);
	//hash(HashTables[z], g, g, i);
	list_node *temp;
	temp = HashTables[g];
	double max;
	int maxPosition;
	

	while(temp != NULL) {//traverse list bucket
		if(temp->g == ID){
			//calculate distance between query and bucket vectors
			//store N nearest neighbors in sorted array
			flag = 0;
			//distance = euclideanDistance(coordsCounter, vectors[temp->vec_pos].coords, query.coords);
			
			/*vector vectorCurved;
			vectorCurved.coords = malloc((coordsCounter)*sizeof(float));
			for(j = 0 ; j < 2*coordsCounter ; j+=2){
				vectorCurved.coords[j] = (j/2);
				vectorCurved.coords[j+1] = vectors[temp->vec_pos].coords[j/2];
			}*/

			//distance = distance(queryCurved, vectorCurved, coordsCounter);
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
	
	
	sortArrays(distanceLSH, lsh_results, N);

}

void frechetContinuous(vector *vectors, vector *queries, int vectorsCounter, int queriesCounter, int coordsCounter, int k, int L, int w, int N, int radius, char* queryFile, int gotQueries, char* outputFile, int gotOutput, double delta){
	//first we save data set vertexes to hash tables and then we search with queries
	
	
	int i, j, z, M, hash_pos, TableSize;
	struct h_func *h; 
	list_node **HashTables;
	clock_t start, stop;

	h = malloc(k*sizeof(struct h_func));

	srand(time(0));
	h[0].tFrechet = malloc(sizeof(float)* coordsCounter);
	uniformDistr(h[0].tFrechet, coordsCounter, 0.0, delta); // uniform
	for(j=0; j<k; j++){
		h[j].t = w*(rand() / (RAND_MAX +1.0));
	}

	for(j=0; j<k; j++){
		h[j].r = rand();
		for(i = 1 ; i < L ; i++){
			h[j].r = h[j].r;
		}
	}
	
	TableSize = vectorsCounter/4;	
	HashTables = malloc(TableSize*sizeof(struct list_node *));
	for(j=0; j<TableSize; j++){
		HashTables[j] = NULL;
	}

	M = pow(2, 32/k);

	//frechetContinuousTrain(vectors, h, HashTables, vectorsCounter, coordsCounter, M, k, L, w, TableSize, delta);
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
			//frechetSearch(vectors, queries[i], h, HashTables, vectorsCounter, coordsCounter, M, k, L, w, TableSize, distanceLSH[i], lsh_results[i], N, temp_radius_results, radius, &radiusCounter, delta);
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
		writeResults(outputFile, queriesCounter,queries, vectors, lsh_results, distanceLSH, distanceTrue, tLSH, tTrue,trueResults, N, 2);

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

	free(h);

	for(j = 0; j < TableSize; j++){
		while(HashTables[j] != NULL){
			temp = HashTables[j];
			HashTables[j] = HashTables[j]->next;
			free(temp);
		}
		
	}
	free(HashTables);



	
}