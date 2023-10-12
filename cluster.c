#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "interfaceCluster.h"

void main (int argc, char *argv[]){
	int i, kCube = 3, M = 10, L = 3, probes = 2,  vectorsCounter, coordsCounter, method = 0, numOfClusters, kLSH = 4, complete = 0;
	char ch, inputFile[256], configurationFile[256], outputFile[256];
	vectorCluster *vectors;
	int argCount = 0;

    //read command line arguments

    for(argCount = 0 ; argCount < argc ; argCount++){
        if(strcmp(argv[argCount],"-i") == 0){
            strcpy(inputFile, argv[argCount+1]);
        }
        else if (strcmp(argv[argCount],"-c") == 0){
            strcpy(configurationFile, argv[argCount+1]);
        }
        else if (strcmp(argv[argCount],"-complete") == 0){
		    complete = 1;
        }
        else if (strcmp(argv[argCount],"-m") == 0){
		    if(strcmp(argv[argCount+1],"Classic") == 0) //Classic -> method = 0
                method = 0;
            else if(strcmp(argv[argCount+1],"LSH") == 0) //LSH -> method = 1
                method = 1;
            else //HyperCube -> method = 2
                method = 2;
        }  
        else if (strcmp(argv[argCount],"-o") == 0){
		    strcpy(outputFile, argv[argCount+1]);
        }
    }

	FILE *fp;
    ssize_t line_size;
    char *line_buf = NULL;
    size_t line_buf_size = 0;
    int line_count = 0;
    char* token;

	fp = fopen(configurationFile,"r");
    line_size = getline(&line_buf, &line_buf_size, fp);
    while (line_size > 0){
        line_count++;

        token = strtok(line_buf, " ");
        if(strcmp(token,"number_of_clusters:") == 0){
            numOfClusters = atoi(strtok(NULL," "));

        }
        else if(strcmp(token,"number_of_vector_hash_tables:") == 0)
            L = atoi(strtok(NULL," "));
        else if(strcmp(token,"number_of_vector_hash_functions:") == 0)
            kLSH = atoi(strtok(NULL," "));
        else if(strcmp(token,"max_number_M_hypercube:") == 0)
            M = atoi(strtok(NULL," "));
        else if(strcmp(token,"number_of_hypercube_dimensions:") == 0)
            kCube = atoi(strtok(NULL," "));        
        else if(strcmp(token,"number_of_probes:") == 0)
            probes = atoi(strtok(NULL," "));


        line_size = getline(&line_buf, &line_buf_size, fp);
    }
    free(line_buf);

    fclose(fp);


	countFile(inputFile, &vectorsCounter, &coordsCounter);
	vectors = malloc(vectorsCounter*sizeof(vectorCluster));
	i=0;
	for(i=0; i<vectorsCounter; i++){
		vectors[i].coords = malloc(coordsCounter*sizeof(int));
	}
	saveFile(inputFile, vectors);

	vectorCluster* centroids = malloc(numOfClusters*sizeof(vectorCluster));
	for(i=0; i<numOfClusters; i++){
		centroids[i].coords = malloc(coordsCounter*sizeof(int));
	}

    clock_t start, stop;
    double clusterTime;
    srand(time(NULL));
    start = clock();
    
    //initialize centroids
    kmeansInit(vectors,vectorsCounter,coordsCounter,numOfClusters,centroids);

    int j = 0;
    if(method == 0){ // Llloyd
        lloydAssignment(vectors, centroids, vectorsCounter, coordsCounter, numOfClusters);

        /* UPDATE */
        int updatesCounter = 0;
        while((update(vectors,centroids,vectorsCounter,numOfClusters,coordsCounter) == 1) && (updatesCounter < 100)){
            updatesCounter++;
            lloydAssignment(vectors, centroids, vectorsCounter, coordsCounter, numOfClusters);
        }

    }
    else if(method == 1){ //LSH
        int w = 1000;
	    struct h_func** h = malloc(L*sizeof(struct h_func *));
	    for(i=0; i<L; i++){	
		    h[i] = malloc(kLSH*sizeof(struct h_func));
	    }
        int TableSize = vectorsCounter/8;
	    list_node*** HashTables = malloc(L*sizeof(struct list_node **));
	    for(i=0; i<L; i++){
		    HashTables[i] = malloc(TableSize*sizeof(struct list_node *));
	    }
	    for(i=0; i<L; i++){	
		    for(j=0; j<TableSize; j++){
			    HashTables[i][j] = NULL;
		    }
	    }
        for(i=0; i<L; i++){
            for(j=0; j<kLSH; j++){
                h[i][j].t = w*(rand() / (RAND_MAX +1.0));
            }
        }

	    for(j=0; j<kLSH; j++){
		    h[0][j].r = rand();
		    for(i = 1 ; i < L ; i++){
			    h[i][j].r = h[0][j].r;
		    }
	    }

        M = pow(2, 32/kLSH);
        lshTrain(vectors, h, HashTables, vectorsCounter, coordsCounter, M, kLSH, L, w, TableSize);
        for( i = 0 ; i < vectorsCounter ; i++){
            vectors[i].nearest = -1; // vector is not asigned to a centroid yet
        }
        lshSearch(vectors, centroids, h, HashTables, vectorsCounter, coordsCounter, M, kLSH, L, w, TableSize, numOfClusters);

        /* UPDATE */
        int updatesCounter = 0;
        while((update(vectors,centroids,vectorsCounter,numOfClusters,coordsCounter) == 1) && (updatesCounter < 100)){
            updatesCounter++;
            lshSearch(vectors, centroids, h, HashTables, vectorsCounter, coordsCounter, M, kLSH, L, w, TableSize, numOfClusters);
        }

        list_node* temp;
	    for(i=0; i<L; i++){	
		    free(h[i]);
            for(j = 0; j < TableSize; j++){
                while(HashTables[i][j] != NULL){
                    temp = HashTables[i][j];
                    HashTables[i][j] = HashTables[i][j]->next;
                    free(temp);
                }
                
            }
            free(HashTables[i]);
	    }
        free(h);
        free(HashTables);


    }
    else if(method == 2){ //HYPERCUBE
        int w = 1000;

        struct h_func* h = malloc(kCube*sizeof(struct h_func));

        srand(time(0));
        for(j=0; j<kCube; j++){
            h[j].t = w*(rand() / (RAND_MAX +1.0));
    	}

        int cubeSize = pow(2,kCube);
        list_node** cube = malloc(cubeSize*sizeof(list_node *));
        for(i=0; i<cubeSize; i++){
            cube[i] = NULL;
    	}

        list_node **f = malloc(kCube*sizeof(list_node*));

        for(i = 0 ; i < kCube ; i++){
            f[i] = NULL;
        }

        cubeTrain(cube, vectors,h, f, vectorsCounter, kCube,coordsCounter,w);
        for( i = 0 ; i < vectorsCounter ; i++){
            vectors[i].nearest = -1; // vector is not assigned to a centroid yet
        }
        cubeSearch(vectors, centroids, h, cube, f, vectorsCounter, coordsCounter, M, kCube, probes, w, numOfClusters);

        /* UPDATE */
        int updatesCounter = 0;
        while((update(vectors,centroids,vectorsCounter,numOfClusters,coordsCounter) == 1) && (updatesCounter < 100)){
            updatesCounter++;
            cubeSearch(vectors, centroids, h, cube, f, vectorsCounter, coordsCounter, M, kCube, probes, w, numOfClusters);
        }

        list_node *temp;

        free(h);

        for(i=0; i<cubeSize; i++){
            
            while(cube[i] != NULL){
                temp = cube[i];
                cube[i] = cube[i]->next;
                free(temp);
            }
        }
        free(cube);

        for(i=0; i<kCube; i++){
            
            while(f[i] != NULL){
                temp = f[i];
                f[i] = f[i]->next;
                free(temp);
            }
        }
        
    	free(f);

    }
	stop = clock();
	clusterTime = (double)( stop - start) / CLOCKS_PER_SEC;;
    float stotal = 0.0;
    float* silhouettes = silhouette(numOfClusters,coordsCounter,vectors,centroids,vectorsCounter,&stotal);

    writeResults(outputFile,method,numOfClusters,vectors,vectorsCounter,coordsCounter,centroids,clusterTime,complete,silhouettes,&stotal);

    // frees

    free(silhouettes);
	for(i=0; i<vectorsCounter; i++){
		free(vectors[i].coords);
	}
    free(vectors);
    for(i = 0 ; i < numOfClusters ; i++){
        free(centroids[i].coords);
    }
    free(centroids);
}