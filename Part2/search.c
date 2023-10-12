#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "interface.h"

int main (int argc, char *argv[]){
	int i, k = -1 , L = 5, vectorsCounter, queriesCounter, coordsCounter, N = 1, gotInput = 0, gotQueries = 0, gotOutput = 0, M = 10, probes = 2, algorithmFlag = 0, metricFlag = 0;
	char ch, inputFile[256], queryFile[256], outputFile[256];
	double radius = 10, delta = 0.5;
	vector *vectors, *queries;
	int argCount = 0;


	//read command line parameters
    for(argCount = 0 ; argCount < argc ; argCount++){
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
        else if (strcmp(argv[argCount],"-L") == 0){
		    L = atoi(argv[argCount+1]);
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
        else if (strcmp(argv[argCount],"-algorithm") == 0){
			if(strcmp(argv[argCount+1],"LSH") == 0)
				algorithmFlag = 0;
			else if(strcmp(argv[argCount+1],"Hypercube") == 0)
				algorithmFlag = 1;
			else
				algorithmFlag = 2; //Frechet
        }
        else if (strcmp(argv[argCount],"-metric") == 0){
		    if(strcmp(argv[argCount+1],"discrete") == 0)
				metricFlag = 0;
			else
				metricFlag = 1; //continuous
        }
        else if (strcmp(argv[argCount],"-delta") == 0){
		    delta = atof(argv[argCount+1]);
        }
    }

	if(algorithmFlag == 0){ //LSH
		if(k == -1) k = 4;
		//if not given from command line
		if(gotInput == 0){
			printf("Insert input file path:\n");
			scanf("%s",inputFile);
		}

		countFile(inputFile, &vectorsCounter, &coordsCounter);
		vectors = malloc(vectorsCounter*sizeof(vector));
		i=0;
		for(i=0; i<vectorsCounter; i++){
			vectors[i].coords = malloc(coordsCounter*sizeof(float));
		}
		printf("saving input\n");
		saveFile(inputFile, vectors);

		int w = 100;

		lsh(vectors, queries, vectorsCounter, queriesCounter, coordsCounter, k, L, w, N, radius, queryFile, gotQueries,outputFile,gotOutput);
		//int flag = 1;
		//writeResults(outputFile, queriesCounter,queries, vectors, lsh_results, distanceLSH, distanceTrue, tLSH, tTrue,radius_results, N, flag);

		//frees

		for(i=0; i< vectorsCounter; i++){
			free(vectors[i].coords);
			free(vectors[i].id);
		}
		free(vectors);
	}
	else if(algorithmFlag == 1){ //HyperCube
		if(k == -1) k = 14;
		countFile(inputFile, &vectorsCounter, &coordsCounter);
		vectors = malloc(vectorsCounter*sizeof(vector));
		i=0;
		for(i=0; i<vectorsCounter; i++){
			vectors[i].coords = malloc(coordsCounter*sizeof(int));
		}
		saveFile(inputFile, vectors);

		int w = 100; //w is averagedistance of vectors * 4;

		cube(vectors,queries, vectorsCounter, queriesCounter, coordsCounter, k, M, w, N, radius, probes, queryFile,gotQueries,outputFile,gotOutput);

		//frees

		for(i = 0; i < vectorsCounter; i++){
			free(vectors[i].coords);
			free(vectors[i].id);
		}
		free(vectors);
	}
	else{ //Frechet
		if(metricFlag == 0){ //Discrete
			if(k == -1) k = 4;
			//if not given from command line
			if(gotInput == 0){
				printf("Insert input file path:\n");
				scanf("%s",inputFile);
			}

			countFile(inputFile, &vectorsCounter, &coordsCounter);
			vectors = malloc(vectorsCounter*sizeof(vector));
			i=0;
			for(i=0; i<vectorsCounter; i++){
				vectors[i].coords = malloc(coordsCounter*sizeof(float));
			}
			printf("saving input\n");
			saveFile(inputFile, vectors);


			//floor(x/delta +1/2) * delta;

			int w = 100;

			frechet(vectors, queries, vectorsCounter, queriesCounter, coordsCounter, k, L, w, N, radius, queryFile, gotQueries,outputFile,gotOutput,delta);
			//int flag = 1;
			//writeResults(outputFile, queriesCounter,queries, vectors, lsh_results, distanceLSH, distanceTrue, tLSH, tTrue,radius_results, N, flag);

			//frees

			for(i=0; i< vectorsCounter; i++){
				free(vectors[i].coords);
				free(vectors[i].id);
			}
			free(vectors);		
		}
		else{ //Continuous
			if(k == -1) k = 4;
			//if not given from command line
			if(gotInput == 0){
				printf("Insert input file path:\n");
				scanf("%s",inputFile);
			}

			countFile(inputFile, &vectorsCounter, &coordsCounter);
			vectors = malloc(vectorsCounter*sizeof(vector));
			i=0;
			for(i=0; i<vectorsCounter; i++){
				vectors[i].coords = malloc(coordsCounter*sizeof(float));
			}
			printf("saving input\n");
			saveFile(inputFile, vectors);


			//floor(x/delta +1/2) * delta;

			int w = 100;

			frechet(vectors, queries, vectorsCounter, queriesCounter, coordsCounter, k, L, w, N, radius, queryFile, gotQueries,outputFile,gotOutput,delta);
			//int flag = 1;
			//writeResults(outputFile, queriesCounter,queries, vectors, lsh_results, distanceLSH, distanceTrue, tLSH, tTrue,radius_results, N, flag);

			//frees

			for(i=0; i< vectorsCounter; i++){
				free(vectors[i].coords);
				free(vectors[i].id);
			}
			free(vectors);		
		}
	}

	return 0;
}