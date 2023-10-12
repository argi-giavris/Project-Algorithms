#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "interfaceCluster.h"

vectorCluster* snap(vectorCluster* vectors, double delta,int vectorsCounter, int coordsCounter, float* t){
	int i, j;
	vectorCluster* vectorsSnapped = malloc(vectorsCounter*sizeof(vectorCluster));
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
			vectorsSnapped[i].coords[j] = floor((vectorsSnapped[i].coords[j]-t[j/2])/delta +1/2) * delta + t[j/2]; //x
			vectorsSnapped[i].coords[j+1] = floor((x[j/2]-t[j/2])/delta +1/2) * delta + t[j/2]; //y
		}
	}
	return vectorsSnapped;
}

float* MeanCurve(float* a, float* b, int coordsCounter){
    float* result = malloc(coordsCounter*sizeof(float));
    for(int i = 0 ; i < coordsCounter ; i++){
        if(b != NULL)
            result[i] = (a[i] + b[i]) / 2;
        else
            result[i] = a[i];
    }
    return result;
}

treeNode *getNewNode(float* data){
    treeNode *newNode = malloc(sizeof(treeNode));
    newNode->data = data;
    newNode->left  = NULL;
    newNode->right = NULL;

    return newNode;
}

treeNode* insertLevelOrder(float** arr, treeNode* root,int i, int n)
{
    // Base case for recursion
    if (i < n)
    {
        treeNode* temp = getNewNode(arr[i]);
        root = temp;
 
        // insert left child
        root->left = insertLevelOrder(arr,root->left, 2 * i + 1, n);
 
        // insert right child
        root->right = insertLevelOrder(arr,root->right, 2 * i + 2, n);
    }
    return root;
}

int isLeaf(treeNode* node){
    if((node->left == NULL) && (node->right == NULL)){
        return 1;
    }
    else
        return 0;
}

float* PostOrderTraversal(treeNode* node, int coordsCounter){
    float* leftCurve, *rightCurve;
    if(isLeaf(node) == 1)
        return node->data;
    else
       leftCurve = PostOrderTraversal(node->left,coordsCounter);
    if(node->right != NULL)
        rightCurve = PostOrderTraversal(node->right,coordsCounter);
    else
        rightCurve = NULL;
    return MeanCurve(leftCurve, rightCurve,coordsCounter);
}

void freeTree(treeNode * node){
   //post-order like FatalError hinted at
    if (node != NULL) {
        //if(isLeaf(node) == 0)
        //    free(node->data); //if data was heap allocated, need to free it
        freeTree(node->right);
        freeTree(node->left);
        free(node);
    }
}

int updateMeanCurves(vectorCluster* vectors, vectorCluster* centroids, int vectorsCounter, int numOfClusters, int coordsCounter){ //mean average update
    int i, j, vectorsInCluster, z, changed = 0;
    //printf("Updating Mean Curves\n");
	vectorCluster* previousCentroids = malloc(numOfClusters*sizeof(vectorCluster));
	for(i=0; i<numOfClusters; i++){ //create vectors to save previous centroids
		previousCentroids[i].coords = malloc(coordsCounter*sizeof(int));
        for(j = 0 ; j < coordsCounter ; j++){ // copy all coords and initialize new coords with 0
            previousCentroids[i].coords[j] = centroids[i].coords[j];
            centroids[i].coords[j] = 0;
        }
	}

    for(i = 0 ; i < numOfClusters ; i++){
        vectorsInCluster = 0;
        for(j = 0 ; j < vectorsCounter ; j++){
            if(vectors[j].nearest == i){
                vectorsInCluster++;
                //for(z = 0 ; z < coordsCounter ; z++){
                //    centroids[i].coords[z] += vectors[j].coords[z]; // calculate the sum of all vectors assigned to this centroid
                //}
            }
        }
        if(vectorsInCluster == 1){
            for(j = 0 ; j < coordsCounter ; j++){
                centroids[i].coords[j] = previousCentroids[i].coords[j];
            }
            free(previousCentroids[i].coords);
            continue;
        }
        int height = ceil(log(vectorsInCluster));
        int maxNodesOfTree = pow(2,height+1)-1;
        int emptyNodes = maxNodesOfTree/2;
        float** arrayTree = malloc((vectorsInCluster+emptyNodes) * sizeof(float*));
        int tempCounter = 0;
        for(tempCounter = 0 ; tempCounter < emptyNodes ; tempCounter++){
            arrayTree[tempCounter] = NULL;
        }
        for(j = 0 ; j < vectorsCounter ; j++){
            if(vectors[j].nearest == i){
                arrayTree[tempCounter] = vectors[j].coords;
                tempCounter++;
            }
        }

        treeNode* root = insertLevelOrder(arrayTree,root,0,vectorsInCluster+emptyNodes);
        centroids[i].coords = PostOrderTraversal(root,coordsCounter);
        freeTree(root);
        free(arrayTree);

        for(j = 0 ; j < coordsCounter ; j++){
            if(centroids[i].coords[j] != previousCentroids[i].coords[j]) changed = 1;
        }
        free(previousCentroids[i].coords);
    }
    free(previousCentroids);
    return changed;

}

void frechetTrain(vectorCluster *vectors, struct h_func **h, list_node ***HashTables, int vectorsCounter, int coordsCounter, int M, int k, int L, int w, int TableSize, double delta){
	int i, j, z, t, hash_pos, *a;
	int g;
	float f;
	srand(time(NULL));

	vectorCluster** vectorsSnapped = malloc(L*sizeof(vectorCluster*));
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

			hash(HashTables[z], g, ID, i);
		}
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

float recursiveC(float** ca,vectorCluster a, vectorCluster b, int i, int j){
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
        ca[i/2][j/2] = maxFloat(minFloat(recursiveC(ca, a, b, i-2, j),recursiveC(ca, a, b, i-2, j-2),recursiveC(ca, a, b, i, j-2)),euclideanDistance(2,tempA,tempB));
	}
}

float discreteFrectet(vectorCluster a, vectorCluster b, int coordsCounter){
	float **ca;
	ca = malloc(coordsCounter*sizeof(float*));
    int i;
	for(i = 0 ; i < coordsCounter ; i++){
		ca[i] = malloc(coordsCounter*sizeof(float));
		for( int j = 0 ; j < coordsCounter ; j++)
			ca[i][j] = -1.0;
	}
	float result = recursiveC(ca,a,b,(coordsCounter*2)-2,(coordsCounter*2)-2);
	for(i = 0 ; i < coordsCounter ; i++){
		free(ca[i]);
	}
	free(ca);
    return result;
}

void frechetSearch(vectorCluster* vectors, vectorCluster* centroids, struct h_func **h, list_node ***HashTables, int vectorsCounter, int coordsCounter, int M, int k, int L, int w, int TableSize, int numOfClusters, double delta){

	int i, j, z, t, hash_pos;
	int g;
	float f;

	float distance;
	float minDistance;
	int minPosition;
    float radii;


    for(i = 0 ; i < numOfClusters ; i++){  // calculate radii min distance among centroids
        for(j = 0 ; j < numOfClusters ; j++){
            if( j != i){
                distance = discreteFrectet(centroids[j], centroids[i], coordsCounter);
                if(i == 0 && j == 1){
                    minDistance = distance;
                }
                else if(distance < minDistance){
                    minDistance = distance;
                }
            }
        }
    }

    radii = minDistance/2;
    int stop = 1;  // flag to know when to stop iterations
    int assigned = 0;//to how many centroids at least a vector assigned
    int flag = 0;
    int firstIteration = 1;

    while((stop == 1) || (radii <= minDistance*2)){
    //int flag2 = 0;
    //while(flag2 == 0){
        //flag2 = 1;
        stop = 1;
        assigned = 0;  
	    for(i = 0 ; i < numOfClusters ; i++){
            flag = 0;
    	    for(z=0; z<L; z++){
	    	    g=0;
                vectorCluster centroidSnapped;
                centroidSnapped.coords = malloc(2*coordsCounter*sizeof(float));
                for(j = 0 ; j < 2*coordsCounter ; j+=2){
                    centroidSnapped.coords[j] = (j/2);
                }
                centroidSnapped.id =  malloc(sizeof(char)*(strlen(centroids[i].id)+1));
                strcpy(centroidSnapped.id,centroids[i].id);
                for(j = 0 ; j < 2*coordsCounter ; j+=2){
                    float* x = centroids[i].coords;
                    centroidSnapped.coords[j] = floor((centroidSnapped.coords[j]-h[z][0].tFrechet[j/2])/delta +1/2) * delta + h[z][0].tFrechet[j/2]; //x
                    centroidSnapped.coords[j+1] = floor((x[j/2]-h[z][0].tFrechet[j/2])/delta +1/2) * delta + h[z][0].tFrechet[j/2]; //y
                }
                for(j = 2 ; j < coordsCounter*2 ; j+=2){
                    if((centroidSnapped.coords[j] == centroidSnapped.coords[j-2]) && (centroidSnapped.coords[j+1] == centroidSnapped.coords[j-1])){
                        for(int z = j ; z < coordsCounter*2 ; z+=2){
                            centroidSnapped.coords[z] = centroidSnapped.coords[z+2];
                            centroidSnapped.coords[z+1] = centroidSnapped.coords[z+3];
                        }
                        centroidSnapped.coords[(coordsCounter*2)-2] = 2147483647; //INT_MAX
                        centroidSnapped.coords[(coordsCounter*2)-1] = 2147483647; //INT_MAX
                    }
        		}

		        for(t=0; t<k; t++){							
		    	    h[z][t].h_floor = 0;
    		    	h[z][t].v = malloc(sizeof(float)* coordsCounter*2);
	    	    	normalDistr(h[z][t].v, coordsCounter*2);

    	    		float pv = inner(centroids[i].coords,h[z][t].v,coordsCounter*2);
     	    		h[z][t].h_floor =floor((pv+h[z][t].t)/w);

        			g = eucmod(eucmod(g,M)+eucmod(h[z][t].r*h[z][t].h_floor,M),M);
                    free(h[z][t].v);
    	    	}
    		    int ID = g;
                g = eucmod(ID,TableSize);
                
                list_node *temp;
                temp = HashTables[z][g];
                double max;
                int maxPosition;

                while(temp != NULL) {
                    if(temp->g == ID){
                        distance = discreteFrectet(centroids[i], vectors[temp->vec_pos], coordsCounter);
                        if(distance <= radii){
                            if(vectors[temp->vec_pos].nearest == -1){	// if this vector doesnt have a centroid
                                vectors[temp->vec_pos].nearest = i;
                              
                                flag = 1;
                            }
                            else{
                                if(distance < discreteFrectet(centroids[vectors[temp->vec_pos].nearest], vectors[temp->vec_pos], coordsCounter)){
                                    vectors[temp->vec_pos].nearest = i;	// this centroid is better
                                    flag = 1;
                                }
                            }
                        }
                    }
                
                    temp = temp->next;
                }
		        free(centroidSnapped.coords);
		        free(centroidSnapped.id);
            }
            if(flag == 1) assigned++;
        }
        if(assigned < 0.5*numOfClusters){
            stop = 0;
        }
        radii = radii*2;
    }

	for(i=0; i<vectorsCounter; i++){
		if(vectors[i].nearest == -1){	// if there wasnt a centroid in the same bucket
			for(j=0; j<numOfClusters; j++){
				distance = discreteFrectet(vectors[i], centroids[j], coordsCounter);
                if(j == 0){
                    vectors[i].nearest = j;
                    minDistance = distance;
                }
				if(distance  < minDistance ){
					vectors[i].nearest = j;
					minDistance = distance;
				}
			}
		}
	}

}

void main (int argc, char *argv[]){
	int i, kCube = 3, M = 10, L = 3, probes = 2,  vectorsCounter, coordsCounter, method = 0, numOfClusters, kLSH = 4, complete = 0, silhouetteFlag = 0, updateFlag = 0;
	char ch, inputFile[256], configurationFile[256], outputFile[256];
	vectorCluster *vectors;
	int argCount = 0;
    double delta = 0.1;
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
        else if (strcmp(argv[argCount],"-assignment") == 0){
		    if(strcmp(argv[argCount+1],"Classic") == 0) //Classic -> method = 0
                method = 0;
            else if(strcmp(argv[argCount+1],"LSH") == 0) //LSH -> method = 1
                method = 2;
            else if(strcmp(argv[argCount+1],"Hypercube") == 0) //HyperCube -> method = 2
                method = 3;
            else
                method = 1;
        }
        else if (strcmp(argv[argCount],"-update") == 0){
		    if(strcmp(argv[argCount+2],"Frechet") == 0) //Mean Frechet;
                updateFlag = 1;
            else
                updateFlag = 0; //Mean Vector
        }
        else if (strcmp(argv[argCount],"-o") == 0){
		    strcpy(outputFile, argv[argCount+1]);
        }
        else if(strcmp(argv[argCount],"-silhouette") == 0){
            silhouetteFlag = 1;
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
		vectors[i].coords = malloc(coordsCounter*sizeof(float));
	}
	saveFile(inputFile, vectors);


    clock_t start, stop;
    double clusterTime;
    srand(time(NULL));
    start = clock();

    vectorCluster* centroids;
    vectorCluster* vectorsCurved;
    
    //initialize centroids
    if(method != 1){
    	centroids = malloc(numOfClusters*sizeof(vectorCluster));
    	for(i=0; i<numOfClusters; i++){
    		centroids[i].coords = malloc(coordsCounter*sizeof(float));
    	}
        kmeansInit(vectors,vectorsCounter,coordsCounter,numOfClusters,centroids);
    }

    int j = 0;
    if(method == 0){ // Llloyd
        lloydAssignment(vectors, centroids, vectorsCounter, coordsCounter, numOfClusters);

        /* UPDATE */
        int updatesCounter = 0;
        if(updateFlag == 1){
            while((update(vectors,centroids,vectorsCounter,numOfClusters,coordsCounter) == 1) && (updatesCounter < 100)){
                updatesCounter++;
                lloydAssignment(vectors, centroids, vectorsCounter, coordsCounter, numOfClusters);
            }
        }
        else{
            while((updateMeanCurves(vectors,centroids,vectorsCounter,numOfClusters,coordsCounter) == 1) && (updatesCounter < 100)){
                updatesCounter++;
                lloydAssignment(vectors, centroids, vectorsCounter, coordsCounter, numOfClusters);
            }            
        }

    }
    else if(method == 2){ //LSH
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
    else if(method == 3){ //HYPERCUBE
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
    else if(method == 1){ //Frechet
        int w = 1000;

    	centroids = malloc(numOfClusters*sizeof(vectorCluster));
    	for(i=0; i<numOfClusters; i++){
    		centroids[i].coords = malloc(coordsCounter*2*sizeof(float));
    	}
        
        vectorsCurved = malloc(vectorsCounter*sizeof(vectorCluster));
        for(i=0; i<vectorsCounter; i++){
            vectorsCurved[i].coords = malloc((2*coordsCounter)*sizeof(float));
            vectorsCurved[i].id = malloc(sizeof(char)*(strlen(vectors[i].id)+1));
            strcpy(vectorsCurved[i].id,vectors[i].id);
        }
        for(i = 0 ; i < vectorsCounter ; i++){
            for(j = 0 ; j < 2*coordsCounter ; j+=2){
                vectorsCurved[i].coords[j] = (j/2);
                vectorsCurved[i].coords[j+1] = vectors[i].coords[j/2];
            }
        }
        printf("kmeans init\n");
        kmeansInit(vectorsCurved,vectorsCounter,coordsCounter*2,numOfClusters,centroids);

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
            h[i][0].tFrechet = malloc(sizeof(float)* coordsCounter);
            uniformDistr(h[i][0].tFrechet, coordsCounter, 0.0, delta); // uniform
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
        frechetTrain(vectors, h, HashTables, vectorsCounter, coordsCounter, M, kLSH, L, w, TableSize, delta);
        for( i = 0 ; i < vectorsCounter ; i++){
            vectorsCurved[i].nearest = -1; // vector is not asigned to a centroid yet
        }
        frechetSearch(vectorsCurved, centroids, h, HashTables, vectorsCounter, coordsCounter, M, kLSH, L, w, TableSize, numOfClusters, delta);

        /* UPDATE */
        int updatesCounter = 0;
        while((updateMeanCurves(vectorsCurved,centroids,vectorsCounter,numOfClusters,coordsCounter*2) == 1) && (updatesCounter < 100)){
            printf("Updated\n");
            updatesCounter++;
            frechetSearch(vectorsCurved, centroids, h, HashTables, vectorsCounter, coordsCounter, M, kLSH, L, w, TableSize, numOfClusters, delta);
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
	stop = clock();
	clusterTime = (double)( stop - start) / CLOCKS_PER_SEC;;
    float stotal = 0.0;
    float* silhouettes = NULL;
    if(silhouetteFlag == 1){
        if(method == 1){ //Frechet
            for(i = 0 ; i < vectorsCounter ; i++){
                vectors[i].nearest = vectorsCurved[i].nearest;
            }
        }
        silhouettes = silhouette(numOfClusters,coordsCounter,vectors,centroids,vectorsCounter,&stotal);

    }
    if(method != 1)
        writeResults(outputFile,method,numOfClusters,vectors,vectorsCounter,coordsCounter,centroids,clusterTime,complete,silhouettes,&stotal,updateFlag);
    else{
        writeResults(outputFile,method,numOfClusters,vectorsCurved,vectorsCounter,coordsCounter,centroids,clusterTime,complete,silhouettes,&stotal,updateFlag);
	    for(i=0; i<vectorsCounter; i++){
	    	free(vectorsCurved[i].coords);
            free(vectorsCurved[i].id);
	    }
    }

    // frees

    free(silhouettes);
	for(i=0; i<vectorsCounter; i++){
		free(vectors[i].coords);
        free(vectors[i].id);
	}
    free(vectors);
    for(i = 0 ; i < numOfClusters ; i++){
        free(centroids[i].coords);
    }
    free(centroids);
}