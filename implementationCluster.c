#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "interfaceCluster.h"

int vectorIsCentroid(vectorCluster* centroids, int id, int numOfClusters){
    int i = 0;
    for(i = 0 ; i < numOfClusters ; i ++){
        if(centroids[i].id == id) return 1;
    }
    return 0;
}

void kmeansInit(vectorCluster* vectors, int vectorsCounter, int coordsCounter, int numOfClusters, vectorCluster* centroids){
    double* distances = malloc((vectorsCounter)*sizeof(double));
    double* P = malloc((vectorsCounter)*sizeof(double));
    double tempDistance;
    int i = 0, centroidsCounter = 1, j = 0;
    int random = rand()%vectorsCounter;
    for(i = 0 ; i < coordsCounter ; i++){
        centroids[0].coords[i] = vectors[random].coords[i];
    }
    centroids[0].id = vectors[random].id;
    
    tempDistance = euclideanDistance(coordsCounter, vectors[0].coords, centroids[0].coords);
    P[0] = 0;
    while(centroidsCounter < numOfClusters){//repeat until initialize all centroids
		for(int i = 1; i < vectorsCounter; i++){
			float max;
			float min;                                     
			for(int j = 0; j < centroidsCounter; j++){				
				tempDistance = euclideanDistance(coordsCounter, vectors[i].coords, centroids[j].coords);
				if(j == 0){
                    min = tempDistance;
                    max = tempDistance;
                }
				if(tempDistance < min)
					min = (float)tempDistance;
				if(tempDistance > max)
					max = (float)tempDistance;
			}

			if (centroidsCounter > 1){
				min = min / max;										
			}
			min = pow(min,2);
			min += P[i-1];
			P[i] = min;
		}
		
        int lower = 0;
        int upper = P[vectorsCounter - centroidsCounter];
        float num = (rand() % (upper - lower + 1)) + lower;
		
        int y = 0;
		for(y = 0; y < vectorsCounter; y++){
			if((num) < P[y] && (vectorIsCentroid(centroids, vectors[y].id, centroidsCounter) == 0)){
				break;
			}
		}
        for(i = 0 ; i < coordsCounter ; i++){
            centroids[centroidsCounter].coords[i] = vectors[y].coords[i];
        }
        centroids[centroidsCounter].id = vectors[y].id;
		centroidsCounter++;
	}

    //frees
    free(distances);
    free(P);

}

void lloydAssignment(vectorCluster *vectors, vectorCluster *centroids, int vectorsCounter, int coordsCounter, int numOfClusters){
	int i, j;
	double min = 0.0, distance;
	
	for(i = 0 ; i < vectorsCounter ; i++){
		for(j = 0 ; j < numOfClusters ; j++){
			distance = euclideanDistance(coordsCounter,vectors[i].coords, centroids[j].coords);
            if(j == 0){
                vectors[i].nearest = j;
                min = distance;
            }
			else if(distance < min ){
				vectors[i].nearest = j;
				min = distance;
                
			}
		}
	}
}

void lshTrain(vectorCluster *vectors, struct h_func **h, list_node ***HashTables, int vectorsCounter, int coordsCounter, int M, int k, int L, int w, int TableSize){
	int i, j, z, t, hash_pos, *a;
	int g;
	float f;

	for(i=0; i<vectorsCounter; i++){
		for(z=0; z<L; z++){
			g=0;
			for(t=0; t<k; t++){								
				h[z][t].h_floor = 0;
				h[z][t].v = malloc(sizeof(float)* coordsCounter);
				normalDistr(h[z][t].v, coordsCounter);
				float pv = inner(vectors[i].coords,h[z][t].v,coordsCounter);
				h[z][t].h_floor =floor((pv+h[z][t].t)/w);
				g = g+h[z][t].r*h[z][t].h_floor;
                free(h[z][t].v);
			}
			int ID = eucmod(g,M);
			g = eucmod(ID,TableSize);
			
			hash(HashTables[z], g, ID, i);
		}
	}
}

void lshSearch(vectorCluster* vectors, vectorCluster* centroids, struct h_func **h, list_node ***HashTables, int vectorsCounter, int coordsCounter, int M, int k, int L, int w, int TableSize, int numOfClusters){

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
                distance = euclideanDistance(coordsCounter,centroids[i].coords, centroids[j].coords);
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
        stop = 1;
        assigned = 0;  
	    for(i = 0 ; i < numOfClusters ; i++){
            flag = 0;
    	    for(z=0; z<L; z++){
	    	    g=0;
		        for(t=0; t<k; t++){							
		    	    h[z][t].h_floor = 0;
    		    	h[z][t].v = malloc(sizeof(float)* coordsCounter);
	    	    	normalDistr(h[z][t].v, coordsCounter);

    	    		float pv = inner(centroids[i].coords,h[z][t].v,coordsCounter);
     	    		h[z][t].h_floor =floor((pv+h[z][t].t)/w);

        			g = g+h[z][t].r*h[z][t].h_floor;
                    free(h[z][t].v);
    	    	}
    		    int ID = eucmod(g,M);
                g = eucmod(ID,TableSize);
                
                list_node *temp;
                temp = HashTables[z][g];
                double max;
                int maxPosition;

                while(temp != NULL) {
                    if(temp->g == ID){
                        distance = euclideanDistance(coordsCounter,vectors[temp->vec_pos].coords,centroids[i].coords);
                        if(distance <= radii){
                            if(vectors[temp->vec_pos].nearest == -1){	// if this vector doesnt have a centroid
                                vectors[temp->vec_pos].nearest = i;
                              
                                flag = 1;
                            }
                            else{
                                if(distance < euclideanDistance(coordsCounter,vectors[temp->vec_pos].coords, centroids[vectors[temp->vec_pos].nearest].coords)){
                                    vectors[temp->vec_pos].nearest = i;	// this centroid is better
                                    flag = 1;
                                }
                            }
                        }
                    }
                
                    temp = temp->next;
                }
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
				distance = euclideanDistance(coordsCounter, vectors[i].coords, centroids[j].coords);
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

void cubeTrain(list_node** cube, vectorCluster * vectors, struct h_func * h, list_node** f, int vectorsCounter, int k, int coordsCounter,int w){
	int i = 0, z = 0, t = 0;
	srand(time(NULL));
	int cubePos = 0;
	list_node* temp;

	for(i=0; i<vectorsCounter; i++){
		cubePos = 0;
		for(t=0; t<k; t++){			
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
		if(cube[cubePos] == NULL){
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

void cubeSearch(vectorCluster* vectors, vectorCluster* centroids, struct h_func *h, list_node **cube, list_node** f, int vectorsCounter, int coordsCounter, int M, int k, int probes, int w, int numOfClusters){

	int j, h_floor,t, i;
	int *binaryString = malloc(k*sizeof(int));
	int* cubePos = malloc(probes*sizeof(int));
	list_node* temp;
    float distance, minDistance;

    for(i = 0 ; i < numOfClusters ; i++){  // calculate radii min distance among centroids
        for(j = 0 ; j < numOfClusters ; j++){
            if( j != i){
                distance = euclideanDistance(coordsCounter,centroids[i].coords, centroids[j].coords);
                if(i == 0 && j == 1){
                    minDistance = distance;
                }
                else if(distance < minDistance){
                    minDistance = distance;
                }
            }
        }
    }

    float radii = minDistance/2;
    int stop = 1;  // flag to know when to stop iterations
    int assigned = 0;
    int flag = 0;
    int firstIteration = 1;

    while((stop == 1) || (radii <= minDistance*2)){
        stop = 1;
        assigned = 0;  
	    for(i = 0 ; i < numOfClusters ; i++){
            flag = 0;
            for(t=0; t<k; t++){								
                h[t].h_floor = 0;
                h[t].v = malloc(sizeof(float)* coordsCounter);
                normalDistr(h[t].v, coordsCounter);
                float pv = inner(centroids[i].coords,h[t].v,coordsCounter);
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
                    
                    distance = euclideanDistance(coordsCounter,vectors[temp->vec_pos].coords,centroids[i].coords);
                    if(distance <= radii){
                        if(vectors[temp->vec_pos].nearest == -1){	// if this vector doesnt have a centroid
                            vectors[temp->vec_pos].nearest = i;
                              
                            flag = 1;
                        }
                        else{
                            if(distance < euclideanDistance(coordsCounter,vectors[temp->vec_pos].coords, centroids[vectors[temp->vec_pos].nearest].coords)){
                                vectors[temp->vec_pos].nearest = i;	// this centroid is better
                                flag = 1;
                            }
                        }
                    }
                    counter++;
                    temp = temp->next;
                }
                
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
				distance = euclideanDistance(coordsCounter, vectors[i].coords, centroids[j].coords);
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

    //frees

    free(cubePos);
    free(binaryString);

}

int update(vectorCluster* vectors, vectorCluster* centroids, int vectorsCounter, int numOfClusters, int coordsCounter){ //mean average update
    int i, j, vectorsInCluster, z, changed = 0;
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
                for(z = 0 ; z < coordsCounter ; z++){
                    centroids[i].coords[z] += vectors[j].coords[z]; // calculate the sum of all vectors asiigned to this centroid
                }
            }
        }
        for(j = 0 ; j < coordsCounter ; j++){
            centroids[i].coords[j] = centroids[i].coords[j] / vectorsInCluster; // find the average
            if(centroids[i].coords[j] != previousCentroids[i].coords[j]) changed = 1;
        }
        free(previousCentroids[i].coords);
    }
    free(previousCentroids);
    return changed;

}

float* silhouette(int numOfClusters, int coordsCounter, vectorCluster* vectors, vectorCluster* centroids, int vectorsCounter, float* stotal){
    float min = 0.0, distance;
    int j, i;
    int* secondNearest = malloc(vectorsCounter*sizeof(int));
    for(i = 0 ; i < vectorsCounter ; i++){
        secondNearest[i] = -1;
    }

    for(i = 0 ; i < vectorsCounter ; i++){ //find Second nearest Cluster for each vector
        min = 0.0;
        for(j = 0; j < numOfClusters ; j++){ 
            distance = euclideanDistance(coordsCounter,vectors[i].coords, centroids[j].coords);
            if((min == 0.0) && (vectors[i].nearest != j)){
                secondNearest[i] = j;
                min = distance;              
            }
            else if((distance < min) && (vectors[i].nearest != j)){
                secondNearest[i] = j;
                min = distance;
            }
        }
    }
    /*Calculate silhouette type for each cluster*/
    float* silhouettes = malloc(numOfClusters*sizeof(float));
    float* distancesA = malloc(vectorsCounter*sizeof(float));
    float* distancesB = malloc(vectorsCounter*sizeof(float));
    int* vectorsInCluster = malloc(numOfClusters*sizeof(int));
    for(i = 0 ; i < numOfClusters ; i++){
        silhouettes[i] = 0.0;
        vectorsInCluster[i] = 0;
    }
    int distancesCounterA, distancesCounterB;
    int z;
    for(j = 0 ; j < vectorsCounter ; j++){ // Calculating S(p)
        distancesA[j] = 0.0;
        distancesB[j] = 0.0;
        distancesCounterA = 0;
        distancesCounterB = 0;
        vectorsInCluster[vectors[j].nearest]++;
        for(z = 0 ; z < vectorsCounter ; z++){
            if((j != z) && (vectors[j].nearest == vectors[z].nearest)){ //distance with vectors in same cluster
                distancesA[j] += euclideanDistance(coordsCounter,vectors[j].coords,vectors[z].coords);
                distancesCounterA++;
            }
            ///if((j != z) && (secondNearest[z] == i) && (vectors[j].nearest == i))
            if((j != z) && (vectors[z].nearest == secondNearest[j])){ // distance with vectors in neighbour cluster
                distancesB[j] += euclideanDistance(coordsCounter,vectors[j].coords,vectors[z].coords);
                distancesCounterB++;
            }
        }
        distancesA[j] = distancesA[j]/distancesCounterA; //Average Distances
        distancesB[j] = distancesB[j]/distancesCounterB;
        if(distancesA[j] > distancesB[j])
            distancesA[j] = (distancesB[j]-distancesA[j])/distancesA[j];  //s(p)
        else
            distancesA[j] = (distancesB[j]-distancesA[j])/distancesB[j];  //s(p)
    }
    for(i = 0 ; i < vectorsCounter ; i++){
        silhouettes[vectors[i].nearest] += distancesA[i];
        *stotal += distancesA[i];
    }
    *stotal = *stotal/vectorsCounter;
    for(i = 0 ; i < numOfClusters ; i++){
        silhouettes[i] = silhouettes[i]/vectorsInCluster[i];
    }

    //frees

    free(secondNearest);
    free(distancesA);
    free(distancesB);
    free(vectorsInCluster);

    return silhouettes;
}