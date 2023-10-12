#ifndef INTERFACECLUSTER_H
#define INTERFACECLUSTER_H

typedef struct vc{
    char* id;
    float *coords;
	int nearest;
}vectorCluster;

struct h_func{
	int h_floor;
	float* v;
	int t;
	int r;
    float* tFrechet;
};

typedef struct ln{
	int g;
	int vec_pos;
	struct ln *next;
}list_node;

typedef struct tn{
	float* data;
	int pos;
	struct tn* left;
	struct tn* right;
}treeNode;


void uniformDistr(float*, int, double, double);
float inner(float *, float *, int );
void normalDistr(float* , int );
float euclideanDistance(int ,float* , float* );
void swap1(double* , double* );
void swap2(int* , int* );
void sortArrays(double *, int *, int );
void hash(list_node **, int , int , int );
void countFile(char* , int *, int *);
void saveFile(char* , vectorCluster *);
void writeResults(char*, int, int, vectorCluster *, int, int, vectorCluster*, double, int,float*,float*,int);
int eucmod(int, int);
int hammingDistance(int, int);
int vectorIsCentroid(vectorCluster* , char*, int);
void kmeansInit(vectorCluster*, int, int, int, vectorCluster*);
void lloydAssignment(vectorCluster *, vectorCluster *, int, int, int);
void lshTrain(vectorCluster *, struct h_func **, list_node ***, int, int, int, int, int, int, int);
void lshSearch(vectorCluster* , vectorCluster*, struct h_func **, list_node ***, int, int, int, int, int, int, int, int);
void cubeTrain(list_node**, vectorCluster *, struct h_func *, list_node**, int, int, int,int);
void cubeSearch(vectorCluster*, vectorCluster*, struct h_func *, list_node **, list_node**, int, int, int, int, int, int, int);
int update(vectorCluster*, vectorCluster*, int, int, int);
float* silhouette(int, int, vectorCluster*, vectorCluster*, int, float*);

#endif