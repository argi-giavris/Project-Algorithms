#ifndef INTERFACE_H
#define INTERFACE_H

typedef struct v{
    char* id;
    float *coords;
}vector;

typedef struct vc{
    char* id;
    float *coords;
	int nearest;
	int secondNearest;
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

void frechetContinuous(vector *, vector *, int, int, int, int, int, int, int, int, char*, int, char*, int, double);
void frechetContinuousSearch(vector*, vector, struct h_func *, list_node **, int, int, int, int, int, int, int, double*, int *, int, double);
void frechetSearch(vector*, vector, struct h_func **, list_node ***, int, int, int, int, int, int, int, double*, int *, int, int*, int, int*, double);
void frechet(vector *, vector *, int, int, int, int, int, int, int, int, char*, int, char*, int, double);
float discreteFrectet(vector, vector, int);
float recursiveC(float**,vector, vector, int, int);
float minFloat(float, float, float);
float maxFloat(float, float);
void frechetContinuousTrain(vector *, struct h_func *, list_node **, int, int, int, int, int, int, int, double);
void frechetTrain(vector *, struct h_func **, list_node ***, int, int, int, int, int, int, int, double);
vector* snapContinuous(vector*, double,int, int, float*, int*);
vector* snap(vector*, double,int, int, float*);
int hammingDistance(int, int);
void cubeTrain(list_node**, vector *, struct h_func *, list_node**, int, int, int,int);
void cubeSearch(vector*, vector, struct h_func *, list_node **, list_node**, int, int, int, int, int, int, double*, int *, int, int*, int, int*);
void cube(vector *, vector *, int, int, int, int, int, int, int, int, int, char*, int, char*, int);
void lshTrain(vector *, struct h_func **, list_node ***, int, int, int, int, int, int, int);
void lshSearch(vector*, vector, struct h_func **, list_node ***, int, int, int, int, int, int, int, double*, int *, int, int*, int, int*);
void lsh(vector *, vector *, int, int, int, int, int, int, int, int, char*, int, char*, int);
float inner(float *, float *, int );
void normalDistr(float* , int );
float euclideanDistance(int ,float* , float* );
void swap1(double* , double* );
void swap2(int* , int* );
void sortArrays(double *, int *, int );
void hash(list_node **, int , int , int );
void countFile(char* , int *, int *);
void saveFile(char* , vector *);
void query_knn(int , int, vector *, vector , double *, int, int*);
void writeResults(char*, int, vector *, vector *, int **, double **, double **, double *, double *, int** , int, int);
int eucmod(int, int);
void uniformDistr(float*, int, double, double);

#endif