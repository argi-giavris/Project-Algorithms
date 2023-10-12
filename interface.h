#ifndef INTERFACE_H
#define INTERFACE_H

typedef struct v{
    int id;
    int *coords;
}vector;

typedef struct vc{
    int id;
    int *coords;
	int nearest;
	int secondNearest;
}vectorCluster;

struct h_func{
	int h_floor;
	float* v;
	int t;
	int r;
};

typedef struct ln{
	int g;
	int vec_pos;
	struct ln *next;
}list_node;

float inner(int *, float *, int );
void normalDistr(float* , int );
float euclideanDistance(int ,int* , int* );
void swap1(double* , double* );
void swap2(int* , int* );
void sortArrays(double *, int *, int );
void hash(list_node **, int , int , int );
void countFile(char* , int *, int *);
void saveFile(char* , vector *);
void query_knn(int , int, vector *, vector , double *, int, int*);
void writeResults(char*, int, vector *, vector *, int **, double **, double **, double *, double *, int** , int, int);
int eucmod(int, int);

#endif