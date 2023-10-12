all: cluster search

cluster: cluster.c methodsCluster.c implementationCluster.c interfaceCluster.h
	gcc -o cluster cluster.c methodsCluster.c implementationCluster.c -lm -g -O2

search: search.c interface.c lsh.c hyperCube.c frechet.c interface.h
	gcc -o search search.c lsh.c hyperCube.c interface.c frechet.c -lm -g

clean:
	rm -r search cluster