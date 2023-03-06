/*
 * This is a mpi version of bellman_ford algorithm
 * Compile: mpic++ -std=c++11 -o teste2 teste2.cpp
 * Run: mpiexec -n <number of processes> ./teste2 matrix.txt, you will find the output file 'output.txt'
 * */

#include <string>
#include <cassert>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <cstring>

#include "mpi.h"

using std::string;
using std::cout;
using std::endl;

#define INF 1000000

/**
 * utils is a namespace for utility functions
 * including I/O (read input file and print results) and matrix dimension convert(2D->1D) function
 */
namespace utils {
    int N; //number of vertices
    int *adjacencyMatrix; // the adjacency matrix

    //translate 2-dimension coordinate to 1-dimension
    int convert_dimension_2D_1D(int x, int y, int inputSize) {
        return x * inputSize + y;
    }

    void read_file(string filename) {
        std::ifstream inputf(filename, std::ifstream::in);
        if (inputf.good()) {
            inputf >> N;
            //input matrix should be smaller than 20MB * 20MB (400MB, we don't have too much memory for multi-processors)
            assert(N < (1024 * 1024 * 20));
            adjacencyMatrix = (int *) malloc(N * N * sizeof(int));
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++) {
                    inputf >> adjacencyMatrix[convert_dimension_2D_1D(i, j, N)];
                }
        } else {
            std::cerr << "Erro com o arquivo";
            exit(0);
        }
    }

    int print_result(bool has_negative_cycle, int *distanceArray) {
        std::ofstream outputf("output/output.txt", std::ofstream::out);
        if (!has_negative_cycle) {
            for (int i = 0; i < N; i++) {
                if (distanceArray[i] > INF)
                    distanceArray[i] = INF;
                outputf << distanceArray[i] << '\n';
            }
            outputf.flush();
        } else {
            outputf << "FOUND NEGATIVE CYCLE!" << endl;
        }
        outputf.close();
        return 0;
    }
}//namespace utils

void free_Mem(int *loc_mat, int *loc_dist, int *loc_enqueue_counter, bool *queue1, bool *queue2);

/**
 * Bellman-Ford algorithm. Find the shortest path from vertex 0 to other vertices.
*/
void bellman_ford(int my_rank, int numberProcesses, MPI_Comm comm, int inputSize, int *adjacencyMatrix, int *distanceArray, bool *has_negative_cycle) {
    int loc_n; // need a local copy for N
    int loc_start, loc_end;
    int *loc_mat; //local matrix
    int *loc_dist; //local distance

    bool *queue1;
    bool *queue2;
    int *loc_enqueue_counter;

    //step 1: broadcast N
    if (my_rank == 0) {
        loc_n = inputSize;
    }
    MPI_Bcast(&loc_n, 1, MPI_INT, 0, comm);

    //step 2: find local task range
    int ave = loc_n / numberProcesses;
    loc_start = ave * my_rank;
    loc_end = ave * (my_rank + 1);
    if (my_rank == numberProcesses - 1) {
        loc_end = loc_n;
    }

    //step 3: allocate local memory
    loc_mat = (int *) malloc(loc_n * loc_n * sizeof(int));
    loc_dist = (int *) malloc(loc_n * sizeof(int));
    queue1 = (bool *) malloc(loc_n * sizeof(bool));
    queue2 = (bool *) malloc(loc_n * sizeof(bool));
    loc_enqueue_counter = (int *) calloc(loc_n, sizeof(int));

    //step 4: broadcast matrix adjacencyMatrix
    if (my_rank == 0)
        memcpy(loc_mat, adjacencyMatrix, sizeof(int) * loc_n * loc_n);
    MPI_Bcast(loc_mat, loc_n * loc_n, MPI_INT, 0, comm);

    //step 5: bellman-ford algorithm
    for (int i = 0; i < loc_n; i++) {
        loc_dist[i] = INF;
    }
    loc_dist[0] = 0;
    queue1[0] = true;

    bool loc_has_change;
    bool has_change;
    auto iter = 0;
    while(1){
        loc_has_change = false;
        memset(queue2, 0, sizeof(bool) * loc_n);// 0 as false
        for(int u = 0; u < loc_n; u++){
            if(queue1[u]){ // if u is active
                for(int v = loc_start; v < loc_end; v++){
                    int weight = loc_mat[utils::convert_dimension_2D_1D(u, v, loc_n)];
                    if (weight < INF) {
                        if (loc_dist[u] + weight < loc_dist[v]) {
                            loc_dist[v] = loc_dist[u] + weight;
                            queue2[v] = true;
                            loc_has_change = true;
                        }
                    }
                }
            }
        }
        MPI_Allreduce(queue2, queue1, loc_n, MPI_CXX_BOOL, MPI_LOR, comm);
        for(int u = 0 ; u < loc_n; u++){
            if(queue1[u]){
                loc_enqueue_counter[u]++;
                if(loc_enqueue_counter[u] >= loc_n) {
                    free_Mem(loc_mat,loc_dist,loc_enqueue_counter, queue1,queue2);
                }
            }
        }
        
            MPI_Allreduce(&loc_has_change, &has_change, 1, MPI_CXX_BOOL, MPI_LOR, comm);
            if(!has_change){
                std::cout << iter << std::endl;
                break;
            }

            MPI_Allreduce(MPI_IN_PLACE, loc_dist, loc_n, MPI_INT, MPI_MIN, comm);
            ++iter;
        
    }


    //step 6: retrieve results back
    if(my_rank == 0)
        memcpy(distanceArray, loc_dist, loc_n * sizeof(int));

    //step 7: remember to free memory
    free_Mem(loc_mat,loc_dist,loc_enqueue_counter, queue1,queue2);
}

void free_Mem(int *loc_mat, int *loc_dist, int *loc_enqueue_counter, bool *queue1, bool *queue2){
    free(loc_mat);
    free(loc_dist);
    free(queue1);
    free(queue2);
    free(loc_enqueue_counter); 
}

int timeOutput(double t1, double t2) {  
  std::ofstream outfile;
  outfile.open("output/timeOutput.txt", std::ios_base::app);
  outfile << std::setprecision(6) << (t2 - t1) << endl; 
  return 0;
}

int main(int argc, char **argv) {
    if (argc <= 1) {
        std::cerr << "Nenhum arquivo passo como parametro";
    }
    string filename = argv[1];

    int *distanceArray;
    bool has_negative_cycle = false;

    //MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm comm;

    int numberProcesses;//number of processors
    int my_rank;//my global rank
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &numberProcesses);
    MPI_Comm_rank(comm, &my_rank);

    //only rank 0 process do the I/O
    if (my_rank == 0) {
        utils::read_file(filename);
        distanceArray = (int *) malloc(sizeof(int) * utils::N);
    }

    //time counter
    double t1, t2;
    MPI_Barrier(comm);
    t1 = MPI_Wtime();

    //bellman-ford algorithm
    bellman_ford(my_rank, numberProcesses, comm, utils::N, utils::adjacencyMatrix, distanceArray, &has_negative_cycle);
    MPI_Barrier(comm);

    //end timer
    t2 = MPI_Wtime();

    if (my_rank == 0) {
        std::cout << utils::N << endl;
        std::cerr.setf(std::ios::fixed);
        std::cerr << std::setprecision(6) << "Time(s): " << (t2 - t1) << endl;
        timeOutput(t1,t2);
        utils::print_result(has_negative_cycle, distanceArray);
        free(distanceArray);
        free(utils::adjacencyMatrix);
    }
    MPI_Finalize();
    return 0;
}