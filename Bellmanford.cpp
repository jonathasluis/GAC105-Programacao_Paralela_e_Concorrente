/*
 * This is a mpi version of bellman_ford algorithm
 * Compile: mpic++ Bellmanford.cpp -o Bellmanford
 * Run: mpiexec -n <number of processes> ./Bellmanford matrix.txt
 * Or : mpirun -np <number of processes> ./Bellmanford matrix.txt
 * you will find the output file 'output.txt'
 * */

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "mpi.h"

using std::cout;
using std::endl;
using std::string;

#define INF 1000000

/**
 * utils is a namespace for utility functions
 * including I/O (read input file and print results) and matrix dimension
 * convert(2D->1D) function
 */
namespace utils {
int N;                 // number of vertices
int *adjacencyMatrix;  // the adjacency matrix

// translate 2-dimension coordinate to 1-dimension
int convert_dimension_2D_1D(int x, int y, int inputSize) {
    return x * inputSize + y;
}

void read_file(string filename) {
    std::ifstream inputf(filename, std::ifstream::in);
    if (inputf.good()) {
        inputf >> N;
        // input matrix should be smaller than 20MB * 20MB (400MB, we don't have
        // too much memory for multi-processors)
        adjacencyMatrix = (int *)malloc(N * N * sizeof(int));
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                inputf >> adjacencyMatrix[convert_dimension_2D_1D(i, j, N)];
            }
    } else {
        std::cerr << "Erro com o arquivo";
        abort();
    }
    inputf.close();
}

int print_result(bool has_negative_cycle, int *distanceArray) {
    std::ofstream outputf("output/output.txt", std::ofstream::out);
    if (!has_negative_cycle) {
        for (int i = 0; i < N; i++) {
            if (distanceArray[i] > INF) distanceArray[i] = INF;
            outputf << distanceArray[i] << '\n';
        }
        outputf.flush();
    } else {
        outputf << "FOUND NEGATIVE CYCLE!" << endl;
    }
    outputf.close();
    return 0;
}

int timeOutput(double t1, double t2, int np) {
    std::ofstream outfile;
    string fileName = "output/timeOutput" + std::to_string(np) + ".txt";
    outfile.open(fileName, std::ios_base::app);
    outfile << std::setprecision(6) << (t2 - t1) << endl;
    outfile.close();
    return 0;
}
}  // namespace utils

/**
 * Bellman-Ford algorithm. Find the shortest path from vertex 0 to other
 * vertices.
 */
void bellman_ford(int my_rank, int numberProcesses, MPI_Comm comm,
                  int inputSize, int *adjacencyMatrix, int *distanceArray,
                  bool *has_negative_cycle) {
    int loc_n;  // need a local copy for N
    int loc_start, loc_end;
    int *loc_mat;   // local matrix
    int *loc_dist;  // local distance

    // step 1: broadcast N
    if (my_rank == 0) {
        loc_n = inputSize;
    }
    MPI_Bcast(&loc_n, 1, MPI_INT, 0, comm);

    // step 2: find local task range
    int ave = loc_n / numberProcesses;
    loc_start = ave * my_rank;
    loc_end = ave * (my_rank + 1);
    if (my_rank == numberProcesses - 1) {
        loc_end = loc_n;
    }

    // step 3: allocate local memory
    loc_mat = (int *)malloc(loc_n * loc_n * sizeof(int));
    loc_dist = (int *)malloc(loc_n * sizeof(int));

    // step 4: broadcast matrix adjacencyMatrix
    if (my_rank == 0) {
        memcpy(loc_mat, adjacencyMatrix, sizeof(int) * loc_n * loc_n);
    }
    MPI_Bcast(loc_mat, loc_n * loc_n, MPI_INT, 0, comm);

    // step 5: bellman-ford algorithm
    for (int i = 0; i < loc_n; i++) {
        loc_dist[i] = INF;
    }
    loc_dist[0] = 0;
    MPI_Barrier(comm);

    bool loc_has_change;
    int loc_iter_num = 0;
    for (int iter = 0; iter < loc_n - 1; iter++) {
        loc_has_change = false;
        loc_iter_num++;
        for (int u = loc_start; u < loc_end; u++) {
            for (int v = 0; v < loc_n; v++) {
                int weight =
                    loc_mat[utils::convert_dimension_2D_1D(u, v, loc_n)];
                if (weight < INF) {
                    if (loc_dist[u] + weight < loc_dist[v]) {
                        loc_dist[v] = loc_dist[u] + weight;
                        loc_has_change = true;
                    }
                }
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &loc_has_change, 1, MPI_CXX_BOOL, MPI_LOR,
                      comm);
        if (!loc_has_change) break;
        MPI_Allreduce(MPI_IN_PLACE, loc_dist, loc_n, MPI_INT, MPI_MIN, comm);
    }

    // do one more step
    if (loc_iter_num == loc_n - 1) {
        loc_has_change = false;
        for (int u = loc_start; u < loc_end; u++) {
            for (int v = 0; v < loc_n; v++) {
                int weight =
                    loc_mat[utils::convert_dimension_2D_1D(u, v, loc_n)];
                if (weight < INF) {
                    if (loc_dist[u] + weight < loc_dist[v]) {
                        loc_dist[v] = loc_dist[u] + weight;
                        loc_has_change = true;
                        break;
                    }
                }
            }
        }
        MPI_Allreduce(&loc_has_change, has_negative_cycle, 1, MPI_CXX_BOOL,
                      MPI_LOR, comm);
    }

    // step 6: retrieve results back
    if (my_rank == 0) memcpy(distanceArray, loc_dist, loc_n * sizeof(int));

    // step 7: remember to free memory
    free(loc_mat);
    free(loc_dist);
}

int main(int argc, char **argv) {
    if (argc <= 1) {
        std::cerr << "Nenhum arquivo passado como parametro";
    }
    string filename = argv[1];

    int *distanceArray;
    bool has_negative_cycle = false;

    // MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm comm;

    int numberProcesses;  // number of processors
    int my_rank;          // my global rank
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &numberProcesses);
    MPI_Comm_rank(comm, &my_rank);

    // only rank 0 process do the I/O
    if (my_rank == 0) {
        utils::read_file(filename);
        distanceArray = (int *)malloc(sizeof(int) * utils::N);
    }

    // time counter
    double t1, t2;
    MPI_Barrier(comm);
    t1 = MPI_Wtime();

    // bellman-ford algorithm
    bellman_ford(my_rank, numberProcesses, comm, utils::N,
                 utils::adjacencyMatrix, distanceArray, &has_negative_cycle);
    MPI_Barrier(comm);

    // end timer
    t2 = MPI_Wtime();

    if (my_rank == 0) {
        std::cout << utils::N << endl;
        std::cerr.setf(std::ios::fixed);
        std::cerr << std::setprecision(6) << "Time(s): " << (t2 - t1) << endl;
        utils::timeOutput(t1, t2, numberProcesses);
        utils::print_result(has_negative_cycle, distanceArray);
        free(distanceArray);
        free(utils::adjacencyMatrix);
    }
    MPI_Finalize();
    return 0;
}