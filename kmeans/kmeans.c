// mpicc -o kmeans kmeans.c

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <locale.h>

int loadVertices(char *arq, double *v[2], int *n) {
    FILE *f = fopen(arq, "r");
    if (f == NULL) {
        return 1;
    }

    int lines = 1;
    while (1) {
        int c = fgetc(f);
        if (c == EOF) break;
        if (c == '\n') lines += 1;
    }

    v[0] = malloc(lines * sizeof(double));
    v[1] = malloc(lines * sizeof(double));

    rewind(f);
    
    for (*n = 0; ; *n += 1) {
        int res = fscanf(f, "%lf %lf", &v[0][*n], &v[1][*n]);
        if (res == EOF) {
            break;
        } else if (res != 2) {
            return 2;
        }
    }

    v[0] = realloc(v[0], *n * sizeof(double));
    v[1] = realloc(v[1], *n * sizeof(double));

    fseek(f, 0, SEEK_SET);
    return 0;
}

double square(double a) {
    return a * a;
}

double calcJ(int n, double *v[2], int k, double mu[2][k]) {
    double J = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            J += square(v[0][i] - mu[0][j]) + square(v[1][i] - mu[1][j]);
        }
    }
    return J;
}

int closestCentroid(int i, int n, double *v[2], int k, double mu[2][k]) {
    double dist = INFINITY;
    int vk = 0;
    for (int j = 0; j < k; j++) {
        double d = square(v[0][i] - mu[0][j]) + square(v[1][i] - mu[1][j]);
        if (d < dist) {
            vk = j;
            dist = d;
        }
    }
    return vk;
}

int nextPermutation(int max, int n, int permutation[n]) {
    int c = 0;
    int finish = 0;
    permutation[0] += 1;
    for (int i = 0; i < n; i++) {
        permutation[i] += c;
        if (permutation[i] == max) {
            permutation[i] = 0;
            c = 1;
            if (i == n-1) finish = 1;
        }
    }
    return !finish;
}

int main(int argc, char *argv[]) {
    setlocale(LC_ALL, "C");

    if (argc < 3) {
        fprintf(stderr, "USO: %s ARQUIVO K\n", argv[0]);
        return 1;
    }

    char *arq = argv[1];
    int k = strtol(argv[2], NULL, 10);

    double *allv[2];
    int n;

    if (k <= 0) {
        fprintf(stderr, "valor inválido para K. Garanta que K >= 0 (valor dado %d)\n", k);
    }

    if (loadVertices(arq, allv, &n)) {
        fprintf(stderr, "Falha ao carregar %s\n", arq);
        return 2;
    }

    int size;
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double t_start = MPI_Wtime();

    srand(1);

    int nprime = n / size;

    double *v[2];
    v[0] = malloc(nprime * sizeof(double));
    v[1] = malloc(nprime * sizeof(double));

    MPI_Scatter(allv[0], nprime, MPI_DOUBLE, v[0], nprime, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(allv[1], nprime, MPI_DOUBLE, v[1], nprime, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double threshold = 1;
    double mu[2][k];

    for (int j = 0; j < k; j++) {
        int i = rand() % nprime;
        mu[0][j] = v[0][i];
        mu[1][j] = v[1][i];
    }

    double J = calcJ(nprime, v, k, mu);
    double Jprime = 0;
    int iter = 0;
    int maxiter = 1000;
    while (fabs(Jprime - J) > threshold && iter < maxiter) {
        Jprime = J;

        double muprime[2][k];
        int musize[k];
        for (int j = 0; j < k; j++) {
            muprime[0][j] = 0;
            muprime[1][j] = 0;
            musize[j] = 0;
        }

        for (int i = 0; i < nprime; i++) {
            int j = closestCentroid(i, nprime, v, k, mu);
            muprime[0][j] += v[0][i];
            muprime[1][j] += v[1][i];
            musize[j] += 1;
        }

        for (int j = 0; j < k; j++) {
            mu[0][j] = muprime[0][j] / musize[j];
            mu[1][j] = muprime[1][j] / musize[j];
        }

        J = calcJ(nprime, v, k, mu);
        iter++;
    }

    double mus[size][2][k];
    MPI_Gather(mu, 2*k, MPI_DOUBLE, mus, 2*k, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double finalmu[2][k];
        for (int j = 0; j < k; j++) {
            finalmu[0][j] = INFINITY;
            finalmu[1][j] = INFINITY;
        }

        for (int muidx = 0; muidx < k; muidx++) {
            int permutation[size];
            for (int i = 0; i < size; i++) {
                permutation[i] = 0;
            }

            int minpermutation[size];
            double mindistance = INFINITY;
            double mincenter[2];

            do {
                double center[2] = {0, 0};
                for (int i = 0; i < size; i++) {
                    center[0] += mus[i][0][permutation[i]];
                    center[1] += mus[i][1][permutation[i]];
                }
                if (center[0] == INFINITY || center[1] == INFINITY) continue;

                center[0] /= size;
                center[1] /= size;

                double D = 0;
                for (int i = 0; i < size; i++) {
                    D += square(mus[i][0][permutation[i]] - center[0]) + square(mus[i][1][permutation[i]] - center[1]);
                }

                if (D < mindistance) {
                    for (int i = 0; i < size; i++) {
                        minpermutation[i] = permutation[i];
                    }
                    mindistance = D;
                    mincenter[0] = center[0];
                    mincenter[1] = center[1];
                }
            } while (nextPermutation(k, size, permutation));

            for (int i = 0; i < size; i++) {
                mus[i][0][minpermutation[i]] = INFINITY;
                mus[i][1][minpermutation[i]] = INFINITY;
            }

            finalmu[0][muidx] = mincenter[0];
            finalmu[1][muidx] = mincenter[1];
        }

        for (int j = 0; j < k; j++) {
            printf("{ %lf, %lf }", finalmu[0][j], finalmu[1][j]);
            if (j < k - 1) {
                printf(",");
            }
            printf("\n");
        }
        
        double t_end = MPI_Wtime();
        fprintf(stderr, "%lf\n", t_end - t_start);
    }

    MPI_Finalize();
    return 0;
}
