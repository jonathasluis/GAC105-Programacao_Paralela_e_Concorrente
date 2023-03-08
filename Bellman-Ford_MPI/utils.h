using std::cout;
using std::endl;
using std::string;
#define INF 1000000

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
