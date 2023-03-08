#include <iostream>
#include <fstream>
#include <string>
#include <random>

using namespace std;


double speedUp(double tempSeq, double tempPara) { return tempSeq / tempPara; }

double efficiency(double speedUp, int np) { return speedUp / np; }

double karpFlat(double speedUp, int np) { return ((1 / speedUp) - (1.0 / np)) / (1 - (1.0 / np)); }

double media_tempo(int np) {
    string command = "mpirun -np " + to_string(np) + " ./kmeans input/entrada.txt 2";
    time_t start, end;
    string fileName = "output/times_np_" + to_string(np) + ".txt";
    fstream outputf(fileName, ofstream::out);
    double media = 0.0;
    for (int i = 0; i < 10; i++) {
        start = time(NULL);
        system(command.c_str());
        end = time(NULL);
        outputf << difftime(end, start) << endl;
        media += difftime(end, start);
    }
    outputf << "Media: " << media / 10 << endl;
    outputf.close();
    return media / 10;
}

void calcula_metricas(double tempo_sequencial, double tempo_npx, int np){
    string fileName = "output/metricas_np_" + to_string(np) + ".txt";
    ofstream outputf(fileName, ofstream::out);
    outputf.seekp(0, ios::end);
    if (outputf.good()) {
        outputf << "tempo;speedup;eficiencia;kf" << endl;
        double speedup = speedUp(tempo_sequencial, tempo_npx);
        double eficiencia = efficiency(speedup, np);
        double e = karpFlat(speedup, np);
        outputf << tempo_npx << ";" << speedup << ";" << eficiencia << ";" << e << endl;
        outputf.close();
    }
}

int main() {

    // gera arquivo de entrada com 100000 pontos aleatorios
    ofstream inputf("input/entrada.txt", ofstream::out);
    if (inputf.good()) {
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(0, 100000);
        for (int i = 0; i < 100000; i++) {
            inputf << dis(gen) << " " << dis(gen) << endl;
        }
        inputf.close();
    }

    double tempo_sequencial = media_tempo(1);
    double tempo_np2 = media_tempo(2);
    double tempo_np4 = media_tempo(4);

    calcula_metricas(tempo_sequencial, tempo_np2, 2);
    calcula_metricas(tempo_sequencial, tempo_np4, 4); 

    return 0;
}
