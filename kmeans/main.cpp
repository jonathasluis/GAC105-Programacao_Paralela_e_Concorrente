#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <time.h>
#include <unistd.h>
#include <sys/wait.h>

using namespace std;


double speedUp(double tempSeq, double tempPara) { return tempSeq / tempPara; }

double efficiency(double speedUp, int np) { return speedUp / np; }

double karpFlat(double speedUp, int np) { return ((1.0 / speedUp) - (1.0 / np)) / (1.0 - (1.0 / np)); }

double runGetTime(string cmd) {
    double t = 0;

    int fd_err[2];
    pipe(fd_err);

    pid_t p;
    if ((p = fork()) == 0) {
        close(fd_err[0]);
        dup2(fd_err[1], 2);

        cout << cmd << endl;
        system(cmd.c_str());
        exit(0);
    } else {
        int status;
        close(fd_err[1]);
        waitpid(p, &status, 0);
        char buf[64];
        read(fd_err[0], buf, 64);
        t = strtod(buf, NULL);
    }
    return t;
}

double media_tempo(int np) {
    string command = "mpirun -np " + to_string(np) + " ./kmeans input/entrada.txt 4";
    time_t start, end;
    string fileName = "output/times_np_" + to_string(np) + ".txt";
    fstream outputf(fileName, ofstream::out);
    double media = 0.0;
    for (int i = 0; i < 10; i++) {
        double t = runGetTime(command);
        cout << t << endl;
        outputf << t << ';';
        media += t;
    }
    outputf << endl << "Media: " << media / 10 << endl;
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

    // gera arquivo de entrada com 1000000 pontos aleatorios
    ofstream inputf("input/entrada.txt", ofstream::out);
    if (inputf.good()) {
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(0, 1000000);
        for (int i = 0; i < 1000000; i++) {
            inputf << dis(gen) << " " << dis(gen) << endl;
        }
        inputf.close();
    }

    double tempo_sequencial = media_tempo(1);
    double tempo_np2 = media_tempo(2);
    double tempo_np4 = media_tempo(4);
    double tempo_np8 = media_tempo(8);

    calcula_metricas(tempo_sequencial, tempo_np2, 2);
    calcula_metricas(tempo_sequencial, tempo_np4, 4); 
    calcula_metricas(tempo_sequencial, tempo_np8, 8);

    return 0;
}
