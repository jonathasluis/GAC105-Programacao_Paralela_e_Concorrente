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
        close(1);
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

void estastisticas(int n, int *nps, int clusters) {
    fstream out("output/estatistica_" + to_string(clusters) + ".tsv", ofstream::out);
    out << "processos\t";
    for (int i = 0; i < 10; i++) {
        out << "t" << i << "\t";
    }
    out << "media" << "\t";
    out << "speedup" << "\t";
    out << "efficiency" << "\t";
    out << "karpflatt" << endl;
    double t0;
    for (int i = 0; i < n; i++) {
        int np = nps[i];
        out << np << "\t";
        string command = "mpirun -np " + to_string(np) + " ./kmeans input/entrada.txt " + to_string(clusters);
        cout << command << endl;
        double media = 0.0;
        for (int j = 0; j < 10; j++) {
            double t = runGetTime(command);
            out << t << "\t";
            media += t;
        }
        if (i == 0) {
            t0 = media;
        }
        out << media / 10 << "\t";
        double speed = speedUp(t0, media);
        out << speed << "\t";
        out << efficiency(speed, np) << "\t";
        double kf = np == 1 ? 0 : karpFlat(speed, np); 
        out << kf << endl;
    }
    out.close();
}

int main() {
    system("mpicc -o kmeans kmeans.c");
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

    int nps[] = {1, 2, 4, 6, 8};
    estastisticas(5, nps, 2);
    estastisticas(5, nps, 3);
    estastisticas(5, nps, 4);
    return 0;
}
