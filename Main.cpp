#include <iostream>
#include <stdlib.h>
#include <fstream>

using namespace std;

double media(int np)
{
    string fileName = "output/timeOutput" + to_string(np) + ".txt";
    ifstream inputf(fileName, ifstream::in);
    double soma = 0.0;
    if (inputf.good())
    {
        for (int i = 0; i < 10; i++)
        {
            double valor = 0.0;
            inputf >> valor;
            soma += valor;
        }
    }
    else
    {
        cerr << "Erro com o arquivo";
        exit(0);
    }
    return soma / 10;
}

double speedUp(double tempSeq, double tempPara)
{
    return tempSeq / tempPara;
}

double efficiency(double speedUp, int np)
{
    return speedUp / np;
}

double karpFlat(double speedUp, int np)
{
    return ((1 / speedUp) - (1 / np)) / (1 - (1 / np));
}

void gravaMedias()
{
    string fileName = "output/timeMedias.txt";
    ofstream outfile(fileName, ios_base::out);
    double tempoSequencial = media(1);
    for (int i = 1; i <= 4; i++)
    {
        if (i == 3)
        {
            continue;
        }
        
        double speedup = speedUp(tempoSequencial,media(i));
        double eficiencia = efficiency(speedup,i);
        double e = karpFlat(speedup,i);
        outfile << i << " : " << media(i) << " - " << speedup<< " - " << eficiencia<< " - " << e << endl;
    }
}

int main(int argc, char const *argv[])
{

    for (int i = 0; i < 10; i++)
    {
        system("mpirun -np 1 ./teste input/matrix.txt");
    }

    for (int i = 0; i < 10; i++)
    {
        system("mpirun -np 2 ./teste input/matrix.txt");
    }

    for (int i = 0; i < 10; i++)
    {
        system("mpirun -np 4 ./teste input/matrix.txt");
    }

    gravaMedias();
}
