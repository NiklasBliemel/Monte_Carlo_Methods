#include <random>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <chrono>

#include "Grid.h"
#include "Ising.h"
#include "Routines.h"

using namespace std;

int main(int argc, char const *argv[])
{
    // make a unique seed based on the current time
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    printf("Seed = %u\n", seed);
    printf("\n");

    // define a mersenne twister generator
    mt19937 gen;

    // initialize generator with the seed
    gen.seed(seed);

    // configure Ising model
    double J = 1;
    double B = 0;
    double L = 64;
    int N_mc = 100000;
    int N_therm = 0;
    int chains_per_T = 3;
    double test = 0;

    vector<double> temps(3);
    temps[0] = 2.0;
    temps[1] = 2.3;
    temps[2] = 2.6;

    vector<int> shape(2);
    shape[0] = L;
    shape[1] = L;
    printf("Shape = [%d, %d]\n\n", shape[0], shape[1]);

    Grid grid(shape);
    Ising ising(&grid, J, B);

    vector<double> energies;
    vector<double> magnetizations;
    double c;
    double chi;

    ofstream wf;
    string filename;
    stringstream stream;

    if (test)
    {
        ising.metropolis(gen, energies, magnetizations, 2.6, N_mc, N_therm);
    }
    else
    {
        for (double temp : temps)
        {
            stream.str("");
            stream << fixed << setprecision(1) << temp;
            printf("\nStarting run with temp : %.1lfK\n\n", temp);
            for (int i = 0; i < chains_per_T; i++)
            {
                chrono::steady_clock::time_point begin = chrono::steady_clock::now();
                printf("\nStarting subrun %d\n", i);
                ising.metropolis(gen, energies, magnetizations, temp, N_mc, N_therm);
                c = vec_var(energies, 10000) / sqr(temp) * sqr(L);
                chi = vec_var(magnetizations, 10000) / temp * sqr(L);
                printf("specific heat c = %.8lf\n", c);
                printf("susceptibility chi = %.8lf\n", chi);

                filename = "/Users/niklasbliemel/Study/Monte_Carlo_Methods/plot_lab/energies/energies_t" + stream.str() + "_" + to_string(i);
                wf.open(filename, ios::binary);
                for (double energy : energies)
                {
                    wf.write(reinterpret_cast<char *>(&energy), sizeof(double));
                }
                wf.write(reinterpret_cast<char *>(&c), sizeof(double));
                wf.close();
                filename = "/Users/niklasbliemel/Study/Monte_Carlo_Methods/plot_lab/magnetizations/magnetizations_t" + stream.str() + "_" + to_string(i);
                wf.open(filename, ios::binary);
                for (double magnetization : magnetizations)
                {
                    wf.write(reinterpret_cast<char *>(&magnetization), sizeof(double));
                }
                wf.write(reinterpret_cast<char *>(&chi), sizeof(double));
                wf.close();
                chrono::steady_clock::time_point end = chrono::steady_clock::now();
                cout << "Time difference = " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << "[µs]" << endl;
                cout << "Time difference = " << chrono::duration_cast<chrono::nanoseconds>(end - begin).count() << "[ns]" << endl;
            }
        }
    }
    return 0;
}
