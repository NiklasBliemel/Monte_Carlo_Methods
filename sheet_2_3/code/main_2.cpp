#include <random>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

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
    double L = 8;
    int N_mc = 100000;
    int N_therm = 1000;
    int chains_per_T = 3;

    vector<int> shape(2);
    shape[0] = L;
    shape[1] = L;
    printf("Shape = [%d, %d]\n\n", shape[0], shape[1]);

    Grid grid(shape);
    Ising ising(&grid, J, B);

    vector<double> energies;
    vector<double> magnetizations;

    ising.metropolis(gen, energies, magnetizations, 2.0, N_mc, N_therm);

    double tau_int = auto_correlation_time(energies);

    printf("auto correlation time test: %.4d\n", tau_int);

    return 0;
}
