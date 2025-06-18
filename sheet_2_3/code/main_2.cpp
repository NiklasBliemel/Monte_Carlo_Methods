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
    double L = 64;
    int N_mc = 90000;
    int N_therm = 10000;

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
    double var_auto_correlation;
    double int_auto_correlation_time;
    double var_auto_corr_std;

    ofstream wf;
    string filename;
    stringstream stream;

    for (double temp : temps)
    {
        stream.str("");
        stream << fixed << setprecision(1) << temp;
        printf("\nStarting run with temp : %.1lfK\n\n", temp);

        ising.metropolis(gen, energies, magnetizations, temp, N_mc, N_therm);

        filename = "../../plot_lab/auto_correlations/auto_corr_e_t" + stream.str();
        wf.open(filename, ios::binary);
        int_auto_correlation_time = auto_correlation_time(energies);
        var_auto_corr_std = auto_corr_std(energies) * ising.get_volume();
        wf.write(reinterpret_cast<char *>(&int_auto_correlation_time), sizeof(double));
        wf.write(reinterpret_cast<char *>(&var_auto_corr_std), sizeof(double));
        for (int t = 0; t < energies.size(); t++)
        {
            var_auto_correlation = auto_correlation(energies, t);
            if (var_auto_correlation < 0)
            {
                break;
            }
            wf.write(reinterpret_cast<char *>(&var_auto_correlation), sizeof(double));
        }
        wf.close();

        filename = "../../plot_lab/auto_correlations/auto_corr_m_t" + stream.str();
        wf.open(filename, ios::binary);
        int_auto_correlation_time = auto_correlation_time(magnetizations);
        var_auto_corr_std = auto_corr_std(magnetizations) * ising.get_volume();
        wf.write(reinterpret_cast<char *>(&int_auto_correlation_time), sizeof(double));
        wf.write(reinterpret_cast<char *>(&var_auto_corr_std), sizeof(double));
        for (int t = 0; t < energies.size(); t++)
        {
            var_auto_correlation = auto_correlation(magnetizations, t);
            if (var_auto_correlation < 0)
            {
                break;
            }
            wf.write(reinterpret_cast<char *>(&var_auto_correlation), sizeof(double));
        }
        wf.close();
    }

    printf("\nDone!!\n");

    return 0;
}
