#include <random>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

#include "Grid.h"
#include "XY.h"
#include "Routines.h"

using namespace std;
int main(int argc, char const *argv[])
{
    // make a unique seed based on the current time
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    printf("Seed = %u\n", seed);

    // define a mersenne twister generator
    mt19937 gen;

    // initialize generator with the seed
    gen.seed(seed);

    // configure XY model
    double J = 1;
    double L = 16;
    int N_mc = 100000;
    int N_therm_default = 10000;

    vector<double> temps(3);
    double T_0 = 0.1;
    double del_T = 0.1;
    int N_steps = 24;

    vector<int> shape(2);
    shape[0] = L;
    shape[1] = L;
    printf("Shape = [%d, %d]\n", shape[0], shape[1]);

    Grid grid(shape);
    XY xy(&grid, J);

    // -------- Start data production -------- //

    vector<double> energies;
    vector<double> magnetizations;
    double del_t = 0.1;
    int n = 4;

    printf("-------- Starting run with HMC from T_0 = %lf to T_N = %lf with Delta T = %lf --------\n", T_0, T_0 + del_T * N_steps, del_T);
    printf("T\t\te\t\tsigma_e\t\tm\t\tsigma_m\n");
    string filename = "/Users/niklasbliemel/Study/Monte_Carlo_Methods/plot_lab/sheet_5_data/ex_1.o";
    ofstream wf;
    wf.open(filename);
    for (size_t i = 0; i < N_steps + 1; i++)
    {
        // calculate temperature
        double T = T_0 + i * del_T;
        printf("%lf\t", T);
        wf << fixed << setprecision(10)
           << T << "\t";

        xy.hmc(gen, energies, magnetizations, T, N_mc, del_t, n);

        // autocorrelation for energies
        double int_auto_correlation_time = auto_correlation_time(energies, N_therm_default);
        int N_therm = (int)(20 * int_auto_correlation_time);

        // calculate observables on energies
        double e = vec_mean(energies, N_therm);
        double sigma_e = auto_corr_std(energies, N_therm);

        // autocorrelation for magnetizations
        int_auto_correlation_time = auto_correlation_time(magnetizations, N_therm_default);
        N_therm = 20 * (int)int_auto_correlation_time;

        // calculate observables on magnetizations
        double m = vec_mean(magnetizations, N_therm);
        double sigma_m = auto_corr_std(magnetizations, N_therm);

        printf("%.10lf\t%.10lf\t%.10lf\t%.10lf\n", e, sigma_e, m, sigma_m);
        wf << fixed << setprecision(10)
           << e << "\t"
           << sigma_e << "\t"
           << m << "\t"
           << sigma_m << "\n";
    }
    return 0;
}