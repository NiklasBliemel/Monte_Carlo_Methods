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

    // define a mersenne twister generator
    mt19937 gen;

    // initialize generator with the seed
    gen.seed(seed);

    // configure Ising model
    double J = 1;
    double B = 0;
    double L = 64;
    int N_mc = 100000;
    int N_therm_default = 10000;

    vector<double> temps(3);
    temps[0] = 2.0;
    temps[1] = 2.3;
    temps[2] = 2.6;

    vector<int> shape(2);
    shape[0] = L;
    shape[1] = L;
    printf("Shape = [%d, %d]\n", shape[0], shape[1]);

    Grid grid(shape);
    Ising ising(&grid, J, B);

    // -------- Start data production -------- //

    vector<double> energies;
    vector<double> magnetizations;
    double T_0 = 1;
    double del_T = 0.1;
    int N_steps = 30;

    printf("-------- Starting run with Metropolis from T_0 = %lf to T_N = %lf with Delta T = %lf --------\n", T_0, T_0 + del_T * N_steps, del_T);
    printf("T\t\te\t\tsigma_e\t\tc\t\tsigma_c\t\tm\t\tsigma_m\t\txi\t\tsigma_xi\n");

    string filename = "/Users/niklasbliemel/Study/Monte_Carlo_Methods/plot_lab/sheet_4_data/metropolis.o";
    ofstream wf;
    wf.open(filename);
    for (size_t i = 0; i < N_steps + 1; i++)
    {
        double tau_g;

        // calculate temperature
        double T = T_0 + i * del_T;

        printf("%lf\t", T);
        wf << T << "\t";

        // generate sample
        ising.metropolis(gen, energies, magnetizations, T, N_mc, 0);

        // autocorrelation for energies
        double int_auto_correlation_time = auto_correlation_time(energies, N_therm_default);
        int N_therm = 20 * (int)int_auto_correlation_time;

        // calculate observables on energies
        double e = vec_mean(energies, N_therm);
        double sigma_e = auto_corr_std(energies, N_therm);
        double c = vec_var(energies, N_therm) * ising.get_volume() / sqr(T);
        double sigma_c = error_propagation(energies, N_therm, tau_g) * ising.get_volume() / sqr(T);

        // autocorrelation for magnetizations
        int_auto_correlation_time = auto_correlation_time(magnetizations, N_therm_default);
        N_therm = 20 * (int)int_auto_correlation_time;

        // calculate observables on magnetizations
        double m = vec_mean(magnetizations, N_therm);
        double sigma_m = auto_corr_std(magnetizations, N_therm);
        double xi = vec_var(magnetizations, N_therm) * ising.get_volume() / T;
        double sigma_xi = error_propagation(magnetizations, N_therm, tau_g) * ising.get_volume() / T;

        printf("%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", e, sigma_e, c, sigma_c, m, sigma_m, xi, sigma_xi);
        wf << fixed << setprecision(10)
           << e << "\t"
           << sigma_e << "\t"
           << c << "\t"
           << sigma_c << "\t"
           << m << "\t"
           << sigma_m << "\t"
           << xi << "\t"
           << sigma_xi << "\n";
    }
    wf.close();
    return 0;
}