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
    vector<double> cluster_densities;
    double T_0 = 2.1;
    double del_T = 0.01;
    int N_steps = 30;
    printf("-------- Starting run with Metropolis and Wolff from T_0 = %lf to T_N = %lf with Delta T = %lf --------\n", T_0, T_0 + del_T * N_steps, del_T);
    printf("T\t\t<n>/V\t\tsigma_n\t\ttau_e_wolff\ttau_m_wolff\ttau_e_metro\ttau_m_metro\n");
    string filename = "../../plot_lab/sheet_3_data/ex_3_4.o";
    ofstream wf;
    wf.open(filename);
    for (size_t i = 0; i < N_steps + 1; i++)
    {
        double tau_g;

        // calculate temperature
        double T = T_0 + i * del_T;
        printf("%.10lf\t", T);
        wf << T << "\t";

        // generate sample with wolff
        ising.wolff_cluster(gen, energies, magnetizations, cluster_densities, T, N_mc);

        // autocorrelation for energies
        double int_auto_correlation_time = auto_correlation_time(cluster_densities, N_therm_default);
        int N_therm = (int)(20 * int_auto_correlation_time);

        // cluster density
        double n = vec_mean(cluster_densities, N_therm);
        printf("%.10lf\t", n);
        wf << n << "\t";

        // cluster density error
        double sigma_n = auto_corr_std(cluster_densities, N_therm);
        printf("%.10lf\t", sigma_n);
        wf << sigma_n << "\t";

        // autocorrelation for energies
        int_auto_correlation_time = auto_correlation_time(energies, N_therm_default) * n;
        printf("%.10lf\t", int_auto_correlation_time);
        wf << int_auto_correlation_time << "\t";

        // autocorrelation for magnetizations
        int_auto_correlation_time = auto_correlation_time(magnetizations, N_therm_default) * n;
        printf("%.10lf\t", int_auto_correlation_time);
        wf << int_auto_correlation_time << "\t";

        // generate sample with metropolis
        ising.metropolis(gen, energies, magnetizations, T, N_mc, 0);

        // autocorrelation for energies
        int_auto_correlation_time = auto_correlation_time(energies, N_therm_default);
        printf("%.10lf\t", int_auto_correlation_time);
        wf << int_auto_correlation_time << "\t";

        // autocorrelation for magnetizations
        int_auto_correlation_time = auto_correlation_time(magnetizations, N_therm_default);
        printf("%.10lf\n", int_auto_correlation_time);
        wf << int_auto_correlation_time << "\n";
    }
    wf.close();
    return 0;
}