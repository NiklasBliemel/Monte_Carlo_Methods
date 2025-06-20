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

double g_c(double beta, int volume, double e_mean, double e)
{
    return sqr(beta) * volume * e * (e - 2 * e_mean);
}

double g_xi(double beta, int volume, double m_mean, double m)
{
    return beta * volume * m * (m - 2 * m_mean);
}

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
    int N_therm_default = 10000;

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

    vector<double> g_cs;
    vector<double> g_xis;

    ofstream wf;
    string filename;
    stringstream stream;

    for (double temp : temps)
    {
        g_cs.clear();
        g_xis.clear();
        stream.str("");
        stream << fixed << setprecision(1) << temp;
        printf("\nStarting run with temp : %.1lfK\n\n", temp);
        ising.metropolis(gen, energies, magnetizations, temp, N_mc, 0);

        // Energie and Specific Heat
        filename = "../../plot_lab/auto_correlations/auto_corr_e_t" + stream.str();
        wf.open(filename, ios::binary);

        double int_auto_correlation_time = auto_correlation_time(energies, N_therm_default);
        int N_therm = 20 * (int)int_auto_correlation_time;
        double var_auto_corr_std = auto_corr_std(energies, N_therm) * ising.get_volume();

        wf.write(reinterpret_cast<char *>(&int_auto_correlation_time), sizeof(double));
        wf.write(reinterpret_cast<char *>(&var_auto_corr_std), sizeof(double));
        for (int t = 0; t < energies.size(); t++)
        {
            double var_auto_correlation = auto_correlation(energies, t, N_therm);
            if (var_auto_correlation < 0)
            {
                break;
            }
            wf.write(reinterpret_cast<char *>(&var_auto_correlation), sizeof(double));
        }
        double mean_e = vec_mean(energies, N_therm);
        for (double e : energies)
        {
            g_cs.push_back(g_c(1 / temp, ising.get_volume(), mean_e, e));
        }
        double c = ising.get_volume() * vec_var(energies, N_therm) / sqr(temp);
        double sigma_c = auto_corr_std(g_cs, N_therm);
        printf("c = %lf\n", c);
        printf("sigma_c = %lf\n", sigma_c);
        wf.close();

        // Magnetization and Susceptibility
        filename = "../../plot_lab/auto_correlations/auto_corr_m_t" + stream.str();
        wf.open(filename, ios::binary);

        int_auto_correlation_time = auto_correlation_time(magnetizations, N_therm_default);
        N_therm = 20 * (int)int_auto_correlation_time;
        var_auto_corr_std = auto_corr_std(magnetizations, N_therm) * ising.get_volume();

        wf.write(reinterpret_cast<char *>(&int_auto_correlation_time), sizeof(double));
        wf.write(reinterpret_cast<char *>(&var_auto_corr_std), sizeof(double));
        for (int t = 0; t < energies.size(); t++)
        {
            double var_auto_correlation = auto_correlation(magnetizations, t, N_therm);
            if (var_auto_correlation < 0)
            {
                break;
            }
            wf.write(reinterpret_cast<char *>(&var_auto_correlation), sizeof(double));
        }
        double mean_m = vec_mean(magnetizations, N_therm);
        for (double m : magnetizations)
        {
            g_xis.push_back(g_c(1 / temp, ising.get_volume(), mean_m, m));
        }
        double xi = ising.get_volume() * vec_var(magnetizations, N_therm) / temp;
        double sigma_xi = auto_corr_std(g_xis, N_therm);
        printf("Xi = %lf\n", xi);
        printf("sigma_Xi = %lf\n", sigma_xi);
        wf.close();
    }

    printf("\nDone!!\n");

    return 0;
}