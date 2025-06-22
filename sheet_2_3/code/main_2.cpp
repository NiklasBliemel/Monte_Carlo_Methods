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

void blocking_c(vector<double> &source_energies, double T, int volume, vector<double> &block_c, int n_b, int N_therm)
{
    block_c.resize(n_b);
    int block_size = (source_energies.size() - N_therm) / n_b;
    for (int i = 0; i < n_b; i++)
    {
        block_c[i] = vec_var(source_energies, i * block_size + N_therm, block_size) * volume / sqr(T);
    }
}

void blocking_xi(vector<double> &source_magnetizations, double T, int volume, vector<double> &block_xi, int n_b, int N_therm)
{
    block_xi.resize(n_b);
    int block_size = (source_magnetizations.size() - N_therm) / n_b;
    for (int i = 0; i < n_b; i++)
    {
        block_xi[i] = volume * vec_var(source_magnetizations, i * block_size + N_therm, block_size) / T;
    }
}

void bootstrap_e(vector<double> &source_vector, double T, int volume, vector<double> &boots, int M, int N_therm, mt19937 gen)
{
    boots.resize(M);
    uniform_int_distribution<int> unidist(N_therm, volume - 1);
    vector<double> boot_sample(volume - N_therm);
    for (size_t i = 0; i < M; i++)
    {
        for (size_t j = 0; j < volume - N_therm; j++)
        {
            int rand_index = unidist(gen);
            boot_sample[j] = source_vector[rand_index];
        }
        boots[i] = volume * vec_var(boot_sample, 0) / sqr(T);
    }
}

void bootstrap_m(vector<double> &source_vector, double T, int volume, vector<double> &boots, int M, int N_therm, mt19937 gen)
{
    boots.resize(M);
    uniform_int_distribution<int> unidist(N_therm, volume - 1);
    vector<double> boot_sample(volume - N_therm);
    for (size_t i = 0; i < M; i++)
    {
        for (size_t j = 0; j < volume - N_therm; j++)
        {
            int rand_index = unidist(gen);
            boot_sample[j] = source_vector[rand_index];
        }
        boots[i] = volume * vec_var(boot_sample, 0) / T;
    }
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

    // -------- Start data production -------- //

    vector<double> energies;
    vector<double> magnetizations;

    vector<double> gs;
    vector<double> blocks;
    int n_b = 20;
    vector<double> boots;
    int M = 1000;

    ofstream wf;
    stringstream stream;

    for (double temp : temps)
    {
        stream.str("");
        stream << fixed << setprecision(1) << temp;
        printf("\n-------- Starting run with temp : %.1lfK --------\n\n", temp);

        ising.metropolis(gen, energies, magnetizations, temp, N_mc, 0);

        // ----- energies and specific heat ----- //

        // auto_correlations
        double int_auto_correlation_time = auto_correlation_time(energies, N_therm_default);
        printf("Auto correlation time %lf\n", int_auto_correlation_time);
        int N_therm = 20 * (int)int_auto_correlation_time;
        double var_auto_corr_std = auto_corr_std(energies, N_therm) * ising.get_volume();
        printf("N_therm = %d\n\n", N_therm);

        /* string filename = "../../plot_lab/auto_correlations/auto_corr_e_t" + stream.str();
        wf.open(filename, ios::binary);
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
        wf.close(); */

        // error propagation
        gs.clear();
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        double mean_e = vec_mean(energies, N_therm);
        for (double e : energies)
        {
            gs.push_back(g_c(1 / temp, ising.get_volume(), mean_e, e));
        }
        double c = ising.get_volume() * vec_var(energies, N_therm) / sqr(temp);
        double sigma_c = auto_corr_std(gs, N_therm);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        printf("Error propagation in %lld µs:\n", chrono::duration_cast<chrono::microseconds>(end - begin).count());
        printf("c = %lf\n", c);
        printf("sigma_c = %lf\n\n", sigma_c);

        // blocking
        blocks.clear();
        begin = chrono::steady_clock::now();
        blocking_c(energies, temp, ising.get_volume(), blocks, n_b, N_therm);
        c = vec_mean(blocks, 0);
        sigma_c = vec_std(blocks, 0);
        end = chrono::steady_clock::now();
        printf("Blocking in %lld µs:\n", chrono::duration_cast<chrono::microseconds>(end - begin).count());
        printf("c = %lf\n", c);
        printf("sigma_c = %lf\n\n", sigma_c);

        // -------- bootstrap (c and xi) -------- //
        boots.clear();
        begin = chrono::steady_clock::now();
        bootstrap_e(energies, temp, ising.get_volume(), boots, M, N_therm, gen);
        c = vec_mean(boots, 0);
        sigma_c = vec_std(boots, 0);
        end = chrono::steady_clock::now();
        printf("Bootstrap in %lld µs:\n", chrono::duration_cast<chrono::microseconds>(end - begin).count());
        printf("c = %lf\n", c);
        printf("sigma_c = %lf\n\n", sigma_c);

        // ----- magnetization and susceptibility ----- //

        // auto_correlations
        int_auto_correlation_time = auto_correlation_time(magnetizations, N_therm_default);
        N_therm = 20 * (int)int_auto_correlation_time;
        var_auto_corr_std = auto_corr_std(magnetizations, N_therm) * ising.get_volume();

        /* filename = "../../plot_lab/auto_correlations/auto_corr_m_t" + stream.str();
        wf.open(filename, ios::binary);
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
        wf.close(); */

        // error propagation
        gs.clear();
        begin = chrono::steady_clock::now();
        double mean_m = vec_mean(magnetizations, N_therm);
        for (double m : magnetizations)
        {
            gs.push_back(g_c(1 / temp, ising.get_volume(), mean_m, m));
        }
        double xi = ising.get_volume() * vec_var(magnetizations, N_therm) / temp;
        double sigma_xi = auto_corr_std(gs, N_therm);
        end = chrono::steady_clock::now();
        printf("Error propagation in %lld µs:\n", chrono::duration_cast<chrono::microseconds>(end - begin).count());
        printf("Xi = %lf\n", xi);
        printf("sigma_Xi = %lf\n\n", sigma_xi);

        // blocking
        blocks.clear();
        begin = chrono::steady_clock::now();
        blocking_xi(magnetizations, temp, ising.get_volume(), blocks, n_b, N_therm);
        xi = vec_mean(blocks, 0);
        sigma_xi = vec_std(blocks, 0);
        end = chrono::steady_clock::now();
        printf("Blocking in %lld µs:\n", chrono::duration_cast<chrono::microseconds>(end - begin).count());
        printf("xi = %lf\n", xi);
        printf("sigma_xi = %lf\n\n", sigma_xi);

        // -------- bootstrap (c and xi) -------- //
        boots.clear();
        begin = chrono::steady_clock::now();
        bootstrap_m(magnetizations, temp, ising.get_volume(), boots, M, N_therm, gen);
        xi = vec_mean(boots, 0);
        sigma_xi = vec_std(boots, 0);
        end = chrono::steady_clock::now();
        printf("Bootstrap in %lld µs:\n", chrono::duration_cast<chrono::microseconds>(end - begin).count());
        printf("xi = %lf\n", xi);
        printf("sigma_xi = %lf\n\n", sigma_xi);
    }

    // ---------- True standard error ---------- //
    /* M = 1000;
    printf("\n\n-------- True Std. for T = 2.3K --------\n\n");
    vector<double> true_cs(M);
    vector<double> true_xis(M);
    for (size_t i = 0; i < M; i++)
    {
        printf("%lf%%\n", (double)i / (double)M * 100);
        ising.metropolis(gen, energies, magnetizations, 2.3, N_mc, 0);
        double int_auto_correlation_time = auto_correlation_time(energies, N_therm_default);
        int N_therm = 20 * (int)int_auto_correlation_time;
        true_cs[i] = ising.get_volume() * vec_var(energies, N_therm) / sqr(2.3);
        int_auto_correlation_time = auto_correlation_time(magnetizations, N_therm_default);
        N_therm = 20 * (int)int_auto_correlation_time;
        true_xis[i] = ising.get_volume() * vec_var(magnetizations, N_therm) / 2.3;
    }
    double true_c = vec_mean(true_cs, 0);
    double true_xi = vec_mean(true_xis, 0);
    double true_sigma_c = vec_var(true_cs, 0);
    double true_sigma_xi = vec_var(true_xis, 0);

    printf("\ntrue c = %lf\n", true_c);
    printf("true sigma_c = %lf\n", true_sigma_c);
    printf("true xi = %lf\n", true_xi);
    printf("true sigma_xi = %lf\n\n", true_sigma_xi); */
    return 0;
}