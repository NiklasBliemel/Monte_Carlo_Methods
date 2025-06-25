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

double error_propagation(vector<double> &sample, int N_therm, double &tau_g)
{
    vector<double> gs(sample.size());
    vector<double> gs_alt(sample.size());
    double mean = vec_mean(sample, N_therm);
    for (size_t i = 0; i < sample.size(); i++)
    {
        gs[i] = sample[i] * (sample[i] - 2 * mean);
    }
    tau_g = auto_correlation_time(gs, N_therm);
    int N_therm_g = (int)(20 * tau_g);
    return sqrt(2 * tau_g / (gs.size() - N_therm_g)) * vec_std(gs, N_therm_g);
}

double blocking(vector<double> &sample, int n_b, int N_therm)
{
    vector<double> blocks(n_b);
    int block_size = (sample.size() - N_therm) / n_b;
    for (int i = 0; i < n_b; i++)
    {
        blocks[i] = vec_var(sample, i * block_size + N_therm, block_size);
    }
    return blocking_std(blocks);
}

double bootstrap(vector<double> &sample, int M, int N_therm, double tau_g, mt19937 gen)
{
    uniform_int_distribution<int> unidist(N_therm, sample.size() - 1);
    vector<double> boots(M);
    vector<double> pseudo_sample((sample.size() - N_therm) / (int)(2 * tau_g));
    for (size_t i = 0; i < M; i++)
    {
        for (size_t j = 0; j < (sample.size() - N_therm) / (int)(2 * tau_g); j++)
        {
            int rand_index = unidist(gen);
            pseudo_sample[j] = sample[rand_index];
        }
        boots[i] = vec_var(pseudo_sample, 0);
    }
    return vec_std(boots, 0);
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

    int n_b = 20;
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
        printf("\nenergies and specific heat:\n\n");

        // auto_correlations
        double int_auto_correlation_time = auto_correlation_time(energies, N_therm_default);
        int N_therm = 20 * (int)int_auto_correlation_time;
        double sigma_e = auto_corr_std(energies, N_therm) * ising.get_volume() / sqr(temp);
        double c = vec_var(energies, N_therm) * ising.get_volume() / sqr(temp);

        printf("Auto correlation time %lf\n", int_auto_correlation_time);
        printf("N_therm = %d\n\n", N_therm);
        printf("c = %lf\n\n", c);

        string filename = "../../plot_lab/auto_correlations/auto_corr_e_t" + stream.str();
        wf.open(filename, ios::binary);
        wf.write(reinterpret_cast<char *>(&int_auto_correlation_time), sizeof(double));
        wf.write(reinterpret_cast<char *>(&sigma_e), sizeof(double));
        for (int t = 0; t < energies.size(); t++)
        {
            double var_auto_correlation = auto_correlation(energies, t, N_therm);
            if (var_auto_correlation < 0)
            {
                break;
            }
            wf.write(reinterpret_cast<char *>(&var_auto_correlation), sizeof(double));
        }
        wf.close();

        // error propagation
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();

        double tau_g;
        double sigma_c = error_propagation(energies, N_therm, tau_g) * ising.get_volume() / sqr(temp);

        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        printf("Error propagation in %lld µs:\n", chrono::duration_cast<chrono::microseconds>(end - begin).count());
        printf("sigma_c = %lf\n\n", sigma_c);

        // blocking
        begin = chrono::steady_clock::now();

        sigma_c = blocking(energies, n_b, N_therm) * ising.get_volume() / sqr(temp);

        end = chrono::steady_clock::now();
        printf("Blocking in %lld µs:\n", chrono::duration_cast<chrono::microseconds>(end - begin).count());
        printf("sigma_c = %lf\n\n", sigma_c);

        // -------- bootstrap (c and xi) -------- //
        begin = chrono::steady_clock::now();

        sigma_c = bootstrap(energies, M, N_therm, tau_g, gen) * ising.get_volume() / sqr(temp);

        end = chrono::steady_clock::now();
        printf("Bootstrap in %lld µs:\n", chrono::duration_cast<chrono::microseconds>(end - begin).count());
        printf("sigma_c = %lf\n\n", sigma_c);

        // ----- magnetization and susceptibility ----- //
        printf("\nenergies and specific heat:\n\n");

        // auto_correlations
        int_auto_correlation_time = auto_correlation_time(magnetizations, N_therm_default);
        N_therm = 20 * (int)int_auto_correlation_time;
        double sigma_m = auto_corr_std(magnetizations, N_therm) * ising.get_volume() / temp;
        double xi = vec_var(magnetizations, N_therm) * ising.get_volume() / temp;

        printf("Auto correlation time %lf\n", int_auto_correlation_time);
        printf("N_therm = %d\n\n", N_therm);
        printf("xi = %lf\n\n", xi);

        filename = "../../plot_lab/auto_correlations/auto_corr_m_t" + stream.str();
        wf.open(filename, ios::binary);
        wf.write(reinterpret_cast<char *>(&int_auto_correlation_time), sizeof(double));
        wf.write(reinterpret_cast<char *>(&sigma_m), sizeof(double));
        for (int t = 0; t < energies.size(); t++)
        {
            double var_auto_correlation = auto_correlation(magnetizations, t, N_therm);
            if (var_auto_correlation < 0)
            {
                break;
            }
            wf.write(reinterpret_cast<char *>(&var_auto_correlation), sizeof(double));
        }
        wf.close();

        // error propagation
        begin = chrono::steady_clock::now();

        double sigma_xi = error_propagation(magnetizations, N_therm, tau_g) * ising.get_volume() / temp;

        end = chrono::steady_clock::now();
        printf("Error propagation in %lld µs:\n", chrono::duration_cast<chrono::microseconds>(end - begin).count());
        printf("sigma_Xi = %lf\n\n", sigma_xi);

        // blocking
        begin = chrono::steady_clock::now();

        sigma_xi = blocking(magnetizations, n_b, N_therm) * ising.get_volume() / temp;

        end = chrono::steady_clock::now();
        printf("Blocking in %lld µs:\n", chrono::duration_cast<chrono::microseconds>(end - begin).count());
        printf("sigma_xi = %lf\n\n", sigma_xi);

        // -------- bootstrap (c and xi) -------- //
        begin = chrono::steady_clock::now();

        sigma_xi = bootstrap(magnetizations, M, N_therm, tau_g, gen) * ising.get_volume() / temp;

        end = chrono::steady_clock::now();
        printf("Bootstrap in %lld µs:\n", chrono::duration_cast<chrono::microseconds>(end - begin).count());
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
    double true_sigma_c = vec_std(true_cs, 0);
    double true_sigma_xi = vec_std(true_xis, 0);

    printf("\ntrue c = %lf\n", true_c);
    printf("true sigma_c = %lf\n", true_sigma_c);
    printf("true xi = %lf\n", true_xi);
    printf("true sigma_xi = %lf\n\n", true_sigma_xi); */
    return 0;
}