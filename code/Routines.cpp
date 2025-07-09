#include "Routines.h"
#include <vector>
#include <functional>
#include <random>

using namespace std;

double sqr(double x)
{
    return x * x;
}

double vec_mean(vector<double> &vec, int N_therm)
{
    double out = 0;
    for (int i = N_therm; i < vec.size(); i++)
    {
        out += vec[i];
    }
    return out / (vec.size() - N_therm);
}

double vec_var(vector<double> &vec, int N_therm)
{
    double mean_value = vec_mean(vec, N_therm);
    double out = 0;
    for (int i = N_therm; i < vec.size(); i++)
    {
        out += sqr(mean_value - vec[i]);
    }
    return out / (vec.size() - N_therm);
}

double vec_std(vector<double> &vec, int N_therm)
{
    return sqrt(vec_var(vec, N_therm));
}

double vec_mean(vector<double> &vec, int start, int size)
{
    double out = 0;
    for (int i = start; i < size + start; i++)
    {
        out += vec[i];
    }
    return out / size;
}

double vec_var(vector<double> &vec, int start, int size)
{
    double mean_value = vec_mean(vec, start, size);
    double out = 0;
    for (int i = start; i < size + start; i++)
    {
        out += sqr(mean_value - vec[i]);
    }
    return out / size;
}

double vec_std(vector<double> &vec, int start, int size)
{
    return sqrt(vec_var(vec, start, size));
}

double blocking_std(vector<double> &blocks)
{
    return sqrt(vec_var(blocks, 0)) / sqrt(blocks.size());
}

double auto_covariance(vector<double> &vec, int t, int N_therm)
{
    int N = vec.size() - N_therm;
    double y_minus = 0;
    double y_plus = 0;
    for (size_t i = N_therm; i < N + N_therm - t; i++)
    {
        y_minus += vec[i];
        y_plus += vec[i + t];
    }
    y_minus /= N - t;
    y_plus /= N - t;

    double out = 0;
    for (size_t i = N_therm; i < N + N_therm - t; i++)
    {
        out += (vec[i] - y_plus) * (vec[i + t] - y_minus);
    }
    return out;
}

double auto_correlation(vector<double> &vec, int t, int N_therm)
{
    return auto_covariance(vec, t, N_therm) / auto_covariance(vec, 0, N_therm);
}

double auto_correlation_time(vector<double> &vec, int N_therm)
{
    double out = 0;
    double C_0 = auto_covariance(vec, 0, N_therm);
    double Ct = 0;
    for (int t = 1; t < vec.size() - N_therm; t++)
    {
        Ct = auto_covariance(vec, t, N_therm);
        if (C_0 * Ct < 0)
        {
            break;
        }
        out += Ct;
    }
    return 0.5 + out / C_0;
}

double auto_corr_std(vector<double> &vec, int N_therm)
{
    double tau = auto_correlation_time(vec, N_therm);
    int N_therm_new = (int)(20 * tau);
    return sqrt(2 * tau / (vec.size() - N_therm_new)) * vec_std(vec, N_therm_new);
}

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