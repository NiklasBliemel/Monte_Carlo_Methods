#include "Routines.h"
#include <vector>
#include <functional>

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
    return sqrt(2 * auto_correlation_time(vec, N_therm) / vec.size()) * vec_std(vec, N_therm);
}
