#include "Routines.h"
#include <vector>

using namespace std;

double sqr(double x)
{
    return x * x;
}

double vec_mean(vector<double> &vec, int offset)
{
    double out = 0;
    for (int i = offset; i < vec.size(); i++)
    {
        out += vec[i];
    }
    return out / (vec.size() - offset);
}

double vec_var(vector<double> &vec, int offset)
{
    double mean_value = vec_mean(vec, offset);
    double out = 0;
    for (int i = offset; i < vec.size(); i++)
    {
        out += sqr(mean_value - vec[i]);
    }
    return out / (vec.size() - 1 - offset);
}

double vec_std(vector<double> &vec, int offset)
{
    return sqrt(vec_var(vec, offset));
}

double auto_covariance(vector<double> &vec, int t)
{
    int N = vec.size();
    double y_minus = 0;
    double y_plus = 0;
    for (size_t i = 0; i < N - t; i++)
    {
        y_minus += vec[i];
        y_plus += vec[i + t]
    }
    y_minus /= N - t;
    y_plus /= N - t;

    double out = 0;
    for (size_t i = 0; i < N - t; i++)
    {
        out += (vec[i] - y_plus) * (vec[i + t] - y_minus);
    }
    return out;
}

double auto_correlation(vector<double> &vec, double t)
{
    return auto_covariance(t) / auto_covariance(0);
}

double auto_correlation_time(vector<double> &vec) // TODO
{
    double out = 0.5;
}