#include "routines.h"
#include <vector>

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

double auto_correlation_time(vector<double> &vec) // TODO
{
    return 0;
}