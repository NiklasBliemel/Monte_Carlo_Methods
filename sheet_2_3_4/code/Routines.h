#ifndef Routines_H
#define Routines_H
#include <vector>
#include <functional>
using namespace std;

double sqr(double x);

double vec_mean(vector<double> &vec, int N_therm);

double vec_var(vector<double> &vec, int N_therm);

double vec_std(vector<double> &vec, int N_therm);

double vec_mean(vector<double> &vec, int start, int size);

double vec_var(vector<double> &vec, int start, int size);

double vec_std(vector<double> &vec, int start, int size);

double blocking_std(vector<double> &blocks);

double auto_covariance(vector<double> &vec, int t, int N_therm);

double auto_correlation(vector<double> &vec, int t, int N_therm);

double auto_correlation_time(vector<double> &vec, int N_therm);

double auto_corr_std(vector<double> &vec, int N_therm);

#endif