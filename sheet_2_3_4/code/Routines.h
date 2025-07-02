#ifndef Routines_H
#define Routines_H
#include <vector>
#include <functional>
#include <random>

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

double error_propagation(vector<double> &sample, int N_therm, double &tau_g);

double blocking(vector<double> &sample, int n_b, int N_therm);

double bootstrap(vector<double> &sample, int M, int N_therm, double tau_g, mt19937 gen);

#endif