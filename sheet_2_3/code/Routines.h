#ifndef Routines_H
#define Routines_H
#include <vector>
using namespace std;

double sqr(double x);

double vec_mean(vector<double> &vec, int offset);

double vec_var(vector<double> &vec, int offset);

double vec_std(vector<double> &vec, int offset);

double auto_correlation_time(vector<double> &vec);

#endif