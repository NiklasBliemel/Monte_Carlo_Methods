#ifndef XY_H
#define XY_H
#include <random>
#include "Grid.h"
using namespace std;

class XY
{
private:
    double J = 1;
    Grid *grid;
    vector<double> data;
    double dV_dqj(int j, double T);
    double H_hmc(vector<double> &momenta, double T);

public:
    XY(Grid *g_grid);
    XY(Grid *g_grid, double g_J);
    double energy();
    double magnetization();
    void org(double state);
    void randomize(mt19937 &gen);
    int get_volume();
    void hmc(mt19937 &gen, vector<double> &energies, vector<double> &magnetizations, double T, int N_mc, double del_t, int n);
};

#endif