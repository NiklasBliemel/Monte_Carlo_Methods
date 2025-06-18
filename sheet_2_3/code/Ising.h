#ifndef Ising_H
#define Ising_H
#include <random>
#include "Grid.h"
using namespace std;

class Ising
{
private:
    double J = 1;
    double B = 0;
    Grid *grid;
    vector<int> data;
    vector<double> exp_delta_E_list;
    int get_exp_list_index(int index);

public:
    Ising(Grid *g_grid);
    Ising(Grid *g_grid, double g_J, double g_B);
    double energy();
    double energy(vector<int> g_data);
    double magnetization();
    double magnetization(vector<int> g_data);
    void org(int state);
    void randomize(mt19937 &gen);
    void metropolis(mt19937 &gen, vector<double> &energies, vector<double> &magnetizations, double T, int N_mc, int N_therm);
    int get_volume();
};

#endif