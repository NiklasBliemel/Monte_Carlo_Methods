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

public:
    Ising(Grid *grid);
    double energy();
    double magnetization();
    void org(int state);
    void randomize(mt19937 &gen);
};

#endif