#include <random>
#include "Ising.h"
using namespace std;

Ising::Ising(Grid *g_grid) : data(g_grid->volume)
{
    grid = g_grid;
    for (int i = 0; i < grid->volume; i++)
    {
        data[i] = 1;
    }
}

double Ising::energy()
{
    double h_1 = 0;
    double dummy;
    double h_2 = 0;
    for (int i = 0; i < grid->volume; i++)
    {
        h_1 += data[i];
        dummy = 0;
        for (int j = 0; j < grid->nn_batch_size; j += 2)
        {
            h_2 += data[grid->nearest_neighbors[i * grid->nn_batch_size + j]];
        }
        h_2 += data[i] * dummy;
    }
    return (-J * h_2 - B * h_1) / grid->volume;
}

double Ising::magnetization()
{
    double out = 0;
    for (int i = 0; i < grid->volume; i++)
    {
        out += data[i];
    }
    return out / grid->volume;
}

void Ising::org(int state)
{
    for (int i = 0; i < grid->volume; i++)
    {
        data[i] = state;
    }
}

void Ising::randomize(mt19937 &gen)
{
    // define a uniform distribution for double between 0 and 1
    uniform_real_distribution<double> unidist(-1, 1);
    double random_double;
    for (int i = 0; i < grid->volume; i++)
    {
        random_double = unidist(gen);
        if (random_double >= 0)
        {
            data[i] = 1;
        }
        else
        {
            data[i] = -1;
        }
    }
}
