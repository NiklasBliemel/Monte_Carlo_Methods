#include <random>
#include <map>
#include "XY.h"
#include "Routines.h"
using namespace std;

XY::XY(Grid *g_grid) : data(g_grid->volume)
{
    grid = g_grid;
    for (int i = 0; i < grid->volume; i++)
    {
        data[i] = 0;
    }
    printf("XY model initialized!\n");
}

XY::XY(Grid *g_grid, double g_J) : data(g_grid->volume)
{
    grid = g_grid;
    J = g_J;
    for (int i = 0; i < grid->volume; i++)
    {
        data[i] = 0;
    }
    printf("XY model initialized!\n");
}

double XY::energy()
{
    double h = 0;
    for (int i = 0; i < grid->volume; i++)
    {
        for (int j = 0; j < grid->nn_batch_size; j += 2)
        {
            int nn = grid->nearest_neighbors[i * grid->nn_batch_size + j];
            h += -J * cos(data[i] - data[nn]);
        }
    }
    return h / grid->volume;
}

double XY::magnetization()
{
    double out_x = 0;
    double out_y = 0;
    for (int i = 0; i < grid->volume; i++)
    {
        out_x += cos(data[i]);
        out_y += sin(data[i]);
    }
    return sqrt(sqr(out_x) + sqr(out_y)) / grid->volume;
}

void XY::org(double state)
{
    for (int i = 0; i < grid->volume; i++)
    {
        data[i] = state;
    }
}

void XY::randomize(mt19937 &gen)
{
    // define a uniform distribution for double between 0 and 1
    uniform_real_distribution<double> unidist(0, 2 * M_PI);
    double random_double;
    for (int i = 0; i < grid->volume; i++)
    {
        data[i] = unidist(gen);
    }
}

int XY::get_volume()
{
    return grid->volume;
}

double XY::dV_dqj(int j, double T)
{
    double out = 0;
    for (size_t k = 0; k < grid->nn_batch_size; k++)
    {
        int nn = grid->nearest_neighbors[j * grid->nn_batch_size + k];
        out += -J * sin(data[j] - data[nn]) / T;
    }
    return out;
}

double XY::H_hmc(vector<double> &momenta, double T)
{
    double out = 0;
    for (size_t i = 0; i < momenta.size(); i++)
    {
        out += sqr(momenta[i]);
    }
    return out / 2 + energy() * grid->volume / T;
}

void XY::hmc(mt19937 &gen, vector<double> &energies, vector<double> &magnetizations, double T, int N_mc, double del_t, int n)
{
    uniform_real_distribution<double> unidist(0, 1);
    energies.clear();
    magnetizations.clear();

    vector<double> old_data = data;
    normal_distribution<double> normaldist(0.0, 1.0);
    vector<double> momenta(grid->volume);

    for (int k = 0; k < N_mc; k++)
    {
        for (size_t j = 0; j < momenta.size(); j++) // t = 0
        {
            momenta[j] = normaldist(gen);
        }

        double del_H = H_hmc(momenta, T);

        for (size_t j = 0; j < momenta.size(); j++) // t = del_t / 2
        {
            momenta[j] += -dV_dqj(j, T) * del_t / 2;
        }

        for (size_t i = 1; i < n; i++) // t = del_t, 3/2 del_t, ..., (n - 1) del_t, (n - 1/2) del_t
        {
            for (size_t j = 0; j < data.size(); j++)
            {
                data[j] += momenta[j] * del_t;
            }

            for (size_t j = 0; j < momenta.size(); j++)
            {
                momenta[j] += -dV_dqj(j, T) * del_t;
            }
        }

        for (size_t j = 0; j < data.size(); j++) // t = n del_t
        {
            data[j] += momenta[j] * del_t;
        }

        for (size_t j = 0; j < momenta.size(); j++) // t = n del_t
        {
            momenta[j] += -dV_dqj(j, T) * del_t / 2;
        }

        del_H -= H_hmc(momenta, T);

        if (-del_H <= 0)
        {
            double r = unidist(gen);
            if (r > exp(-del_H))
            {
                data = old_data;
            }
        }
        energies.push_back(energy());
        magnetizations.push_back(magnetization());
    }
}
