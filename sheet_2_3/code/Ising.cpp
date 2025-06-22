#include <random>
#include <map>
#include "Ising.h"
using namespace std;

Ising::Ising(Grid *g_grid) : data(g_grid->volume), exp_delta_E_list(g_grid->nn_batch_size + 1)
{
    grid = g_grid;
    for (int i = 0; i < grid->volume; i++)
    {
        data[i] = 1;
    }
    printf("Ising model initialized!\n");
}

Ising::Ising(Grid *g_grid, double g_J, double g_B) : data(g_grid->volume), exp_delta_E_list(g_grid->nn_batch_size + 1)
{
    grid = g_grid;
    J = g_J;
    B = g_B;
    for (int i = 0; i < grid->volume; i++)
    {
        data[i] = 1;
    }
    printf("Ising model initialized!\n");
}

double Ising::energy()
{
    double h_1 = 0;
    double h_2 = 0;
    for (int i = 0; i < grid->volume; i++)
    {
        h_1 += data[i];
        double dummy = 0;
        for (int j = 0; j < grid->nn_batch_size; j += 2)
        {
            dummy += data[grid->nearest_neighbors[i * grid->nn_batch_size + j]];
        }
        h_2 += data[i] * dummy;
    }
    return (-J * h_2 - B * h_1) / grid->volume;
}

double Ising::energy(vector<int> g_data)
{
    if (g_data.size() != grid->volume)
    {
        printf("Given data size %zu does not match grid volume %d\n", g_data.size(), grid->volume);
        return 0;
    }

    double h_1 = 0;
    double dummy;
    double h_2 = 0;
    for (int i = 0; i < grid->volume; i++)
    {
        h_1 += g_data[i];
        dummy = 0;
        for (int j = 0; j < grid->nn_batch_size; j += 2)
        {
            dummy += g_data[grid->nearest_neighbors[i * grid->nn_batch_size + j]];
        }
        h_2 += g_data[i] * dummy;
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
    return abs(out / grid->volume);
}

double Ising::magnetization(vector<int> g_data)
{
    if (g_data.size() != grid->volume)
    {
        printf("Given data size %zu does not match grid volume %d\n", g_data.size(), grid->volume);
        return 0;
    }

    double out = 0;
    for (int i = 0; i < grid->volume; i++)
    {
        out += g_data[i];
    }
    return abs(out / grid->volume);
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

int Ising::get_exp_list_index(int index)
{
    int out = 0;
    for (int i = 0; i < grid->nn_batch_size; i++)
    {
        out += data[grid->nearest_neighbors[index * grid->nn_batch_size + i]];
    }
    return (out * data[index] + grid->nn_batch_size) / 2;
}

void Ising::metropolis(mt19937 &gen, vector<double> &energies, vector<double> &magnetizations, double T, int N_mc, int N_therm)
{
    uniform_real_distribution<double> unidist_real(0, 1);
    double accept_rate = 0;
    double beta = 1 / T;
    printf("Start Metropolis T = %.2lf; N_cm = %d; N_therm = %d\n", T, N_mc, N_therm);

    for (unsigned i = 0; i < exp_delta_E_list.size(); i++)
    {
        exp_delta_E_list[i] = exp(-beta * (-2 * J * grid->nn_batch_size + 4 * J * i));
    }

    energies.clear();
    magnetizations.clear();

    for (unsigned i = 0; i < N_therm + N_mc; i++)
    {
        for (unsigned j = 0; j < grid->volume; j++)
        {
            double proposal_prob = exp_delta_E_list[get_exp_list_index(j)];
            if (proposal_prob > 1)
            {
                accept_rate++;
                data[j] *= -1;
            }
            else if (unidist_real(gen) < proposal_prob)
            {
                accept_rate++;
                data[j] *= -1;
            }
        }
        if (i > N_therm)
        {
            energies.push_back(energy());
            magnetizations.push_back(magnetization());
        }
    }
    printf("Metropolis accept-rate = %.3lf%%\n", accept_rate / ((N_therm + N_mc) * grid->volume) * 100);
}

int Ising::get_volume()
{
    return grid->volume;
}
