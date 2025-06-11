#include <random>
#include "Ising.h"
using namespace std;

Ising::Ising(Grid *g_grid) : data(g_grid->volume), exp_delta_E_list(2 * (2 * g_grid->nn_batch_size + 1))
{
    grid = g_grid;
    for (int i = 0; i < grid->volume; i++)
    {
        data[i] = 1;
    }
    printf("Ising model initialized!\n");
}

Ising::Ising(Grid *g_grid, double g_J, double g_B) : data(g_grid->volume), exp_delta_E_list(2 * (2 * g_grid->nn_batch_size + 1))
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
            h_2 += g_data[grid->nearest_neighbors[i * grid->nn_batch_size + j]];
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

void Ising::get_proposal_prob(double &proposal_prob, int rand_index)
{
    int index;
    int dummy = 0;
    index = grid->nn_batch_size * (data[rand_index] + 1);
    for (size_t i = 0; i < grid->nn_batch_size; i++)
    {
        dummy += data[grid->nearest_neighbors[grid->nn_batch_size * rand_index + i]];
    }
    index += data[rand_index] * dummy + grid->nn_batch_size;
    proposal_prob = exp_delta_E_list[index];
}

void Ising::metropolis(mt19937 &gen, vector<double> &energies, vector<double> &magnetizations, double T, int N_mc, int N_therm)
{
    randomize(gen); // random inital state
    uniform_int_distribution<int> unidist(0, grid->volume - 1);
    uniform_real_distribution<double> unidist_real(0, 1);
    int rand_index;
    double proposal_prob;
    double rand_real;
    double accept_rate = 0;
    double beta = 1 / T;

    int contrib_J;
    int contrib_B;
    for (size_t i = 0; i < exp_delta_E_list.size(); i++)
    {
        contrib_J = i % (2 * grid->nn_batch_size + 1) - grid->nn_batch_size;
        contrib_B = i / (2 * grid->nn_batch_size + 1) * 2 - 1;
        exp_delta_E_list[i] = exp(-beta * (2 * J * contrib_J + 2 * B * contrib_B));
    }
    energies.resize(N_mc);
    magnetizations.resize(N_mc);

    for (size_t i = 0; i < N_therm + N_mc * grid->volume; i++)
    {
        rand_index = unidist(gen);
        get_proposal_prob(proposal_prob, rand_index);
        rand_real = unidist_real(gen);
        if (rand_real < proposal_prob)
        {
            accept_rate++;
            data[rand_index] *= -1;
        }
        if (i >= N_therm && i % grid->volume == 0)
        {
            energies[(i - N_therm) / grid->volume] = energy();
            magnetizations[(i - N_therm) / grid->volume] = magnetization();
        }
    }
    printf("Metropolis accept-rate = %.3lf%%\n", accept_rate / (N_therm + N_mc * grid->volume) * 100);
}
