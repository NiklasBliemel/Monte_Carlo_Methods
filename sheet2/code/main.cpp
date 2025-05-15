#include <iostream>
#include <chrono>
#include <random>
#include <vector>
#include "Grid.h"
#include "Ising.h"
using namespace std;

int main(int argc, char const *argv[])
{
    // make a unique seed based on the current time
    unsigned seed = 123;
    cout << "Seed: " << seed << std::endl;
    printf("\n");

    // define a mersenne twister generator
    mt19937 gen;

    // initialize generator with the seed
    gen.seed(seed);

    vector<int> shape(2);
    shape[0] = 4;
    shape[1] = 5;
    printf("Shape = [%d, %d]\n\n", shape[0], shape[1]);

    Grid grid(shape);
    Ising ising(&grid);
    int test_size = 1;
    for (int i = 0; i < test_size; i++)
    {
        ising.org(1);
        printf("Energy = %lf\n", ising.energy());
        printf("Magnetization = %lf\n", ising.magnetization());
    }
    return 0;
}
