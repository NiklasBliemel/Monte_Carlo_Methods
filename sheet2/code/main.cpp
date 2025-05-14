#include "grid.h"
#include <vector>
using namespace std;
#define mod(a, b) ((((a) % (b)) + b) % b)

void print_vec(vector<int> vec)
{
    printf("[%d", vec[0]);
    for (int i = 1; i < vec.size(); i++)
    {
        printf(", %d", vec[i]);
    }
    printf("]\n");
}

int main(int argc, char const *argv[])
{
    vector<int> data(20);
    vector<int> shape(2);
    for (int i = 0; i < data.size(); i++)
    {
        data[i] = i;
    }
    shape[0] = 4;
    shape[1] = 5;

    Grid grid(data, shape);
    print_vec(grid.nn(12));
    return 0;
}
