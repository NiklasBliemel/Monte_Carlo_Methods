#include "grid.h"
#include <vector>
using namespace std;
#define mod(a, b) ((((a) % (b)) + b) % b)

Grid::Grid(vector<int> &g_data, vector<int> &g_shape) : data(g_data)
{
    data = g_data;
    shape = g_shape;
    stride.resize(g_shape.size());
    stride[stride.size() - 1] = 1;
    for (int i = stride.size() - 2; i >= 0; i--)
    {
        stride[i] = stride[i + 1] * shape[i + 1];
    }
}

void Grid::print()
{
    printf("[%d", data[0]);
    for (int i = 1; i < data.size(); i++)
    {
        printf(", %d", data[i]);
    }
    printf("]\n");
}

vector<int> Grid::nn(int index)
{
    vector<int> out(2 * shape.size());
    int temp;
    for (int i = 0; i < shape.size(); i++)
    {
        temp = stride[i] * shape[i];
        out[2 * i] = mod((index + stride[i]), temp) + index / temp * temp;
        out[2 * i + 1] = mod((index - stride[i]), temp) + index / temp * temp;
    }
    return out;
}
