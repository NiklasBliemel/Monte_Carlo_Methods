#include "Grid.h"
#include <vector>
using namespace std;
#define mod(a, b) ((((a) % (b)) + b) % b)

Grid::Grid(vector<int> &g_shape) : stride(g_shape.size())
{
    shape = g_shape;
    rank = g_shape.size();
    volume = shape[rank - 1];
    stride[rank - 1] = 1;
    for (int i = stride.size() - 2; i >= 0; i--)
    {
        stride[i] = stride[i + 1] * shape[i + 1];
        volume *= shape[i];
    }
    nn_batch_size = 2 * rank;
    nearest_neighbors.resize(nn_batch_size * volume);
    int temp;
    for (int index = 0; index < volume; index++)
    {
        for (int j = 0; j < rank; j++)
        {
            temp = stride[j] * shape[j];
            nearest_neighbors[2 * j + index * nn_batch_size] = mod((index + stride[j]), temp) + index / temp * temp;
            nearest_neighbors[2 * j + 1 + index * nn_batch_size] = mod((index - stride[j]), temp) + index / temp * temp;
        }
    }
}

void Grid::print()
{
    int extra_print; // used to decide if and how many "]", "[", "\n" and " ," have to be printed
    for (unsigned flat_index = 0; flat_index < volume - 1; flat_index++)
    {
        // special rule for the first element
        if (flat_index == 0)
        {
            extra_print = 0;
        }
        else
        {
            extra_print = -1;
        }

        // count extra prints
        for (int i = 0; i < rank; i++)
        {
            if (flat_index % stride[i] == 0)
            {
                extra_print++;
            }
        }

        // print according to extraprint counting results
        if (extra_print == 0)
        {
            printf(" ,");
        }
        else
        {
            if (flat_index != 0)
            {
                for (int i = 0; i < extra_print; i++)
                {
                    printf("]");
                }
                printf(",");
                if (flat_index != volume - 1)
                {
                    for (int i = 0; i < extra_print; i++)
                    {
                        printf("\n");
                    }
                }
            }
            for (int i = 0; i < rank; i++)
            {
                if (i < rank - extra_print)
                {
                    printf(" ");
                }
                else
                {
                    printf("[");
                }
            }
        }
        printf("%4d", flat_index * stride[rank - 1] % (volume - 1));
    }

    // special print for the last element
    printf(", %4d", volume - 1);
    for (int i = 0; i < rank; i++)
    {
        printf("]");
    }
    printf("\n");
}
