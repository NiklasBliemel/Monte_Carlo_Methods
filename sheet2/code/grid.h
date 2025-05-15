#ifndef Grid_H
#define Grid_H
#include <vector>
using namespace std;

class Grid
{
public:
    vector<int> shape;
    vector<int> stride;
    vector<int> nearest_neighbors;
    int volume;
    int rank;
    int nn_batch_size;
    Grid(vector<int> &shape);
    void print();
};

#endif