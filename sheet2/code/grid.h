#ifndef grid_H
#define grid_H
#include <vector>
using namespace std;

class Grid
{
public:
    vector<int> data;
    vector<int> shape;
    vector<int> stride;
    Grid(vector<int> &data, vector<int> &shape);
    vector<int> nn(int index);
    void print();
};

#endif