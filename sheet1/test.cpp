#include <vector>
#include <iostream>
#include <chrono>
#include <random>
using namespace std;
#define sqr(x) ((x) * (x))

double vec_mean(vector<double> &vec)
{
    double out = 0;
    for (int i = 0; i < vec.size(); i++)
    {
        out += vec[i];
    }
    return out / vec.size();
}

double vec_var(vector<double> &vec)
{
    double mean_value = vec_mean(vec);
    double out = 0;
    for (int i = 0; i < vec.size(); i++)
    {
        out += sqr(mean_value - vec[i]);
    }
    return out / (vec.size() - 1);
}

double vec_std(vector<double> &vec)
{
    return sqrt(vec_var(vec));
}

double f(double x)
{
    return cos(x);
}

int main()
{
    // make unique seed based on current time
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    cout << "Seed: " << seed << endl;

    // define and initialize rng
    mt19937 gen;
    gen.seed(seed);

    // define uniform distribution for int name(int start, int end);
    uniform_int_distribution<int> unidist(0, 100);

    // define unifrom distribution for double name(double start, double end);
    uniform_real_distribution<double> unireal(0, 1);

    // define normal distribution name(double mean, double std);
    normal_distribution<double> normaldist(0, 1);

    vector<double> vec(100);
    for (int i = 0; i < vec.size(); i++)
    {
        vec[i] = unireal(gen);
    }

    cout << "Mean: " << vec_mean(vec) << endl;
    cout << "Std: " << vec_std(vec) << endl;
}
