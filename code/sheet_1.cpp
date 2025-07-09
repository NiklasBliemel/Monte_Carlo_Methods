#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <string>

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

void normal_N(vector<double> &vec, int N, mt19937 &gen)
{
    // define normal distribution name(double mean, double std);
    normal_distribution<double> normaldist(0, 1);

    vec.resize(N);
    for (size_t i = 0; i < N; i++)
    {
        vec[i] = f(normaldist(gen));
    }
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

    int N_size = 101;
    vector<int> Ns(N_size);
    for (int i = 0; i < Ns.size(); i++)
    {
        Ns[i] = (int)pow(10, 5.0 / 100.0 * i + 1);
    }

    int M = 100;
    double mean;
    string filename;
    string Ns_filename;
    vector<double> sample;
    filename = "/Users/niklasbliemel/Study/Monte_Carlo_Methods/plot_lab/sheet_1_data/data.bin";
    Ns_filename = "/Users/niklasbliemel/Study/Monte_Carlo_Methods/plot_lab/sheet_1_data/data_Ns.bin";
    ofstream wf(filename, ios::binary);
    ofstream wfNs(Ns_filename, ios::binary);
    if (!wf && !wfNs)
    {
        cout << "Cannot open files!" << endl;
        return 1;
    }
    wf.write(reinterpret_cast<char *>(&N_size), sizeof(int));
    wf.write(reinterpret_cast<char *>(&M), sizeof(int));
    for (int i = 0; i < Ns.size(); i++)
    {
        wfNs.write(reinterpret_cast<char *>(&Ns[i]), sizeof(int));
        cout << "run " << i + 1 << "/" << Ns.size() << endl;
        for (int j = 0; j < M; j++)
        {
            normal_N(sample, Ns[i], gen);
            mean = vec_mean(sample);
            wf.write(reinterpret_cast<char *>(&mean), sizeof(double));
        }
    }
    wf.close();
    wfNs.close();
}
