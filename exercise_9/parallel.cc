//
// 1D molecular dynamics simulation code for course
// MATR326 Tools for high performance computing
// Antti Kuronen, University of Helsinki
//
// Try e.g.
// ./a.out 10000 0.001 100000 1 100 0
//

#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <string>
#include <time.h>
#include <vector>

#define d 1.0
#define k1 1.0
#define k2 0.1
#define xsc 2.35

using namespace std;

// --------------------------------------------------------------------------------------

int main(int argc, char **argv)
{

    vector<double> x;  // atom positions
    vector<double> v;  //      velocities
    vector<double> v0; //      previous veloocities (leap frog needs them)
    vector<double> a;  //      accelerations

    double epsum, eksum; // system energies
    double dt;           // time step
    double vsc;          // mean initial velocity
    double box;          // system size
    int nat;             // number of atoms
    int maxt;            // number of time steps simulated
    int eout;            // energy output interval
    int coout;           // coordinate output interval (lot of data, beware!)

    int n;
    double vsum, rn;

    // Get number of atoms, time step and simulation length from command line
    if (argc < 5 || argc > 7)
    {
        cerr << "usage: " << argv[0] << " nat dt maxt vsc [eout [coout]]\n";
        cerr << "    nat   = number of atoms\n";
        cerr << "    dt    = time step\n";
        cerr << "    maxt  = number of time steps in simulation\n";
        cerr << "    vsc   = mean velocity of atoms in the beginning ('temperature')\n";
        cerr << "    eout  = interval for printing energies to stdout\n";
        return (128);
    }

    coout = 0;
    eout = 1;
    nat = atoi(*++argv);
    dt = atof(*++argv);
    maxt = atoi(*++argv);
    vsc = atof(*++argv);

    if (argc > 5)
        eout = atoi(*++argv);
    if (argc > 6)
        coout = atoi(*++argv);

    x = vector<double>(nat);
    v = vector<double>(nat);
    v0 = vector<double>(nat);
    a = vector<double>(nat);

    // Initialize atoms positions and give them random velocities
    box = nat;
    srand(time(NULL));
    for (int i = 0; i < nat; i++)
    {
        x[i] = i;
        rn = (double)rand() / RAND_MAX;
        v[i] = vsc * (rn - 0.5); // Scale the velocities to vsc*[-½,½]
    }

    n = 0;

    // Simulation proper

    int nThreads = 4;
    if (nat % nThreads != 0)
    {
        cerr << "Num. atoms is not divisible by num. threads.";
        exit(1);
    }
    int nPerThread = nat / nThreads;

    omp_set_num_threads(nThreads);

    vector<double> ep = vector<double>(nat);
    vector<double> ek = vector<double>(nat);

    auto t0 = std::chrono::system_clock::now();

    for (n = 0; n < maxt; n++)
    {
        for (int i = 0; i < nat; i++)
            v0[i] = v[i];
        //  shared (a)  reduction(+ : ep)
        // #pragma omp parallel for
        for (int i = 0; i < nat; i++)
        {
            // New potential energy and acceleration
            int j, k;
            double dxl, dxr;

            j = i - 1;
            if (j < 0)
                j = nat - 1;
            k = i + 1;
            if (k >= nat)
                k = 0;

            dxl = x[i] - x[j];
            dxr = x[k] - x[i];
            if (dxl < -box / 2.0)
                dxl += box;
            if (dxl >= box / 2.0)
                dxl -= box;
            if (dxr < -box / 2.0)
                dxr += box;
            if (dxr >= box / 2.0)
                dxr -= box;
            dxl -= d;
            dxr -= d;

            ep[i] = (k1 * (dxl * dxl + dxr * dxr) + k2 * (dxl * dxl * dxl + dxr * dxr * dxr)) / 2.0;
            a[i] = -(2.0 * k1 * (dxl - dxr) + 3.0 * k2 * (dxl * dxl - dxr * dxr));
        }

        // double ek = 0;

#pragma omp parallel for
        for (int iThread = 0; iThread < nThreads; ++iThread)
            for (int i = iThread * nPerThread; i < (iThread + 1) * nPerThread; ++i)
            {
                double vave;
                // Leap frog integration algorithm: update position and velocity
                v[i] = v[i] + dt * a[i];
                x[i] = x[i] + dt * v[i];
                // Check periodic boundary conditions
                if (x[i] < 0.0)
                    x[i] = x[i] + box;
                if (x[i] >= box)
                    x[i] = x[i] - box;
                // Calculate kinetic energy (note: mass=1)
                vave = (v0[i] + v[i]) / 2.0;
                ek[i] = 1.0 / 2.0 * vave * vave;
                // ek = 1.0 / 2.0 * vave * vave;
            }

        // Calculate and print total potential end kinetic energies
        // and their sum (= total energy that should be conserved).
        if (n % eout == 0)
        {
            epsum = accumulate(ep.begin(), ep.end(), 0.0);
            eksum = accumulate(ek.begin(), ek.end(), 0.0);
            // eksum = ek;
            printf("%20.10g %20.10g %20.10g %20.10g\n", dt * n, epsum + eksum, epsum, eksum);
        }
    }

    auto t1 = std::chrono::system_clock::now();
    auto wct = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    cerr << "Wall clock time: " << wct.count() / 1000.0 << " seconds\n";

    return (0);
}
