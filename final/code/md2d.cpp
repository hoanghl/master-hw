#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <time.h>
#include <vector>

#define d 1.0
#define k1 1.0
#define k2 1.0
#define k3 0.4
#define xsc 2.35

const int TAG_POPULATE = 100;

typedef struct
{
    double x, y;
} Coordinate;
typedef struct
{
    int top, left, right, bottom;
} Neightbor;

using namespace std;

ofstream fcoord;
string coord_file = "coords.dat";

// --------------------------------------------------------------------------------------
double getRand()
{
    return (double)rand() / RAND_MAX;
}

void getNeighbor(int parID, int direction, int nuc, int &parIDNb)
{
    int i = parID / nuc, j = parID % nuc;
    int iNb = i, jNb = j;

    switch (direction)
    {
    case 1:
        iNb--;
        if (iNb < 0)
            iNb = nuc - 1;
        break;
    case 3:
        iNb++;
        if (iNb >= nuc)
            iNb = 0;
        break;
    case 0:
        jNb--;
        if (jNb < 0)
            jNb = nuc - 1;
        break;
    case 2:
        jNb++;
        if (jNb >= nuc)
            jNb = 0;
        break;
    }

    parIDNb = iNb * nuc + jNb;
}

void accel(
    int nuc,
    int i,
    vector<double> &u,
    vector<double> &a_x,
    vector<double> &a_y,
    double box,
    vector<double> &x,
    vector<double> &y)
{

    double u2 = 0, u3 = 0;
    a_x[i] = a_y[i] = 0;

    double dx, dy, r, fx2, fy2, fx3, fy3;

    for (int j = 0; j < 4; ++j)
    {
        int in;
        getNeighbor(i, j, nuc, in);
        dx = x[in] - x[i];
        if (dx < -box / 2.0)
            dx += box;
        if (dx >= box / 2.0)
            dx -= box;

        dy = y[in] - y[i];
        if (dy < -box / 2.0)
            dy += box;
        if (dy >= box / 2.0)
            dy -= box;

        r = sqrt(dx * dx + dy * dy);

        u2 = u2 + 1.0 / 2.0 * k2 * pow(r - d, 2);
        u3 = u3 + 1.0 / 3.0 * k3 * pow(r - d, 3);
        fx2 = k2 * (r - d) * dx / r;
        fx3 = k3 * pow(r - d, 2) * dx / r;
        fy2 = k2 * (r - d) * dy / r;
        fy3 = k3 * pow(r - d, 2) * dy / r;
        a_x[i] += fx2 + fx3;
        a_y[i] += fy2 + fy3;
    }

    // Remember the factor 1 / 2 due to 'double counting
    u[i] = (u2 + u3) / 2.0;
}

// --------------------------------------------------------------------------------------

int main(int argc, char **argv)
{

    vector<double> x;   //      atom positions
    vector<double> y;   //      atom positions
    vector<double> v_x; //      velocities
    vector<double> v_y; //      velocities
    vector<double> a_x; //      accelerations
    vector<double> a_y; //      accelerations
    vector<double> ep;  //      potential energies
    vector<double> ek;  //      kinetic energies

    double epsum, eksum; // system energies
    double dt;           // time step
    double vsc;          // mean initial velocity
    double box;          // system size
    int nuc;             // number of unit cells
    int nat;             // number of atoms
    int maxt;            // number of time steps simulated
    int eout;            // energy output interval
    int coout;           // coordinate output interval (lot of data, beware!)

    int i;
    double vsum, vx, vy;

    // Get number of atoms, time step and simulation length from command line
    if (argc < 5 || argc > 7)
    {
        cerr << "usage: " << argv[0] << " nuc dt maxt vsc [eout [coout]]\n";
        cerr << "    nuc   = number of unit cells\n";
        cerr << "    dt    = time step\n";
        cerr << "    maxt  = number of time steps in simulation\n";
        cerr << "    vsc   = mean velocity of atoms in the beginning ('temperature')\n";
        cerr << "    eout  = interval for printing energies to stdout\n";
        cerr << "    coout = interval for printing coordinates to '" << coord_file << "'\n";
        return (128);
    }

    coout = 0;
    eout = 1;
    nuc = atoi(*++argv);
    dt = atof(*++argv);
    maxt = atoi(*++argv);
    vsc = atof(*++argv);

    if (argc > 5)
        eout = atoi(*++argv);
    if (argc > 6)
        coout = atoi(*++argv);

    nat = nuc * nuc;
    x = vector<double>(nat);
    y = vector<double>(nat);
    v_x = vector<double>(nat);
    v_y = vector<double>(nat);
    a_x = vector<double>(nat);
    a_y = vector<double>(nat);

    ep = vector<double>(nat);
    ek = vector<double>(nat);

    // Initialize atoms positions and give them random velocities
    box = nuc;
    srand(time(NULL));
    for (i = 0; i < nat; i++)
    {
        x[i] = (i / nuc) * 1.0;
        y[i] = (i % nuc) * 1.0;

        // Scale the velocities to vsc*[-½,½]
        v_x[i] = vsc * (getRand() - 0.5);
        v_y[i] = vsc * (getRand() - 0.5);
    }

    // Remove center of mass velocity
    double sum_vx = 0, sum_vy = 0;
    for (int i = 0; i < nat; ++i)
    {
        sum_vx += v_x[i];
        sum_vy += v_y[i];
    }
    sum_vx /= nat;
    sum_vy /= nat;
    for (int i = 0; i < nat; ++i)
    {
        v_x[i] -= sum_vx;
        v_y[i] -= sum_vy;
    }

    // Simulation proper
    auto t0 = std::chrono::system_clock::now();

    double v0_x, v0_y;
    for (int n = 0; n < maxt; n++)
    {
        for (i = 0; i < nat; i++)
        {
            accel(nuc, i, ep, a_x, a_y, box, x, y);
            v0_x = v_x[i];
            v0_y = v_y[i];

            // Leap frog integration algorithm: update position and velocity
            v_x[i] += dt * a_x[i];
            v_y[i] += dt * a_y[i];
            x[i] += dt * v_x[i];
            y[i] += dt * v_y[i];

            // Check periodic boundary conditions
            if (x[i] < 0.0)
                x[i] += box;
            else if (x[i] >= box)
                x[i] -= box;
            if (y[i] < 0.0)
                y[i] += box;
            else if (y[i] >= box)
                y[i] -= box;

            vx = (v0_x + v_x[i]) / 2.0;
            vy = (v0_y + v_y[i]) / 2.0;
            ek[i] = 1.0 / 2.0 * (pow(vx, 2) + pow(vy, 2));
        }

        // Calculate and print total potential end kinetic energies
        // and their sum (= total energy that should be conserved).

        epsum = accumulate(ep.begin(), ep.end(), 0.0);
        eksum = accumulate(ek.begin(), ek.end(), 0.0);

        if (eout > 0 && n % eout == 0)
            printf("%20.10g %20.10g %20.10g %20.10g\n", dt * n, epsum + eksum, epsum, eksum);
    }

    auto t1 = std::chrono::system_clock::now();
    auto wct = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    double wtime = wct.count() / 1000.0;
    if (eout > 0)
        printf("Wall clock time: %.4f seconds\n", wtime);
    else
        printf("%f\n", wtime);

    if (coout > 0)
        fcoord.close();
    return (0);
}
