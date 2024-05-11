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

void getNeighList(vector<vector<int>> &neigh, int nat, int nuc)
{
    int row, col;
    int left_row, left_col;
    int right_row, right_col;
    int top_row, top_col;
    int bottom_row, bottom_col;

    for (int i = 0; i < nat; ++i)
    {
        row = i / nuc;
        col = i % nuc;

        left_row = right_row = row;
        top_col = bottom_col = col;

        right_col = col + 1;
        if (right_col >= nuc)
            right_col = 0;
        left_col = col - 1;
        if (left_col < 0)
            left_col = nuc - 1;
        top_row = row - 1;
        if (top_row < 0)
            top_row = nuc - 1;
        bottom_row = row + 1;
        if (bottom_row >= nuc)
            bottom_row = 0;

        vector<int> v1 = {
            left_col + (left_row)*nuc,
            top_col + (top_row)*nuc,
            right_col + (right_row)*nuc,
            bottom_col + (bottom_row)*nuc,
        };
        neigh.push_back(v1);
    }
}

void accel(int nat, int i, double *u, Coordinate *a, double box, vector<Coordinate> *x, vector<vector<int>> *neighbors)
{

    double u2 = 0, u3 = 0;
    a->x = 0;
    a->y = 0;
    double dx, dy, r, fx2, fy2, fx3, fy3;

    for (int j = 0; j < 4; ++j)
    {
        int in = (*neighbors)[i][j];
        dx = (*x)[in].x - (*x)[i].x;
        if (dx < -box / 2.0)
            dx += box;
        if (dx >= box / 2.0)
            dx -= box;

        dy = (*x)[in].y - (*x)[i].y;
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
        a->x += fx2 + fx3;
        a->y += fy2 + fy3;
    }

    // Remember the factor 1 / 2 due to 'double counting
    *u = (u2 + u3) / 2.0;
}

// --------------------------------------------------------------------------------------

int main(int argc, char **argv)
{

    vector<Coordinate> x;  //      atom positions
    vector<Coordinate> v;  //      velocities
    vector<Coordinate> a;  //      accelerations
    vector<Coordinate> v0; //      previous veloocities (leap frog needs them)
    vector<double> ep;     //      potential energies
    vector<double> ek;     //      kinetic energies

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
    x = vector<Coordinate>(nat);
    v = vector<Coordinate>(nat);
    a = vector<Coordinate>(nat);
    v0 = vector<Coordinate>(nat);

    ep = vector<double>(nat);
    ek = vector<double>(nat);

    // Initialize atoms positions and give them random velocities
    box = nuc;
    srand(time(NULL));
    for (i = 0; i < nat; i++)
    {
        x[i] = {
            .x = (i / nuc) * 1.0,
            .y = (i % nuc) * 1.0,
        };
        v[0] = {.x = 0, .y = 0};

        // Scale the velocities to vsc*[-½,½]
        v[i] = {
            .x = vsc * (getRand() - 0.5),
            .y = vsc * (getRand() - 0.5),
        };
    }

    // Remove center of mass velocity
    double sum_vx = 0, sum_vy = 0;
    for (int i = 0; i < nat; ++i)
    {
        sum_vx += v[i].x;
        sum_vy += v[i].y;
    }
    sum_vx /= nat;
    sum_vy /= nat;
    for (int i = 0; i < nat; ++i)
    {
        v[i].x -= sum_vx;
        v[i].y -= sum_vy;
    }

    // Find neighbors for each atam
    vector<vector<int>> neighbors;
    getNeighList(neighbors, nat, nuc);

    // Simulation proper
    auto t0 = std::chrono::system_clock::now();

    double v0_x, v0_y;
    for (int n = 0; n < maxt; n++)
    {
        for (i = 0; i < nat; i++)
        {
            accel(nat, i, &ep[i], &a[i], box, &x, &neighbors);
            v0_x = v[i].x;
            v0_y = v[i].y;

            // Leap frog integration algorithm: update position and velocity
            v[i].x = v[i].x + dt * a[i].x;
            v[i].y = v[i].y + dt * a[i].y;
            x[i].x = x[i].x + dt * v[i].x;
            x[i].y = x[i].y + dt * v[i].y;

            // Check periodic boundary conditions
            if (x[i].x < 0.0)
                x[i].x += box;
            if (x[i].y < 0.0)
                x[i].y += box;
            if (x[i].x >= box)
                x[i].x -= box;
            if (x[i].y >= box)
                x[i].y -= box;

            vx = (v0_x + v[i].x) / 2.0;
            vy = (v0_y + v[i].y) / 2.0;
            ek[i] = 1.0 / 2.0 * (pow(vx, 2) + pow(vy, 2));
        }

        // Calculate and print total potential end kinetic energies
        // and their sum (= total energy that should be conserved).

        epsum = accumulate(ep.begin(), ep.end(), 0.0);
        eksum = accumulate(ek.begin(), ek.end(), 0.0);
        if (eout > 0)
            if (n % eout == 0)
                printf("%20.10g %20.10g %20.10g %20.10g\n", dt * n, epsum + eksum, epsum, eksum);
    }

    auto t1 = std::chrono::system_clock::now();
    auto wct = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    cerr << "Wall clock time: " << wct.count() / 1000.0 << " seconds\n";

    if (coout > 0)
        fcoord.close();
    return (0);
}
