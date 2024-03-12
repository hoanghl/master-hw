#include <chrono>
#include <iostream>
#include <numeric>
#include <time.h>
#include <vector>

#include <mpi.h>

using namespace std;

const double D = 1.0;
const double K1 = 1.0;
const double K2 = 0.1;
const double XSC = 2.35;
const int TAG_POPULATE = 100;

enum Energy
{
    potential,
    kinetic
};

void accel(int nAtoms, int i, double *u, double *a, double box, vector<double> &x, double xLeftMost, double xRightMost)
{
    // Calculate the potential energy u and acceleration a of atom i
    int j, k;
    double dxl, dxr;

    // NOTE: HoangLe [Mar-12]: This is the part where a particle interacts with the neighbors

    double xLeft = i - 1 < 0 ? xLeftMost : x[i - 1];
    double xRight = i + 1 >= nAtoms ? xRightMost : x[i + 1];

    dxl = x[i] - xLeft;
    dxr = xRight - x[i];
    if (dxl < -box / 2.0)
        dxl += box;
    if (dxl >= box / 2.0)
        dxl -= box;
    if (dxr < -box / 2.0)
        dxr += box;
    if (dxr >= box / 2.0)
        dxr -= box;
    dxl -= D;
    dxr -= D;

    *u = (K1 * (dxl * dxl + dxr * dxr) + K2 * (dxl * dxl * dxl + dxr * dxr * dxr)) / 2.0;
    *a = -(2.0 * K1 * (dxl - dxr) + 3.0 * K2 * (dxl * dxl - dxr * dxr));
}

int getTag(int rank, Energy energy)
{
    int tag = 0;

    switch (energy)
    {
    case potential:
        tag = rank * 2;
        break;
    case kinetic:
        tag = rank * 2 + 1;
        break;

    default:
        cout << "Invalid energy type: " << energy << endl;
        exit(1);
    }

    return tag;
}

double getRand()
{
    return static_cast<double>(rand()) / RAND_MAX;
}

int main(int argc, char *argv[])
{
    // Parse the arguments
    if (argc != 5)
    {
        cerr << "usage: " << argv[0] << " nAtoms dt maxt vsc \n";
        cerr << "    nAtoms   = number of atoms\n";
        cerr << "    dt    = time step\n";
        cerr << "    maxt  = number of time steps in simulation\n";
        cerr << "    vsc   = mean velocity of atoms in the beginning ('temperature')\n";
        return (128);
    }

    int nAtoms;     // number of atoms
    int nAtomsEach; // number of atoms each process
    int maxt;       // number of time steps simulated
    double dt;      // time step
    double vsc;     // mean initial velocity

    nAtoms = atoi(*++argv);
    dt = atof(*++argv);
    maxt = atoi(*++argv);
    vsc = atof(*++argv);

    // Init things:
    // 1. Init MPI_comm
    int rank, nProcs;
    MPI_Comm comm;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

    int dim[] = {nProcs}, periods[] = {1}, reorder = false;
    MPI_Cart_create(MPI_COMM_WORLD, 1, dim, periods, reorder, &comm);

    // assert num. particles is divisible by the num. processes
    if (nAtoms % nProcs != 0)
    {
        cerr << "Error: Number of atoms not divisible by number of processes: nAtoms = "
             << nAtoms << " ; nProcs = " << nProcs << endl;
        return 1;
    }
    nAtomsEach = nAtoms / nProcs;

    // 2. Init atom position in each processes
    vector<double> x = vector<double>(nAtomsEach);
    vector<double> v = vector<double>(nAtomsEach);
    vector<double> v0 = vector<double>(nAtomsEach);
    vector<double> a = vector<double>(nAtomsEach);
    vector<double> ep = vector<double>(nAtomsEach);
    vector<double> ek = vector<double>(nAtomsEach);

    double box = nAtoms;
    srand(time(NULL));
    for (int i = 0; i < nAtomsEach; i++)
    {
        x[i] = rank * nAtomsEach + i;
        v[i] = vsc * (getRand() - 0.5); // Scale the velocities to vsc*[-½,½]
    }

    // Determine the rank of the left and right neighbor
    int direction = 0;
    int disp = 1;
    int procLeft = 0, procRight = 0;
    MPI_Cart_shift(comm, direction, disp, &procLeft, &procRight);

    MPI_Request requestSent;
    vector<MPI_Request> vReqRecvKinetic(nProcs);
    vector<MPI_Request> vReqRecvPotential(nProcs);
    vector<MPI_Request> vReqSend(2);
    vector<double> ekProcs(nProcs);
    vector<double> epProcs(nProcs);

    auto t0 = std::chrono::system_clock::now();
    // Start simulation
    for (int step = 0; step < maxt; ++step)
    {

        // cout << "rank = " << rank << ": start step: " << step << endl;

        // Each process calculates the following:
        // Get initial velocity
        for (int i = 0; i < nAtomsEach; i++)
            v0[i] = v[i];

        // cout << "rank = " << rank << ": procLeft: " << procLeft << " -- procRight = " << procRight << endl;

        // Populate the leftmost particle's coordinate to the left process and receive back the coordinate
        double xLeftMost = 0;
        double xRightMost = 0;

        MPI_Sendrecv(&x[0], 1, MPI_DOUBLE, procLeft, TAG_POPULATE,
                     &xRightMost, 1, MPI_DOUBLE, procRight, TAG_POPULATE, comm,
                     &status);

        // cout << "rank = " << rank << ": finish leftmost" << endl;

        MPI_Sendrecv(&x[x.size() - 1], 1, MPI_DOUBLE, procRight, TAG_POPULATE,
                     &xLeftMost, 1, MPI_DOUBLE, procLeft, TAG_POPULATE, comm,
                     &status);

        // cout << "rank = " << rank << ": finish rightmost" << endl;

        // Calculate new acceleration and ep for each particle
        for (int i = 0; i < nAtomsEach; ++i)
        {
            accel(nAtomsEach, i, &ep[i], &a[i], box, x, xLeftMost, xRightMost);
        }

        // Calculate the kinetic energy
        for (int i = 0; i < nAtomsEach; ++i)
        {
            // Leap frog integration algorithm: update position and velocity
            v[i] = v[i] + dt * a[i];
            x[i] = x[i] + dt * v[i];
            // Check periodic boundary conditions
            if (x[i] < 0.0)
                x[i] = x[i] + box;
            if (x[i] >= box)
                x[i] = x[i] - box;
            // Calculate kinetic energy (note: mass=1)
            double vave = (v0[i] + v[i]) / 2.0;
            ek[i] = 1.0 / 2.0 * vave * vave;
        }

        // Send all accumulated ep and ek to rank 0
        double epSum = accumulate(ep.begin(), ep.end(), 0.0);
        double ekSum = accumulate(ek.begin(), ek.end(), 0.0);

        if (rank != 0)
        {
            MPI_Isend(&epSum, 1, MPI_DOUBLE, 0, getTag(rank, potential), MPI_COMM_WORLD, &vReqSend[0]);
            MPI_Isend(&ekSum, 1, MPI_DOUBLE, 0, getTag(rank, kinetic), MPI_COMM_WORLD, &vReqSend[1]);
            MPI_Wait(&vReqSend[0], MPI_STATUS_IGNORE);
            MPI_Wait(&vReqSend[1], MPI_STATUS_IGNORE);
        }
        else
        {
            ekProcs[0] = ekSum;
            epProcs[0] = epSum;

            for (int i = 1; i < nProcs; ++i)
            {
                MPI_Irecv(&ekProcs[i], 1, MPI_DOUBLE, i, getTag(i, kinetic), MPI_COMM_WORLD, &vReqRecvKinetic[i]);
                MPI_Irecv(&epProcs[i], 1, MPI_DOUBLE, i, getTag(i, potential), MPI_COMM_WORLD, &vReqRecvPotential[i]);
            }

            for (int i = 1; i < nProcs; ++i)
            {
                MPI_Wait(&vReqRecvKinetic[i], MPI_STATUS_IGNORE);
                MPI_Wait(&vReqRecvPotential[i], MPI_STATUS_IGNORE);
            }

            double eksumFinal = accumulate(ekProcs.begin(), ekProcs.end(), 0.0);
            double epsumFinal = accumulate(epProcs.begin(), epProcs.end(), 0.0);

            printf("%20.10g %20.10g %20.10g %20.10g\n", dt * step, epsumFinal + eksumFinal, epsumFinal, eksumFinal);
        }
    }

    if (rank == 0)
    {
        auto t1 = std::chrono::system_clock::now();
        auto wct = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
        cout << "Wall clock time: " << wct.count() / 1000.0 << " seconds\n";
    }

    MPI_Finalize();

    return 0;
}
