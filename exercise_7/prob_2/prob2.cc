#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <random>
#include <string>
#include <time.h>

using namespace std;

const int TAG = 0;
const int LEN_MSG = 1;
const unsigned int N = 1000000000;

double getRand()
{
    return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}

void log_send(int fromRoutine, int dest_id, int msg)
{
    cout << "= Routine " << fromRoutine << ": Sent to " << dest_id << ": " << msg << endl;
}

void log_recv(int fromRoutine, MPI_Status status, int msg)
{
    cout << "= Routine " << fromRoutine << ": Received from " << status.MPI_SOURCE << ": " << msg << endl;
}

int main(int argc, char *argv[])
{
    // Initialize MP
    int id, ntasks, dest_id, rc;
    MPI_Status status;
    int pnlen;
    char pname[MPI_MAX_PROCESSOR_NAME];

    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS)
    {
        printf("MPI initialization failed\n");
        exit(1);
    }
    rc = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &id);
    rc = MPI_Get_processor_name(pname, &pnlen);

    if (id > 0)
    {
        // 1. Receive seed from process 0
        int seed;
        MPI_Recv(&seed, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        log_recv(id, status, seed);

        srand(seed);

        // 2. Start looping
        int countInside = 0;
        double x, y;
        for (int i = 0; i < N; ++i)
        {
            x = getRand();
            y = getRand();

            if (x * x + y * y < 1)
                ++countInside;
        }

        // 3. Send back result to process 0
        dest_id = 0;
        MPI_Send(&countInside, 1, MPI_INT, dest_id, TAG, MPI_COMM_WORLD);
        log_send(id, dest_id, countInside);

        MPI_Barrier(MPI_COMM_WORLD);
    }
    else
    {
        // 1. Send seed to other processes
        srand(time(0));

        int seed;
        for (int dest_id = 1; dest_id < ntasks; ++dest_id)
        {
            seed = rand();
            MPI_Send(&seed, 1, MPI_INT, dest_id, TAG, MPI_COMM_WORLD);
            log_send(id, dest_id, seed);
        }

        // 2. Receive resuls from other processes

        int countProc = 0;
        double countPoint;
        int msgRecv;

        while (countProc < ntasks - 1)
        {
            MPI_Recv(&msgRecv, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            log_recv(id, status, msgRecv);

            countPoint += static_cast<double>(msgRecv);
            countProc += 1;
        }

        // 3. Calculate the pi
        double pi = 4 * countPoint / ((ntasks - 1) * static_cast<double>(N));
        cout << "pi = " << setprecision(10) << pi << endl;

        MPI_Barrier(MPI_COMM_WORLD);
    }
}
