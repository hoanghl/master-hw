#include <iostream>
#include <mpi.h>
#include <random>

using namespace std;

int *genMessage(int size)
{
    int *message = new int[size];
    for (int i = 0; i < size; ++i)
        message[i] = rand();

    return message;
}

int main(int argc, char *argv[])
{
    // Initialize MP
    const int tag = 50;
    int id, ntasks, source_id, dest_id, rc, i;
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

    // Start
    int msgSizes[] = {32, 654, 5324, 6343, 232, 52423, 52492, 12391, 4532, 1237194, 4242342, 23423425, 342342357, 3264824, 312, 54423};
    constexpr int N = 16;

    if (id == 0)
    {
        dest_id = id + 1;
        double start, stop, wallclock;

        for (int i = 0; i < N; ++i)
        {
            int *msgSend = genMessage(msgSizes[i]);
            int *msgRecv = new int[msgSizes[i]];

            start = MPI_Wtime();
            rc = MPI_Send(msgSend, msgSizes[i], MPI_INT, dest_id, tag, MPI_COMM_WORLD);
            rc = MPI_Recv(msgRecv, msgSizes[i], MPI_INT, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
            stop = MPI_Wtime();
            wallclock = stop - start;

            cout << "= Message size: " << msgSizes[i] << " -> time: " << wallclock << endl;

            delete[] msgSend;
            delete[] msgRecv;
        }
    }
    else if (id == 1)
    {
        for (int i = 0; i < N; ++i)
        {
            int *msgRecv = new int[msgSizes[i]];

            rc = MPI_Recv(msgRecv, msgSizes[i], MPI_INT, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
            rc = MPI_Send(msgRecv, msgSizes[i], MPI_INT, status.MPI_SOURCE, tag, MPI_COMM_WORLD);

            delete[] msgRecv;
        }
    }

    return 0;
}
