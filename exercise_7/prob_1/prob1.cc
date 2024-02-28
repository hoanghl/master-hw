#include <fstream>
#include <iostream>
#include <mpi.h>
#include <random>
#include <string>

using namespace std;

constexpr int NUM_PACKETS = 100;
constexpr int NUM_PROCS = 5;
const string FILENAME = "prob1_ids.txt";

void log_send(int fromRoutine, int dest_id, int msg)
{
    cout << "= Routine " << fromRoutine << ": Sent to " << dest_id << ": " << msg << endl;
}

void log_recv(int fromRoutine, MPI_Status status, int msg)
{
    cout << "= Routine " << fromRoutine << ": Received from " << status.MPI_SOURCE << ": " << msg << endl;
}

void log_wait(int fromRoutine)
{
    cout << "= Routine " << fromRoutine << ": Waiting" << endl;
    // printf("= Routine %d: Waiting\n", fromRoutine);
}

void log_terminate(int fromRoutine)
{
    cout << "= Routine " << fromRoutine << ": Terminating" << endl;
    // printf("= Routine %d: Terminating\n", fromRoutine);
}

int main(int argc, char *argv[])
{
    // Initialize MP
    const int tag = 0;
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

    // Start
    const int LEN_MSG = 1;

    // 1. Send
    if (id > 0)
    {
        dest_id = 0;

        int msgSend;
        for (int i = 0; i < NUM_PACKETS / NUM_PROCS; ++i)
        {
            msgSend = rand();
            rc = MPI_Send(&msgSend, LEN_MSG, MPI_INT, dest_id, tag, MPI_COMM_WORLD);

            log_send(id, dest_id, msgSend);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // 2. Wait & 3. Receive
    if (id == 0)
    {
        log_wait(id);

        int count = 0;
        int msgRecv;
        int ids[NUM_PACKETS];

        while (count < NUM_PACKETS)
        {
            rc = MPI_Recv(&msgRecv, LEN_MSG, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            log_recv(id, status, msgRecv);

            ids[count] = status.MPI_SOURCE;
            count += 1;

            cout << "count: " << count << endl;
        }

                // Write list of ids of processes sending to process 0
        fstream file(FILENAME, ios::out);
        if (file.is_open())
        {
            for (int i = 0; i < NUM_PACKETS; ++i)
                file << ids[i] << endl;

            file.close();
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // 4. Terminate
    log_terminate(id);
}
