#include <iostream>
#include <mpi.h>

using namespace std;

void log_send(int fromRoutine, int *msg)
{
    cout << "= Routine " << fromRoutine << ": Sent [" << msg[0] << ", " << msg[1] << "]" << endl;
    // printf("= Routine %d: Sent [%d, %d]\n", fromRoutine, msg[0], msg[1]);
}

void log_recv(int fromRoutine, MPI_Status status, int *msg)
{
    cout << "= Routine " << fromRoutine << ": Received from " << status.MPI_SOURCE << ": [" << msg[0] << ", " << msg[1] << "]" << endl;
    // printf("= Routine %d: Received from %d: [%d, %d]\n", fromRoutine, status.MPI_SOURCE, msg[0], msg[1]);
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
    const int LEN_MSG = 2;

    // 1. Send
    if (id < ntasks - 1)
    {
        dest_id = id + 1;

        int msgSend[] = {id, (id + 1) * 10};
        msgSend[0] = id;

        rc = MPI_Send(msgSend, LEN_MSG, MPI_INT, dest_id, tag, MPI_COMM_WORLD);
        log_send(id, msgSend);
    }

    // 2. Wait
    if (id > 0)
    {
        log_wait(id);
    }

    // 3. Receive
    if (id > 0)
    {
        int msgRecv[LEN_MSG];
        rc = MPI_Recv(&msgRecv, LEN_MSG, MPI_INT, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
        log_recv(id, status, msgRecv);
    }

    log_terminate(id);

    return 0;
}
