#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

void log_send(int msg, int fromRoutine)
{
    printf("= Routine %d: Sent '%d'\n", fromRoutine, msg);
}

void log_recv(int msg, int fromRoutine, MPI_Status status)
{
    printf("= Routine %d: Received from %d: '%d'\n", fromRoutine, status.MPI_SOURCE, msg);
}

int main(int argc, char *argv[])
{
    const int tag = 50;
    int id, ntasks, source_id, dest_id, rc, i;
    MPI_Status status;
    int msgSend, msgRecv, pnlen;
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

    if (id == 0)
    {
        msgSend = 1000;
        dest_id = 1;
    }
    else
    {
        msgSend = 100;
        dest_id = 0;
    }

    rc = MPI_Send(&msgSend, 1, MPI_INT, dest_id, tag, MPI_COMM_WORLD);
    log_send(msgSend, id);

    rc = MPI_Recv(&msgRecv, 1, MPI_INT, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
    log_recv(msgRecv, id, status);

    msgSend = msgRecv + 1;
    rc = MPI_Send(&msgSend, 1, MPI_INT, dest_id, tag, MPI_COMM_WORLD);
    log_send(msgSend, id);

    rc = MPI_Recv(&msgRecv, 1, MPI_INT, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
    log_recv(msgRecv, id, status);

    rc = MPI_Finalize();
    exit(0);
}
