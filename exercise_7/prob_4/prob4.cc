#include <iostream>
#include <queue>
#include <vector>

#include <mpi.h>

using namespace std;

const int ID_LEADER = 0;
const int TAG = 0;

struct ProcInfo
{
    queue<int> waitFrom;
    int sendTo;
};

void log_send(int fromRoutine, int destId, int msg)
{
    cout << "= Routine " << fromRoutine << ": Sent to " << destId << ": " << msg << endl;
}

void log_recv(int fromRoutine, MPI_Status status, int msg)
{
    cout << "= Routine " << fromRoutine << ": Received from " << status.MPI_SOURCE << ": " << msg << endl;
}

int main(int argc, char *argv[])
{
    // Initialize MP
    int id, nTasks, rc;
    MPI_Status status;
    int pnlen;
    char pname[MPI_MAX_PROCESSOR_NAME];

    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS)
    {
        printf("MPI initialization failed\n");
        exit(1);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Get_processor_name(pname, &pnlen);

    vector<ProcInfo> confgs(nTasks);

    // 1. Process 0 run algorithm to find the configuration for each process
    if (id == 0)
    {
        queue<int> processes;

        // 1.1. Init
        for (int i = 0; i < nTasks; ++i)
            processes.push(i);

        // 1.2. Start algorithm
        while (processes.size() > 1)
        {
            int nLeft = processes.size();

            while (nLeft > 0)
            {
                nLeft -= 2;
                int x1 = processes.front();
                processes.pop();
                int x2 = processes.front();
                processes.pop();

                processes.push(x1);
                confgs[x1].waitFrom.push(x2);
                confgs[x2].sendTo = x1;
            }
        }
    }

    // 2. Process 0 populate configurations to
    if (id == 0)
    {
        for (int destID = 1; destID < nTasks; ++destID)
        {
            // 2.1. Send 'sendTo'
            MPI_Send(&confgs[destID].sendTo, 1, MPI_INT, destID, TAG, MPI_COMM_WORLD);
            // log_send(id, destID, confgs[destID].sendTo);

            // 2.2. Send length of 'waitFrom' and 'waitFrom', only processes with ID is even number need 'waitFrom'
            if (destID % 2 == 0)
            {
                int nWaitFrom = confgs[destID].waitFrom.size();
                // 2.2.1. Send length of 'waitFrom'
                MPI_Send(&nWaitFrom, 1, MPI_INT, destID, TAG, MPI_COMM_WORLD);
                // log_send(id, destID, nWaitFrom);

                // 2.2.2. Create an integer array to store values in 'waitFrom'
                int *waitFrom = new int[nWaitFrom];
                for (int i = 0; i < nWaitFrom; ++i)
                {
                    waitFrom[i] = confgs[destID].waitFrom.front();
                    confgs[destID].waitFrom.pop();
                }

                // 2.2.3. Send array 'waitFrom' to sub-process
                MPI_Send(waitFrom, nWaitFrom, MPI_INT, destID, TAG, MPI_COMM_WORLD);
            }
        }
    }

    // 3. Each sub-process receives configuration from process 0
    int sendTo = 0;
    int nWaitFrom = 0;
    int *waitFrom = NULL;
    if (id == 0)
    {
        nWaitFrom = confgs[0].waitFrom.size();
        waitFrom = new int[nWaitFrom];
        for (int i = 0; i < nWaitFrom; ++i)
        {
            waitFrom[i] = confgs[0].waitFrom.front();
            confgs[0].waitFrom.pop();
        }
    }
    else
    {
        MPI_Recv(&sendTo, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        // log_recv(id, status, nWaitFrom);

        if (id % 2 == 0)
        {
            MPI_Recv(&nWaitFrom, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            // log_recv(id, status, nWaitFrom);

            waitFrom = new int[nWaitFrom];
            MPI_Recv(waitFrom, nWaitFrom, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        }
    }

    // 4. Each sub-process calculates and send result to designated process
    int sum = id;
    for (int i = 0; i < nWaitFrom; ++i)
    {
        int waitFromID = waitFrom[i];
        int sumComponent = 0;

        MPI_Recv(&sumComponent, 1, MPI_INT, MPI_ANY_SOURCE, waitFromID, MPI_COMM_WORLD, &status);

        sum += sumComponent;
    }

    if (id != 0)
    {
        MPI_Send(&sum, 1, MPI_INT, sendTo, id, MPI_COMM_WORLD);
        cout << "At " << id << " send to " << sendTo << " value: " << sum << endl;

        // MPI_Barrier(MPI_COMM_WORLD);
    }

    // 5. Print result at process 0
    if (id == 0)
    {
        cout << "= At process " << id << " -> sum = " << sum << endl;
        // MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
