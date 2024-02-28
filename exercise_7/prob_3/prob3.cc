#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <random>
#include <sstream>
#include <string>
#include <time.h>

using namespace std;

const int TAG = 0;
const string FILENAME = "ex7p3.dat";
const int SIZE = 500000;
const int ID_LEADER = 0;
const int N_MIN_PER_PROC = 4;

double getRand()
{
    return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}

void log_send(int fromRoutine, int destId, int msg)
{
    cout << "= Routine " << fromRoutine << ": Sent to " << destId << ": " << msg << endl;
}

void log_recv(int fromRoutine, MPI_Status status, int msg)
{
    cout << "= Routine " << fromRoutine << ": Received from " << status.MPI_SOURCE << ": " << msg << endl;
}

void log_recv(int fromRoutine, MPI_Status status, double msg)
{
    cout << "= Routine " << fromRoutine << ": Received from " << status.MPI_SOURCE << ": " << setprecision(6) << msg << endl;
}

int getNumProc(int n, int nTasks, int idx)
{
    int nTasksUsable = nTasks - 1;

    int nPerProc = std::max<int>(N_MIN_PER_PROC, n / nTasksUsable);

    int num;
    if (idx == 0)
        num = 0;
    else if (n > (idx - 1) * nPerProc)
    {
        if (idx == nTasksUsable)
            num = n - (idx - 1) * nPerProc;
        else
            num = std::min<int>(nPerProc, n - (idx - 1) * nPerProc);
    }

    else
        num = 0;
    return num;
}

int main(int argc, char *argv[])
{
    // Initialize MP
    int id, nTasks, destId, rc;
    MPI_Status status;
    int pnlen;
    char pname[MPI_MAX_PROCESSOR_NAME];

    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS)
    {
        printf("MPI initialization failed\n");
        exit(1);
    }
    rc = MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &id);
    rc = MPI_Get_processor_name(pname, &pnlen);

    int nTasksUsable = nTasks - 1; // Process 0 wont process data, it is as the leader

    if (id == 0)
    {
        // 1. Read data
        fstream file(FILENAME, ios::in);
        float arr[SIZE];
        int n = 0;
        if (file.is_open())
        {
            string line;
            while (getline(file, line))
            {
                istringstream iss(line);
                iss >> arr[n];
                ++n;
            }
            file.close();
        }

        // 2. Send array size to other
        for (int destId = 1; destId < nTasks; ++destId)
            MPI_Send(&n, 1, MPI_INT, destId, TAG, MPI_COMM_WORLD);

        // 2. Send to other
        int *displs = new int[nTasks];
        int *sendcounts = new int[nTasks];
        int displacement = 0;
        for (int i = 1; i <= nTasksUsable; ++i)
        {
            int num = getNumProc(n, nTasks, i);
            if (num == 0)
                break;

            sendcounts[i] = num;
            displs[i] = displacement;

            displacement += num;

            cout << "disp[" << i << "] = " << displs[i] << " & sendcounts[" << i << "] = " << sendcounts[i] << endl;
        }
        MPI_Scatterv(arr, sendcounts, displs, MPI_FLOAT,
                     NULL, 0, MPI_FLOAT, 0, MPI_COMM_WORLD);

        // 3. Receive sum from each sub-process and calculate mean
        int countSentProc = 0;
        double mean = 0, sumRecv = 0;
        while (countSentProc < nTasksUsable)
        {
            MPI_Recv(&sumRecv, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            log_recv(id, status, sumRecv);

            mean += sumRecv;
            ++countSentProc;
        }
        mean /= n;

        // 4. Populate 'mean' back to sub-processes to calculate variance
        for (int destId = 1; destId < nTasks; ++destId)
            MPI_Send(&mean, 1, MPI_DOUBLE, destId, TAG, MPI_COMM_WORLD);

        // 5. Receive variance from sub-processes and calculate the final variance
        countSentProc = 0;
        double variance = 0, varianceRecv = 0;
        while (countSentProc < nTasksUsable)
        {
            MPI_Recv(&varianceRecv, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            // log_recv(id, status, sumRecv);

            variance += varianceRecv;
            ++countSentProc;
        }
        variance /= n;

        cout << "mean = " << setprecision(8) << mean << "  -  variance = " << setprecision(8) << variance << endl;

        MPI_Barrier(MPI_COMM_WORLD);
    }
    else
    {
        // 1. Receive configurations including array size and nTasksUsable
        int n = 0;
        MPI_Recv(&n, 1, MPI_INT, ID_LEADER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        // log_recv(id, status, n);

        // 2. Calculate recvcount for each sub-process
        int recvcount = getNumProc(n, nTasks, id);

        // 3. Receive data from process 0
        float *recvbuf = new float[recvcount];
        MPI_Scatterv(NULL, NULL, NULL, MPI_FLOAT,
                     recvbuf, recvcount, MPI_FLOAT, ID_LEADER, MPI_COMM_WORLD);

        // 4. Calculate sum and send back to process 0
        double sumProc = 0;
        for (int i = 0; i < recvcount; ++i)
        {
            sumProc += static_cast<double>(recvbuf[i]);
        }

        MPI_Send(&sumProc, 1, MPI_DOUBLE, ID_LEADER, TAG, MPI_COMM_WORLD);

        // 5. Receive mean from process 0, calculate and send back variance to process 0
        double mean = 0;
        MPI_Recv(&mean, 1, MPI_DOUBLE, ID_LEADER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        double variance = 0;
        for (int i = 0; i < recvcount; ++i)
        {
            variance += pow(static_cast<double>(recvbuf[i]) - mean, 2);
        }

        MPI_Send(&variance, 1, MPI_DOUBLE, ID_LEADER, TAG, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
    }
}
