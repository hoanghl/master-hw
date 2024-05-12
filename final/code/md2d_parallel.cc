#include <iostream>
#include <map>
#include <math.h>
#include <vector>

#include <mpi.h>

using namespace std;

const int N_DIRECTIONS = 4;
const double k2 = 1.0;
const double k3 = 0.4;
const double d = 1.0;

const int TAG_X = 10;
const int TAG_Y = 20;

const int PROC_LEADER = 0;

enum DIRECTION
{
    LEFT,
    UP,
    RIGHT,
    BOTTOM
};

//- Neighbor ---------------------------------------------------------------------------------

class Neighbor
{
public:
    // ----- Attributes -----
    int proc_id;
    int minParIDSend = -1, maxParIDSend = -1, minParIDRecv = -1, maxParIDRecv = -1;

    // ----- Methods -----

    void signup_recv(int parID)
    {
        if (this->minParIDRecv == -1)
            this->minParIDRecv = this->maxParIDRecv = parID;
        else
        {
            if (parID < this->minParIDRecv)
                this->minParIDRecv = parID;
            else if (parID > this->maxParIDRecv)
                this->maxParIDRecv = parID;
        }
    }
    void signup_send(int parID)
    {
        if (this->minParIDSend == -1)
            this->minParIDSend = this->maxParIDSend = parID;
        else
        {
            if (parID < this->minParIDSend)
                this->minParIDSend = parID;
            else if (parID > this->maxParIDSend)
                this->maxParIDSend = parID;
        }
    }
};

//- Proces -----------------------------------------------------------------------------------

class Process
{

public:
    // ----- Attributes -----
    int proc_id;
    int maxt, num_procs, nuc;
    int eout;
    double box, dt, vsc;

    double ep_total, ek_total;
    int num_each, partStart, partEnd;

    vector<double> x, y, v_x, v_y, a_x, a_y, ep, ek;

    unordered_map<int, Neighbor> procsNeighbor;

    MPI_Comm comm;
    MPI_Status status;

    // ----- Methods -----
    Process(
        int proc_id,
        int maxt,
        int num_procs,
        int nuc,
        double box,
        double dt,
        double vsc,
        int eout = 100) : proc_id(proc_id), maxt(maxt), num_procs(num_procs), nuc(nuc), box(box), dt(dt), vsc(vsc), eout(eout)
    {
        // Some checkings
        if (proc_id < 0 || proc_id >= num_procs)
        {
            printf("Err: proc_id not in range [0, %d)\n", num_procs);
            exit(1);
        }

        // 1. Determine which particle is handled by current proc
        int nat = this->nuc * this->nuc;
        this->num_each = nat / this->num_procs;

        this->partStart = proc_id * num_each;
        this->partEnd = (proc_id + 1) * this->num_each - 1;
        if (proc_id == num_procs - 1)
            this->partEnd = nat - 1;

        // 2. Initialize particle: coordinate, velocity
        this->x = vector<double>(nat, 0.0);
        this->y = vector<double>(nat, 0.0);
        this->v_x = vector<double>(nat, 0.0);
        this->v_y = vector<double>(nat, 0.0);
        this->a_x = vector<double>(nat, 0.0);
        this->a_y = vector<double>(nat, 0.0);
        this->ep = vector<double>(nat, 0.0);
        this->ek = vector<double>(nat, 0.0);

        // 2.1. Init coordinate, velocity
        for (int parID = this->partStart; parID <= this->partEnd; ++parID)
        {
            int i, j;
            this->parID2ij(parID, i, j);

            this->x[parID] = i * 1.0;
            this->y[parID] = j * 1.0;

            this->v_x[parID] = this->vsc * (getRand() - 0.5);
            this->v_y[parID] = this->vsc * (getRand() - 0.5);
        }

        // 2.2. Update particles' coordinate
        double sum_vx = 0, sum_vy = 0;
        for (int parID = this->partStart; parID <= this->partEnd; ++parID)
        {
            sum_vx += this->v_x[parID];
            sum_vy += this->v_y[parID];
        }
        sum_vx /= nat;
        sum_vy /= nat;
        for (int parID = this->partStart; parID <= this->partEnd; ++parID)
        {
            this->v_x[parID] -= sum_vx;
            this->v_y[parID] -= sum_vy;
        }

        // 2. Define neighbor (neigbor particles, processes) in 4 directions
        for (int parID = this->partStart; parID <= this->partEnd; ++parID)
        {
            int i, j;
            parID2ij(parID, i, j);

            for (int direction = DIRECTION::LEFT; direction <= DIRECTION::BOTTOM; ++direction)
            {

                // 2.1. Determine the neighbor particle and corresponding process in each direction

                int parIDNb;
                this->getNeighbor(parID, direction, parIDNb);
                int procIDNb = parID2proc(parIDNb);

                // If the process handling the neighbor particle is not the current process
                if (procIDNb != this->proc_id)
                {
                    // Get neigbor process (create new one if not existed)
                    Neighbor *procNeighbor = nullptr;
                    if (this->procsNeighbor.find(procIDNb) == this->procsNeighbor.end())
                    {
                        this->procsNeighbor[procIDNb] = {.proc_id = procIDNb};
                    }
                    procNeighbor = &this->procsNeighbor[procIDNb];

                    // Sign up pNeighbor and p to procNeighbor
                    procNeighbor->signup_recv(parIDNb);
                    procNeighbor->signup_send(parID);
                }
            }
        }
    }

    double getRand()
    {
        return (double)rand() / RAND_MAX;
    }

    int ij2parID(int i, int j)
    {
        return i * this->nuc + j;
    }

    void parID2ij(int parID, int &i, int &j)
    {
        i = parID / this->nuc;
        j = parID % this->nuc;
    }

    void getNeighbor(int i, int j, int direction, int &iNb, int &jNb)
    {
        iNb = i;
        jNb = j;

        switch (direction)
        {
        case DIRECTION::UP:
            iNb--;
            if (iNb < 0)
                iNb = this->nuc - 1;
            break;
        case DIRECTION::BOTTOM:
            iNb++;
            if (iNb >= this->nuc)
                iNb = 0;
            break;
        case DIRECTION::LEFT:
            jNb--;
            if (jNb < 0)
                jNb = this->nuc - 1;
            break;
        case DIRECTION::RIGHT:
            jNb++;
            if (jNb >= this->nuc)
                jNb = 0;
            break;
        }
    }

    void getNeighbor(int parID, int direction, int &parIDNb)
    {
        int i, j;
        parID2ij(parID, i, j);
        int iNb = i, jNb = j;

        switch (direction)
        {
        case DIRECTION::UP:
            iNb--;
            if (iNb < 0)
                iNb = this->nuc - 1;
            break;
        case DIRECTION::BOTTOM:
            iNb++;
            if (iNb >= nuc)
                iNb = 0;
            break;
        case DIRECTION::LEFT:
            jNb--;
            if (jNb < 0)
                jNb = nuc - 1;
            break;
        case DIRECTION::RIGHT:
            jNb++;
            if (jNb >= nuc)
                jNb = 0;
            break;
        }

        parIDNb = ij2parID(iNb, jNb);
    }

    int ij2proc(int i, int j)
    {
        int nat = this->nuc * this->nuc;
        int idParticle = ij2parID(i, j);
        if (idParticle < 0 || idParticle >= nat)
        {
            // TODO: HoangLe [Apr-28]: Yield error
            return -1;
        }

        int proc = idParticle / this->num_each;
        if (proc >= this->num_procs)
        {
            proc = this->num_procs - 1;
        }

        return proc;
    }

    int parID2proc(int parID)
    {
        int nat = this->nuc * this->nuc;

        if (parID < 0 || parID >= nat)
        {
            // TODO: HoangLe [Apr-28]: Yield error
            return -1;
        }

        int proc = parID / this->num_each;
        if (proc >= this->num_procs)
        {
            proc = this->num_procs - 1;
        }

        return proc;
    }

    // --- Particle update ---------
    void update_a(int parID)
    {
        double dx, dy, r, fx2, fy2, fx3, fy3;
        double u2 = 0, u3 = 0;

        int i, j;
        this->parID2ij(parID, i, j);
        for (int direction = LEFT; direction <= BOTTOM; ++direction)
        {
            int parIDNb;
            this->getNeighbor(parID, direction, parIDNb);

            dx = this->x[parIDNb] - this->x[parID];
            if (dx < -box / 2.0)
                dx += box;
            if (dx >= box / 2.0)
                dx -= box;

            dy = this->y[parIDNb] - this->y[parID];
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

            this->a_x[parID] += fx2 + fx3;
            this->a_y[parID] += fy2 + fy3;
        }

        this->ep[parID] = (u2 + u3) / 2.0;
    }

    void update_v(int parID)
    {
        this->v_x[parID] += this->dt * this->a_x[parID];
        this->v_y[parID] += this->dt * this->a_y[parID];
    }

    void update_xy(int parID)
    {
        this->x[parID] += this->dt * this->v_x[parID];
        this->y[parID] += this->dt * this->v_y[parID];

        if (this->x[parID] < 0.0)
            this->x[parID] += box;
        else if (this->x[parID] >= box)
            this->x[parID] -= box;
        if (this->y[parID] < 0.0)
            this->y[parID] += box;
        else if (this->y[parID] >= box)
            this->y[parID] -= box;
    }

    void calc_ek(int parID, double v_x_prev, double v_y_prev)
    {
        double vx = (v_x_prev + this->v_x[parID]) / 2.0;
        double vy = (v_y_prev + this->v_y[parID]) / 2.0;
        this->ek[parID] = 1.0 / 2.0 * (pow(vx, 2) + pow(vy, 2));
    }

    void loop()
    {
        for (int n = 0; n < this->maxt; ++n)
        {
            // if (this->proc_id == PROC_LEADER)
            //     printf("n = %d\n", n);

            this->ek_total = this->ep_total = 0;
            double epFinal = 0, ekFinal = 0;

            fill(this->a_x.begin(), this->a_x.end(), 0);
            fill(this->a_y.begin(), this->a_y.end(), 0);

            // Exchange the coordinate to the neighbors
            this->exchange_neighbors();

            for (int parID = this->partStart; parID <= this->partEnd; ++parID)
            {

                // Calulcate acceleration, velocity, coordinate and energies
                double v_x_prev = this->v_x[parID], v_y_prev = this->v_y[parID];

                update_a(parID);
                update_v(parID);
                update_xy(parID);
                calc_ek(parID, v_x_prev, v_y_prev);

                // Accumulate total energy
                this->ek_total += this->ek[parID];
                this->ep_total += this->ep[parID];
            }

            // Process 0 () accumulates the total_ep, total_ek sent from followers'
            reduce_energies(epFinal, ekFinal);

            // Print energies
            print_energies(n, epFinal, ekFinal);
        }
    }

    void exchange_neighbors()
    {

        // Create thins
        int N = 4;
        int num_reqs = N * this->procsNeighbor.size();
        MPI_Request *requests = new MPI_Request[num_reqs];

        // Start iterating procNeighbor and exchanging
        int i = 0;
        for (auto &[_, p] : this->procsNeighbor)
        {
            // Use MPI API
            MPI_Isend(&this->x[p.minParIDSend], p.maxParIDSend - p.minParIDSend + 1, MPI_DOUBLE, p.proc_id, TAG_X, MPI_COMM_WORLD, &requests[i + 0]);
            MPI_Irecv(&this->x[p.minParIDRecv], p.maxParIDRecv - p.minParIDRecv + 1, MPI_DOUBLE, p.proc_id, TAG_X, MPI_COMM_WORLD, &requests[i + 1]);

            MPI_Isend(&this->y[p.minParIDSend], p.maxParIDSend - p.minParIDSend + 1, MPI_DOUBLE, p.proc_id, TAG_Y, MPI_COMM_WORLD, &requests[i + 2]);
            MPI_Irecv(&this->y[p.minParIDRecv], p.maxParIDRecv - p.minParIDRecv + 1, MPI_DOUBLE, p.proc_id, TAG_Y, MPI_COMM_WORLD, &requests[i + 3]);

            i += N;
        }
        MPI_Waitall(num_reqs, requests, MPI_STATUSES_IGNORE);

        // Remove pointers
        delete[] requests;
    }

    void reduce_energies(double &epFinal, double &ekFinal)
    {
        MPI_Reduce(&this->ep_total, &epFinal, 1, MPI_DOUBLE, MPI_SUM, PROC_LEADER, MPI_COMM_WORLD);
        MPI_Reduce(&this->ek_total, &ekFinal, 1, MPI_DOUBLE, MPI_SUM, PROC_LEADER, MPI_COMM_WORLD);
    }

    void print_energies(int n, double epFinal, double ekFinal)
    {
        if (this->proc_id == PROC_LEADER)
        {
            if (this->eout > 0 && n % this->eout == 0)
                printf("%20.10g %20.10g %20.10g %20.10g\n", dt * n, epFinal + ekFinal, epFinal, ekFinal);
        }
    }

    void check()
    {
        int i, j, iNb, jNb;

        // Check handling particles
        printf("== handling:\n");
        for (int parID = this->partStart; parID <= this->partEnd; ++parID)
        {
            parID2ij(parID, i, j);

            printf("particle %2d: \n", parID);
            printf("   - i = %2d - j = %2d - x = %.2f - y = %.2f\n", i, j, this->x[parID], this->y[parID]);
            for (int d = LEFT; d <= BOTTOM; ++d)
            {
                this->getNeighbor(i, j, d, iNb, jNb);
                printf("   - %d: i = %2d - j = %2d\n", d, iNb, jNb);
            }
        }

        // Check neighbor processes
        printf("== Neighbor processes:\n");
        for (auto &[k, v] : this->procsNeighbor)
        {
            printf("neigh_proc %d:\n", k);

            printf("--- Receive: min: %d - max: %d\n", v.minParIDRecv, v.maxParIDRecv);
            printf("--- Send:    min: %d - max: %d\n", v.minParIDSend, v.maxParIDSend);
        }
    }
};

int main(int argc, char *argv[])
{
    double dt;  // time step
    double vsc; // mean initial velocity
    double box; // system size
    int nuc;    // number of unit cells
    int nat;    // number of atoms
    int maxt;   // number of time steps simulated
    int eout;   // energy output eout

    int num_procs, proc_id;

    if (argc < 5 || argc > 6)
    {
        cerr << "usage: " << argv[0] << " nuc dt maxt vsc [eout [coout]]\n";
        cerr << "    nuc   = number of unit cells\n";
        cerr << "    dt    = time step\n";
        cerr << "    maxt  = number of time steps in simulation\n";
        cerr << "    vsc   = mean velocity of atoms in the beginning ('temperature')\n";
        cerr << "    eout  = eout for printing energies to stdout\n";
        return (128);
    }

    nuc = atoi(*++argv);
    dt = atof(*++argv);
    maxt = atoi(*++argv);
    vsc = atof(*++argv);
    eout = atof(*++argv);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    box = nuc;

    Process process(proc_id, maxt, num_procs, nuc, box, dt, vsc, eout);
    // process.check();

    auto t0 = std::chrono::system_clock::now();

    process.loop();

    auto t1 = std::chrono::system_clock::now();
    auto wct = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    if (proc_id == PROC_LEADER)
    {
        double wtime = wct.count() / 1000.0;
        if (eout > 0)
            printf("Wall clock time: %.4f seconds\n", wtime);
        else
            printf("%f\n", wtime);
    }

    MPI_Finalize();

    return 0;
}
