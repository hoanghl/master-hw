#include <map>
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

//- Particle ---------------------------------------------------------------------------------

class Particle
{
public:
    double x, y;
    double v_x, v_y, v_x_prev, v_y_prev, a_x, a_y, ep, ek;
    int i, j;

    double box, dt, vsc;

    vector<Particle *> neighbors = vector<Particle *>(N_DIRECTIONS);

    Particle(double i = 0, double j = 0, double box = 0, double dt = 1, double vsc = 1) : i(i), j(j), box(box), dt(dt), vsc(vsc)
    {
        this->x = i * 1.0;
        this->y = j * 1.0;

        this->v_x = vsc * (getRand() - 0.5);
        this->v_y = vsc * (getRand() - 0.5);
    }

    void update_xy()
    {
        this->x += dt * this->v_x;
        this->y += dt * this->v_y;

        if (this->x < 0.0)
            this->x += box;
        else if (this->x >= box)
            this->x -= box;
        if (this->y < 0.0)
            this->y += box;
        else if (this->y >= box)
            this->y -= box;
    }
    void update_v()
    {
        this->v_x_prev = this->v_x;
        this->v_y_prev = this->v_y;

        this->v_x += this->dt * this->a_x;
        this->v_y += this->dt * this->a_y;
    }
    void update_a()
    {
        this->a_x = this->a_y = 0;

        double dx, dy, r, fx2, fy2, fx3, fy3;
        double u2 = 0, u3 = 0;

        Particle *neighbor;
        for (int i = LEFT; i <= BOTTOM; ++i)
        {
            neighbor = neighbors[i];

            dx = neighbor->x - this->x;
            if (dx < -box / 2.0)
                dx += box;
            if (dx >= box / 2.0)
                dx -= box;

            dy = neighbor->y - this->y;
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

            this->a_x += fx2 + fx3;
            this->a_y += fy2 + fy3;
        }

        this->ep = (u2 + u3) / 2.0;
    }

    void calc_ek()
    {
        double vx = (this->v_x_prev + this->v_x) / 2.0;
        double vy = (this->v_y_prev + this->v_y) / 2.0;
        this->ek = 1.0 / 2.0 * (pow(vx, 2) + pow(vy, 2));
    }

    double getRand()
    {
        return (double)rand() / RAND_MAX;
    }
};

//- Neighbor ---------------------------------------------------------------------------------

class Neighbor
{
public:
    // ----- Attributes -----
    int proc_id;
    map<int, Particle *> particlesSend, particlesRecv;

    vector<double> buffSendX, buffSendY, buffRecvX, buffRecvY;

    // ----- Methods -----

    void signup_recv(int particleID, Particle *particle)
    {
        if (this->particlesRecv.find(particleID) == this->particlesRecv.end())
        {
            this->particlesRecv[particleID] = particle;

            buffRecvX.push_back(0.0);
            buffRecvY.push_back(0.0);
        }
    }
    void signup_send(int particleID, Particle *particle)
    {
        if (this->particlesSend.find(particleID) == this->particlesSend.end())
        {
            this->particlesSend[particleID] = particle;

            buffSendX.push_back(0.0);
            buffSendY.push_back(0.0);
        }
    }

    void harvest()
    {
        int i = 0;
        for (auto &[_, p] : this->particlesSend)
        {
            this->buffSendX[i] = p->x;
            this->buffSendY[i] = p->y;

            ++i;
        }
    }
    void redistribute()
    {
        int i = 0;
        for (auto &[_, p] : this->particlesRecv)
        {
            this->buffRecvX[i] = p->x;
            this->buffRecvY[i] = p->y;

            ++i;
        }
    }
};

//- Proces -----------------------------------------------------------------------------------

class Process
{
public:
    // ----- Attributes -----
    int proc_id;
    int num_iters, num_procs, nuc, num_each;
    int interval;
    double box, dt, vsc;

    double ep_total, ek_total;
    vector<Particle *> particles, particlesNeighbor;
    unordered_map<int, Particle *> dict_id2par;
    unordered_map<int, Neighbor> procsNeighbor;

    // MPI_Comm comm;
    // MPI_Status status;

    // ----- Methods -----
    Process(
        int proc_id,
        int num_iters,
        int num_procs,
        int nuc,
        double box,
        double dt,
        double vsc,
        int interval = 100) : proc_id(proc_id), num_iters(num_iters), num_procs(num_procs), nuc(nuc), box(box), dt(dt), vsc(vsc), interval(interval)
    {
        // Some checkings
        if (proc_id < 0 || proc_id >= num_procs)
        {
            printf("Err: proc_id not in range [0, %d)\n", num_procs);
            exit(1);
        }

        // Determine which particle is handled by current proc
        int nat = this->nuc * this->nuc;
        this->num_each = nat / this->num_procs;

        int partStart = proc_id * num_each;
        int partEnd = (proc_id + 1) * this->num_each - 1;
        if (proc_id == num_procs - 1)
            partEnd = nat - 1;

        // 1. Initialize particle

        // 1.1. Create new particles
        for (int p = partStart; p <= partEnd; ++p)
        {
            // 1.1. Declare particle
            int i = p / this->nuc;
            int j = p % this->nuc;

            Particle *particle = new Particle(i, j, box, dt, vsc);

            // 1.2. Push particle to vector and dictionary
            this->particles.push_back(particle);
            this->dict_id2par[this->ij2parID(particle)] = particle;
        }

        // 1.2. Update particles' velocity
        double sum_vx = 0, sum_vy = 0;
        for (Particle *p : this->particles)
        {
            sum_vx += p->v_x;
            sum_vy += p->v_y;
        }
        sum_vx /= nat;
        sum_vy /= nat;
        for (Particle *p : this->particles)
        {
            p->x -= sum_vx;
            p->y -= sum_vy;
        }

        // 2. Define neighbor (neigbor particles, processes) in 4 directions
        for (Particle *p : this->particles)
        {
            for (int direction = DIRECTION::LEFT; direction <= DIRECTION::BOTTOM; ++direction)
            {

                // 2.1. Determine the neighbor particle and corresponding in each direction

                int iNeighbor = 0, jNeighbor = 0;
                this->getNeighbor(p->i, p->j, direction, iNeighbor, jNeighbor);
                int idNeighbor = ij2parID(iNeighbor, jNeighbor);

                // Get neighbor particle (create new one if not existed)

                Particle *pNeighbor = nullptr;
                if (this->dict_id2par.find(idNeighbor) == this->dict_id2par.end())
                {
                    pNeighbor = new Particle(iNeighbor, jNeighbor);

                    this->particlesNeighbor.push_back(pNeighbor);
                    this->dict_id2par[idNeighbor] = pNeighbor;
                }
                else
                {
                    pNeighbor = this->dict_id2par[idNeighbor];
                }

                int idProcNeighbor = ij2proc(iNeighbor, jNeighbor);

                // If the process handling the neighbor particle is not the current process
                if (idProcNeighbor != this->proc_id)
                {
                    // Get neigbor process (create new one if not existed)
                    Neighbor *procNeighbor = nullptr;
                    if (this->procsNeighbor.find(idProcNeighbor) == this->procsNeighbor.end())
                    {
                        this->procsNeighbor[idProcNeighbor] = {.proc_id = idProcNeighbor};
                    }
                    procNeighbor = &this->procsNeighbor[idProcNeighbor];

                    // Sign up pNeighbor and p to procNeighbor
                    procNeighbor->signup_recv(idNeighbor, pNeighbor);
                    procNeighbor->signup_send(ij2parID(p), p);
                }

                // 2.2. Assign neighbor to each particle
                p->neighbors[direction] = pNeighbor;
            }
        }

        // 3. Initialize MPI
    }

    ~Process()
    {
        for (int i = 0; i < this->particlesNeighbor.size(); ++i)
        {
            delete this->particlesNeighbor[i];
            delete this->particles[i];
        }
    }

    int ij2parID(int i, int j)
    {
        return i * this->nuc + j;
    }
    int ij2parID(Particle *p)
    {
        return p->i * this->nuc + p->j;
    }

    void getNeighbor(int i, int j, int direction, int &iNeighbor, int &jNeighbor)
    {
        iNeighbor = i;
        jNeighbor = j;

        switch (direction)
        {
        case DIRECTION::UP:
            iNeighbor--;
            if (iNeighbor < 0)
                iNeighbor = this->nuc - 1;
            break;
        case DIRECTION::BOTTOM:
            iNeighbor++;
            if (iNeighbor >= nuc)
                iNeighbor = 0;
            break;
        case DIRECTION::LEFT:
            jNeighbor--;
            if (jNeighbor < 0)
                jNeighbor = nuc - 1;
            break;
        case DIRECTION::RIGHT:
            jNeighbor++;
            if (jNeighbor >= nuc)
                jNeighbor = 0;
            break;
        }
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

    void
    loop()
    {
        printf("Start looping with num_iters = %d...\n", this->num_iters);

        for (int n = 0; n < this->num_iters; ++n)
        {
            if (this->proc_id == PROC_LEADER)
                printf("n = %d\n", n);

            this->ek_total = this->ep_total = 0;
            double epFinal = 0, ekFinal = 0;

            for (Particle *p : this->particles)
            {
                // Exchange the coordinate to the neighbors
                this->exchange_neighbors();

                // Calulcate acceleration, velocity, coordinate and energies
                p->update_a();
                p->update_v();
                p->update_xy();
                p->calc_ek();

                // Accumulate total energy
                this->ek_total += p->ek;
                this->ep_total += p->ep;

                // Process 0 () accumulates the total_ep, total_ek sent from followers'
                reduce_energies(epFinal, ekFinal);
            }

            // Print energies
            print_energies(n, epFinal, ekFinal);
        }
    }

    void exchange_neighbors()
    {
        // Create thins
        int num_reqs = 4 * this->procsNeighbor.size();
        MPI_Request *requests = new MPI_Request[num_reqs];

        // Start iterating procNeighbor and exchanging
        int i = 0;
        for (auto &[_, p] : this->procsNeighbor)
        {
            // p.harvest();

            // Use MPI API
            MPI_Isend(&p.buffSendX, p.buffSendX.size(), MPI_DOUBLE, p.proc_id, TAG_X, MPI_COMM_WORLD, &requests[i]);
            MPI_Isend(&p.buffSendY, p.buffSendY.size(), MPI_DOUBLE, p.proc_id, TAG_Y, MPI_COMM_WORLD, &requests[i + 1]);

            MPI_Irecv(&p.buffRecvX, p.buffRecvX.size(), MPI_DOUBLE, p.proc_id, TAG_X, MPI_COMM_WORLD, &requests[i + 2]);
            MPI_Irecv(&p.buffRecvY, p.buffRecvY.size(), MPI_DOUBLE, p.proc_id, TAG_Y, MPI_COMM_WORLD, &requests[i + 3]);

            i += 4;
        }

        MPI_Waitall(num_reqs, requests, MPI_STATUSES_IGNORE);

        // After exchanging, each neighbor redistributes the received coordinates
        // for (auto &[_, p] : this->procsNeighbor)
        // {
        //     p.redistribute();
        // }

        // Remove pointers
        delete[] requests;
    }

    void reduce_energies(double &epFinal, double &ekFinal)
    {
        // MPI_Reduce(&this->ep_total, &epFinal, 1, MPI_DOUBLE, MPI_SUM, PROC_LEADER, MPI_COMM_WORLD);
        // MPI_Reduce(&this->ek_total, &ekFinal, 1, MPI_DOUBLE, MPI_SUM, PROC_LEADER, MPI_COMM_WORLD);
    }

    void print_energies(int n, double epFinal, double ekFinal)
    {
        if (this->proc_id == PROC_LEADER)
        {
            if (n % this->interval == 0)
                printf("%20.10g %20.10g %20.10g %20.10g\n", dt * n, epFinal + ekFinal, epFinal, ekFinal);
        }
    }

    void check()
    {
        // Check handling particles
        printf("== handling:\n");
        for (Particle *p : this->particles)
        {
            printf("particle %2d: \n", ij2parID(p->i, p->j));
            printf("   - i = %2d - j = %2d - x = %.2f - y = %.2f\n", p->i, p->j, p->x, p->y);
            for (int d = LEFT; d <= BOTTOM; ++d)
            {
                printf("   - %d: i = %2d - j = %2d\n", d, p->neighbors[d]->i, p->neighbors[d]->j);
            }
        }

        // Check neighbor particles
        printf("== Neighbor particles:\n");
        for (Particle *p : this->particlesNeighbor)
        {
            printf("particle %2d: \n", ij2parID(p->i, p->j));
            printf("   - i = %2d - j = %2d - x = %.2f - y = %.2f\n", p->i, p->j, p->x, p->y);
        }

        // Check neighbor processes
        printf("== Neighbor processes:\n");
        for (auto &[k, v] : this->procsNeighbor)
        {
            printf("neigh_proc %d\n: ", k);

            printf("--- Receive:\n");
            for (auto &[pID, p] : v.particlesRecv)
            {
                printf("particle %2d: \n", pID);
                printf("   - i = %2d - j = %2d - x = %.2f - y = %.2f\n", p->i, p->j, p->x, p->y);
            }

            printf("--- Send:\n");
            for (auto &[pID, p] : v.particlesSend)
            {
                printf("particle %2d: \n", pID);
                printf("   - i = %2d - j = %2d - x = %.2f - y = %.2f\n", p->i, p->j, p->x, p->y);
            }
        }
    }
};

int main(int argc, char *argv[])
{
    // TODO: HoangLe [May-11]: Write code to receive arguments

    int num_procs, proc_id;

    int nuc = 100, num_iters = 100000;
    double box = nuc, dt = 0.001, vsc = 0.5;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm comm;
    MPI_Status status;

    int dim[] = {num_procs}, periods[] = {1}, reorder = false;
    MPI_Cart_create(MPI_COMM_WORLD, 1, dim, periods, reorder, &comm);
    Process process(proc_id, num_iters, num_procs, nuc, box, dt, vsc);

    printf("proc_id = %d , num_procs = %d\n", proc_id, num_procs);

    process.loop();

    // process.check();

    return 0;
}
