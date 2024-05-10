#include <map>
#include <vector>

using namespace std;

const int N_DIRECTIONS = 4;
const double k2 = 1.0;
const double k3 = 0.4;
const double d = 1.0;

const int PROC_MAIN = 0;

enum DIRECTION
{
    LEFT,
    RIGHT,
    UP,
    BOTTOM
};

//- Particle ---------------------------------------------------------------------------------

class Particle
{
public:
    double x, y;
    double v_x, v_y, v_x_prev, v_y_prev, a_x, a_y, ep, ek;
    double box;
    int i, j;

    double dt;

    vector<Particle *> neighbors = vector<Particle *>(N_DIRECTIONS);

    Particle(double i, double j, double box = 0) : i(i), j(j), box(box)
    {
        this->x = i * 1.0;
        this->y = j * 1.0;
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

        this->v_x += this->dt * this->a_x;
        this->v_y += this->dt * this->a_y;
    }
    void update_a()
    {
        this->a_x = this->a_y = 0;

        double dx, dy, r, fx2, fy2, fx3, fy3;
        double u2 = 0, u3 = 0;

        for (int i = LEFT; i <= BOTTOM; ++i)
        {
            Particle *neighbor = neighbors[i];

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
};

//- Neighbor ---------------------------------------------------------------------------------

class Neighbor
{
public:
    // ----- Attributes -----
    int proc_id;
    map<int, Particle *> particlesSend, particlesRecv;

    // ----- Methods -----

    void signup_recv(int particleID, Particle *particle)
    {
        if (this->particlesRecv.find(particleID) == this->particlesRecv.end())
            this->particlesRecv[particleID] = particle;
    }
    void signup_send(int particleID, Particle *particle)
    {
        if (this->particlesSend.find(particleID) == this->particlesSend.end())
            this->particlesSend[particleID] = particle;
    }
};

//- Proces -----------------------------------------------------------------------------------

class Process
{
public:
    // ----- Attributes -----
    int proc_id;
    int num_iters, num_procs, nuc, num_each;
    double box;

    double ep_total, ek_total;
    vector<Particle *> particles, particlesNeighbor;
    unordered_map<int, Particle *> dict_id2par;
    unordered_map<int, Neighbor> procsNeighbor;

    // ----- Methods -----
    Process(
        int proc_id, int num_iters, int num_procs, int nuc,
        double box) : proc_id(proc_id), num_iters(num_iters), num_procs(num_procs), nuc(nuc), box(box)
    {
        // Determine which particle is handled by current proc
        int nat = this->nuc * this->nuc;
        this->num_each = (nat + this->num_procs - 1) / this->num_procs;

        int partStart = this->proc_id * this->num_each;
        int partEnd = min((this->proc_id + 1) * this->num_each, nat) - 1;

        // 1.Initialize particle
        for (int p = partStart; p <= partEnd; ++p)
        {
            // 1.1. Declare particle
            int i = p / this->nuc;
            int j = p % this->nuc;

            Particle *particle = new Particle(i, j, this->box);

            // 1.2. Push particle to vector and dictionary
            this->particles.push_back(particle);
            this->dict_id2par[this->ij2parID(particle)] = particle;
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
    }

    ~Process()
    {
        // TODO: HoangLe [May-10]: Remove all elements in particlesNeighbors
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

        return proc;
    }

    void
    loop()
    {
        for (int n = 0; n < this->num_iters; ++n)
        {
            this->ek_total = this->ep_total = 0;

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

                // Send energy to the process 0
                if (this->proc_id != PROC_MAIN)
                {
                    send_energies();
                }
                else
                {
                    receive_energies();
                }

                // Print energies
                if (this->proc_id == PROC_MAIN)
                {
                    print_energies();
                }
            }
        }
    }

    void exchange_neighbors()
    {
    }

    void send_energies()
    {
        //
    }

    void receive_energies()
    {
    }

    void print_energies()
    {
        // NOTE: HoangLe [May-09]: Run by process 0 only
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

int main(int argc, char const *argv[])
{
    int nuc = 4, num_procs = 3, proc_id = 1, num_iters = 10;
    double box = 1;
    Process process(proc_id, num_iters, num_procs, nuc, box);

    process.check();

    return 0;
}
