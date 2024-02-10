#include <cmath>
#include <vector>

using namespace std;

const int N_PER_PACK = 4;

typedef double double4 __attribute__((vector_size(N_PER_PACK * sizeof(double))));

const double4 vZeros = {0., 0., 0., 0.};

static inline double sumInternal(double4 x)
{
    double sum = 0;
    for (int i = 0; i < N_PER_PACK; ++i)
        sum += x[i];

    return sum;
}

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
void correlate(int ny, int nx, const float *data, float *result)
{
    vector<double> norm(ny * nx, 0.);

    // 1. row-wise 0-mean normalization
    // cout << "Step 1" << endl;

    for (int i = 0; i < ny; ++i)
    {
        double mean = 0;
        for (int j = 0; j < nx; ++j)
            mean += data[i * nx + j];

        mean /= nx;

        for (int j = 0; j < nx; ++j)
            norm[i * nx + j] = data[i * nx + j] - mean;
    }

    // 2. row-wise square-sum normalization
    // cout << "Step 2" << endl;

    for (int i = 0; i < ny; ++i)
    {
        double sq_sum = 0;
        for (int j = 0; j < nx; ++j)
        {
            sq_sum += norm[i * nx + j] * norm[i * nx + j];
        }
        sq_sum = sqrt(sq_sum);
        for (int j = 0; j < nx; ++j)
            norm[i * nx + j] = norm[i * nx + j] / sq_sum;
    }

    // 3. upper-triangular matmul
    // 3.1. Expand vector 'norm'
    int nPadded = nx;
    if (nPadded % N_PER_PACK != 0)
        nPadded = static_cast<int>(ceil(nx * 1.0 / N_PER_PACK)) * N_PER_PACK;
    int nPacks = nPadded / N_PER_PACK;

    // 3.2. Move 'norm' to new padded
    vector<double4> vnorm(ny * nPacks, vZeros);
    for (int i = 0; i < ny; ++i)
        for (int j = 0; j < nx; ++j)
            vnorm[nPacks * i + j / N_PER_PACK][j % N_PER_PACK] = norm[nx * i + j];

    // 3.3. Start calculating
    for (int i = 0; i < ny; ++i)
        for (int j = i; j < ny; ++j)
        {
            double4 vcor = vZeros;
            for (int k = 0; k < nPacks; ++k)
                vcor += vnorm[nPacks * i + k] * vnorm[nPacks * j + k];

            result[i * ny + j] = (float)sumInternal(vcor);
        }
}
