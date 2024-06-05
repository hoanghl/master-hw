#include <algorithm>
#include <vector>

using namespace std;

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in in[x + y*nx]
- for each pixel (x, y), store the median of the pixels (a, b) which satisfy
  max(x-hx, 0) <= a < min(x+hx+1, nx), max(y-hy, 0) <= b < min(y+hy+1, ny)
  in out[x + y*nx].
*/
void mf(int ny, int nx, int hy, int hx, const float *in, float *out)
{
    const int MAX = (2 * hx + 1) * (2 * hy + 1);

    for (int i = 0; i < ny; ++i)
    {
        for (int j = 0; j < nx; ++j)
        {
            int count = 0;
            vector<float> v = vector<float>(MAX);
            for (int y = max(0, i - hy); y <= min(i + hy, ny - 1); ++y)
            {
                for (int x = max(0, j - hx); x <= min(j + hx, nx - 1); ++x)
                {
                    v[count] = in[x + nx * y];
                    ++count;
                }
            }

            float median = 0;
            if (count % 2 == 0)
            {
                auto m = v.begin() + count / 2 - 1;
                nth_element(v.begin(), m, v.begin() + count);
                median += v[count / 2 - 1];

                m = v.begin() + count / 2;
                nth_element(v.begin(), m, v.begin() + count);
                median += v[count / 2];

                median *= 0.5;
            }
            else
            {
                auto m = v.begin() + count / 2;
                nth_element(v.begin(), m, v.begin() + count);
                median += v[count / 2];
            }

            out[j + nx * i] = median;
        }
    }
}
