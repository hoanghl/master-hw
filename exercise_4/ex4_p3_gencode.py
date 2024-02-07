import os
import re
import subprocess
from typing import Literal

re_pat = r"\:\s((\d|\.)*)"

CC = "clang++"
LOWER = """
}
    clock_t end = clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    cout << "CPU time (second): " << setprecision(12) << elapsed << endl;

    return 0;
}
"""


def get_upper(n: int):
    return f"""
    #include <iomanip>
    #include <iostream>
    #include <random>
    #include <time.h>

    #include "ex4.hpp"

    using namespace std;

    int main(int argc, char const *argv[])
    {{
        const int N = 16777216;
        int *arr = new int[N];
        for (int i = 0; i < N; ++i)
            arr[i] = rand();
        int *sub = new int[N];

        clock_t start = clock();
        for (int i = 0; i < N; i = i + {n})
        {{
    """


def get_filename(n: int, filetype: Literal["code", "exe"]):
    extension = ".cc" if filetype == "code" else ""
    return f"ex4_p3_{n}{extension}"


if __name__ == "__main__":
    results = []

    for i in range(1):
        n = 2**i

        # Create code
        filename_code = get_filename(n, "code")
        upper = get_upper(n)

        with open(filename_code, "w+") as file:
            file.write(upper)
            file.write("\n")
            for j in range(n):
                file.write(f"sub[i + {j}] = arr[i + {j}] * 2;")
                file.write("\n")
            file.write(LOWER)

        # Compile & Run
        filename_exe = get_filename(n, "exe")
        for flag_opt in ["-O0", "-O3"]:
            ## Compile
            subprocess.run([CC, flag_opt, "-o", filename_exe, filename_code])

            # Run
            result = subprocess.run([f"./{filename_exe}"], stdout=subprocess.PIPE)

            # Extract result
            time = float(re.findall(re_pat, result.stdout.decode("utf-8"))[0][0])

            results.append({"n": n, "flag": flag_opt, "time": time})

            os.remove(filename_code)
            os.remove(filename_exe)

        break

    print(results)
