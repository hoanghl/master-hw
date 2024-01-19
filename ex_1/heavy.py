import multiprocessing

import numpy as np
from tqdm.contrib.itertools import product


def f():
    len_a = np.random.randint(1000, 10000)
    len_b = np.random.randint(1000, 10000)
    len_c = np.random.randint(1000, 10000)

    a = np.random.randint(10, 1000, (len_a))
    b = np.random.randint(10, 1000, (len_b))
    c = np.random.randint(10, 1000, (len_c))

    results = []

    for xa, xb, xc in product(a, b, c):
        results.append(xa * xb - np.sqrt(xc))

    np.save("results.npy", np.array(results, dtype=np.float32))


if __name__ == "__main__":
    p1 = multiprocessing.Process(target=f)
    p2 = multiprocessing.Process(target=f)

    p1.start()
    p2.start()

    p1.join()
    p2.join()
