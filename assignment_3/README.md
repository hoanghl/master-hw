### Run experiment

In file `hw3_p2.py`, the function `experiment()` runs experiments. To run experiments, run as follow:

```python
    list_p = [0.1, 0.5]
    results = experiment(G, list_p, r=1, k=5, num_init_cand=7, num_cand_set=3, sigma=5)
```

where

| argument        | meaning                                                                  | value                    |
| --------------- | ------------------------------------------------------------------------ | ------------------------ |
| `G`             | The graph                                                                | Graph read by `networkx` |
| `list_p`        | List of infection probability                                            | list of float            |
| `r`             | Recovery probability                                                     | single float value       |
| `k`             | Number of greedy algorithm runs                                          | single integer value     |
| `num_init_cand` | Number of initial candidates (chosen without replacement from `G.nodes`) | single integer value     |
| `num_cand_set`  | Number of candidate sets                                                 | single integer value     |
| `sigma`         | Number of runs of calculating expected spread                            | single integer value     |
