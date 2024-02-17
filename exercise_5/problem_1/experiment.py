import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from tqdm import trange

plt.style.use("seaborn-v0_8")
plt.rcParams.update({"font.size": 8})

re_pat = r"\:\s((\d|\.){2,})"
FILENAME = "./p1"
N = 10

if __name__ == "__main__":
    results = []

    assert Path(FILENAME).exists()

    for i in trange(N):
        # Run
        result = subprocess.run([FILENAME], stdout=subprocess.PIPE)

        # Extract result
        out = re.findall(re_pat, result.stdout.decode("utf-8"))
        results.append({"time": float(out[0][0]), "type": "formatted"})
        results.append({"time": float(out[1][0]), "type": "unformatted"})

    # Plot
    df = pd.DataFrame.from_records(results)
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    sns.violinplot(df, x="type", y="time", ax=ax)

    fig.savefig("p1_cputime.png", bbox_inches="tight")
