import os
import json
import pandas as pd
import matplotlib.pyplot as plt

# Load stats JSON files
stats = []
for fname in ["global_stats.json", "local_stats.json"]:
    path = f"../build/{fname}"
    if os.path.exists(path):
        with open(path) as f:
            data = json.load(f)
            data["method"] = data.get("method", fname.split("_")[0])
            stats.append(data)

# Create DataFrame of stats
df_stats = pd.DataFrame(stats)

# Bar charts for each metric
metrics = ["score", "time_ms", "identity", "coverage"]
for metric in metrics:
    plt.figure()
    df_stats.plot(x="method", y=metric, kind="bar", legend=False)
    plt.title(f"{metric.capitalize()} by Method")
    plt.xlabel("Method")
    plt.ylabel(metric.replace("_", " ").capitalize())
    plt.tight_layout()
    plt.show()

# Parse per-block identity from global_alignment.txt
blocks = []
align_path = "../build/global_alignment.txt"
with open(align_path) as f:
    lines = f.readlines()

for idx, line in enumerate(lines):
    if line.startswith("Position:"):
        parts = line.strip().split()
        start, end = int(parts[1]), int(parts[3])
        seq1 = lines[idx+1].strip()
        seq2 = lines[idx+2].strip()
        total = len(seq1)
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        blocks.append({
            "start": start,
            "end": end,
            "identity": matches / total
        })

df_blocks = pd.DataFrame(blocks)

# Line plot of per-block identity
plt.figure()
x = (df_blocks["start"] + df_blocks["end"]) / 2
y = df_blocks["identity"]
plt.plot(x, y, marker="o", linewidth=1)
plt.title("Per-block Identity Along Sequence")
plt.xlabel("Position")
plt.ylabel("Identity")
plt.tight_layout()
plt.show()
