#!/usr/bin/env python3
from os import path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

label_rot = 90

plt.style.use("seaborn-v0_8-dark")


def scale_int(x):
    def fmt(val, suffix=""):
        d_val = float(f"{val:.2f}")
        if float(d_val).is_integer():
            return f"{int(val)}{suffix}"
        else:
            return f"{d_val}{suffix}"

    suffixes = ["", "K", "M", "B", "T"]
    divisor = 1
    idx = 0

    while idx < len(suffixes) - 1 and x >= divisor * 1000:
        divisor *= 1000
        idx += 1

    return fmt(x / divisor, suffixes[idx])


output_dir = "results"
results_file = f"{output_dir}/benchmark_results.csv"
df = pd.read_csv(results_file)

labels = [
    f"{path.splitext(path.basename(g))[0]} ({scale_int(v)},{scale_int(e)})"
    for g, v, e in zip(df["Graph"], df["Vertices"], df["Edges"])
]

# Phase Throughput Comparison
phase_columns = [
    "ReadThr(MB/s)",
    "ComputeDegreesThr(M/s)",
    "ConstructThr(M/s)",
]
phase_colours = plt.cm.tab10.colors


def colour(i, colours=plt.cm.tab10.colors):
    return colours[i % len(colours)]


fig, ax = plt.subplots(figsize=(14, 6))
read_thr = df["ReadThr(MB/s)"]
comp_deg_thr = df["ComputeDegreesThr(M/s)"]
# ratio = read_thr / np.maximum(parse_thr, comp_deg_thr)

for i, col in enumerate(phase_columns):
    ax.plot(labels, df[col], marker="o", label=col, color=colour(i))

# for i, r in enumerate(ratio):
#     if r < 0.5:
#         ax.plot(labels[i], read_thr[i], marker='v', markersize=8)
# ax.annotate(df['Graph'][i], ("", read_thr[i]), textcoords="offset points", xytext=(0,10), ha='center', fontsize=6)

ax.set_ylabel("Throughput")
ax.set_title("Phase Throughput Comparison Across Graphs (Anomalies Marked)")
ax.legend(loc="best")
plt.setp(ax.get_xticklabels(), rotation=label_rot, ha="right")
plt.tight_layout()
plt.savefig(f"{output_dir}/phase_throughput_comparison.png")
print(f"wrote {output_dir}/phase_throughput_comparison.png")
plt.close()

# Normalised Times
phase_times = []
phase_labels = []

if "ReadTime(ns)" in df:
    phase_times.append(df["ReadTime(ns)"])
    phase_labels.append("Read")
if "ParseTime(ns)" in df:
    phase_times.append(df["ParseTime(ns)"])
    phase_labels.append("Parse")
if "ConstructTime(ns)" in df:
    phase_times.append(df["ConstructTime(ns)"])
    phase_labels.append("Construct")
if "MergeVerticesTime(ns)" in df:
    phase_times.append(df["MergeVerticesTime(ns)"])
    phase_labels.append("MergeVertices")
if "MergeEdgesTime(ns)" in df:
    phase_times.append(df["MergeEdgesTime(ns)"])
    phase_labels.append("MergeEdges")
if "ComputeDegreesTime(ns)" in df:
    phase_times.append(df["ComputeDegreesTime(ns)"])
    phase_labels.append("ComputeDegrees")

phase_pcts = [(p / df["TotalTime(ns)"]).fillna(0) * 100 for p in phase_times]

fig, ax = plt.subplots(figsize=(12, 6))
bottom = None

for pct, label in zip(phase_pcts, phase_labels):
    bars = ax.bar(labels, pct, bottom=bottom, label=f"{label} (%)")
    # if label == 'Read':
    #     for i, bar in enumerate(bars):
    #         if df['ReadTime(ns)'][i] > 2 * df['ParseTime(ns)'][i]:
    #             bar.set_edgecolor('black')
    #             bar.set_linewidth(1.5)
    bottom = pct if bottom is None else bottom + pct

ax.set_ylabel("Percentage of Total Time (%)")
ax.set_title("Normalised Benchmarking Times per Graph (Anomalies Outlined)")
ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
plt.setp(ax.get_xticklabels(), rotation=label_rot, ha="right")
plt.tight_layout()
plt.savefig(f"{output_dir}/benchmark_times_normalised.png")
print(f"wrote {output_dir}/benchmark_times_normalised.png")
plt.close()

# Total Time vs Number of Edges
df_sorted = df.sort_values(by="Edges")
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(
    df_sorted["Edges"],
    df_sorted["TotalTime(ns)"],
    c="darkblue",
    marker="o",
    label="Edges",
)
ax.plot(df_sorted["Edges"], df_sorted["TotalTime(ns)"], c="darkblue", linewidth=1)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Edges (log scale)")
ax.set_ylabel("Total Time (s, log scale)")
ax.set_title("Total Time vs Number of Edges")
ax.legend(loc="best")
plt.tight_layout()
plt.savefig(f"{output_dir}/total_time_vs_edges.png")
print(f"wrote {output_dir}/total_time_vs_edges.png")
plt.close()

# Stream Size vs Estimated Chunks and Edges per Chunk
df_sorted = df.sort_values(by="StreamSize(B)")
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(
    df_sorted["StreamSize(B)"],
    df_sorted["EstimatedChunks"],
    c="steelblue",
    marker="o",
    label="Estimated Chunks",
)
ax.plot(
    df_sorted["StreamSize(B)"], df_sorted["EstimatedChunks"], c="steelblue", linewidth=1
)
ax.scatter(
    df_sorted["StreamSize(B)"],
    df_sorted["EdgesPerChunk"],
    c="darkorange",
    marker="s",
    label="Edges per Chunk",
)
ax.plot(
    df_sorted["StreamSize(B)"], df_sorted["EdgesPerChunk"], c="darkorange", linewidth=1
)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Stream Size (B, log scale)")
ax.set_ylabel("Count (log scale)")
ax.set_title("Stream Size vs Estimated Chunks and Edges per Chunk")
ax.legend(loc="best")
plt.tight_layout()
plt.savefig(f"{output_dir}/streamsize_vs_chunks_edges.png")
print(f"wrote {output_dir}/streamsize_vs_chunks_edges.png")
plt.close()

# CPU and Memory Usage per Graph
peak_mem_mb = df["PeakMem(B)"] / (1024 * 1024)
fig, ax1 = plt.subplots(figsize=(14, 6))

colour_cpu = "tab:red"
ax1.set_ylabel("Average CPU (%)", color=colour_cpu)
cpu_bars = ax1.bar(
    labels, df["AvgCPU(%)"], color=colour_cpu, alpha=0.6, label="Avg CPU (%)"
)
ax1.tick_params(axis="y", labelcolor=colour_cpu)

colour_mem = "tab:blue"
ax2 = ax1.twinx()
ax2.set_ylabel("Peak Memory (MB)", color=colour_mem)
(mem_line,) = ax2.plot(
    labels, peak_mem_mb, color=colour_mem, marker="o", label="Peak Memory (MB)"
)
ax2.set_yscale("log")
ax2.tick_params(axis="y", labelcolor=colour_mem)

plt.setp(ax1.get_xticklabels(), rotation=label_rot, ha="right")
fig.legend(
    [cpu_bars, mem_line],
    ["Avg CPU (%)", "Peak Memory (MB)"],
    loc="center left",
    bbox_to_anchor=(1.25, 0.5),
)
plt.title("CPU and Memory Usage per Graph")
plt.tight_layout()
plt.savefig(f"{output_dir}/cpu_mem_usage_per_graph.png")
print(f"wrote {output_dir}/cpu_mem_usage_per_graph.png")
plt.close()

fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

# CPU usage
for phase in ["Read", "Parse", "Construct"]:
    axes[0].plot(labels, df[f"AvgCPU(%)/{phase}"], marker="o", label=phase)
axes[0].set_ylabel("Average CPU (%)")
axes[0].set_title("CPU Usage per Phase")
axes[0].legend()

# Memory usage
for phase in ["Read", "Parse", "Construct"]:
    axes[1].plot(
        labels, df[f"PeakMem(B)/{phase}"] / (1024 * 1024), marker="o", label=phase
    )
axes[1].set_ylabel("Peak Memory (MB)")
axes[1].set_title("Memory Usage per Phase")
axes[1].legend()

plt.xticks(rotation=label_rot, ha="right")
plt.tight_layout()
plt.savefig(f"{output_dir}/phase_resource_usage.png")
print(f"wrote {output_dir}/phase_resource_usage.png")
plt.close()

# Heatmaps
vertex_bins = [0, 1e3, 1e4, 1e5, 1e6, 1e7, np.inf]
vertex_labels = ["<1K", "1K–10K", "10K–100K", "100K–1M", "1M–10M", ">10M"]

edge_bins = [0, 1e4, 1e5, 1e6, 1e7, 1e8, np.inf]
edge_labels = ["<10K", "10K–100K", "100K–1M", "1M–10M", "10M–100M", ">100M"]

df["VertexBin"] = pd.cut(df["Vertices"], bins=vertex_bins, labels=vertex_labels)
df["EdgeBin"] = pd.cut(df["Edges"], bins=edge_bins, labels=edge_labels)


def create_heatmap(metric, title, output_file, cm, compare_metric=None):
    pivot_table = df.pivot_table(
        index="EdgeBin",
        columns="VertexBin",
        values=metric,
        aggfunc="mean",
        observed=False,
    ).reindex(index=edge_labels, columns=vertex_labels)
    fig, ax = plt.subplots(figsize=(10, 6))
    cax = ax.imshow(pivot_table, cmap=cm, aspect="auto")

    for i in range(pivot_table.shape[0]):
        for j in range(pivot_table.shape[1]):
            value = pivot_table.iloc[i, j]
            if pd.notna(value):
                ax.text(
                    j,
                    i,
                    f"{scale_int(value)}",
                    ha="center",
                    va="center",
                    color="black",
                    fontsize=15,
                ).set_path_effects(
                    [
                        path_effects.Stroke(linewidth=2, foreground="white"),
                        path_effects.Normal(),
                    ]
                )

    if compare_metric:
        comp_table = df.pivot_table(
            index="EdgeBin",
            columns="VertexBin",
            values=compare_metric,
            aggfunc="mean",
            observed=False,
        ).reindex(index=edge_labels, columns=vertex_labels)
        residuals = pivot_table - comp_table
        threshold = np.nanstd(residuals) * 2
        for i in range(pivot_table.shape[0]):
            for j in range(pivot_table.shape[1]):
                if (
                    pd.notna(residuals.iloc[i, j])
                    and abs(residuals.iloc[i, j]) > threshold
                ):
                    ax.add_patch(
                        plt.Rectangle(
                            (j - 0.5, i - 0.5),
                            0.99,
                            0.99,
                            fill=False,
                            edgecolor="black",
                            linewidth=1.5,
                        )
                    )

    ax.set_xticks(np.arange(len(vertex_labels)))
    ax.set_xticklabels(vertex_labels)
    ax.set_yticks(np.arange(len(edge_labels)))
    ax.set_yticklabels(edge_labels)
    ax.set_title(title)
    ax.set_xlabel("Vertex Count Bin")
    ax.set_ylabel("Edge Count Bin")
    fig.colorbar(cax, ax=ax, label=f"Average {metric}")
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"wrote {output_file}")
    plt.close()


create_heatmap(
    "ReadThr(MB/s)",
    "Average Read Throughput (MB/s)",
    f"{output_dir}/heatmap_read_thr.png",
    plt.cm.RdYlGn,
)
create_heatmap(
    "ReadTime(ns)",
    "Average Read Time (ns)",
    f"{output_dir}/heatmap_read_time.png",
    plt.cm.RdYlGn_r,
    compare_metric="ParseTime(ns)",
)
create_heatmap(
    "EdgesPerChunk",
    "Average Edges Per Chunk",
    f"{output_dir}/heatmap_edges_per_chunk.png",
    plt.cm.RdYlGn,
)
