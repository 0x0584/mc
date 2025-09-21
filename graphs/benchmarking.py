#!/usr/bin/env python3
import os, sys, subprocess, psutil, re, signal, time

import pandas as pd
import matplotlib.pyplot as plt

plt.style.use("seaborn-v0_8")

output_dir = "results"
os.makedirs(output_dir, exist_ok=True)

threads = 4
runs_per_graph = 5
binary = "max-clique"
graph_groups = {
    "small": {
        # <~1K vertices, <~20K edges — fits in cache, dominated by fixed overheads
        "low density": [
            "johnson8-2-4.mtx",  # v=28, e=210. Very small, symmetric; correctness checks.
            "c-fat200-1.mtx",  # v=200, e=1534. Sparse planted clique.
            "c-fat200-2.mtx",  # v=200, e=3235. Slightly denser planted clique.
            "c-fat500-1.mtx",  # v=500, e=4459. Sparse, larger vertex set.
            "johnson16-2-4.mtx",  # v=120, e=5460. Symmetric, validation target.
            "c-fat500-2.mtx",  # v=500, e=9139. Medium sparse planted clique.
            "c-fat500-5.mtx",  # v=500, e=23191. Denser planted clique.
        ],
        "medium-high density": [
            "c125.9.clq",  # v=125, e=6963. Medium density, large planted clique.
            "keller4.clq",  # v=171, e=9435. Structured, challenging.
            "brock200_2.clq",  # v=200, e=9876. Many large near-cliques.
            "brock200-3.mtx",  # v=200, e=12048. Dense, hard instance.
            "brock200-1.mtx",  # v=200, e=14834. Densest brock200 variant.
            "brock200-4.clq",  # v=200, e=13089. Another brock200 variant.
            "gen200-p0-9-44.mtx",  # v=200, e=17910. Random dense.
            "p_hat300-1.clq",  # v=300, e=10933. Sparse random.
            "p_hat300-2.clq",  # v=300, e=21928. Denser random.
        ],
    },
    "medium": {
        # ~1K–10K vertices, 10K–500K edges — parse cost visible, memory still modest
        "low density": [
            "p-hat1000-1.mtx",  # v=1000, e=122253. Large sparse random.
            "p-hat1500-1.mtx",  # v=1500, e=284923. Sparse random.
            "p-hat500-1.mtx",  # v=500, e=31569. Sparse random.
            "p-hat700-1.mtx",  # v=700, e=60999. Sparse random.
        ],
        "medium-high density": [
            "brock400-2.mtx",  # v=400, e=59786. Larger brock, hard.
            "brock400-4.mtx",  # v=400, e=59765. Variant of brock400-2.
            "brock800-4.mtx",  # v=800, e=207643. Very large brock.
            "brock800-2.mtx",  # v=800, e=208166. Variant of brock800-4.
            "p-hat1000-3.mtx",  # v=1000, e=371746. Large dense random.
            "p-hat1500-3.mtx",  # v=1500, e=1681696. Dense random.
            "p-hat500-3.mtx",  # v=500, e=121925. Dense random.
            "p-hat700-3.mtx",  # v=700, e=244799. Dense random.
            "hamming10-4.mtx",  # v=1024, e=434176. Large structured Hamming.
            "hamming10-2.mtx",  # v=1024, e=5120. Sparse structured Hamming.
            "MANN-a9.mtx",  # v=45, e=918. Small structured MANN.
            "MANN-a27.mtx",  # v=378, e=7020. Structured MANN.
            "MANN-a45.mtx",  # v=1035, e=19890. Structured MANN.
            "MANN-a81.mtx",  # v=3321, e=63720. Structured MANN.
            "adaptive.mtx",  # Synthetic adaptive graph.
            "m14b.mtx",  # Synthetic benchmark graph.
            "packing-500x100x100-b050.mtx",  # Packing problem graph.
            "san1000.mtx",  # v=1000, e=250500. Dense synthetic.
        ],
    },
    "large": {
        # 10K–1M vertices, 200K–5M edges — I/O and parse both significant
        "dense and structured": [
            "keller5.mtx",  # v=776, e=225990. Dense structured.
            "keller6.mtx",  # v=3361, e=4619898. Dense, structured, high complexity.
            "C1000-9.mtx",  # v=1000, e=450450. Dense planted clique.
            "C2000-9.mtx",  # v=2000, e=1799532. Large planted clique.
            "C4000-5.mtx",  # v=4000, e=4000268. Massive, small clique.
            "kron_g500-logn20.mtx",  # Kronecker synthetic, skewed degree.
            "kron_g500-logn21.mtx",  # Larger Kronecker synthetic.
        ],
        "real-world and geometric": [
            "ca-AstroPh.mtx",  # v=17903, e=196972. Collaboration network.
            "soc-slashdot.mtx",  # v=70068, e=358647. Social network, skewed degree.
            "hugetrace-00020.mtx",  # Large trace graph.
            "inf-asia_osm.mtx",  # Road network, planar.
            "inf-belgium_osm.mtx",  # Road network, planar.
            "inf-europe_osm.mtx",  # Road network, planar.
            "inf-germany_osm.mtx",  # Road network, planar.
            "inf-great-britain_osm.mtx",  # Road network, planar.
            "inf-luxembourg_osm.mtx",  # Road network, planar.
            "inf-netherlands_osm.mtx",  # Road network, planar.
            "inf-road_central.mtx",  # Road network, planar.
        ],
        "synthetic geometric": [
            "rgg_n_2_15_s0.mtx",  # Random geometric, small.
            "rgg_n_2_16_s0.mtx",
            "rgg_n_2_17_s0.mtx",
            "rgg_n_2_18_s0.mtx",
            "rgg_n_2_19_s0.mtx",
            "rgg_n_2_20_s0.mtx",
            "rgg_n_2_22_s0.mtx",
            "rgg_n_2_23_s0.mtx",
        ],
    },
    "huge": {
        # >10M vertices, tens to hundreds of millions of edges — I/O-bound, memory locality critical
        "extreme scale": [
            "delaunay_n24.mtx",  # v=16777216, e=50331601. Massive geometric mesh.
            "inf-road_usa.mtx",  # v=23947347, e=28854312. Sparse road network.
            "rgg_n_2_24_s0.mtx",  # v=16777216, e=132557200. Largest random geometric.
        ],
    },
}

selected_groups = ["small", "medium", "large", "huge"]

graphs = []
for group in selected_groups:
    for subgroup in graph_groups[group]:
        graphs.extend(graph_groups[group][subgroup])

graph_name_padding = max(len(graph) for graph in graphs) + 2

bin_path = os.path.abspath(binary)
if not os.path.exists(bin_path):
    print(f"FATAL: Binary not found: {bin_path}")
    sys.exit(1)
if not os.access(bin_path, os.X_OK):
    print(f"FATAL: Binary is not executable: {bin_path}")
    sys.exit(1)


def handle_sigint(signum, frame):
    print("")
    sys.exit(signum)


signal.signal(signal.SIGINT, handle_sigint)

# ----------------------------
# Unified scaled number regex & converters
# ----------------------------
SCALED_NUMBER = r"([\d.]+)\s*(ns|µs|ms|s|hour|minute|minutes|hours|day|days|B|KiB|MiB|GiB|TiB|/s|K/s|M/s|B/s|T/s|K|M|B|T)?"

size_units = {"B": 1, "KiB": 1024, "MiB": 1024**2, "GiB": 1024**3, "TiB": 1024**4}
count_units = {
    "": 1,
    "K": 1_000,
    "M": 1_000_000,
    "B": 1_000_000_000,
    "T": 1_000_000_000_000,
}
throughput_units = {"/s": 1, "K/s": 1e3, "M/s": 1e6, "B/s": 1e9, "T/s": 1e12}
time_units = {
    "ns": 1,
    "µs": 1_000,
    "ms": 1_000_000,
    "s": 1_000_000_000,
    "minutes": 600_00_000_000,
    "hours": 3_600_000_000_000,
    "days": 86_400_000_000_000,
    "minute": 600_00_000_000,
    "hour": 3_600_000_000_000,
    "day": 86_400_000_000_000,
}


def convert_scaled(val, unit, unit_map):
    return float(val) * unit_map.get(unit or "", 1)


def format_value(x):
    x_str = f"{x:.3f}"
    y = float(x_str)
    return f"{int(x)}" if y.is_integer() else x_str


def scale_value(x):
    if x < 1_000:
        return f"{format_value(x)}"
    elif x < 1_000_000:
        return f"{format_value(x/1_000)}K"
    elif x < 1_000_000_000:
        return f"{format_value(x/1_000_000)}M"
    else:
        return f"{format_value(x/1_000_000_000)}B"


time_scale_units = [
    (86400000000000, "day"),
    (3600000000000, "hour"),
    (60000000000, "minute"),
    (1000000000, "s"),
    (1000000, "ms"),
    (1000, "μs"),
    (1, "ns"),
]


def format_duration(stamp):
    for unit_ns, unit_str in time_scale_units:
        if stamp >= unit_ns:
            value = stamp / unit_ns
            plural = (
                "s" if value >= 2.0 and unit_str in ["minute", "hour", "day"] else ""
            )
            return f"{format_value(value)} {unit_str}{plural}"
    return f"{format_value(stamp)} ns"


# ----------------------------
# Parsing helpers
# ----------------------------
def extract_time(label, text):
    m = re.search(rf"{label}.*in\s+{SCALED_NUMBER}", text, re.IGNORECASE)
    if m:
        val, unit = m.groups()
        return convert_scaled(val, unit, time_units)
    return 0.0


def extract_total_time(text):
    m = re.search(
        rf"Graph with\s+\d+\s+Vertices\s+and\s+\d+\s+Edges\s+was\s+Loaded\s+in\s+{SCALED_NUMBER}",
        text,
        re.IGNORECASE,
    )
    if m:
        val, unit = m.groups()
        return convert_scaled(val, unit, time_units)
    return 0.0


def extract_vertices_edges(text):
    m = re.search(r"with\s+(\d+)\s+vertices\s+and\s+(\d+)\s+edges", text, re.IGNORECASE)
    return (float(m.group(1)), float(m.group(2))) if m else (0, 0)


def extract_first_scaled(pattern, text, unit_map):
    m = re.search(pattern.replace("NUM", SCALED_NUMBER), text, re.IGNORECASE)
    if m:
        val, unit = m.groups()
        return convert_scaled(val, unit, unit_map)
    return 0.0


def extract_stream_and_chunk(text):
    m = re.search(
        rf"Stream Size\s+{SCALED_NUMBER}\s+and\s+Chunk Size\s+{SCALED_NUMBER}",
        text,
        re.IGNORECASE,
    )
    if m:
        sv, su, cv, cu = m.groups()
        return float(convert_scaled(sv, su, size_units)), float(
            convert_scaled(cv, cu, size_units)
        )
    return 0, 0


def extract_chunks_and_edges_per_chunk(text):
    m = re.search(
        rf"Estimating\s+{SCALED_NUMBER}\s+Chunks\s+with\s+{SCALED_NUMBER}\s+Edges per Chunk",
        text,
        re.IGNORECASE,
    )
    if m:
        cv, cu = m.groups()[0:2]
        ev, eu = m.groups()[2:4]
        return (
            float(convert_scaled(cv, cu, count_units)),
            float(convert_scaled(ev, eu, count_units)),
        )
    return 0, 0


# ----------------------------
# Patterns for throughput
# ----------------------------

READ_THR_RE = rf"Reading\s+{SCALED_NUMBER}"
MERGE_VERTICES_THR_RE = rf"Merging\s+Vertices\s+{SCALED_NUMBER}"
MERGE_EDGES_THR_RE = rf"Merging\s+Edges\s+{SCALED_NUMBER}"
COMPUTE_DEGREES_THR_RE = rf"Computing\s+Degrees\s+{SCALED_NUMBER}"
CONSTR_THR_RE = rf"Constructing\s+Edges\s+{SCALED_NUMBER}"


def run_once(graph_path):
    proc = subprocess.Popen(
        [bin_path, "-t", str(threads), "-i", graph_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    p = psutil.Process(proc.pid)
    proc_start_time = time.time_ns()

    cpu_samples = []
    mem_samples = []
    sample_times = []

    while proc.poll() is None:
        try:
            cpu = p.cpu_percent(interval=0.001)
            mem = p.memory_info().rss
            ts = time.time_ns()
            cpu_samples.append(cpu)
            mem_samples.append(mem)
            sample_times.append(ts)
        except psutil.NoSuchProcess:
            break

    elapsed = time.time_ns() - proc_start_time
    stdout, stderr = proc.communicate()
    output = (stdout or "") + (stderr or "")

    if proc.returncode != 0:
        print(f"ERROR: {proc.returncode} exit for {graph_path}")
        print(output)
        return None

    pattern = r"^.*\b(?:WARN|ERROR)\b.*$"
    matches = re.findall(pattern, output, flags=re.MULTILINE)
    warnings = []
    for m in matches:
        warnings = m

    stats = {"graph": graph_path}

    stats["reading"] = extract_time("Finished Reading Graph", output)
    stats["parsing"] = extract_time("Finished Parsing Graph", output)
    stats["total_time"] = extract_total_time(output)
    stats["constructing"] = max(
        0.0, stats["total_time"] - stats["reading"] - stats["parsing"]
    )

    stats["read_thr"] = extract_first_scaled(READ_THR_RE, output, throughput_units)
    stats["merge_vertices_thr"] = extract_first_scaled(
        MERGE_VERTICES_THR_RE, output, throughput_units
    )
    stats["merge_edges_thr"] = extract_first_scaled(
        MERGE_EDGES_THR_RE, output, throughput_units
    )
    stats["compute_degrees_thr"] = extract_first_scaled(
        COMPUTE_DEGREES_THR_RE, output, throughput_units
    )
    stats["construct_thr"] = extract_first_scaled(
        CONSTR_THR_RE, output, throughput_units
    )

    stream_size, chunk_size = extract_stream_and_chunk(output)
    chunks, edges_per_chunk = extract_chunks_and_edges_per_chunk(output)
    stats["stream_size_bytes"] = stream_size
    stats["chunk_size_bytes"] = chunk_size
    stats["estimated_chunks"] = chunks
    stats["edges_per_chunk"] = edges_per_chunk

    v, e = extract_vertices_edges(output)
    stats["vertices"], stats["edges"] = v, e

    read_time = stats["reading"]
    parse_time = stats["parsing"]
    construct_time = stats["constructing"]
    total_time = stats["total_time"]

    read_end = proc_start_time + read_time
    parse_end = read_end + parse_time
    construct_end = parse_end + construct_time

    def slice_phase(start, end):
        idx = [i for i, t in enumerate(sample_times) if start <= t < end]
        if not idx:
            return 0.0, 0
        avg_cpu = sum(cpu_samples[i] for i in idx) / len(idx)
        peak_mem = max(mem_samples[i] for i in idx)
        return avg_cpu, peak_mem

    cpu_read, mem_read = slice_phase(proc_start_time, read_end)
    cpu_parse, mem_parse = slice_phase(read_end, parse_end)
    cpu_construct, mem_construct = slice_phase(parse_end, construct_end)

    stats["cpu_read"] = cpu_read
    stats["mem_read"] = mem_read
    stats["cpu_parse"] = cpu_parse
    stats["mem_parse"] = mem_parse
    stats["cpu_construct"] = cpu_construct
    stats["mem_construct"] = mem_construct

    stats["avg_cpu_percent"] = (
        sum(cpu_samples) / len(cpu_samples) if cpu_samples else 0.0
    )
    stats["peak_mem_bytes"] = max(mem_samples) if mem_samples else 0

    return warnings, stats, elapsed


results = []
worst_results = []  # TODO: output worst

overall_start = time.time_ns()
for i, graph in enumerate(graphs):
    if not os.path.exists(graph):
        print(f"ERROR: Missing graph file: {graph}")
        continue
    best_run, best_time = None, float("inf")
    worst_run, worst_time = None, float("-inf")
    graph_elapsed = 0
    warn_lst = []

    # print(f"running [{i+1:2}/{len(graphs):2}] {graph:<{graph_name_padding}}", end=" ", flush=True)
    print(f"running [{i+1:2}/{len(graphs):2}] {graph}", end=" ", flush=True)
    for run_idx in range(1, runs_per_graph + 1):
        print(".", end="", flush=True)
        warnings, stats, elapsed = run_once(graph)
        if len(warnings):
            warn_lst.append(warnings)
        graph_elapsed += elapsed
        if best_time is None or elapsed < best_time:
            best_time, best_run = elapsed, stats
        if worst_time is None or elapsed > worst_time:
            worst_time, worst_run = elapsed, stats
    print(
        f" (done in \033[1m{format_duration(graph_elapsed)}\033[0m"
        f" worst \033[35m{format_duration(worst_time)}\033[0m"
        f" best \033[36m{format_duration(best_time)}\033[0m)",
        flush=True,
    )
    for warn in warn_lst:
        print(warn, file=sys.stderr)
    if best_run:
        results.append(best_run)
    proc = subprocess.run(
        [
            "osascript",
            "-e",
            f'display notification "[{i+1:2}/{len(graphs):2}] {graph} done in {format_duration(graph_elapsed)}" with title "Benchmarking"',
        ]
    )
overall_elapsed = time.time_ns() - overall_start

if len(results) == 0:
    sys.exit(1)
print(f"completed in {format_duration(overall_elapsed)}.")

header = [
    # Graph information
    "Graph",
    "Vertices",
    "Edges",
    # File size information
    "StreamSize(B)",
    "ChunkSize(B)",
    # Estimates on reading
    "EstimatedChunks",
    "EdgesPerChunk",
    # Telemetry of Phases
    "ReadTime(ns)",
    "ReadThr(MB/s)",
    "MergeVerticesThr(M/s)",
    "MergeEdgesThr(M/s)",
    "ComputeDegreesThr(M/s)",
    "ParseTime(ns)",
    "ConstructTime(ns)",
    "ConstructThr(M/s)",
    "TotalTime(ns)",
    # Telemetry of resources
    "AvgCPU(%)",
    "PeakMem(B)",
    "AvgCPU(%)/Read",
    "PeakMem(B)/Read",
    "AvgCPU(%)/Parse",
    "PeakMem(B)/Parse",
    "AvgCPU(%)/Construct",
    "PeakMem(B)/Construct",
]

key_map = {
    # Graph information
    "Graph": "graph",
    "Vertices": "vertices",
    "Edges": "edges",
    # File size information
    "StreamSize(B)": "stream_size_bytes",
    "ChunkSize(B)": "chunk_size_bytes",
    # Estimates on reading
    "EstimatedChunks": "estimated_chunks",
    "EdgesPerChunk": "edges_per_chunk",
    # Telemetry of Phases
    "ReadTime(ns)": "reading",
    "ReadThr(MB/s)": "read_thr",
    "MergeVerticesThr(M/s)": "merge_vertices_thr",
    "MergeEdgesThr(M/s)": "merge_edges_thr",
    "ComputeDegreesThr(M/s)": "compute_degrees_thr",
    "ParseTime(ns)": "parsing",
    "ConstructTime(ns)": "constructing",
    "ConstructThr(M/s)": "construct_thr",
    "TotalTime(ns)": "total_time",
    # Telemetry of resources
    "AvgCPU(%)": "avg_cpu_percent",
    "PeakMem(B)": "peak_mem_bytes",
    "AvgCPU(%)/Read": "cpu_read",
    "PeakMem(B)/Read": "mem_read",
    "AvgCPU(%)/Parse": "cpu_parse",
    "PeakMem(B)/Parse": "mem_parse",
    "AvgCPU(%)/Construct": "cpu_construct",
    "PeakMem(B)/Construct": "mem_construct",
}

rows = [[r.get(key_map[col], 0) for col in header] for r in results]
df = pd.DataFrame(rows, columns=header)
df.to_csv(f"{output_dir}/benchmark_results.csv", index=False, float_format="%.6f")
print(f"wrote {output_dir}/benchmark_results.csv")
