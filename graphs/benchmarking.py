#!/usr/bin/env python3
import os, sys, subprocess, psutil, re, signal, time

import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('seaborn-v0_8')

output_dir = 'results'
os.makedirs(output_dir, exist_ok=True)

threads = 4
runs_per_graph = 5
binary = "max-clique-Release"

graphs = [
    "johnson8-2-4.mtx",       # kept, dropped johnson8-4-4.mtx
    "MANN-a9.mtx",
    "c-fat200-1.mtx",
    "c-fat200-2.mtx",
    "c-fat500-1.mtx",
    "johnson16-2-4.mtx",
    "C125.9.clq",
    "c-fat200-5.mtx",
    "c-fat500-2.mtx",
    "keller4.clq",            # kept, dropped keller4.mtx
    "brock200_2.clq",
    "p_hat300-1.clq",
    "brock200-3.mtx",
    "brock200_4.clq",
    "brock200-1.mtx",
    "gen200-p0-9-44.mtx",
    "hamming8-4.clq",
    "p_hat300-2.clq",
    "c-fat500-5.mtx",
    "C250-9.mtx",
    "p-hat500-1.mtx",
    "p-hat300-3.mtx",
    "c-fat500-10.mtx",
    "brock400-4.mtx",         # kept, dropped brock400-2.mtx
    "p-hat700-1.mtx",
    "p-hat500-2.mtx",
    "MANN-a27.mtx",
    "p-hat500-3.mtx",
    "johnson32-2-4.mtx",
    "inf-luxembourg_osm.mtx",
    "p-hat700-2.mtx",
    "p-hat1000-1.mtx",
    "DSJC500-5.mtx",
    "rgg_n_2_15_s0.mtx",      # kept, dropped rgg_n_2_16_s0.mtx
    "p-hat700-3.mtx",
    "ca-AstroPh.mtx",
    "brock800-4.mtx",         # kept, dropped brock800-2.mtx
    "keller5.mtx",
    "p-hat1000-2.mtx",
    "san1000.mtx",
    "p-hat1500-1.mtx",
    "soc-slashdot.mtx",
    "p-hat1000-3.mtx",
    "hamming10-4.mtx",        # kept, dropped hamming10-2.mtx
    "C1000-9.mtx",
    "MANN-a45.mtx",
    "p-hat1500-2.mtx",
    "rgg_n_2_17_s0.mtx",      # kept, dropped rgg_n_2_18_s0.mtx
    "p-hat1500-3.mtx",
    "inf-belgium_osm.mtx",
    "m14b.mtx",
    "C2000-9.mtx",
    "inf-netherlands_osm.mtx",
    "rgg_n_2_19_s0.mtx",
    "C4000-5.mtx",
    "keller6.mtx",
    "MANN-a81.mtx",
    "rgg_n_2_20_s0.mtx",      # kept, dropped rgg_n_2_22_s0.mtx
    "inf-great-britain_osm.mtx",
    "inf-germany_osm.mtx",
    "inf-asia_osm.mtx",
    "adaptive.mtx",
    "inf-road_central.mtx",
    "packing-500x100x100-b050.mtx",
    "hugetrace-00020.mtx",
    "inf-road_usa.mtx",
    "kron_g500-logn20.mtx",
    "delaunay_n24.mtx",
    "inf-europe_osm.mtx",
    "rgg_n_2_23_s0.mtx",
    "kron_g500-logn21.mtx",
    "rgg_n_2_24_s0.mtx"
]

bin_path = os.path.abspath(binary)
if not os.path.exists(bin_path):
    print(f"FATAL: Binary not found: {bin_path}"); sys.exit(1)
if not os.access(bin_path, os.X_OK):
    print(f"FATAL: Binary is not executable: {bin_path}"); sys.exit(1)

def handle_sigint(signum, frame):
    print("")
    sys.exit(signum)
signal.signal(signal.SIGINT, handle_sigint)

# ----------------------------
# Unified scaled number regex & converters
# ----------------------------
SCALED_NUMBER = r"([\d.]+)\s*(ns|µs|ms|s|B|KB|MB|GB|TB|/s|K/s|M/s|B/s|T/s|K|M|B|T)?"

size_units = {"B":1, "KB":1024, "MB":1024**2, "GB":1024**3, "TB":1024**4}
count_units = {"":1, "K":1_000, "M":1_000_000, "B":1_000_000_000, "T":1_000_000_000_000}
throughput_units = {"/s":1, "K/s":1e3, "M/s":1e6, "B/s":1e9, "T/s":1e12}
time_units = {"ns":1e-9, "µs":1e-6, "ms":1e-3, "s":1}

def convert_scaled(val, unit, unit_map):
    return float(val) * unit_map.get(unit or "", 1)

def scale_int(x):
    if x < 1_000:
        return f'{x}'
    elif x < 1_000_000:
        return f"{x//1_000}K"
    elif x < 1_000_000_000:
        return f"{x//1_000_000}M"
    else:
        return f"{x//1_000_000_000}B"

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
    m = re.search(rf"Graph with\s+\d+\s+Vertices\s+and\s+\d+\s+Edges\s+was\s+Loaded\s+in\s+{SCALED_NUMBER}", text, re.IGNORECASE)
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
    m = re.search(rf"Stream Size\s+{SCALED_NUMBER}\s+and\s+Chunk Size\s+{SCALED_NUMBER}", text, re.IGNORECASE)
    if m:
        sv, su, cv, cu = m.groups()
        return float(convert_scaled(sv, su, size_units)), float(convert_scaled(cv, cu, size_units))
    return 0, 0

def extract_chunks_and_edges_per_chunk(text):
    m = re.search(rf"Estimating\s+{SCALED_NUMBER}\s+Chunks\s+with\s+{SCALED_NUMBER}\s+Edges per Chunk", text, re.IGNORECASE)
    if m:
        cv, cu = m.groups()[0:2]
        ev, eu = m.groups()[2:4]
        return float(convert_scaled(cv, cu, count_units)), float(convert_scaled(ev, eu, count_units))
    return 0, 0

# ----------------------------
# Patterns for throughput
# ----------------------------
READ_THR_RE = r"Reading\s+NUM"
MERGE_VERTICES_THR_RE = r"Merging\s+Vertices\s+NUM"
MERGE_EDGES_THR_RE = r"Merging\s+Edges\s+NUM"
COMPUTE_DEGREES_THR_RE = r"Computing\s+Degrees\s+NUM"
PARSE_THR_RE = r"Parsing\s+Graph\s+NUM"
CONSTR_THR_RE = r"Constructing\s+Edges\s+NUM"

results = []

def run_once(graph_path):
    proc = subprocess.Popen(
        [bin_path, "-t", str(threads), "-i", graph_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    p = psutil.Process(proc.pid)
    proc_start_time = p.create_time()

    cpu_samples = []
    mem_samples = []
    sample_times = []

    while proc.poll() is None:
        try:
            cpu = p.cpu_percent(interval=0.001)
            mem = p.memory_info().rss
            ts = time.time()
            cpu_samples.append(cpu)
            mem_samples.append(mem)
            sample_times.append(ts)
        except psutil.NoSuchProcess:
            break
    elapsed = time.time() - proc_start_time
    stdout, stderr = proc.communicate()
    output = (stdout or "") + (stderr or "")

    if proc.returncode != 0:
        print(f"ERROR: Non-zero exit for {graph_path}")
        print(output)
        return None

    pattern = r'^.*\b(?:WARN|ERROR)\b.*$'
    matches = re.findall(pattern, output, flags=re.MULTILINE)
    warnings = []
    for m in matches:
        warnings = m

    stats = {"graph": graph_path}

    stats["reading"] = extract_time("Finished Reading Graph", output)
    stats["parsing"] = extract_time("Finished Parsing Graph", output)
    stats["total_time"] = extract_total_time(output)
    stats["constructing"] = max(0.0, stats["total_time"] - stats["reading"] - stats["parsing"])

    stats["read_thr"] = extract_first_scaled(READ_THR_RE, output, throughput_units)
    stats["merge_vertices_thr"] = extract_first_scaled(MERGE_VERTICES_THR_RE, output, throughput_units)
    stats["merge_edges_thr"] = extract_first_scaled(MERGE_EDGES_THR_RE, output, throughput_units)
    stats["compute_degrees_thr"] = extract_first_scaled(COMPUTE_DEGREES_THR_RE, output, throughput_units)
    stats["parse_thr"] = extract_first_scaled(PARSE_THR_RE, output, throughput_units)
    stats["construct_thr"] = extract_first_scaled(CONSTR_THR_RE, output, throughput_units)

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

    stats["avg_cpu_percent"] = sum(cpu_samples) / len(cpu_samples) if cpu_samples else 0.0
    stats["peak_mem_bytes"] = max(mem_samples) if mem_samples else 0

    return warnings, stats, elapsed

overall_start = time.time()
for i, graph in enumerate(graphs):
    if not os.path.exists(graph):
        print(f"ERROR: Missing graph file: {graph}"); continue
    best_run, best_time = None, float("inf")
    print(f"running [{i+1:2}/{len(graphs):2}] {graph}", end=" ", flush=True)
    graph_elapsed = 0
    warn_lst = []
    for run_idx in range(1, runs_per_graph + 1):
        print(".", end="", flush=True)
        warnings, stats, elapsed = run_once(graph)
        warn_lst.append(warnings)
        graph_elapsed += elapsed
        if stats and stats["total_time"] < best_time:
            best_time, best_run = stats["total_time"], stats
    print(f" (done in {graph_elapsed:.2f}s)", flush=True)
    for warn in warn_lst:
        print(warn, file=sys.stderr)
    proc = subprocess.run(["osascript", "-e", f'display notification "[{i+1:2}/{len(graphs):2}] {graph} done in {graph_elapsed:.2f}s" with title "Benchmarking"'])
    if best_run:
        results.append(best_run)
overall_elapsed = time.time() - overall_start
if len(results) == 0:
    sys.exit(1)
print(f"completed in {overall_elapsed:.2f} seconds.")

header = [
    "Graph","Vertices","Edges","StreamSize(B)","ChunkSize(B)","EstimatedChunks","EdgesPerChunk",
    "ReadTime(s)","ReadThr(MB/s)","MergeVerticesThr(M/s)","MergeEdgesThr(M/s)","ComputeDegreesThr(M/s)",
    "ParseTime(s)","ParseThr(M/s)","ConstructTime(s)","ConstructThr(M/s)","TotalTime(s)",
    "AvgCPU(%)","PeakMem(B)","AvgCPU(%)/Read","PeakMem(B)/Read","AvgCPU(%)/Parse","PeakMem(B)/Parse","AvgCPU(%)/Construct","PeakMem(B)/Construct",
]

key_map = {
    "Graph": "graph",
    "Vertices": "vertices",
    "Edges": "edges",
    "StreamSize(B)": "stream_size_bytes",
    "ChunkSize(B)": "chunk_size_bytes",
    "EstimatedChunks": "estimated_chunks",
    "EdgesPerChunk": "edges_per_chunk",
    "ReadTime(s)": "reading",
    "ReadThr(MB/s)": "read_thr",
    "MergeVerticesThr(M/s)": "merge_vertices_thr",
    "MergeEdgesThr(M/s)": "merge_edges_thr",
    "ComputeDegreesThr(M/s)": "compute_degrees_thr",
    "ParseTime(s)": "parsing",
    "ParseThr(M/s)": "parse_thr",
    "ConstructTime(s)": "constructing",
    "ConstructThr(M/s)": "construct_thr",
    "TotalTime(s)": "total_time",
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
