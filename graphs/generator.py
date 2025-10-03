#!/usr/bin/env python3
import os
import argparse
import random
import sys
import time
import math
import threading
import queue
from itertools import combinations

BATCH_SIZE = 10000


def edge_index_to_pair(idx, n):
    u = idx // (n - 1) + 1
    v = idx % (n - 1) + 1
    if v >= u:
        v += 1
    return (u, v)


def lcg_params(N, rng):
    if N <= 1:
        return 1, 0
    a = 2 * rng.randrange(1, max(1, N // 2)) + 1
    while math.gcd(a, N) != 1:
        a = 2 * rng.randrange(1, max(1, N // 2)) + 1
    b = rng.randrange(N)
    return a, b


def planted_clique_edges(m):
    clique = list(range(1, m + 1))
    planted = []
    for u, v in combinations(clique, 2):
        planted.append((u, v))
        planted.append((v, u))
    return clique, planted


def build_edges(n, m, num_edges=None, density=None, seed=None):
    max_possible = n * (n - 1)
    if density is not None:
        num_edges = int(round(density * max_possible))
    if num_edges is None:
        num_edges = min(3 * n, max_possible)
    if not (0 <= num_edges <= max_possible):
        raise ValueError("invalid num_edges for the given graph")

    clique, planted = planted_clique_edges(m)
    planted_set = set(planted)
    if len(planted) > num_edges:
        raise ValueError("num_edges too small for planted clique")

    rng = random.Random(seed)
    a, b = lcg_params(max_possible, rng)

    edges = set(planted_set)
    target_remaining = num_edges - len(planted)

    collected = 0
    for i in range(max_possible):
        if collected >= target_remaining:
            break
        j = (a * i + b) % max_possible
        e = edge_index_to_pair(j, n)
        if e in edges:
            continue
        edges.add(e)
        collected += 1

    if len(edges) < num_edges:
        raise RuntimeError("Could not collect enough unique edges")

    remaining_sorted = sorted(edges - planted_set)
    final = planted + remaining_sorted
    assert len(final) == num_edges
    return clique, final


def producer_thread(edges, q):
    batch = []
    for (u, v) in edges:
        batch.append((u, v))
        if len(batch) >= BATCH_SIZE:
            q.put(batch)
            batch = []
    if batch:
        q.put(batch)
    q.put(None)


def consumer_thread(path, header, q, expected_total):
    print("Phase (writing) started")
    out = open(path, "w") if path else sys.stdout
    out.write(header)
    written = 0
    start = time.time()
    next_report = start + 1

    while True:
        batch = q.get()
        if batch is None:
            break
        lines = [f"{u} {v}\n" for u, v in batch]
        out.writelines(lines)
        written += len(batch)

        now = time.time()
        if now >= next_report:
            percent = (written / expected_total) * 100
            elapsed = int(now - start)
            print(f"{elapsed}s: Writing {percent:.1f}%")
            next_report = now + 1

    if path:
        out.close()
    elapsed = int(time.time() - start)
    print(f"Phase (writing) finished in {elapsed}s")

    if written != expected_total:
        raise RuntimeError(
            f"Written edge count mismatch: expected {expected_total}, wrote {written}")


def write_clique_info(clique, output_path=None):
    msg = f"Maximum planted clique (size {len(clique)}): {clique}\n"
    if output_path:
        base, _ = os.path.splitext(output_path)
        with open(base + ".clique", "w") as f:
            f.write(msg)
    else:
        sys.stderr.write(msg)


parser = argparse.ArgumentParser(
    description="Directed graph generator with planted clique (exact counts).")
parser.add_argument("-n", "--num-vertices", type=int, required=True)
parser.add_argument("-m", "--max-clique-size", type=int, required=True)
parser.add_argument("-e", "--num-edges", type=int, default=None)
parser.add_argument("-d", "--density", type=float, default=None)
parser.add_argument("-s", "--seed", type=int, default=None)
parser.add_argument("-o", "--output", type=str, default=None)
args = parser.parse_args()

clique, edges = build_edges(
    n=args.num_vertices,
    m=args.max_clique_size,
    num_edges=args.num_edges,
    density=args.density,
    seed=args.seed
)

total_edges = len(edges)
header = f"{args.num_vertices} {total_edges} {args.max_clique_size}\n"

q = queue.Queue(maxsize=50)
prod = threading.Thread(target=producer_thread, args=(edges, q))
cons = threading.Thread(target=consumer_thread,
                        args=(args.output, header, q, total_edges))
prod.start()
cons.start()
prod.join()
cons.join()

write_clique_info(clique, args.output)
