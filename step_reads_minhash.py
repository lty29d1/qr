from __future__ import annotations

import argparse
import gzip
import json
import zlib
import heapq
from pathlib import Path

COMP = str.maketrans("ACGT", "TGCA")

def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]

def splitmix64(x: int) -> int:
    x = (x + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    x ^= (x >> 30)
    x = (x * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    x ^= (x >> 27)
    x = (x * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    x ^= (x >> 31)
    return x & 0xFFFFFFFFFFFFFFFF

def hash64_bytes(b: bytes, seed: int) -> int:
    hi = zlib.crc32(b, seed) & 0xFFFFFFFF
    lo = zlib.crc32(b, seed ^ 0xA5A5A5A5) & 0xFFFFFFFF
    x = ((hi << 32) | lo) & 0xFFFFFFFFFFFFFFFF
    return splitmix64(x ^ (seed & 0xFFFFFFFFFFFFFFFF))

def iter_kmers(seq: str, k: int):
    L = len(seq)
    for i in range(L - k + 1):
        kmer = seq[i:i + k]
        # skip kmers with any non-ACGT (handles N etc.)
        if any(c not in "ACGT" for c in kmer):
            continue
        yield kmer

def fastq_iter_seqs_gz(path: str):
    """
    Robust FASTQ parser:
    - header line starts with '@'
    - sequence may span multiple lines until '+' line
    - quality may span multiple lines until it matches seq length
    """
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore", newline="") as f:
        while True:
            header = f.readline()
            if not header:
                return
            header = header.strip()
            if not header:
                continue
            if not header.startswith("@"):
                # try to resync: skip until next header
                continue

            # read sequence lines until '+'
            seq_parts = []
            while True:
                line = f.readline()
                if not line:
                    return
                line = line.rstrip("\r\n")
                if line.startswith("+"):
                    break
                seq_parts.append(line)

            seq = "".join(seq_parts).upper()

            # read quality lines until we have >= len(seq)
            qlen = 0
            while qlen < len(seq):
                qline = f.readline()
                if not qline:
                    return
                qlen += len(qline.rstrip("\r\n"))

            if seq:
                yield seq

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("r1")
    ap.add_argument("r2")
    ap.add_argument("--out", required=True)
    ap.add_argument("-k", type=int, default=21)
    ap.add_argument("-n", type=int, default=256)
    ap.add_argument("--bits", type=int, default=32)
    ap.add_argument("--canonical", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    heap = []

    def push_hash(h: int):
        if len(heap) < args.n:
            heapq.heappush(heap, -h)
        else:
            if h < -heap[0]:
                heapq.heapreplace(heap, -h)

    def process_seq(seq: str):
        for kmer in iter_kmers(seq, args.k):
            if args.canonical:
                rc = revcomp(kmer)
                if rc < kmer:
                    kmer = rc
            h = hash64_bytes(kmer.encode("ascii"), args.seed)
            if args.bits == 32:
                h &= 0xFFFFFFFF
            else:
                h &= 0xFFFFFFFFFFFFFFFF
            push_hash(h)

    for seq in fastq_iter_seqs_gz(args.r1):
        process_seq(seq)
    for seq in fastq_iter_seqs_gz(args.r2):
        process_seq(seq)

    sketch = sorted(set(-x for x in heap))

    meta = {
        "name": f"reads:{Path(args.r1).stem}",
        "algo": "minhash-bottom-n",
        "k": args.k,
        "n": args.n,
        "hash_bits": args.bits,
        "canonical": args.canonical,
        "seed": args.seed,
        "sketch": sketch,
    }

    Path(args.out).write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"Wrote {args.out} with {len(sketch)} hashes")

if __name__ == "__main__":
    main()