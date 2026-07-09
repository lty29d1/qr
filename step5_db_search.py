#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import concurrent.futures as cf
import json
import os
import struct
import zlib
import heapq
from pathlib import Path
from typing import Dict, Any, Iterable, List, Tuple, Optional

DNA_OK = set("ACGTNRYSWKMBDHV")
COMP = str.maketrans("ACGTN", "TGCAN")

# -------------------------
# FASTA + k-mers
# -------------------------
def read_fasta(path: str) -> List[Tuple[str, str]]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"FASTA not found: {p.resolve()}")

    records: List[Tuple[str, str]] = []
    header: str | None = None
    chunks: List[str] = []

    with p.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(chunks).upper()
                    _validate_seq(seq, header)
                    records.append((header, seq))
                header = line[1:].strip()
                chunks = []
            else:
                if header is None:
                    raise ValueError("FASTA format error: sequence line before any header ('>').")
                chunks.append(line.replace(" ", "").replace("\t", ""))

    if header is not None:
        seq = "".join(chunks).upper()
        _validate_seq(seq, header)
        records.append((header, seq))

    if not records:
        raise ValueError("No FASTA records found (empty file?).")
    return records

def _validate_seq(seq: str, header: str) -> None:
    bad = sorted(set(seq) - DNA_OK)
    if bad:
        raise ValueError(f"Invalid characters in FASTA record '{header}': {bad[:20]}")
    if len(seq) == 0:
        raise ValueError(f"Empty sequence for FASTA record '{header}'.")

def iter_kmers(seq: str, k: int) -> Iterable[str]:
    if k <= 0:
        raise ValueError("k must be positive")
    L = len(seq)
    if L < k:
        return
    if "N" not in seq:
        for i in range(0, L - k + 1):
            yield seq[i:i+k]
        return
    for i in range(0, L - k + 1):
        kmer = seq[i:i+k]
        if "N" in kmer:
            continue
        yield kmer

def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]

# -------------------------
# Hashing + MinHash bottom-n
# -------------------------
def splitmix64(x: int) -> int:
    x = (x + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    x ^= (x >> 30) & 0xFFFFFFFFFFFFFFFF
    x = (x * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    x ^= (x >> 27) & 0xFFFFFFFFFFFFFFFF
    x = (x * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    x ^= (x >> 31) & 0xFFFFFFFFFFFFFFFF
    return x & 0xFFFFFFFFFFFFFFFF

def hash64_bytes(b: bytes, seed: int = 0) -> int:
    hi = zlib.crc32(b, seed) & 0xFFFFFFFF
    lo = zlib.crc32(b, seed ^ 0xA5A5A5A5) & 0xFFFFFFFF
    x = ((hi << 32) | lo) & 0xFFFFFFFFFFFFFFFF
    return splitmix64(x ^ (seed & 0xFFFFFFFFFFFFFFFF))

def to_uint(val: int, bits: int) -> int:
    if bits == 64:
        return val & 0xFFFFFFFFFFFFFFFF
    if bits == 32:
        return val & 0xFFFFFFFF
    raise ValueError("bits must be 32 or 64")

def minhash_bottom_n(records: List[Tuple[str, str]], k: int, n: int, bits: int,
                    canonical: bool, seed: int) -> List[int]:
    if n <= 0:
        raise ValueError("n must be positive")
    heap: List[int] = []  # negative values => max-heap
    for _, seq in records:
        for kmer in iter_kmers(seq, k):
            if canonical:
                rc = revcomp(kmer)
                if rc < kmer:
                    kmer = rc
            h = to_uint(hash64_bytes(kmer.encode("ascii"), seed=seed), bits)
            if len(heap) < n:
                heapq.heappush(heap, -h)
            else:
                if h < -heap[0]:
                    heapq.heapreplace(heap, -h)
    return sorted((-x for x in heap))

# -------------------------
# Payload decode (VGQR1.*)
# -------------------------
def unpack_sketch_uint32(b: bytes) -> List[int]:
    if len(b) % 4 != 0:
        raise ValueError("Corrupt sketch bytes (not divisible by 4)")
    n = len(b)//4
    return list(struct.unpack("<" + "I"*n, b))

def decode_payload(payload: str) -> Dict[str, Any]:
    payload = payload.strip()
    if payload.lower().endswith(".txt") and Path(payload).exists():
        payload = Path(payload).read_text(encoding="utf-8").strip()

    if not payload.startswith("VGQR1."):
        raise ValueError("Payload must start with 'VGQR1.'")

    b64 = payload.split(".", 1)[1]
    pad = "=" * ((4 - (len(b64) % 4)) % 4)
    raw = base64.urlsafe_b64decode((b64 + pad).encode("ascii"))

    if len(raw) < 14:
        raise ValueError("Payload too short / corrupt")

    magic, k, n, bits, flags, seed = struct.unpack("<4s H H B B I", raw[:14])
    if magic != b"VGQ1":
        raise ValueError("Bad magic (not a VGQR v1 payload)")
    canonical = bool(flags & 1)

    comp = raw[14:]
    body = zlib.decompress(comp)

    if bits != 32:
        raise ValueError("This script currently supports bits=32 payloads (matches your step4 encoder).")

    sketch = unpack_sketch_uint32(body)
    if len(sketch) != n:
        raise ValueError("Decoded sketch length mismatch")

    return {
        "algo": "minhash-bottom-n",
        "k": int(k),
        "n": int(n),
        "hash_bits": int(bits),
        "canonical": canonical,
        "seed": int(seed),
        "sketch": sorted(int(x) for x in sketch),
    }

# -------------------------
# Similarity (fast)
# -------------------------
def intersection_size_sorted(a: List[int], b: List[int]) -> int:
    i = j = 0
    inter = 0
    la, lb = len(a), len(b)
    while i < la and j < lb:
        av, bv = a[i], b[j]
        if av == bv:
            inter += 1
            i += 1
            j += 1
        elif av < bv:
            i += 1
        else:
            j += 1
    return inter

def minhash_jaccard_est(a: List[int], b: List[int]) -> float:
    # For equal-size bottom-n sketches with same hash function,
    # a common estimator is |A ∩ B| / n (since |A|=|B|=n).
    if len(a) == 0 or len(b) == 0:
        return 0.0
    if len(a) != len(b):
        # still compute overlap / min(n)
        denom = min(len(a), len(b))
    else:
        denom = len(a)
    inter = intersection_size_sorted(a, b)
    return inter / float(denom)

def confidence_label(top_score: float, gap: float) -> str:
    """
    Shared confidence labeling logic for QR search and benchmark scripts.

    Current thresholds are heuristic and should be refined after benchmark testing.
    """
    # Absolute similarity threshold first:
    # prevents unrelated samples from being called confident matches
    if top_score < 0.15:
        return "Out of DB / No confident match"

    # Related to something in the DB, but not a strong exact match
    if top_score < 0.80:
        return "Near DB / Related"

    # Strong similarity, but top two matches are too close
    if gap < 0.05:
        return "In DB / Ambiguous close match"

    # Strong similarity and clear separation from second-best
    return "In DB / Strong Match"

# -------------------------
# DB build + search
# -------------------------
def sketch_from_fasta_file(fasta_path: str, k: int, n: int, bits: int, canonical: bool, seed: int) -> Dict[str, Any]:
    recs = read_fasta(fasta_path)
    # If multi-contig, we sketch across all contigs (concatenated kmers) — fine for now.
    sketch = minhash_bottom_n(recs, k=k, n=n, bits=bits, canonical=canonical, seed=seed)
    return {
        "name": Path(fasta_path).stem,
        "source_fasta": Path(fasta_path).name,
        "algo": "minhash-bottom-n",
        "k": k,
        "n": n,
        "hash_bits": bits,
        "canonical": canonical,
        "seed": seed,
        "sketch": sketch,
    }

def load_db_jsonl(db_path: str) -> List[Dict[str, Any]]:
    items = []
    with Path(db_path).open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            items.append(json.loads(line))
    if not items:
        raise ValueError("DB is empty")
    return items

def check_compat(query: Dict[str, Any], ref: Dict[str, Any]) -> Optional[str]:
    keys = ["k", "n", "hash_bits", "canonical", "seed"]
    for k in keys:
        if int(query[k]) != int(ref[k]):
            return f"Mismatch on {k}: query={query[k]} ref={ref[k]}"
    return None

def main():
    ap = argparse.ArgumentParser(
        description="Step 5 (more complex): build a MinHash reference DB and search it using a query (FASTA, sketch.json, or VGQR payload)."
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build-db", help="Build refs.jsonl from a folder of FASTA files")
    b.add_argument("fasta_dir", help="Folder containing .fa/.fasta/.fna files")
    b.add_argument("--out", default="refs.jsonl", help="Output DB file (default: refs.jsonl)")
    b.add_argument("-k", type=int, default=21)
    b.add_argument("-n", type=int, default=256)
    b.add_argument("--bits", type=int, choices=[32, 64], default=32)
    b.add_argument("--canonical", action="store_true")
    b.add_argument("--seed", type=int, default=1)
    b.add_argument("--workers", type=int, default=0, help="0=auto, 1=sequential, >1=processes")

    s = sub.add_parser("search", help="Search DB with a query")
    s.add_argument("db", help="refs.jsonl produced by build-db")
    s.add_argument("query", help="One of: query.fasta OR sketch.json OR payload.txt OR VGQR1.<...> string")
    s.add_argument("--top", type=int, default=10, help="Top N hits to show (default: 10)")
    s.add_argument("--require-match", action="store_true", help="Error out if k/n/bits/seed/canonical mismatch with DB")
    s.add_argument("--json", dest="json_out", default=None, help="Optional: write ranked results to JSON file")

    # If query is FASTA, we need to know sketch params
    s.add_argument("-k", type=int, default=None, help="(Only for FASTA query) k-mer length")
    s.add_argument("-n", type=int, default=None, help="(Only for FASTA query) sketch size")
    s.add_argument("--bits", type=int, choices=[32, 64], default=None, help="(Only for FASTA query) hash bits")
    s.add_argument("--canonical", action="store_true", help="(Only for FASTA query) canonical k-mers")
    s.add_argument("--seed", type=int, default=None, help="(Only for FASTA query) hash seed")

    args = ap.parse_args()

    if args.cmd == "build-db":
        d = Path(args.fasta_dir)
        if not d.exists() or not d.is_dir():
            raise SystemExit(f"Not a directory: {d}")

        exts = {".fa", ".fasta", ".fna", ".ffn", ".frn"}
        fasta_files = [str(p) for p in sorted(d.iterdir()) if p.suffix.lower() in exts]
        if not fasta_files:
            raise SystemExit("No FASTA files found in directory (expected .fa/.fasta/.fna/...)")

        workers = args.workers
        if workers == 0:
            workers = max(1, (os.cpu_count() or 2) - 1)

        outp = Path(args.out)
        if workers == 1:
            items = [sketch_from_fasta_file(fp, args.k, args.n, args.bits, args.canonical, args.seed) for fp in fasta_files]
        else:
            items = []
            with cf.ProcessPoolExecutor(max_workers=workers) as ex:
                futs = [ex.submit(sketch_from_fasta_file, fp, args.k, args.n, args.bits, args.canonical, args.seed) for fp in fasta_files]
                for fut in cf.as_completed(futs):
                    items.append(fut.result())
            items.sort(key=lambda x: x["name"])

        with outp.open("w", encoding="utf-8") as f:
            for it in items:
                f.write(json.dumps(it) + "\n")

        print(f"Wrote DB: {outp}  (n_refs={len(items)})")
        print(f"Params: k={args.k} n={args.n} bits={args.bits} canonical={args.canonical} seed={args.seed}")
        return

    if args.cmd == "search":
        db_items = load_db_jsonl(args.db)

        # Load query
        qsrc = args.query.strip()
        query_meta: Dict[str, Any]
        if qsrc.startswith("VGQR1.") or (qsrc.lower().endswith(".txt") and Path(qsrc).exists()):
            query_meta = decode_payload(qsrc)
            query_meta["name"] = "query_from_payload"
        elif qsrc.lower().endswith(".json") and Path(qsrc).exists():
            query_meta = json.loads(Path(qsrc).read_text(encoding="utf-8"))
            query_meta["name"] = "query_from_json"
        else:
            # assume FASTA query
            # Infer params from DB unless user overrides
            ref0 = db_items[0]
            k = args.k if args.k is not None else int(ref0["k"])
            n = args.n if args.n is not None else int(ref0["n"])
            bits = args.bits if args.bits is not None else int(ref0["hash_bits"])
            seed = args.seed if args.seed is not None else int(ref0["seed"])
            canonical = True if args.canonical else bool(ref0.get("canonical", False))
            query_meta = sketch_from_fasta_file(qsrc, k=k, n=n, bits=bits, canonical=canonical, seed=seed)
            query_meta["name"] = Path(qsrc).stem

        # Score
        results = []
        for ref in db_items:
            mismatch = check_compat(query_meta, ref)
            if mismatch and args.require_match:
                raise SystemExit("Query/DB parameter mismatch: " + mismatch)
            if mismatch:
                # Still compute, but warn by tagging
                tag = "MISMATCH"
            else:
                tag = ""

            sim = minhash_jaccard_est(query_meta["sketch"], ref["sketch"])
            results.append({
                "ref": ref.get("name", ref.get("source_fasta", "ref")),
                "source_fasta": ref.get("source_fasta", ""),
                "jaccard_est": sim,
                "tag": tag,
            })

        results.sort(key=lambda r: r["jaccard_est"], reverse=True)

        top = max(1, int(args.top))
        print(f"Query: {query_meta.get('name','query')}  k={query_meta['k']} n={query_meta['n']} bits={query_meta['hash_bits']} canonical={query_meta.get('canonical', False)} seed={query_meta.get('seed', 0)}")
        print(f"DB refs: {len(db_items)}  Showing top {top}\n")
        for i, r in enumerate(results[:top], start=1):
            print(f"{i:>2}. jaccard_est={r['jaccard_est']:.4f}  ref={r['ref']}  file={r['source_fasta']} {r['tag']}")

        if args.json_out:
            Path(args.json_out).write_text(json.dumps({
                "query": {k: query_meta.get(k) for k in ["name","k","n","hash_bits","canonical","seed"]},
                "results": results,
            }, indent=2), encoding="utf-8")
            print(f"\nWrote ranked results JSON to: {args.json_out}")

if __name__ == "__main__":
    main()
