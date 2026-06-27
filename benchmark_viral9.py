import csv
from pathlib import Path
import step5_db_search as db

DB_PATH = "viral9_refs.jsonl"
QUERY_DIR = Path("viral9_benchmark_refs")
OUT_CSV = "viral9_benchmark_results.csv"

refs = db.load_db_jsonl(DB_PATH)

rows = []

for fasta in sorted(QUERY_DIR.glob("*.fasta")):
    if fasta.name.startswith("._"):
        continue

    ref0 = refs[0]
    q = db.sketch_from_fasta_file(
        str(fasta),
        k=int(ref0["k"]),
        n=int(ref0["n"]),
        bits=int(ref0["hash_bits"]),
        canonical=bool(ref0["canonical"]),
        seed=int(ref0["seed"]),
    )

    results = []
    for r in refs:
        sim = db.minhash_jaccard_est(q["sketch"], r["sketch"])
        results.append((r["name"], sim))

    results.sort(key=lambda x: x[1], reverse=True)

    top_name, top_score = results[0]
    second_name, second_score = results[1]
    gap = top_score - second_score

    if gap >= 0.20:
        confidence = "High"
    elif gap >= 0.05:
        confidence = "Medium"
    else:
        confidence = "Low"

    query_name = fasta.stem
    correct = "YES" if query_name == top_name else "NO"

    rows.append([
        query_name,
        top_name,
        f"{top_score:.4f}",
        second_name,
        f"{second_score:.4f}",
        f"{gap:.4f}",
        confidence,
        correct,
    ])

with open(OUT_CSV, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "Query genome",
        "Top match",
        "Top similarity",
        "Second match",
        "Second similarity",
        "Confidence gap",
        "Confidence",
        "Correct?"
    ])
    writer.writerows(rows)

print(f"Wrote {OUT_CSV}")
for row in rows:
    print(row)