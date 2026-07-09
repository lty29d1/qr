import csv
from pathlib import Path
import step5_db_search as db

DB_PATH = "viral9_refs.jsonl"

# Put out-of-database FASTA files here later:
# examples: non_sars_coronavirus.fasta, hepC.fasta, unrelated_virus.fasta
QUERY_DIR = Path("ood_queries")

OUT_CSV = "ood_confidence_results.csv"




def main():
    refs = db.load_db_jsonl(DB_PATH)

    if not QUERY_DIR.exists():
        QUERY_DIR.mkdir()
        print(f"Created folder: {QUERY_DIR}")
        print("Put out-of-database FASTA files in this folder, then rerun.")
        return

    fasta_files = sorted(
        list(QUERY_DIR.glob("*.fasta")) +
        list(QUERY_DIR.glob("*.fa")) +
        list(QUERY_DIR.glob("*.fna"))
    )

    if not fasta_files:
        print(f"No FASTA files found in {QUERY_DIR}")
        return

    rows = []

    ref0 = refs[0]

    for fasta in fasta_files:
        if fasta.name.startswith("._"):
            continue

        query = db.sketch_from_fasta_file(
            str(fasta),
            k=int(ref0["k"]),
            n=int(ref0["n"]),
            bits=int(ref0["hash_bits"]),
            canonical=bool(ref0["canonical"]),
            seed=int(ref0["seed"]),
        )

        results = []
        for ref in refs:
            sim = db.minhash_jaccard_est(query["sketch"], ref["sketch"])
            results.append((ref["name"], sim))

        results.sort(key=lambda x: x[1], reverse=True)

        top_name, top_score = results[0]
        second_name, second_score = results[1]
        gap = top_score - second_score
        label = db.confidence_label(top_score, gap)

        rows.append([
            fasta.stem,
            top_name,
            f"{top_score:.4f}",
            second_name,
            f"{second_score:.4f}",
            f"{gap:.4f}",
            label
        ])

    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Query",
            "Top match",
            "Top similarity",
            "Second match",
            "Second similarity",
            "Gap",
            "Confidence category"
        ])
        writer.writerows(rows)

    print(f"Wrote {OUT_CSV}")
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
