**# QR Genome Identification Workflow

This repository contains the current QR-based genome identification workflow. The workflow uses MinHash sketches to represent genome/read data in a compact form, then performs database search, similarity scoring, and confidence classification.

The current confidence thresholds are heuristic and should be refined using additional benchmark testing.

## Main workflow

The basic workflow is:

1. Build a MinHash reference database from FASTA genomes.
2. Generate a MinHash sketch from FASTQ read files or a FASTA query.
3. Search the query sketch against the reference database.
4. Decode a QR payload or QR image and search it against the database.
5. Run benchmark scripts to evaluate top-match accuracy and confidence scoring.

## Main scripts

### `step5_db_search.py`

Main database/search script.

This script can:

* build a reference database from FASTA files
* search with a FASTA file, sketch JSON file, text payload, or direct `VGQR1...` payload
* compute estimated Jaccard similarity
* return ranked database matches

Build a reference database:

```bash
python step5_db_search.py build-db reference_folder --out refs.jsonl --canonical
```

Search using a query sketch:

```bash
python step5_db_search.py search refs.jsonl query_sketch.json --top 10
```

Search using a FASTA query:

```bash
python step5_db_search.py search refs.jsonl query.fasta --top 10
```

---

### `step_reads_minhash.py`

Creates a MinHash sketch from paired FASTQ read files.

Example:

```bash
python step_reads_minhash.py reads_R1.fastq.gz reads_R2.fastq.gz --out reads_sketch.json --canonical
```

The output is a JSON sketch file that can be searched with `step5_db_search.py`.

---

### `qr_decode_search.py`

Decodes a QR image, payload text file, or direct `VGQR1...` payload, then searches it against the viral reference database.

This script currently expects the database file:

```text
viral9_refs.jsonl
```

Example with a QR image:

```bash
python qr_decode_search.py qr_image.png
```

Example with a payload text file:

```bash
python qr_decode_search.py payload.txt
```

Example with a direct payload string:

```bash
python qr_decode_search.py VGQR1....
```

The script reports:

* top match
* top similarity score
* second-best match
* confidence gap
* confidence label
* top 5 matches

---

### `benchmark_viral9.py`

Runs an in-database benchmark using FASTA files in:

```text
viral9_benchmark_refs/
```

It searches each benchmark genome against:

```text
viral9_refs.jsonl
```

Run:

```bash
python benchmark_viral9.py
```

Output:

```text
viral9_benchmark_results.csv
```

---

### `ood_confidence_benchmark.py`

Runs out-of-database/confidence testing using FASTA files in:

```text
ood_queries/
```

It searches each query against:

```text
viral9_refs.jsonl
```

Run:

```bash
python ood_confidence_benchmark.py
```

Output:

```text
ood_confidence_results.csv
```

This is used to test whether the confidence logic can separate:

* in-database matches
* near-database or related matches
* sibling/related genomes
* far or out-of-database genomes

---

### `unzip_ood.py`

Helper script for extracting `.fasta.gz` files inside:

```text
ood_queries/genomes_to_send/
```

Run:

```bash
python unzip_ood.py
```

## Suggested order to run

### 1. Build the reference database

```bash
python step5_db_search.py build-db reference_folder --out viral9_refs.jsonl --canonical
```

### 2. Generate a sketch from read files

```bash
python step_reads_minhash.py reads_R1.fastq.gz reads_R2.fastq.gz --out reads_sketch.json --canonical
```

### 3. Search the read sketch against the database

```bash
python step5_db_search.py search viral9_refs.jsonl reads_sketch.json --top 10
```

### 4. Decode and search a QR payload

```bash
python qr_decode_search.py qr_image.png
```

or

```bash
python qr_decode_search.py payload.txt
```

### 5. Run benchmark tests

```bash
python benchmark_viral9.py
python ood_confidence_benchmark.py
```

## Confidence scoring

The current confidence logic uses:

1. the similarity score of the top match
2. the gap between the top match and the second-best match

The current confidence categories include:

* strong in-database match
* near-database or related match
* ambiguous close match
* out-of-database / no confident match

These thresholds are currently hand-set and should be validated with a larger benchmark panel.

## Notes

Large raw FASTQ/FASTA files are not included in this repository. The repository focuses on the core scripts needed to inspect the workflow and scoring logic.

Some scripts currently expect specific local file or folder names, such as:

```text
viral9_refs.jsonl
viral9_benchmark_refs/
ood_queries/
ood_queries/genomes_to_send/
```

These paths can be changed later to make the workflow more flexible.
******
