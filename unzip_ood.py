from pathlib import Path
import gzip
import shutil

src = Path("ood_queries/genomes_to_send")

for gz_file in src.glob("*.fasta.gz"):
    out_file = gz_file.with_suffix("")  # removes .gz

    with gzip.open(gz_file, "rb") as fin:
        with open(out_file, "wb") as fout:
            shutil.copyfileobj(fin, fout)

    print("Extracted:", out_file.name)

print("Done.")