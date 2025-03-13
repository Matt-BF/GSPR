#!/usr/bin/env python

import gzip
from pathlib import Path
import sys
import pyhmmer


def run_phhmer(query_seqs, db_seqs):
    phmmer_results = []

    search = pyhmmer.phmmer(query_seqs, db_seqs, E=1e-8)

    for hits in search:
        for hit in hits:
            phmmer_results.append(
                (
                    hits.query_name.decode(),
                    hit.name.decode(),
                    hit.evalue,
                    hit.score,
                    hit.description.decode(),
                )
            )
    return phmmer_results


p = Path(sys.argv[1])
prefix = p.stem.replace("_proteins", "")
print(f"\nrunning for {prefix}")
with open(f"{prefix}_phmmer.tsv", "w") as fout:
    fout.write("Query_name\tHit_name\tHit_score\tHit_evalue\tHit_description\n")

    with pyhmmer.easel.SequenceFile(p, digital=True) as seq_file:
        seqs = seq_file.read_block()

    with pyhmmer.easel.SequenceFile(
        "../dbs/iceberg/all_ICEs.faa", digital=True
    ) as db_file:
        db_seqs = db_file.read_block()

    phmmer_results = run_phhmer(seqs, db_seqs)
    print('Writing results')
    for result in phmmer_results:
        fout.write(
            f"{result[0]}\t{result[1]}\t{result[2]:.2E}\t{result[3]:.2f}\t{result[4]}\n"
        )
