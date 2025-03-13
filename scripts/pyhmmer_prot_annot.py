#!/usr/bin/env python

from pathlib import Path
import pyhmmer
import gzip

hmm_lengths = {}
for p in Path('../dbs').glob("*.hmm.gz"):
    with pyhmmer.plan7.HMMFile(gzip.open(p)) as hmm_file:
        for hmm in hmm_file:
            if hmm.accession:
                hmm_lengths[hmm.accession] = len(hmm.consensus)
            else:
                hmm_lengths[hmm.name] = len(hmm.consensus)

for p in Path("./").glob("soil_isolate_plasmids/soil_isolate_plasmid.faa"):
    prefix = p.stem.replace("_proteins", "")
    with open(f"{prefix}.tsv", "w") as fout:
        with pyhmmer.easel.SequenceFile(p, digital=True) as seq_file:
            seqs = seq_file.read_block()

        with pyhmmer.plan7.HMMFile(gzip.open("../dbs/Pfam-A.hmm.gz")) as hmm_file:
            for hits in pyhmmer.hmmsearch(hmm_file, seqs, bit_cutoffs="gathering"):
                for hit in hits:
                    n_aligned_positions = len(
                        hit.best_domain.alignment.hmm_sequence
                    ) - hit.best_domain.alignment.hmm_sequence.count(".")
                    hmm_coverage = (
                        n_aligned_positions
                        / hmm_lengths[hit.best_domain.alignment.hmm_accession]
                    )
                    fout.write(
                        f"{hits.query_accession.decode()}\t{hit.name.decode()}\t{hit.evalue:.2E}\t{hit.score:.2f}\t{hmm_coverage:.2f}\n"
                    )

        with pyhmmer.plan7.HMMFile(gzip.open("../dbs/NCBIFAM.hmm.gz")) as hmm_file:
            for hits in pyhmmer.hmmsearch(hmm_file, seqs, bit_cutoffs="gathering"):
                for hit in hits:
                    n_aligned_positions = len(
                        hit.best_domain.alignment.hmm_sequence
                    ) - hit.best_domain.alignment.hmm_sequence.count(".")
                    hmm_coverage = (
                        n_aligned_positions
                        / hmm_lengths[hit.best_domain.alignment.hmm_accession]
                    )
                    fout.write(
                        f"{hits.query_accession.decode()}\t{hit.name.decode()}\t{hit.evalue:.2E}\t{hit.score:.2f}\t{hmm_coverage:.2f}\n"
                    )


        with pyhmmer.plan7.HMMFile(gzip.open("../dbs/kofam.hmm.gz")) as hmm_file:
            for hits in pyhmmer.hmmsearch(hmm_file, seqs, Z=10_000, E=1e-8):
                for hit in hits:
                    n_aligned_positions = len(
                        hit.best_domain.alignment.hmm_sequence
                    ) - hit.best_domain.alignment.hmm_sequence.count(".")
                    hmm_coverage = (
                        n_aligned_positions
                        / hmm_lengths[hit.best_domain.alignment.hmm_name]
                    )
                    if hmm_coverage >= 0.6:
                        fout.write(
                            f"{hits.query_name.decode()}\t{hit.name.decode()}\t{hit.evalue:.2E}\t{hit.score:.2f}\t{hmm_coverage:.2f}\n"
                        )
