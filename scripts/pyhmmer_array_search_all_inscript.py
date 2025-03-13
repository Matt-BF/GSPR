#!/usr/bin/env python

import gzip
from pathlib import Path
import sys
import pyhmmer


hmm_lengths = {}
hmm_list = []

def run_hmmscan(hmm, seqs, has_gathering=False):
    hmm_results = {}
    hmm_results[hmm.parts[-1].replace('.hmm.gz','')] = []

    with pyhmmer.plan7.HMMFile(gzip.open(hmm)) as hmm_file:
            if has_gathering:
                search = pyhmmer.hmmsearch(hmm_file, seqs, bit_cutoffs="gathering")
            else:
                search = pyhmmer.hmmsearch(hmm_file, seqs, Z=10_000, E=1e-8)
                
            for hits in search:
                for hit in hits:
                    n_aligned_positions = len(
                        hit.best_domain.alignment.hmm_sequence
                    ) - hit.best_domain.alignment.hmm_sequence.count(".")
                    try:
                        hmm_coverage = (
                            n_aligned_positions
                            / hmm_lengths[hit.best_domain.alignment.hmm_accession.decode()]
                        )
                    except KeyError:
                        hmm_coverage = (
                            n_aligned_positions
                            / hmm_lengths[hit.best_domain.alignment.hmm_name.decode()])

                    if hmm_coverage > 0.6:
                        try:
                            hmm_results[hmm.parts[-1].replace('.hmm.gz','')].append([hits.query_accession.decode(), hit.name.decode(), hit.evalue, hit.score, hmm_coverage])
                        except AttributeError:
                            hmm_results[hmm.parts[-1].replace('.hmm.gz','')].append([hits.query_name.decode(), hit.name.decode(), hit.evalue, hit.score, hmm_coverage])
    return hmm_results                    

for p in Path("../dbs/").glob("*.hmm.gz"):
    hmm_list.append(p)
    with pyhmmer.plan7.HMMFile(gzip.open(p)) as hmm_file:
        print(f"Opening {p.parts[-1]}")
        for hmm in hmm_file:
            if hmm.accession:
                hmm_lengths[hmm.accession.decode()] = len(hmm.consensus)
            else:
                hmm_lengths[hmm.name.decode()] = len(hmm.consensus)

with open('remaining_paths.txt') as f:
    for p in f:
        p = p.strip()
        prefix = Path(p).stem.replace("_proteins", "")
        print(f"\nrunning for {prefix}")
        with open(f"hmmsearch_outputs_5kb/{prefix}.tsv", "w") as fout:
            fout.write("Query_name\tHit_name\tHit_evalue\tHit_score\tHMM_coverage\n")

            with pyhmmer.easel.SequenceFile(p, digital=True) as seq_file:
                seqs = seq_file.read_block()
            
            print("Running HMM Searches")
            for hmm in hmm_list:
                print(f"Doing {hmm.parts[-1]}")
                if "PFAM" in hmm.parts[-1] or "NCBIFAM" in hmm.parts[-1]:
                    hmmscan = run_hmmscan(hmm, seqs, has_gathering=True)
                else:
                    hmmscan = run_hmmscan(hmm, seqs)
                
                print("Writing results")
                for hit_list in hmmscan.values():
                    for hit_sublist in hit_list:
                        fout.write(f"{hit_sublist[0]}\t{hit_sublist[1]}\t{hit_sublist[2]:.2E}\t{hit_sublist[3]:.2f}\t{hit_sublist[4]:.2f}\n")

