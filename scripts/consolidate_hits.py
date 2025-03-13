#!/usr/bin/env python
#this is just the first part of stephen's script to assign host per crispr spacer

import os, csv, gzip, argparse
from collections import defaultdict
from operator import itemgetter

def yield_hits(input):
    
    myid = None
    hits = []
    handle = gzip.open(input, "rt") if input.endswith(".gz") else open(input)
    
    for line in handle:
        rec = line.split()
        virus_id, spacer_id = rec[0:2]
        if int(rec[4]) + int(rec[5]) > 1: continue # mismatches + gaps > 1
        if float(rec[3])/float(rec[-1]) < 0.95: continue # aln coverage < 95%
        if int(rec[3]) < 25: continue # aln < 25 bp
        myid = virus_id
        hits.append(spacer_id)
        break

    for line in handle:
        rec = line.split()
        virus_id, spacer_id = rec[0:2]
        if int(rec[4]) + int(rec[5]) > 1: continue
        if float(rec[3])/float(rec[-1]) < 0.95: continue
        if int(rec[3]) < 25: continue
        if virus_id == myid:
            hits.append(spacer_id)
        else:
            yield myid, hits
            myid = virus_id
            hits = [spacer_id]
            
    yield myid, hits



with open('blastn_consolidated.tsv', 'w') as output:
    for plasmid_id, hits in yield_hits('blastn.tsv'):
        output.write(f'{plasmid_id}\t{",".join([hit for hit in hits])}\n')
