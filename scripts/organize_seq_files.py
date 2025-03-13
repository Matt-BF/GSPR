#!/usr/bin/env python

import os
import shutil
import tempfile
import sys
from natsort import natsorted
from pathlib import Path
from tqdm import tqdm


# Global parameters
INPUT_DIRECTORY = Path(sys.argv[1])
OUTPUT_DIRECTORY = Path(sys.argv[2])
MIN_N_SEQS_FILE = int(sys.argv[3])


def read_fasta(fp, uppercase=False, strip_n=False):
    with open(fp) as fin:
        last = None
        while True:
            if not last:
                for l in fin:
                    if l[0] == ">":
                        last = l[:-1]
                        break
            if not last: break
            name, seqs, last = last[1:].partition(" ")[0], [], None
            for l in fin:
                if l[0] == ">":
                    last = l[:-1]
                    break
                seqs.append(l[:-1].upper())
            seqs = "".join(seqs)
            if uppercase:
                seqs = seqs.upper()
            if strip_n:
                seqs = seqs.strip("nN")
            yield name, seqs
            if not last: break


input_file_list = [i for i in INPUT_DIRECTORY.iterdir() if os.path.getsize(i) > 0]
input_file_list = natsorted(input_file_list)

if OUTPUT_DIRECTORY.exists():
    shutil.rmtree(OUTPUT_DIRECTORY)
OUTPUT_DIRECTORY.mkdir()

with tempfile.TemporaryDirectory() as tmpd:
    tmpf = open(Path(tmpd).joinpath("tmpfile"), "w")
    first_file = input_file_list[0].stem
    counter = 0
    for input_file in tqdm(input_file_list, ncols=100, smoothing=0):
        if not first_file:
            first_file = input_file.stem
        for header, seq in read_fasta(input_file, uppercase=True, strip_n=False):
            counter += 1
            tmpf.write(f">{header}\n{seq}\n")
        if counter >= MIN_N_SEQS_FILE:
            last_file = input_file.stem
            output_file = OUTPUT_DIRECTORY.joinpath(f"{first_file}_{last_file}.fna")
            tmpf.close()
            shutil.copy(tmpf.name, output_file)
            Path(tmpd).joinpath("tmpfile").unlink()
            tmpf = open(Path(tmpd).joinpath("tmpfile"), "w")
            first_file = None
            counter = 0
    last_file = input_file.stem
    output_file = OUTPUT_DIRECTORY.joinpath(f"{first_file}_{last_file}.fna")
    tmpf.close()
    shutil.copy(tmpf.name, output_file)
