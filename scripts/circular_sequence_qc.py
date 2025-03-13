#!/usr/bin/env python

import bz2
import fileinput
import gzip
import lzma
import multiprocessing.pool
import subprocess
import textwrap
from contextlib import contextmanager
from copy import deepcopy
from enum import Enum, auto
from pathlib import Path

import kcounter
import numpy as np
import pyrodigal_gv


class Compression(Enum):
    bzip2 = auto()
    gzip = auto()
    xz = auto()
    zstd = auto()
    uncompressed = auto()


def is_compressed(filepath: Path):
    with open(filepath, "rb") as fin:
        signature = fin.peek(8)[:8]
        if tuple(signature[:2]) == (0x1F, 0x8B):
            return Compression.gzip
        elif tuple(signature[:3]) == (0x42, 0x5A, 0x68):
            return Compression.bzip2
        elif tuple(signature[:7]) == (0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00, 0x00):
            return Compression.xz
        elif tuple(signature[:4]) == (0x28, 0xB5, 0x2F, 0xFD):
            return Compression.zstd
        else:
            return Compression.uncompressed


@contextmanager
def open_file(filepath):
    filepath_compression = is_compressed(filepath)
    if filepath_compression == Compression.gzip:
        fin = gzip.open(filepath, "rt")
    elif filepath_compression == Compression.bzip2:
        fin = bz2.open(filepath, "rt")
    elif filepath_compression == Compression.xz:
        fin = lzma.open(filepath, "rt")
    else:
        fin = open(filepath, "r")
    try:
        yield fin
    finally:
        fin.close()


class Sequence:
    def __init__(self, header: str, seq: str, compress: bool = False):
        self._compress = compress
        self._header = header
        if self._compress:
            self._seq = gzip.compress(seq.encode("ascii"), 1)
        else:
            self._seq = seq.encode("ascii")

    @property
    def header(self):
        return self._header

    @property
    def accession(self):
        return self._header.split()[0]

    @property
    def seq(self):
        if self._compress:
            return gzip.decompress(self._seq).decode()
        else:
            return self._seq.decode()

    @property
    def seq_ascii(self):
        return self.seq.upper().encode("ascii")

    def count(self, substring: str):
        return self.seq.count(substring)

    def rc(self):
        tab = self.seq.maketrans("ACTGNactgn", "TGACNtgacn")
        return Sequence(self.header, self.seq.translate(tab)[::-1], self._compress)

    def has_dtr(self, min_length: int = 21):
        substring = self.seq.casefold()[:min_length]
        pos = self.seq.casefold().rfind(substring)
        if pos < len(self) / 2:
            return (False, 0)
        substring = self.seq.casefold()[pos:]
        return (self.seq.casefold()[: len(substring)] == substring, len(substring))

    def has_itr(self, min_len: int = 21, max_len: int = 1000):
        rev = self.rc()
        if self.seq[:min_len] == rev.seq[:min_len]:
            i = min_len + 1
            while self.seq[:i] == rev.seq[:i] and i <= max_len:
                i += 1
            return (True, len(self.seq[: i - 1]))
        else:
            return (False, 0)

    def __str__(self):
        return f">{self.header}\n{textwrap.fill(self.seq, 60)}\n"

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, k: int):
        return Sequence(self.header, self.seq[k], self._compress)

    def __eq__(self, other: object):
        if other.__class__ is self.__class__:
            return self.seq.casefold() == other.seq.casefold()
        elif other.__class__ is str:
            return self.seq.casefold() == other.casefold()
        return NotImplemented

    def __hash__(self):
        return hash(self.seq.casefold())

    def __add__(self, other: object):
        if other.__class__ is not self.__class__:
            return NotImplemented
        compress = other._compress or self._compress
        return Sequence(
            f"{self.accession}+{other.accession}", f"{self.seq}{other.seq}", compress
        )


def read_fasta(uppercase=False, strip_n=False, compress=False):
    with fileinput.input() as fin:
        last = None
        while True:
            if not last:
                for l in fin:
                    if l[0] == ">":
                        last = l[:-1]
                        break
            if not last:
                break
            name, seqs, last = last[1:], [], None
            for l in fin:
                if l[0] == ">":
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            seqs = "".join(seqs)
            if uppercase:
                seqs = seqs.upper()
            if strip_n:
                seqs = seqs.strip("nN")
            if len(seqs):
                yield Sequence(name, seqs, compress)
            if not last:
                break


def longest_consecutive_substring_len(string, char):
    if not string:
        return 0
    if char not in string:
        return 0
    max_length = 1
    current_length = 1
    for i in range(1, len(string)):
        if string[i] == string[i - 1] and string[i] == char:
            current_length += 1
        else:
            current_length = 1
        max_length = max(max_length, current_length)
    return max_length


def dustmasker(record, level=40):
    p = subprocess.Popen(
        [
            "dustmasker",
            "-in",
            "-",
            "-outfmt",
            "fasta",
            "-level",
            str(level),
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    stdout, stderr = p.communicate(str(record))
    if p.returncode:
        raise RuntimeError(stderr)
    header, seq = stdout.strip(">").split("\n", 1)
    seq = seq.replace("\n", "")
    return Sequence(header, seq)


def repeat_match(record, n=20):
    p = subprocess.Popen(
        ["repeat-match", "-n", str(n), "/dev/stdin"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    stdout, stderr = p.communicate(str(record))
    if stdout.startswith("***"):
        return ()
    elif p.returncode:
        raise RuntimeError(stderr)
    output = stdout.strip("\n").split("\n")[2:]
    repeats = []
    for line in output:
        start_1, start_2, length = line.split()
        strand = -1 if start_2.endswith("r") else 1
        start_1, start_2, length = int(start_1), int(start_2.strip("r")), int(length)
        repeats.append((start_1, start_2, length, strand))
    return tuple(repeats)


def merge_intervals(intervals):
    intervals = deepcopy(intervals)
    intervals.sort()
    stack = [intervals[0]]
    for i in intervals[1:]:
        if stack[-1][0] <= i[0] <= stack[-1][-1]:
            stack[-1][-1] = max(stack[-1][-1], i[-1])
        else:
            stack.append(i)
    return stack


def compute_coding_density(record):
    gene_finder = pyrodigal_gv.ViralGeneFinder(meta=True)
    gene_intervals = [
        [pred.begin, pred.end] for pred in gene_finder.find_genes(record.seq_ascii)
    ]
    if len(gene_intervals):
        gene_intervals = merge_intervals(gene_intervals)
        c = sum(interval[1] - interval[0] + 1 for interval in gene_intervals)
        return c / len(record)
    else:
        return 0.0


def process_sequence(record):
    dtr, itr, concatemer = False, False, False
    if repeats := repeat_match(record, 1000):
        longest_repeat = max(i[2] for i in repeats)
        if longest_repeat < len(record) * 0.9:
            concatemer = True
    dtr, dtr_length = record.has_dtr()
    if dtr:
        masked_record = dustmasker(record)
        if (
            masked_record.seq.upper()[:dtr_length].count("N") >= dtr_length / 2
            or dtr_length >= 1_000
        ):
            dtr = False
    else:
        itr, itr_length = record.has_itr()
        if itr:
            masked_record = dustmasker(record)
            if (
                masked_record.seq.upper()[:itr_length].count("N") >= itr_length / 2
                or itr_length >= 1_000
            ):
                itr = False
    coding_density = compute_coding_density(record)
    avg_kmer = kcounter.count_kmers(record.seq.upper(), 21, canonical_kmers=True)
    avg_kmer = np.mean(list(avg_kmer.values()))
    return (
        record.accession,
        len(record),
        record.seq.upper().count("N"),
        longest_consecutive_substring_len(record.seq.upper(), "N"),
        coding_density,
        avg_kmer,
        concatemer,
        dtr,
        itr,
    )


print(
    "seq_name\tseq_length\tn_ambiguous\tlongest_ambiguous\t"
    "coding_density\tavg_kmer\tconcatemer\tdtr\titr"
)
with multiprocessing.pool.Pool() as pool:
    for (
        accession,
        length,
        n_ambiguous,
        longest_ambiguous,
        coding_density,
        avg_kmer,
        concatemer,
        dtr,
        itr,
    ) in pool.imap(process_sequence, read_fasta()):
        print(
            f"{accession}\t{length}\t{n_ambiguous}\t{longest_ambiguous}\t"
            f"{coding_density:.4f}\t{avg_kmer:.4f}\t{concatemer}\t{dtr}\t{itr}"
        )
