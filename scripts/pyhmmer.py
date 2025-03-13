import bz2
import gzip
import lzma
import textwrap
from contextlib import contextmanager
from enum import Enum, auto
from pathlib import Path
import pyhmmer


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
            return False
        substring = self.seq.casefold()[pos:]
        return self.seq.casefold()[: len(substring)] == substring

    def has_itr(self, min_len: int = 21):
        rev = self.rc().seq
        return self.seq.casefold()[:min_len] == rev.casefold()[:min_len]

    def __str__(self):
        return f">{self.header}\n{textwrap.fill(self.seq, 60)}\n"

    def __repr__(self):
        if len(self) > 40:
            start = self.seq[:34]
            end = self.seq[-3:]
            seq = f"{start}...{end}"
        else:
            seq = self.seq
        return f"Sequence({self.accession}, {seq})"

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


def read_fasta(filepath, uppercase=False, strip_n=False, compress=False):
    with open_file(filepath) as fin:
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


with pyhmmer.easel.SequenceFile("test.faa", digital=True) as seq_file:
    seqs = seq_file.read_block()

with pyhmmer.plan7.HMMFile("Pfam-A.h3m") as hmm_file:
    for hits in pyhmmer.hmmsearch(hmm_file, seqs, bit_cutoffs="gathering"):
        if len(hits) >= 40:
            print(f"HMM {hits.query_name.decode()} found {len(hits)} hits in the target sequences")
            break