#!/usr/bin/env python3

import argparse
import re
import subprocess as sp
import tempfile
from collections import OrderedDict
from pathlib import Path

N_THREADS = 4
HIT_COVERAGE = 0.3
HIT_COVERAGE_TRUNCATED = 0.8
FEATURE_END_UNKNOWN = "?"
FEATURE_END_5_PRIME = "5-prime"
FEATURE_END_3_PRIME = "3-prime"
RE_MULTIWHITESPACE = re.compile(r"\s+")
FEATURE_R_RNA = "rRNA"
STRAND_FORWARD = "+"
STRAND_REVERSE = "-"
TMP_DIR = tempfile.TemporaryDirectory()


def predict_r_rnas(cm_path: Path, contigs_path: Path):
    """Search for ribosomal RNA sequences."""

    output_path = Path(TMP_DIR.name) / "rrna.tsv"
    cmd = [
        "cmscan",
        "--noali",
        "--cut_ga", # use CM's GA gathering cutoffs as reporting thresholds
        "-g", # activate glocal mode
        "--nohmmonly", # strictly use CM models
        "--rfam",
        "--cpu",
        str(N_THREADS),
        "--tblout",
        str(output_path),
        str(cm_path),
        str(contigs_path),
    ]
    proc = sp.run(
        cmd,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True,
    )
    if proc.returncode != 0:
        print("stdout='%s', stderr='%s'", proc.stdout, proc.stderr)
        print("rRNAs failed! cmscan-error-code=%d", proc.returncode)
        raise Exception(f"cmscan error! error code: {proc.returncode}")

    rrnas = []
    with output_path.open() as fh:
        for line in fh:
            if line[0] != "#":
                (
                    subject,
                    accession,
                    contig_id,
                    contig_acc,
                    mdl,
                    mdl_from,
                    mdl_to,
                    start,
                    stop,
                    strand,
                    trunc,
                    passed,
                    gc,
                    bias,
                    score,
                    evalue,
                    inc,
                    description,
                ) = RE_MULTIWHITESPACE.split(line.strip(), maxsplit=17)

                if strand == "-":
                    (start, stop) = (stop, start)
                (start, stop) = (int(start), int(stop))
                evalue = float(evalue)
                score = float(score)
                length = stop - start + 1
                if trunc == "5'":
                    truncated = FEATURE_END_5_PRIME
                elif trunc == "3'":
                    truncated = FEATURE_END_3_PRIME
                else:
                    truncated = None

                if accession == "RF00001":
                    rrna_tag = "5S"
                    consensus_length = 119
                elif accession == "RF00177":
                    rrna_tag = "16S"
                    consensus_length = 1533
                elif accession == "RF02541":
                    rrna_tag = "23S"
                    consensus_length = 2925
                elif accession == "RF01959":
                    rrna_tag = "16S archaeal"
                    consensus_length = 1477
                elif accession == "RF02540":
                    rrna_tag = "23S archaeal"
                    consensus_length = 2987
                else:
                    continue

                coverage = length / consensus_length
                if coverage < HIT_COVERAGE_TRUNCATED:
                    truncated = FEATURE_END_UNKNOWN

                if coverage >= HIT_COVERAGE:
                    rrna = OrderedDict()
                    rrna["type"] = FEATURE_R_RNA
                    rrna["contig"] = contig_id
                    rrna["start"] = start
                    rrna["stop"] = stop
                    rrna["strand"] = STRAND_FORWARD if strand == "+" else STRAND_REVERSE
                    rrna["accession"] = accession
                    if accession == "RF00001":
                        rrna["gene"] = "rrf"
                    elif accession in {"RF00177", "RF01959"}:
                        rrna["gene"] = "rrs"
                    elif accession in {"RF02541", "RF02540"}:
                        rrna["gene"] = "rrl"

                    if truncated is None:
                        rrna["product"] = f"{rrna_tag} ribosomal RNA"
                    elif truncated == FEATURE_END_UNKNOWN:
                        rrna["product"] = f"(partial) {rrna_tag} ribosomal RNA"
                    elif truncated == FEATURE_END_5_PRIME:
                        rrna["product"] = f"(5' truncated) {rrna_tag} ribosomal RNA"
                    elif truncated == FEATURE_END_3_PRIME:
                        rrna["product"] = f"(3' truncated) {rrna_tag} ribosomal RNA"

                    if truncated:
                        rrna["truncated"] = truncated

                    rrna["coverage"] = coverage
                    rrna["score"] = score
                    rrna["evalue"] = evalue

                    rrnas.append(rrna)

    return rrnas


def parse_cli():
    argp = argparse.ArgumentParser(description="Find rRNAs")
    argp.add_argument(
        "cm_path",
        type=Path,
        help="Covariation models that will be used by cmscan.",
    )
    argp.add_argument(
        "contigs_path",
        type=Path,
        help="FASTA file containing the contigs that will be processed with cmscan.",
    )
    argp.add_argument(
        "output_path",
        type=Path,
        help="Path to the output file.",
    )
    opts = argp.parse_args()
    return argp, opts


def main():
    argp, opts = parse_cli()
    rrnas = predict_r_rnas(
        opts.cm_path,
        opts.contigs_path,
    )
    with open(opts.output_path, "w") as fout:
        for i in rrnas:
            fout.write(
                f"{i['contig']}\t{i['accession']}\t{i['product']}\t{i['start']}\t"
                f"{i['stop']}\t{i['coverage']:.4f}\t{i['score']}\t{i['evalue']}\n"
            )


if __name__ == "__main__":
    main()
