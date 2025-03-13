import requests
from Bio import Entrez, SeqIO

with open('refsoil_plasmids.txt') as f:
    plasmids = [i.strip() for i in f]


Entrez.email = "mbfiamenghi@lbl.gov"  # Always tell NCBI who you are
with open('refsoil_plasmid_origin.tsv', 'w') as f:
    for plasmid in plasmids:
        stream = Entrez.efetch(db="nucleotide", id=plasmid, rettype="gb", retmode="text")
        record = SeqIO.read(stream, "genbank")
        stream.close()
        print(f"{record.id}\t{record.description}")

        f.write(f"{record.id}\t{record.description}\n")
