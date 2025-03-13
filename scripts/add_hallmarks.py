from collections import Counter
import pandas as pd

with open("contigs_with_hallmarks.txt") as f:
    contigs = [i.strip() for i in f]
    hallmarks = Counter(contigs)

df = pd.read_csv("all_plasmid_summary.tsv", sep="\t")

print("adding hallmarks")
df["n_hmm_hallmarks"] = df["seq_name"].apply(lambda x: hallmarks.get(x, 0))
# for i in df.index:
#   df.loc[i, "n_hmm_hallmarks"] = hallmarks.get(i, 0)

df.to_csv("all_plasmid_summary_hallmarks.tsv", sep="\t", index=None)
