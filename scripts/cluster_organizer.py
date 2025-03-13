import sys
from typing import Dict,List

with open(sys.argv[1]) as f:
    clusters:Dict[str,List] = {}

    for line in f:
        line = line.strip().split()
        if line[0] not in clusters.keys():
            clusters[line[0]] = []
            clusters[line[0]].append(line[1])
        else:
            clusters[line[0]].append(line[1])
with open(f"{sys.argv[1]}.fixed", "w") as j:
    j.write("cluster_representative\tcluster_members\tcluster_length\n")
    for k in clusters:
        j.write(f"{k}\t{','.join(clusters[k])}\t{len(clusters[k])}\n")
