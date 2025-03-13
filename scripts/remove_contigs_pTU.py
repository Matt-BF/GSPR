with open('exclude_contig_list.txt') as f:
    exclude = [i.strip() for i in f]

with open('pTUs_clusters.tsv') as f:
    with open('pTUs_clusters_fixed.tsv','w') as fout:
        fout.write('cluster_representative\tcluster_members\t\n')
        next(f)
        for line in f:
            line = line.strip().split('\t')

            if line[0] in exclude:
                print(f'{line[0]} in exclude')
                continue 
            else:
                cluster = []
                for member in line[1].split(','):
                    if member not in exclude:
                        cluster.append(member)
                fout.write(f'{line[0]}\t{",".join(cluster)}\n')
