from collections import defaultdict

gtdb_taxons = defaultdict(set)
soil_plasmid_taxons = defaultdict(set)

with open('gtdb_ncbi_isolation_soil_taxonomy.txt') as f:
    for line in f:
        line = line.strip()
        for tax in line.split(';'):
            gtdb_taxons[tax[0:3]].add(tax[3:])

with open('soil_plasmid_hosts.txt') as f:
    for line in f:
        line = line.strip()
        for tax in line.split(';'):
            soil_plasmid_taxons[tax[0:3]].add(tax[3:])

print('GTDB-Plasmid')
for k, v in gtdb_taxons.items():
    print(k)
    print(f'taxon_size {len(v)}')
    print(f'intersection size gtdb-plasmids: {len(v.intersection(soil_plasmid_taxons[k]))}')


print('Plasmid-GTDB')
for k,v in soil_plasmid_taxons.items():
    print(k)
    print(f'taxon_size {len(v)}')
    print(f'intersection size plasmids-gtdb: {len(v.intersection(gtdb_taxons[k]))}')
    print(f'{v-gtdb_taxons[k]}')

