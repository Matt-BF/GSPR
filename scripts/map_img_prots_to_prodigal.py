from collections import defaultdict

img_contig_to_prot = defaultdict(list)

with open('../helper_tables/IMG_prots.txt') as f:
    for line in f:
        protein = line.strip()
        contig = protein.split('|')[0]
        img_contig_to_prot[contig].append(protein)


old_prot_to_new = {}

for contig,prots in img_contig_to_prot.items():
    for idx, prot in enumerate(prots):
        old_prot_to_new[prot] = f"{prot.split('|')[0]}_{idx+1}"

with open ('IMG_prots_original_to_prodigal.tsv','w') as j:
    for k,v in old_prot_to_new.items():
        j.write(f"{k}\t{v}\n")
