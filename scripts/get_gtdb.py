import requests

with open('NCBI_missing_isolates.txt') as f:
    accessions = {i.strip().split()[0]:i.strip().split()[1] for i in f}


results = {}

with open('isolate_host_fix.txt', 'w') as j:
    for k,v in accessions.items():
        response = requests.get(f'https://gtdb-api.ecogenomic.org/search/gtdb?search={v}&page=1&itemsPerPage=100&searchField=all&gtdbSpeciesRepOnly=false&ncbiTypeMaterialOnly=false')
        j.write(f"{k}\t{response.json()['rows'][0]['gtdbTaxonomy']}\n")

