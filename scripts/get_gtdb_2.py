import requests
from tqdm import tqdm

with open("isolates_assembly_accessions.txt") as f:
    accessions = {i.strip().split()[0]: i.strip().split()[1] for i in f}


results = {}
errors = {}

with open("isolate_host_update.txt", "w") as j:
    with tqdm(total=len(accessions)) as pbar:
        for k, v in accessions.items():
            response = requests.get(
                f"https://gtdb-api.ecogenomic.org/search/gtdb?search={v}&page=1&itemsPerPage=100&searchField=all&gtdbSpeciesRepOnly=false&ncbiTypeMaterialOnly=false"
            )
            try:
                j.write(f"{k}\t{response.json()['rows'][0]['gtdbTaxonomy']}\n")
            except Exception as e:
                errors[k] = e
            pbar.update(1)

with open("errors.txt", "w") as f:
    for k, v in errors.items():
        f.write(f"{k}\t{v}\n")
