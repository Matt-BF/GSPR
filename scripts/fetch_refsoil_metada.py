from Bio import Entrez
from Bio import SeqIO

with open("../novelty/refsoil_plasmids_no_dups.txt") as f:
    ids = f.read().splitlines()

Entrez.email = "mbfiamenghi@lbl.gov"  # Replace with your email address
handle = Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="gb", retmode="text")
records = list(SeqIO.parse(handle, "genbank"))

metadata = {}
for record in records:
    metadata[record.id] = {"Country": "Unknown", "Latitude/Longitude": "Unknown"}
    for feature in record.features:
        if feature.type == "source":
            qualifiers = feature.qualifiers
            country = qualifiers.get("geo_loc_name")
            if country:
                metadata[record.id]["Country"] = country[0]
            if "lat_lon" in qualifiers:
                metadata[record.id]["Latitude/Longitude"] = qualifiers["lat_lon"][0]

with open("../novelty/refsoil_efetch_metadata.tsv", "w") as h:
    h.write("Accession\tCountry\tLatitude/Longitude\n")
    for k, v in metadata.items():
        h.write(f"{k}\t{v['Country']}\t{v['Latitude/Longitude']}\n")


#     metadata = {}
#     for r in record:
#         print(r)
#         metadata[r.id] = {"Country": "Unknown", "Latitude/Longitude": "Unknown"}
#         for feature in r.features:
#             if feature.type == "source":
#                 qualifiers = feature.qualifiers
#                 country = qualifiers.get("geo_loc_name")
#                 if country:
#                     metadata["Country"] = country[0].split(":")[0]
#                 if "lat_lon" in qualifiers:
#                     metadata["Latitude/Longitude"] = qualifiers["lat_lon"][0]

#     return metadata
