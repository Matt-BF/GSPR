from Bio import Entrez

Entrez.email = "mbfiamenghi@lbl.gov"  # Always tell NCBI who you are

id_list = []
with open("plsdb_all_soil_plasmids.txt") as f:
    next(f)
    for line in f:
        id_list.append(line.strip())
database = "nuccore"  # Change this to the appropriate database

# Fetch records for the IDs
print("fetching...")
handle = Entrez.efetch(db=database, id=id_list, retmode="xml")
records = Entrez.read(handle)

# Extract taxIDs
with open("plsdb_all_soil_plasmids_taxid.txt", "w") as f:
    for record in records:
        # The path to the taxID can vary by database and record format.
        # You may need to adjust the code based on your specific requirements.
        taxID = record["TaxId"]
        f.write(f"{record['Id']}\t{taxID}\n")
