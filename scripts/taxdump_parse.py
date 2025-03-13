import taxopy
import sys
plsdb_taxons = {}
taxdb_ncbi = taxopy.TaxDb(keep_files=True)
taxdb_gtdb = taxopy.TaxDb(
    names_dmp="../phylo/mob_analysis/gtdb_taxdump/R214.1/names.dmp",
    nodes_dmp="../phylo/mob_analysis/gtdb_taxdump/R214.1/nodes.dmp",
    merged_dmp="../phylo/mob_analysis/gtdb_taxdump/R214.1/merged.dmp",
    keep_files=True,
)

taxon = sys.argv[1]
ncbi_taxon = taxopy.Taxon(int(taxon), taxdb_ncbi)
for rank, name in ncbi_taxon.rank_name_dictionary.items():
    if taxopy.taxid_from_name(name, taxdb_gtdb):
        gtdb_taxid = taxopy.taxid_from_name(name, taxdb_gtdb)
        print(f"{taxopy.Taxon(int(gtdb_taxid[0]), taxdb_gtdb)}")
        break
    elif taxopy.taxid_from_name(name + "_B", taxdb_gtdb):
        gtdb_taxid = taxopy.taxid_from_name(name + "_B", taxdb_gtdb)
        print(f"{k}\t{taxopy.Taxon(int(gtdb_taxid[0]), taxdb_gtdb)}")
        break
