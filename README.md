# GSPR
This repository contains notebooks, scripts and other helper files used in the work published at [link]

Please note many paths are still hardcoded. I am slowly fixing this to point mainly to database file (see below) or otherwise indicate how to (re)generate the original file for analysis. Please open an issue if you have any inquiries.


## Data
Data generated in this work is available as a [DuckDB](https://duckdb.org/) database, available at the IMG website or in Zenodo.

The database has the following tables:

`duckdb soil_plasmids -c "SHOW TABLES"`
- **contig_qc** &rarr; QC results for the metaG/metaT contigs
- **taxon_metadata** &rarr; various metadata associated with a taxon_oid, for instance ecosystem, Lat/Lon, number of plasmids, soil classification, soil attributes, etc.
- **genomad_results_unfiltered** &rarr; plasmids detected via geNomad before any kind of filtering described in the paper.
- **genomad_results** &rarr; filtered geNomad results as described in the paper.
- **ptu_derep** &rarr; plasmids dereplicated and clustered into pOTUs, with additional ecosystem and origin (Meta vs Isolate) information.
- **hmm_outputs** &rarr; plasmid protein annotations using HMMs from different databases (PFam, KEGG, NCBIFam, TIGRFam).
- **host** &rarr; plasmids with detected host, the host taxonomy and prediction method.
- **iphop_blast_sources** &rarr; for hosts predicted via iphop-BLAST, where is the source genome located (IMG, GTDB, etc).
- **plasmid_num_genes** &rarr; for each plasmid, how many genes were found and how many had prediction/annotations.
- **conj_types** &rarr; predicted conjugation machinery using CONJScan.
- **plasmid_cripr** &rarr; CRISPR detected in plasmids with cctyper.
- **crispr_hits** &rarr; Consolidated results based on BLAST hits of plasmid CRISPR spacers against GSPR plasmids and viruses.
- **crispr_hits_assemblies** &rarr; Consolidated results based on BLAST hits of soil contig spacers against GSPR plasmids and viruses.
- **eggnog** &rarr; annotations using EggNOG.
- **amr** &rarr; Antimicrobial Resistance Genes predicted with AMRFinderPlus
- **bgc** &rarr; Biosynthetic Gene Clusters predicted via AntiSMASH and consolidated with BiG-SCAPE.
- **cog** &rarr; COG categories and COG genes extracted from EggNOG results.
- **cazy** &rarr; Cazymes in plasmids, predicted via dbCAN.
- **amp** &rarr; Antimicrobial Peptides (AMPs) in plasmids, predicted via Macrel.


## Notebooks
Analysis was done in several notebooks
- **amp** &rarr; Analysis of the Antimicrobial peptides results.
- **amr_results** &rarr; Analysis of the Antimicrobial resistance results.
- **bgc** &rarr; Analysis of the BGC results.
- **cluster_analysis** &rarr; Analysis of the pOTU clustering.
- **cog** &rarr; Analysis of the COG results.
- **complete_size_host** &rarr; Analysis of the correlations between host origin (Meta or Isolate), host size, plasmid origin and plasmid sizes.
- **conjscan** &rarr; Analysis of the CONJScan results (conjugation machinery).
- **consolidate_metadata** &rarr; Notebook that consolidates metadata from IMG/PR, MGnify, PLSDB and Refsoil, which results in the `taxon_metadata` table in the db.
- **contig_qc** &rarr; Analysis of the Contig QC step, with exploratory analysis to select thresholds for some contig filters (such as `avg_kmer` and `concatemer`) as described in the manuscript.
- **dbcan** &rarr; Analysis of the Cazymes results.
- **eggnog** &rarr; Analysis of the EggNOG results.
- **filter_soil_rhizo_IMG** &rarr; Notebook demonstrating the filters used to keep soil and rhizosphere samples from IMG metadata according to the exclusion terms mentioned in the manuscript
- **filter_soil_rhizo_MGnify** &rarr; Same as above, but for MGnify.
- **gtdb_plasmid_intersection** &rarr; Analysis of soil species intersection between GTDB soil taxa and GSPR hosts.
- **hosts** &rarr; Analysis of the GSPR predicted hosts.
- **hosts_environment** &rarr; Analysis of which soil habitats GSPR predicted hosts are found, and comparison with which habitats the pOTUs are found.







