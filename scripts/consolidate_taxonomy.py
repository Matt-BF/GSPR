import pandas as pd

df_hosts = pd.read_csv(
    "Host_taxonomy_results_all_methods_new.tsv", sep="\t"
).drop_duplicates(subset="Plasmid")
df_hosts.set_index("Plasmid", inplace=True)

# Apparently the crispr host is not adding any new information from what we had from iphop
# df_hosts[(pd.notna(df_hosts["crispr_blast_host"])) & (pd.isna(df_hosts["iphop_host"]))]


final_hosts = {}
for idx in df_hosts.index:
    if pd.isna(df_hosts.loc[idx, "isolate_host"]):
        if pd.isna(df_hosts.loc[idx, "iphop_host"]) and pd.notna(
            df_hosts.loc[idx, "MMseqs_taxonomy"]
        ):
            final_hosts[idx] = {
                "host": df_hosts.loc[idx, "MMseqs_taxonomy"].replace("_", "__"),
                "method": "mmseqs-taxonomy",
            }
        elif pd.notna(df_hosts.loc[idx, "iphop_host"]):
            final_hosts[idx] = {
                "host": df_hosts.loc[idx, "iphop_host"],
                "method": f'iphop-{df_hosts.loc[idx,"List of methods"].split(";")[0]}',
            }
    elif pd.notna(df_hosts.loc[idx, "isolate_host"]):
        final_hosts[idx] = {
            "host": df_hosts.loc[idx, "isolate_host"],
            "method": "isolate",
        }

df_final_hosts = pd.DataFrame.from_dict(final_hosts, orient="index")

df_final_hosts.reset_index().to_csv(
    "plasmid_host_taxonomy_consolidated_new.tsv", sep="\t", index=False
)
