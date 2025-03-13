rule get_from_db:
    input:
        db_file = "../../soil_plasmid.db",
        plasmid_prot_ids = "ids.txt"
    params:
        pfam_arg = lambda wildcards: wildcards.pfam_id
    output:
        "{pfam_id}.txt"
    shell:
        """
        duckdb {input.db_file} -c "COPY (SELECT Hit_name FROM hmm_outputs WHERE Query_name LIKE '{params.pfam_arg}%') TO '{output}'"
        rg -Ff {input.plasmid_prot_ids} {output} | sponge {output}
        """


rule separate_meta_and_isolates:
    input:
        pfam_file = "{pfam_id}.txt"
    output:
        meta = "{pfam_id}_meta.txt.sep",
        isolates = "{pfam_id}_isolates.txt.sep"
    shell:
        """
        rg 'IMGPR|Refsoil|PLSDB' {input.pfam_file} > {output.isolates}
        rg -v 'IMGPR|Refsoil|PLSDB' {input.pfam_file} > {output.meta}
        """

rule fetch_fasta:
    input:
        meta = "{pfam_id}_meta.txt.sep",
        isolates = "{pfam_id}_isolates.txt.sep"
    output:
        meta_fasta = "{pfam_id}_meta.faa",
        isolates_fasta = "{pfam_id}_isolates.faa"
    shell:
        """
        seqkit grep -f {input.meta} ../../ptu_derep/all_soil_plasmids_derep.faa.gz > {output.meta_fasta}
        seqkit grep -f {input.isolates} ../../ptu_derep/all_soil_plasmids_derep.faa.gz > {output.isolates_fasta}
        """

rule run_cdhit_and_concat:
    input:
        meta_fasta = "{pfam_id}_meta.faa",
        isolates_fasta = "{pfam_id}_isolates.faa"
    output:
        cdhit_results = "{pfam_id}_cdhit.faa"
    shell:
        """
        cd-hit -c 1 -aS 1 -i {input.meta_fasta} -o {input.meta_fasta}.cdhit.faa
        cd-hit -c 1 -aS 1 -i {input.isolates_fasta} -o {input.isolates_fasta}.cdhit.faa

        cat {input.meta_fasta}.cdhit.faa {input.isolates_fasta}.cdhit.faa > {output.cdhit_results}
        rm {input.meta_fasta}.cdhit.faa {input.isolates_fasta}.cdhit.faa
        """

rule run_hmmalign:
    input:
        cdhit_file = "{pfam_id}_cdhit.faa",
        hmm_file = "{pfam_id}.hmm"
    output:
        hmmalign_file = "{pfam_id}_hmmalign.faa",
    shell:
        """
        hmmalign --trim -o {input.cdhit_file}.sto {input.hmm_file} {input.cdhit_file}
        seqmagick convert {input.cdhit_file}.sto {output.hmmalign_file}
        rm {input.cdhit_file}.sto
        """

rule clipkit:
    input:
        hmmalign_file = "{pfam_id}_hmmalign.faa"
    output:
        trimmed_file = "{pfam_id}_hmmalign.faa.clipkit"
    shell:
        """
        clipkit {input.hmmalign_file} -m kpic-smart-gap
        """

rule fasttree:
    input:
        trimmed_file = "{pfam_id}_hmmalign.faa.clipkit"
    output:
        tree = "{pfam_id}.tree"
    shell:
        """
        FastTreeMP {input.trimmed_file} > {output.tree}
        """
