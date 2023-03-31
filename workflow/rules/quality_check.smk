# rule checkv_database:
#     output:
#         directory("./"),
#     threads: 1
#     log:
#         "logs/checkv_database/checkv_database.log"
#     envmodules:
#         "plgrid/tools/checkv/0.9.0"
#     resources:
#         partition="plgrid-short",
#         nodes=1,
#         ntasks=2,
#         mem_mb="64GB",
#         time="01:00:00",
#     shell: """
#         checkv download_database {output}
#         """

rule checkv:
    input:
        ancient_viruses_fasta="../results/04-aDNA_authentication/ancient_viruses_fasta/{sample}_ancient_viruses.fasta",
    output:
        outdir=directory("../results/05-quality_check/{sample}_checkv"),
        summary="../results/05-quality_check/{sample}_checkv/quality_summary.tsv"
    params:
        db="../resources/checkv-db-v1.2"
    threads: 1 # it is a problem in checkv, if sample have few viral sequences resulted "Error: 1(and other numbers as well) prodigal tasks failed. Program should be rerun."
    log:
        "logs/checkv/{sample}.log"
    envmodules:
        "plgrid/tools/checkv/0.9.0"
    resources:
        partition="plgrid-short",
        nodes=1,
        ntasks=1,
        mem_mb="64GB",
        disk_mb="32GB",
        time="01:00:00",
    shell: """
        checkv end_to_end {input.ancient_viruses_fasta} {output.outdir} -d {params.db} -t {threads}

        """
rule checkv2fasta:
    input:
        checkv_summary="../results/05-quality_check/{sample}_checkv/quality_summary.tsv",
        viruses_fasta="../results/04-aDNA_authentication/ancient_viruses_fasta/{sample}_ancient_viruses.fasta"
    output:
        checkv_filtered_csv="../results/05-quality_check/checkv_filtered_csv/{sample}_checkv_filtered.csv",
        checkv_filtered_fasta="../results/05-quality_check/checkv_filtered_fasta/{sample}_checkv_filtered.fasta",
    params:
        sample_name="{sample}"
    threads: 1
    log:
        "logs/checkv2fasta/{sample}.log"
    resources:
        partition="plgrid-short",
        nodes=1,
        ntasks=1,
        mem_mb="4GB",
        time="00:05:00",
    script:
        "../scripts/checkv2fasta.py"

rule concat:
    input:
        all_csv="../results/05-quality_check/checkv_filtered_csv/",
        all_fasta="../results/05-quality_check/checkv_filtered_fasta/"
    output:
        out_all_csv="../results/05-quality_check/checkv_filtered_csv/all_ancient_viruses_checkv_filtered.csv",
        out_all_fasta="../results/05-quality_check/checkv_filtered_fasta/all_ancient_viruses_checkv_filtered.fasta"
    threads: 1
    log:
        "logs/concat/concat.log"
    resources:
        partition="plgrid-short",
        nodes=1,
        ntasks=1,
        mem_mb="4GB",
        time="00:05:00",
    script:
        "../scripts/concat.py"
