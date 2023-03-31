rule all_vs_all_blast:
    input:
        checkv_summary="../results/05-quality_check/{sample}_checkv/quality_summary.tsv",
        viruses_fasta="../results/04-aDNA_authentication/ancient_viruses_fasta/{sample}_ancient_viruses.fasta"
    output:
        checkv_filtered_csv="../results/05-quality_check/checkv_filtered_csv/{sample}_checkv_filtered.csv",
        checkv_filtered_fasta="../results/05-quality_check/checkv_filtered_fasta/{sample}_checkv_filtered.fasta"
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



rule blast_ani:


rule cluster:
