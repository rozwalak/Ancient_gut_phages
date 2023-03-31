# rule jaeger:
#     input:
#         ""
#     output:
#         ""
#     shell:
#---------------------------------------------------------------------------------------------------------------------------------------------------
rule jaeger2fasta:
    input:
        predicted="../results/03-viral_prediction/jaeger/{sample}_scaffolds.jaeger.txt",
        fasta="../results/02-assembly/output_scaffolds/{sample}_scaffolds.fasta",
    output:
        viruses="../results/03-viral_prediction/viruses/{sample}_viruses.txt",
        viruses_fasta="../results/03-viral_prediction/viruses_fasta/{sample}_viruses.fasta",
    params:
        name="{sample}"
    threads: 1
    resources:
        partition="plgrid-short",
        nodes=1,
        ntasks=1,
        mem_mb="4GB",
        time="00:05:00"
    script:
        "../scripts/jaeger2fasta.py"
#---------------------------------------------------------------------------------------------------------------------------------------------------
