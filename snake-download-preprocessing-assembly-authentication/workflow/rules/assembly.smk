rule metaspades:
    input:
        reads1="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.1_trimmed_kneaddata_paired_1.fastq",
        reads2="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.1_trimmed_kneaddata_paired_2.fastq",
    output:
        filtered="../results/02-assembly/{sample}/{sample}_scaffolds_filtered_4k_20cov.fasta"
    log:
        "logs/metaspades/{sample}_metaspades.log",
    threads: 16
    params:
        meta_dir=directory("../results/02-assembly/{sample}"),
        scaffolds="../results/02-assembly/{sample}/scaffolds.fasta",
        new_dir="../results/02-assembly/output_scaffolds",
    conda:
        "../envs/metaspades.yaml"
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=16,
        mem_mb="120GB",
        disk_mb="80GB",
        time="01:00:00",
    shell: """
        metaspades.py \
	    --meta \
	    -1 {input.reads1} \
	    -2 {input.reads2} \
	    -o {params.meta_dir} \
	    -t {threads} \
	    -m 120

        python scripts/assembly_filter.py {params.scaffolds} {output.filtered}
        mkdir -p {params.new_dir}
        cp {output.filtered} {params.new_dir}
        """
