rule metaspades:
    input:
        reads1="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.1_trimmed_kneaddata_paired_1.fastq",
        reads2="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.1_trimmed_kneaddata_paired_2.fastq",
    output:
        meta_dir=directory("../results/02-assembly/{sample}"),
        new_name="../results/02-assembly/{sample}/{sample}_scaffolds.fasta"
    benchmark:
        "logs/benchmarks/metaspades/{sample}_metaspades.txt"
    log:
        "logs/metaspades/{sample}_metaspades.log",
    threads: 16
    params:
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
	    -o {output.meta_dir} \
	    -t {threads} \
	    -m 120

        mkdir -p {params.new_dir}
        mv {params.scaffolds} {output.new_name}
        cp {output.new_name} {params.new_dir}
        """
rule assembly_filter:
    input:
        "../results/02-assembly/{sample}/{sample}_scaffolds.fasta",
    output:
        "../results/02-assembly/{sample}/{sample}_filtered_scaffolds.fasta",
    threads: 1
    log:
        "logs/pydamage2fasta/{sample}.log"
    conda:
        "../envs/aDNA_authentication_pydamage.yaml"
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=1,
        mem_mb="4GB",
        time="01:00:00",
    script:
        "../scripts/pydamage2fasta.py"
