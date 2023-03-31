rule metaspades:
    input:
        reads1="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.1_trimmed_kneaddata_paired_1.fastq",
        reads2="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.1_trimmed_kneaddata_paired_2.fastq",
    output:
        meta_dir=directory("../results/02-assembly/{sample}"),
    benchmark:
        "logs/benchmarks/metaspades/{sample}_metaspades.txt"
    log:
        "logs/metaspades/{sample}_metaspades.log",
    threads: 16
    params:
        scaffolds="../results/02-assembly/{sample}/scaffolds.fasta",
        new_name="../results/02-assembly/{sample}/{sample}_scaffolds.fasta",
        new_dir="../results/02-assembly/output_scaffolds",
    envmodules:
        "plgrid/tools/spades/3.15.4"
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=16,
        mem_mb="120GB",
        disk_mb="80GB",
        time="10:00:00",
    shell: """
        metaspades.py \
	    --meta \
	    -1 {input.reads1} \
	    -2 {input.reads2} \
	    -o {output.meta_dir} \
	    -t {threads} \
	    -m 120

        mv {params.scaffolds} {params.new_name}
        mkdir -p {params.new_dir}
        cp {params.new_name} {params.new_dir}
        """
