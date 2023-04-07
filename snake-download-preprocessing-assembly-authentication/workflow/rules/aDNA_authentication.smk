rule bowtie2_build:
    input:
        "../results/02-assembly/{sample}/{sample}_scaffolds_filtered_4k_20cov.fasta"
    params:
        basename="../results/04-aDNA_authentication/indexes/{sample}_index/{sample}.index"
    output:
        output1="../results/04-aDNA_authentication/indexes/{sample}_index/{sample}.index.1.bt2",
        output2="../results/04-aDNA_authentication/indexes/{sample}_index/{sample}.index.2.bt2",
        output3="../results/04-aDNA_authentication/indexes/{sample}_index/{sample}.index.3.bt2",
        output4="../results/04-aDNA_authentication/indexes/{sample}_index/{sample}.index.4.bt2",
        outputrev1="../results/04-aDNA_authentication/indexes/{sample}_index/{sample}.index.rev.1.bt2",
        outputrev2="../results/04-aDNA_authentication/indexes/{sample}_index/{sample}.index.rev.2.bt2"
    log:
        "logs/bowtie2_build/{sample}.log",
    conda:
        "../envs/aDNA_authentication_bowtie2_samtools.yaml"
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=2,
        mem_mb="64GB",
        time="01:00:00",
    shell: "bowtie2-build {input} {params.basename}"

#---------------------------------------------------------------------------------------------------------------------------------------------------
rule bowtie2:
    input:
        r1="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.1_trimmed_kneaddata_paired_1.fastq",
        r2="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.1_trimmed_kneaddata_paired_2.fastq",
        input1="../results/04-aDNA_authentication/indexes/{sample}_index/{sample}.index.1.bt2",
        input2="../results/04-aDNA_authentication/indexes/{sample}_index/{sample}.index.2.bt2",
        input3="../results/04-aDNA_authentication/indexes/{sample}_index/{sample}.index.3.bt2",
        input4="../results/04-aDNA_authentication/indexes/{sample}_index/{sample}.index.4.bt2",
        inputrev1="../results/04-aDNA_authentication/indexes/{sample}_index/{sample}.index.rev.1.bt2",
        inputrev2="../results/04-aDNA_authentication/indexes/{sample}_index/{sample}.index.rev.2.bt2"
    output:
        bam="../results/04-aDNA_authentication/bowtie2_out/{sample}.bam",
        bam_bai="../results/04-aDNA_authentication/bowtie2_out/{sample}.bam.bai"
    threads: 24
    log:
        "logs/bowtie2/{sample}.log"
    conda:
        "../envs/aDNA_authentication_bowtie2_samtools.yaml"
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=24,
        mem_mb="64GB",
        time="01:00:00",
    shell: """
        bowtie2 -1 {input.r1} -2 {input.r2} -x ../results/04-aDNA_authentication/indexes/{wildcards.sample}_index/{wildcards.sample}.index -p {threads} | samtools view -@ {threads} -Sb | samtools sort -@ {threads} > {output.bam}
        samtools index {output.bam} {output.bam_bai}
        """
#---------------------------------------------------------------------------------------------------------------------------------------------------

rule pydamage:
    input:
        "../results/04-aDNA_authentication/bowtie2_out/{sample}.bam"
    output:
        outdir=directory("../results/04-aDNA_authentication/pydamage/{sample}_pydamage"),
        pydamage="../results/04-aDNA_authentication/pydamage/{sample}_pydamage/pydamage_results.csv",
        pydamage_filtered="../results/04-aDNA_authentication/pydamage/{sample}_pydamage/{sample}_pydamage_filtered_results.csv",
    threads: 24
    log:
        "logs/pydamage/{sample}.log"
    conda:
        "../envs/aDNA_authentication_pydamage.yaml"
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=24,
        mem_mb="64GB",
        time="01:00:00",
    shell: """

        pydamage \
        --outdir {output.outdir} \
        analyze {input} \
        --force \
        --process={threads}

        pydamage filter -t 0 {output.pydamage}

        mv pydamage_results/pydamage_filtered_results.csv {output.pydamage_filtered}
        rm -r pydamage_results
        """
#---------------------------------------------------------------------------------------------------------------------------------------------------
rule pydamage2fasta:
    input:
        pydamage_csv="../results/04-aDNA_authentication/pydamage/{sample}_pydamage/{sample}_pydamage_filtered_results.csv",
        viruses_fasta="../results/02-assembly/{sample}/{sample}_scaffolds_filtered_4k_20cov.fasta",
    output:
        ancient_viruses_fasta="../results/04-aDNA_authentication/ancient_contigs_fasta/{sample}_ancient.fasta",
    params:
        all_pydamage_outputs=directory("../results/04-aDNA_authentication/pydamage/pydamage_all_outputs"),
        sample_name="{sample}"
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
