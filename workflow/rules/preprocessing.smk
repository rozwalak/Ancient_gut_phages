rule fastqc_raw_data_1:
    input:
        "../results/00-raw_data/{sample}_1.fastq.gz"
    output:
        html="../results/01-preprocessing/01-raw_data-fastqc/{sample}.1_raw_data.html",
        zip="../results/01-preprocessing/01-raw_data-fastqc/{sample}.1_raw_data_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        "logs/fastqc_raw_data_1/{sample}.log"
    threads: 2
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=2,
        time="00:30:00"
    wrapper:
        "v1.25.0/bio/fastqc"

rule fastqc_raw_data_2:
    input:
        "../results/00-raw_data/{sample}_2.fastq.gz"
    output:
        html="../results/01-preprocessing/01-raw_data-fastqc/{sample}.2_raw_data.html",
        zip="../results/01-preprocessing/01-raw_data-fastqc/{sample}.2_raw_data_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        "logs/fastqc_raw_data_2/{sample}.log"
    threads: 1
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=1,
        time="00:30:00"
    wrapper:
        "v1.25.0/bio/fastqc"

#---------------------------------------------------------------------------------------------------------------------------------------------------
rule cutadapt:
    input:
        ["../results/00-raw_data/{sample}_1.fastq.gz", "../results/00-raw_data/{sample}_2.fastq.gz"]
    output:
        fastq1="../results/01-preprocessing/02-cutadapt/{sample}.1_trimmed.fastq.gz",
        fastq2="../results/01-preprocessing/02-cutadapt/{sample}.2_trimmed.fastq.gz",
        qc="../results/01-preprocessing/02-cutadapt/{sample}.qc.txt"
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters="-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        extra="-q 25 -m 30 -O 1 -j 0"
    log:
        "logs/cutadapt/{sample}.log"
    threads: 14
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=14,
        time="01:00:00"
    envmodules:
        "plgrid/tools/cutadapt/3.7"
    wrapper:
        "v1.25.0/bio/cutadapt/pe"

#---------------------------------------------------------------------------------------------------------------------------------------------------
rule fastqc_cutadapt_1:
    input:
        "../results/01-preprocessing/02-cutadapt/{sample}.1_trimmed.fastq.gz"
    output:
        html="../results/01-preprocessing/03-cutadapt_output-fastqc/{sample}.1_cutadapt_output.html",
        zip="../results/01-preprocessing/03-cutadapt_output-fastqc/{sample}.1_cutadapt_output_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        "logs/fastqc_cutadapt_1/{sample}.log"
    threads: 1
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=1,
        time="00:30:00"
    wrapper:
        "v1.25.0/bio/fastqc"

rule fastqc_cutadapt_2:
    input:
        "../results/01-preprocessing/02-cutadapt/{sample}.2_trimmed.fastq.gz"
    output:
        html="../results/01-preprocessing/03-cutadapt_output-fastqc/{sample}.2_cutadapt_output.html",
        zip="../results/01-preprocessing/03-cutadapt_output-fastqc/{sample}.2_cutadapt_output_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        "logs/fastqc_cutadapt_2/{sample}.log"
    threads: 1
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=1,
        time="00:30:00"
    wrapper:
        "v1.25.0/bio/fastqc"

#---------------------------------------------------------------------------------------------------------------------------------------------------
rule kneaddata_database:
    output:
        directory("../resources/human_database")
    conda:
        "../envs/kneaddata.yaml"
    log:
        "logs/kneaddata_database/human_database.log"
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=1,
        time="01:00:00"
    shell:
        "kneaddata_database --download human_genome bowtie2 {output}"

#---------------------------------------------------------------------------------------------------------------------------------------------------
rule seqtk:
    output:
        tmp_fwd2="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.tmp2_1.fastq",
        tmp_rev2="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.tmp2_2.fastq",
    input:
        fwd="../results/01-preprocessing/02-cutadapt/{sample}.1_trimmed.fastq.gz",
        rev="../results/01-preprocessing/02-cutadapt/{sample}.2_trimmed.fastq.gz",
    params:
        tmp_fwd="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.tmp_1.fastq",
        tmp_rev="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.tmp_2.fastq",
    threads: 2
    conda:
        "../envs/seqtk_bbmap.yaml"
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=2,
        mem_mb="12GB",
        time="00:30:00"
    log:
        "logs/seqtk/{sample}.log"
    shell: """
        seqtk seq -C {input.fwd} > {params.tmp_fwd}
        seqtk seq -C {input.rev} > {params.tmp_rev}
        reformat.sh in={params.tmp_fwd} in2={params.tmp_rev} out1={output.tmp_fwd2} out2={output.tmp_rev2} addslash spaceslash=f
        """

#---------------------------------------------------------------------------------------------------------------------------------------------------
rule kneaddata:
    output:
        fwd="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.tmp2_1_kneaddata_paired_1.fastq",
        rev="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.tmp2_1_kneaddata_paired_2.fastq",
    input:
        tmp_fwd2="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.tmp2_1.fastq",
        tmp_rev2="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.tmp2_2.fastq",
    # params:
        # out_dir=directory("../results/01-preprocessing/04-kneaddata/{sample}")
    params:
        outdir=directory("../results/01-preprocessing/04-kneaddata/{sample}"),
        indx="../../resources/human_database",
    conda:
        "../envs/kneaddata.yaml"
    threads: 24
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=24,
        mem_mb="94GB",
        time="06:00:00" #06:00:00 for big samples
    log:
        "logs/kneaddata/{sample}.log"

    shell: """
        kneaddata --threads {threads} --processes 24 \
        -i1 {input.tmp_fwd2} -i2 {input.tmp_rev2}\
        --output {params.outdir} \
        -db {params.indx} \
        --bypass-trim
            """
#---------------------------------------------------------------------------------------------------------------------------------------------------
rule repair:
    output:
        fwd="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.1_trimmed_kneaddata_paired_1.fastq",
        rev="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.1_trimmed_kneaddata_paired_2.fastq",
    input:
        fwd="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.tmp2_1_kneaddata_paired_1.fastq",
        rev="../results/01-preprocessing/04-kneaddata/{sample}/{sample}.tmp2_1_kneaddata_paired_2.fastq",
    threads: 12
    conda:
        "../envs/seqtk_bbmap.yaml"
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=12,
        mem_mb="32GB",
        time="00:30:00"
    shell: """
        repair.sh in={input.fwd} in2={input.rev} out={output.fwd} out2={output.rev} repair
        """

#---------------------------------------------------------------------------------------------------------------------------------------------------
rule fastqc_kneaddata_1:
    input:
        "../results/01-preprocessing/04-kneaddata/{sample}/{sample}.1_trimmed_kneaddata_paired_1.fastq"
    output:
        html="../results/01-preprocessing/05-kneaddata_output-fastqc/{sample}.1_kneaddata_output.html",
        zip="../results/01-preprocessing/05-kneaddata_output-fastqc/{sample}.1_kneaddata_output_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        "logs/fastqc_kneaddata_1/{sample}.log"
    threads: 2
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=2,
        time="00:30:00"
    wrapper:
        "v1.25.0/bio/fastqc"

rule fastqc_kneaddata_2:
    input:
        "../results/01-preprocessing/04-kneaddata/{sample}/{sample}.1_trimmed_kneaddata_paired_2.fastq"
    output:
        html="../results/01-preprocessing/05-kneaddata_output-fastqc/{sample}.2_kneaddata_output.html",
        zip="../results/01-preprocessing/05-kneaddata_output-fastqc/{sample}.2_kneaddata_output_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        "logs/fastqc_kneaddata_2/{sample}.log"
    threads: 1
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=1,
        time="00:30:00"
    wrapper:
        "v1.25.0/bio/fastqc"
