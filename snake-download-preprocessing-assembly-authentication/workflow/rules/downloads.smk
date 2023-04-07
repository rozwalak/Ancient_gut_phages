rule download_raw_data:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        "../results/00-raw_data/{accession}_1.fastq.gz",
        "../results/00-raw_data/{accession}_2.fastq.gz",
    log:
        "logs/download_raw_data/{accession}.gz.log"
    params:
        extra="--skip-technical"
    threads: 2
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=2,
        time="03:00:00"
    wrapper:
        "v1.7.0/bio/sra-tools/fasterq-dump"
