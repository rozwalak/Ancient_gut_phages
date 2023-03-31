rule phagcn2:
    input:
        "../results/05-quality_check/checkv_filtered_fasta/all_ancient_viruses_checkv_filtered.fasta",
    output:
        "../results/06-taxonomic_classification/all_ancient_viruses_taxonomy.csv"
    params:
        input="../../results/05-quality_check/checkv_filtered_fasta/all_ancient_viruses.fasta",
        output="../../results/06-taxonomic_classification/"
    threads: 24
    log:
        "logs/phagcn2/phagcn2.log"  #naprawić te logi! bo terazjest problem z tym, że nie ma job indywidualnej nazwy
    conda:
        "../envs/taxonomic_classification.yaml"
    resources:
        partition="plgrid",
        nodes=1,
        ntasks=24,
        mem_mb="64GB",
        time="01:30:00",
    shell: """
        export MKL_SERVICE_FORCE_INTEL=1

        cd ../resources/PhaGCN2.0
        cd database
        tar -zxvf ALL_protein.tar.gz
        cd ..

        python run_Speed_up.py --contigs {params.input} --len 2000

        mv final_prediction.csv all_ancient_viruses_taxonomy.csv
        cp all_ancient_viruses_taxonomy.csv {params.output}
        cp -r network {params.output}

        """
