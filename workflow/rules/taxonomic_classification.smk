rule phagcn2:
    input:
        "../results/05-quality_check/checkv_filtered_fasta/{sample}_checkv_filtered.fasta",
    output:
        out="../results/06-taxonomic_classification/{sample}_taxonomy.csv"
    params:
        dir="../resources/PhaGCN2.0/*",
        new_dir="../results/06-taxonomic_classification/{sample}_PhaGCN2.0",
        db="../results/06-taxonomic_classification/{sample}_PhaGCN2.0/database/ALL_protein.tar.gz",
        untar_db="../results/06-taxonomic_classification/{sample}_PhaGCN2.0/database",
        pred="../results/06-taxonomic_classification/{sample}_PhaGCN2.0/final_prediction.csv",
        new_pred="../results/06-taxonomic_classification/{sample}_PhaGCN2.0/{sample}_taxonomy.csv",
        tmp_input="../../05-quality_check/checkv_filtered_fasta/{sample}_checkv_filtered.fasta",
        network="../results/06-taxonomic_classification/{sample}_PhaGCN2.0/network/*",
        new_network="../results/06-taxonomic_classification/{sample}_PhaGCN2.0_network",

    threads: 24
    log:
        "logs/phagcn2/{sample}.log"
    conda:
        "../envs/taxonomic_classification.yaml"
    resources:
        partition="plgrid-short",
        nodes=1,
        ntasks=24,
        mem_mb="64GB",
        time="01:00:00",
    shell: """
        export MKL_SERVICE_FORCE_INTEL=1
        mkdir -p {params.new_dir}
        cp -r -f {params.dir} {params.new_dir}

        tar -zxvf {params.db} -C {params.untar_db}
        cd {params.new_dir}
        python run_Speed_up.py --contigs {params.tmp_input} --len 2000
        cd ..
        cd ..
        cd ..
        cd workflow
        mv {params.pred} {params.new_pred}
        cp -f {params.new_pred} {output.out}
        mkdir -p {params.new_network}
        cp -f {params.network} {params.new_network}
        rm -r -f {params.new_dir}
        """
