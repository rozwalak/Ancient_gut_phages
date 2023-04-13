# Reconstruction of ancient phage genomes from the human gut

## Table of contents
- [Introduction](#introduction)
- [Downloading, preprocessing, assembly and authentication](#downloading-preprocessing-assembly-and-authentication)
- [Viral contigs identification](#viral-contigs-identification)
- [Viral contigs clustering](#viral-contigs-clustering)
- [Genomes quality assessment](#genomes-quality-assessment)
- [Gene-sharing network](#gene-sharing-network)
- [Host prediction](#host-prediction)
- [Taxonomy assignment, clustering and phylogenetic analysis](#taxonomy-assignment-clustering-and-phylogenetic-analysis)
- [Analyses of Mushuvirus mushu genome](#analyses-of-mushuvirus-mushu-genome)

### Introduction
This repository is a guide for analyses from publication about de novo assembly of ancient phage genomes.

### Downloading, preprocessing, assembly and authentication
We prepared snakemake pipeline (snake-download-preprocessing-assembly-authentication) to:
- downloading 72 Illumina metagenomic libraries of palaeofeces and human gut content samples (see: supplementary table S1)
- trimming paired-end reads using [Cutadapt](https://github.com/marcelm/cutadapt) (v.4.1)
- filtering out human DNA using [KneadData](https://github.com/biobakery/kneaddata) (v.0.12.0)
- assembly reads from samples into contigs using [Metaspades](https://github.com/ablab/spades) (v.3.15.5)
- filtering out contigs shorter than 4kb and with coverage <20 using custom python script
- mapping filtered reads to contigs using [Bowtie2](https://github.com/BenLangmead/bowtie2) (v.2.4.4)
- sorting and indexing with [SAMtools](https://github.com/samtools/samtools) (v.1.14)
- authenticating ancient contigs using [Pydamage](https://github.com/maxibor/pydamage) (v.0.70)
- quality control in FastQC before and after downloading, Cutadapt, KneadData


Metadata in suplementary table S1 is coming from [AncientMetagenomeDir](https://github.com/SPAAM-community/AncientMetagenomeDir), a community curated resource of lists of all published shotgun-sequenced ancient metagenomes.

The first step is installation snakemake following instructions from [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).<br>  
Then you can download this repository and run our pipeline which automatically create specific environments using conda and download all 72 libraries. 
You can modify config.yaml from config folder to select different samples to analyze. 
<br>  

```
conda activate snakemake

git clone https://github.com/rozwalak/Ancient_gut_phages.git

cd Ancient_gut_phages/snake-download-preprocessing-assembly-authentication/workflow

snakemake --profile ../config/snakemake/slurm --use-conda 
```

### Viral contigs identification
We identified viral sequences in collection of ancient contigs using three different tools:
- [Jaeger](https://github.com/Yasas1994/Jaeger)
- [VIBRANT](https://github.com/AnantharamanLab/VIBRANT) (v.1.2.1)
- [VirSorter2](https://github.com/jiarong/VirSorter2) (v.2.2.3)

After installation of tools we performed predictions using following commands for combined fasta file with all ancient contigs:

```
Jaeger:
python inference.py -i input_file.fasta -o jaeger_output.fasta

VIBRANT:
python VIBRANT_run.py -i input_file.fasta -t 37 -folder vibrant_output

VirSorter2:
virsorter run --prep-for-dramv -w virsorter2_output -i input_file.fasta -j 37 --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae all
```
Finally, we selected contigs classified as viral by at least 2 methods to further analyses.
### Genomes quality assessment

We assessed quality of ancient viral contigs using [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) (v.1.0.1) <br>
After installation we run:
```
checkv download_database ./

checkv end_to_end input_file.fasta output_directory -t 8 -d PATH_to_DB
```
We considered ancient viral genomes classified as complete, high-quality, medium-quality, or fragments longer than 20kb. Additionally, we filtered out sequences with viral genes <= 1 and host genes >= 1 to clean putative contaminants. 
### Viral contigs clustering
We clustered selected ancient viral genomes on the basis of 95% average nucleotide identity (ANI) and 85% alignment fraction of the shorter sequence, as recommended in [MIUViG](https://www.nature.com/articles/nbt.4306) (Minimum information about an uncultivated virus genome)
<br>
<br>
For this purpose, we used custom scripts published in the [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/scripts/) repository
<br>
```
First, create a blast+ database:
makeblastdb -in <my_seqs.fna> -dbtype nucl -out <my_db>

Next, use megablast from blast+ package to perform all-vs-all blastn of sequences:
blastn -query <my_seqs.fna> -db <my_db> -outfmt '6 std qlen slen' -max_target_seqs 10000 -o <my_blast.tsv> -num_threads 32

Next, calculate pairwise ANI by combining local alignments between sequence pairs:
anicalc.py -i <my_blast.tsv> -o <my_ani.tsv>

Finally, perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):
aniclust.py --fna <my_seqs.fna> --ani <my_ani.tsv> --out <my_clusters.tsv> --min_ani 95 --min_tcov 85 --min_qcov 0
```

### Gene-sharing network

### Host prediction

### Taxonomy assignment, clustering and phylogenetic analysis

### Analyses of Mushuvirus mushu genome
