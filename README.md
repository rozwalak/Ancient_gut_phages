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


### Introduction
This repository is a guide for analyses performed in publication about de novo assembly of ancient phage genomes.

### Downloading, preprocessing, assembly and authentication
We prepared snakemake pipeline (snake-download-preprocessing-assembly-authentication) to:
- downloading 72 Illumina metagenomic libraries of palaeofeces and human gut content samples (see: Supplementary Table S1)
- trimming paired-end reads using [Cutadapt](https://github.com/marcelm/cutadapt) (v.4.1)
- filtering out human DNA using [KneadData](https://github.com/biobakery/kneaddata) (v.0.12.0)
- assembly reads from samples into contigs using [Metaspades](https://github.com/ablab/spades) (v.3.15.5)
- filtering out contigs shorter than 4kb and with coverage <20 using custom python script
- mapping filtered reads to contigs using [Bowtie2](https://github.com/BenLangmead/bowtie2) (v.2.4.4)
- sorting and indexing with [SAMtools](https://github.com/samtools/samtools) (v.1.14)
- authenticating ancient contigs using [Pydamage](https://github.com/maxibor/pydamage) (v.0.70)
- quality control in FastQC before and after downloading, Cutadapt, KneadData


Metadata in Suplementary Table S1 is coming from [AncientMetagenomeDir](https://github.com/SPAAM-community/AncientMetagenomeDir), a community curated resource of lists of all published shotgun-sequenced ancient metagenomes.

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
We identified viral sequences in the collection of ancient contigs using three different tools:
- [Jaeger](https://github.com/Yasas1994/Jaeger) (v.1.1.0)
- [VIBRANT](https://github.com/AnantharamanLab/VIBRANT) (v.1.2.1)
- [VirSorter2](https://github.com/jiarong/VirSorter2) (v.2.2.3)

In order to predict viral sequences after installing the necessary tools, we used the combined fasta file with all of the ancient contigs and the following commands:

```
Jaeger:
python inference.py -i input_file.fasta -o jaeger_output.fasta

VIBRANT:
python VIBRANT_run.py -i input_file.fasta -t 37 -folder vibrant_output

VirSorter2:
virsorter run --prep-for-dramv -w virsorter2_output -i input_file.fasta -j 37 --include-groups dsDNAphage,NCLDV,ssDNA,lavidaviridae all
```
Finally, we selected contigs classified as viral by at least 2 methods to downstream analyses.
### Genomes quality assessment

We assessed quality of ancient viral contigs using [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) (v.1.0.1) <br>
After installation we run:
```
checkv download_database ./

checkv end_to_end input_file.fasta output_directory -t 8 -d PATH_to_DB
```
We considered ancient viral genomes classified as complete, high-quality, medium-quality, or fragments longer than 20kb. Additionally, we filtered out sequences with viral genes <= 1 and host genes >= 1 assigned by CheckV to clean putative contaminants. 
```
python checkv2fasta.py
```

### Viral contigs clustering
We clustered selected ancient viral genomes on the basis of 95% average nucleotide identity (ANI) and 85% alignment fraction of the shorter sequence, as recommended in [MIUViG](https://www.nature.com/articles/nbt.4306) (Minimum information about an uncultivated virus genome)
<br>
<br>
For this purpose, we used custom scripts published in the [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/scripts/) repository
<br>
```
#First, create a blast+ database:
makeblastdb -in <my_seqs.fna> -dbtype nucl -out <my_db>

#Next, use megablast from blast+ package to perform all-vs-all blastn of sequences:
blastn -query <my_seqs.fna> -db <my_db> -outfmt '6 std qlen slen' -max_target_seqs 10000 -o <my_blast.tsv> -num_threads 32

#Next, calculate pairwise ANI by combining local alignments between sequence pairs:
anicalc.py -i <my_blast.tsv> -o <my_ani.tsv>

#Finally, perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):
aniclust.py --fna <my_seqs.fna> --ani <my_ani.tsv> --out <my_clusters.tsv> --min_ani 95 --min_tcov 85 --min_qcov 0
```

### Gene-sharing network
At first, we predicted protein-coding genes in ancient viral genomes using [Prodigal](https://github.com/hyattpd/Prodigal) (v.2.6.3):

```
prodigal -i ancient_viruses.fasta -a ancient_viruses_proteins.faa -d ancient_viruses_proteins.ffn -p meta -f gff > ancient_viruses_proteins.gff
```

Next, we downloaded IMG_VR_2022-09-20_6.1/IMGVR_all_proteins-high_confidence.faa.gz from [IMG/VR](https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html) (v.4 high-confidence genomes only) and searched with modern viral proteins from IMG/VR as a query in [DIAMOND](https://github.com/bbuchfink/diamond) (v2.0.15):
```
gunzip IMGVR_all_proteins-high_confidence.faa.gz

diamond makedb --in ancient_viruses_proteins.faa -d ancient_viruses_proteins

diamond blastp -d ancient_viruses_proteins -q IMGVR_all_proteins-high_confidence.faa -o IMGVR_vs_all_ancient_viruses_proteins.tsv -p 96 -f 6 qseqid sseqid scovhsp pident length mismatch gapopen qstart qend sstart send evalue bitscore
```
We selected 10 modern genomes related to every ancient virus genome based on the highest number of shared proteins at a minimum of 50% query coverage
and 50% identity in Diamond output using custom script:
```
python extractTOP10.py
```
To the collection of ancient phages and related genomes from IMG/VR, we added sequences downloaded from NCBI RefSeq or GenBank using an accession number from [ICTV Virus Metadata Resources](https://ictv.global/vmr) (VMR_20-190822_MSL37.2, created 08/31/2022) and finally we made gene-sharing network in [vContact2](https://bitbucket.org/MAVERICLab/vcontact2/src/master/) (v.0.11.3)

```
#predict genes for all viruses in the network
prodigal -i SelectedIMGVR_ICTV_AncientViruses.fasta -a SelectedIMGVR_ICTV_AncientViruses.faa -p meta

#prepare file gene_2_genome
vcontact2_gene2genome -p SelectedIMGVR_ICTV_AncientViruses.faa -o SelectedIMGVR_ICTV_AncientViruses_g2g.csv -s 'Prodigal-FAA'

#make gene-sharing network
vcontact2 --raw-proteins SelectedIMGVR_ICTV_AncientViruses.faa --proteins-fp SelectedIMGVR_ICTV_AncientViruses_g2g.csv --db 'None' --output-dir SelectedIMGVR_ICTV_AncientViruses_vContact2
```
We visualised the network (Fig. 2A) from vContact2 in [Cytoscape](https://cytoscape.org/) (v.3.9.0) and refined it in [Inkscape](https://inkscape.org/) (v.1.2.2.)

### Host prediction
We used four computational tools to predict hosts of ancient viral genomes: 
- [BLASTn](https://github.com/enormandeau/ncbi_blast_tutorial)
- [PHIST](https://github.com/refresh-bio/PHIST)
- [VirHostMatcher-Net](https://github.com/WeiliWw/VirHostMatcher-Net)
- [RaFAH](https://sourceforge.net/projects/rafah/)
```
#BLAST
cat hosts/*.fna > blastdb.fna

makeblastdb -in blastdb.fna -dbtype nucl

ls aMGVs/*.fna | xargs -I {} blastn -task blastn -query {} -db blastdb.fna -outfmt 6 -num_threads 16 -out blast/{}.txt

#PHIST
./phist.py --t 16 aMGVs/ hosts/ phist/

#VirHostMatcher-Net
./VirHostMatcher-Net.py -t 16 -q aMGVs -o vhm-net/

#RaFAH
perl RaFAH.pl --predict --genomes_dir aMGVs/ --extension .fasta --file_prefix RaFAH_aMGVs
```
### Taxonomy assignment, clustering and phylogenetic analysis
We performed taxonomic assignment of ancient viral genomes using [geNomad](https://portal.nersc.gov/genomad/) (v.1.3.3)

```
genomad download-database .

genomad annotate --cleanup ancient_viruses.fasta ancient_viruses_genomad genomad_db
```
To identify genus- and family-level ancient viral genomes, we clustered genomes using a combination of gene sharing and AAI following scripts from [here](https://github.com/snayfach/MGV/tree/master/aai_cluster) (Nayfach et al., 2021)

Firstly, we selected ancient viral genomes assessed as high-quality or complete by CheckV, and we then added prokaryotic viruses from RefSeq (n = 4703; access: 30.01.2023) and IMG/VR sequences (n = 265) forming clusters (VC) with ancient viruses (see: Gene-sharing network) = ancient_viruses_RefSeq_selectedIMGVR.fasta

```
#Predict protein-coding genes
prodigal -i ancient_viruses_RefSeq_selectedIMGVR.fasta -a ancient_viruses_RefSeq_selectedIMGVR.faa -p meta

#Make DIAMOND database
diamond makedb --in ancient_viruses_RefSeq_selectedIMGVR.faa --db viral_proteins --threads 10

#Perform all-vs-all BLASTP
diamond blastp --query ancient_viruses_RefSeq_selectedIMGVR.faa --db viral_proteins --out blastp.tsv --outfmt 6 --evalue 1e-5 --max-target-seqs 10000 --query-cover 50 --subject-cover 50

#Compute AAI from BLAST results
python amino_acid_identity.py --in_faa ancient_viruses_RefSeq_selectedIMGVR.faa --in_blast blastp.tsv --out_tsv aai.tsv

#Amino acid identity is computed based on the average BLAST percent identity between all genes shared between each pair of genomes (E-value <1e-5)

#Filter edges and prepare MCL input
python filter_aai.py --in_aai aai.tsv --min_percent_shared 20 --min_num_shared 16 --min_aai 40 --out_tsv genus_edges.tsv
python filter_aai.py --in_aai aai.tsv --min_percent_shared 10 --min_num_shared 8 --min_aai 20 --out_tsv family_edges.tsv

#Here we're keeping edges between genomes with >=20% AAI and genomes with either 8 shared genes or at least 20% of shared genes (relative to both genomes)

#Perform MCL-based clustering
mcl genus_edges.tsv -te 8 -I 2.0 --abc -o genus_clusters.txt
mcl family_edges.tsv -te 8 -I 1.2 --abc -o family_clusters.txt

#In the output each row indicates the members belonging to each cluster (including singletons)
```
To visualise phylogenetic relationships of genus- and family-level groups, we generated a proteomic tree (Fig. 3) of 5017 viral sequences using [ViPTreeGen](https://github.com/yosuken/ViPTreeGen) (v1.1.3) and [GraPhlAn](https://github.com/biobakery/graphlan) (v1.1.3).
```
ViPTreeGen --ncpus 24 ancient_viruses_RefSeq_selectedIMGVR.fasta ancient_viruses_RefSeq_selectedIMGVR_VipTree
```
