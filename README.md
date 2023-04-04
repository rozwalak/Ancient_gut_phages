# Reconstruction of ancient phage genomes from the human gut

## Table of contents
- [Introduction](#introduction)
- [Downloading, preprocessing and assembly](#downloading-preprocessing-and-assembly)
- [aDNA authentication](#adna-authentication)
- [Viral contigs identification](#viral-contigs-identification)
- [Viral contigs clustering](#viral-contigs-clustering)
- [Genomes quality assessment](#genomes-quality-assessment)
- [Gene-sharing network](#gene-sharing-network)
- [Host prediction](#host-prediction)
- [Taxonomy assignment, clustering and phylogenetic analysis](#taxonomy-assignment-clustering-and-phylogenetic-analysis)
- [Analyses of Mushuvirus mushu genome](#analyses-of-mushuvirus-mushu-genome)

### Introduction
This repository is a guide for analyses from publication about de novo assembly of ancient phage genomes.

### Downloading, preprocessing and assembly
We prepared snakemake pipeline (snake-download-preprocessing-assembly) to:
- downloading 72 Illumina metagenomic libraries of palaeofeces and human gut content samples (see: supplementary table S1)
- trimming paired-end reads using Cutadapt (v.4.1)
- filtering out human DNA using KneadData (v.0.12.0)
- assembly reads from samples into contigs using Metaspades (v.3.15.5)

The first step is installation snakemake environment following instructions from [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)<br>  
Then you can download this repository and run our pipeline which automatically create specific environments using conda and download all 72 libraries. 
You can modify config.yaml from config folder to select different samples to analyze. 
<br>  

Metadata in Suplementary Table S1 coming from [AncientMetagenomeDir](https://github.com/SPAAM-community/AncientMetagenomeDir), a community curated resource of lists of all published shotgun-sequenced ancient metagenomes.

### aDNA authentication

### Viral contigs identification

### Viral contigs clustering

### Genomes quality assessment

### Gene-sharing network

### Host prediction

### Taxonomy assignment, clustering and phylogenetic analysis

### Analyses of Mushuvirus mushu genome
