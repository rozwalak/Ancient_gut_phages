import pandas as pd
from Bio import SeqIO

# Load output from Diamond
print("Loading output from Diamond...")
col_headers_list = ["qseqid", "sseqid", "scovhsp", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
df = pd.read_csv("IMGVR_vs_all_ancient_viruses_proteins.tsv", sep='\t', header=None, names=col_headers_list)

# Remove '_id' for each ancient genome name
print("Removing '_id' for each ancient genome name...")
df['sseqid'] = df['sseqid'].apply(lambda r: '_'.join(r.split('_')[:-1]))

# Filter sequences with at minimum qcovs >= 50 and pident >= 50
print("Filtering sequences with at minimum qcovs >= 50 and pident >= 50...")
df_filtered = df[(df['scovhsp'] >= 50)  & (df['pident'] >= 50)]
df_filtered = df_filtered.groupby(["qseqid", "sseqid"]).size().reset_index(name="frequency")

# Groupby and sort shared proteins between ancient and modern viral genomes
print("Grouping by and sorting shared proteins between ancient and modern viral genomes...")
df_sorted = df_filtered.groupby(["sseqid"]).apply(lambda x: x.sort_values(["frequency"], ascending = False)).reset_index(drop=True)

# Select top 10 modern genomes for each ancient genome
print("Selecting top 10 modern genomes for each ancient genome...")
df_sorted = df_sorted.groupby('sseqid').head(10)

# Drop duplicates in modern genomes
print("Dropping duplicates in modern genomes...")
df_sorted = df_sorted.drop_duplicates(subset=['qseqid'])

# Save to csv
print("Saving to CSV...")
df_sorted.to_csv("TOP10_IMGVR_vs_all_ancient_viruses_proteins.csv")

# Read fasta with modern viral genomes
print("Reading FASTA with modern viral genomes...")
db_fasta = list(SeqIO.parse('IMGVR_v4.fasta', 'fasta'))
db_fasta = pd.DataFrame({'qseqid': [r.id for r in db_fasta], 'f_seq': db_fasta})

# Read fasta with ancient viral genomes
print("Reading FASTA with ancient viral genomes...")
viruses_fasta = list(SeqIO.parse('ancient_viruses.fasta', 'fasta'))
viruses_fasta = pd.DataFrame({'qseqid': [r.id for r in viruses_fasta], 'f_seq': viruses_fasta})

# Create fasta file with selected modern and ancient genomes
print("Creating FASTA file with selected modern and ancient genomes...")
db_fasta_selected = pd.merge(db_fasta, df_sorted, on='qseqid')
final = pd.concat([viruses_fasta, db_fasta_selected])
final = list(final.f_seq.values)

# Write final fasta file
print("Writing final FASTA file...")
with open('TOP10_IMGVR_ancient_viruses.fasta', "w") as output_handle:
    SeqIO.write(final, output_handle, "fasta")