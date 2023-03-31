import pandas as pd
from Bio import SeqIO


df = pd.read_csv(snakemake.input.predicted, sep="\t")
df_viruses = df[df['prediction'] == "Virus"]
df_viruses.to_csv(snakemake.output.viruses, sep="\t", index=False)

fasta = list(SeqIO.parse(snakemake.input.fasta, "fasta"))
fasta = pd.DataFrame({"contig_id": [r.id for r in fasta], "f_seq": fasta})

fasta_viruses = pd.merge(fasta, df_viruses, on='contig_id')
fasta_viruses = list(fasta_viruses.f_seq.values)

with open(snakemake.output.viruses_fasta, "w") as output_handle:
        SeqIO.write(fasta_viruses, output_handle, "fasta")
