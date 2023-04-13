import pandas as pd
from Bio import SeqIO

#input_file = change to your name

checkv_tsv = pd.read_csv('quality_summary.tsv', sep="\t")

checkv_filtered = pydamage_csv[(checkv_tsv['checkv_quality'] == "Complete") | (checkv_tsv['checkv_quality'] == "High-quality") | (checkv_tsv["checkv_quality"] == "Medium-quality") | (checkv_tsv["contig_length"] >= 20000)]
checkv_filtered = checkv_filtered.drop(checkv_filtered[checkv_filtered['viral_genes'] <= 1 & (checkv_filtered['host_genes'] >= 1)].index)


checkv_filtered.to_csv("checkv_filtered.csv", sep=",", header=True, index=False)

fasta = list(SeqIO.parse("input_file.fasta", 'fasta'))
fasta = pd.DataFrame({'contig_id': [r.id for r in fasta], 'f_seq': fasta})

final = pd.merge(fasta, checkv_filtered, on='contig_id')
final = list(final.f_seq.values)

with open("input_file_checkv_filtered", "w") as output_handle:
	    SeqIO.write(final, output_handle, "fasta")
