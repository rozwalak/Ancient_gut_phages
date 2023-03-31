import pandas as pd
from Bio import SeqIO


pydamage_csv = pd.read_csv(snakemake.input.checkv_summary, sep="\t")

checkv_filtered = pydamage_csv[(pydamage_csv['checkv_quality'] == "Complete") | (pydamage_csv['checkv_quality'] == "High-quality") | (pydamage_csv["checkv_quality"] == "Medium-quality") | (pydamage_csv["contig_length"] >= 20000)]
checkv_filtered = checkv_filtered.drop(checkv_filtered[checkv_filtered['viral_genes'] <= 1 & (checkv_filtered['host_genes'] >= 1)].index)
checkv_filtered["run_accession"] = snakemake.params.sample_name

checkv_filtered.to_csv(snakemake.output.checkv_filtered_csv, sep=",", header=True, index=False)

fasta = list(SeqIO.parse(snakemake.input.viruses_fasta, 'fasta'))
fasta = pd.DataFrame({'contig_id': [r.id for r in fasta], 'f_seq': fasta})

final = pd.merge(fasta, checkv_filtered, on='contig_id')
final = list(final.f_seq.values)

with open(snakemake.output.checkv_filtered_fasta, "w") as output_handle:
	    SeqIO.write(final, output_handle, "fasta")
