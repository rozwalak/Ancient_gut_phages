import pandas as pd
from Bio import SeqIO


pydamage_csv = pd.read_csv(snakemake.input.pydamage_csv
	)
filtered_pydamage_csv = pydamage_csv[pydamage_csv['CtoT-0'] >= 0.01]

fasta = list(SeqIO.parse(snakemake.input.viruses_fasta, 'fasta')
	)
fasta = pd.DataFrame({'reference': [r.id for r in fasta], 'f_seq': fasta}
	)

final = pd.merge(fasta, filtered_pydamage_csv, on='reference'
	)
final = list(final.f_seq.values
	)

with open(snakemake.output.ancient_viruses_fasta, "w") as output_handle:
	    SeqIO.write(final, output_handle, "fasta")
