from Bio import SeqIO

input_file = snakemake.params.scaffolds
output_file = snakemake.output.filtered
min_length = 4000
min_coverage = 20

with open(output_file, "w") as output_handle:
    for record in SeqIO.parse(input_file, "fasta"):
        node_id, length, cov = record.id.split("_")
        if int(length) >= min_length and float(cov) >= min_coverage:
            SeqIO.write(record, output_handle, "fasta")